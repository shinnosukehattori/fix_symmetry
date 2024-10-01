/*
  * fix_symmetry.cpp
  *
  *  Created on: 2024/09/29 shattori
*/

#include "fix_symmetry.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "update.h"
#include "comm.h"
#include "memory.h"
#include "math_extra.h"
#include <cmath>
#include <iostream>

using namespace LAMMPS_NS;


FixSymmetry::FixSymmetry(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (narg < 4 || narg > 8)
    error->all(FLERR, "Illegal fix symmetry command. Must specify space group number and optionals, symprec, symforce,  symstress and debug-mode.");

  // Get the space group number
  spacegroup_number = atoi(arg[3]);

  // Get the tolerance (if not specified, set default value)
  if (narg >= 5) {
    symprec = utils::numeric(FLERR, arg[4], false, lmp);
  } else {
    symprec = 1e-5;  // Default tolerance
  }
  if (narg >= 6) {
    symforce = utils::logical(FLERR, arg[5], false, lmp);
  } else {
    symforce = false;
  }
  if (narg >= 7) {
    symforce = utils::logical(FLERR, arg[6], false, lmp);
  } else {
    symstress = false;
  }
  if (narg >= 8) {
    debug = utils::logical(FLERR, arg[7], false, lmp);
  } else {
    debug = false;
  }


  std_cell[0][0] = std_cell[0][1] = std_cell[0][2] = 0.0;
  std_cell[1][0] = std_cell[1][1] = std_cell[1][2] = 0.0;
  std_cell[2][0] = std_cell[2][1] = std_cell[2][2] = 0.0;
  prim_cell[0][0] = prim_cell[0][1] = prim_cell[0][2] = 0.0;
  prim_cell[1][0] = prim_cell[1][1] = prim_cell[1][2] = 0.0;
  prim_cell[2][0] = prim_cell[2][1] = prim_cell[2][2] = 0.0;


  // Initialize previous positions array
  dataset = nullptr;
  prev_positions = nullptr;
  all_positions = nullptr;
  all_types = nullptr;

}

FixSymmetry::~FixSymmetry() {
  // Free memory
  if (prev_positions) {
    memory->destroy(prev_positions);
  }
  if (all_positions) {
    memory->destroy(all_positions);
  }
  if (all_types) {
    memory->destroy(all_types);
  }
  if (dataset) {
    spg_free_dataset(dataset);
  }
}

int FixSymmetry::setmask() {
  int mask = 0;
  mask |= FixConst::POST_FORCE;
  mask |= FixConst::MIN_POST_FORCE;
  return mask;
}

void FixSymmetry::init() {
  // Save previous atom positions

  initial_prep();
  refine_symmetry();
  prep_symmetry();

  memory->create(prev_positions, atom->nlocal, 3, "fix_symmetry:prev_positions");
  for (int i = 0; i < atom->nlocal; i++) {
    prev_positions[i][0] = atom->x[i][0];
    prev_positions[i][1] = atom->x[i][1];
    prev_positions[i][2] = atom->x[i][2];
  }
}

void FixSymmetry::post_force(int vflag) {
  // Monitor atom displacement and update symmetry operations if needed
  // if (need_to_update_symmetry()) {
  if (1) {
    adjust_cell();
    adjust_positions();
    if (symforce) {
      adjust_forces();
    }
    if (symstress) {
      adjust_stress();
    }
  }
}

void FixSymmetry::min_post_force(int vflag) {
  post_force(vflag);
}

bool FixSymmetry::need_to_update_symmetry() {
  double max_displacement = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    double dx = atom->x[i][0] - prev_positions[i][0];
    double dy = atom->x[i][1] - prev_positions[i][1];
    double dz = atom->x[i][2] - prev_positions[i][2];
    double disp = sqrt(dx*dx + dy*dy + dz*dz);
    if (disp > max_displacement) {
      max_displacement = disp;
    }
    // Update previous positions
    prev_positions[i][0] = atom->x[i][0];
    prev_positions[i][1] = atom->x[i][1];
    prev_positions[i][2] = atom->x[i][2];
  }

  // Compare maximum displacement across all processes
  double global_max_disp;
  MPI_Allreduce(&max_displacement, &global_max_disp, 1, MPI_DOUBLE, MPI_MAX, world);

  // Determine if symmetry operations need to be updated
  printf("global_max_disp = %f\n", global_max_disp);
  if (global_max_disp > symprec) {
    return true;  // Need to recalculate symmetry operations
  }
  return false;  // No need to recalculate
}

int FixSymmetry::get_index(std::vector<int> &vec, int val) {
  for (int i = 0; i < vec.size(); i++) {
    if (vec[i] == val) {
      return i;
    }
  }
  return -1;
}

void FixSymmetry::store_std_cell() {
  std_cell[0][1] = std_cell[0][2] = std_cell[1][2] = 0.0;
  std_cell[0][0] = domain->h[0];
  std_cell[1][1] = domain->h[1];
  std_cell[2][2] = domain->h[2];
  if (domain->triclinic) {
    std_cell[1][0] = std_cell[2][0] = std_cell[2][1] = 0.0;
  } else {
    std_cell[1][0] = domain->xy;
    std_cell[2][0] = domain->xz;
    std_cell[2][1] = domain->yz;
  }
}

void FixSymmetry::restore_std_cell() {

  // Update cell matrix
  domain->boxlo[0] = domain->boxlo[1] = domain->boxlo[2] = 0.0;
  domain->boxhi[0] = MathExtra::len3(std_cell[0]);
  domain->boxhi[1] = MathExtra::len3(std_cell[1]);
  domain->boxhi[1] = MathExtra::len3(std_cell[2]);
  domain->xy = std_cell[1][0];
  domain->xz = std_cell[2][0];
  domain->yz = std_cell[2][1];

  // Update related parameters
  domain->set_global_box();
  domain->set_local_box();
}


void FixSymmetry::adjust_cell() {
  store_std_cell();
  symmetrize_rank2(std_cell);
  restore_std_cell();
}

void FixSymmetry::adjust_positions() {
  double **x = atom->x;
  int nlocal = atom->nlocal;

  std::vector<double[3]> pos(nlocal);
  for (int i = 0; i < nlocal; i++) {
    pos[i][0] = x[i][0];
    pos[i][1] = x[i][1];
    pos[i][2] = x[i][2];
  }

  // Symmetrize Positions
  symmetrize_rank1(pos);

  // Update 
  for (int i = 0; i < nlocal; i++) {
    x[i][0] = pos[i][0];
    x[i][1] = pos[i][1];
    x[i][2] = pos[i][2];
  }
}

void FixSymmetry::adjust_forces() {
  double **f = atom->f;
  int nlocal = atom->nlocal;

  // Convert forces to vector format
  std::vector<double[3]> forces(nlocal);
  for (int i = 0; i < nlocal; i++) {
    forces[i][0] = f[i][0];
    forces[i][1] = f[i][1];
    forces[i][2] = f[i][2];
  }

  // Symmetrize forces
  symmetrize_rank1(forces);

  // Update forces
  for (int i = 0; i < nlocal; i++) {
    f[i][0] = forces[i][0];
    f[i][1] = forces[i][1];
    f[i][2] = forces[i][2];
  }
}

void FixSymmetry::adjust_stress() {
  // Get global stress tensor
  double *stress = virial;
  double stress_tensor[3][3];
  stress_tensor[0][0] = stress[0];
  stress_tensor[1][1] = stress[1];
  stress_tensor[2][2] = stress[2];
  stress_tensor[1][2] = stress_tensor[2][1] = stress[3];
  stress_tensor[0][2] = stress_tensor[2][0] = stress[4];
  stress_tensor[0][1] = stress_tensor[1][0] = stress[5];

  // Symmetrize stress tensor
  symmetrize_rank2(stress_tensor);

  // Update stress tensor in LAMMPS
  stress[0] = stress_tensor[0][0];
  stress[1] = stress_tensor[1][1];
  stress[2] = stress_tensor[2][2];
  stress[3] = stress_tensor[1][2];
  stress[4] = stress_tensor[0][2];
  stress[5] = stress_tensor[0][1];
}


void FixSymmetry::x2lambda(const double pos[3], double lambda[3]) {
  double inv_std_cell[3][3];
  MathExtra::invert3(std_cell, inv_std_cell);
  MathExtra::matvec(inv_std_cell, pos, lambda);
}

void FixSymmetry::print_symmetry() {
  std::cout << "Precision: " << symprec
            << ", Space group number: " << dataset->spacegroup_number
            << ", Sym N operations: " << dataset->spacegroup_number
            << ", International symbol: " << dataset->n_operations
            << ", Hall symbol: " << dataset->hall_symbol << std::endl;
}

void FixSymmetry::refine_symmetry() {
  check_and_symmetrize_cell();
  check_and_symmetrize_positions();
  return check_symmetry(false);
}

void FixSymmetry::check_and_symmetrize_cell(){
  check_symmetry(false);
  symmetrize_cell();
}

void FixSymmetry::symmetrize_cell() {

  // Get standard cell matrix
  double trans_std_cell[3][3];
  MathExtra::transpose_times3(dataset->transformation_matrix,  dataset->std_lattice, trans_std_cell);
  MathExtra::times3(trans_std_cell,  dataset->std_rotation_matrix, std_cell);

  //dump cell matrix
  printf("Symmetrized Cell Matrix = ");
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      printf("%f,", std_cell[i][j]);
    }
  }
  printf("\n");

  restore_std_cell();
}

void FixSymmetry::check_and_symmetrize_positions(){
  check_symmetry(true);
  symmetrize_positions();
}

void FixSymmetry::symmetrize_positions() {
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double rot_std_cell[3][3];
  double rot_prim_cell[3][3];
  double inv_rot_prim_cell[3][3];

  double **rot_std_pos;

  MathExtra::times3(std_cell,  dataset->std_rotation_matrix, rot_std_cell);

  int n_std = dataset->n_std_atoms;
  memory->create(rot_std_pos, n_std, 3, "fix_symmetry:rot_std_pos");
  for (int i=0; i < n_std; i++) {
    MathExtra::matvec(rot_std_cell, dataset->std_positions[i], rot_std_pos[i]);
  }
  int id0 = get_index(mapping_to_primitive, 0);
  int id1 = get_index(std_mapping_to_primitive, 0);

  double dp0[3];
  dp0[0] = x[id0][0] - rot_std_pos[id1][0];
  dp0[1] = x[id0][1] - rot_std_pos[id1][1];
  dp0[2] = x[id0][2] - rot_std_pos[id1][2];

  MathExtra::times3(prim_cell,  dataset->std_rotation_matrix, rot_prim_cell);
  MathExtra::invert3(rot_prim_cell, inv_rot_prim_cell);

  for (int i=0; i < nlocal; i++) {
    if (debug) {
      printf("mapping_to_primitive[%d] = %d\n", i, mapping_to_primitive[i]);
      printf("std_mapping_to_primitive[%d] = %d\n", i, std_mapping_to_primitive[i]);
    }

    int std_i = get_index(std_mapping_to_primitive, mapping_to_primitive[i]);

    double dp[3];
    double dp_s[3];
    double dp_sr[3];
    double dp_sr_mul[3];

    dp[0] = rot_std_pos[std_i][0] + dp0[0] - x[i][0];
    dp[1] = rot_std_pos[std_i][1] + dp0[1] - x[i][1];
    dp[2] = rot_std_pos[std_i][2] + dp0[2] - x[i][2];

    MathExtra::matvec(inv_rot_prim_cell, dp, dp_s);
    dp_sr[0] = std::round(dp_s[0]);
    dp_sr[1] = std::round(dp_s[1]);
    dp_sr[2] = std::round(dp_s[2]);
    MathExtra::matvec(rot_prim_cell, dp_sr, dp_sr_mul);


    x[i][0] = rot_std_pos[std_i][0] + dp0[0] - dp_sr_mul[0]; 
    x[i][1] = rot_std_pos[std_i][1] + dp0[1] - dp_sr_mul[1]; 
    x[i][2] = rot_std_pos[std_i][2] + dp0[2] - dp_sr_mul[2]; 
  }

  memory->destroy(rot_std_pos);
}

void FixSymmetry::initial_prep() {
  // Prepare data structures for spglib
  int num_atoms = atom->natoms;  // Total number of atoms

  // Collect positions and types from all processes
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;

  // Total number of atoms
  int natoms = atom->natoms;

  memory->create(all_positions, natoms, 3, "fix_symmetry:all_positions");
  memory->create(all_types, natoms, "fix_symmetry:all_types");

  // Store local positions and types
  for (int i = 0; i < nlocal; i++) {
    int id = tag[i] - 1;  // Assuming atom IDs start from 1
    all_positions[id][0] = x[i][0];
    all_positions[id][1] = x[i][1];
    all_positions[id][2] = x[i][2];
    all_types[id] = type[i];
  }

  // Gather data from all processes
  MPI_Allreduce(MPI_IN_PLACE, all_positions, natoms * 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, all_types, natoms, MPI_INT, MPI_SUM, world);
}

void FixSymmetry::check_symmetry(bool do_find_prim) {

  int natoms = atom->natoms;
  double (*my_all_positions)[3] = new double[natoms][3];

  // copy all_positions to my_all_positions
  for (int i = 0; i < natoms; i++) {
    my_all_positions[i][0] = all_positions[i][0];
    my_all_positions[i][1] = all_positions[i][1];
    my_all_positions[i][2] = all_positions[i][2];
  }

  store_std_cell();
  dataset = spg_get_dataset(std_cell, my_all_positions, all_types, natoms, symprec);

  if (dataset == nullptr) {
    std::string spg_error_msg = spg_get_error_message(spg_get_error_code());
    error->all(FLERR, spg_error_msg.c_str());
  }
  //add debug code to display all_positions and all_types
  if (debug) {
    for (int i = 0; i < natoms; i++) {
      printf("pos[%d] = %f %f %f, ", i, all_positions[i][0], all_positions[i][1], all_positions[i][2]);
      printf("type[%d] = %d\n", i, all_types[i]);
    }
  }
  print_symmetry();

  if (do_find_prim) {
    int find_prim = spg_find_primitive(std_cell, my_all_positions, all_types, natoms, symprec);
    printf("find_prim = %d, ", find_prim);
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        if (find_prim > 0) {
          prim_cell[j][k] = dataset->primitive_lattice[j][k];
        } else {
          prim_cell[j][k] =std_cell[j][k];
        }
      }
    }
    printf("prim_cell = ");
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        printf("%f,", prim_cell[j][k]);
      }
    }
    printf("\n");

    // Atom mapping
    mapping_to_primitive.clear();
    std_mapping_to_primitive.clear();
    for (int i = 0; i < atom->nlocal; i++) {
      int id = atom->tag[i] - 1;
      mapping_to_primitive.push_back(dataset->mapping_to_primitive[id]);
      std_mapping_to_primitive.push_back(dataset->std_mapping_to_primitive[id]);
    }
  }
}

void FixSymmetry::prep_symmetry() {
  // Get symmetry operations
  int nsym = dataset->n_operations;
  rotation_matrices.clear();
  translation_vectors.clear();

  std::vector<std::vector<double>> scaled_pos;
  //eval scaled pos
  for (int i = 0; i < atom->nlocal; i++) {
    double pos[3];
    pos[0] = atom->x[i][0];
    pos[1] = atom->x[i][1];
    pos[2] = atom->x[i][2];
    x2lambda(pos, pos);
    scaled_pos.push_back(std::vector<double>(pos, pos + 3));
  }

  for (int i = 0; i < nsym; i++) {
    std::array<std::array<double, 3>, 3> R;
    std::array<double, 3> t;
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        R[j][k] = dataset->rotations[i][j][k];
      }
      t[j] = dataset->translations[i][j];
    }
    rotation_matrices.push_back(R);
    translation_vectors.push_back(t);
  }
}

void FixSymmetry::symmetrize_rank1(std::vector<double[3]> &vec) {
  int nlocal = atom->nlocal;
  int nsym = rotation_matrices.size();

  // Initialize symmetrized vectors
  std::vector<double[3]> sym_vec(nlocal);
  for (int i = 0; i < nlocal; i++) {
    sym_vec[i][0] = 0.0;
    sym_vec[i][1] = 0.0;
    sym_vec[i][2] = 0.0;
  }

  for (int k = 0; k < nsym; k++) {
    double R[3][3];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        R[i][j] = rotation_matrices[k][i][j];
      }
    }

    for (int i = 0; i < nlocal; i++) {
      int j = mapping_to_primitive[i];
      if (j >= 0 && j < nlocal) {
        double tmp[3];
        MathExtra::matvec(R, vec[j], tmp);
        sym_vec[i][0] += tmp[0];
        sym_vec[i][1] += tmp[1];
        sym_vec[i][2] += tmp[2];
      }
    }
  }

  // Average
  for (int i = 0; i < nlocal; i++) {
    sym_vec[i][0] /= nsym;
    sym_vec[i][1] /= nsym;
    sym_vec[i][2] /= nsym;
    vec[i][0] = sym_vec[i][0];
    vec[i][1] = sym_vec[i][1];
    vec[i][2] = sym_vec[i][2];
  }
}

void FixSymmetry::symmetrize_rank2(double tensor[3][3]) {
  int nsym = rotation_matrices.size();

  double sym_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    sym_tensor[i][0] = 0.0;
    sym_tensor[i][1] = 0.0;
    sym_tensor[i][2] = 0.0;
  }

  for (int k = 0; k < nsym; k++) {
    double R[3][3];
    double Rt[3][3];
    double tmp1[3][3];
    double tmp2[3][3];

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
      R[i][j] = rotation_matrices[k][i][j];
      }
    }

    MathExtra::transpose3(R, Rt);
    MathExtra::times3(tensor, Rt, tmp1);
    MathExtra::times3(R, tmp1, tmp2);

    for (int i = 0; i < 3; i++) {
      sym_tensor[i][0] += tmp2[i][0];
      sym_tensor[i][1] += tmp2[i][1];
      sym_tensor[i][2] += tmp2[i][2];
    }
  }

  // Average
  for (int i = 0; i < 3; i++) {
    tensor[i][0] = sym_tensor[i][0] / nsym;
    tensor[i][1] = sym_tensor[i][1] / nsym;
    tensor[i][2] = sym_tensor[i][2] / nsym;
  }
}


/*
  if(spacegroup_number != dataset->spacegroup_number) {
    std::string str = "Space group number mismatch. Expected " + std::to_string(spacegroup_number) + " but got " + std::to_string(dataset->spacegroup_number);
    //error->all(FLERR, str.c_str());
    printf("%s", str.c_str());
  }
  */


