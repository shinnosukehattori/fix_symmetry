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
    double cell[3][3];
    double inv_cell[3][3];

    get_cell(cell);
    MathExtra::invert3(cell, inv_cell);
    //MathExtra::transpose3(inv_cell, inv_cell);
    adjust_cell(cell, inv_cell);

    get_cell(cell);
    MathExtra::invert3(cell, inv_cell);
    //adjust_positions(cell, inv_cell);
    if (symforce) {
      adjust_forces(cell, inv_cell);
    }
    if (symstress) {
      adjust_stress(cell, inv_cell);
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

void FixSymmetry::get_cell(double cell[3][3]) {
  cell[0][1] = cell[0][2] = cell[1][2] = 0.0;
  cell[0][0] = domain->h[0];
  cell[1][1] = domain->h[1];
  cell[2][2] = domain->h[2];
  if (domain->triclinic) {
    cell[1][0] = cell[2][0] = cell[2][1] = 0.0;
  } else {
    cell[1][0] = domain->xy;
    cell[2][0] = domain->xz;
    cell[2][1] = domain->yz;
  }
}

void FixSymmetry::set_cell(double cell[3][3]) {

  // Update cell matrix
  domain->boxlo[0] = domain->boxlo[1] = domain->boxlo[2] = 0.0;
  domain->boxhi[0] = MathExtra::len3(cell[0]);
  domain->boxhi[1] = MathExtra::len3(cell[1]);
  domain->boxhi[1] = MathExtra::len3(cell[2]);
  domain->xy = cell[1][0];
  domain->xz = cell[2][0];
  domain->yz = cell[2][1];

  // Update related parameters
  domain->set_global_box();
  domain->set_local_box();
}


void FixSymmetry::adjust_cell(double cell[3][3], double inv_cell[3][3]) {

  double delta_defrom_grad[3][3];

  MathExtra::transpose_times3(inv_cell, sym_cell, delta_defrom_grad);
  delta_defrom_grad[0][0] -= 1.0;
  delta_defrom_grad[1][1] -= 1.0;
  delta_defrom_grad[2][2] -= 1.0;

  double max_delta_grad = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double delta = fabs(delta_defrom_grad[i][j]);
      if (delta > max_delta_grad) {
        max_delta_grad = delta;
      }
    }
  }
  if (max_delta_grad > 0.25){
    error->all(FLERR, "Too large deformation gradient in FixSymmetry::adjust_cell");
  } else if (max_delta_grad > 1.0) {
    printf("Warning! max_delta_grad = %f\n", max_delta_grad);
  }

  symmetrize_rank2(delta_defrom_grad, cell, inv_cell);

  delta_defrom_grad[0][0] += 1.0;
  delta_defrom_grad[1][1] += 1.0;
  delta_defrom_grad[2][2] += 1.0;

  //update
  MathExtra::transpose_times3(delta_defrom_grad, cell, sym_cell);
  set_cell(sym_cell);
}

void FixSymmetry::adjust_positions(double cell[3][3], double inv_cell[3][3]) {
  double **x = atom->x;
  int nlocal = atom->nlocal;

  std::vector<double[3]> pos(nlocal);
  for (int i = 0; i < nlocal; i++) {
    pos[i][0] = x[i][0];
    pos[i][1] = x[i][1];
    pos[i][2] = x[i][2];
  }

  // Symmetrize Positions
  symmetrize_rank1(pos, cell, inv_cell);

  // Update 
  for (int i = 0; i < nlocal; i++) {
    x[i][0] = pos[i][0];
    x[i][1] = pos[i][1];
    x[i][2] = pos[i][2];
  }
}

void FixSymmetry::adjust_forces(double cell[3][3], double inv_cell[3][3]) {
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
  symmetrize_rank1(forces, cell, inv_cell);

  // Update forces
  for (int i = 0; i < nlocal; i++) {
    f[i][0] = forces[i][0];
    f[i][1] = forces[i][1];
    f[i][2] = forces[i][2];
  }
}

void FixSymmetry::adjust_stress(double cell[3][3], double inv_cell[3][3]) {
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
  symmetrize_rank2(stress_tensor, cell, inv_cell);

  // Update stress tensor in LAMMPS
  stress[0] = stress_tensor[0][0];
  stress[1] = stress_tensor[1][1];
  stress[2] = stress_tensor[2][2];
  stress[3] = stress_tensor[1][2];
  stress[4] = stress_tensor[0][2];
  stress[5] = stress_tensor[0][1];
}


void FixSymmetry::x2lambda(const double pos[3], double lambda[3]) {
  double inv_cell[3][3];
  double cell[3][3];
  get_cell(cell);
  MathExtra::invert3(cell, inv_cell);
  MathExtra::matvec(inv_cell, pos, lambda);
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
  MathExtra::times3(trans_std_cell,  dataset->std_rotation_matrix, sym_cell);

  //dump cell matrix
  printf("Symmetrized Cell Matrix = ");
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      printf("%f,", sym_cell[i][j]);
    }
  }
  printf("\n");

  set_cell(sym_cell);
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

  if (debug) {
    for (int i = 0; i < nlocal; i++) {
      printf("pos0[%d] = %f %f %f, ", i, x[i][0], x[i][1], x[i][2]);
      printf("type[%d] = %d\n", i, all_types[i]);
    }
  }

  MathExtra::times3(dataset->std_lattice,  dataset->std_rotation_matrix, rot_std_cell);

  int n_std = dataset->n_std_atoms;
  memory->create(rot_std_pos, n_std, 3, "fix_symmetry:rot_std_pos");
  for (int i=0; i < n_std; i++) {
    MathExtra::transpose_matvec(rot_std_cell, dataset->std_positions[i], rot_std_pos[i]);
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
      int j = mapping_to_primitive[i];
      printf("mapping_to_primitive[%d] = %d\n", i, j);
      printf("std_mapping_to_primitive[%d] = %d\n", j, std_mapping_to_primitive[j]);
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
  if (debug) {
    for (int i = 0; i < nlocal; i++) {
      printf("pos[%d] = %f %f %f, ", i, x[i][0], x[i][1], x[i][2]);
      printf("type[%d] = %d\n", i, all_types[i]);
    }
  }
  print_symmetry();

  memory->destroy(rot_std_pos);
}

void FixSymmetry::initial_prep() {
  // Prepare data structures for spglib
  int natoms = atom->natoms;  // Total number of atoms

  // Collect positions and types from all processes
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;

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

  //add debug code to display all_positions and all_types
  if (debug) {
    for (int i = 0; i < natoms; i++) {
      printf("allpos[%d] = %f %f %f, ", i, all_positions[i][0], all_positions[i][1], all_positions[i][2]);
      printf("type[%d] = %d\n", i, all_types[i]);
    }
  }
}

void FixSymmetry::check_symmetry(bool do_find_prim) {

  initial_prep();
  int natoms = atom->natoms;
  double (*scaled_all_positions)[3] = new double[natoms][3];

  // copy all_positions to my_all_positions
  for (int i = 0; i < natoms; i++) {
    x2lambda(all_positions[i],scaled_all_positions[i]);
  }

  double cell[3][3];
  get_cell(cell);
  dataset = spg_get_dataset(cell, scaled_all_positions, all_types, natoms, symprec);

  if (dataset == nullptr) {
    std::string spg_error_msg = spg_get_error_message(spg_get_error_code());
    error->all(FLERR, spg_error_msg.c_str());
  }
  print_symmetry();

  if (do_find_prim) {
    int find_prim = spg_find_primitive(cell, scaled_all_positions, all_types, natoms, symprec);
    printf("find_prim = %d, ", find_prim);
    printf("prim_cell = ");
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        prim_cell[k][j] = cell[j][k]; //transposed
      }
    }
    printf("\n");

    // Atom mapping for primmitive
    mapping_to_primitive.clear();
    std_mapping_to_primitive.clear();
    for (int i = 0; i < atom->nlocal; i++) {
      int id = atom->tag[i] - 1;
      mapping_to_primitive.push_back(dataset->mapping_to_primitive[id]);
      std_mapping_to_primitive.push_back(dataset->std_mapping_to_primitive[id]);
    }
  }
}


std::vector<std::vector<double>> dp_mat_sub(int size, double  **mat, const double vec[3]) {
    std::vector<std::vector<double>> result(size, std::vector<double>(3));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = mat[i][j] - vec[j];
            result[i][j] -= std::round(result[i][j]); // Ensure periodic boundary conditions
        }
    }
    return result;
}

void FixSymmetry::prep_symmetry() {
  // Get symmetry operations
  int natoms = atom->nlocal;
  int nsym = dataset->n_operations;
  //memory->create(all_positions, natoms, 3, "fix_symmetry:all_positions");
  memory->create(rotation_matrices, nsym, 3, 3, "fix_symmetry:rotation_matrices");
  memory->create(translation_vectors, nsym, 3, "fix_symmetry:translation_vectors");
  symm_map.clear();

  double **scaled_pos;
  memory->create(scaled_pos, atom->nlocal, 3, "fix_symmetry:scaled_pos");

  //eval scaled pos
  for (int i = 0; i < atom->nlocal; i++) {
    double pos[3];
    pos[0] = atom->x[i][0];
    pos[1] = atom->x[i][1];
    pos[2] = atom->x[i][2];
    x2lambda(pos, pos);
    scaled_pos[i][0] = pos[0];
    scaled_pos[i][1] = pos[1];
    scaled_pos[i][2] = pos[2];
  }

  printf("nsym = %d\n", nsym);
  for (int i = 0; i < nsym; i++) {
    double R[3][3];
    double t[3];
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        rotation_matrices[i][j][k] = dataset->rotations[i][j][k];
        R[j][k] = dataset->rotations[i][j][k];
      }
      translation_vectors[i][j] = dataset->translations[i][j];
      t[j] = dataset->translations[i][j];
    }
    auto this_op_map = std::vector<int>(natoms, -1);

    for (int j = 0; j < natoms; j++) {
      double new_p[3];
      MathExtra::matvec(R, scaled_pos[j], new_p);
      new_p[0] += t[0];
      new_p[1] += t[1];
      new_p[2] += t[2];

      auto dp = dp_mat_sub(natoms, scaled_pos, new_p);
      int min_index = 0;
      double min_delta = 9999.0;
      for (int k = 0; k < natoms; k++) {
        double delta = sqrt(dp[k][0]*dp[k][0] + dp[k][1]*dp[k][1] + dp[k][2]*dp[k][2]);
        if (delta < min_delta) {
          min_index = k;
          min_delta = delta;
        }
      }
      this_op_map[j] = min_index;
    }
    symm_map.push_back(this_op_map);
  }
}

void FixSymmetry::symmetrize_rank1(std::vector<double[3]> &vec, const double cell[3][3], const double inv_cell[3][3]) {
  int nlocal = atom->nlocal;
  int nsym = dataset->n_operations;

  std::vector<double[3]> scaled_vec(nlocal);

  for (int i = 0; i < nlocal; i++) {
    MathExtra::matvec(inv_cell, vec[i], scaled_vec[i]);
  }


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
      int j = symm_map[k][i];
      if (j >= 0 && j < nlocal) {
        double tmp[3];
        MathExtra::matvec(R, scaled_vec[i], tmp);
        sym_vec[j][0] += tmp[0];
        sym_vec[j][1] += tmp[1];
        sym_vec[j][2] += tmp[2];
      }
    }
  }

  // Average
  for (int i = 0; i < nlocal; i++) {
    sym_vec[i][0] /= nsym;
    sym_vec[i][1] /= nsym;
    sym_vec[i][2] /= nsym;
  }

  for (int i = 0; i < nlocal; i++) {
    MathExtra::matvec(cell, sym_vec[i], vec[i]);
  }
}

void FixSymmetry::symmetrize_rank2(double tensor[3][3], const double cell[3][3], const double inv_cell[3][3]) {
  int nlocal = atom->nlocal;
  int nsym = dataset->n_operations;

  std::vector<double[3]> scaled_vec(nlocal);
  double t_cell[3][3];
  double t_inv_cell[3][3];
  MathExtra::transpose3(cell, t_cell);
  MathExtra::transpose3(inv_cell, t_inv_cell);

  double tmp0[3][3];
  double tmp1[3][3];
  double tmp4[3][3];
  double tmp5[3][3];

  MathExtra::times3(cell, tensor, tmp0);
  MathExtra::times3(tmp0, t_cell, tmp1);

  double sym_tensor[3][3];
  for (int i = 0; i < 3; i++) {
    sym_tensor[i][0] = 0.0;
    sym_tensor[i][1] = 0.0;
    sym_tensor[i][2] = 0.0;
  }

  for (int k = 0; k < nsym; k++) {
    double R[3][3];
    double Rt[3][3];
    double tmp2[3][3];
    double tmp3[3][3];

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
      R[i][j] = rotation_matrices[k][i][j];
      }
    }

    MathExtra::transpose3(R, Rt);
    MathExtra::times3(tmp1, Rt, tmp2);
    MathExtra::times3(R, tmp2, tmp3);

    for (int i = 0; i < 3; i++) {
      sym_tensor[i][0] += tmp3[i][0];
      sym_tensor[i][1] += tmp3[i][1];
      sym_tensor[i][2] += tmp3[i][2];
    }
  }

  // Average
  for (int i = 0; i < 3; i++) {
    tmp4[i][0] = sym_tensor[i][0] / nsym;
    tmp4[i][1] = sym_tensor[i][1] / nsym;
    tmp4[i][2] = sym_tensor[i][2] / nsym;
  }

  MathExtra::times3(inv_cell, tmp4, tmp5);
  MathExtra::times3(tmp5, t_inv_cell, tensor);
}
