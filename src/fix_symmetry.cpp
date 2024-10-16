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
#include <sstream>

using namespace LAMMPS_NS;


FixSymmetry::FixSymmetry(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {

  symprec = 1e-4;  // Default tolerance
  symcell = false;
  symposs = true;
  debug = false;

  int proc_x = comm->procgrid[0];
  int proc_y = comm->procgrid[1];
  int proc_z = comm->procgrid[2];

  if (proc_x >= 2 || proc_y >= 2 || proc_z >= 2) {
    error->all(FLERR, "Processor grid size in any direction must be less than 2 for fix symmetry.");
  }


  if (narg < 3 || narg > 7)
    error->all(FLERR, "Illegal fix symmetry command. optional symprec, symcell, symforce, symstress and debug-mode.");

  // Get the tolerance (if not specified, set default value)
  if (narg >= 4) symprec = utils::numeric(FLERR, arg[3], false, lmp);
  if (narg >= 5) symcell = utils::logical(FLERR, arg[4], false, lmp);
  if (narg >= 6) symposs= utils::logical(FLERR, arg[5], false, lmp);
  if (narg >= 7) debug = utils::logical(FLERR, arg[6], false, lmp);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      prim_cell[i][j] = 0.0;
    }
  }

  spacegroup_number = -1;

  // Initialize previous positions array
  prev_positions = nullptr;
  all_positions = nullptr;
  all_types = nullptr;
  rotation_matrices = nullptr;
  translation_vectors = nullptr;
  dataset = nullptr;
}

FixSymmetry::~FixSymmetry() {
  // Free memory
  memory->destroy(prev_positions);
  memory->destroy(all_positions);
  memory->destroy(all_types);
  memory->destroy(rotation_matrices);
  memory->destroy(translation_vectors);
  if (dataset) {
    spg_free_dataset(dataset);
  }
}

int FixSymmetry::setmask() {
  int mask = 0;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::MIN_PRE_FORCE;
  mask |= FixConst::END_OF_STEP;
  mask |= FixConst::MIN_POST_FORCE;
  mask |= FixConst::POST_RUN;
  return mask;
}

void FixSymmetry::init() {
  refine_symmetry();
  prep_symmetry();
}


void FixSymmetry::setup_pre_force(int vflag) {
  double cell[3][3];
  double inv_cell[3][3];

  get_cell(cell, inv_cell);
  adjust_forces(cell, inv_cell);
  adjust_stress(cell, inv_cell);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sym_cell[i][j] = cell[i][j];
    }
  }
  if(symposs) save_prev_position();
}

void FixSymmetry::end_of_step() {
  double cell[3][3];
  double inv_cell[3][3];


  if(symcell) adjust_cell();

  get_cell(cell, inv_cell);
  if(symposs) adjust_positions(cell, inv_cell);
  adjust_forces(cell, inv_cell);
  adjust_stress(cell, inv_cell);

  if(symposs) save_prev_position();
}

void FixSymmetry::min_post_force(int vflag) {
  end_of_step();
}

void FixSymmetry::post_run() {
  store_all_coordinates();
  check_symmetry(true, false);
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
    cell[2][1] = domain->h[3]; //yz
    cell[2][0] = domain->h[4]; //xz
    cell[1][0] = domain->h[5]; //xy
  } else {
    cell[1][0] = cell[2][0] = cell[2][1] = 0.0;
  }
}
void FixSymmetry::get_cell(double cell[3][3], double inv_cell[3][3]) {
  get_cell(cell);
  inv_cell[0][1] = inv_cell[0][2] = inv_cell[1][2] = 0.0;
  inv_cell[0][0] = domain->h_inv[0];
  inv_cell[1][1] = domain->h_inv[1];
  inv_cell[2][2] = domain->h_inv[2];
  if (domain->triclinic) {
    inv_cell[2][1] = domain->h_inv[3]; //yz
    inv_cell[2][0] = domain->h_inv[4]; //xz
    inv_cell[1][0] = domain->h_inv[5]; //xy
  } else {
    inv_cell[1][0] = inv_cell[2][0] = inv_cell[2][1] = 0.0;
  }
}

void FixSymmetry::set_cell(double cell[3][3]) {

  // Update cell matrix
  domain->boxhi[0] = cell[0][0];
  domain->boxhi[1] = cell[1][1];
  domain->boxhi[2] = cell[2][2];
  if (domain->triclinic) {
    domain->yz = cell[2][1];
    domain->xz = cell[2][0];
    domain->xy = cell[1][0];
  }

  // Update related parameters
  domain->set_global_box();
  domain->set_local_box();
}

void FixSymmetry::save_prev_position() {

  if (prev_positions != nullptr || 3*atom->nlocal != sizeof(prev_positions)) {
    memory->destroy(prev_positions);
  }
  if (prev_positions == nullptr) {
    memory->create(prev_positions, atom->nlocal, 3, "fix_symmetry:prev_positions");
  }
  for (int i = 0; i < atom->nlocal; i++) {
    domain->unmap(atom->x[i], atom->image[i], prev_positions[i]);
  }
}

void FixSymmetry::adjust_cell() {

  double cell[3][3];
  double inv_sym_cell[3][3];
  double tmp[3][3];

  get_cell(cell);
  MathExtra::invert3(sym_cell, inv_sym_cell);
  double delta_deform_grad0[3][3];
  double delta_deform_grad[3][3];
  MathExtra::times3(inv_sym_cell, cell, delta_deform_grad0);
  MathExtra::transpose3(delta_deform_grad0, delta_deform_grad);

  delta_deform_grad[0][0] -= 1.0;
  delta_deform_grad[1][1] -= 1.0;
  delta_deform_grad[2][2] -= 1.0;

  double max_delta_grad = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double delta = fabs(delta_deform_grad[i][j]);
      if (delta > max_delta_grad) {
        max_delta_grad = delta;
      }
    }
  }
  if (max_delta_grad > 0.25){
    error->all(FLERR, "Too large deformation gradient in FixSymmetry::adjust_cell");
  } else if (max_delta_grad > 1.0) {
    std::string message = fmt::format("Warning! max_delta_grad = {:f}\n", max_delta_grad);
    utils::logmesg(lmp, message);
  }

  symmetrize_rank2(delta_deform_grad, sym_cell, inv_sym_cell);

  delta_deform_grad[0][0] += 1.0;
  delta_deform_grad[1][1] += 1.0;
  delta_deform_grad[2][2] += 1.0;

  //update
  MathExtra::transpose3(delta_deform_grad, delta_deform_grad0);
  MathExtra::times3(sym_cell, delta_deform_grad0, tmp);

  if (debug) utils::logmesg(lmp, "Symmetrized Cell Matrix (adj) = ");
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (debug){
        std::string message = fmt::format("{:5.3f}->{:5.3f}, ", cell[i][j], tmp[i][j]);
        utils::logmesg(lmp, message);
      }
      sym_cell[i][j] = tmp[i][j];
    }
  }
  if (debug) utils::logmesg(lmp, "\n");
  domain->x2lamda(atom->nlocal);
  set_cell(tmp);
  domain->lamda2x(atom->nlocal);
}

void FixSymmetry::adjust_positions(double cell[3][3], double inv_cell[3][3]) {
  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  std::vector<double[3]> step(nlocal);
  for (int i = 0; i < nlocal; i++) {
    double pos[3];
    domain->unmap(x[i], image[i], pos);
    step[i][0] = pos[0] - prev_positions[i][0];
    step[i][1] = pos[1] - prev_positions[i][1];
    step[i][2] = pos[2] - prev_positions[i][2];
  }

  // Symmetrize Positions
  symmetrize_rank1(step, cell, inv_cell);

  // Update 
  for (int i = 0; i < nlocal; i++) {
    x[i][0] = prev_positions[i][0] + step[i][0];
    x[i][1] = prev_positions[i][1] + step[i][1];
    x[i][2] = prev_positions[i][2] + step[i][2];
    domain->unmap_inv(x[i], image[i]);
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
  stress_tensor[1][2] = stress_tensor[2][1] = stress[3]*0.5;
  stress_tensor[0][2] = stress_tensor[2][0] = stress[4]*0.5;
  stress_tensor[0][1] = stress_tensor[1][0] = stress[5]*0.5;

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

void FixSymmetry::print_symmetry(int prev) {
  std::ostringstream message;
  message  << "[Sym]"
           << "SG: " << dataset->spacegroup_number;
  if (prev >= 0) {
     message << " (" << prev << ")";
  }
  message << " Prec.: " << symprec
          << " SymNops: " << dataset->n_operations
          << " Int.Symbol: " << dataset->international_symbol
          << " HallSymbol: " << dataset->hall_symbol << std::endl;
  if(comm->me == 0) utils::logmesg(lmp, message.str());
}

void FixSymmetry::refine_symmetry() {
  //check_and_symmetrize_cell
  store_all_coordinates();
  check_symmetry(true, false);
  symmetrize_cell();
  check_symmetry(false, true); //use symmetrized cell and find primitive
  symmetrize_positions();

  store_all_coordinates();
  check_symmetry(false, false);
  return;
}

void FixSymmetry::symmetrize_cell() {

  // Get standard cell matrix
  double trans_std_cell[3][3];
  double std_cell[3][3];
  double cell[3][3];
  MathExtra::transpose3(dataset->std_lattice, std_cell);
  MathExtra::transpose_times3(dataset->transformation_matrix,  std_cell, trans_std_cell);
  MathExtra::times3(trans_std_cell,  dataset->std_rotation_matrix, cell);

  //dump cell matrix
  utils::logmesg(lmp, "Symmetrized Cell Matrix = ");
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      std::string message = fmt::format("{:f}", cell[i][j]);
      utils::logmesg(lmp, message);
      message = fmt::format("({:f}),", std_cell[i][j]);
      utils::logmesg(lmp, message);
    }
  }
  utils::logmesg(lmp, "\n");

  domain->x2lamda(atom->nlocal);
  set_cell(cell);
  domain->lamda2x(atom->nlocal);
}

void FixSymmetry::symmetrize_positions() {
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double std_cell[3][3];
  MathExtra::transpose3(dataset->std_lattice, std_cell);
  double rot_std_cell[3][3];
  double rot_prim_cell[3][3];
  double inv_rot_prim_cell[3][3];


  if (debug) {
    for (int i = 0; i < nlocal; i++) {
      printf("pos0[%d] = %f %f %f, ", i, x[i][0], x[i][1], x[i][2]);
      printf("type[%d] = %d\n", i, all_types[i]);
    }
  }

  MathExtra::times3(std_cell,  dataset->std_rotation_matrix, rot_std_cell);

  int n_std = dataset->n_std_atoms;
  double **rot_std_pos;

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
      printf("std_mapping_to_primitive[%d] = %d\n", j, get_index(std_mapping_to_primitive, j));
    }

    int std_i = get_index(std_mapping_to_primitive, mapping_to_primitive[i]);

    double dp[3];
    double dp_s[3];
    double dp_sr[3];
    double dp_sr_mul[3];

    dp[0] = rot_std_pos[std_i][0] + dp0[0] - x[i][0];
    dp[1] = rot_std_pos[std_i][1] + dp0[1] - x[i][1];
    dp[2] = rot_std_pos[std_i][2] + dp0[2] - x[i][2];

    MathExtra::transpose_matvec(inv_rot_prim_cell, dp, dp_s);
    dp_sr[0] = std::round(dp_s[0]);
    dp_sr[1] = std::round(dp_s[1]);
    dp_sr[2] = std::round(dp_s[2]);
    MathExtra::transpose_matvec(rot_prim_cell, dp_sr, dp_sr_mul);


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

  memory->destroy(rot_std_pos);
}

void FixSymmetry::store_all_coordinates() {
  // Prepare data structures for spglib
  int natoms = atom->natoms;  // Total number of atoms

  // Collect positions and types from all processes
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;

  if (debug) printf("Store coordinates: :natoms = %d, nlocal = %d\n", natoms, nlocal);
  if (all_positions) delete [] all_positions;
  if (all_types) memory->destroy(all_types);

  all_positions  = new double[natoms*3];
  memory->create(all_types, natoms, "fix_symmetry:all_types");

  for (int i = 0; i < natoms; i++) {
    all_positions[3*i] = 0.0;
    all_positions[3*i+1] = 0.0;
    all_positions[3*i+2] = 0.0;
    all_types[i] = 0;
  }

  // Store local positions and types
  domain->x2lamda(nlocal);
  for (int i = 0; i < nlocal; i++) {
    int id = tag[i] - 1;  // Assuming atom IDs start from 1
    all_positions[3*id+0] = x[i][0];
    all_positions[3*id+1] = x[i][1];
    all_positions[3*id+2] = x[i][2];
    all_types[id] = type[i];
  }
  domain->lamda2x(nlocal);

  //!!!!!!! Dose not support space decomposition system!!!!!!!!!!!
  //MPI_Allreduce(MPI_IN_PLACE, all_positions, natoms*3, MPI_DOUBLE, MPI_MAX, world);
  //MPI_Allreduce(MPI_IN_PLACE, all_types, natoms, MPI_INT, MPI_SUM, world);
}

void FixSymmetry::check_symmetry(bool do_print, bool do_find_prim) {

  int natoms = atom->natoms;

  double cell[3][3];
  double t_cell[3][3];
  get_cell(cell);
  MathExtra::transpose3(cell, t_cell);

  if (dataset) spg_free_dataset(dataset);
  dataset = spg_get_dataset(t_cell, (const double (*)[3])all_positions, all_types, natoms, symprec);

  if (dataset == nullptr) {
    std::string spg_error_msg = spg_get_error_message(spg_get_error_code());
    error->all(FLERR, spg_error_msg.c_str());
  }

  if ( do_print || spacegroup_number != dataset->spacegroup_number) {
    print_symmetry(spacegroup_number);
  }

  spacegroup_number = dataset->spacegroup_number;

  if (do_find_prim) {
    int find_prim = spg_find_primitive(t_cell, (double (*)[3])all_positions, all_types, natoms, symprec);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        prim_cell[k][j] = t_cell[j][k]; //transposed
      }
    }

    // Atom mapping for primmitive
    mapping_to_primitive.clear();
    std_mapping_to_primitive.clear();

    for (int i = 0; i < atom->nlocal; i++) {
      int id = atom->tag[i] - 1;
      mapping_to_primitive.push_back(dataset->mapping_to_primitive[id]);
    }

    int n_std = dataset->n_std_atoms;
    for (int i = 0; i < n_std; i++) {
      std_mapping_to_primitive.push_back(dataset->std_mapping_to_primitive[i]);
    }
  }
}


void FixSymmetry::prep_symmetry() {
  // Get symmetry operations
  int natoms = atom->natoms;
  int nsym = dataset->n_operations;
  if (rotation_matrices) memory->destroy(rotation_matrices);
  if (translation_vectors) memory->destroy(translation_vectors);
  memory->create(rotation_matrices, nsym, 3, 3, "fix_symmetry:rotation_matrices");
  memory->create(translation_vectors, nsym, 3, "fix_symmetry:translation_vectors");
  symm_map.clear();


  for (int k = 0; k < nsym; k++) {
    double R[3][3];
    double t[3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        rotation_matrices[k][i][j] = dataset->rotations[k][i][j];
        R[i][j] = dataset->rotations[k][i][j];
      }
      translation_vectors[k][i] = dataset->translations[k][i];
      t[i] = dataset->translations[k][i];
    }
    auto this_op_map = std::vector<int>(natoms, -1);

    for (int i = 0; i < natoms; i++) {
      double p[3] = {all_positions[3*i], all_positions[3*i+1], all_positions[3*i+2]};
      double new_p[3];
      MathExtra::matvec(R, p, new_p);
      new_p[0] += t[0];
      new_p[1] += t[1];
      new_p[2] += t[2];

      std::vector<double[3]> dp(natoms);
      for (int j = 0; j < natoms; j++) {
        dp[j][0] = all_positions[3*j]   - new_p[0];
        dp[j][1] = all_positions[3*j+1] - new_p[1];
        dp[j][2] = all_positions[3*j+2] - new_p[2];
        dp[j][0] -= std::round(dp[j][0]); // Ensure periodic boundary conditions
        dp[j][1] -= std::round(dp[j][1]); // Ensure periodic boundary conditions
        dp[j][2] -= std::round(dp[j][2]); // Ensure periodic boundary conditions
      }
      int min_index = 0;
      double min_delta = 9999.0;
      for (int j = 0; j < natoms; j++) {
        double delta = sqrt(dp[j][0]*dp[j][0] + dp[j][1]*dp[j][1] + dp[j][2]*dp[j][2]);
        if (delta < min_delta) {
          min_index = j;
          min_delta = delta;
        }
      }
      this_op_map[i] = min_index;
    }
    symm_map.push_back(this_op_map);
  }
}

void FixSymmetry::symmetrize_rank1(std::vector<double[3]> &vec, double cell[3][3], double inv_cell[3][3]) {
  int nlocal = atom->nlocal;
  int nsym = dataset->n_operations;

  std::vector<double[3]> scaled_vec(nlocal);

  for (int i = 0; i < nlocal; i++) {
    MathExtra::transpose_matvec(inv_cell, vec[i], scaled_vec[i]);
  }

  // Initialize symmetrized vectors
  std::vector<double[3]> sym_vec(nlocal);
  std::vector<double[3]> transformed_vec(nlocal);
  for (int i = 0; i < nlocal; i++) {
    sym_vec[i][0] = 0.0;
    sym_vec[i][1] = 0.0;
    sym_vec[i][2] = 0.0;
    transformed_vec[i][0] = 0.0;
    transformed_vec[i][1] = 0.0;
    transformed_vec[i][2] = 0.0;
  }

  for (int k = 0; k < nsym; k++) {
    double R[3][3];

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        R[i][j] = rotation_matrices[k][i][j];
      }
    }
    for (int i = 0; i < nlocal; i++) {
      MathExtra::matvec(R, scaled_vec[i], transformed_vec[i]);
    }
    for (int i = 0; i < nlocal; i++) {
      int iglobal = atom->tag[i] - 1;
      int jglobal = symm_map[k][iglobal];
      for (int j = 0; j < nlocal; j++) {
        if (jglobal == atom->tag[j] - 1) {
          sym_vec[j][0] += transformed_vec[i][0];
          sym_vec[j][1] += transformed_vec[i][1];
          sym_vec[j][2] += transformed_vec[i][2];
        }
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

void FixSymmetry::symmetrize_rank2(double tensor[3][3], double cell[3][3], double inv_cell[3][3]) {
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

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
      R[i][j] = rotation_matrices[k][i][j];
      }
    }

    MathExtra::transpose3(R, Rt);
    MathExtra::times3(Rt, tmp1, tmp2);
    MathExtra::times3(tmp2, R, tmp3);

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
