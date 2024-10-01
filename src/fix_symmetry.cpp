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
#include "spglib.h"

using namespace LAMMPS_NS;


FixSymmetry::FixSymmetry(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (narg < 4 || narg > 5)
    error->all(FLERR, "Illegal fix symmetry command. Must specify space group number and optional tolerance.");

  // Get the space group number
  spacegroup_number = atoi(arg[3]);

  // Get the tolerance (if not specified, set default value)
  if (narg == 5) {
    tolerance = utils::numeric(FLERR, arg[4], false, lmp);
  } else {
    tolerance = 1e-5;  // Default tolerance
  }

  // Initialize previous positions array
  prev_positions = nullptr;
}

FixSymmetry::~FixSymmetry() {
  // Free memory
  if (prev_positions) {
    memory->destroy(prev_positions);
  }
}

int FixSymmetry::setmask() {
  int mask = 0;
  mask |= FixConst::POST_FORCE;  // Apply after force calculation
  mask |= FixConst::MIN_POST_FORCE;  // Apply after force calculation at the end of the minimization step
  return mask;
}

void FixSymmetry::init() {
  // Save previous atom positions
  memory->create(prev_positions, atom->nlocal, 3, "fix_symmetry:prev_positions");
  for (int i = 0; i < atom->nlocal; i++) {
    prev_positions[i][0] = atom->x[i][0];
    prev_positions[i][1] = atom->x[i][1];
    prev_positions[i][2] = atom->x[i][2];
  }
}

void FixSymmetry::post_force(int vflag) {
  // Monitor atom displacement and update symmetry operations if needed
  if (need_to_update_symmetry()) {
    generate_symmetry_operations();
  }

  // Apply symmetrization using precomputed symmetry operations and mappings
  symmetrize_positions();
  symmetrize_cell();
  symmetrize_forces();
  symmetrize_stress();
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
  if (global_max_disp > tolerance) {
    return true;  // Need to recalculate symmetry operations
  }
  return false;  // No need to recalculate
}

void FixSymmetry::symmetrize_positions() {
  double **x = atom->x;
  int nlocal = atom->nlocal;

  // Convert positions to vector format
  std::vector<double[3]> positions(nlocal);
  for (int i = 0; i < nlocal; i++) {
    positions[i][0] = x[i][0];
    positions[i][1] = x[i][1];
    positions[i][2] = x[i][2];
  }

  // Symmetrize positions
  symmetrize_rank1(positions);

  // Update positions
  for (int i = 0; i < nlocal; i++) {
    printf("bef_x[%d] = %f %f %f\n", i, x[i][0], x[i][1], x[i][2]);
    x[i][0] = positions[i][0];
    x[i][1] = positions[i][1];
    x[i][2] = positions[i][2];
    printf("aft_x[%d] = %f %f %f\n", i, x[i][0], x[i][1], x[i][2]);
  }
}

void FixSymmetry::symmetrize_cell() {
  // Get cell matrix
  double h[3][3];
  set_lattice_from_domain(h);

  // Symmetrize cell
  symmetrize_rank2(h);

  // Update cell matrix
  domain->boxlo[0] = domain->boxlo[1] = domain->boxlo[2] = 0.0;
  domain->boxhi[0] = MathExtra::len3(h[0]);
  domain->boxhi[1] = MathExtra::len3(h[1]);
  domain->boxhi[1] = MathExtra::len3(h[2]);
  domain->xy = h[1][0];
  domain->xz = h[2][0];
  domain->yz = h[2][1];

  // Update related parameters
  domain->set_global_box();
  domain->set_local_box();
}

void FixSymmetry::symmetrize_forces() {
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

void FixSymmetry::symmetrize_stress() {
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

void FixSymmetry::symmetrize_rank1(std::vector<double[3]> &vec) {
  int nlocal = atom->nlocal;
  int nsym = symmetry_matrices.size();

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
        R[i][j] = symmetry_matrices[k][i][j];
        printf("R[%d][%d] = %f\n", i, j, R[i][j]);
      }
    }
    std::vector<int> &mapping = symm_map[k];

    for (int i = 0; i < nlocal; i++) {
      int j = mapping[i];
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
  int nsym = symmetry_matrices.size();

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
      R[i][j] = symmetry_matrices[k][i][j];
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

void FixSymmetry::set_lattice_from_domain(double lattice[3][3]) {
  lattice[0][1] = lattice[0][2] = lattice[1][2] = 0.0;
  lattice[0][0] = domain->h[0];
  lattice[1][1] = domain->h[1];
  lattice[2][2] = domain->h[2];
  if (domain->triclinic) {
    lattice[1][0] = lattice[2][0] = lattice[2][1] = 0.0;
  } else {
    lattice[1][0] = domain->xy;
    lattice[2][0] = domain->xz;
    lattice[2][1] = domain->yz;
  }
}

void FixSymmetry::generate_symmetry_operations() {
  // Prepare data structures for spglib
  int num_atoms = atom->natoms;  // Total number of atoms

  // Get lattice constants
  double lattice_spg[3][3];
  set_lattice_from_domain(lattice_spg);
  printf("lattice_spg = %f %f %f\n",  lattice_spg[0][0], lattice_spg[1][1], lattice_spg[2][2]);

  // Collect positions and types from all processes
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;

  // Total number of atoms
  int natoms = atom->natoms;

  // Allocate arrays for positions and types
  double (*all_positions)[3] = new double[natoms][3];
  int *all_types = new int[natoms];

  // Initialize arrays
  for (int i = 0; i < natoms; i++) {
    all_positions[i][0] = 0.0;
    all_positions[i][1] = 0.0;
    all_positions[i][2] = 0.0;
    all_types[i] = 0;
  }

  // Store local positions and types
  for (int i = 0; i < nlocal; i++) {
    double lamda[3];
    //domain->x2lamda(x[i], lamda);

    int id = tag[i] - 1;  // Assuming atom IDs start from 1
    all_positions[id][0] = x[i][0];//lamda[0];
    all_positions[id][1] = x[i][1];//lamda[1];
    all_positions[id][2] = x[i][2];//lamda[2];

    all_types[id] = type[i];
  }

  // Gather data from all processes
  MPI_Allreduce(MPI_IN_PLACE, all_positions, natoms * 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, all_types, natoms, MPI_INT, MPI_SUM, world);


  //add debug code to display all_positions and all_types
  for (int i = 0; i < natoms; i++) {
    printf("pos[%d] = %f %f %f, ", i, all_positions[i][0], all_positions[i][1], all_positions[i][2]);
    printf("type[%d] = %d\n", i, all_types[i]);
  }
  printf("torelance %g\n", tolerance);

  // Get symmetry operations
  if (spacegroup_number <= 0) {
    SpglibDataset *dataset = spg_get_dataset(lattice_spg, all_positions, all_types, natoms, tolerance);
    if (dataset == NULL) {
      error->all(FLERR, "Failed to get symmetry operations from spglib.");
    }
    spacegroup_number = dataset->spacegroup_number;
    printf("spacegroup_number = %d\n", spacegroup_number);
  }

  SpglibDataset *dataset = spg_get_dataset(lattice_spg, all_positions, all_types, natoms, tolerance);
  if (dataset == NULL) {
    error->all(FLERR, "Failed to get symmetry operations from spglib.");
  }

  if(spacegroup_number != dataset->spacegroup_number) {
    std::string str = "Space group number mismatch. Expected " + std::to_string(spacegroup_number) + " but got " + std::to_string(dataset->spacegroup_number);
    //error->all(FLERR, str.c_str());
  }

  int nsym = dataset->n_operations;

  symmetry_matrices.clear();
  translation_vectors.clear();
  symm_map.clear();

  for (int i = 0; i < nsym; i++) {
    double S[3][3];
    double t[3];

    // Symmetry rotation matrix
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        S[j][k] = dataset->rotations[i][j][k];
      }
      t[j] = dataset->translations[i][j];
    }

    symmetry_matrices.push_back({});
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        symmetry_matrices.back()[j][k] = S[j][k];
      }
    }
    translation_vectors.push_back({});
    std::copy(t, t + 3, translation_vectors.back().begin());

    // Atom mapping
    std::vector<int> mapping(nlocal);
    for (int j = 0; j < nlocal; j++) {
      int id = tag[j] - 1;
      int mapped_id = dataset->equivalent_atoms[id];
      mapping[j] = mapped_id;
    }
    symm_map.push_back(mapping);
  }

  // Free memory
  delete[] all_positions;
  delete[] all_types;
  spg_free_dataset(dataset);
}

