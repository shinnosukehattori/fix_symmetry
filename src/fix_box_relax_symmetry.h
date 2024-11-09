/*
  * fix_symmetry.h
  *
  *  Created on: 2024/09/29 shattori
*/

#ifdef FIX_CLASS
// clang-format off
FixStyle(box/relax/symmetry,FixBoxRelaxSymmetry);
// clang-format on
#else

#ifndef LMP_FIX_BOX_RELAX_SYMMETRY_H
#define LMP_FIX_BOX_RELAX_SYMMETRY_H

#include "fix_box_relax.h"
#include "lammps.h"
#include <vector>
#include <array>

#include "spglib.h"

namespace LAMMPS_NS {


class ARGRemover {
 public:
  ARGRemover(class LAMMPS *, int, char **);
  ~ARGRemover();
 protected:
  int narg_new;
  char **arg_new;
  double symprec;
  double symcell;
  double symposs;
  bool debug;
};

class FixBoxRelaxSymmetry : protected ARGRemover, public FixBoxRelax {
 public:
  FixBoxRelaxSymmetry(class LAMMPS *, int, char **);
  ~FixBoxRelaxSymmetry();
  int setmask() override;
  void init() override;
  void setup_pre_force(int vlag) override;
  void post_run() override;
  double min_energy(double *fextra) override;

 private:
  int spacegroup_number;

  double prim_cell[3][3];
  double sym_cell[3][3];

  double ***rotation_matrices;
  double **translation_vectors;
  std::vector<std::vector<int>> symm_map;

  SpglibDataset *dataset;
  double **prev_positions;
  imageint *prev_image;

  double *all_positions;
  int *all_types;

  std::vector<int> mapping_to_primitive;
  std::vector<int> std_mapping_to_primitive;

  // Symmetrization functions
  int get_index(std::vector<int> &vec, int val);
  void get_cell(double cell[3][3]);
  void get_cell(double cell[3][3], double inv_cell[3][3]);
  void set_cell(double cell[3][3]);
  void save_prev_info();
  void unmap_inv(double cell[3][3], double x[3], imageint image);

  void adjust_cell();
  void adjust_positions(double cell[3][3], double inv_cell[3][3]);
  void adjust_forces(double cell[3][3], double inv_cell[3][3]);
  void adjust_stress(double cell[3][3], double inv_cell[3][3]);

  void print_symmetry(SpglibDataset*, int, bool);
  void refine_symmetry();
  void symmetrize_cell();
  void symmetrize_positions();

  void store_all_coordinates();
  void prep_symmetry();
  void check_symmetry(bool, bool, bool);

  void symmetrize_rank1(std::vector<double[3]> &vec, double cell[3][3], double inv_cell[3][3]);
  void symmetrize_rank2(double vec[3][3], double cell[3][3], double inv_cell[3][3]);

};

}

#endif
#endif

