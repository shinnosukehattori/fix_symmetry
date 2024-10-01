/*
  * fix_symmetry.h
  *
  *  Created on: 2024/09/29 shattori
*/

#ifdef FIX_CLASS
// clang-format off
FixStyle(symmetry,FixSymmetry);
// clang-format on
#else

#ifndef LMP_FIX_SYMMETRY_H
#define LMP_FIX_SYMMETRY_H

#include "fix.h"
#include <vector>
#include <array>

#include "spglib.h"

namespace LAMMPS_NS {

class FixSymmetry : public Fix {
 public:
  FixSymmetry(class LAMMPS *, int, char **);
  ~FixSymmetry();
  int setmask();
  void init();
  void post_force(int vflag);
  void min_post_force(int vflag);

 private:
  int spacegroup_number;
  double symprec;  // Tolerance for position displacement
  bool symforce;   // Symmetrize forces
  bool symstress;  // Symmetrize stress
  bool debug;

  double std_cell[3][3];
  double prim_cell[3][3];

  std::vector<std::array<std::array<double, 3>, 3>> rotation_matrices;
  std::vector<std::array<double, 3>> translation_vectors;
  std::vector<std::array<double, 3>> symm_map;

  SpglibDataset *dataset;
  double **prev_positions;
  double **all_positions;
  int *all_types;

  std::vector<int> mapping_to_primitive;
  std::vector<int> std_mapping_to_primitive;

  // Symmetrization functions
  bool need_to_update_symmetry();
  int get_index(std::vector<int> &vec, int val);
  void store_std_cell();
  void restore_std_cell();

  void adjust_cell();
  void adjust_positions();
  void adjust_forces();
  void adjust_stress();

  void x2lambda(const double pos[3], double lambda[3]);

  void print_symmetry();
  void refine_symmetry();
  void check_and_symmetrize_cell();
  void symmetrize_cell();
  void check_and_symmetrize_positions();
  void symmetrize_positions();

  void initial_prep();
  void prep_symmetry();
  void check_symmetry(bool do_find_prim);

  void symmetrize_forces();
  void symmetrize_stress();

  void symmetrize_rank1(std::vector<double[3]> &vec);
  void symmetrize_rank2(double vec[3][3]);

};


}

#endif
#endif

