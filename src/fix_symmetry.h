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
  double tolerance;  // Tolerance for position displacement

  std::vector<std::array<std::array<double, 3>, 3>> symmetry_matrices;
  std::vector<std::array<double, 3>> translation_vectors;
  std::vector<std::vector<int>> symm_map;  // Atom mapping for symmetry operations

  // Previous atom positions for displacement monitoring
  double **prev_positions;

  // Symmetrization functions
  void generate_symmetry_operations();
  void symmetrize_positions();
  void symmetrize_forces();
  void symmetrize_stress();
  void symmetrize_cell();

  void symmetrize_rank1(std::vector<double[3]> &vec);
  void symmetrize_rank2(double vec[3][3]);

  bool need_to_update_symmetry();
  void set_lattice_from_domain(double lattice[3][3]);

};


}

#endif
#endif

