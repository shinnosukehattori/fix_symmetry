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
  ~FixSymmetry() override;
  int setmask() override;
  void init() override;
  void setup_pre_force(int vlag) override;
  void end_of_step() override;
  void min_post_force(int vflag) override;
  void post_run() override;

 private:
  int spacegroup_number;
  double symprec;  // Tolerance for position displacement
  bool symcell;   // Symmetrize cell
  bool symposs;  // Symmetrize position
  bool debug;

  double sym_cell[3][3];
  double prim_cell[3][3];

  double ***rotation_matrices;
  double **translation_vectors;
  std::vector<std::vector<int>> symm_map;

  SpglibDataset *dataset;
  double **prev_positions;
  double *all_positions;
  int *all_types;

  std::vector<int> mapping_to_primitive;
  std::vector<int> std_mapping_to_primitive;

  // Symmetrization functions
  int get_index(std::vector<int> &vec, int val);
  void get_cell(double cell[3][3]);
  void get_cell(double cell[3][3], double inv_cell[3][3]);
  void set_cell(double cell[3][3]);
  void save_prev_position();

  void adjust_cell();
  void adjust_positions(double cell[3][3], double inv_cell[3][3]);
  void adjust_forces(double cell[3][3], double inv_cell[3][3]);
  void adjust_stress(double cell[3][3], double inv_cell[3][3]);

  void print_symmetry(int);
  void refine_symmetry();
  void symmetrize_cell();
  void symmetrize_positions();

  void store_all_coordinates();
  void prep_symmetry();
  void check_symmetry(bool, bool);

  void symmetrize_rank1(std::vector<double[3]> &vec, double cell[3][3], double inv_cell[3][3]);
  void symmetrize_rank2(double vec[3][3], double cell[3][3], double inv_cell[3][3]);

};

}

#endif
#endif

