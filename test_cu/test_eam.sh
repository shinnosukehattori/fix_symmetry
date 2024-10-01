#!/bin/bash

module load oneapi
LAMMPS_DIR="../../mylammps"
LAMMPS_SRC_DIR="$LAMMPS_DIR/src"
dirname="build_fix_symmetry"


cp "$LAMMPS_SRC_DIR/../potentials/Cu_u3.eam" .
#LAMMPS_DIR/$dirname/lmp -h
$LAMMPS_DIR/$dirname/lmp -in in.fix_symmetry_test
