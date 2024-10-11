#!/bin/bash

module load oneapi

# Paths (please adjust according to your environment)
LAMMPS_DIR="$HOME/app/lammps_build/mylammps"
SPGLIB_INSTALL_DIR="$HOME/.local"
dirname="build_fix_symmetry"

# Check paths
if [ ! -d "$SPGLIB_INSTALL_DIR" ]; then
  echo "spglib install directory not found: $SPGLIB_INSTALL_DIR"
  exit 1
fi
SPGLIB_INCLUDE_DIR="$SPGLIB_INSTALL_DIR/include"
SPGLIB_LIB_DIR="$SPGLIB_INSTALL_DIR/lib64"
CMAKE_MODULE_PATH="$SPGLIB_LIBRARIES/cmake/Spglib"


cp "src/fix_symmetry.h" "$LAMMPS_DIR/src/"
cp "src/fix_symmetry.cpp" "$LAMMPS_DIR/src/"
cp "src/fix_box_relax_symmetry.h" "$LAMMPS_DIR/src/"
cp "src/fix_box_relax_symmetry.cpp" "$LAMMPS_DIR/src/"

cd $LAMMPS_DIR
rm -rf $dirname
mkdir $dirname
cd $dirname

# Update CMakeLists.txt to include fix_symmetry.cpp and link spglib
echo "Updating CMakeLists.txt to include fix_symmetry.cpp and link spglib..."

# Check if fix_symmetry.cpp is already included
if ! grep -q "fix_symmetry.cpp" "../cmake/CMakeLists.txt"; then
  echo "Adding fix_symmetry.cpp to CMakeLists.txt..."
  echo "set(SOURCE_EXTRA \${SOURCE_EXTRA} fix_symmetry.cpp:fix_box_relax_symmetry.cpp)" >> "../cmake/CMakeLists.txt"
fi

# Add spglib include and library directories
if ! grep -q "find_package(Spglib REQUIRED)" "../cmake/CMakeLists.txt"; then
  echo "Finding and linking Spglib..."
  echo "find_package(Spglib REQUIRED)" >> "../cmake/CMakeLists.txt"
  echo "set(SPGLIB_INCLUDE_DIRS ${SPGLIB_INCLUDE_DIRS} $SPGLIB_INCLUDE_DIR)" >> "../cmake/CMakeLists.txt"
  echo "set(SPGLIB_LIBRARIES ${SPGLIB_LIBRARIES} $SPGLIB_LIB_DIR/libsymspg.so)" >> "../cmake/CMakeLists.txt"
  echo "include_directories(\${SPGLIB_INCLUDE_DIRS})" >> "../cmake/CMakeLists.txt"
  echo "target_link_libraries(lmp PUBLIC \${SPGLIB_LIBRARIES})" >> "../cmake/CMakeLists.txt"
fi


# Configure LAMMPS with CMake
cmake ../cmake \
  -C ../cmake/presets/basic.cmake \
  -D CMAKE_BUILD_TYPE=Debug \
  -D CMAKE_CXX_COMPILER=icpx \
  -D CMAKE_CXX_FLAGS="-qopenmp" \
  -D CMAKE_MODULE_PATH=$CMAKE_MODULE_PATH \
  -D BUILD_MPI=on \
  -D BUILD_SHARED_LIBS=yes \
  -D LAMMPS_EXCEPTIONS=on \
  -D PKG_USER-SYMMETRY=on \
  -D PKG_REPLICA=yes \
  -D PKG_EXTRA-MOLECULE=yes -D PKG_REAXFF=yes \
  -D PKG_USER-MISC=on \
  -D FFT=FFTW3 \
  -D CMAKE_INSTALL_PREFIX=$HOME/.local

# Build LAMMPS
make -j8

# Install LAMMPS (optional)
#make install

