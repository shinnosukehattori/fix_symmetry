# Installation Instructions for FixSymmetry in LAMMPS

## Overview

`FixSymmetry` is a custom LAMMPS fix that enforces crystal symmetry during simulations, porting from ase. It utilizes the `spglib` library to obtain symmetry operations based on the specified space group.

This guide provides step-by-step instructions to install `FixSymmetry` into LAMMPS and run a test simulation to verify its functionality.

## Prerequisites

- **LAMMPS Source Code**: Ensure you have the LAMMPS source code available on your system.
- **spglib Library**: Install `spglib` (version compatible with your system).

## Installation Steps

### 1. Install spglib

Clone the `spglib` repository and install it:

```bash
git clone https://github.com/spglib/spglib.git
cd spglib
mkdir build
cd build
cmake .. --install-prefix=$HOME/.local
make
sudo make install
```

Ensure that the `spglib` header files and libraries are installed in standard system directories or note their installation paths for later use.

### 2. Modify the instalation related files in lammps for FixSymmetry
Basically, you can use install_fix_symmetry.sh with CMake compiliation in lammps.
If you would like to compile manually with Make, please check below.

### 2-1. Copy FixSymmetry Source Files

Copy `fix_symmetry.h` and `fix_symmetry.cpp` into the LAMMPS `src` directory:

```bash
cp fix_symmetry.h /path/to/lammps/src/
cp fix_symmetry.cpp /path/to/lammps/src/
```

### 2-2. Update LAMMPS Build Configuration
Use
Install_?

Navigate to the LAMMPS `src` directory:

```bash
cd /path/to/lammps/src
```

#### For Makefile Build System:

- **Add Source File**: Append `fix_symmetry.cpp` and `fix_box_relax_symmetry.cpp` to the list of source files. Edit `Makefile` or `Makefile.list` and include:

  ```
  EXTRA_SRCS += fix_symmetry.cpp:fix_box_relax_symmtery.cpp
  ```

- **Include spglib Headers and Libraries**: Add the spglib include and library paths. Edit `Makefile.lammps` and include:

  ```
  EXTRA_INC += -I/path/to/spglib/include
  EXTRA_LIB += -L/path/to/spglib/lib -lspglib
  ```

  Replace `/path/to/spglib/include` and `/path/to/spglib/lib` with the actual paths if `spglib` is not in standard directories.

#### For CMake Build System:

- **Modify CMakeLists.txt**: If using CMake, add `fix_symmetry.cpp` to the source list and link `spglib`.

  In `CMakeLists.txt`, add:

  ```cmake
  set(SOURCES ${SOURCES} fix_symmetry.cpp:fix_box_relax_symmetry.cpp)
  find_package(spglib REQUIRED)
  include_directories(${SPGLIB_INCLUDE_DIRS})
  target_link_libraries(lammps ${SPGLIB_LIBRARIES})
  ```

### 3. Build LAMMPS


#### CMake
```bash
mkdir build
cd build
cmake ../cmake ...


make -j8
```

#### Make
Clean previous builds (optional):

```bash
make clean-all
```

Compile LAMMPS:

```bash
make mpi    # Replace 'mpi' with your desired build target
```

Ensure that the compiler can find `spglib` headers and libraries. If not in standard locations, you may need to set environment variables or adjust compiler flags.

### 5. Verify Installation

After successful compilation, verify that `FixSymmetry` is included:

- Run `lmp_mpi -h` and check if `fix symmetry` and `fix box/relax/symmetry` appears in the list of available fixes.
- Alternatively, proceed to run the test simulation.

## Usage Example without box change

In your LAMMPS input script, use `FixSymmetry` as follows:

```lammps
fix sym all symmetry <tolerance>
```
- `<tolerance>`: Optional. The tolerance for symmetry enforcement (default is `1e-5`).

## Usage Example with box change

In your LAMMPS input script, use `FixSymmetry` as follows:

```lammps
fix sym all box/relax/symmetry symprec <tolerance>
```
- `<tolerance>`: Optional. The tolerance for symmetry enforcement (default is `1e-5`).

## Troubleshooting

- **Compilation Errors**: Ensure that `spglib` is correctly installed and that the compiler and linker can find the necessary headers and libraries.
- **Runtime Errors**: Verify that the space group number and tolerance are correctly specified and appropriate for your simulation.
- **Missing Potentials**: Ensure that all required potential files are available in your working directory.

## License

This code is distributed under the GNU General Public License v2. See the `LICENSE` file for more details.

## Acknowledgments

- **spglib**: [https://github.com/spglib/spglib](https://github.com/spglib/spglib)
- **LAMMPS**: [https://lammps.sandia.gov/](https://lammps.sandia.gov/)


