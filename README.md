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

### 2. Copy FixSymmetry Source Files

Copy `fix_symmetry.h` and `fix_symmetry.cpp` into the LAMMPS `src` directory:

```bash
cp fix_symmetry.h /path/to/lammps/src/
cp fix_symmetry.cpp /path/to/lammps/src/
```

### 3. Update LAMMPS Build Configuration
Use
Install_?

Navigate to the LAMMPS `src` directory:

```bash
cd /path/to/lammps/src
```

#### For Makefile Build System:

- **Add Source File**: Append `fix_symmetry.cpp` to the list of source files. Edit `Makefile` or `Makefile.list` and include:

  ```
  EXTRA_SRCS += fix_symmetry.cpp
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
  set(SOURCES ${SOURCES} fix_symmetry.cpp)
  find_package(spglib REQUIRED)
  include_directories(${SPGLIB_INCLUDE_DIRS})
  target_link_libraries(lammps ${SPGLIB_LIBRARIES})
  ```

### 4. Build LAMMPS

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

- Run `lmp_mpi -h` and check if `fix symmetry` appears in the list of available fixes.
- Alternatively, proceed to run the test simulation.

## Running the Test Simulation

### 1. Copy the Test Input Script

Copy the test input script `in.fix_symmetry_test` to a working directory:

```bash
cp in.fix_symmetry_test /path/to/working/directory/
```

### 2. Prepare the Potential File

Ensure that the required EAM potential file `Cu_u3.eam` is available:

- Copy it from the LAMMPS `potentials` directory:

  ```bash
  cp /path/to/lammps/potentials/Cu_u3.eam /path/to/working/directory/
  ```

### 3. Run the Simulation

Navigate to the working directory and execute the simulation:

```bash
cd /path/to/working/directory/
lmp_mpi -in in.fix_symmetry_test
```

Replace `lmp_mpi` with the name of your LAMMPS executable if different.

### 4. Analyze Results

- Check the output to ensure that the simulation runs without errors.
- Verify that the symmetry is enforced by inspecting the dump files or output data.

## Usage Example

In your LAMMPS input script, use `FixSymmetry` as follows:

```lammps
fix sym all symmetry <spacegroup_number> <tolerance>
```

- `<spacegroup_number>`: The international number of the space group (e.g., 225 for Fm-3m).
- `<tolerance>`: Optional. The tolerance for symmetry enforcement (default is `1e-5`).

**Example:**

```lammps
fix sym all symmetry 225 1e-5
```

## Troubleshooting

- **Compilation Errors**: Ensure that `spglib` is correctly installed and that the compiler and linker can find the necessary headers and libraries.
- **Runtime Errors**: Verify that the space group number and tolerance are correctly specified and appropriate for your simulation.
- **Missing Potentials**: Ensure that all required potential files are available in your working directory.

## License

This code is distributed under the GNU General Public License v2. See the `LICENSE` file for more details.

## Acknowledgments

- **spglib**: [https://github.com/spglib/spglib](https://github.com/spglib/spglib)
- **LAMMPS**: [https://lammps.sandia.gov/](https://lammps.sandia.gov/)


