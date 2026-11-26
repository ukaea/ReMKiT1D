<p align="left">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="https://github.com/ukaea/ReMKiT1D/blob/master/logos/thumbnail_ReMKiT1D_white.png">
      <source media="(prefers-color-scheme: light)" srcset="https://github.com/ukaea/ReMKiT1D/blob/master/logos/thumbnail_ReMKiT1D_black.png">
      <img src="https://github.com/ukaea/ReMKiT1D/blob/master/logos/thumbnail_ReMKiT1D_black.png" alt="ReMKiT1D logo" width="60%">
    </picture>
</p>

# ReMKiT1D - **Re**active **M**ultifluid and **Ki**netic **T**ransport in **1D**

[![CI](https://github.com/ukaea/ReMKiT1D/actions/workflows/CI.yml/badge.svg)](https://github.com/ukaea/ReMKiT1D/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ukaea/ReMKiT1D/branch/master/graph/badge.svg?token=I709666D08)](https://codecov.io/gh/ukaea/ReMKiT1D)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)

Contact: stefan.mijin@ukaea.uk

Logo credits: Vera Oberhauser

## Overview

ReMKiT1D is a Fortran framework for building 1D multi-fluid models of the tokamak Scrape-Off Layer with support for kinetic electron effects and collisional-radiative modelling. The framework is configured and analysed through the companion Python package [ReMKiT1D-Python (`RMK_support`)](https://github.com/ukaea/ReMKiT1D-Python), where the high-level user manual, tutorials, and examples live. API documentation and installation notes for the Fortran codebase are hosted at [ukaea.github.io/ReMKiT1D](https://ukaea.github.io/ReMKiT1D/).

For an overview of both the Fortran framework and the Python interface, see the [code paper](https://doi.org/10.1016/j.cpc.2024.109195) (*Computer Physics Communications*, 300, 109195, 2024).

### Version compatibility (`v1.4.0`)

`ReMKiT1D v1.4.0` is paired with [`ReMKiT1D-Python v2.1.0`](https://github.com/ukaea/ReMKiT1D-Python). Key changes affecting input files and dependencies:

- Requires **PETSc >= 3.17** (COO matrix assembly; compatible with PETSc 3.17 through 3.23+) and **SUNDIALS >= 7.0**.
- The deprecated HYPRE `"euclid"` preconditioner default and matrix term `skipPattern` option have been removed.
- Adds support for PETSc preconditioner reuse (`reusePreconditioner`), Anderson acceleration depth (`andersonDepth`), CRM `energyInterpolationGrid`, and the `filterWithin` calculation-tree node.

## Prerequisites

Building and compiling:

- **CMake** >= 3.18 (3.28+ recommended)
- **Fortran & C compilers** with Fortran 2018 support: `gfortran` >= 11 (tested with `gfortran-11` and modern `gfortran-13+`)
- **pkg-config** (used to locate PETSc)

Libraries:

- **MPI** — MPICH >= 3.4 (tested with MPICH 3.4.2 and 4.2.x)
- **PETSc** >= 3.17 (tested with 3.17–3.23+) built with HYPRE
- **HDF5** with Fortran interface (tested with 1.13 and 1.14 LTS)
- **json-fortran** >= 8.2.5 (tested with 8.2.5 and 9.x)
- **SUNDIALS** >= 7.0 (tested with 7.2–7.4+) with MPI, LAPACK, and Fortran 2003 interface
- **pFUnit** — optional, required only if building unit tests (`-DBUILD_TESTING=ON`; tested with 4.2.2 on GCC 11 and 4.12+ on GCC 13/14)

Make sure that all libraries are compiled using the same compiler and MPI distribution. ReMKiT1D has been tested on Linux (including WSL) and macOS.

## Building prerequisites

On Ubuntu 24.04 LTS (which ships GCC 13 and CMake 3.28+), the following commands build a clean dependency stack under a single prefix directory (`RMK_PREFIX`). Adjust versions or compiler flags if working on an older distribution.

### 1. System packages

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential gfortran cmake pkg-config \
    libblas-dev liblapack-dev m4 git wget ca-certificates \
    gcovr lcov python3 python3-dev
```

*(Optional)* To build the FORD API documentation locally, install `ford` into a virtual environment or via `pipx`:

```bash
python3 -m venv ~/.venvs/ford
~/.venvs/ford/bin/pip install "ford>=7.0,<8"
```

### 2. Choose an installation prefix

```bash
export RMK_PREFIX="${HOME}/remkit1d-deps"
mkdir -p "${RMK_PREFIX}/src"
cd "${RMK_PREFIX}/src"
```

### 3. MPICH

```bash
wget -O mpich.tar.gz https://www.mpich.org/static/downloads/4.2.3/mpich-4.2.3.tar.gz
tar xfz mpich.tar.gz
cd mpich-4.2.3
./configure CC=gcc CXX=g++ FC=gfortran \
    --prefix="${RMK_PREFIX}/mpich" \
    FFLAGS="-fallow-argument-mismatch" \
    --with-device=ch3
make -j"$(nproc)"
make install
export PATH="${RMK_PREFIX}/mpich/bin:${PATH}"
export LD_LIBRARY_PATH="${RMK_PREFIX}/mpich/lib:${LD_LIBRARY_PATH:-}"
cd "${RMK_PREFIX}/src"
```

### 4. PETSc (with HYPRE)

```bash
git clone --depth 1 --branch v3.23.7 https://gitlab.com/petsc/petsc.git petsc
cd petsc
./configure \
    --prefix="${RMK_PREFIX}/petsc" \
    --with-mpi-dir="${RMK_PREFIX}/mpich" \
    --download-hypre=1 \
    --download-fblaslapack=1 \
    --with-debugging=0 \
    COPTFLAGS="-O3" CXXOPTFLAGS="-O3" FOPTFLAGS="-O3"
make -j"$(nproc)" all
make install
cd "${RMK_PREFIX}/src"
```

> **Note:** When PETSc is installed via `--prefix`, only `PETSC_DIR` (or `PKG_CONFIG_PATH`) is needed; `PETSC_ARCH` can be left unset. In-place PETSc builds using both `PETSC_DIR` and `PETSC_ARCH` are equally supported.

### 5. pFUnit *(optional, for unit tests)*

```bash
git clone --depth 1 --branch v4.12.0 --recursive https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
cmake -S pFUnit -B pFUnit/build \
    -DSKIP_OPENMP=YES \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_INSTALL_PREFIX="${RMK_PREFIX}/pFUnit"
cmake --build pFUnit/build -j"$(nproc)"
cmake --install pFUnit/build
```

### 6. json-fortran

```bash
git clone --depth 1 --branch 9.0.5 https://github.com/jacobwilliams/json-fortran.git
cmake -S json-fortran -B json-fortran/build \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DSKIP_DOC_GEN=TRUE \
    -DCMAKE_INSTALL_PREFIX="${RMK_PREFIX}/json-fortran"
cmake --build json-fortran/build -j"$(nproc)"
cmake --install json-fortran/build
```

### 7. HDF5

```bash
wget -O hdf5.tar.gz https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz
tar xfz hdf5.tar.gz
cmake -S hdf5-1.14.6 -B hdf5-1.14.6/build \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_COMPILER=gcc \
    -DHDF5_BUILD_FORTRAN=ON \
    -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_TOOLS=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DBUILD_TESTING=OFF \
    -DCMAKE_INSTALL_PREFIX="${RMK_PREFIX}/hdf5"
cmake --build hdf5-1.14.6/build -j"$(nproc)"
cmake --install hdf5-1.14.6/build
```

### 8. SUNDIALS

```bash
git clone --depth 1 --branch v7.4.0 https://github.com/LLNL/sundials.git
cmake -S sundials -B sundials/build \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DENABLE_MPI=ON \
    -DENABLE_LAPACK=ON \
    -DBUILD_FORTRAN_MODULE_INTERFACE=ON \
    -DBUILD_TESTING=OFF \
    -DEXAMPLES_ENABLE_C=OFF \
    -DEXAMPLES_ENABLE_F2003=OFF \
    -DCMAKE_INSTALL_PREFIX="${RMK_PREFIX}/sundials"
cmake --build sundials/build -j"$(nproc)"
cmake --install sundials/build
```

### 9. Environment variables

```bash
export RMK_PREFIX="${HOME}/remkit1d-deps"
export PETSC_DIR="${RMK_PREFIX}/petsc"
# Set PETSC_ARCH only if using an in-place (non-prefix) PETSc build:
# export PETSC_ARCH="arch-linux-c-opt"
export PFUNIT_DIR="${RMK_PREFIX}/pFUnit"
export HDF5_ROOT="${RMK_PREFIX}/hdf5"
export SUNDIALS_DIR="${RMK_PREFIX}/sundials/lib/cmake/sundials"
export CMAKE_PREFIX_PATH="${RMK_PREFIX}/hdf5:${RMK_PREFIX}/json-fortran/jsonfortran-gnu-9.0.5:${RMK_PREFIX}/pFUnit:${RMK_PREFIX}/sundials:${CMAKE_PREFIX_PATH:-}"
export PKG_CONFIG_PATH="${PETSC_DIR}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export PATH="${RMK_PREFIX}/mpich/bin:${PATH}"
export LD_LIBRARY_PATH="${RMK_PREFIX}/mpich/lib:${PETSC_DIR}/lib:${RMK_PREFIX}/hdf5/lib:${RMK_PREFIX}/json-fortran/jsonfortran-gnu-9.0.5/lib:${RMK_PREFIX}/sundials/lib:${LD_LIBRARY_PATH:-}"
```

> **Note:** On some Linux distributions, CMake installs libraries into `lib64/` rather than `lib/`. Adjust `LD_LIBRARY_PATH` and `SUNDIALS_DIR` accordingly if needed.

## Building, testing, and installing ReMKiT1D

With the environment variables configured, run an out-of-source CMake build from the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build -j"$(nproc)"
```

Run the unit test suite (requires pFUnit):

```bash
(cd build && ctest --output-on-failure)
```

Install the executable, libraries, Fortran module files, and CMake package config (`find_package(ReMKiT1D)`):

```bash
cmake --install build --prefix "${HOME}/remkit1d-install"
```

### Build options

- **`CMAKE_BUILD_TYPE`**:
  - `RELEASE` — optimised production build (default if omitted)
  - `DEBUG` — runtime assertions, bounds/floating-point checks, and debug symbols
  - `TESTING` — debug symbols with `gcov` instrumentation (enables `gcov` and `gcovr` build targets)
  - `RelWithDebInfo` — optimised build with debug symbols
- **`-DBUILD_TESTING=OFF`**: skip pFUnit discovery and unit test targets.
- **Optional targets**:
  - `cmake --build build --target docs` — generate FORD API documentation (when `ford` is on `PATH`)
  - `cmake --build build --target gcov` / `gcovr` — generate coverage data in `TESTING` builds

## Docker

To avoid building dependencies manually, a multi-architecture (`linux/amd64` and `linux/arm64`) Ubuntu 24.04 Dockerfile is provided in [docker/Dockerfile](docker/Dockerfile):

```bash
docker build -t remkit1d-deps:latest docker/
docker run --rm -it -v "$(pwd):/workspace/ReMKiT1D" -w /workspace/ReMKiT1D remkit1d-deps:latest bash
```

A prebuilt CI image is available on Docker Hub as `smijin/remkit1d-ci:latest`. For HPC cluster build advice, please contact the developers.

## Running ReMKiT1D

After building, the executable is located at `build/src/executables/ReMKiT1D/ReMKiT1D` (or `<prefix>/bin/ReMKiT1D` after `cmake --install`).

ReMKiT1D reads a `config.json` file (typically generated via `ReMKiT1D-Python`) from the current working directory:

```bash
mpirun -np [num_procs] ./build/src/executables/ReMKiT1D/ReMKiT1D
```

Supported command-line options:

- `-with_config_path=/path/to/config.json` — read the configuration from an explicit path
- `--version` or `-v` — print the ReMKiT1D version string and exit

## Citation

If you use ReMKiT1D in your research, please cite the code paper (see [CITATION.cff](CITATION.cff)):

> S. Mijin, D. Power, R. Holden, W. Hornsby, D. Moulton, F. Militello, *"ReMKiT1D - A framework for building reactive multi-fluid models of the tokamak scrape-off layer with coupled electron kinetics in 1D"*, Computer Physics Communications **300** (2024) 109195, [doi:10.1016/j.cpc.2024.109195](https://doi.org/10.1016/j.cpc.2024.109195).

## Licence

ReMKiT1D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

ReMKiT1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with ReMKiT1D. If not, see <https://www.gnu.org/licenses/>.

Copyright 2023-2026 United Kingdom Atomic Energy Authority
