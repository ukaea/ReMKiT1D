Title: Building and Installing
Author: Stefan Mijin
Date: 2026

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

On Ubuntu 24.04 LTS (which ships GCC 13 and CMake 3.28+), the following commands build all prerequisites under a single prefix (`RMK_PREFIX`). Keep in mind that versions or compiler flags may need adjustment on older distributions.

First install system packages:

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential gfortran cmake pkg-config \
    libblas-dev liblapack-dev m4 git wget ca-certificates \
    gcovr lcov python3 python3-dev
```

*(Optional)* Install `ford` into a virtual environment to build API documentation locally:

```bash
python3 -m venv ~/.venvs/ford
~/.venvs/ford/bin/pip install "ford>=7.0,<8"
```

Choose a common installation prefix:

```bash
export RMK_PREFIX="${HOME}/remkit1d-deps"
mkdir -p "${RMK_PREFIX}/src"
cd "${RMK_PREFIX}/src"
```

Install MPICH:

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

Install PETSc (with HYPRE):

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

Install pFUnit *(optional, for unit tests)*:

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

Install json-fortran:

```bash
git clone --depth 1 --branch 9.0.5 https://github.com/jacobwilliams/json-fortran.git
cmake -S json-fortran -B json-fortran/build \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DSKIP_DOC_GEN=TRUE \
    -DCMAKE_INSTALL_PREFIX="${RMK_PREFIX}/json-fortran"
cmake --build json-fortran/build -j"$(nproc)"
cmake --install json-fortran/build
```

Install HDF5:

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

Install SUNDIALS:

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

Set the environment variables (note that `PETSC_ARCH` is only needed for in-place PETSc builds without `--prefix`):

```bash
export RMK_PREFIX="${HOME}/remkit1d-deps"
export PETSC_DIR="${RMK_PREFIX}/petsc"
# export PETSC_ARCH="arch-linux-c-opt"
export PFUNIT_DIR="${RMK_PREFIX}/pFUnit"
export HDF5_ROOT="${RMK_PREFIX}/hdf5"
export SUNDIALS_DIR="${RMK_PREFIX}/sundials/lib/cmake/sundials"
export CMAKE_PREFIX_PATH="${RMK_PREFIX}/hdf5:${RMK_PREFIX}/json-fortran/jsonfortran-gnu-9.0.5:${RMK_PREFIX}/pFUnit:${RMK_PREFIX}/sundials:${CMAKE_PREFIX_PATH:-}"
export PKG_CONFIG_PATH="${PETSC_DIR}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export PATH="${RMK_PREFIX}/mpich/bin:${PATH}"
export LD_LIBRARY_PATH="${RMK_PREFIX}/mpich/lib:${PETSC_DIR}/lib:${RMK_PREFIX}/hdf5/lib:${RMK_PREFIX}/json-fortran/jsonfortran-gnu-9.0.5/lib:${RMK_PREFIX}/sundials/lib:${LD_LIBRARY_PATH:-}"
```

## Building ReMKiT1D

1. Make sure all required library paths are exported (`PETSC_DIR`, `PKG_CONFIG_PATH`, `CMAKE_PREFIX_PATH`, `SUNDIALS_DIR`, and `PFUNIT_DIR` if testing).
2. Configure and build out-of-source from the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build -j"$(nproc)"
```

Available build types are `RELEASE`, `DEBUG`, `TESTING`, and `RelWithDebInfo`. Pass `-DBUILD_TESTING=OFF` to disable unit tests when pFUnit is not installed.

3. Run the unit test suite (if pFUnit is enabled):

```bash
ctest --test-dir build --output-on-failure
```

4. *(Optional)* Install the executable, libraries, module files, and CMake package config:

```bash
cmake --install build --prefix "${HOME}/remkit1d-install"
```

5. *(Optional)* Generate FORD documentation or test coverage reports (coverage requires `-DCMAKE_BUILD_TYPE=TESTING`):

```bash
cmake --build build --target docs
cmake --build build --target gcovr
```

## Setup with Docker

A multi-architecture (`linux/amd64` and `linux/arm64`) `docker/Dockerfile` is included in the repository that builds an Ubuntu 24.04 image with all ReMKiT1D prerequisites preinstalled.

Build the image from the repository root:

```bash
docker build --tag remkit1d:latest docker/
```

Mount your local ReMKiT1D checkout into the container and build inside it:

```bash
docker run --rm -it \
    -v /path/to/ReMKiT1D:/workspace/ReMKiT1D \
    -w /workspace/ReMKiT1D \
    remkit1d:latest bash
```

Inside the container, all prerequisite environment variables are already configured, so you can run the `cmake` and `ctest` commands above directly.
