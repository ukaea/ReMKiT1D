# Multi-architecture (linux/amd64, linux/arm64) dependency image for ReMKiT1D v1.4.0
FROM ubuntu:24.04

LABEL org.opencontainers.image.title="ReMKiT1D CI / Dependency Image" \
      org.opencontainers.image.description="Multi-architecture Ubuntu 24.04 environment with MPICH, PETSc+HYPRE, pFUnit, json-fortran, HDF5, and SUNDIALS for building and testing ReMKiT1D" \
      org.opencontainers.image.source="https://github.com/ukaea/ReMKiT1D" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.version="1.4.0" \
      org.opencontainers.image.authors="Stefan Mijin <stefan.mijin@ukaea.uk>"

ARG DEBIAN_FRONTEND=noninteractive

ARG INSTALL_PREFIX=/opt/remkit1d
ARG MPICH_VERSION=4.2.3
ARG PETSC_VERSION=v3.23.7
ARG PFUNIT_VERSION=v4.12.0
ARG JSON_FORTRAN_VERSION=9.0.5
ARG HDF5_VERSION=1.14.6
ARG SUNDIALS_VERSION=v7.4.0

# Install compilers, build tools, and coverage utilities
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        gcc \
        g++ \
        cmake \
        pkg-config \
        python3 \
        python3-dev \
        libblas-dev \
        liblapack-dev \
        m4 \
        git \
        wget \
        ca-certificates \
        gcovr \
        lcov \
    && rm -rf /var/lib/apt/lists/*

ENV FC=gfortran \
    CC=gcc \
    CXX=g++

# 1. MPICH (portable ch3 device works reliably across amd64 and arm64 containers)
RUN set -eux; \
    mkdir -p /tmp/build-mpich "${INSTALL_PREFIX}"; \
    cd /tmp/build-mpich; \
    wget -qO mpich.tar.gz "https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz"; \
    tar xfz mpich.tar.gz; \
    cd "mpich-${MPICH_VERSION}"; \
    ./configure \
        CC="${CC}" CXX="${CXX}" FC="${FC}" \
        FFLAGS="-fallow-argument-mismatch" \
        --prefix="${INSTALL_PREFIX}/mpich" \
        --with-device=ch3; \
    make -j"$(nproc)"; \
    make install; \
    rm -rf /tmp/build-mpich

ENV PATH="${INSTALL_PREFIX}/mpich/bin:${PATH}" \
    LD_LIBRARY_PATH="${INSTALL_PREFIX}/mpich/lib"

# 2. PETSc with HYPRE (prefix install — no PETSC_ARCH required)
RUN set -eux; \
    git clone --depth 1 --branch "${PETSC_VERSION}" https://gitlab.com/petsc/petsc.git /tmp/build-petsc; \
    cd /tmp/build-petsc; \
    ./configure \
        --prefix="${INSTALL_PREFIX}/petsc" \
        --with-mpi-dir="${INSTALL_PREFIX}/mpich" \
        --download-hypre=1 \
        --download-fblaslapack=1 \
        --with-debugging=0 \
        COPTFLAGS="-O3" \
        CXXOPTFLAGS="-O3" \
        FOPTFLAGS="-O3"; \
    make -j"$(nproc)" all; \
    make install; \
    rm -rf /tmp/build-petsc

# 3. pFUnit (for unit testing)
RUN set -eux; \
    git clone --depth 1 --branch "${PFUNIT_VERSION}" --recursive https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git /tmp/build-pfunit; \
    cmake -S /tmp/build-pfunit -B /tmp/build-pfunit/build \
        -DSKIP_OPENMP=YES \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}/pFUnit"; \
    cmake --build /tmp/build-pfunit/build -j"$(nproc)"; \
    cmake --install /tmp/build-pfunit/build; \
    rm -rf /tmp/build-pfunit

# 4. json-fortran
RUN set -eux; \
    git clone --depth 1 --branch "${JSON_FORTRAN_VERSION}" https://github.com/jacobwilliams/json-fortran.git /tmp/build-jsonfortran; \
    cmake -S /tmp/build-jsonfortran -B /tmp/build-jsonfortran/build \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DSKIP_DOC_GEN=TRUE \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}/json-fortran"; \
    cmake --build /tmp/build-jsonfortran/build -j"$(nproc)"; \
    cmake --install /tmp/build-jsonfortran/build; \
    rm -rf /tmp/build-jsonfortran

# 5. HDF5 (Fortran interface)
RUN set -eux; \
    mkdir -p /tmp/build-hdf5; \
    cd /tmp/build-hdf5; \
    wget -qO hdf5.tar.gz "https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz"; \
    tar xfz hdf5.tar.gz; \
    cmake -S "/tmp/build-hdf5/hdf5-${HDF5_VERSION}" -B /tmp/build-hdf5/build \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DHDF5_BUILD_FORTRAN=ON \
        -DHDF5_BUILD_CPP_LIB=OFF \
        -DHDF5_BUILD_TOOLS=OFF \
        -DHDF5_BUILD_EXAMPLES=OFF \
        -DBUILD_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}/hdf5"; \
    cmake --build /tmp/build-hdf5/build -j"$(nproc)"; \
    cmake --install /tmp/build-hdf5/build; \
    rm -rf /tmp/build-hdf5

# 6. SUNDIALS (MPI + LAPACK + Fortran 2003 interface)
RUN set -eux; \
    git clone --depth 1 --branch "${SUNDIALS_VERSION}" https://github.com/LLNL/sundials.git /tmp/build-sundials; \
    cmake -S /tmp/build-sundials -B /tmp/build-sundials/build \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DENABLE_MPI=ON \
        -DENABLE_LAPACK=ON \
        -DBUILD_FORTRAN_MODULE_INTERFACE=ON \
        -DBUILD_TESTING=OFF \
        -DEXAMPLES_ENABLE_C=OFF \
        -DEXAMPLES_ENABLE_F2003=OFF \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}/sundials"; \
    cmake --build /tmp/build-sundials/build -j"$(nproc)"; \
    cmake --install /tmp/build-sundials/build; \
    rm -rf /tmp/build-sundials

# Export library discovery and runtime linker paths for ReMKiT1D CMake builds
ENV PETSC_DIR="${INSTALL_PREFIX}/petsc" \
    PFUNIT_DIR="${INSTALL_PREFIX}/pFUnit" \
    HDF5_ROOT="${INSTALL_PREFIX}/hdf5" \
    SUNDIALS_DIR="${INSTALL_PREFIX}/sundials/lib/cmake/sundials" \
    PKG_CONFIG_PATH="${INSTALL_PREFIX}/petsc/lib/pkgconfig" \
    CMAKE_PREFIX_PATH="${INSTALL_PREFIX}/hdf5:${INSTALL_PREFIX}/json-fortran/jsonfortran-gnu-${JSON_FORTRAN_VERSION}:${INSTALL_PREFIX}/pFUnit:${INSTALL_PREFIX}/sundials" \
    LD_LIBRARY_PATH="${INSTALL_PREFIX}/mpich/lib:${INSTALL_PREFIX}/petsc/lib:${INSTALL_PREFIX}/hdf5/lib:${INSTALL_PREFIX}/json-fortran/jsonfortran-gnu-${JSON_FORTRAN_VERSION}/lib:${INSTALL_PREFIX}/sundials/lib"

WORKDIR /workspace
