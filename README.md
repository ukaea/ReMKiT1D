# ReMKiT1D
![CI](https://github.com/ukaea/ReMKiT1D/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/ukaea/ReMKiT1D/branch/master/graph/badge.svg?token=I709666D08)](https://codecov.io/gh/ukaea/ReMKiT1D)

Contact: stefan.mijin@ukaea.uk


## **Re**active **M**ultifluid and **Ki**netic **T**ransport in **1D**

## Overview 

ReMKiT1D is a framework for building 1D multi-fluid models of the tokamak Scrape-Off Layer with support for kinetic electron effects and collisional-radiative modelling. The framework is controlled through a Python interface available [here](https://github.com/ukaea/ReMKiT1D-Python). The global user manual, covering both the framework and the Python wrappers, is available on the Python repository. Requirements and installation instructions for ReMKiT1D can also be found [here](https://ukaea.github.io/ReMKiT1D/)

For an overview of both the Fortran framework and the Python interface see the code paper (available soon).
## Prerequisites

Building and compiling: 
- CMake 3.14 or higher (CMake 3.19 found to have issues with HDF5)
- gfortran-11 (fully tested, current default)
- ifort (should work, but not used in current version)

Libraries: 
- MPI - tested with mpich-3.4.1
- HDF5 
- pFUnit - only if unit testing required 
- PETSc - tested with versions 3.16 and 3.17 (3.18 not compatible) - requires hypre
- json-fortran - currently using version 8.2.5

Make sure that all libraries have been compiled using the same compiler/mpi distribution. 

The code has so far been tested on Ubuntu and MacOS. 

## Building prerequisites 

On Ubuntu, the following commands can be used to build all of the prerequisites. Keep in mind conflicts might occur depending on the state of your system. First run

```
sudo apt-get install g++-11 gcc-11 gfortran-11 libblas-dev liblapack-dev m4 gcovr python-pip
sudo pip install ford
```

Make sure that the correct CMake version is installed. For example, to install CMake 3.18.0:

```
wget https://github.com/Kitware/CMake/releases/download/v3.18.0/cmake-3.18.0.tar.gz
tar -zvxf cmake-3.18.0.tar.gz
cd ./cmake-3.18.0
./bootstrap -- -DCMAKE_USE_OPENSSL=OFF
make
make install
```

Then, to install a fresh build of MPI:

```
wget -O mpich.tar.gz https://www.mpich.org/static/downloads/3.4.2/mpich-3.4.2.tar.gz
tar xfz mpich.tar.gz
mkdir mpich-install
cd ./mpich-3.4.2
./configure CC=gcc-11 CXX=g++-11 FC=gfortran-11 --prefix=/home/mpich-install  FFLAGS=-fallow-argument-mismatch --with-device=ch3
make 
make install
```

To install a compatible version of PETSc (make sure to set correct paths!):

```
git clone -b release-3.17 https://gitlab.com/petsc/petsc.git petsc
cd petsc
./configure --with-mpi-dir=/path/to/mpich-install --download-hypre=1 --download-fblaslapack=1 --with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3
make PETSC_DIR=/home/petsc PETSC_ARCH=arch-linux-c-opt all check
```

To install pFUnit (make sure to set paths):
```
git clone -b v4.2.2 https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
cd pFUnit
mkdir build
cd build
cmake .. -DSKIP_OPENMP=yes -DCMAKE_Fortran_COMPILER=gfortran-11 -DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_INSTALL_PREFIX=/path/to/installs/pFUnit
make install
```

The json-fortran library version used by ReMKiT1D can be installed using:

```
git clone -b 8.2.5 https://github.com/jacobwilliams/json-fortran.git
cd json-fortran
mkdir build
cd build
cmake .. -DCMAKE_Fortran_COMPILER=gfortran-11 -DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_INSTALL_PREFIX=/path/to/installs/json-fortran
make install
```

Install the hdf5 library:

```
git clone -b hdf5-1_13_0 https://github.com/HDFGroup/hdf5.git
cd hdf5
mkdir build
cd build
cmake .. -DCMAKE_Fortran_COMPILER=gfortran-11 -DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11 -DHDF5_BUILD_FORTRAN=ON -DCMAKE_INSTALL_PREFIX=/path/to/installs/hdf5
make install
```

Set the environmental variables

```
export PFUNIT_DIR=/home/installs/pFUnit
export PETSC_DIR=/home/petsc
export PETSC_ARCH=arch-linux-c-opt
export PATH=$PATH:/home/installs/petsc
export PATH=$PATH:/home/installs/hdf5
export PATH=$PATH:/home/installs/json-fortran/jsonfortran-gnu-8.2.5
export PATH=$PATH:/home/mpich-install/bin:$PATH
export LD_LIBRARY_PATH=/home/mpich-install/lib:$LD_LIBRARY_PATH
```
NOTE: Some of the above paths might have to be set differently depending on your system setup and you might have to define the ones you need after the corresponding installation. Another option is to export the paths before the start of the installation process. 

## Building ReMKiT1D 

1. Make sure all used libraries have their paths set (set PETSC_DIR and PETSC_ARCH, for example)
2. Create a debug or build folder. On Linux:
``` 
mkdir build
cd build
```
3. Run cmake (see [here](https://cmake.org/cmake/help/latest/manual/cmake.1.html) for cmake documentation)
```
cmake ..
```
4. Run make and make test (if pFUnit enabled)
```
make
make test
```

In order to avoid the above installation process, a dockerfile is made available for those who prefer to use it. For installation instructions and support on HPC platforms, please contact developers.
## Running ReMKiT1D 

If the above build process is successful, the executable will be in build/src/executables/ReMKiT1D.
ReMKiT1D requires a config.json file as input. Assuming a valid config file is available, the code can be run from the executable directory using, for example:
```
mpirun -np [num_procs] ./ReMKiT1D
```
where [num_procs] should be set to the desired number of processes. An alternative path to the config file can be specified using the command line option 

```
-with_config_path=/path/to/config/file.json
```

## Licence

ReMKiT1D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

ReMKiT1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with ReMKiT1D. If not, see <https://www.gnu.org/licenses/>. 

Copyright 2023 United Kingdom Atomic Energy Authority