Title: Building and Installing
Author: Stefan Mijin 
Date: 16.05.2023.

## Building prerequisites

On Ubuntu, the following commands can be used to build all of the prerequisites. Keep in mind conflicts might occur depending on the state of your system. First run

```
sudo apt-get install g++-11 gcc-11 gfortran-11 libblas-dev liblapack-dev m4 gcovr python-pip
sudo pip install ford
```

Make sure that the correct CMake version is installed (there are possible conflicts with newer versions). For example, to install CMake 3.18.0:

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
export LD_LIBRARY_PATH=/path/to/mpich-install/lib:$LD_LIBRARY_PATH
export PFUNIT_DIR=/home/installs/pFUnit
export PETSC_DIR=/home/petsc
export PETSC_ARCH=arch-linux-c-debug
export PATH=$PATH:/home/installs/petsc
export PATH=$PATH:/home/installs/hdf5
export PATH=$PATH:/home/installs/json-fortran/jsonfortran-gnu-8.2.5
export PATH=/home/mpich-install/bin:$PATH
```
NOTE: Some of the above paths might have to be set differently depending on your system setup. 

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
For test coverage enabled run with `-DCMAKE_BUILD_TYPE=TESTING`

4. Run make and make test (if pFUnit enabled)
```
make
make test
``` 
5. Optional: To generate documentation and get gcovr test coverage data
```
make docs
make gcovr
```

## Setup with Docker

A Dockerfile is included in the repository that sets up an image with all ReMKiT1D prerequisites. 

To build the image with docker run the following command in the docker folder of the repository

```
docker build --tag remkit1d:latest .
``` 
Note that ReMKiT1D needs to be copied into the docker image. This can be done either by [mounting](https://docs.docker.com/storage/volumes/) the local repository folder or using the [docker cp](https://docs.docker.com/engine/reference/commandline/cp/) command. 

The docker cp approach can be done by running 

```
docker run -it --name [CONTAINER_NAME] remkit1d
```
where [CONTAINER_NAME] should be replaced by the name you wish to give the container for use with the docker cp command. Then, in another terminal, run

```
docker cp /path/to/ReMKiT1D_repo [CONTAINER_NAME]:/home/
```
which will copy the ReMKiT1D repository folder into the docker container. After this is done follow the above ReMKiT1D build instructions within the docker image. 