#
# ReMKiT1DDependencies.cmake
# This file is part of ReMKiT1D.
#
# ReMKiT1D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# ReMKiT1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with ReMKiT1D. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2023 United Kingdom Atomic Energy Authority (stefan.mijin@ukaea.uk)
#

# Shared dependency discovery used both during the ReMKiT1D build and by
# installed ReMKiT1DConfig.cmake so downstream consumers see the exact same
# imported targets (MPI::MPI_Fortran, PkgConfig::PETSC, ReMKiT1D::HDF5,
# ReMKiT1D::jsonfortran, SUNDIALS::*).

MACRO(remkit1d_find_dependencies)
    CMAKE_PARSE_ARGUMENTS(_RMK_DEP "AS_DEPENDENCY" "" "" ${ARGN})

    IF(_RMK_DEP_AS_DEPENDENCY)
        INCLUDE(CMakeFindDependencyMacro)
        find_dependency(MPI)
        find_dependency(PkgConfig)
    ELSE()
        FIND_PACKAGE(MPI REQUIRED)
        FIND_PACKAGE(PkgConfig REQUIRED)
    ENDIF()

    IF(DEFINED ENV{PETSC_DIR} AND NOT "$ENV{PETSC_DIR}" STREQUAL "")
        IF(DEFINED ENV{PETSC_ARCH} AND NOT "$ENV{PETSC_ARCH}" STREQUAL "")
            SET(_petsc_prefix "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")
        ELSE()
            SET(_petsc_prefix "$ENV{PETSC_DIR}")
        ENDIF()
        IF(DEFINED ENV{PKG_CONFIG_PATH} AND NOT "$ENV{PKG_CONFIG_PATH}" STREQUAL "")
            SET(ENV{PKG_CONFIG_PATH} "${_petsc_prefix}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
        ELSE()
            SET(ENV{PKG_CONFIG_PATH} "${_petsc_prefix}/lib/pkgconfig")
        ENDIF()
    ENDIF()
    IF(NOT TARGET PkgConfig::PETSC)
        pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc>=3.17)
    ENDIF()

    IF(_RMK_DEP_AS_DEPENDENCY)
        find_dependency(HDF5 COMPONENTS Fortran)
    ELSE()
        FIND_PACKAGE(HDF5 REQUIRED COMPONENTS Fortran)
    ENDIF()
    IF(NOT TARGET ReMKiT1D::HDF5)
        ADD_LIBRARY(ReMKiT1D::HDF5 INTERFACE IMPORTED)
        IF(TARGET HDF5::HDF5)
            TARGET_LINK_LIBRARIES(ReMKiT1D::HDF5 INTERFACE HDF5::HDF5)
        ELSE()
            TARGET_LINK_LIBRARIES(ReMKiT1D::HDF5 INTERFACE ${HDF5_Fortran_LIBRARIES})
            IF(DEFINED HDF5_Fortran_INCLUDE_DIRS)
                TARGET_INCLUDE_DIRECTORIES(ReMKiT1D::HDF5 INTERFACE ${HDF5_Fortran_INCLUDE_DIRS})
            ELSE()
                TARGET_INCLUDE_DIRECTORIES(ReMKiT1D::HDF5 INTERFACE ${HDF5_INCLUDE_DIRS})
            ENDIF()
        ENDIF()
    ENDIF()

    # Support json-fortran >= 8.2.5 (including 9.x) across static, shared, and namespaced targets
    SET(_jsonfortran_pkg "jsonfortran-${CMAKE_Fortran_COMPILER_ID}")
    FIND_PACKAGE(${_jsonfortran_pkg} QUIET)
    IF(NOT ${_jsonfortran_pkg}_FOUND)
        STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" _fc_id_lower)
        SET(_jsonfortran_pkg "jsonfortran-${_fc_id_lower}")
        IF(_RMK_DEP_AS_DEPENDENCY)
            find_dependency(${_jsonfortran_pkg})
        ELSE()
            FIND_PACKAGE(${_jsonfortran_pkg} REQUIRED)
        ENDIF()
    ENDIF()
    IF(DEFINED ${_jsonfortran_pkg}_VERSION AND ${_jsonfortran_pkg}_VERSION VERSION_LESS 8.2.5)
        MESSAGE(FATAL_ERROR "json-fortran >= 8.2.5 is required (found ${${_jsonfortran_pkg}_VERSION})")
    ENDIF()
    IF(NOT TARGET ReMKiT1D::jsonfortran)
        ADD_LIBRARY(ReMKiT1D::jsonfortran INTERFACE IMPORTED)
        IF(TARGET jsonfortran::jsonfortran)
            TARGET_LINK_LIBRARIES(ReMKiT1D::jsonfortran INTERFACE jsonfortran::jsonfortran)
        ELSEIF(TARGET jsonfortran-static)
            TARGET_LINK_LIBRARIES(ReMKiT1D::jsonfortran INTERFACE jsonfortran-static)
        ELSEIF(TARGET jsonfortran)
            TARGET_LINK_LIBRARIES(ReMKiT1D::jsonfortran INTERFACE jsonfortran)
        ELSE()
            MESSAGE(FATAL_ERROR "Found ${_jsonfortran_pkg} but no supported jsonfortran target was exported")
        ENDIF()
        IF(DEFINED jsonfortran_INCLUDE_DIRS)
            TARGET_INCLUDE_DIRECTORIES(ReMKiT1D::jsonfortran INTERFACE ${jsonfortran_INCLUDE_DIRS})
        ENDIF()
    ENDIF()

    IF(_RMK_DEP_AS_DEPENDENCY)
        find_dependency(SUNDIALS 7.0)
    ELSE()
        FIND_PACKAGE(SUNDIALS 7.0 REQUIRED)
    ENDIF()
ENDMACRO()
