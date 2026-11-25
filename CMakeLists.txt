#
# CMakeLists.txt
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
# This file is adapted from cmake_fortran_template
# <https://github.com/SethMMorton/cmake_fortran_template>
# Copyright (c) 2018 Seth M. Morton
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.18)

PROJECT(ReMKiT1D VERSION 1.4.0
  DESCRIPTION "A Fortran framework for 1D multifluid and kinetic simulations geared towards Scrape-Off Layer plasmas"
  LANGUAGES Fortran C)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ReMKiT1D dependencies
FIND_PACKAGE(MPI REQUIRED)

FIND_PACKAGE(PkgConfig REQUIRED)

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
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc>=3.17)

INCLUDE(CTest)
IF(BUILD_TESTING)
    FIND_PACKAGE(PFUNIT)
ENDIF()

FIND_PACKAGE(HDF5 REQUIRED COMPONENTS Fortran)
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

# Support json-fortran >= 8.2.5 (including 9.x) across static, shared, and namespaced targets
SET(_jsonfortran_pkg "jsonfortran-${CMAKE_Fortran_COMPILER_ID}")
FIND_PACKAGE(${_jsonfortran_pkg} QUIET)
IF(NOT ${_jsonfortran_pkg}_FOUND)
    STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" _fc_id_lower)
    SET(_jsonfortran_pkg "jsonfortran-${_fc_id_lower}")
    FIND_PACKAGE(${_jsonfortran_pkg} REQUIRED)
ENDIF()
IF(DEFINED ${_jsonfortran_pkg}_VERSION AND ${_jsonfortran_pkg}_VERSION VERSION_LESS 8.2.5)
    MESSAGE(FATAL_ERROR "json-fortran >= 8.2.5 is required (found ${${_jsonfortran_pkg}_VERSION})")
ENDIF()
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

FIND_PACKAGE(SUNDIALS 7.0 REQUIRED)

# Add our local modules to the module path
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

# This INCLUDE statement executes code that sets the compile flags for
# DEBUG, RELEASE, and TESTING, stored in the variables
# DEFAULT_Fortran_FLAGS_{DEBUG,RELEASE,TESTING}.  The flags for your
# build type are available in the variable DEFAULT_Fortran_FLAGS.  This
# will allow for better granularity of control over compiler options.
INCLUDE("${PROJECT_SOURCE_DIR}/cmake/Modules/SetFortranFlags.cmake")
MESSAGE(STATUS "Used compiler flags: ${DEFAULT_Fortran_FLAGS}")
# This INCLUDE statement provides a function to set up automatic
# documentation
INCLUDE("${PROJECT_SOURCE_DIR}/cmake/Modules/MakeFordDocs.cmake")

############################################################
# Define the actual files and folders that make up the build
############################################################

INCLUDE(GNUInstallDirs)

# Define the name of the unit test binary
SET(LIBDIR ${CMAKE_INSTALL_FULL_LIBDIR}/${PROJECT_NAME})
SET(INCDIR ${CMAKE_INSTALL_FULL_INCLUDEDIR}/${PROJECT_NAME})

# The source 
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src)

# Add the source for the unit tests
IF(BUILD_TESTING)
    IF(PFUNIT_FOUND)
        ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/tests)
    ELSE()
        MESSAGE(WARNING "pFUnit not installed, testing disabled")
    ENDIF()
ENDIF()


# Build the documentation
IF(FORD_PATH)
    ADD_FORD_DOCUMENTATION(docs ford-docs.md)
ENDIF()
# Export this project for use by others
SET(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
INSTALL(EXPORT ${PROJECT_NAME}_installation
  FILE
    ${PROJECT_NAME}Targets.cmake
  NAMESPACE
    ${PROJECT_NAME}::
  DESTINATION
    ${INSTALL_CONFIGDIR}
  )

INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

CONFIGURE_PACKAGE_CONFIG_FILE(${CMAKE_CURRENT_LIST_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    PATH_VARS INSTALL_CONFIGDIR
)

# Install the config and configversion files
INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

# Install the documentation
INSTALL(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/docs/
  DESTINATION ${CMAKE_INSTALL_DOCDIR}
  OPTIONAL MESSAGE_NEVER
)

# Generate coverage reports for unit tests when building in TESTING mode
IF(CMAKE_BUILD_TYPE STREQUAL "TESTING")
  FIND_PACKAGE(Python3 COMPONENTS Interpreter)
  find_program(GCOVR_PATH gcovr)

  IF(GCOVR_PATH)
    add_custom_target(gcovr
      # Create folder
      COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/gcovr
      # Running gcovr
      COMMAND ${GCOVR_PATH} --html --html-details -s
              -r ${PROJECT_SOURCE_DIR} --object-directory=${PROJECT_BINARY_DIR}
              -o gcovr/index.html
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      COMMENT "Running gcovr to produce HTML code coverage report."
    )
    # Show info where to find the report
    add_custom_command(TARGET gcovr POST_BUILD
      COMMAND ;
      COMMENT "Open ./gcovr/index.html in your browser to view the full coverage report."
    )
  ENDIF()

  add_custom_target(gcov
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/coverage
    COMMAND echo "=================== GCOV ===================="
    COMMAND bash ${PROJECT_SOURCE_DIR}/get_code_cov.sh ${CMAKE_BINARY_DIR}/coverage ${PROJECT_SOURCE_DIR}/src
    COMMAND echo "-- Coverage files have been output to ${CMAKE_BINARY_DIR}/coverage"
  )
ENDIF()
