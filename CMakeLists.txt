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

CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

PROJECT(ReMKiT1D VERSION 1.1.0
  DESCRIPTION "A Fortran framework for 1D multifluid and kinetic simulations geared towards Scrape-Off Layer plasmas"
  LANGUAGES Fortran)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ReMKiT1D dependencies
FIND_PACKAGE(MPI REQUIRED)

FIND_PACKAGE(PkgConfig REQUIRED) 

SET(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
SET(ENV{PKG_CONFIG_PATH} $ENV{PKG_CONFIG_PATH}:${PETSC}/lib/pkgconfig)
pkg_check_modules(PETSC REQUIRED petsc)
INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})
LINK_DIRECTORIES(${PETSC_LINK_LIBRARIES})

FIND_PACKAGE(PFUNIT)

FIND_PACKAGE(HDF5 REQUIRED COMPONENTS Fortran HL)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})

FIND_PACKAGE(jsonfortran-${CMAKE_Fortran_COMPILER_ID} 8.2.5 REQUIRED)
INCLUDE_DIRECTORIES(${jsonfortran_INCLUDE_DIRS})

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

# This INCLUDE statement executes code that sets the compile flags for
# DEBUG, RELEASE, and TESTING, stored in the variables
# DEFAULT_Fortran_FLAGS_{DEBUG,RELEASE,TESTING}.  The flags for your
# build type are available in the variable DEFAULT_Fortran_FLAGS.  This
# will allow for better granularity of control over compiler options.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
MESSAGE(STATUS "Used compiler flags: ${DEFAULT_Fortran_FLAGS}")
# This INCLUDE statement provides a functions to set up automatic
# documentation
INCLUDE(${CMAKE_MODULE_PATH}/MakeFordDocs.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

INCLUDE(GNUInstallDirs)

# Define the name of the unit test binary
SET(LIBDIR ${CMAKE_INSTALL_FULL_LIBDIR}/${PROJECT_NAME})
SET(INCDIR ${CMAKE_INSTALL_FULL_INCLUDEDIR}/${PROJECT_NAME})

# The source 
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src)

# Add the source for the unit tests
IF(${PFUNIT_FOUND})

    ENABLE_TESTING()
    ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/tests)

ELSE()
    MESSAGE(WARNING "pFUnit not installed, testing disabled")
ENDIF()


# Build the documentation
IF (FORD_PATH)
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

# Set up to generate coverage reports
FIND_PACKAGE(PythonInterp REQUIRED)
find_program(GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/scripts/test)

# Generate a coverage report for unit tests
IF(CMAKE_BUILD_TYPE STREQUAL "TESTING")
  add_custom_target(gcovr
    # Create folder
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/gcovr
    
    # Running gcovr
    COMMAND ${Python_EXECUTABLE} ${GCOVR_PATH} --html --html-details -s
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

  add_custom_target(gcov
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/coverage
    COMMAND echo "=================== GCOV ===================="
    COMMAND bash ${CMAKE_SOURCE_DIR}/get_code_cov.sh ${CMAKE_BINARY_DIR}/coverage ${CMAKE_SOURCE_DIR}/src
    COMMAND echo "-- Coverage files have been output to ${CMAKE_BINARY_DIR}/coverage"
    )
  
ENDIF()
