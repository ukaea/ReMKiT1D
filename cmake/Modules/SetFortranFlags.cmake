#
# SetFortranFlags.cmake
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
# Adapted from cmake_fortran_template
# <https://github.com/SethMMorton/cmake_fortran_template>
# Copyright (c) 2018 Seth M. Morton
#

######################################################
# Determine and set the Fortran compiler & link flags
######################################################

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/SetCompileFlag.cmake")

STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

SET(_rmk_bt_doc "Choose the type of build: DEBUG, RELEASE, TESTING, RelWithDebInfo, or MinSizeRel.")
IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "${_rmk_bt_doc}" FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "${_rmk_bt_doc}" FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET(CMAKE_BUILD_TYPE "TESTING" CACHE STRING "${_rmk_bt_doc}" FORCE)
ELSEIF(BT STREQUAL "RELWITHDEBINFO")
    SET(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "${_rmk_bt_doc}" FORCE)
ELSEIF(BT STREQUAL "MINSIZEREL")
    SET(CMAKE_BUILD_TYPE "MinSizeRel" CACHE STRING "${_rmk_bt_doc}" FORCE)
ELSEIF(NOT BT)
    MESSAGE(STATUS "Build type not explicitly set. Trying to guess...")
    GET_FILENAME_COMPONENT(BINDIR "${PROJECT_BINARY_DIR}" NAME)
    IF(BINDIR MATCHES "debug")
        SET(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "${_rmk_bt_doc}" FORCE)
        SET(BT "DEBUG")
        MESSAGE(STATUS "Set to ${CMAKE_BUILD_TYPE}")
    ELSE()
        SET(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "${_rmk_bt_doc}" FORCE)
        SET(BT "RELEASE")
        MESSAGE(STATUS "Defaulting CMAKE_BUILD_TYPE to ${CMAKE_BUILD_TYPE}")
    ENDIF()
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}' not valid. Choices are DEBUG, RELEASE, TESTING, RelWithDebInfo, or MinSizeRel.")
ENDIF()

SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "RELEASE" "DEBUG" "TESTING" "RelWithDebInfo" "MinSizeRel")

# Individual flag probes are cached by CheckFortranCompilerFlag; compute the
# flag lists for the selected build type directly on each configure.
SET(_rmk_basic_flags "")
SET(BUILD_TYPE_Fortran_FLAGS "")
SET(BUILD_TYPE_Fortran_LINK_OPTIONS "")

#####################
### GENERAL FLAGS ###
#####################

SET_COMPILE_FLAG(_rmk_basic_flags Fortran "-warn all" "/warn:all" "-Wall")
SET_COMPILE_FLAG(_rmk_basic_flags Fortran "-Wextra")
SET_COMPILE_FLAG(_rmk_basic_flags Fortran "-Wsurprising")
SET_COMPILE_FLAG(_rmk_basic_flags Fortran "-Wpedantic")
SET_COMPILE_FLAG(_rmk_basic_flags Fortran "-stand=f18" "-std=f2018" "-std=f2008ts")

########################
### BUILD-TYPE FLAGS ###
########################

IF(BT STREQUAL "RELEASE")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran REQUIRED "-O3" "/O3")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-funroll-loops" "-unroll" "/unroll")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-inline" "/Qinline" "-finline-functions")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-ip" "/Qip")
ELSEIF(BT STREQUAL "DEBUG")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran REQUIRED "-O0" "/Od")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-g" "/debug:full")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-traceback" "/traceback" "-fbacktrace")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-check all" "/check:all" "-fcheck=all")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-fpe-all=0" "/fpe-all=0" "-ffpe-trap=invalid,zero,overflow")
ELSEIF(BT STREQUAL "TESTING")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran REQUIRED "-O0" "/Od")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran REQUIRED "-g" "/debug:full")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "--coverage" "/Qcov-gen")
    IF("--coverage" IN_LIST BUILD_TYPE_Fortran_FLAGS)
        LIST(APPEND BUILD_TYPE_Fortran_LINK_OPTIONS "--coverage")
    ENDIF()
ELSEIF(BT STREQUAL "RELWITHDEBINFO")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran REQUIRED "-O2" "/O2")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran "-g" "/debug:full")
ELSEIF(BT STREQUAL "MINSIZEREL")
    SET_COMPILE_FLAG(BUILD_TYPE_Fortran_FLAGS Fortran REQUIRED "-Os" "-O2" "/Os")
ENDIF()

SET(DEFAULT_Fortran_FLAGS ${_rmk_basic_flags} ${BUILD_TYPE_Fortran_FLAGS})
