#
# SetCompileFlag.cmake
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

#############################################################################
# Try each candidate compiler flag in order and append the first supported
# flag (as list tokens) to FLAGVAR in the caller's scope. Fails with
# FATAL_ERROR if REQUIRED is specified and none of the candidates work.
#
# Usage:
#   SET_COMPILE_FLAG(FLAGVAR (Fortran|C|CXX) [REQUIRED] flag1 [flag2 ...])
#############################################################################

INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckCXXCompilerFlag)
INCLUDE(CheckFortranCompilerFlag)

FUNCTION(SET_COMPILE_FLAG FLAGVAR LANG)
    SET(FLAG_REQUIRED FALSE)
    SET(FLAG_FOUND FALSE)
    SET(FLAGLIST "")
    FOREACH(var ${ARGN})
        STRING(TOUPPER "${var}" UP)
        IF(UP STREQUAL "REQUIRED")
            SET(FLAG_REQUIRED TRUE)
        ELSE()
            LIST(APPEND FLAGLIST "${var}")
        ENDIF()
    ENDFOREACH()

    SET(_current "${${FLAGVAR}}")
    FOREACH(flag ${FLAGLIST})
        STRING(MAKE_C_IDENTIFIER "RMK_${LANG}_FLAG_${flag}" _flag_var)

        IF(LANG STREQUAL "C")
            CHECK_C_COMPILER_FLAG("${flag}" ${_flag_var})
        ELSEIF(LANG STREQUAL "CXX")
            CHECK_CXX_COMPILER_FLAG("${flag}" ${_flag_var})
        ELSEIF(LANG STREQUAL "Fortran")
            IF(flag STREQUAL "--coverage")
                SET(_saved_link_opts "${CMAKE_REQUIRED_LINK_OPTIONS}")
                LIST(APPEND CMAKE_REQUIRED_LINK_OPTIONS "--coverage")
                CHECK_FORTRAN_COMPILER_FLAG("${flag}" ${_flag_var})
                SET(CMAKE_REQUIRED_LINK_OPTIONS "${_saved_link_opts}")
            ELSE()
                CHECK_FORTRAN_COMPILER_FLAG("${flag}" ${_flag_var})
            ENDIF()
        ELSE()
            MESSAGE(FATAL_ERROR "Unknown language in SET_COMPILE_FLAG: ${LANG}")
        ENDIF()

        IF(${${_flag_var}})
            SEPARATE_ARGUMENTS(_flag_tokens UNIX_COMMAND "${flag}")
            LIST(APPEND _current ${_flag_tokens})
            SET(${FLAGVAR} "${_current}" PARENT_SCOPE)
            SET(FLAG_FOUND TRUE)
            BREAK()
        ENDIF()
    ENDFOREACH()

    IF(FLAG_REQUIRED AND NOT FLAG_FOUND)
        MESSAGE(FATAL_ERROR "No supported ${LANG} compile flag found in: ${FLAGLIST}")
    ENDIF()
ENDFUNCTION()
