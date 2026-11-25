#
# ReMKiT1DModule.cmake
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

# Define, configure, and install a ReMKiT1D Fortran module library target.
#
# Usage:
#   remkit1d_add_module(<name>
#       [PREPROCESS]
#       [NO_STRICT_FLAGS]
#       [EXTRA_SOURCES <file> ...]
#       [LINK_LIBRARIES <target> ...]
#       [COMPILE_DEFINITIONS <def> ...]
#   )
FUNCTION(remkit1d_add_module MODULE_NAME)
    CMAKE_PARSE_ARGUMENTS(
        _RMK_MOD
        "PREPROCESS;NO_STRICT_FLAGS"
        ""
        "EXTRA_SOURCES;LINK_LIBRARIES;COMPILE_DEFINITIONS"
        ${ARGN}
    )

    FILE(GLOB_RECURSE _mod_sources CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.f90"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
    )
    IF(_RMK_MOD_EXTRA_SOURCES)
        LIST(APPEND _mod_sources ${_RMK_MOD_EXTRA_SOURCES})
    ENDIF()

    ADD_LIBRARY(${MODULE_NAME} STATIC ${_mod_sources})
    ADD_LIBRARY(ReMKiT1D::${MODULE_NAME} ALIAS ${MODULE_NAME})

    IF(_RMK_MOD_PREPROCESS)
        SET_TARGET_PROPERTIES(${MODULE_NAME} PROPERTIES Fortran_PREPROCESS ON)
    ENDIF()

    IF(_RMK_MOD_NO_STRICT_FLAGS)
        TARGET_COMPILE_OPTIONS(${MODULE_NAME} PRIVATE ${BUILD_TYPE_Fortran_FLAGS})
    ELSE()
        TARGET_COMPILE_OPTIONS(${MODULE_NAME} PRIVATE ${DEFAULT_Fortran_FLAGS})
    ENDIF()

    IF(BUILD_TYPE_Fortran_LINK_OPTIONS)
        TARGET_LINK_OPTIONS(${MODULE_NAME} PUBLIC ${BUILD_TYPE_Fortran_LINK_OPTIONS})
    ENDIF()

    IF(_RMK_MOD_COMPILE_DEFINITIONS)
        TARGET_COMPILE_DEFINITIONS(${MODULE_NAME} PUBLIC ${_RMK_MOD_COMPILE_DEFINITIONS})
    ENDIF()

    IF(_RMK_MOD_LINK_LIBRARIES)
        TARGET_LINK_LIBRARIES(${MODULE_NAME} PUBLIC ${_RMK_MOD_LINK_LIBRARIES})
    ENDIF()

    SET(_build_moddir "${CMAKE_CURRENT_BINARY_DIR}/include")
    SET(_install_moddir "${CMAKE_INSTALL_INCLUDEDIR}/ReMKiT1D/${MODULE_NAME}")

    SET_TARGET_PROPERTIES(${MODULE_NAME} PROPERTIES
        Fortran_MODULE_DIRECTORY "${_build_moddir}"
    )

    TARGET_INCLUDE_DIRECTORIES(${MODULE_NAME}
        PUBLIC
            "$<BUILD_INTERFACE:${_build_moddir}>"
            "$<INSTALL_INTERFACE:${_install_moddir}>"
    )

    INSTALL(TARGETS ${MODULE_NAME}
        EXPORT ReMKiT1D_installation
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/ReMKiT1D"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/ReMKiT1D"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
    INSTALL(DIRECTORY "${_build_moddir}/"
        DESTINATION "${_install_moddir}"
    )
ENDFUNCTION()
