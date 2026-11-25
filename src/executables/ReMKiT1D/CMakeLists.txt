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

FILE(GLOB_RECURSE ReMKiT1D_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
)

ADD_EXECUTABLE(ReMKiT1D ${ReMKiT1D_sources})

TARGET_COMPILE_OPTIONS(ReMKiT1D PRIVATE ${DEFAULT_Fortran_FLAGS})

IF(BUILD_TYPE_Fortran_LINK_OPTIONS)
    TARGET_LINK_OPTIONS(ReMKiT1D PRIVATE ${BUILD_TYPE_Fortran_LINK_OPTIONS})
ENDIF()

TARGET_LINK_LIBRARIES(ReMKiT1D
    PRIVATE
        parameters
        basic_support
        grid
        variables
        solver_support
        modeller
        mpi_support
        IO_support
        extended_support
        timeloop
        modeller_assembly
)

INSTALL(TARGETS ReMKiT1D
    EXPORT ReMKiT1D_installation
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
