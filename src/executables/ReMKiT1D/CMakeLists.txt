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
PROJECT(ReMKiT1D Fortran)

FILE(GLOB_RECURSE ${PROJECT_NAME}_sources *.f90 *.F90)

MESSAGE(STATUS "Sources for project ${PROJECT_NAME}: ${${PROJECT_NAME}_sources}")

# Define the executable in terms of the source files
ADD_EXECUTABLE(${PROJECT_NAME} ${${PROJECT_NAME}_sources})

TARGET_COMPILE_OPTIONS(${PROJECT_NAME} PRIVATE ${DEFAULT_Fortran_FLAGS})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} 
  PUBLIC
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
  PRIVATE
    ${BUILD_TYPE_Fortran_FLAGS}
)

INSTALL(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}_installation
  ARCHIVE DESTINATION ${LIBDIR}
  LIBRARY DESTINATION ${LIBDIR}
  PUBLIC_HEADER DESTINATION ${INCDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
