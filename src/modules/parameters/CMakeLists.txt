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
PROJECT(parameters Fortran)

FILE(GLOB_RECURSE ${PROJECT_NAME}_sources *.f90 *.F90)

MESSAGE(STATUS "Sources for project ${PROJECT_NAME}: ${${PROJECT_NAME}_sources}")

ADD_LIBRARY(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_sources})

TARGET_COMPILE_OPTIONS(${PROJECT_NAME}  PRIVATE -cpp)

IF(${CMAKE_BUILD_TYPE} MATCHES "DEBUG")
TARGET_COMPILE_DEFINITIONS(${PROJECT_NAME} PUBLIC DEBUG)
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "DEBUG")

SET(MODDIR ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(INCDIR ${INCDIR}/${PROJECT_NAME})
SET_PROPERTY(TARGET ${PROJECT_NAME} PROPERTY Fortran_MODULE_DIRECTORY ${MODDIR})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME}
  PUBLIC
    "$<INSTALL_INTERFACE:${INCDIR}>"
    "$<BUILD_INTERFACE:${MODDIR}>")
ADD_LIBRARY(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

INSTALL(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}_installation
  ARCHIVE DESTINATION ${LIBDIR}
  LIBRARY DESTINATION ${LIBDIR}
  PUBLIC_HEADER DESTINATION ${INCDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
INSTALL(DIRECTORY ${MODDIR}/
  DESTINATION ${INCDIR})