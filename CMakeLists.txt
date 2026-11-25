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
# Adapted from cmake_fortran_template
# <https://github.com/SethMMorton/cmake_fortran_template>
# Copyright (c) 2018 Seth M. Morton
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.18)

PROJECT(ReMKiT1D VERSION 1.4.0
    DESCRIPTION "A Fortran framework for 1D multifluid and kinetic simulations geared towards Scrape-Off Layer plasmas"
    LANGUAGES Fortran C)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

# ReMKiT1D dependencies
INCLUDE(ReMKiT1DDependencies)
remkit1d_find_dependencies()

INCLUDE(CTest)
IF(BUILD_TESTING)
    FIND_PACKAGE(PFUNIT)
ENDIF()

# Compiler/link flags and build helpers
INCLUDE(SetFortranFlags)
MESSAGE(STATUS "Used compiler flags: ${DEFAULT_Fortran_FLAGS}")
INCLUDE(ReMKiT1DModule)
INCLUDE(MakeFordDocs)

############################################################
# Define the actual files and folders that make up the build
############################################################

INCLUDE(GNUInstallDirs)

ADD_SUBDIRECTORY(src)

IF(BUILD_TESTING)
    IF(PFUNIT_FOUND)
        ADD_SUBDIRECTORY(tests)
    ELSE()
        MESSAGE(WARNING "pFUnit not installed, testing disabled")
    ENDIF()
ENDIF()

# Build the documentation
IF(FORD_PATH)
    ADD_FORD_DOCUMENTATION(docs ford-docs.md)
ENDIF()

# Export this project for downstream find_package(ReMKiT1D)
SET(INSTALL_CONFIGDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
INSTALL(EXPORT ${PROJECT_NAME}_installation
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION "${INSTALL_CONFIGDIR}"
)

INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

CONFIGURE_PACKAGE_CONFIG_FILE(
    "${CMAKE_CURRENT_LIST_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${INSTALL_CONFIGDIR}"
)

INSTALL(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/Modules/ReMKiT1DDependencies.cmake"
    DESTINATION "${INSTALL_CONFIGDIR}"
)

# Install the documentation
INSTALL(DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/docs/"
    DESTINATION "${CMAKE_INSTALL_DOCDIR}"
    OPTIONAL MESSAGE_NEVER
)

# Generate coverage reports for unit tests when building in TESTING mode
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" _rmk_bt_upper)
IF(_rmk_bt_upper STREQUAL "TESTING")
    FIND_PACKAGE(Python3 COMPONENTS Interpreter)
    FIND_PROGRAM(GCOVR_PATH gcovr)

    IF(GCOVR_PATH)
        ADD_CUSTOM_TARGET(gcovr
            COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/gcovr"
            COMMAND ${GCOVR_PATH} --html --html-details -s
                    -r "${PROJECT_SOURCE_DIR}" --object-directory="${PROJECT_BINARY_DIR}"
                    -o gcovr/index.html
            WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
            COMMENT "Running gcovr to produce HTML code coverage report."
        )
        ADD_CUSTOM_COMMAND(TARGET gcovr POST_BUILD
            COMMAND ;
            COMMENT "Open ./gcovr/index.html in your browser to view the full coverage report."
        )
    ENDIF()

    ADD_CUSTOM_TARGET(gcov
        COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/coverage"
        COMMAND echo "=================== GCOV ===================="
        COMMAND bash "${PROJECT_SOURCE_DIR}/get_code_cov.sh" "${CMAKE_BINARY_DIR}/coverage" "${PROJECT_SOURCE_DIR}/src"
        COMMAND echo "-- Coverage files have been output to ${CMAKE_BINARY_DIR}/coverage"
    )
ENDIF()
