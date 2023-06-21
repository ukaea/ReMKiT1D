#
# MakeFordDocs.cmake
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

# Defines a target for running the FORD automatic documentation
# generator. Optional arguments determine the command line arguments
# when executing FORD and will override any settings provided in the
# project file.
#
# SETUP_TARGET_FOR_COVERAGE_GCOVR_XML(
#     doc_name                        - New target name
#     project_file                    - Path to the FORD project file
#     [ WARN ]                        - Run FORD with the --warn flag
#     [ NOSEARCH ]                    - Don't generate a search index
#                                       for the documentation (leading
#                                       to faster build)
#     [ PAGE_DIR pdir ]               - A directory containing source
#                                       files for pages describing the
#                                       project
#     [ CSS cssfile ]                 - User-supplied CSS file for docs
#     [ REVISION rev_name ]           - Name of particular revision
#                                       being documented
#     [ SRC_DIR directory [...] ]     - Directories containing
#                                       source files to be documented
#     [ EXCLUDE filename [...] ]      - Files to be excluded from
#                                       documentation
#     [ EXCLUDE_DIR directory [...] ] - Directories for which source
#                                       files will not be parsed for
#                                       documentation
#     [ EXTENSIONS ext [...] ]        - Extensions of files containing
#                                       Fortran code to parse for
#                                       documentation
#     [ MACROS macro[=val] [...] ]    - Preprocessor macros/values to
#                                       use when processing source files
# )

include(CMakeParseArguments)

# Check prereqs
find_package(PythonInterp)
find_program(FORD_PATH ford)
find_program(GRAPHVIZ_PATH NAMES dot neato twopi circo gdp sfdp patchwork osage)

if(NOT FORD_PATH)
  message(STATUS "FORD not found! Cannot build documenation!")
else()
  function(list_to_absolute_paths var base_dir)
    set(listvar "")
    foreach(f ${ARGN})
      get_filename_component(absf ${f} ABSOLUTE BASE_DIR ${base_dir})
      list(APPEND listvar ${absf})
    endforeach()
    set(${var} "${listvar}" PARENT_SCOPE)
  endfunction()

  function(add_ford_documentation doc_name project_file)
    set(options WARN NOSEARCH)
    set(oneValueArgs PAGE_DIR CSS REVISION)
    set(multiValueArgs SRC_DIR EXCLUDE EXCLUDE_DIR EXTENSIONS MACROS)
    cmake_parse_arguments(FORD "${options}" "${oneValueArgs}" "${multiValueArgs}" "${ARGN}")

    get_filename_component(project_file_abs ${project_file} ABSOLUTE)
    get_filename_component(projectfile_dir_rel ${project_file} DIRECTORY)

    get_filename_component(projectfile_dir "${projectfile_dir_rel}" ABSOLUTE)
    set(fordArgs "")
    if(FORD_WARN)
      list(APPEND fordArgs "--warn")
    else()
      list(APPEND fordArgs "--quiet")
    endif()
    if(FORD_NOSEARCH)
      list(APPEND fordArgs "--no-search")
    endif()
    if(FORD_PAGE_DIR)
      get_filename_component(page_dir_abs ${FORD_PAGE_DIR} ABSOLUTE)
      list(APPEND fordArgs --page_dir ${page_dir_abs})
    endif()
    if(FORD_CSS)
      get_filename_component(css_abs ${FORD_CSS} ABSOLUTE)
      list(APPEND fordArgs --css ${css_abs})
    endif()
    if(FORD_REVISION)
      list(APPEND fordArgs --revision ${FORD_REVISION})
    endif()
    if(FORD_SRC_DIR)
      list_to_absolute_paths(src_dirs ${CMAKE_CURRENT_SOURCE_DIR} ${FORD_SRC_DIR})
      list(append fordArgs --src_dir ${src_dirs})
    endif()
    if(FORD_EXCLUDE)
      list(append fordArgs --exclude ${FORD_EXCLUDE})
    endif()
    if(FORD_EXCLUDE_DIR)
      list_to_absolute_paths(excl_dirs ${CMAKE_CURRENT_SOURCE_DIR} ${FORD_EXCLUDE_DIR})
      list(append fordArgs --exclude_dir ${excl_dirs})
    endif()
    if(FORD_EXTENSIONS)
      list(append fordArgs --extensions ${FORD_EXTENSIONS})
    endif()
    if(FORD_MACROS)
      foreach(m ${FORD_MACROS})
        list(append fordArgs -m ${m})
      endforeach()
    endif()

    file(READ ${project_file_abs} projfile_contents)
    if(NOT GRAPHVIZ_PATH AND ${projfile_contents}
        MATCHES "(^|\n|\r)[ \t]*[Gg][Rr][Aa][Pp][Hh][ \t]*:[ \t]*[Tt][Rr][Uu][Ee][ \t]*($|\n|\r)")
      message(FATAL_ERROR "FORD documentation requies Graphviz, which is not installed! Aborting...")
    endif()

    if(NOT FORD_PAGE_DIR)
      if(${projfile_contents} MATCHES "(^|\n|\r)[ \t]*[Pp][Aa][Gg][Ee]_[Dd][Ii][Rr][ \t]*:[ \t]*([-A-Za-z0-9_]+)[ \t]*($|\n|\r)")
        get_filename_component(page_dir_abs ${CMAKE_MATCH_2} ABSOLUTE
    BASE_DIR ${projectfile_dir})
      endif()
    endif()

    if(page_dir_abs)
      file(GLOB_RECURSE page_files FOLLOW_SYMLINKS ${page_dir_abs}/*)
    endif()

    add_custom_target(${doc_name}
      ${FORD_PATH} ${fordArgs} -o ${CMAKE_CURRENT_BINARY_DIR}/${doc_name} ${project_file_abs}
      DEPENDS ${project_file_abs}
      SOURCES ${page_files}
      COMMENT "Building documentation with FORD. This may take some time."
    )

    add_custom_command(TARGET ${doc_name}
      POST_BUILD
      COMMAND ;
      COMMENT "Open ./${doc_name}/index.html in your browser to view documentation."
    )

  endfunction()
endif()