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
project(tests Fortran)

add_subdirectory(test_basic_support)
add_subdirectory(test_grid)
add_subdirectory(test_variables)
add_subdirectory(test_solver_support)
add_subdirectory(test_mpi)
add_subdirectory(test_terms)
add_subdirectory(test_model)
add_subdirectory(test_modeller)
add_subdirectory(test_hdf5)
add_subdirectory(test_json)
add_subdirectory(test_fluid_terms)
add_subdirectory(test_inel_mapping)
add_subdirectory(test_extended_support)
add_subdirectory(test_kinetic_terms)
add_subdirectory(test_fluid_models)
add_subdirectory(test_kinetic_models)
add_subdirectory(test_common_derivations)
add_subdirectory(test_timeloop)
