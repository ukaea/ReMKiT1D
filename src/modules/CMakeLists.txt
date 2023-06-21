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
project(modules Fortran)

add_subdirectory(parameters)
add_subdirectory(basic_support)
add_subdirectory(grid)
add_subdirectory(variables)
add_subdirectory(solver_support)
add_subdirectory(mpi_support)
add_subdirectory(terms)
add_subdirectory(modelbound_data)
add_subdirectory(model)
add_subdirectory(testing_support)
add_subdirectory(abstract_manipulators)
add_subdirectory(modeller)
add_subdirectory(integrator_abstractions)
add_subdirectory(basic_integrators)
add_subdirectory(IO_support)
add_subdirectory(abstract_model_builder)
add_subdirectory(fluid_terms)
add_subdirectory(normalization)
add_subdirectory(collisional_radiative_modelling)
add_subdirectory(kinetic_terms)
add_subdirectory(extended_support)
add_subdirectory(fluid_models)
add_subdirectory(common_derivations)
add_subdirectory(modeller_assembly)
add_subdirectory(timeloop)
add_subdirectory(signals)
add_subdirectory(common_manipulators)
add_subdirectory(general_terms)
add_subdirectory(general_models)
add_subdirectory(kinetic_models)
add_subdirectory(extrapolation_support)


