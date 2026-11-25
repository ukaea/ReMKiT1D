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

remkit1d_add_module(extended_support
    LINK_LIBRARIES
        parameters
        basic_support
        grid
        variables
        IO_support
        mpi_support
        solver_support
        normalization
        common_derivations
        basic_integrators
        integrator_abstractions
        signals
        modelbound_data
        collisional_radiative_modelling
        model
        common_manipulators
        abstract_manipulators
        extrapolation_support
        kinetic_terms
)
