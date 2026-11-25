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

remkit1d_add_module(kinetic_models
    LINK_LIBRARIES
        parameters
        basic_support
        grid
        variables
        solver_support
        terms
        kinetic_terms
        modelbound_data
        abstract_model_builder
        model
        modeller
        IO_support
        extended_support
        common_derivations
        signals
        collisional_radiative_modelling
        general_terms
        normalization
        fluid_terms
)
