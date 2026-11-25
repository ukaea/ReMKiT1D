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

remkit1d_add_module(basic_integrators
    NO_STRICT_FLAGS
    LINK_LIBRARIES
        parameters
        basic_support
        variables
        abstract_manipulators
        integrator_abstractions
        modeller
        mpi_support
        IO_support
        MPI::MPI_Fortran
        SUNDIALS::cvode
        SUNDIALS::fcvode_mod
        SUNDIALS::core
        SUNDIALS::fcore_mod
        SUNDIALS::fnvecparallel_mod
        SUNDIALS::nvecparallel
        SUNDIALS::nvecserial
        SUNDIALS::fnvecserial_mod
        SUNDIALS::fsunlinsolspgmr_mod
)

# SUNDIALS F2003 matrix-free Krylov solvers pass an unassociated SUNMatrix pointer
# to FCVodeSetLinearSolver, which triggers gfortran's -fcheck=pointer runtime check.
SET_SOURCE_FILES_PROPERTIES(cvode_integrator_procedures.f90 PROPERTIES
    COMPILE_OPTIONS "$<$<AND:$<CONFIG:DEBUG>,$<Fortran_COMPILER_ID:GNU>>:-fcheck=no-pointer>"
)
