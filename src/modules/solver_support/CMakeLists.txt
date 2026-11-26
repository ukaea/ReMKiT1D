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

remkit1d_add_module(solver_support
    PREPROCESS
    NO_STRICT_FLAGS
    LINK_LIBRARIES
        PkgConfig::PETSC
        parameters
        basic_support
        grid
        variables
        mpi_support
)

# PETSc >= 3.23 expands CHKERRQ(ierr) with __FILE__, which can exceed Fortran's
# default 132-column free-form line limit on long source paths.
SET(_solver_support_extra_flags "")
SET_COMPILE_FLAG(_solver_support_extra_flags Fortran "-ffree-line-length-none")
IF(_solver_support_extra_flags)
    TARGET_COMPILE_OPTIONS(solver_support PRIVATE ${_solver_support_extra_flags})
ENDIF()
