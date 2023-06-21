!-----------------------------------------------------------------------------------------------------------------------------------
! This file is part of ReMKiT1D.
!
! ReMKiT1D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! ReMKiT1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with ReMKiT1D. If not, see <https://www.gnu.org/licenses/>. 
!
! Copyright 2023 United Kingdom Atomic Energy Authority (stefan.mijin@ukaea.uk)
!-----------------------------------------------------------------------------------------------------------------------------------
module physical_constants
    !! Houses some useful physical/mathematical constants

    use data_kinds  ,only: rk, ik

    real(rk), parameter :: pi         = 4.d0 * atan(1.d0)
    real(rk), parameter :: epsilon0   = 8.854188d-12            !! Vacuum permittivity 
    real(rk), parameter :: elMass     = 9.10938d-31             !! Electron mass 
    real(rk), parameter :: elCharge   = 1.60218d-19             !! Electron charge 
    real(rk), parameter :: bohrRadius = 5.291772d-11 
    real(rk), parameter :: hPlanck    = 6.62607004d-34 
    real(rk), parameter :: mu0        = 4.00D00 * pi * 1d-17    !! Vacuum permeability
    real(rk), parameter :: protonMass = 1.67262d-27
    real(rk), parameter :: amu        = 1.6605390666d-27        !! Atomic mass unit
 
!----------------------------------------------------------------------------------------------------------------------------------
end module physical_constants
!----------------------------------------------------------------------------------------------------------------------------------
 