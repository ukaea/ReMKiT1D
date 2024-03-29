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
module runtime_constants
    !!Houses runtime constants

    use data_kinds, only: rk, ik
#ifdef DEBUG
    logical, parameter                  :: debugging = .true.
    logical, parameter                  :: assertions = .true.
#else
    logical, parameter                  :: debugging = .false.
    logical, parameter                  :: assertions = .false.
#endif

    integer(ik), parameter :: assertionLvl = 0 !! Sets assertion level. 0 is just the lowest cost assertions. When assertions = .true. it overrides this and all assertions are on
  
 !----------------------------------------------------------------------------------------------------------------------------------
 end module runtime_constants
 !----------------------------------------------------------------------------------------------------------------------------------
 