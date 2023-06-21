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
module modeller_surrogate_class
    !! author: Stefan Mijin
    !!
    !! Contains empty Object extension type for use as a workaround for Fortran circular referencing - to be used for Modeller 
    !! objects
    
    use god_objects                 ,only: Object

    implicit none
    private 

    type, public, extends(Object) :: ModellerSurrogate
    end type ModellerSurrogate
!-----------------------------------------------------------------------------------------------------------------------------------
end module modeller_surrogate_class
!-----------------------------------------------------------------------------------------------------------------------------------

