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
submodule (god_objects) god_objects_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains basic definition routines for base Object class

    implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
    contains
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isDefinedObject(this) result(defined)
        !! Getter for userDefined

        class(Object) ,intent(in)  :: this 
        logical                    :: defined

        defined = this%userDefined

    end function isDefinedObject
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine makeDefinedObject(this) 
        !! Set userDefined to .true.

        class(Object) ,intent(inout)  :: this 

        this%userDefined = .true.

    end subroutine makeDefinedObject
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine makeUndefinedObject(this) 
        !! Set userDefined to .false.

        class(Object) ,intent(inout)  :: this 

        this%userDefined = .false.

    end subroutine makeUndefinedObject
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule god_objects_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
