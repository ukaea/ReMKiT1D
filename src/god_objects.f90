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
module god_objects
    !! author:: Stefan Mijin
    !!
    !! Houses base Object used to determine whether complex objects have been properly initialized

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging

    implicit none
    private

    type ,public ,abstract :: Object
        !! Base object class used to determine whether complex objects have been properly initialized

        logical                 :: userDefined = .false. !! True only if user explicitly sets it to true

        contains

            procedure ,public :: isDefined   => isDefinedObject
            procedure ,public :: makeDefined => makeDefinedObject
            procedure ,public :: makeUndefined => makeUndefinedObject

    end type Object
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isDefinedObject(this) result(defined)
        !! Getter for userDefined

        class(Object) ,intent(in)  :: this 
        logical                    :: defined

    end function isDefinedObject
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine makeDefinedObject(this) 
        !! Set userDefined to .true.

        class(Object) ,intent(inout)  :: this 

    end subroutine makeDefinedObject
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine makeUndefinedObject(this) 
        !! Set userDefined to .false.

        class(Object) ,intent(inout)  :: this 

    end subroutine makeUndefinedObject
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module god_objects
!-----------------------------------------------------------------------------------------------------------------------------------
 