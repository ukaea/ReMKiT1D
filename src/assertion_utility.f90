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
module assertion_utility
    !! author: Stefan Mijin
    !! 
    !! Contains a toggleable assertion utility - based on Figure 10.8 of Scientific Software Design by Rouson et al. 

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging
    use iso_fortran_env             ,only: error_unit

    implicit none
    public

!-----------------------------------------------------------------------------------------------------------------------------------
    interface assert

        module subroutine assertSingle(assertion,error)

        logical                ,intent(in) :: assertion
        character(*) ,optional ,intent(in) :: error

        end subroutine assertSingle

        module subroutine assertVector(assertion,error)
    
            logical  ,dimension(:) ,intent(in) :: assertion
            character(*) ,optional ,intent(in) :: error

        end subroutine assertVector

    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    interface assertPure

        pure module subroutine assertPureSingle(assertion,error)

        logical                ,intent(in) :: assertion
        character(*) ,optional ,intent(in) :: error

        end subroutine assertPureSingle

        pure module subroutine assertPureVector(assertion,error)
    
            logical  ,dimension(:) ,intent(in) :: assertion
            character(*) ,optional ,intent(in) :: error

        end subroutine assertPureVector

    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    contains
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine assertSingle(assertion,error)
        !! Checks whether assertion is .false. and stops the program excecution with error text

        implicit none 

        logical                ,intent(in) :: assertion
        character(*) ,optional ,intent(in) :: error

        if (.not. assertion) then 

            write(error_unit,*) "Assertion failed with message: "
            if (present(error)) then
                write(error_unit,*) error
            else
                write(error_unit,*) "(no message provided)."
            end if
            error stop "Execution halted on failed assertion(s)!"

        end if
    end subroutine assertSingle
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine assertVector(assertion,error)
        !!  Checks whether any assertion is .false. and stops the program excecution with error text

        implicit none 

        logical  ,dimension(:) ,intent(in) :: assertion
        character(*) ,optional ,intent(in) :: error

        if (any(.not. assertion)) then 

            write(error_unit,*) "Assertion failed with message: "
            if (present(error)) then
                write(error_unit,*) error
            else
                write(error_unit,*) "(no message provided)."
            end if
            error stop "Execution halted on failed assertion(s)!"

        end if
    end subroutine assertVector
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine assertIdentical(integers, error)
        !! Checks all integers are identical and if not stops the program excecution with error text

        implicit none 

        integer ,dimension(:)  ,intent(in) :: integers
        character(*) ,optional ,intent(in) :: error

        if (any(integers /= integers(1))) then 

            write(error_unit,*) "Assertion failed with message: "
            if (present(error)) then
                write(error_unit,*) error
            else
                write(error_unit,*) "(no message provided)."
            end if

            error stop "Execution halted on failed assertion!"
        end if

    end subroutine assertIdentical
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine assertPureVector(assertion,error)
        !! Pure version of assertVector routine, does not write separately to std out

        implicit none 

        logical  ,dimension(:) ,intent(in) :: assertion
        character(*) ,optional ,intent(in) :: error
        character(:) ,allocatable          :: outputError

        if (any(.not. assertion)) then 
            outputError = "Assertion failed with message: "
            if (present(error)) then
                outputError = outputError // error 
                error stop outputError
            else
                error stop "Assertion failed with message: (no message provided)."
            end if

        end if
    end subroutine assertPureVector
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine assertPureSingle(assertion,error)
        !! Pure version of assertSingle routine, does not write separately to std out

        implicit none 

        logical                ,intent(in) :: assertion
        character(*) ,optional ,intent(in) :: error
        character(:) ,allocatable          :: outputError

        if (.not. assertion) then 
            outputError = "Assertion failed with message: "
            if (present(error)) then
                outputError = outputError // error 
                error stop outputError
            else
                error stop "Assertion failed with message: (no message provided)."
            end if

        end if
    end subroutine assertPureSingle
!-----------------------------------------------------------------------------------------------------------------------------------
 end module assertion_utility
!-----------------------------------------------------------------------------------------------------------------------------------
 