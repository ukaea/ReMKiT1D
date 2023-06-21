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
submodule (composite_manipulator_class) composite_manipulator_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the composite manipulator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initCompositeManipulator(this,numManipulators) 
    !! Compostite manipulator initialization routine

    class(CompositeManipulator)         ,intent(inout)  :: this
    integer(ik)                         ,intent(in)     :: numManipulators !! Number of manipulators expected

    allocate(this%manipulators(numManipulators))

    this%numManipulatorsAdded = 0
    this%allManipulatorsAdded = .false.

    allocate(this%manipulatorPriority(numManipulators))

    this%manipulatorPriority = 0

    call this%makeDefined()

end subroutine initCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addManipulator(this,manip,priority) 
    !! Add Manipulator object to composite manipulator

    class(CompositeManipulator)          ,intent(inout)  :: this
    class(Manipulator)                   ,intent(in)     :: manip
    integer(ik)                          ,intent(in)     :: priority

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add manipulator component to undefined composite manipulator")
        call assertPure(manip%isDefined(),"Attempted to add undefined manipulator component to composite manipulator")
        call assertPure(.not. this%allManipulatorsAdded,&
        "Attempted to add manipulator component to composite manipulator with no free manipulator slots")
    end if

    this%numManipulatorsAdded = this%numManipulatorsAdded + 1

    allocate(this%manipulators(this%numManipulatorsAdded)%entry,source=manip)

    this%manipulatorPriority(this%numManipulatorsAdded) = priority

    if (this%numManipulatorsAdded == size(this%manipulators)) this%allManipulatorsAdded = .true.

end subroutine addManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine manipulate(this,manipulatedModeller,outputVars,inputVars,priority) 
    !! Call affect routines of all manipulators whose manipulatorPriority <= priority

    class(CompositeManipulator)           ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller object used in callback
    class(VariableContainer)              ,intent(inout) :: outputVars !! Container for manipulation output
    class(VariableContainer)              ,intent(in)    :: inputVars !! Manipulation input variables
    integer(ik)                           ,intent(in)    :: priority !! Priority for this call

    integer(ik) :: i 

    if (assertions) then 

        call assertPure(this%isDefined(),"Attempted to manipulate modeller using undefined composite manipulator")
        call assertPure(manipulatedModeller%isDefined(),"Attempted to manipulate undefined modeller using composite manipulator")
        call assertPure(outputVars%isDefined(),"Undefined outputVars passed to manipulate all routine of composite manipulator")
        call assertPure(inputVars%isDefined(),"Undefined inputVars passed to manipulate all routine of composite manipulator")
        call assertPure(this%allManipulatorsAdded,"Attempted to manipulate modeller when not all manipulators have been added to &
        &composite manipulator object")

    end if

    do i = 1,size(this%manipulators)
        if (this%manipulatorPriority(i) <= priority) &
        call this%manipulators(i)%entry%affect(manipulatedModeller,outputVars,inputVars)
    end do

end subroutine manipulate
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule composite_manipulator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
