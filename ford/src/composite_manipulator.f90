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
module composite_manipulator_class
    !! Houses container/controller for multiple manipulator objects

    use data_kinds                           ,only: rk, ik
    use runtime_constants                    ,only: debugging, assertions
    use god_objects                          ,only: object
    use assertion_utility                    ,only: assert, assertIdentical, assertPure
    use variable_container_class             ,only: VariableContainer
    use modeller_surrogate_class             ,only: ModellerSurrogate
    use manipulator_abstract_class           ,only: Manipulator ,ManipulatorContainer

    implicit none
    private

    type ,public ,extends(Object) :: CompositeManipulator
        !! Composite manipulator object allowing for application of multiple manipulators in series based on their priority

        type(ManipulatorContainer) ,allocatable ,dimension(:) ,private :: manipulators !! Manipulators contained in this composite manipulator

        integer(ik)                ,allocatable ,dimension(:) ,private :: manipulatorPriority !! Manipulator priority for each manipulator object

        integer(ik)                                           ,private :: numManipulatorsAdded !! Counter keeping track of how many integrators have been added

        logical                                               ,private :: allManipulatorsAdded !! True when all manipulators have been allocated

        contains

        procedure ,public :: manipulate
        procedure ,public :: addManipulator

        procedure ,public :: init => initCompositeManipulator

    end type CompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initCompositeManipulator(this,numManipulators) 
        !! Compostite manipulator initialization routine

        class(CompositeManipulator)         ,intent(inout)  :: this
        integer(ik)                         ,intent(in)     :: numManipulators !! Number of manipulators expected

    end subroutine initCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addManipulator(this,manip,priority) 
        !! Add Manipulator object to composite manipulator

        class(CompositeManipulator)          ,intent(inout)  :: this
        class(Manipulator)                   ,intent(in)     :: manip
        integer(ik)                          ,intent(in)     :: priority

    end subroutine addManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine manipulate(this,manipulatedModeller,outputVars,inputVars,priority) 
        !! Call affect routines of all manipulators whose manipulatorPriority <= priority

        class(CompositeManipulator)           ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller object used in callback
        class(VariableContainer)              ,intent(inout) :: outputVars !! Container for manipulation output
        class(VariableContainer)              ,intent(in)    :: inputVars !! Manipulation input variables
        integer(ik)                           ,intent(in)    :: priority !! Priority for this call
        
    end subroutine manipulate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module composite_manipulator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 