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
module simple_transition_class
    !! author: Stefan Mijin
    !!
    !! Houses simple fixed rate and energy transition object

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use modelbound_data_abstract_class        ,only: ModelboundData
    use variable_container_class              ,only: VariableContainer
    use transition_abstract_class             ,only: Transition

    implicit none
    private

    type ,public ,extends(Transition) :: SimpleTransition
        !! Simple transition object with a fixed rate and energy, and single ingoing and outgoing state. The default energy rate is
        !! set to rate*energy, and the momentum rate is 0

        real(rk) ,private :: transitionEnergy 

        contains

        procedure ,public :: init => initSimpleTransition

        procedure ,public :: getTransitionEnergy => getEnergy

    end type SimpleTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initSimpleTransition(this,locNumX,inState,outState,energy,fixedRate) 
        !! Initialization routine for SimpleTransition object

        class(SimpleTransition)            ,intent(inout)  :: this
        integer(ik)                        ,intent(in)     :: locNumX !! Local number of spatial cells
        integer(ik)                        ,intent(in)     :: inState !! Pre-transition state
        integer(ik)                        ,intent(in)     :: outState !! Post-transition state
        real(rk)                           ,intent(in)     :: energy !! Transition energy
        real(rk)                           ,intent(in)     :: fixedRate !! Fixed transition rate

    end subroutine initSimpleTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getEnergy(this) result(energyCost)
        !! Returns array representing energy cost of this transition

        class(SimpleTransition)            ,intent(in) :: this 
        real(rk) ,allocatable ,dimension(:)            :: energyCost

    end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module simple_transition_class
!-----------------------------------------------------------------------------------------------------------------------------------
 