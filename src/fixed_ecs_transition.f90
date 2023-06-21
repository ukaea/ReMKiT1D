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
module fixed_ecs_transition_class
    !! author: Stefan Mijin
    !!
    !! Houses fixed transition energy and fixed cross-section transition object

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use modelbound_data_abstract_class        ,only: ModelboundData
    use variable_container_class              ,only: VariableContainer
    use transition_abstract_class             ,only: Transition
    use v_space_class                         ,only: VSpace
    use modelbound_CRM_data_class             ,only: ModelboundCRMData

    implicit none
    private

    type ,public ,extends(Transition) :: FixedECSTransition
        !! Transition object with a fixed transition energy and cross-section, which takes moments of a distribution function to
        !! calculate rates

        real(rk)              ,private :: transitionEnergy 
        type(VSpace) ,pointer ,private :: vSpaceRef !! Reference to velocity space object used in moment-taking

        logical               ,private :: takeMomentumMoment !! True if the momentum rate should be allocated and updated by taking 
                                                             !! the first moment of the l = 1 harmonic
        integer(ik)           ,private :: l1HarmonicInd !! Index of the l=1 m = 0 harmonic on the harmonic grid 
        integer(ik)           ,private :: distFunVarIndex !! Variable index of the distribution function to be used in moment-taking
        integer(ik)           ,private :: locNumX !! Size of local rate vectors 

        integer(ik)           ,private :: fixedWIndex !! Index of the inelastic weigth array corresponding to this transition
        
        contains

        procedure ,public :: init => initFixedECSTransition

        procedure ,public :: update => updateRates

        procedure ,public :: getTransitionEnergy => getEnergy

        final             :: finalizeFixedECSTransition

    end type FixedECSTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initFixedECSTransition(this,locNumX,inStates,outStates,energy,crossSection,distVarIndex,refVSpace, &
                                                fixedWIndex,momentumMoment, l1Index) 
        !! Initialization routine for FixedECSTransition object

        class(FixedECSTransition)          ,intent(inout)  :: this
        integer(ik)                        ,intent(in)     :: locNumX !! Local number of spatial cells
        integer(ik)  ,dimension(:)         ,intent(in)     :: inStates !! Pre-transition states
        integer(ik)  ,dimension(:)         ,intent(in)     :: outStates !! Post-transition states
        real(rk)                           ,intent(in)     :: energy !! Transition energy
        real(rk) ,dimension(:,:)           ,intent(in)     :: crossSection !! Cross-section harmonics
        integer(ik)                        ,intent(in)     :: distVarIndex !! Distribution function variable index
        type(VSpace) ,target               ,intent(inout)  :: refVSpace !! Target for the reference pointer
        integer(ik)                        ,intent(in)     :: fixedWIndex !! Index of the inelastic weight array corresponding to this transition
        logical ,optional                  ,intent(in)     :: momentumMoment !! Set to true if the momentum rate should be calculated
        integer(ik) ,optional              ,intent(in)     :: l1Index !! Index of the l=1 harmonic - must be provided if calculating momentum rate

    end subroutine initFixedECSTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getEnergy(this) result(energyCost)
        !! Returns array representing energy cost of this transition

        class(FixedECSTransition)          ,intent(in) :: this 
        real(rk) ,allocatable ,dimension(:)            :: energyCost

    end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateRates(this,varCont,hostModel,hostData,updatePriority)
        !! Update transition and moment exchange rate (if applicable)

        class(FixedECSTransition)       ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
        class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data used to determine emitting cells for discrete cross-section
        integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call - unused

    end subroutine updateRates 
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeFixedECSTransition(this) 
        !! Deallocate pointer component

        type(FixedECSTransition)                    ,intent(inout) :: this

    end subroutine finalizeFixedECSTransition 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module fixed_ecs_transition_class
!-----------------------------------------------------------------------------------------------------------------------------------
 