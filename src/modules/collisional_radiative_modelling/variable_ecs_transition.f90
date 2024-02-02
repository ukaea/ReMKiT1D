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
module variable_ecs_transition_class
    !! author: Stefan Mijin
    !!
    !! Houses variable energy and cross section transition

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions, assertionLvl
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use modelbound_data_abstract_class        ,only: ModelboundData
    use variable_container_class              ,only: VariableContainer
    use transition_abstract_class             ,only: Transition
    use derivation_abstract_class             ,only: Derivation ,DerivationContainer
    use v_space_class                         ,only: VSpace
    use modelbound_CRM_data_class             ,only: ModelboundCRMData
    use support_types

    implicit none
    private

    type ,public ,extends(Transition) :: VariableECSTransition
        !! Transition where the transition energy and cross-section are derived 

        real(rk) ,allocatable ,dimension(:) ,private :: energyCost
        real(rk) ,allocatable ,dimension(:) ,private :: flattenedEmissionVector

        class(Derivation)         ,allocatable ,private :: energyDeriv !! Derivation object used in transition energy calculation
        integer(ik) ,allocatable ,dimension(:) ,private :: energyDerivIndices !! Indices for transition energy derivation

        type(DerivationContainer) ,allocatable ,dimension(:) ,private :: csDerivs !! Derivation objects for the various cross section data columns
        type(IntArray)            ,allocatable ,dimension(:) ,private :: csDerivsReqIndices !! Required indices for the various cross section data column derivations
        integer(ik)               ,allocatable ,dimension(:) ,private :: csCols !! Cross section columns corresponding to each of the derivations

        integer(ik)                            ,private :: locNumX !! Size of local rate vectors 
        
        type(VSpace) ,pointer ,private :: vSpaceRef !! Reference to velocity space object used in moment-taking
        logical               ,private :: takeMomentumMoment !! True if the momentum rate should be allocated and updated by taking 
                                                             !! the first moment of the l = 1 harmonic
        integer(ik)           ,private :: l1HarmonicInd !! Index of the l=1 m = 0 harmonic on the harmonic grid 
        integer(ik)           ,private :: distFunVarIndex !! Variable index of the distribution function to be used in moment-taking

        contains

        procedure ,public :: init => initVariableECSTransition

        procedure ,public :: update => updateCSRates

        procedure ,public :: getTransitionEnergy => getEnergy

        final             :: finalizeVariableECSTransition

    end type VariableECSTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initVariableECSTransition(this,locNumX,inStates,outStates,energyDeriv,energyDerivIndices,csDerivs,&
                                                     csDerivsReqIndices,csCols,distVarIndex,refVSpace,momentumMoment,l1Index) 
        !! Initialization routine for VariableECSTransition object

        class(VariableECSTransition)            ,intent(inout)  :: this
        integer(ik)                             ,intent(in)     :: locNumX !! Local number of spatial cells
        integer(ik) ,dimension(:)               ,intent(in)     :: inStates !! Pre-transition states
        integer(ik) ,dimension(:)               ,intent(in)     :: outStates !! Post-transition states
        class(Derivation)                       ,intent(in)     :: energyDeriv !! Derivation object used in rate calculation
        integer(ik) ,dimension(:)               ,intent(in)     :: energyDerivIndices !! Indices for rate derivation
        type(DerivationContainer) ,dimension(:) ,intent(in)     :: csDerivs !! Derivation objects for the various cross section data columns
        type(IntArray)            ,dimension(:) ,intent(in)     :: csDerivsReqIndices !! Required indices for the various cross section data column derivations
        integer(ik) ,dimension(:)               ,intent(in)     :: csCols !! Cross section columns corresponding to each of the derivations
        integer(ik)                             ,intent(in)     :: distVarIndex !! Distribution function variable index
        type(VSpace) ,target                    ,intent(inout)  :: refVSpace !! Target for the reference pointer
        logical ,optional                       ,intent(in)     :: momentumMoment !! Set to true if the momentum rate should be calculated
        integer(ik) ,optional                   ,intent(in)     :: l1Index !! Index of the l=1 harmonic - must be provided if calculating momentum rate

    end subroutine initVariableECSTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateCSRates(this,varCont,hostModel,hostData,updatePriority)
        !! Update transition properties

        class(VariableECSTransition)    ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
        class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data - used to get emission vectors
        integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call - unused

    end subroutine updateCSRates 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getEnergy(this) result(energyCost)
        !! Returns array representing energy cost of this transition

        class(VariableECSTransition)       ,intent(in) :: this 
        real(rk) ,allocatable ,dimension(:)            :: energyCost

    end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeVariableECSTransition(this) 
        !! Deallocate pointer component

        type(VariableECSTransition)                    ,intent(inout) :: this

    end subroutine finalizeVariableECSTransition 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module variable_ecs_transition_class
!-----------------------------------------------------------------------------------------------------------------------------------
 