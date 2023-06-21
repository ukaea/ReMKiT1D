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
module db_transition_class
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
    use sparse_row_data_class                 ,only: SparseRowData
    use derivation_abstract_class             ,only: Derivation
    use modelbound_CRM_data_class             ,only: ModelboundCRMData
    use support_types                 

    implicit none
    private

    type ,public ,extends(Transition) :: DBTransition
        !! Transition with fixed energy where the cross-section is calculated using SOL-KiT's version of discrete detailed balance

        real(rk)              ,private :: transitionEnergy 
        type(VSpace) ,pointer ,private :: vSpaceRef !! Reference to velocity space object used in moment-taking

        logical               ,private :: takeMomentumMoment !! True if the momentum rate should be allocated and updated by taking 
                                                             !! the first moment of the l = 1 harmonic

        logical ,private :: strictDB !! Set to false if strict detailed balance should not be enforced by scaling cross-sections. Defaults to true.

        integer(ik)           ,private :: l1HarmonicInd !! Index of the l=1 m = 0 harmonic on the harmonic grid 
        integer(ik)           ,private :: distFunVarIndex !! Variable index of the distribution function to be used in moment-taking
        integer(ik)           ,private :: locNumX !! Size of local rate vectors 
        integer(ik)           ,private :: directTransitionIndex !! Index of the direct transition in the host model data
        integer(ik)           ,private :: fixedWIndexDirect !! Index of the direct transition inelastic weight matrix in the host model inelastic data object
        integer(ik)           ,private :: temperatureVarIndex !! Index of the temperature variable used to calculate the detailed balance cross-section
        integer(ik)           ,private :: fixedWIndex !! Index of this transition's inelastic weight matrix in the host model inelastic data object
        
        integer(ik)           ,private :: maxCSl !! Highest harmonic of the cross-section to calculate

        real(rk)              ,private :: degeneracyRatio !! Ratio of the degeneracy of the initial and final states of the transition

        class(Derivation) ,allocatable ,private :: degeneracyFun !! Optional derivation object when the degeneracy is a function of variables in the variable container
        integer(ik) ,allocatable ,dimension(:) ,private :: degeneracyFunIndices !! Variable indices needed for the degeneracy function calculation 
        
        integer(ik) ,private :: csUpdatePriority !! Update priority for cross-section data

        type(SparseRowData) ,allocatable ,private :: directFixedW !! Local copy to avoid expensive copy routines
        type(RealArray) ,allocatable ,dimension(:) ,private :: discreteEnergyErrors !! Difference between discrete and analytic transition energies for each row/column of directFixedW
        contains

        procedure ,public :: init => initDBTransition

        procedure ,public :: update => updateCSRates

        procedure ,public :: getTransitionEnergy => getEnergy

        final             :: finalizeDBTransition

    end type DBTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initDBTransition(this,locNumX,inStates,outStates,energy,distVarIndex,refVSpace, &
        directTransitionIndex, fixedWIndexDirect,fixedWIndex,&
        temperatureVarIndex, maxCSl,degeneracyRatio,&
        degeneracyFun,degeneracyFunIndices, momentumMoment ,l1Index ,csUpdatePriority,strictDB) 
        !! Initialization routine for DBTransition object

        class(DBTransition)                 ,intent(inout)  :: this
        integer(ik)                         ,intent(in)     :: locNumX !! Local number of spatial cells
        integer(ik)  ,dimension(:)          ,intent(in)     :: inStates !! Pre-transition states
        integer(ik)  ,dimension(:)          ,intent(in)     :: outStates !! Post-transition states
        real(rk)                            ,intent(in)     :: energy !! Transition energy
        integer(ik)                         ,intent(in)     :: distVarIndex !! Distribution function variable index
        type(VSpace) ,target                ,intent(inout)  :: refVSpace !! Target for the reference pointer
        integer(ik)                         ,intent(in)     :: directTransitionIndex !! Index of the direct transition in the host model data
        integer(ik)                         ,intent(in)     :: fixedWIndexDirect !! Index of the direct transition inelastic weight matrix in the host model inelastic data object
        integer(ik)                         ,intent(in)     :: fixedWIndex !! Index of this transition's inelastic weight matrix in the host model inelastic data object
        integer(ik)                         ,intent(in)     :: temperatureVarIndex !! Index of the temperature variable used to calculate the detailed balance cross-section
        integer(ik)                         ,intent(in)     :: maxCSl !! Highest harmonic of the cross-section to calculate
        real(rk)                            ,intent(in)     :: degeneracyRatio !! Ratio of the degeneracy of the initial and final states of the transition
        class(Derivation) ,optional         ,intent(in)     :: degeneracyFun !! Optional derivation object when the degeneracy is a function of variables in the variable container
        integer(ik) ,dimension(:) ,optional ,intent(in)     :: degeneracyFunIndices !! Variable indices needed for the degeneracy function calculation 
        logical ,optional                   ,intent(in)     :: momentumMoment !! Set to true if the momentum rate should be calculated
        integer(ik) ,optional               ,intent(in)     :: l1Index !! Index of the l=1 harmonic - must be provided if calculating momentum rate
        integer(ik) ,optional               ,intent(in)     :: csUpdatePriority !! Update priority for cross-section data. Defaults to highest priority (0)
        logical ,optional                   ,intent(in)     :: strictDB !! Set to false if strict detailed balance should not be enforced by scaling cross-sections. Defaults to true.

    end subroutine initDBTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getEnergy(this) result(energyCost)
        !! Returns array representing energy cost of this transition

        class(DBTransition)                ,intent(in) :: this 
        real(rk) ,allocatable ,dimension(:)            :: energyCost

    end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateCSRates(this,varCont,hostModel,hostData,updatePriority)
        !! Update cross-section and transition and moment exchange rates (if applicable)
    
        class(DBTransition)             ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
        class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data used to access direct transition data 
        integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call used to determine if cross-sections should be updated

    end subroutine updateCSRates 
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeDBTransition(this) 
        !! Deallocate pointer component

        type(DBTransition)                    ,intent(inout) :: this

    end subroutine finalizeDBTransition 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module db_transition_class
!-----------------------------------------------------------------------------------------------------------------------------------
 