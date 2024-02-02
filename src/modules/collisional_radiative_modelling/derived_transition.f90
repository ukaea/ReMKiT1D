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
module derived_transition_class
    !! author: Stefan Mijin
    !!
    !! Houses fixed energy transition with a rate calculated using a derivation object

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions, assertionLvl
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use modelbound_data_abstract_class        ,only: ModelboundData
    use variable_container_class              ,only: VariableContainer
    use transition_abstract_class             ,only: Transition
    use derivation_abstract_class             ,only: Derivation

    implicit none
    private

    type ,public ,extends(Transition) :: DerivedTransition
        !! Transition where the rate is calculated using a derivation object and with a fixed transition energy.
        !! Optionally, momentum and energy rates can be calculated using derivation objects

        real(rk) ,private :: transitionEnergy 

        class(Derivation)         ,allocatable ,private :: rateDeriv !! Derivation object used in rate calculation
        integer(ik) ,allocatable ,dimension(:) ,private :: rateDerivIndices !! Indices for rate derivation

        class(Derivation)         ,allocatable ,private :: momentumRateDeriv !! Derivation object used in momentum rate calculation
        integer(ik) ,allocatable ,dimension(:) ,private :: momentumRateDerivIndices !! Indices for momentum rate derivation

        class(Derivation)         ,allocatable ,private :: energyRateDeriv !! Derivation object used in energy rate calculation
        integer(ik) ,allocatable ,dimension(:) ,private :: energyRateDerivIndices !! Indices for energy rate derivation

        integer(ik)                            ,private :: locNumX !! Size of local rate vectors 

        contains

        procedure ,public :: init => initDerivedTransition

        procedure ,public :: update => updateRates

        procedure ,public :: getTransitionEnergy => getEnergy

    end type DerivedTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initDerivedTransition(this,locNumX,inStates,outStates,energy,rateDeriv,rateDerivIndices,&
                                                momentumRateDeriv,momentumRateDerivIndices,energyRateDeriv,energyRateDerivIndices) 
        !! Initialization routine for DerivedTransition object

        class(DerivedTransition)            ,intent(inout)  :: this
        integer(ik)                         ,intent(in)     :: locNumX !! Local number of spatial cells
        integer(ik) ,dimension(:)           ,intent(in)     :: inStates !! Pre-transition states
        integer(ik) ,dimension(:)           ,intent(in)     :: outStates !! Post-transition states
        real(rk)                            ,intent(in)     :: energy !! Transition energy
        class(Derivation)                   ,intent(in)     :: rateDeriv !! Derivation object used in rate calculation
        integer(ik) ,dimension(:)           ,intent(in)     :: rateDerivIndices !! Indices for rate derivation
        class(Derivation) ,optional         ,intent(in)     :: momentumRateDeriv !! Derivation object used in momentum rate calculation
        integer(ik) ,optional ,dimension(:) ,intent(in)     :: momentumRateDerivIndices !! Indices for momentum rate derivation
        class(Derivation) ,optional         ,intent(in)     :: energyRateDeriv !! Derivation object used in energy rate calculation
        integer(ik) ,optional ,dimension(:) ,intent(in)     :: energyRateDerivIndices !! Indices for energy rate derivation

    end subroutine initDerivedTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateRates(this,varCont,hostModel,hostData,updatePriority)
        !! Update transition rates using derivation object

        class(DerivedTransition)        ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
        class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data - unused
        integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call - unused

    end subroutine updateRates 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getEnergy(this) result(energyCost)
        !! Returns array representing energy cost of this transition

        class(DerivedTransition)           ,intent(in) :: this 
        real(rk) ,allocatable ,dimension(:)            :: energyCost

    end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module derived_transition_class
!-----------------------------------------------------------------------------------------------------------------------------------
 