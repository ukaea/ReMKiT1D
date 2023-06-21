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
module transition_abstract_class
    !! author: Stefan Mijin
    !!
    !! Houses abstract object responsible for handling rates and cross-sections, as well as characterization of various transitions

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use modelbound_data_abstract_class        ,only: ModelboundData
    use variable_container_class              ,only: VariableContainer

    implicit none
    private

    type ,public ,extends(Object), abstract :: Transition
        !! Represents a process in which an ingoing set of particles(states) is transformed into an outgoing set and to which 
        !! quantities representing particle/momentum/energy exchange rates can be meaningfully assigned. Optionally, a Transition can have
        !! a cross-section object of rank 2 to accommodate anisotropic differential cross-sections.

        integer(ik) ,allocatable ,dimension(:) ,private :: ingoingStates  !! Array of IDs associated with the ingoing states
        integer(ik) ,allocatable ,dimension(:) ,private :: outgoingStates !! Array of IDs associated with the outgoing states

        real(rk)    ,allocatable ,dimension(:) ,private :: rate         !! Spatial array associated with the (particle) Transition rate
        real(rk)    ,allocatable ,dimension(:) ,private :: rateMomentum !! Spatial array associated with the momentum exchange rate
        real(rk)    ,allocatable ,dimension(:) ,private :: rateEnergy   !! Spatial array associated with the energy rate

        real(rk)    ,allocatable ,dimension(:,:) ,private :: crossSection !! Differential cross-section associated with the process - optional

        logical                                 ,private :: rateContainsElDensity = .false. !! Set to true if rate calculations already include one electron 
                                                                                            !! density factor from taking the rate moment

        integer                                 ,private :: csDim = -1 !! Cross-section dimensionality - used to properly expose cs data 

        contains

        procedure ,public :: getIngoingStates
        procedure ,public :: getOutgoingStates
        procedure ,public :: setStates

        procedure ,public :: setRate 
        procedure ,public :: getRate 

        procedure ,public :: setRateMomentum 
        procedure ,public :: getRateMomentum

        procedure ,public :: setRateEnergy 
        procedure ,public :: getRateEnergy 

        procedure ,public :: setCrossSection
        procedure ,public :: setCrossSectionCol
        procedure ,public :: getCrossSectionCol

        procedure ,public :: includesElDensity
        procedure ,public :: setIncludeElectronDensity

        procedure ,public :: setCSDim
        procedure ,public :: getCSDim

        procedure ,public :: getRateSize

        procedure ,public :: update => noUpdate

        procedure(returnEnergy) ,deferred :: getTransitionEnergy

    end type Transition
!-----------------------------------------------------------------------------------------------------------------------------------
!Abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function returnEnergy(this) result(energyCost)
            !! Returns array representing energy cost of this transition

            import :: Transition ,rk

            class(Transition)                  ,intent(in) :: this 
            real(rk) ,allocatable ,dimension(:)            :: energyCost

        end function returnEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getIngoingStates(this) result(inStates)
        !! Getter for ingoingStates

        class(Transition)          ,intent(in)  :: this
        integer(ik) ,allocatable ,dimension(:)  :: inStates

    end function getIngoingStates  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getOutgoingStates(this) result(outStates)
        !! Getter for outgoingStates

        class(Transition)          ,intent(in)  :: this
        integer(ik) ,allocatable ,dimension(:)  :: outStates

    end function getOutgoingStates  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setStates(this,inStates,outStates)
        !! Setter for both ingoing and outgoing states

        class(Transition)          ,intent(inout)  :: this
        integer(ik)  ,dimension(:) ,intent(in)     :: inStates !! Ingoing state IDs
        integer(ik)  ,dimension(:) ,intent(in)     :: outStates !! Outgoing state IDs

    end subroutine setStates  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setRate(this,rate)
        !! Setter for rate values

        class(Transition)        ,intent(inout)  :: this
        real(rk)  ,dimension(:)  ,intent(in)     :: rate

    end subroutine setRate  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getRate(this) result(rate)
        !! Getter for rate values

        class(Transition)          ,intent(in)  :: this
        real(rk)    ,allocatable ,dimension(:)  :: rate

    end function getRate  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setRateMomentum(this,rate)
        !! Setter for rateMomentum values

        class(Transition)        ,intent(inout)  :: this
        real(rk)  ,dimension(:)  ,intent(in)     :: rate

    end subroutine setRateMomentum  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getRateMomentum(this) result(rate)
        !! Getter for rateMomentum values

        class(Transition)          ,intent(in)  :: this
        real(rk)    ,allocatable ,dimension(:)  :: rate

    end function getRateMomentum  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setRateEnergy(this,rate)
        !! Setter for rateEnergy

        class(Transition)        ,intent(inout)  :: this
        real(rk)  ,dimension(:)  ,intent(in)     :: rate

    end subroutine setRateEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getRateEnergy(this) result(rate)
        !! Getter for rateEnergy

        class(Transition)          ,intent(in)  :: this
        real(rk)    ,allocatable ,dimension(:)  :: rate

    end function getRateEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setCrossSectionCol(this,crossSection,col)
        !! Set cross-section values in column col

        class(Transition)        ,intent(inout)  :: this
        real(rk)  ,dimension(:)  ,intent(in)     :: crossSection
        integer(ik)              ,intent(in)     :: col

    end subroutine setCrossSectionCol  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setCrossSection(this,crossSection)
        !! Setter for crossSection

        class(Transition)         ,intent(inout)  :: this
        real(rk)  ,dimension(:,:) ,intent(in)     :: crossSection

    end subroutine setCrossSection  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCrossSectionCol(this,col) result(crossSection)
        !! Get cross-section values from column col

        class(Transition)          ,intent(in)  :: this
        integer(ik)              ,intent(in)    :: col
        real(rk)    ,allocatable ,dimension(:)  :: crossSection

    end function getCrossSectionCol  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function includesElDensity(this) result(includesDens)
        !! Check whether rates in this transition include an electron density factor

        class(Transition)        ,intent(in)  :: this
        logical                               :: includesDens

    end function includesElDensity  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setIncludeElectronDensity(this,includeDens)
        !! Setter for rateContainsElDensity

        class(Transition)    ,intent(inout)  :: this
        logical              ,intent(in)     :: includeDens

    end subroutine setIncludeElectronDensity  
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine noUpdate(this,varCont,hostModel,hostData,updatePriority)
        !! Default update routine - does nothing

        class(Transition)               ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Optional host model reference for callbacks during update
        class(ModelboundData) ,optional ,intent(in)     :: hostData !! Optional host data reference for callbacks during update
        integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call 

    end subroutine noUpdate 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCSDim(this) result(csDim)
        !! Getter for csDim

        class(Transition)        ,intent(in)  :: this
        integer(ik)                           :: csDim

    end function getCSDim  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setCSDim(this,csDim)
        !! Setter for csDim

        class(Transition)    ,intent(inout)  :: this
        integer(ik)          ,intent(in)     :: csDim

    end subroutine setCSDim   
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getRateSize(this) result(rateSize)
        !! Getter for rate array length

        class(Transition)        ,intent(in)  :: this
        integer(ik)                           :: rateSize

    end function getRateSize  
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module transition_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 