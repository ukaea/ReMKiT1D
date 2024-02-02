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
submodule (derived_transition_class) derived_transition_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the fixed energy derived rate transition class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
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

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(rateDeriv%isDefined(),&
        "Derivation object passed to DerivedTransition constructor is not defined")

        if (present(momentumRateDeriv)) then
            call assertPure(present(momentumRateDerivIndices),"If momentumRateDeriv is passed to DerivedTransition constructor &
            &so must be momentumRateDerivIndices")
            call assertPure(momentumRateDeriv%isDefined(),"Undefined momentum rate derivation object passed to DerivedTransition&
            & constructor")
        end if

        if (present(energyRateDeriv)) then
            call assertPure(present(energyRateDerivIndices),"If energyRateDeriv is passed to DerivedTransition constructor &
            &so must be energyRateDerivIndices")
            call assertPure(energyRateDeriv%isDefined(),"Undefined momentum rate derivation object passed to DerivedTransition&
            & constructor")
        end if

    end if

    call this%setStates(inStates,outStates)

    this%transitionEnergy = energy
    this%locNumX = locNumX

    allocate(this%rateDeriv,source=rateDeriv)
    this%rateDerivIndices = rateDerivIndices 

    if (present(momentumRateDeriv)) then
        allocate(this%momentumRateDeriv,source=momentumRateDeriv)
        this%momentumRateDerivIndices = momentumRateDerivIndices 
    end if

    if (present(energyRateDeriv)) then
        allocate(this%energyRateDeriv,source=energyRateDeriv)
        this%energyRateDerivIndices = energyRateDerivIndices 
    end if

    call this%makeDefined()

end subroutine initDerivedTransition
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getEnergy(this) result(energyCost)
    !! Returns array representing energy cost of this transition 

    class(DerivedTransition)           ,intent(in) :: this 
    real(rk) ,allocatable ,dimension(:)            :: energyCost

    if (assertions) call assertPure(this%isDefined(),"getEnergy called on undefined DerivedTransition object")

    if (allocated(this%energyRateDeriv)) then 
        energyCost = this%getRateEnergy()/this%getRate()
    else
        allocate(energyCost(1))
        energyCost = this%transitionEnergy    
    end if

end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateRates(this,varCont,hostModel,hostData,updatePriority)
    !! Update transition rates using derivation object

    class(DerivedTransition)        ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
    class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data - unused
    integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call - unused

    real(rk) ,allocatable ,dimension(:)   :: rateBuffer

    integer(ik) :: inferredHaloWidth

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update undefined DerivedTransition object")
        call assert(varCont%isDefined(),"Attempted to update DerivedTransition object using undefined variable container")
    end if

    rateBuffer = this%rateDeriv%calculate(varCont%variables,this%rateDerivIndices)
    inferredHaloWidth = (size(rateBuffer) - this%locNumX)/2
    call this%setRate(rateBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth)) !Remove halo before saving

    if (allocated(this%energyRateDeriv)) then
        rateBuffer = this%energyRateDeriv%calculate(varCont%variables,this%energyRateDerivIndices)
        inferredHaloWidth = (size(rateBuffer) - this%locNumX)/2
        call this%setRateEnergy(rateBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth)) !Remove halo before saving
    else
        call this%setRateEnergy(rateBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth)*this%transitionEnergy)
    end if

    if (allocated(this%momentumRateDeriv)) then
        rateBuffer = this%momentumRateDeriv%calculate(varCont%variables,this%momentumRateDerivIndices)
        inferredHaloWidth = (size(rateBuffer) - this%locNumX)/2
        call this%setRateMomentum(rateBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth)) !Remove halo before saving
    end if
    
end subroutine updateRates 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule derived_transition_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
