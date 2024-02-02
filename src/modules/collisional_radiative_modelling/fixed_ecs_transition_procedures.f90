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
submodule (fixed_ecs_transition_class) fixed_ecs_transition_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the fixed energy and cross section transition class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
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
    integer(ik)                        ,intent(in)     :: fixedWIndex !! Index of the inelastic weigth array corresponding to this transition
    logical ,optional                  ,intent(in)     :: momentumMoment !! Set to true if the momentum rate should be calculated
    integer(ik) ,optional              ,intent(in)     :: l1Index !! Index of the l=1 harmonic - must be provided if calculating momentum rate

    real(rk) ,allocatable ,dimension(:) :: rateVec

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(refVSpace%isDefined(),&
        "Undefined target VSpace object passed to fixed energy/cross-section transition constructor")
        if (present(momentumMoment)) then 
            if (momentumMoment) call assertPure(present(l1Index),"When momentumMoment is set to true in FixedECSTransition &
            &cosntructor the l1Index must be provided")
        end if
    end if

    this%locNumX = locNumX 
    allocate(rateVec(locNumX))
    rateVec = 0 
    call this%setRate(rateVec)
    this%takeMomentumMoment = .false.
    this%distFunVarIndex = distVarIndex
    this%transitionEnergy = energy 
    if (present(momentumMoment)) then 
        if (momentumMoment) call this%setRateMomentum(rateVec)
        this%takeMomentumMoment = .true. 
        this%l1HarmonicInd = l1Index
    end if 

    call this%setStates(inStates,outStates)

    this%vSpaceRef => refVSpace 
    this%fixedWIndex = fixedWIndex

    call this%setCrossSection(crossSection)
    call this%setCSDim(1)

    call this%setIncludeElectronDensity(.true.)

    call this%makeDefined()


end subroutine initFixedECSTransition
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getEnergy(this) result(energyCost)
    !! Returns array representing energy cost of this transition 

    class(FixedECSTransition)            ,intent(in) :: this 
    real(rk) ,allocatable ,dimension(:)            :: energyCost

    if (assertions) call assertPure(this%isDefined(),"getEnergy called on undefined FixedECSTransition object")

    allocate(energyCost(1))
    energyCost = this%transitionEnergy

end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateRates(this,varCont,hostModel,hostData,updatePriority)
    !! Update transition and moment exchange rate (if applicable)

    class(FixedECSTransition)       ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
    class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data used to determine emitting cells for discrete cross-section
    integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call - unused

    real(rk) ,allocatable ,dimension(:) :: rateBuffer , emissionVec
    integer(ik) :: inferredHaloWidth

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update undefined FixedECSTransition object")
        call assert(varCont%isDefined(),"Attempted to update FixedECSTransition object using undefined variable container")
    end if

    select type (hostData)
        type is (ModelboundCRMData)
            emissionVec = hostData%getFixedEmissionVector(this%fixedWIndex)
        class default 
            error stop "incompatible modelbound data passed to FixedECSTransition update routine"
        end select

    rateBuffer = this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,1,1,&
                                                this%getCrossSectionCol(1)*emissionVec)
    
    inferredHaloWidth = (size(rateBuffer) - this%locNumX)/2
    call this%setRate(rateBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth)) !Remove halo before saving
    call this%setRateEnergy(rateBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth)*this%transitionEnergy) 

    if (this%takeMomentumMoment) then

        rateBuffer = this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,&
                                                    this%l1HarmonicInd,2,this%getCrossSectionCol(2)*emissionVec)

        call this%setRateMomentum(rateBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth)) !Remove halo before saving
        
    end if

end subroutine updateRates 
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeFixedECSTransition(this) 
    !! Deallocate pointer component

    type(FixedECSTransition)                    ,intent(inout) :: this

    nullify(this%vSpaceRef)

end subroutine finalizeFixedECSTransition 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule fixed_ecs_transition_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
