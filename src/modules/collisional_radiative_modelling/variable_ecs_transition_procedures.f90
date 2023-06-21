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
submodule (variable_ecs_transition_class) variable_ecs_transition_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the variable energy and cross section transition class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
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

    real(rk) ,allocatable ,dimension(:) :: rateVec
    real(rk) ,allocatable ,dimension(:,:) :: crossSection

    if (assertions) then 
        call assertPure(refVSpace%isDefined(),&
        "Undefined target VSpace object passed to variable energy/cross-section transition constructor")
        if (present(momentumMoment)) then 
            if (momentumMoment) call assertPure(present(l1Index),"When momentumMoment is set to true in VariableECSTransition &
            &constructor the l1Index must be provided")
        end if

        call assertPure(size(csDerivs) == size(csDerivsReqIndices),&
        "csDerivs and csDerivsReqIndices passed to initVariableECSTransition must be of the same size")
        call assertPure(size(csDerivs) == size(csCols),&
        "csDerivs and csCols passed to initVariableECSTransition must be of the same size")
    end if

    this%locNumX = locNumX 
    allocate(rateVec(locNumX))
    rateVec = 0 
    call this%setRate(rateVec)
    this%takeMomentumMoment = .false.
    this%distFunVarIndex = distVarIndex
    if (present(momentumMoment)) then 
        if (momentumMoment) call this%setRateMomentum(rateVec)
        this%takeMomentumMoment = .true. 
        this%l1HarmonicInd = l1Index
    end if 

    call this%setStates(inStates,outStates)

    this%vSpaceRef => refVSpace 

    allocate(crossSection(refVSpace%getNumV()*this%locNumX, maxval(csCols)))
    crossSection = 0
    call this%setCrossSection(crossSection)
    call this%setCSDim(2)

    call this%setIncludeElectronDensity(.true.)

    allocate(this%energyDeriv,source=energyDeriv)
    this%energyDerivIndices = energyDerivIndices

    allocate(this%csDerivs,source=csDerivs)
    this%csDerivsReqIndices = csDerivsReqIndices
    this%csCols = csCols

    allocate(this%energyCost(this%locNumX))
    this%energyCost = 0

    allocate(this%flattenedEmissionVector(this%locNumX*refVSpace%getNumV()))
    this%flattenedEmissionVector = 0

    call this%makeDefined()

end subroutine initVariableECSTransition
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateCSRates(this,varCont,hostModel,hostData,updatePriority)
    !! Update transition properties

    class(VariableECSTransition)    ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
    class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data - used to get emission vectors
    integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call - unused

    real(rk) ,allocatable ,dimension(:)   :: calcBuffer

    integer(ik) :: i, inferredHaloWidth ,numV

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update undefined VariableECSTransition object")
        call assert(varCont%isDefined(),"Attempted to update VariableECSTransition object using undefined variable container")
    end if

    calcBuffer = this%energyDeriv%calculate(varCont%variables,this%energyDerivIndices)
    inferredHaloWidth = (size(calcBuffer) - this%locNumX)/2
    this%energyCost = calcBuffer(1+inferredHaloWidth:this%locNumX+inferredHaloWidth) !Remove halo before saving

    numV = this%vSpaceRef%getNumV()
    select type (hostData)
    type is (ModelboundCRMData)
        do i = 1,this%locNumX
            this%flattenedEmissionVector((i-1)*numV+1:i*numV) = hostData%getInterpolatedEmissionVector(this%energyCost(i))
        end do
    class default 
        error stop "incompatible modelbound data passed to VariableECSTransition update routine"
    end select

    deallocate(calcBuffer)

    do i = 1,size(this%csCols)
        calcBuffer = this%csDerivs(i)%entry%calculate(varCont%variables,this%csDerivsReqIndices(i)%entry)
        inferredHaloWidth = (size(calcBuffer) - this%locNumX*numV)/2
        call this%setCrossSectionCol(calcBuffer(1+inferredHaloWidth:this%locNumX*numV+inferredHaloWidth),this%csCols(i))
    end do

    deallocate(calcBuffer)

    calcBuffer = this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,1,1,&
                                                g=this%getCrossSectionCol(1)*this%flattenedEmissionVector,gDependsOnX=.true.)
    
    call this%setRate(calcBuffer) 
    call this%setRateEnergy(calcBuffer*this%energyCost) 

    if (this%takeMomentumMoment) then

        calcBuffer = this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,&
                            this%l1HarmonicInd,2,this%getCrossSectionCol(2)*this%flattenedEmissionVector,gDependsOnX=.true.)

        call this%setRateMomentum(calcBuffer) 
        
    end if

end subroutine updateCSRates 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getEnergy(this) result(energyCost)
    !! Returns array representing energy cost of this transition

    class(VariableECSTransition)       ,intent(in) :: this 
    real(rk) ,allocatable ,dimension(:)            :: energyCost

    if (assertions) call assertPure(this%isDefined(),"getEnergy called on undefined VariableECSTransition object")

    energyCost = this%energyCost

end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeVariableECSTransition(this) 
    !! Deallocate pointer component

    type(VariableECSTransition)                    ,intent(inout) :: this

    nullify(this%vSpaceRef)

end subroutine finalizeVariableECSTransition 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule variable_ecs_transition_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
