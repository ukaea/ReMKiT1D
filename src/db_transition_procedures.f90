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
submodule (db_transition_class) db_transition_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the fixed energy detailed balance transition class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
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

    real(rk) ,allocatable ,dimension(:) :: rateVec

    if (assertions) then 
        call assertPure(refVSpace%isDefined(),&
        "Undefined target VSpace object passed to fixed energy detailed balance transition constructor")
        if (present(momentumMoment)) then 
            if (momentumMoment) call assertPure(present(l1Index),"When momentumMoment is set to true in DBTransition &
            &cosntructor the l1Index must be provided")
        end if
        call assertPure(present(degeneracyFun) .eqv. present(degeneracyFunIndices),"If a degeneracy function derivation object &
        &is passed to the DBTransition so must the corresponding variable indices and vice versa")
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

    this%directTransitionIndex = directTransitionIndex
    this%fixedWIndexDirect = fixedWIndexDirect 
    this%temperatureVarIndex = temperatureVarIndex
    this%maxCSl = maxCSl
    call this%setStates(inStates,outStates)

    this%vSpaceRef => refVSpace 

    this%degeneracyRatio = degeneracyRatio 
    if (present(degeneracyFun)) then 
        allocate(this%degeneracyFun,source=degeneracyFun)
        this%degeneracyFunIndices = degeneracyFunIndices 
    end if

    this%csUpdatePriority = 0

    if (present(csUpdatePriority)) this%csUpdatePriority = csUpdatePriority

    call this%setCSDim(2)

    call this%setIncludeElectronDensity(.true.)

    this%fixedWIndex = fixedWIndex 

    this%strictDB = .true.

    if (present(strictDB)) this%strictDB = strictDB

    call this%makeDefined()

end subroutine initDBTransition
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getEnergy(this) result(energyCost)
    !! Returns array representing energy cost of this transition 

    class(DBTransition)                ,intent(in) :: this 
    real(rk) ,allocatable ,dimension(:)            :: energyCost

    if (assertions) call assertPure(this%isDefined(),"getEnergy called on undefined DBTransition object")

    allocate(energyCost(1))
    energyCost = this%transitionEnergy

end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateCSRates(this,varCont,hostModel,hostData,updatePriority)
    !! Update cross-section and transition and moment exchange rates (if applicable)

    class(DBTransition)             ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused
    class(ModelboundData) ,optional ,intent(in)     :: hostData !! Host data used to access direct transition data 
    integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call used to determine if cross-sections should be updated

    real(rk) ,allocatable ,dimension(:)   :: rateBuffer ,emissionVec ,dbCorrection
    real(rk) ,allocatable ,dimension(:)   :: directCrossSection , tempVar ,degenFun ,vGrid ,correctedCS
    real(rk) ,allocatable ,dimension(:,:) :: thisCrossSection 

    integer(ik) :: i, j ,l ,csIndex ,usedPriority

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update undefined DBTransition object")
        call assert(varCont%isDefined(),"Attempted to update DBTransition object using undefined variable container")
    end if

    ! Calculate cross-section harmonic using detailed balance and direct process data (see equation 4.41 in Mijin thesis)

    usedPriority = this%csUpdatePriority
    if (present(updatePriority)) usedPriority = updatePriority
    if (usedPriority >= this%csUpdatePriority) then 

        allocate(tempVar,source = varCont%variables(this%temperatureVarIndex)%entry)
        allocate(degenFun,mold=tempVar)
        degenFun = real(1,kind=rk)
        if (allocated(this%degeneracyFun)) degenFun = this%degeneracyFun%calculate(varCont%variables,this%degeneracyFunIndices)
        vGrid = this%VSpaceRef%getVGrid()
        select type (hostData)
        type is (ModelboundCRMData)
                    emissionVec = hostData%getFixedEmissionVector(this%fixedWIndex)
                class default 
                    error stop "incompatible modelbound data passed to DBTransition update routine"
        end select

        if (.not. allocated(this%directFixedW)) then 

            select type (hostData)
                type is (ModelboundCRMData)
                    allocate(this%directFixedW)
                    this%directFixedW = hostData%getFixedW(this%fixedWIndexDirect) ! Fixed weights of direct process
                class default 
                    error stop "incompatible modelbound data passed to DBTransition update routine"
            end select

            allocate(this%discreteEnergyErrors(size(this%directFixedW%rowIndex)))
            do j = 1,size(this%directFixedW%rowIndex)
                this%discreteEnergyErrors(j)%entry=vGrid(this%directFixedW%columnVector(j)%entry)**2&
                -vGrid(this%directFixedW%rowIndex(j))**2+this%transitionEnergy
            end do
        end if

        do l = 1, this%maxCSl+1

            select type (hostData)
            type is (ModelboundCRMData)
                directCrossSection = hostData%getTransitionCrossSection(l,this%directTransitionIndex)
            class default 
                error stop "incompatible modelbound data passed to DBTransition update routine"
            end select

            if (.not. allocated(thisCrossSection)) then 
                allocate(thisCrossSection(this%locNumX*size(directCrossSection),1+this%maxCSl))
            end if
            thisCrossSection = 0

            do i = 1, this%locNumX
                do j = 1,size(this%directFixedW%rowIndex)
                    csIndex = (i-1)*size(vGrid) + this%directFixedW%rowIndex(j)
                    thisCrossSection(csIndex,l) = dot_product(this%directFixedW%values(j)%entry,&
                    vGrid(this%directFixedW%columnVector(j)%entry)*directCrossSection(this%directFixedW%columnVector(j)%entry)&
                    *exp(-(this%discreteEnergyErrors(j)%entry)/tempVar(i)))
                end do
                thisCrossSection((i-1)*size(vGrid)+1:i*size(vGrid),l) = &
                thisCrossSection((i-1)*size(vGrid)+1:i*size(vGrid),l) * this%degeneracyRatio &
                * degenFun(i)  /vGrid 
            end do

        end do

        !Handle error in detailed balance due to some cells potentially not emitting by rescaling the cross-section locally
        correctedCS = thisCrossSection(:,1)
        do i = 1, this%locNumX
            correctedCS((i-1)*size(vGrid)+1:i*size(vGrid)) = &
            correctedCS((i-1)*size(vGrid)+1:i*size(vGrid)) * emissionVec
        end do  

        dbCorrection = this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,1,1,&
            thisCrossSection(:,1),.true.)&
            /this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,1,1,&
            correctedCS,.true.)
        
        if (.not. this%strictDB) dbCorrection = real(1,kind=rk)
        do l = 1, this%maxCSl+1

            do i = 1, this%locNumX
                thisCrossSection((i-1)*size(vGrid)+1:i*size(vGrid),l) = &
                thisCrossSection((i-1)*size(vGrid)+1:i*size(vGrid),l) * dbCorrection(i) * emissionVec
            end do
        end do
        call this%setCrossSection(thisCrossSection)

    end if
    
    rateBuffer = this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,1,1,&
                                                this%getCrossSectionCol(1),.true.)
    call this%setRate(rateBuffer)
    call this%setRateEnergy(rateBuffer*this%transitionEnergy) 

    if (this%takeMomentumMoment) then

        rateBuffer = this%vSpaceRef%calculateMoment(varCont%variables(this%distFunVarIndex)%entry,&
                                                    this%l1HarmonicInd,2,this%getCrossSectionCol(2),.true.)

        call this%setRateMomentum(rateBuffer)
        
    end if

end subroutine updateCSRates 
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeDBTransition(this) 
    !! Deallocate pointer component

    type(DBTransition)                    ,intent(inout) :: this

    nullify(this%vSpaceRef)

end subroutine finalizeDBTransition 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule db_transition_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
