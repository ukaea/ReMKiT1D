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
submodule(janev_crm_data_support) janev_crm_data_support_procedures
!! author: Stefan Mijin 
!! 
!! Contains implementation of Janev CRM support routines 

implicit none

contains
!-----------------------------------------------------------------------------------------------------------------------------------  
module subroutine addJanevRadRecombTransition(crmData,locNumX,endState,temperatureVarIndex,eVTempNorm,densNorm,timeNorm)
    !! Adds Janev radiative recombination transition to endState of hydrogen. Normalizes the rate to densNorm/timeNorm for 
    !! easier use.
    
    type(ModelboundCRMData) ,intent(inout) :: crmData 
    integer(ik)             ,intent(in)    :: locNumX 
    integer(ik)             ,intent(in)    :: endState 
    integer(ik)             ,intent(in)    :: temperatureVarIndex
    real(rk)                ,intent(in)    :: eVTempNorm
    real(rk)                ,intent(in)    :: densNorm
    real(rk)                ,intent(in)    :: timeNorm
    
    type(FunWrapperDerivation1I1) :: funWrapper 
    type(DerivedTransition) :: transitionObj

    real(rk) :: multConstNorm ,transitionEnergy

    if (assertions .or. assertionLvl >= 0) call assert(crmData%isDefined(),&
    "Undefined CRM Data passed to addJanevRadRecombTransition")

    multConstNorm = real(1.0d-20,kind=rk) * timeNorm*densNorm !Prefactor comes from fact Janev routine returns rate in 10^-14 cm^3/s 

    call funWrapper%init(radRecombRateHydrogen,endState,multConst=multConstNorm,multConstInner=eVTempNorm)

    transitionEnergy = real(13.6d00,kind=rk)/real(endState**2,kind=rk)
    transitionEnergy = transitionEnergy/eVTempNorm

    call transitionObj%init(locNumX,[0,-1],[endState],transitionEnergy,funWrapper,[temperatureVarIndex])

    call crmData%addTransition(transitionObj)

end subroutine addJanevRadRecombTransition
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addJanevCollExIonTransition(crmData,locNumX,startState,endState,distVarIndex,&
                                              eVTempNorm,csNorm,refVSpace,fixedWIndex)
    !! Adds Janev collisional excitation/ionization transition from startState to endState of hydrogen.
    !! Normalizes the cross-section for easier later use. Cross-section only isotropic. endState = 0 interpreted as ionization. 
    !! Assumes velocity is normalized to electron thermal velocity so that me*v_0^2/2 = e * eVTempNorm 
    
    type(ModelboundCRMData) ,intent(inout) :: crmData 
    integer(ik)             ,intent(in)    :: locNumX 
    integer(ik)             ,intent(in)    :: startState 
    integer(ik)             ,intent(in)    :: endState 
    integer(ik)             ,intent(in)    :: distVarIndex
    real(rk)                ,intent(in)    :: eVTempNorm
    real(rk)                ,intent(in)    :: csNorm
    type(VSpace) ,target    ,intent(inout) :: refVSpace !! Target for the transition reference pointer
    integer(ik)             ,intent(in)    :: fixedWIndex !! Index of the inelastic weight array corresponding to this transition

    type (FixedECSTransition) :: transitionObj

    real(rk) ,allocatable ,dimension(:,:) :: crossSection 

    real(rk) ,allocatable ,dimension(:) :: eGrid

    integer(ik) ,allocatable ,dimension(:) :: inStates, outStates

    real(rk) :: transitionEnergy
    integer(ik) :: numV

    if (assertions .or. assertionLvl >= 0) call assert(crmData%isDefined(),&
    "Undefined CRM Data passed to addJanevCollExIonTransition")

    eGrid = refVSpace%getVGrid() ** 2 * eVTempNorm
    numV = refVSpace%getNumV()

    allocate(crossSection(numV,1))

    if (endState == 0) then 

        crossSection(:,1) = real(1.0d-20,kind=rk)/csNorm * ionizationCrossSectionHydrogen(eGrid,startState)
        inStates = [0,startState]
        outStates = [0,0,-1]

        transitionEnergy = real(13.6d00,kind=rk)/real(startState**2,kind=rk)
        transitionEnergy = transitionEnergy/eVTempNorm
    else

        if (assertions .or. assertionLvl >= 0) call assert(endState>startState,&
        "endState must be greater than startState in addJanevCollExIonTransition")

        crossSection(:,1) = real(1.0d-20,kind=rk)/csNorm * excitationCrossSectionHydrogen(eGrid,startState,endState)

        inStates = [0,startState]
        outStates = [0,endState]

        transitionEnergy=real(13.6d00,kind=rk) * (real(startState,kind=rk)**(-2) - real(endState,kind=rk)** (-2))
        transitionEnergy = transitionEnergy/eVTempNorm

    end if
    
    call transitionObj%init(locNumX,inStates,outStates,transitionEnergy,crossSection,distVarIndex,refVSpace,fixedWIndex)

    call crmData%addTransition(transitionObj)

end subroutine addJanevCollExIonTransition
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addJanevCollDeexRecombTransition(crmData,locNumX,startState,endState,distVarIndex,&
    tempVarIndex,eVTempNorm,densNorm,refVSpace,&
    directFixedWIndex,fixedWIndex,directTransitionIndex,csUpdatePriority)
    !! Adds Janev collisional deexcitation/recombination transition from startState to endState of hydrogen.
    !! Normalizes the cross-section for easier later use. Cross-section only isotropic.

    type(ModelboundCRMData) ,intent(inout) :: crmData 
    integer(ik)             ,intent(in)    :: locNumX 
    integer(ik)             ,intent(in)    :: startState 
    integer(ik)             ,intent(in)    :: endState 
    integer(ik)             ,intent(in)    :: distVarIndex
    integer(ik)             ,intent(in)    :: tempVarIndex
    real(rk)                ,intent(in)    :: eVTempNorm
    real(rk)                ,intent(in)    :: densNorm
    type(VSpace) ,target    ,intent(inout) :: refVSpace !! Target for the transition reference pointer
    integer(ik)             ,intent(in)    :: directFixedWIndex !! Index of the inelastic weight array corresponding to 
                    !! the direct counterpart of this transition
    integer(ik)             ,intent(in)    :: fixedWIndex !! Index of the inelastic weight array corresponding to this transition
    integer(ik)             ,intent(in)    :: directTransitionIndex !! Index of the direct transition in crmData
    integer(ik)             ,intent(in)    :: csUpdatePriority !! Update priority for detailed balance cross-section

    type (DBTransition) :: transitionObj

    integer(ik) ,allocatable ,dimension(:) :: inStates, outStates

    real(rk) :: transitionEnergy ,degeneracyRatio ,degenFunMultConst

    type(SimpleDerivation) :: degeneracyFun

    if (assertions .or. assertionLvl >= 0) call assert(crmData%isDefined(),&
    "Undefined CRM Data passed to addJanevCollDeexRecombTransition")

    if (startState == 0) then 

        inStates = [0,0,-1]
        outStates = [0,endState]

        transitionEnergy = -real(13.6d00,kind=rk)/real(endState**2,kind=rk)
        transitionEnergy = transitionEnergy/eVTempNorm

        degeneracyRatio = real(endState**2,kind=rk)

        degenFunMultConst = hPlanck**2/(2*pi*elMass*elCharge*eVTempNorm)
        degenFunMultConst = densNorm * degenFunMultConst ** real(1.5d0,kind=rk) !Included density normalization here so that units agree but left 
                                                                                ! out normalized secondary electron density dependence
        call degeneracyFun%init([real(-1.5d0,kind=rk)],multConst=degenFunMultConst)

        call transitionObj%init(locNumX,inStates,outStates,transitionEnergy,distVarIndex,refVSpace, &
                                directTransitionIndex, directFixedWIndex,fixedWIndex,&
                                tempVarIndex, 0,degeneracyRatio,&
                                degeneracyFun,degeneracyFunIndices=[tempVarIndex],csUpdatePriority=csUpdatePriority)

    else

        if (assertions .or. assertionLvl >= 0) call assert(endState<startState,&
        "endState must be lower than startState in addJanevCollExIonTransition")

        inStates = [0,startState]
        outStates = [0,endState]

        transitionEnergy= real(13.6d00,kind=rk) * (real(startState,kind=rk)**(-2) - real(endState,kind=rk)** (-2))
        transitionEnergy = transitionEnergy/eVTempNorm

        degeneracyRatio = real(endState**2,kind=rk)/real(startState**2,kind=rk)

        call transitionObj%init(locNumX,inStates,outStates,transitionEnergy,distVarIndex,refVSpace, &
                                directTransitionIndex, directFixedWIndex,fixedWIndex,&
                                tempVarIndex, 0,degeneracyRatio,csUpdatePriority=csUpdatePriority)

    end if

    call crmData%addTransition(transitionObj)

end subroutine addJanevCollDeexRecombTransition
!-----------------------------------------------------------------------------------------------------------------------------------  
end submodule janev_crm_data_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------  
