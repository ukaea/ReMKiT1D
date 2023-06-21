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
submodule (modelbound_data_support) modelbound_data_support_procedures
!! author: Stefan Mijin 
!! 
!! Contains modelbound data support procedures

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addModelboundDataToModel(modelObj,modelTag,envObj,normObj)
    !! Initialize modelbound data and add to corresponding model object

    type(Model)                ,intent(inout) :: modelObj
    character(*)               ,intent(in)    :: modelTag 
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    class(Normalization)        ,intent(in)    :: normObj
    
    type(NamedString) ,dimension(1) :: dataType 

    dataType(1) = NamedString(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyModelboundDataType,"")

    call envObj%jsonCont%load(dataType)
    call envObj%jsonCont%output(dataType)
    select case(dataType(1)%value)
    case (keyVarLikeMB)
        call addCustomVarlikeMBDataToModel(modelObj,modelTag,envObj)
    case (keyCRMData)
        call addCustomCRMDataToModel(modelObj,modelTag,envObj,normObj)
    case (keyLBCData)
        call addCustomLBCDataToModel(modelObj,modelTag,envObj,normObj)
    case default
        error stop "Unsupported modelbound data type detected by addModelboundDataToModel"
    end select

end subroutine addModelboundDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addCustomVarlikeMBDataToModel(modelObj,modelTag,envObj)
    !! Initialize varlike modelbound data and add to corresponding model object

    type(Model)                ,intent(inout) :: modelObj
    character(*)               ,intent(in)    :: modelTag 
    type(EnvironmentWrapper)   ,intent(inout) :: envObj

    type(NamedStringArray) ,dimension(1) :: dataNames ,reqVars
    type(NamedLogical)     ,dimension(4) :: logicalParams 
    type(NamedString)      ,dimension(1) :: ruleName

    type(NamedInteger)     ,dimension(1) :: varPriority
    
    type(ModelboundDataVarlike) :: mbData 
    type(CalculationRule) ,allocatable ,dimension(:) :: cRules
    type(VariableList)                               :: varList 
    class(Derivation) ,allocatable                   :: method
    integer(ik) :: i 
    integer(ik) ,allocatable ,dimension(:) :: dataDerivIndices

    dataNames(1)%name = keyModels//"."//modelTag//"."//keyModelboundData//"."//keyDataNames
    allocate(dataNames(1)%values(0))

    call envObj%jsonCont%load(dataNames)
    call envObj%jsonCont%output(dataNames)

    call varList%init()
    allocate(cRules(size(dataNames(1)%values)))
    allocate(dataDerivIndices(0))

    do i = 1, size(dataNames(1)%values) 
        
        logicalParams(1) = &
        NamedLogical(keyModels//"."//modelTag//"."//keyModelboundData//"."//dataNames(1)%values(i)%string//"."//keyIsDist,.false.)

        logicalParams(2) = &
        NamedLogical(keyModels//"."//modelTag//"."//keyModelboundData//"."//dataNames(1)%values(i)%string//"."//keyIsScalar,.false.)

        logicalParams(3) = NamedLogical(keyModels//"."//modelTag//"."//keyModelboundData//"."//dataNames(1)%values(i)%string//"."//&
                                        keyIsSingleHarmonic,.false.)

        logicalParams(4) = NamedLogical(keyModels//"."//modelTag//"."//keyModelboundData//"."//dataNames(1)%values(i)%string//"."//&
                                        keyIsDerivedFromOtherData,.false.)
        
        ruleName(1) = NamedString(keyModels//"."//modelTag//"."//keyModelboundData//"."//dataNames(1)%values(i)%string//"."//&
                                  keyRuleName,"")

        reqVars(1)%name = keyModels//"."//modelTag//"."//keyModelboundData//"."//dataNames(1)%values(i)%string//"."//keyReqVarNames
        allocate(reqVars(1)%values(0))

        varPriority(1) = NamedInteger(keyModels//"."//modelTag//"."//keyModelboundData//&
                                     "."//dataNames(1)%values(i)%string//"."//keyDerivationPriority,0)

        call envObj%jsonCont%load(ruleName)
        call envObj%jsonCont%output(ruleName)

        call envObj%jsonCont%load(logicalParams)
        call envObj%jsonCont%output(logicalParams)

        call envObj%jsonCont%load(reqVars)
        call envObj%jsonCont%output(reqVars)

        call envObj%jsonCont%load(varPriority)
        call envObj%jsonCont%output(varPriority)

        call varList%addVar(dataNames(1)%values(i)%string,isDist=logicalParams(1)%value,&
                            isScalar=logicalParams(2)%value,isSingleHarmonic=logicalParams(3)%value,priority=varPriority(1)%value)

        if (logicalParams(4)%value) dataDerivIndices = [dataDerivIndices,i]
        call envObj%textbookObj%copyDerivation(ruleName(1)%value,method)
        call cRules(i)%init(method,reqVars(1)%values)
        deallocate(reqVars(1)%values)

    end do

    call mbData%init(varList,cRules,envObj%partitionObj,envObj%indexingObj,envObj%mpiCont%getXHaloWidth(),&
                    envObj%externalVars,envObj%mpiCont%getWorldRank(),dataDerivIndices=dataDerivIndices)

    call modelObj%setModelData(mbData)

end subroutine addCustomVarlikeMBDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addCustomCRMDataToModel(modelObj,modelTag,envObj,normObj)
    !! Initialize custom CRM modelbound data and add to corresponding model object

    type(Model)                ,intent(inout) :: modelObj
    character(*)               ,intent(in)    :: modelTag 
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    class(Normalization)       ,intent(in)    :: normObj

    type(ModelboundCRMData) :: crmData
    type(InelasticGridData) :: inelData

    type(NamedStringArray) ,dimension(1) :: transitionTags
    type(NamedString) ,allocatable ,dimension(:) :: transitionType
    type(NamedRealArray) ,dimension(1) :: fixedTransitionEnergies
    type(NamedLogical) ,dimension(1) :: activeInelData

    integer(ik) :: i

    transitionTags(1)%name = keyModels//"."//modelTag//"."//keyModelboundData//"."//keyTransitionTags
    allocate(transitionTags(1)%values(0))

    call envObj%jsonCont%load(transitionTags)
    call envObj%jsonCont%output(transitionTags)

    allocate(transitionType(size(transitionTags(1)%values)))

    activeInelData(1) = NamedLogical(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyInelGrid//"."//keyActive,.false.)

    fixedTransitionEnergies(1)%name = &
    keyModels//"."//modelTag//"."//keyModelboundData//"."//keyInelGrid//"."//keyFixedTransitionEnergies
    allocate(fixedTransitionEnergies(1)%values(0))

    call envObj%jsonCont%load(activeInelData)
    call envObj%jsonCont%output(activeInelData)

    call envObj%jsonCont%load(fixedTransitionEnergies)
    call envObj%jsonCont%output(fixedTransitionEnergies)

    if (activeInelData(1)%value) then 
        call inelData%init(envObj%vSpaceObj,fixedEnergies=fixedTransitionEnergies(1)%values)
        call crmData%setInelData(inelData)
    end if

    do i  = 1, size(transitionTags(1)%values)
        transitionType(i) = NamedString(keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                        keyTransitions//"."//transitionTags(1)%values(i)%string//"."//keyType,"")
    end do

    call envObj%jsonCont%load(transitionType)
    call envObj%jsonCont%output(transitionType)

    call crmData%init(size(transitionTags(1)%values))

    do i = 1,size(transitionTags(1)%values)

        select case (transitionType(i)%value)
        case (keySimpleTransition)
            call addSimpleTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                             keyTransitions//"."//transitionTags(1)%values(i)%string,&
                                                             fixedTransitionEnergies(1)%values)
        case (keyDerivedTransition)
            call addDerivedTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                             keyTransitions//"."//transitionTags(1)%values(i)%string,&
                                                             fixedTransitionEnergies(1)%values)
        case (keyFixedECSTransition)
            call addFixedECSTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                             keyTransitions//"."//transitionTags(1)%values(i)%string,&
                                                             fixedTransitionEnergies(1)%values)
        case (keyVariableECSTransition)
            call addVariableECSTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                                keyTransitions//"."//transitionTags(1)%values(i)%string)
        case (keyDBTransition)
            call addDBTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                                keyTransitions//"."//transitionTags(1)%values(i)%string,&
                                                                fixedTransitionEnergies(1)%values)
        case (keyJanevRadRecomb)
            call addJanevRadRecombTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                      keyTransitions//"."//transitionTags(1)%values(i)%string,normObj)
        case (keyJanevCollExIon)
            call addJanevCollExIonTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                      keyTransitions//"."//transitionTags(1)%values(i)%string,normObj)
        case (keyJanevCollDeexRecomb)
            call addJanevCollDeexRecombTransitionToCRMData(crmData,envObj,keyModels//"."//modelTag//"."//keyModelboundData//"."//&
                                                            keyTransitions//"."//transitionTags(1)%values(i)%string,normObj)
        case default 
            error stop "unsupported transition type detected in custom modelbound data"
        end select

    end do

    call modelObj%setModelData(crmData)

end subroutine addCustomCRMDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addCustomLBCDataToModel(modelObj,modelTag,envObj,normObj)
    !! Initialize custom LBC modelbound data and add to corresponding model object

    type(Model)                ,intent(inout) :: modelObj
    character(*)               ,intent(in)    :: modelTag 
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    class(Normalization)       ,intent(in)    :: normObj

    type(NamedLogical) ,dimension(1) :: leftBoundary
    type(NamedStringArray) ,dimension(1) :: derivReqVars
    type(NamedString) ,dimension(1) :: derivName ,ionCurrentVar, totalCurrentVar
    type(NamedReal) ,dimension(1) :: bisTol

    integer(ik) ,allocatable ,dimension(:) :: derivReqIndices

    integer(ik) :: i

    class(MatDerivation) ,allocatable :: derivObj

    type(ModelboundLBCData) :: mbData

    logical :: isActive


    derivName(1) = NamedString(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyRuleName,"")
    leftBoundary(1) = NamedLogical(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyLeftBoundary,.false.)
    derivReqVars(1)%name = keyModels//"."//modelTag//"."//keyModelboundData//"."//keyReqVarNames
    allocate(derivReqVars(1)%values(0))

    ionCurrentVar(1) = NamedString(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyIonCurrentVar,"")
    totalCurrentVar(1) = NamedString(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyTotalCurrentVar,keyNone)
    bisTol(1) = NamedReal(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyBisTol,real(1.0d-12,kind=rk))

    call envObj%jsonCont%load(leftBoundary)
    call envObj%jsonCont%output(leftBoundary)
    call envObj%jsonCont%load(derivReqVars)
    call envObj%jsonCont%output(derivReqVars)
    call envObj%jsonCont%load(derivName)
    call envObj%jsonCont%output(derivName)
    call envObj%jsonCont%load(ionCurrentVar)
    call envObj%jsonCont%output(ionCurrentVar)
    call envObj%jsonCont%load(totalCurrentVar)
    call envObj%jsonCont%output(totalCurrentVar)
    call envObj%jsonCont%load(bisTol)
    call envObj%jsonCont%output(bisTol)

    allocate(derivReqIndices(size(derivReqVars(1)%values)))
    do i = 1,size(derivReqIndices)
        derivReqIndices(i) = envObj%externalVars%getVarIndex(derivReqVars(1)%values(i)%string)
    end do

    call envObj%textbookObj%copyMatDerivation(derivName(1)%value,derivObj)

    isActive = envObj%partitionObj%getMaxXAtInd(envObj%mpiCont%getWorldRank()+1) == envObj%gridObj%getNumX()
    if (leftBoundary(1)%value) isActive = envObj%partitionObj%getMinXAtInd(envObj%mpiCont%getWorldRank()+1) == 1
    select type (derivObj)
    type is (FScalingDerivation)

        if (totalCurrentVar(1)%value /= keyNone) then 
            call mbData%init(envObj%vSpaceObj,derivObj,derivReqIndices,envObj%externalVars%getVarIndex(ionCurrentVar(1)%value),&
                             isActive ,totalCurrentVarIndex= envObj%externalVars%getVarIndex(totalCurrentVar(1)%value),&
                             bisTol=bisTol(1)%value,isLeftBoundary=leftBoundary(1)%value)
        else
            call mbData%init(envObj%vSpaceObj,derivObj,derivReqIndices,envObj%externalVars%getVarIndex(ionCurrentVar(1)%value),&
                             isActive ,bisTol=bisTol(1)%value,isLeftBoundary=leftBoundary(1)%value)
        end if
    class default 
        error stop "Derivation object detected by addCustomLBCData to model must be of type FScalingDerivation"
    end select

    call modelObj%setModelData(mbData)
    
end subroutine addCustomLBCDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addSimpleTransitionToCRMData(crmData,envObj,tTag,fixedEnergies)
    !! Add custom simple transition to CRM data object based on JSON file

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag 
    real(rk) ,dimension(:)     ,intent(in)    :: fixedEnergies

    type(SimpleTransition) :: transitionObj 

    type(NamedInteger) ,dimension(1) :: inState, outState ,fixedEnergyIndex
    type(NamedReal) ,dimension(1) :: fixedEnergy ,fixedRate

    real(rk) :: usedEnergy
    
    integer(ik) :: locNumX ,locRank

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 

    inState(1) = NamedInteger(tTag//"."//keyIngoingState,1)
    outState(1) = NamedInteger(tTag//"."//keyOutgoingState,1)
    fixedEnergyIndex(1) = NamedInteger(tTag//"."//keyFixedEnergyIndex,0)

    fixedEnergy(1) = NamedReal(tTag//"."//keyFixedEnergy,real(0,kind=rk))
    fixedRate(1) = NamedReal(tTag//"."//keyRate,real(0,kind=rk))

    call envObj%jsonCont%load(inState)
    call envObj%jsonCont%output(inState)
    call envObj%jsonCont%load(outState)
    call envObj%jsonCont%output(outState)
    call envObj%jsonCont%load(fixedEnergyIndex)
    call envObj%jsonCont%output(fixedEnergyIndex)

    call envObj%jsonCont%load(fixedEnergy)
    call envObj%jsonCont%output(fixedEnergy)
    call envObj%jsonCont%load(fixedRate)
    call envObj%jsonCont%output(fixedRate)

    usedEnergy = fixedEnergy(1)%value 

    if (fixedEnergyIndex(1)%value > 0) usedEnergy = fixedEnergies(fixedEnergyIndex(1)%value)

    call transitionObj%init(locNumX,inState(1)%value,outState(1)%value,usedEnergy,fixedRate(1)%value)

    call crmData%addTransition(transitionObj)

end subroutine addSimpleTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addDerivedTransitionToCRMData(crmData,envObj,tTag,fixedEnergies)
    !! Add custom derived transition to CRM data object based on JSON file

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag 
    real(rk) ,dimension(:)     ,intent(in)    :: fixedEnergies

    type(DerivedTransition) :: transitionObj 

    type(NamedIntegerArray) ,dimension(1) :: inStates, outStates
    type(NamedInteger) ,dimension(1) :: fixedEnergyIndex
    type(NamedReal) ,dimension(1) :: fixedEnergy 
    type(NamedString) ,dimension(3) :: ruleName 
    type(NamedStringArray) ,dimension(3) :: reqVars

    class(Derivation) ,allocatable :: deriv ,deriv2, deriv3

    real(rk) :: usedEnergy
    
    integer(ik) :: locNumX ,locRank, i

    integer(ik) ,allocatable ,dimension(:) :: reqVarIndices,reqVarIndices2,reqVarIndices3

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 

    inStates(1) = NamedIntegerArray(tTag//"."//keyIngoingStatePlural,[1])
    outStates(1) = NamedIntegerArray(tTag//"."//keyOutgoingStatePlural,[1])
    fixedEnergyIndex(1) = NamedInteger(tTag//"."//keyFixedEnergyIndex,0)

    fixedEnergy(1) = NamedReal(tTag//"."//keyFixedEnergy,real(0,kind=rk))

    ruleName(1) = NamedString(tTag//"."//keyRuleName,"")

    reqVars(1)%name = tTag//"."//keyReqVarNames
    allocate(reqVars(1)%values(0))

    ruleName(2) = NamedString(tTag//"."//keyMomentumRule,keyNone)
    reqVars(2)%name = tTag//"."//keyMomentumReqVarNames
    allocate(reqVars(2)%values(0))

    ruleName(3) = NamedString(tTag//"."//keyEnergyRule,keyNone)
    reqVars(3)%name = tTag//"."//keyEnergyReqVarNames
    allocate(reqVars(3)%values(0))

    call envObj%jsonCont%load(inStates)
    call envObj%jsonCont%output(inStates)
    call envObj%jsonCont%load(outStates)
    call envObj%jsonCont%output(outStates)
    call envObj%jsonCont%load(fixedEnergyIndex)
    call envObj%jsonCont%output(fixedEnergyIndex)

    call envObj%jsonCont%load(fixedEnergy)
    call envObj%jsonCont%output(fixedEnergy)

    call envObj%jsonCont%load(ruleName)
    call envObj%jsonCont%output(ruleName)

    call envObj%jsonCont%load(reqVars)
    call envObj%jsonCont%output(reqVars)

    usedEnergy = fixedEnergy(1)%value 

    if (fixedEnergyIndex(1)%value > 0) usedEnergy = fixedEnergies(fixedEnergyIndex(1)%value)

    allocate(reqVarIndices(size(reqVars(1)%values)))

    do i = 1,size(reqVars(1)%values)
        reqVarIndices(i) = envObj%externalVars%getVarIndex(reqVars(1)%values(i)%string)
    end do

    allocate(reqVarIndices2(size(reqVars(2)%values)))

    do i = 1,size(reqVars(2)%values)
        reqVarIndices2(i) = envObj%externalVars%getVarIndex(reqVars(2)%values(i)%string)
    end do

    allocate(reqVarIndices3(size(reqVars(3)%values)))

    do i = 1,size(reqVars(3)%values)
        reqVarIndices3(i) = envObj%externalVars%getVarIndex(reqVars(3)%values(i)%string)
    end do

    call envObj%textbookObj%copyDerivation(ruleName(1)%value,deriv)
    if (ruleName(2)%value /= keyNone) call envObj%textbookObj%copyDerivation(ruleName(2)%value,deriv2)
    if (ruleName(3)%value /= keyNone) call envObj%textbookObj%copyDerivation(ruleName(3)%value,deriv3)

    if (ruleName(2)%value /= keyNone) then
        if (ruleName(3)%value /= keyNone) then
            
            call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,usedEnergy,deriv,reqVarIndices,&
            momentumRateDeriv=deriv2,momentumRateDerivIndices=reqVarIndices2,&
            energyRateDeriv=deriv3,energyRateDerivIndices=reqVarIndices3)
        else
            call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,usedEnergy,deriv,reqVarIndices,&
            momentumRateDeriv=deriv2,momentumRateDerivIndices=reqVarIndices2)
        end if

    else

        if (ruleName(3)%value /= keyNone) then
            call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,usedEnergy,deriv,reqVarIndices,&
            energyRateDeriv=deriv3,energyRateDerivIndices=reqVarIndices3)
        else
            call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,usedEnergy,deriv,reqVarIndices)
        end if
    end if

    call crmData%addTransition(transitionObj)

end subroutine addDerivedTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addFixedECSTransitionToCRMData(crmData,envObj,tTag,fixedEnergies)
    !! Add custom fixed energy/cross-section transition to CRM data object based on JSON file

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag 
    real(rk) ,dimension(:)     ,intent(in)    :: fixedEnergies

    type(FixedECSTransition) :: transitionObj 

    type(NamedIntegerArray) ,dimension(1) :: inStates, outStates ,csLIndices
    type(NamedInteger) ,dimension(1) :: fixedEnergyIndex ,maxCSl
    type(NamedString) ,dimension(1) :: distributionVarName
    type(NamedRealArray) ,allocatable ,dimension(:) :: csValues 
    type(NamedLogical) ,dimension(1) :: takeMomentumMoment

    integer(ik) :: locNumX ,locRank ,distVarIndex ,i ,numV ,maxL ,l1Index

    character(len=30) :: intToStrBuffer

    real(rk) ,allocatable ,dimension(:,:) :: crossSection

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 
    numV = envObj%gridObj%getNumV()
    maxL = envObj%gridObj%getMaxL()

    inStates(1) = NamedIntegerArray(tTag//"."//keyIngoingStatePlural,[1])
    outStates(1) = NamedIntegerArray(tTag//"."//keyOutgoingStatePlural,[1])
    csLIndices(1)%name = tTag//"."//keyCrossSectionData//"."//keyPresentHarmonics
    allocate(csLIndices(1)%values(0))

    fixedEnergyIndex(1) = NamedInteger(tTag//"."//keyFixedEnergyIndex,0)
    maxCSl(1) = NamedInteger(tTag//"."//keyCrossSectionData//"."//keyMaxCrossSectionL,0)

    call envObj%jsonCont%load(inStates)
    call envObj%jsonCont%output(inStates)
    call envObj%jsonCont%load(outStates)
    call envObj%jsonCont%output(outStates)
    call envObj%jsonCont%load(fixedEnergyIndex)
    call envObj%jsonCont%output(fixedEnergyIndex)

    call envObj%jsonCont%load(maxCSl)

    call envObj%jsonCont%load(csLIndices)

    if (size(csLIndices(1)%values)==0) csLIndices(1)%values=[(i,i=0,maxCSl(1)%value)]

    maxCSl(1)%value = maxval(csLIndices(1)%values)

    call envObj%jsonCont%output(maxCSl)
    call envObj%jsonCont%output(csLIndices)

    allocate(csValues(size(csLIndices(1)%values)))

    do i = 1,size(csLIndices(1)%values)
        intToStrBuffer=""
        write(intToStrBuffer,'(I0)') csLIndices(1)%values(i)
        csValues(i)%name = tTag//"."//keyCrossSectionData//"."//keyLPrefix//trim(intToStrBuffer)
        csValues(i)%values = real([(0,i=1,numV)],kind=rk)
    end do

    call envObj%jsonCont%load(csValues)
    call envObj%jsoncont%output(csValues)

    allocate(crossSection(numV,maxCSl(1)%value))

    crossSection = real(0,kind=rk)

    do i =1,size(csLIndices(1)%values)
        crossSection(:,csLIndices(1)%values(i)+1) = csValues(i)%values
    end do

    distributionVarName(1) = NamedString(tTag//"."//keyDistributionVarName,"")

    call envObj%jsonCont%load(distributionVarName)

    if (assertions) call assert(envObj%externalVars%isVarNAmeRegistered(distributionVarName(1)%value),distributionVarName(1)%value&
                                //" not registered in environment wrapper")

    call envObj%jsonCont%output(distributionVarName)

    distVarIndex = envObj%externalVars%getVarIndex(distributionVarName(1)%value)

    takeMomentumMoment(1) = NamedLogical(tTag//"."//keyTakeMomentumMoment,.false.)

    call envObj%jsonCont%load(takeMomentumMoment)

    if (assertions) then 
        if (takeMomentumMoment(1)%value) call assert(maxCSl(1)%value>0,"Fixed energy/cross-section transition set to take momentum&
                                                    & moment when maximum cross-section harmonic is not 1 or higher")
    end if

    call envObj%jsonCont%output(takeMomentumMoment)

    l1Index = 1
    if (maxL > 0) l1Index = envObj%gridObj%getH(1,0,.false.)

    call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,fixedEnergies(fixedEnergyIndex(1)%value),&
                            crossSection,distVarIndex,envObj%vSpaceObj, &
                            fixedEnergyIndex(1)%value,momentumMoment=takeMomentumMoment(1)%value, l1Index=l1Index)

    call crmData%addTransition(transitionObj)

end subroutine addFixedECSTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addVariableECSTransitionToCRMData(crmData,envObj,tTag)
    !! Add custom variable energy/cross-section transition to CRM data object based on JSON file

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag 

    type(VariableECSTransition) :: transitionObj 

    type(NamedIntegerArray) ,dimension(1) :: inStates, outStates ,csLIndices 
    type(NamedString) ,dimension(1) :: distributionVarName ,energyDerivName
    type(NamedStringArray) ,dimension(1) :: energyDerivationReqVars
    type(NamedString) ,allocatable ,dimension(:) :: csDerivs 
    type(NamedStringArray) ,allocatable ,dimension(:) :: csDerivsReqVars 
    type(NamedLogical) ,dimension(1) :: takeMomentumMoment

    integer(ik) :: locNumX ,locRank ,distVarIndex ,i ,j ,maxL ,l1Index

    type(DerivationContainer) ,allocatable ,dimension(:) :: csDerivCont
    class(Derivation) ,allocatable :: energyDeriv

    integer(ik) ,allocatable ,dimension(:) :: reqVarIndicesEnergy 
    type(IntArray) ,allocatable ,dimension(:) :: reqVarIndicesCS

    character(len=30) :: intToStrBuffer

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 
    maxL = envObj%gridObj%getMaxL()

    inStates(1) = NamedIntegerArray(tTag//"."//keyIngoingStatePlural,[1])
    outStates(1) = NamedIntegerArray(tTag//"."//keyOutgoingStatePlural,[1])
    csLIndices(1)%name = tTag//"."//keyCrossSectionDerivations//"."//keyCSDerivationHarmonics
    allocate(csLIndices(1)%values(0))


    call envObj%jsonCont%load(inStates)
    call envObj%jsonCont%output(inStates)
    call envObj%jsonCont%load(outStates)
    call envObj%jsonCont%output(outStates)

    call envObj%jsonCont%load(csLIndices)
    call envObj%jsonCont%output(csLIndices)

    allocate(csDerivsReqVars(size(csLIndices(1)%values)))
    allocate(csDerivs(size(csLIndices(1)%values)))

    do i = 1,size(csLIndices(1)%values)
        intToStrBuffer=""
        write(intToStrBuffer,'(I0)') csLIndices(1)%values(i)
        csDerivs(i)%name = tTag//"."//keyCrossSectionDerivations//"."//keyLPrefix//trim(intToStrBuffer)//"."//keyRuleName
        csDerivs(i)%value = ""
        csDerivsReqVars(i)%name = tTag//"."//keyCrossSectionDerivations//"."//keyLPrefix//trim(intToStrBuffer)//"."//keyReqVarNames
        allocate(csDerivsReqVars(i)%values(0))
    end do

    distributionVarName(1) = NamedString(tTag//"."//keyDistributionVarName,"")

    call envObj%jsonCont%load(distributionVarName)

    if (assertions) call assert(envObj%externalVars%isVarNameRegistered(distributionVarName(1)%value),distributionVarName(1)%value&
                                //" not registered in environment wrapper")

    call envObj%jsonCont%output(distributionVarName)

    distVarIndex = envObj%externalVars%getVarIndex(distributionVarName(1)%value)

    call envObj%jsonCont%load(csDerivs)
    call envObj%jsonCont%output(csDerivs)

    call envObj%jsonCont%load(csDerivsReqVars)
    call envObj%jsonCont%output(csDerivsReqVars)

    energyDerivName(1) = NamedString(tTag//"."//keyEnergyDerivationName,"")
    energyDerivationReqVars(1)%name = tTag//"."//keyEnergyDerivationReqVars
    allocate(energyDerivationReqVars(1)%values(0))

    call envObj%jsonCont%load(energyDerivName)
    call envObj%jsonCont%output(energyDerivName)

    call envObj%jsonCont%load(energyDerivationReqVars)
    call envObj%jsonCont%output(energyDerivationReqVars)

    allocate(reqVarIndicesEnergy(size(energyDerivationReqVars(1)%values)))

    do i = 1,size(energyDerivationReqVars(1)%values)
        reqVarIndicesEnergy(i) = envObj%externalVars%getVarIndex(energyDerivationReqVars(1)%values(i)%string)
    end do

    call envObj%textbookObj%copyDerivation(energyDerivName(1)%value,energyDeriv)

    allocate(csDerivCont(size(csDerivs)))
    allocate(reqVarIndicesCS(size(csDerivs)))

    do i = 1,size(csDerivs)
        call envObj%textbookObj%copyDerivation(energyDerivName(1)%value,csDerivCont(i)%entry)
        allocate(reqVarIndicesCS(size(csDerivsReqVars(i)%values)))
        do j =1,size(csDerivsReqVars(i)%values)
            reqVarIndicesCS(i)%entry(j) = envObj%externalVars%getVarIndex(csDerivsReqVars(i)%values(j)%string)
        end do
    end do

    takeMomentumMoment(1) = NamedLogical(tTag//"."//keyTakeMomentumMoment,.false.)

    call envObj%jsonCont%load(takeMomentumMoment)

    if (assertions) then 
        if (takeMomentumMoment(1)%value) call assert(maxval(csLIndices(1)%values)>0,&
                                                    "Fixed energy/cross-section transition set to take momentum&
                                                    & moment when maximum cross-section harmonic is not 1 or higher")
    end if

    call envObj%jsonCont%output(takeMomentumMoment)

    l1Index = 1
    if (maxL > 0) l1Index = envObj%gridObj%getH(1,0,.false.)

    call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,energyDeriv,reqVarIndicesEnergy,csDerivCont,&
                            reqVarIndicesCS,csLIndices(1)%values+1,distVarIndex,envObj%vSpaceObj, &
                            momentumMoment=takeMomentumMoment(1)%value, l1Index=l1Index) 

    call crmData%addTransition(transitionObj)

end subroutine addVariableECSTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addDBTransitionToCRMData(crmData,envObj,tTag,fixedEnergies)
    !! Add custom detailed balance transition to CRM data object based on JSON file.

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag 
    real(rk) ,dimension(:)     ,intent(in)    :: fixedEnergies

    type(DBTransition) :: transitionObj 

    type(NamedIntegerArray) ,dimension(1) :: inStates, outStates 
    type(NamedInteger) ,dimension(1) :: fixedEnergyIndex ,maxCSl ,directTransitionIndex &
                                       ,csUpdatePriority,directTransitionFixedEnergyIndex
    type(NamedString) ,dimension(1) :: distributionVarName, temperatureVarName ,degenRule
    type(NamedStringArray) ,dimension(1) :: degenRuleVars
    type(NamedLogical) ,dimension(1) :: takeMomentumMoment
    type(NamedReal) ,dimension(1) :: degeneracyRatio 

    class(Derivation) ,allocatable :: deriv

    integer(ik) :: locNumX ,locRank ,distVarIndex ,temperatureVarIndex,i ,maxL ,l1Index

    integer(ik) ,allocatable ,dimension(:) :: reqVarIndices

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 
    maxL = envObj%gridObj%getMaxL()

    inStates(1) = NamedIntegerArray(tTag//"."//keyIngoingStatePlural,[1])
    outStates(1) = NamedIntegerArray(tTag//"."//keyOutgoingStatePlural,[1])

    fixedEnergyIndex(1) = NamedInteger(tTag//"."//keyDirectTransitionEnergyIndex,0) 
    directTransitionFixedEnergyIndex(1) = NamedInteger(tTag//"."//keyDirectTransitionEnergyIndex,1)
    directTransitionIndex(1) = NamedInteger(tTag//"."//keyDirectTransitionIndex,0)
    maxCSl(1) = NamedInteger(tTag//"."//keyMaxCrossSectionL,0)

    csUpdatePriority(1) = NamedInteger(tTag//"."//keyCSUpdatePriority,0)

    call envObj%jsonCont%load(inStates)
    call envObj%jsonCont%output(inStates)
    call envObj%jsonCont%load(outStates)
    call envObj%jsonCont%output(outStates)
    call envObj%jsonCont%load(fixedEnergyIndex)
    call envObj%jsonCont%output(fixedEnergyIndex)
    call envObj%jsonCont%load(directTransitionFixedEnergyIndex)
    call envObj%jsonCont%output(directTransitionFixedEnergyIndex)
    call envObj%jsonCont%load(csUpdatePriority)
    call envObj%jsonCont%output(csUpdatePriority)

    call envObj%jsonCont%load(directTransitionIndex)
    call envObj%jsonCont%output(directTransitionIndex)

    call envObj%jsonCont%load(maxCSl)

    call envObj%jsonCont%output(maxCSl)


    distributionVarName(1) = NamedString(tTag//"."//keyDistributionVarName,"")
    temperatureVarName(1) = NamedString(tTag//"."//keyElectronTemperatureVar,"")

    call envObj%jsonCont%load(distributionVarName)
    call envObj%jsonCont%load(temperatureVarName)

    if (assertions) then 
        call assert(envObj%externalVars%isVarNameRegistered(distributionVarName(1)%value),distributionVarName(1)%value&
                                //" not registered in environment wrapper")
        call assert(envObj%externalVars%isVarNameRegistered(temperatureVarName(1)%value),temperatureVarName(1)%value&
                                //" not registered in environment wrapper")
    end if
    call envObj%jsonCont%output(distributionVarName)

    distVarIndex = envObj%externalVars%getVarIndex(distributionVarName(1)%value)
    temperatureVarIndex = envObj%externalVars%getVarIndex(temperatureVarName(1)%value)

    takeMomentumMoment(1) = NamedLogical(tTag//"."//keyTakeMomentumMoment,.false.)

    call envObj%jsonCont%load(takeMomentumMoment)

    if (assertions) then 
        if (takeMomentumMoment(1)%value) call assert(maxCSl(1)%value>0,"Fixed energy/cross-section transition set to take momentum&
                                                    & moment when maximum cross-section harmonic is not 1 or higher")
    end if

    call envObj%jsonCont%output(takeMomentumMoment)

    l1Index = 1
    if (maxL > 0) l1Index = envObj%gridObj%getH(1,0,.false.)

    degenRule(1) = NamedString(tTag//"."//keyDegeneracyRule,keyNone)

    degenRuleVars(1)%name = tTag//"."//keyDegeneracyRuleReqVars
    allocate(degenRuleVars(1)%values(0))

    degeneracyRatio = NamedReal(tTag//"."//keyFixedDegenRatio,real(1,kind=rk))

    call envObj%jsonCont%load(degeneracyRatio)
    call envObj%jsonCont%output(degeneracyRatio)

    call envObj%jsonCont%load(degenRule)
    call envObj%jsonCont%output(degenRule)

    call envObj%jsonCont%load(degenRuleVars)
    call envObj%jsonCont%output(degenRuleVars)

    allocate(reqVarIndices(size(degenRuleVars(1)%values)))

    do i = 1,size(degenRuleVars(1)%values)
        reqVarIndices(i) = envObj%externalVars%getVarIndex(degenRuleVars(1)%values(i)%string)
    end do

    if (degenRule(1)%value /= keyNone) then

        call envObj%textbookObj%copyDerivation(degenRule(1)%value,deriv)

        call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,fixedEnergies(fixedEnergyIndex(1)%value),&
                                distVarIndex,envObj%vSpaceObj,directTransitionIndex(1)%value,&
                                directTransitionFixedEnergyIndex(1)%value,fixedEnergyIndex(1)%value,&
                                temperatureVarIndex,maxCSl(1)%value,degeneracyRatio(1)%value,&
                                momentumMoment=takeMomentumMoment(1)%value, l1Index=l1Index,&
                                degeneracyFun=deriv,degeneracyFunIndices=reqVarIndices,csUpdatePriority=csUpdatePriority(1)%value)
    else
        call transitionObj%init(locNumX,inStates(1)%values,outStates(1)%values,fixedEnergies(fixedEnergyIndex(1)%value),&
                                distVarIndex,envObj%vSpaceObj,directTransitionIndex(1)%value,&
                                directTransitionFixedEnergyIndex(1)%value,fixedEnergyIndex(1)%value,&
                                temperatureVarIndex,maxCSl(1)%value,degeneracyRatio(1)%value,&
                                momentumMoment=takeMomentumMoment(1)%value, l1Index=l1Index,&
                                csUpdatePriority=csUpdatePriority(1)%value)
    end if

    call crmData%addTransition(transitionObj)

end subroutine addDBTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addJanevRadRecombTransitionToCRMData(crmData,envObj,tTag,normObj)
    !! Add Janev data radiative recombination transition to CRM data object based on JSON file

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag 
    class(Normalization)       ,intent(in)    :: normObj

    type(NamedInteger) ,dimension(1) :: endState

    type(NamedString) ,dimension(1) :: temperatureVarName
    
    integer(ik) :: locNumX ,locRank ,temperatureVarIndex

    real(rk) :: tempNorm ,densNorm ,timeNorm

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 

    endState(1) = NamedInteger(tTag//"."//keyEndHState,1)

    call envObj%jsonCont%load(endState)
    call envObj%jsonCont%output(endState)

    temperatureVarName(1) = NamedString(tTag//"."//keyElectronTemperatureVar,"")

    call envObj%jsonCont%load(temperatureVarName)
    call envObj%jsonCont%output(temperatureVarName)

    if (assertions) call assert(envObj%externalVars%isVarNameRegistered(temperatureVarName(1)%value),temperatureVarName(1)%value//&
                                " not registered in environment wrapper")

    temperatureVarIndex = envObj%externalVars%getVarIndex(temperatureVarName(1)%value)

    tempNorm = normObj%getNormalizationValue(keyTempEVNorm)
    densNorm = normObj%getNormalizationValue(keyDensNorm)
    timeNorm = normObj%getNormalizationValue(keyTimeNorm)

    call addJanevRadRecombTransition(crmData,locNumX,endState(1)%value,temperatureVarIndex,tempNorm,densNorm,timeNorm)

end subroutine addJanevRadRecombTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addJanevCollExIonTransitionToCRMData(crmData,envObj,tTag,normObj)
    !! Add Janev data collisional excitation/ionization transition to CRM data object based on JSON file

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag
    class(Normalization)       ,intent(in)    :: normObj

    type(NamedInteger) ,dimension(1) :: startState, endState ,fixedEnergyIndex

    type(NamedString) ,dimension(1) :: distributionVarName
    
    integer(ik) :: locNumX ,locRank ,distVarIndex

    real(rk) :: tempNorm ,csNorm

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 

    startState(1) = NamedInteger(tTag//"."//keyStartHState,1)
    endState(1) = NamedInteger(tTag//"."//keyEndHState,1)
    fixedEnergyIndex(1) = NamedInteger(tTag//"."//keyFixedEnergyIndex,1)
    
    call envObj%jsonCont%load(startState)
    call envObj%jsonCont%output(startState)

    call envObj%jsonCont%load(endState)
    call envObj%jsonCont%output(endState)

    call envObj%jsonCont%load(fixedEnergyIndex)
    call envObj%jsonCont%output(fixedEnergyIndex)

    distributionVarName(1) = NamedString(tTag//"."//keyDistributionVarName,"")

    call envObj%jsonCont%load(distributionVarName)
    call envObj%jsonCont%output(distributionVarName)

    if (assertions) call assert(envObj%externalVars%isVarNameRegistered(distributionVarName(1)%value),&
                                distributionVarName(1)%value//" not registered in environment wrapper")

    distVarIndex = envObj%externalVars%getVarIndex(distributionVarName(1)%value)

    tempNorm = normObj%getNormalizationValue(keyTempEVNorm)
    csNorm = normObj%getNormalizationValue(keyCrossSectionNorm)

    call addJanevCollExIonTransition(crmData,locNumX,startState(1)%value,endState(1)%value,distVarIndex,&
                                    tempNorm,csNorm,envObj%vSpaceObj,fixedEnergyIndex(1)%value)

end subroutine addJanevCollExIonTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addJanevCollDeexRecombTransitionToCRMData(crmData,envObj,tTag,normObj)
    !! Add Janev data collisional deexcitation/recombination transition to CRM data object based on JSON file

    type(ModelboundCRMData)    ,intent(inout) :: crmData
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: tTag 
    class(Normalization)       ,intent(in)    :: normObj

    type(NamedInteger) ,dimension(1) :: startState, endState &
                                       ,directTransitionFixedEnergyIndex ,directTransitionIndex ,csUpdatePriority,fixedEnergyIndex

    type(NamedString) ,dimension(1) :: distributionVarName ,temperatureVarName
    
    integer(ik) :: locNumX ,locRank ,distVarIndex ,tempVarIndex

    real(rk) :: tempNorm ,densNorm

    locRank = envObj%mpiCont%getWorldRank()
    locNumX = envObj%partitionObj%getMaxXAtInd(locRank+1) - envObj%partitionObj%getMinXAtInd(locRank+1) + 1 

    startState(1) = NamedInteger(tTag//"."//keyStartHState,1)
    endState(1) = NamedInteger(tTag//"."//keyEndHState,1)
    directTransitionFixedEnergyIndex(1) = NamedInteger(tTag//"."//keyDirectTransitionEnergyIndex,1)
    directTransitionIndex(1) = NamedInteger(tTag//"."//keyDirectTransitionIndex,1)
    fixedEnergyIndex(1) = NamedInteger(tTag//"."//keyFixedEnergyIndex,1)
    csUpdatePriority(1) = NamedInteger(tTag//"."//keyCSUpdatePriority,0)
    
    call envObj%jsonCont%load(startState)
    call envObj%jsonCont%output(startState)

    call envObj%jsonCont%load(endState)
    call envObj%jsonCont%output(endState)

    call envObj%jsonCont%load(directTransitionFixedEnergyIndex)
    call envObj%jsonCont%output(directTransitionFixedEnergyIndex)

    call envObj%jsonCont%load(fixedEnergyIndex)
    call envObj%jsonCont%output(fixedEnergyIndex)

    call envObj%jsonCont%load(csUpdatePriority)
    call envObj%jsonCont%output(csUpdatePriority)

    call envObj%jsonCont%load(directTransitionIndex)
    call envObj%jsonCont%output(directTransitionIndex)

    distributionVarName(1) = NamedString(tTag//"."//keyDistributionVarName,"")
    temperatureVarName(1) = NamedString(tTag//"."//keyElectronTemperatureVar,"")

    call envObj%jsonCont%load(distributionVarName)
    call envObj%jsonCont%output(distributionVarName)

    call envObj%jsonCont%load(temperatureVarName)
    call envObj%jsonCont%output(temperatureVarName)

    if (assertions) then 
        call assert(envObj%externalVars%isVarNameRegistered(distributionVarName(1)%value),distributionVarName(1)%value//&
                                " not registered in environment wrapper")
        call assert(envObj%externalVars%isVarNameRegistered(temperatureVarName(1)%value),temperatureVarName(1)%value//&
                                " not registered in environment wrapper")
    end if

    distVarIndex = envObj%externalVars%getVarIndex(distributionVarName(1)%value)
    tempVarIndex = envObj%externalVars%getVarIndex(temperatureVarName(1)%value)

    tempNorm = normObj%getNormalizationValue(keyTempEVNorm)
    densNorm = normObj%getNormalizationValue(keyDensNorm)

    call addJanevCollDeexRecombTransition(crmData,locNumX,startState(1)%value,endState(1)%value,distVarIndex,tempVarIndex,&
                                          tempNorm,densNorm,envObj%vSpaceObj,directTransitionFixedEnergyIndex(1)%value,&
                                          fixedEnergyIndex(1)%value,&
                                          directTransitionIndex(1)%value,csUpdatePriority(1)%value)

end subroutine addJanevCollDeexRecombTransitionToCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule modelbound_data_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
