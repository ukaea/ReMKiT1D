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
submodule (manipulator_support) manipulator_support_procedures
!! author: Stefan Mijin 
!! 
!! Contains implementation of manipulator initialization support procedures

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCompositeManipulatorFromJSON(manip,envObj,normObj)
    !! Initialize modelbound data and add to corresponding model object

    type(CompositeManipulator) ,allocatable ,intent(inout) :: manip
    type(EnvironmentWrapper)                ,intent(inout) :: envObj
    class(Normalization)                    ,intent(in)    :: normObj

    type(NamedStringArray)     ,dimension(1) :: manipulatorTags
    type(NamedString)          ,dimension(1) :: manipulatorType

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(envObj%isDefined(),"Undefined environment wrapper passed to initCompositeManipulatorFromJSON")
        call assert(normObj%isDefined(),"Undefined normalization object passed to initCompositeManipulatorFromJSON")
    end if
    
    manipulatorTags(1)%name = keyManipulators//"."//keyTags
    allocate(manipulatorTags(1)%values(0))

    call envObj%jsonCont%load(manipulatorTags)

    if (size(manipulatorTags(1)%values) > 0) then 
        allocate(manip)
        call manip%init(size(manipulatorTags(1)%values))
    end if

    do i = 1,size(manipulatorTags(1)%values)
        manipulatorType(1) = NamedString(keyManipulators//"."//manipulatorTags(1)%values(i)%string//"."//keyType,"")

        call envObj%jsonCont%load(manipulatorType)

        select case(manipulatorType(1)%value)
        case (keyGroupEvaluator)
            call addGroupEvaluatorToCompositeManipulator(manip,envObj,normObj,&
                                                         keyManipulators//"."//manipulatorTags(1)%values(i)%string)
        case (keyTermEvaluator)
            call addTermEvaluatorToCompositeManipulator(manip,envObj,normObj,&
                                                            keyManipulators//"."//manipulatorTags(1)%values(i)%string)
        case (keyExtractor)
            call addMBDataExtractorToCompositeManipulator(manip,envObj,normObj,&
                                                          keyManipulators//"."//manipulatorTags(1)%values(i)%string)
        case default 
            error stop "Unsupported manipulator type detected by initCompositeManipulatorFromJSON"
        end select
    end do

end subroutine initCompositeManipulatorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addGroupEvaluatorToCompositeManipulator(manip,envObj,normObj,jsonPrefix)
    !! Add GroupEvaluator type manipulator based json file

    type(CompositeManipulator)              ,intent(inout) :: manip
    type(EnvironmentWrapper)                ,intent(inout) :: envObj
    class(Normalization)                    ,intent(in)    :: normObj
    character(*)                            ,intent(in)    :: jsonPrefix

    type(GroupEvaluator) :: evalManip 

    type(NamedInteger) ,dimension(2) :: namedIntParams 
    type(NamedString)  ,dimension(2) :: namedStringParams 
    type(NamedStringArray) ,dimension(1) :: modelTags

    integer(ik) :: i ,modelIndex ,varIndex

    logical :: modelFound

    namedStringParams(1) = NamedString(jsonPrefix//"."//keyResultVarName,"")

    namedStringParams(2) = NamedString(jsonPrefix//"."//keyModelTag,"")

    namedIntParams(1) = NamedInteger(jsonPrefix//"."//keyEvaluatedGroup,1)
    namedIntParams(2) = NamedInteger(jsonPrefix//"."//keyPriority,0)

    modelTags(1)%name = keyModels//"."//keyTags
    allocate(modelTags(1)%values(0))

    call envObj%jsonCont%load(modelTags)
    call envObj%jsonCont%load(namedStringParams)
    call envObj%jsonCont%load(namedIntParams)

    call envObj%jsonCont%output(namedStringParams)
    call envObj%jsonCont%output(namedIntParams)

    modelFound = .false. 
    do i = 1, size(modelTags(1)%values)
        if (modelTags(1)%values(i)%string == namedStringParams(2)%value) then 
            modelIndex = i 
            modelFound = .true. 
            exit 
        end if
    end do

    if (.not. modelFound) error stop "model requested by GroupEvaluator not found"

    if (assertions .or. assertionLvl >= 0) call assert(envObj%externalVars%isVarNameRegistered(namedStringParams(1)%value),&
    namedStringParams(1)%value//&
    " variable requested by GroupEvaluator not found in environment wrapper")

    varIndex = envObj%externalVars%getVarIndex(namedStringParams(1)%value)

    call evalManip%init(varIndex,modelIndex,namedIntParams(1)%value)

    call manip%addManipulator(evalManip,namedIntParams(2)%value)
    
end subroutine addGroupEvaluatorToCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addTermEvaluatorToCompositeManipulator(manip,envObj,normObj,jsonPrefix)
    !! Add TermEvaluator type manipulator based json file

    type(CompositeManipulator)              ,intent(inout) :: manip
    type(EnvironmentWrapper)                ,intent(inout) :: envObj
    class(Normalization)                    ,intent(in)    :: normObj
    character(*)                            ,intent(in)    :: jsonPrefix

    type(TermEvaluator) :: evalManip 

    type(NamedInteger) ,dimension(1) :: manipPriority 
    type(NamedString)  ,dimension(1) :: resultVarName 
    type(NamedStringArray) ,dimension(1) :: modelEvalTags, termNames 
    type(NamedStringArray) ,dimension(1) :: modelTags
    type(NamedLogical) ,dimension(2) :: evalOptions

    integer(ik) :: i ,j ,varIndex

    integer(ik) ,allocatable ,dimension(:) :: modelIndices
    type(StringArray) ,allocatable ,dimension(:) :: termNameStrings

    logical :: modelFound

    modelTags(1)%name = keyModels//"."//keyTags
    allocate(modelTags(1)%values(0))

    resultVarName(1) = NamedString(jsonPrefix//"."//keyResultVarName,"")
    modelEvalTags(1)%name = jsonPrefix//"."//keyEvaluatedModelNames
    allocate(modelEvalTags(1)%values(0))
    termNames(1)%name = jsonPrefix//"."//keyEvaluatedTermNames
    allocate(termNames(1)%values(0))
    
    manipPriority(1) = NamedInteger(jsonPrefix//"."//keyPriority,0)

    evalOptions(1) = NamedLogical(jsonPrefix//"."//keyUpdate,.false.)
    evalOptions(2) = NamedLogical(jsonPrefix//"."//keyAccumulate,.false.)

    call envObj%jsonCont%load(modelTags)
    call envObj%jsonCont%load(modelEvalTags)
    call envObj%jsonCont%load(termNames)
    call envObj%jsonCont%load(resultVarName)
    call envObj%jsonCont%load(manipPriority)
    call envObj%jsonCont%load(evalOptions)

    call envObj%jsonCont%output(modelEvalTags)
    call envObj%jsonCont%output(termNames)
    call envObj%jsonCont%output(resultVarName)
    call envObj%jsonCont%output(manipPriority)
    call envObj%jsonCont%output(evalOptions)

    if (assertions .or. assertionLvl >= 0) call assert(size(termNames(1)%values) == size(modelEvalTags(1)%values),&
                                     modelEvalTags(1)%name//" and "//termNames(1)%name//" must be of the same size")
    
    allocate(modelIndices(size(modelEvalTags(1)%values)))
    allocate(termNameStrings(size(modelEvalTags(1)%values)))
    do i = 1,size(modelEvalTags(1)%values)
        modelFound = .false.
        do j = 1, size(modelTags(1)%values)
            if (modelTags(1)%values(j)%string == modelEvalTags(1)%values(i)%string) then 
                modelIndices(i) = j 
                modelFound = .true. 
                exit 
            end if
        end do
        if (.not. modelFound) error stop "model requested by TermEvaluator not found"

        termNameStrings(i)%string = termNames(1)%values(i)%string
    end do

    if (assertions .or. assertionLvl >= 0) call assert(envObj%externalVars%isVarNameRegistered(resultVarName(1)%value),&
    resultVarName(1)%value//&
    " variable requested by TermEvaluator not found in environment wrapper")

    varIndex = envObj%externalVars%getVarIndex(resultVarName(1)%value)

    call evalManip%init(varIndex,modelIndices,termNameStrings,&
        accumulate = evalOptions(1)%value, update = evalOptions(2)%value)

    call manip%addManipulator(evalManip,manipPriority(1)%value)
    
end subroutine addTermEvaluatorToCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addMBDataExtractorToCompositeManipulator(manip,envObj,normObj,jsonPrefix)
    !! Add ModelboundDataExtractor type manipulator based json file

    type(CompositeManipulator)              ,intent(inout) :: manip
    type(EnvironmentWrapper)                ,intent(inout) :: envObj
    class(Normalization)                    ,intent(in)    :: normObj
    character(*)                            ,intent(in)    :: jsonPrefix

    type(ModelboundDataExtractor) :: extrManip 

    type(NamedInteger) ,dimension(1) :: namedIntParams 
    type(NamedString)  ,dimension(3) :: namedStringParams 
    type(NamedStringArray) ,dimension(1) :: modelTags

    integer(ik) :: i ,modelIndex ,varIndex

    logical :: modelFound

    namedStringParams(1) = NamedString(jsonPrefix//"."//keyResultVarName,"")

    namedStringParams(2) = NamedString(jsonPrefix//"."//keyModelTag,"")

    namedStringParams(3) = NamedString(jsonPrefix//"."//keyModelboundDataName,"")

    namedIntParams(1) = NamedInteger(jsonPrefix//"."//keyPriority,0)

    modelTags(1)%name = keyModels//"."//keyTags
    allocate(modelTags(1)%values(0))

    call envObj%jsonCont%load(modelTags)
    call envObj%jsonCont%load(namedStringParams)
    call envObj%jsonCont%load(namedIntParams)

    call envObj%jsonCont%output(namedStringParams)
    call envObj%jsonCont%output(namedIntParams)

    modelFound = .false. 
    do i = 1, size(modelTags(1)%values)
        if (modelTags(1)%values(i)%string == namedStringParams(2)%value) then 
            modelIndex = i 
            modelFound = .true. 
            exit 
        end if
    end do

    if (.not. modelFound) error stop "model requested by ModelboundDataExtractor not found"

    if (assertions .or. assertionLvl >= 0) call assert(envObj%externalVars%isVarNameRegistered(namedStringParams(1)%value),&
    namedStringParams(1)%value//&
    " variable requested by ModelboundDataExtractor not found in environment wrapper")

    varIndex = envObj%externalVars%getVarIndex(namedStringParams(1)%value)

    call extrManip%init(varIndex,modelIndex,namedStringParams(3)%value)

    call manip%addManipulator(extrManip,namedIntParams(1)%value)
    
end subroutine addMBDataExtractorToCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule manipulator_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
