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
submodule (custom_model_builder_class) custom_model_builder_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the CustomModelBuilder class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCustomBuilder(this,env,normObject,modelTag) 
    !! Constructs the model for this builder and sets it into defined state

    class(CustomModelBuilder)         ,intent(inout)  :: this 
    class(EnvironmentWrapper)         ,intent(inout)  :: env 
    class(Normalization)              ,intent(in)     :: normObject !! Reference normalization object 
    character(*)                      ,intent(in)     :: modelTag !! Tag of this model

    integer(ik) :: currentProc

    logical :: modelActive ! True if fluid variable are evolved on current processor

    type(NamedScalarContainer) :: scalarParams 
    type(NamedArrayContainer)  :: arrayParams
    integer(ik) :: i 

    class(ModelboundData) ,allocatable :: mbData

    if (assertions) then 

        call assert(.not. this%isDefined(),"CustomModelBuilder constructor called for already defined object")
        call assert(env%isDefined(),"Undefined environment wrapper passed to CustomModelBuilder constructor")

    end if

    currentProc = env%mpiCont%getWorldRank()

    modelActive = env%partitionObj%getMinHAtInd(currentProc+1) == 1

    ! Set up default parameters
    allocate(scalarParams%intData(2))
    allocate(scalarParams%realData(0))
    allocate(scalarParams%logicalData(0))
    allocate(scalarParams%stringData(1))

    allocate(arrayParams%intData(0))
    allocate(arrayParams%realData(0))
    allocate(arrayParams%logicalData(0))
    allocate(arrayParams%stringData(1))

    scalarParams%intData(1) = NamedInteger(keyIntegrator//"."//keyNumImplicitGroups,1)
    scalarParams%intData(2) = NamedInteger(keyIntegrator//"."//keyNumGenGroups,1)
    scalarParams%stringData(1) = NamedString(keyModels//"."//modelTag//"."//keyModelboundData//"."//keyModelboundDataType,keyNone)

    arrayParams%stringData(1)%name = keyModels//"."//modelTag//"."//keyTermTags

    allocate(arrayParams%stringData(1)%values(0))

    call this%setScalarParams(scalarParams)
    call this%setArrayParams(arrayParams)

    call this%loadParams(env%jsonCont,env%mpiCont)
    scalarParams = this%getScalarParams()
    arrayParams = this%getArrayParams()

    call this%makeDefined()
    allocate(this%modelBuffer)

    call this%outputUsedParams(env%jsonCont,env%mpiCont)

    call this%modelBuffer%init(numGeneralTerms=0,numImplicitGroups=scalarParams%intData(1)%value,&
                               numGeneralGroups=scalarParams%intData(2)%value)

    ! Add modelbound data

    if (scalarParams%stringData(1)%value /= keyNone) then 
        call addModelboundDataToModel(this%modelBuffer,modelTag,env,normObject)
        call this%modelBuffer%copyModelData(mbData)
    end if

    if (allocated(mbData)) then
        ! Add any generated terms and finish modelBuffer initialization
        call this%applyTermGenerator(env,normObject,modelTag,size(arrayParams%stringData(1)%values),mbData=mbData)
        ! Adding each custom term
        do i = 1, size(arrayParams%stringData(1)%values)
            call printMessage("Adding term: "//arrayParams%stringData(1)%values(i)%string//" to "//modelTag)
            call this%addTermToModel(keyModels//"."//modelTag,arrayParams%stringData(1)%values(i)%string,env,normObject&
                                     ,mbData)
        end do
    else
        ! Add any generated terms and finish modelBuffer initialization
        call this%applyTermGenerator(env,normObject,modelTag,size(arrayParams%stringData(1)%values))
        ! Adding each custom term
        do i = 1, size(arrayParams%stringData(1)%values)
            call printMessage("Adding term: "//arrayParams%stringData(1)%values(i)%string//" to "//modelTag)
            call this%addTermToModel(keyModels//"."//modelTag,arrayParams%stringData(1)%values(i)%string,env,normObject)
        end do
    end if

end subroutine initCustomBuilder
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addCustomModel(this,modellerObj) 
    !! Adds the model built by the builder and resets the builder to become undefined for further use

    class(CustomModelBuilder)           ,intent(inout) :: this 
    class(Modeller)                     ,intent(inout) :: modellerObj

    if (assertions) then 

        call assert(this%isDefined(),"Attempted to add advection model using undefined builder")
        call assert(modellerObj%isDefined(),"Attempted to add advection model to undefined modeller object")

    end if

    call modellerObj%addModel(this%modelBuffer)
    call this%makeUndefined()

end subroutine addCustomModel
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addTermToModel(this,termJSONPrefix,termTag,env,normObject,mbData) 
    !! Adds individual term to model buffer based on JSON file data

    class(CustomModelBuilder)                 ,intent(inout)  :: this 
    character(*)                              ,intent(in)     :: termJSONPrefix 
    character(*)                              ,intent(in)     :: termTag 
    class(EnvironmentWrapper)                 ,intent(inout)  :: env 
    class(Normalization)                      ,intent(in)     :: normObject
    class(ModelboundData) ,optional           ,intent(in)     :: mbData

    class(MatrixTerm) ,allocatable :: termBuffer

    type(GeneralMatrixTerm)             :: genTerm
    type(StencilTemplate)               :: stencilData
    type(CoordProfiles)                 :: coordProfs

    type(VarData)  :: vData
    type(TimeSignalData) :: tData

    type(NamedString)       ,allocatable ,dimension(:) :: stringParams ! evolvedVar, implicitVar ,tSignalName

    type(NamedStringArray)  ,allocatable ,dimension(:) :: stringArrayParams ! requiredRowVars ,requiredColVars, requiredRowMBVars ,requiredColMBVars ,normVars
    type(NamedRealArray)    ,allocatable ,dimension(:) :: realArrayParams ! rowVarPowers ,colVarPowers ,rowMBVarPowers ,colMBVarPowers ,normVarPowers ,xProfile, hProfile,vProfile ,tParams
    type(NamedReal)         ,allocatable ,dimension(:) :: realParams !normMultConst , tPeriod 
    type(NamedInteger)      ,allocatable ,dimension(:) :: intParams !evaluatedGroupIndex 
    type(NamedIntegerArray) ,allocatable ,dimension(:) :: intArrayParams !implicitTermGroups ,generalTermGroups
    type(NamedLogical)      ,allocatable ,dimension(:) :: logicalParams ! realTimePeriod

    integer(ik) :: i ,j

    real(rk) :: normConst ,timeNorm

    logical :: pGrid ,kineticTerm

    character(:) ,allocatable  :: copyTermNameBuffer
    
    allocate(stringParams(4))

    stringParams(1) = NamedString(termJSONPrefix//"."//termTag//"."//keyEvolvedVar,"")
    stringParams(2) = NamedString(termJSONPrefix//"."//termTag//"."//keyImplicitVar,"")
    stringParams(3) = NamedString(termJSONPrefix//"."//termTag//"."//keyTimeSignalData//"."//keySignal,keyNone)
    stringParams(4) = NamedString(termJSONPrefix//"."//termTag//"."//keyMultCopyTermName,keyNone)

    allocate(stringArrayParams(5))

    stringArrayParams(1)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqRowVarNames
    allocate(stringArrayParams(1)%values(0))
    stringArrayParams(2)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqColVarNames
    allocate(stringArrayParams(2)%values(0))
    stringArrayParams(3)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqMBRowVarNames
    allocate(stringArrayParams(3)%values(0))
    stringArrayParams(4)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqMBColVarNames
    allocate(stringArrayParams(4)%values(0))
    stringArrayParams(5)%name = termJSONPrefix//"."//termTag//"."//keyCustomNormConst//"."//keyNormNames
    allocate(stringArrayParams(5)%values(0))

    allocate(realArrayParams(9))

    realArrayParams(1)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqRowVarPowers
    allocate(realArrayParams(1)%values(0))
    realArrayParams(2)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqColVarPowers
    allocate(realArrayParams(2)%values(0))
    realArrayParams(3)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqMBRowVarPowers
    allocate(realArrayParams(3)%values(0))
    realArrayParams(4)%name = termJSONPrefix//"."//termTag//"."//keyVarData//"."//keyReqMBColVarPowers
    allocate(realArrayParams(4)%values(0))
    realArrayParams(5)%name = termJSONPrefix//"."//termTag//"."//keyCustomNormConst//"."//keyNormPowers
    allocate(realArrayParams(5)%values(0))
    realArrayParams(6)%name = termJSONPrefix//"."//termTag//"."//keySpatialProfile
    allocate(realArrayParams(6)%values(0))
    realArrayParams(7)%name = termJSONPrefix//"."//termTag//"."//keyHarmonicProfile
    allocate(realArrayParams(7)%values(0))
    realArrayParams(8)%name = termJSONPrefix//"."//termTag//"."//keyVelocityProfile
    allocate(realArrayParams(8)%values(0))
    realArrayParams(9)%name = termJSONPrefix//"."//termTag//"."//keyTimeSignalData//"."//keyTimeSignalParams
    allocate(realArrayParams(9)%values(0))

    allocate(realParams(2))

    realParams(1) = NamedReal(termJSONPrefix//"."//termTag//"."//keyCustomNormConst//"."//keyMultConst,real(1,kind=rk))
    realParams(2) = NamedReal(termJSONPrefix//"."//termTag//"."//keyTimeSignalData//"."//keyTimeSignalPeriod,real(1,kind=rk))

    allocate(intParams(1))

    intParams(1) = NamedInteger(termJSONPrefix//"."//termTag//"."//keyEvaluatedGroup,0)

    allocate(intArrayParams(2))

    intArrayParams(1) = NamedIntegerArray(termJSONPrefix//"."//termTag//"."//keyImplicitTermGroups,[1])
    intArrayParams(2) = NamedIntegerArray(termJSONPrefix//"."//termTag//"."//keyGeneralTermGroups,[1])

    allocate(logicalParams(3))

    logicalParams(1) = NamedLogical(termJSONPrefix//"."//termTag//"."//keyTimeSignalData//"."//keyRealTimePeriod,.false.)
    logicalParams(2) = NamedLogical(termJSONPrefix//"."//termTag//"."//keySkipPattern,.false.)
    logicalParams(3) = NamedLogical(termJSONPrefix//"."//termTag//"."//keyFixedMatrix,.false.)

    call env%jsonCont%load(stringParams)
    call env%jsonCont%load(stringArrayParams)
    call env%jsonCont%load(realArrayParams)
    call env%jsonCont%load(realParams)
    call env%jsonCont%load(intParams)
    call env%jsonCont%load(intArrayParams)
    call env%jsonCont%load(logicalParams)

    call env%jsonCont%output(stringParams)
    call env%jsonCont%output(stringArrayParams)
    call env%jsonCont%output(realArrayParams)
    call env%jsonCont%output(realParams)
    call env%jsonCont%output(intParams)
    call env%jsonCont%output(intArrayParams)
    call env%jsonCont%output(logicalParams)

    if (size(realArrayParams(6)%values) == 0) realArrayParams(6)%values = [(real(1,kind=rk),i=1,env%gridObj%getNumX())]
    if (size(realArrayParams(7)%values) == 0) realArrayParams(7)%values = [(real(1,kind=rk),i=1,env%gridObj%getNumH())]
    if (size(realArrayParams(8)%values) == 0) realArrayParams(8)%values = [(real(1,kind=rk),i=1,env%gridObj%getNumV())]


    do i = 1,5

        if (size(stringArrayParams(i)%values) > 0 .and. size(realArrayParams(i)%values)==0)&
        realArrayParams(i)%values = [(real(1,kind=rk),j=1,size(stringArrayParams(i)%values))]

    end do

    if (assertions) then 

        call assert(env%externalVars%isVarNameRegistered(stringParams(1)%value),stringParams(1)%name//&
        " not registered in environment wrapper")

        call assert(env%externalVars%isVarNameRegistered(stringParams(2)%value),stringParams(2)%name//&
        " not registered in environment wrapper")

        do i = 1,4
            call assert(size(stringArrayParams(i)%values) == size(realarrayParams(i)%values),stringArrayParams(i)%name//" and "//&
            realArrayParams(i)%name//" must be of the same size")

            if (i<3) then
                do j = 1,size(stringArrayParams(i)%values)
                    call assert(env%externalVars%isVarNameRegistered(stringArrayParams(i)%values(j)%string),&
                    "Variable with name "//stringArrayParams(i)%values(j)%string//" not registered in environment wrapper")
                end do
            end if
        end do

        call assert(size(stringArrayParams(5)%values) == size(realarrayParams(5)%values),stringArrayParams(5)%name//" and "//&
            realArrayParams(5)%name//" must be of the same size")

        call assert(size(realArrayParams(6)%values) == env%gridObj%getNumX(),"Spatial profile "//realArrayParams(6)%name//&
        " does not conform to grid size")

        call assert(size(realArrayParams(7)%values) == env%gridObj%getNumH(),"Harmonic profile "//realArrayParams(7)%name//&
        " does not conform to grid size")

        call assert(size(realArrayParams(8)%values) == env%gridObj%getNumV(),"Velocity profile "//realArrayParams(8)%name//&
        " does not conform to grid size")

        if (stringParams(3)%value /= keyNone) call assert(env%signalCollectionObj%isSignalRegistered(stringParams(3)%value),&
        stringParams(3)%name//" not registered in the signal collection object in environment wrapper")
    end if

    vData%rowVars = stringArrayParams(1)%values
    vData%colVars = stringArrayParams(2)%values
    vData%modelboundRowVars = stringArrayParams(3)%values
    vData%modelboundColVars = stringArrayParams(4)%values

    vData%rowVarPowers = realArrayParams(1)%values
    vData%colVarPowers = realArrayParams(2)%values
    vData%modelboundRowVarPowers = realArrayParams(3)%values
    vData%modelboundColVarPowers = realArrayParams(4)%values

    normConst = normObject%getCustomNormalization(stringArrayParams(5)%values,realArrayParams(5)%values,realParams(1)%value)

    pGrid = env%geometryObj%isPeriodic()

    kineticTerm = env%externalVars%isVarDist(env%externalVars%getVarIndex(stringParams(1)%value)) &
             .or. env%externalVars%isVarDist(env%externalVars%getVarIndex(stringParams(2)%value))

    if (kineticTerm) then 
        if (present(mbData)) then
            call initKineticStencilTemplate(stencilData,env,termJSONPrefix//"."//termTag,&
            stringParams(1)%value,stringParams(2)%value,mbData)
        else
            call initKineticStencilTemplate(stencilData,env,termJSONPrefix//"."//termTag,&
            stringParams(1)%value,stringParams(2)%value)
        end if
    else
        call initFluidStencilTemplate(stencilData,env,termJSONPrefix//"."//termTag,stringParams(1)%value,stringParams(2)%value)
    end if

    coordProfs%xProfile = realArrayParams(6)%values
    coordProfs%hProfile = realArrayParams(7)%values
    coordProfs%vProfile = realArrayParams(8)%values

    if (stringParams(4)%value /= keyNone) copyTermNameBuffer = stringParams(4)%value

    if (stringParams(3)%value /= keyNone) then 

        timeNorm = normObject%getNormalizationValue(keyTimeNorm)

        tData%tParams = realArrayParams(9)%values
        tData%tPeriod = realParams(2)%value 
        if (logicalParams(1)%value) tData%tPeriod = tData%tPeriod/timeNorm
        call env%signalCollectionObj%copySignal(stringParams(3)%value,tData%tSignal)

        call genTerm%init(env%gridObj,env%partitionObj,env%indexingObj,env%mpiCont%getWorldRank(),&
        stringParams(1)%value ,stringParams(2)%value,env%externalVars,stencilData,&
        vData=vData,normConst=normConst,tData=tData,coordProfile=coordProfs,evalTermGroup=intParams(1)%value,&
        copyTermName=copyTermNameBuffer)

    else

        call genTerm%init(env%gridObj,env%partitionObj,env%indexingObj,env%mpiCont%getWorldRank(),&
        stringParams(1)%value ,stringParams(2)%value,env%externalVars,stencilData,&
        vData=vData,normConst=normConst,coordProfile=coordProfs,evalTermGroup=intParams(1)%value,&
        copyTermName=copyTermNameBuffer)

    end if

    if (logicalParams(3)%value) call genTerm%setFixedMatrix(.true.)
    allocate(termBuffer,source=genTerm)
    call this%modelBuffer%addImplicitTerm(termBuffer,intArrayParams(1)%values,intArrayParams(2)%values,termTag,&
                                          skipPattern=logicalParams(2)%value)

end subroutine addTermToModel
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine applyTermGenerator(this,env,normObject,modelTag,currentNumTerms,mbData) 
    !! Checks for associated term generator in JSON file and applies to custom model

    class(CustomModelBuilder)                 ,intent(inout)  :: this 
    class(EnvironmentWrapper)                 ,intent(inout)  :: env 
    class(Normalization)                      ,intent(in)     :: normObject !! Reference normalization object
    character(*)                              ,intent(in)     :: modelTag !! Tag of this model
    integer(ik)                               ,intent(in)     :: currentNumTerms
    class(ModelboundData) ,optional           ,intent(in)     :: mbData

    type(NamedString) ,allocatable ,dimension(:) :: generatorType
    type(NamedStringArray) ,dimension(1) :: generatorTags

    type(NamedIntegerArray) ,allocatable ,dimension(:) :: generatedTermImplicitGroups

    integer(ik) :: numGeneratedImpTerms ,i ,j

    class(TermGeneratorContainer) ,allocatable ,dimension(:) :: tGenerators 

    type(IntArray) ,allocatable ,dimension(:) :: iTermIGroups

    generatorTags(1)%name = keyModels//"."//modelTag//"."//keyTermGenerators//"."//keyTags
    allocate(generatorTags(1)%values(0))

    call env%jsonCont%load(generatorTags)
    call env%jsonCont%output(generatorTags)

    allocate(generatorType(size(generatorTags(1)%values)))
    allocate(generatedTermImplicitGroups(size(generatorTags(1)%values)))

    do i = 1,size(generatorTags(1)%values)
        generatorType(i) = NamedString(keyModels//"."//modelTag//"."//keyTermGenerators//&
                                       "."//generatorTags(1)%values(i)%string//"."//keyType,keyNone)

        generatedTermImplicitGroups(i) = NamedIntegerArray(keyModels//"."//modelTag//"."//keyTermGenerators//&
                                                           "."//generatorTags(1)%values(i)%string//"."//keyImplicitTermGroups,[1])
    end do

    call env%jsonCont%load(generatorType)
    call env%jsonCont%output(generatorType)
    call env%jsonCont%load(generatedTermImplicitGroups)
    call env%jsonCont%output(generatedTermImplicitGroups)

    allocate(tGenerators(size(generatorTags(1)%values)))

    numGeneratedImpTerms = 0 

    do i = 1,size(generatorTags(1)%values)
        if (generatorType(i)%value /= keyNone) then 

            call printMessage("Applying term generator: "//generatorTags(1)%values(i)%string//" to "//modelTag)
            !Add supported generator calls here
            select case (generatorType(i)%value)
            case (keyCRMDensTermGen)
                call initCRMDensTermGeneratorFromJSON(tGenerators(i)%entry,this%modelBuffer,env,&
                                                      keyModels//"."//modelTag//"."//keyTermGenerators//&
                                                      "."//generatorTags(1)%values(i)%string,generatorTags(1)%values(i)%string)
            case (keyCRMElEnergyTermGen)
                call initCRMElEnergyTermGeneratorFromJSON(tGenerators(i)%entry,this%modelBuffer,env,&
                                                        keyModels//"."//modelTag//"."//keyTermGenerators//&
                                                        "."//generatorTags(1)%values(i)%string,generatorTags(1)%values(i)%string)
            case (keyCRMBoltzTermGen)
                call initCRMFixedBoltzTermGeneratorFromJSON(tGenerators(i)%entry,normObject,this%modelBuffer,env,&
                                                        keyModels//"."//modelTag//"."//keyTermGenerators//&
                                                        "."//generatorTags(1)%values(i)%string,generatorTags(1)%values(i)%string)
            case (keyCRMVarBoltzTermGen)
                call initCRMVarBoltzTermGeneratorFromJSON(tGenerators(i)%entry,normObject,this%modelBuffer,env,&
                                                        keyModels//"."//modelTag//"."//keyTermGenerators//&
                                                        "."//generatorTags(1)%values(i)%string,generatorTags(1)%values(i)%string)
            case (keyCRMSecElTermGen)
                call initCRMSecElTermGeneratorFromJSON(tGenerators(i)%entry,this%modelBuffer,env,&
                                                        keyModels//"."//modelTag//"."//keyTermGenerators//&
                                                        "."//generatorTags(1)%values(i)%string,generatorTags(1)%values(i)%string)
            case default 
                error stop "applyTermGenerator in CustomModelBuilder detected unsupported term generator type"
            end select
            
            if (present(mbData)) then 
                call tGenerators(i)%entry%generate(mbData)
            else
                call tGenerators(i)%entry%generate()
            end if
            numGeneratedImpTerms = numGeneratedImpTerms + tGenerators(i)%entry%getNumImplicitTerms()

        end if

    end do

    call this%modelBuffer%setNumImplicitTerms(numGeneratedImpTerms+currentNumTerms)

    do i = 1,size(generatorTags(1)%values)
        if (generatorType(i)%value /= keyNone) then 

            if (tGenerators(i)%entry%getNumImplicitTerms() > 0) then 
                if (allocated(iTermIGroups)) deallocate(iTermIGroups)
                allocate(iTermIGroups(tGenerators(i)%entry%getNumImplicitTerms()))
                do j = 1, tGenerators(i)%entry%getNumImplicitTerms()
                    iTermIGroups(j)%entry = generatedTermImplicitGroups(i)%values
                end do
                call tGenerators(i)%entry%moveTerms(this%modelBuffer,impTermImpGroups=iTermIGroups)
            end if

        end if
    end do
end subroutine applyTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule custom_model_builder_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
