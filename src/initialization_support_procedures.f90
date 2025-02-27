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
submodule (initialization_support) initialization_support_procedures
!! author: Stefan Mijin  
!!
!! Contains the implementations of various initialization support routines

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVarListFromJSON(varList,jsonCont,mpiCont,isDerivedList)
    !! Initialize variable list based on variable data from a JSON file. If isDerivedList = .true. the list is initialized
    !! based on derived list data, otherwise it uses the implicit list.

    type(VariableList)      ,intent(inout) :: varList
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
    logical ,optional       ,intent(in)    :: isDerivedList

    type(NamedStringArray) ,dimension(1) :: varNames 
    integer(ik) :: i
    logical :: derivList

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. varList%isDefined(),"initVarListFromJSON called on already defined variable list")
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to initVarListFromJSON")
    end if

    derivList = .false. 
    if (present(isDerivedList)) derivList = isDerivedList 

    varNames(1)%name = keyVariables//"."//keyImplicitVars//"."//keyNames
    if (derivList) varNames(1)%name = keyVariables//"."//keyDerivedVars//"."//keyNames

    allocate(varNames(1)%values(0))

    call jsonCont%load(varNames)
    call jsonCont%output(varNames)

    call varList%init()

    do i = 1,size(varNames(1)%values)
        call addVarToListFromJSON(varList,jsonCont,mpiCont,derivList,varNames(1)%values(i)%string)
    end do

end subroutine initVarListFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine addVarToListFromJSON(varList,jsonCont,mpiCont,isDerivedList,varName)
    !! Add variable to variable list based on varName 

    type(VariableList)      ,intent(inout) :: varList
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
    logical                 ,intent(in)    :: isDerivedList
    character(*)            ,intent(in)    :: varName

    type(NamedLogical) ,allocatable ,dimension(:) :: logicalParams
    type(NamedInteger) ,dimension(1) :: varPriority
    character(:) ,allocatable :: jsonPrefix

    if (assertions .or. assertionLvl >= 0) then 
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to addVarToListFromJSON")
    end if

    jsonPrefix = keyVariables//"."//keyImplicitVars//"."//varName
    if (isDerivedList) jsonPrefix = keyVariables//"."//keyDerivedVars//"."//varName
    allocate(logicalParams(5))

    logicalParams(1) = NamedLogical(jsonPrefix//"."//keyIsDist,.false.)
    logicalParams(2) = NamedLogical(jsonPrefix//"."//keyIsSingleHarmonic,.false.)
    logicalParams(3) = NamedLogical(jsonPrefix//"."//keyIsScalar,.false.)
    logicalParams(4) = NamedLogical(jsonPrefix//"."//keyIsOnDualGrid,.false.)
    logicalParams(5) = NamedLogical(jsonPrefix//"."//keyIsStationary,.false.)

    varPriority(1) = NamedInteger(jsonPrefix//"."//keyPriority,0)
    
    call jsonCont%load(logicalParams)
    call jsonCont%load(varPriority)

    call jsonCont%output(logicalParams)
    call jsonCont%output(varPriority)

    call varList%addVar(varName,isDist=logicalParams(1)%value,isSingleHarmonic=logicalParams(2)%value,&
                                isScalar=logicalParams(3)%value,isOnDualGrid=logicalParams(4)%value,&
                                isStationary=logicalParams(5)%value, priority=varPriority(1)%value)

end subroutine addVarToListFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initGridFromJSON(gridObj,jsonCont,mpiCont,lengthNorm)
    !! Initialize grid object based on data from a JSON file. 

    type(Grid)              ,intent(inout) :: gridObj
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
    real(rk) ,optional      ,intent(in)    :: lengthNorm  !! Length normalization if the supplied grid is in meters, defaults to 1

    type(NamedRealArray) ,allocatable ,dimension(:) :: gridPoints ,cellWidth
    type(NamedInteger) ,allocatable ,dimension(:) :: intParams
    type(NamedLogical) ,dimension(2) :: constructGridFromWidths
    type(NamedLogical) ,dimension(1) :: gridInMeters

    real(rk) ,allocatable ,dimension(:) :: usedVGrid ,usedXGrid 
    
    real(rk) :: usedNorm 

    integer(ik) :: i 

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. gridObj%isDefined(),"initGridFromJSON called on already defined grid object")
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to initGridFromJSON")
    end if
    
    allocate(gridPoints(2))

    usedNorm = real(1,kind=rk)
    if (present(lengthNorm)) usedNorm = lengthNorm

    gridPoints(1)%name = keyXGrid//"."//keyCellCoords
    gridPoints(2)%name = keyVGrid//"."//keyCellCoords

    allocate(gridPoints(1)%values(0))
    allocate(gridPoints(2)%values(0))

    allocate(cellWidth(2))

    cellWidth(1)%name = keyXGrid//"."//keyCellWidths
    cellWidth(2)%name = keyVGrid//"."//keyCellWidths

    allocate(cellWidth(1)%values(0))
    allocate(cellWidth(2)%values(0))


    constructGridFromWidths(1) = NamedLogical(keyXGrid//"."//keyBuildFromWidths,.false.)
    constructGridFromWidths(2) = NamedLogical(keyVGrid//"."//keyBuildFromWidths,.false.)

    allocate(intParams(2))

    intParams(1) = NamedInteger(keyVGrid//"."//keyMaxL,0)
    intParams(2) = NamedInteger(keyVGrid//"."//keyMaxM,0)

    gridInMeters(1) = NamedLogical(keyXGrid//"."//keyLengthInMeters,.false.)

    call jsonCont%load(gridPoints)
    call jsonCont%load(cellWidth)
    call jsonCont%load(constructGridFromWidths)
    call jsonCont%load(gridInMeters)
    call jsonCont%load(intParams)

    if (assertions .or. assertionLvl >= 0) then 
        if (constructGridFromWidths(1)%value) then
            call assert(size(cellWidth(1)%values) > 0,"No "//cellWidth(1)%name//" entry found in JSON file by initGridFromJSON &
            &when "//constructGridFromWidths(1)%name//" was read as true")
        else
            call assert(size(gridPoints(1)%values) > 0,"No "//gridPoints(1)%name//" entry found in JSON file by initGridFromJSON")
        end if 
        if (constructGridFromWidths(2)%value) then 
            call assert(size(cellWidth(2)%values) > 0,"No "//cellWidth(2)%name//" entry found in JSON file by initGridFromJSON &
            &when "//constructGridFromWidths(2)%name//" was read as true")
        else
            call assert(size(gridPoints(2)%values) > 0,"No "//gridPoints(2)%name//" entry found in JSON file by initGridFromJSON")
        end if
    end if

    if (constructGridFromWidths(1)%value) then
        allocate(usedXGrid(size(cellWidth(1)%values)))

        usedXGrid(1) = cellWidth(1)%values(1)/2

        do i = 2, size(usedXGrid)
            usedXGrid(i) = usedXGrid(i-1) + cellWidth(1)%values(i-1)/2 + cellWidth(1)%values(i)/2
        end do
    else
        usedXGrid = gridPoints(1)%values
    end if 
    if (constructGridFromWidths(2)%value) then 
        allocate(usedVGrid(size(cellWidth(2)%values)))

        usedVGrid(1) = cellWidth(2)%values(1)/2

        do i = 2, size(usedVGrid)
            usedVGrid(i) = usedVGrid(i-1) + cellWidth(2)%values(i-1)/2 + cellWidth(2)%values(i)/2
        end do
    else
        usedVGrid = gridPoints(2)%values
    end if

    if (.not. gridInMeters(1)%value) usedNorm = real(1,kind=rk)
    call gridObj%init(usedXGrid/usedNorm,usedVGrid,intParams(1)%value,intParams(2)%value)

    gridPoints(1)%values = usedXGrid 
    gridPoints(2)%values = usedVGrid

    call jsonCont%output(gridPoints)
    call jsonCont%output(cellWidth)
    call jsonCont%output(constructGridFromWidths)
    call jsonCont%output(gridInMeters)
    call jsonCont%output(intParams)


end subroutine initGridFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initPartFromJSON(partObj,gridObj,jsonCont,mpiCont)
    !! Initialize simple partition object based on data from a JSON file and grid object

    type(Partition)         ,intent(inout) :: partObj
    type(Grid)              ,intent(in)    :: gridObj     !! Grid object used to get numX and numH
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedInteger) ,allocatable ,dimension(:) :: intParams

    integer(ik) :: numX, numH ,numProcs

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. partObj%isDefined(),"initPartFromJSON called on already defined partition object")
        call assert(gridObj%isDefined(),"Undefined grid object passed to initPartFromJSON")
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to initPartFromJSON")
    end if

    numX = gridObj%getNumX()
    numH = gridObj%getNumH()

    allocate(intParams(2))

    numProcs = mpiCont%getWorldSize()
    intParams(1) = NamedInteger(keyMPI//"."//keyNumPX,numProcs)
    intParams(2) = NamedInteger(keyMPI//"."//keyNumPH,1)

    call jsonCont%load(intParams)
    
    if (assertions .or. assertionLvl >= 0) call assert(intParams(1)%value*intParams(2)%value == numProcs,&
    intParams(1)%name//" and "//intParams(2)%name//" in config.json do not conform to total number of MPI processes")
    
    call partObj%initSimplePartition(intParams(1)%value,intParams(2)%value,numX,numH)

    call jsonCont%output(intParams)

end subroutine initPartFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initGeometryFromJSON(geometryObj,gridObj,jsonCont,mpiCont)
    !! Initialize geometry object based on data from a JSON file and grid object

    type(Geometry)          ,intent(inout) :: geometryObj
    type(Grid)              ,intent(in)    :: gridObj     !! Grid object used to infer cell widths for consistency
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    real(rk) ,allocatable ,dimension(:) :: xGrid, inferredWidths
    type(NamedRealArray)  ,dimension(1) :: cellFaceJacobians 
    type(NamedLogical)    ,dimension(1) :: periodicGrid

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. geometryObj%isDefined(),"initGeometryFromJSON called on already defined geometry object")
        call assert(gridObj%isDefined(),"Undefined grid object passed to initGeometryFromJSON")
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to initGeometryFromJSON")
    end if

    xGrid = gridObj%getXGrid()
    !Infer x grid widths assuming that left grid boundary is at x=0 
    allocate(inferredWidths(size(xGrid)))

    inferredWidths(1) = 2*xGrid(1)

    do i = 2,size(xGrid)
        inferredWidths(i) = 2*(xGrid(i) - sum(inferredWidths(1:i-1)))
    end do

    cellFaceJacobians(1)%name = keyXGrid//"."//keyFaceJacobians
    allocate(cellFaceJacobians(1)%values(size(xGrid)+1))
    cellFaceJacobians(1)%values = real(1,kind=rk)

    periodicGrid = NamedLogical(keyXGrid//"."//keyPeriodic,.false.)

    call jsonCont%load(cellFaceJacobians)
    call jsonCont%load(periodicGrid)
    if (assertions .or. assertionLvl >= 0) &
    call assert(size(cellFaceJacobians(1)%values) == size(xGrid)+1,&
    cellFaceJacobians(1)%name//" from config.json do not conform to xGrid")

    call geometryObj%init(inferredWidths,cellFaceJacobians(1)%values(1:size(xGrid)),cellFaceJacobians(1)%values(2:size(xGrid)+1),&
                         periodicGrid(1)%value) 

    call jsonCont%output(cellFaceJacobians)
    call jsonCont%output(periodicGrid)

end subroutine initGeometryFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initPETScContFromJSON(petscCont,indexingObj,jsonCont,mpiCont)
    !! Initialize PETSc controller object based on data from a JSON file and indexing object

    type(PETScController)   ,intent(inout) :: petscCont
    type(Indexing)          ,intent(in)    :: indexingObj 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController and to initialize PETSc

    type(NamedInteger) ,dimension(2) :: intParams 
    type(NamedReal)    ,dimension(3) :: realParams 
    type(NamedString)  ,dimension(3) :: charParams

    type(SolverOptions) :: solOptions

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. petscCont%isDefined(),"initPETScContFromJSON called on already defined controller object")
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to initPETScContFromJSON")
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to initPETScContFromJSON")
    end if

    intParams(1) = NamedInteger(keyPETSc//"."//keySolOptions//"."//keyMaxSolIters,10000)
    intParams(2) = NamedInteger(keyPETSc//"."//keyObjGroups,1)
    realParams(1) = NamedReal(keyPETSc//"."//keySolOptions//"."//keySolRelTol,real(1.0d-17,kind=rk))
    realParams(2) = NamedReal(keyPETSc//"."//keySolOptions//"."//keySolAbsTol,real(1.0d-20,kind=rk))
    realParams(3) = NamedReal(keyPETSc//"."//keySolOptions//"."//keySolDivTol,real(1.0d07,kind=rk))
    charParams(1) = NamedString(keyPETSc//"."//keySolOptions//"."//keySolKSP,"bcgs")
    charParams(2) = NamedString(keyPETSc//"."//keySolOptions//"."//keyHyprePC,"euclid")
    charParams(3) = NamedString(keyPETSc//"."//keySolOptions//"."//keyPETScOpts,"")

    call jsonCont%load(realParams)
    call jsonCont%load(intParams)
    call jsonCont%load(charParams)

    solOptions%kspSolverType = charParams(1)%value
    solOptions%solverToleranceRel =  realParams(1)%value
    solOptions%solverToleranceAbs = realParams(2)%value
    solOptions%solverToleranceDiv = realParams(3)%value
    solOptions%maxSolverIters =  intParams(1)%value 
    solOptions%hyprePC = charParams(2)%value
    solOptions%petscOptions = charParams(3)%value
    call petscCont%init(indexingObj,mpiCont,solOptions,numObjs=intParams(2)%value)

    call jsonCont%output(realParams)
    call jsonCont%output(intParams)
    call jsonCont%output(charParams)

end subroutine initPETScContFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initHDF5ContFromJSON(hdf5Cont,varCont,jsonCont,mpiCont)
    !! Initialize HDF5 controller object based on data from a JSON file and variable container object

    type(HDF5Controller)    ,intent(inout) :: hdf5Cont
    type(VariableContainer) ,intent(in)    :: varCont 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController and to initialize PETSc

    type(NamedString)       ,dimension(1) :: filepath
    type(NamedStringArray)  ,dimension(1) :: varNames

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. hdf5Cont%isDefined(),"initHDF5ContFromJSON called on already defined controller object")
        call assert(varCont%isDefined(),"Undefined variable container object passed to initHDF5ContFromJSON")
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to initHDF5ContFromJSON")
    end if

    filepath(1) = NamedString(keyHDF5//"."//keyFilepath,"./")

    varNames(1)%name = keyHDF5//"."//keyOutputVars
    allocate(varNames(1)%values(0))


    call jsonCont%load(filepath)
    call jsonCont%load(varNames)

    call hdf5Cont%init(varCont,varNames(1)%values,filepath(1)%value)
    
    call jsonCont%output(filepath)
    call jsonCont%output(varNames)
    
end subroutine initHDF5ContFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVarContFromJSON(varCont,indexingObj,partObj,textbookObj,jsonCont,mpiCont)
    !! Initialize variable container based on variable data from a JSON file. 

    type(VariableContainer) ,intent(inout) :: varCont
    type(Indexing)          ,intent(in)    :: indexingObj !! Indexing object used in variable container initialization 
    type(Partition)         ,intent(in)    :: partObj     !! Partition object used in variable container initialization 
    type(Textbook)          ,intent(in)    :: textbookObj !! Textbook object used to retrieve derivation rules for variable container 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(StringArray) ,allocatable ,dimension(:) :: implicitVarNames ,derivedVarNames
    type(NamedStringArray) ,allocatable ,dimension(:) :: derivationReqVarNames 
    type(NamedString) ,allocatable ,dimension(:) :: derivationNames
    type(NamedRealArray) ,allocatable ,dimension(:) :: initValsImplicit, initValsDerived
    type(NamedInteger) ,dimension(1) :: haloWidth
    type(CalculationRule) ,allocatable ,dimension(:) :: cRules

    type(VariableList) :: implicitVars ,derivedVars
    class(Derivation) ,allocatable :: derivBuffer 
    integer(ik) :: numImplicitVars,numDerivedVars

    integer(ik) :: i ,numX, numH, numV, minX, maxX

    if (assertions .or. assertionLvl >= 0) then 

        call assert(.not. varCont%isDefined(),"initVarContFromJSON called on already defined variable container")
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to initVarContFromJSON")
        call assert(partObj%isDefined(),"Undefined partition object passed to initVarContFromJSON")
        call assert(textbookObj%isDefined(),"Undefined textbook object passed to initVarContFromJSON")
        call assert(mpiCont%isDefined(),"Undefined MPIController passed to initVarContFromJSON")

    end if

    call initVarListFromJSON(implicitVars,jsonCont,mpiCont)
    call initVarListFromJSON(derivedVars,jsonCont,mpiCont,isDerivedList=.true.)

    haloWidth(1) = NamedInteger(keyMPI//"."//keyXHaloWidth,1)

    implicitVarNames = implicitVars%getVarNames()
    derivedVarNames = derivedVars%getVarNames()
    
    numImplicitVars = size(implicitVarNames) 
    numDerivedVars = size(derivedVarNames)

    if (assertions .or. assertionLvl >= 0) &
    call assert(numImplicitVars + numDerivedVars > 0,"No variables detected in config.json file by initVarContFromJSON")

    allocate(derivationReqVarNames(numDerivedVars))
    allocate(derivationNames(numDerivedVars))
    if (numDerivedVars > 0) then  

        do i = 1, numDerivedVars

            derivationNames(i) = NamedString(keyVariables//"."//keyDerivedVars//"."//&
            derivedVarNames(i)%string//"."//keyDerivRule//"."//keyRuleName,keyNone)

            derivationReqVarNames(i)%name = keyVariables//"."//keyDerivedVars//"."//&
            derivedVarNames(i)%string//"."//keyDerivRule//"."//keyReqVarNames
            
            allocate(derivationReqVarNames(i)%values(0))
        end do

    end if

    call jsonCont%load(derivationNames)
    call jsonCont%load(derivationReqVarNames)

    allocate(cRules(numDerivedVars))

    if (numDerivedVars > 0) then 
        
        do i = 1,numDerivedVars
            if (derivationNames(i)%value /= keyNone) then 
                call textbookObj%copyDerivation(derivationNames(i)%value,derivBuffer)
                call cRules(i)%init(derivBuffer,derivationReqVarNames(i)%values)
            else 
                call cRules(i)%init()
            end if
        end do
    end if

    call jsonCont%load(haloWidth)

    call varCont%init(implicitVars,derivedVars,cRules,indexingObj,partObj,haloWidth(1)%value,mpiCont%getWorldRank())

    numX = indexingObj%getNumX()
    numH = indexingObj%getNumH()
    numV = indexingObj%getNumV()
    minX = partObj%getMinXAtInd(mpiCont%getWorldRank()+1)
    maxX = partObj%getMaxXAtInd(mpiCont%getWorldRank()+1)

    allocate(initValsImplicit(numImplicitVars))
    allocate(initValsDerived(numDerivedVars))

    !Initialize implicit variables
    do i = 1,numImplicitVars
        initValsImplicit(i)%name=keyVariables//"."//keyImplicitVars//"."//implicitVarNames(i)%string//"."//keyInitVals
        allocate(initValsImplicit(i)%values(0))
    end do

    call jsonCont%load(initValsImplicit)

    do i = 1,numImplicitVars
        if (size(initValsImplicit(i)%values) > 0) then
            if (implicitVars%isVarDist(i)) then 
                if (assertions .or. assertionLvl >= 0) call assert(size(initValsImplicit(i)%values) == numX*numH*numV,&
                "Initial value vector "//initValsImplicit(i)%name//" does not conform to grid")
                varCont%variables(varCont%getVarIndex(implicitVarNames(i)%string))%entry(1:(maxX-minX+1)*numH*numV) = &
                initValsImplicit(i)%values((minX-1)*numH*numV+1:maxX*numH*numV)
                if (haloWidth(1)%value > 0) then
                    if (minX > 1) varCont%variables(varCont%getVarIndex(implicitVarNames(i)%string))%entry&
                    (1-haloWidth(1)%value*numH*numV:0) &
                    = initValsImplicit(i)%values((minX-1-haloWidth(1)%value)*numH*numV:(minX-1)*numH*numV-1)

                    if (maxX < numX) &
                    varCont%variables(varCont%getVarIndex(implicitVarNames(i)%string))%entry&
                    ((maxX-minX+1)*numH*numV + 1:(maxX-minX+1+haloWidth(1)%value)*numH*numV) &
                    = initValsImplicit(i)%values(maxX*numH*numV + 1:(maxX+haloWidth(1)%value)*numH*numV)
                end if

            else if (implicitVars%isVarScalar(i)) then 

                if (assertions .or. assertionLvl >= 0) call assert(size(initValsImplicit(i)%values) == 1,&
                "Initial value vector "//initValsImplicit(i)%name//" should be size 1")

                varCont%variables(varCont%getVarIndex(implicitVarNames(i)%string))%entry = &
                initValsImplicit(i)%values

            else

                if (assertions .or. assertionLvl >= 0) call assert(size(initValsImplicit(i)%values) == numX,&
                "Initial value vector "//initValsImplicit(i)%name//" does not conform to grid")
                varCont%variables(varCont%getVarIndex(implicitVarNames(i)%string))%entry(1:maxX-minX+1) = &
                initValsImplicit(i)%values(minX:maxX)

                if (haloWidth(1)%value > 0) then
                    if (minX > 1) varCont%variables(varCont%getVarIndex(implicitVarNames(i)%string))%entry(1-haloWidth(1)%value:0) &
                    = initValsImplicit(i)%values(minX-haloWidth(1)%value:minX-1)

                    if (maxX < numX) &
                    varCont%variables(varCont%getVarIndex(implicitVarNames(i)%string))%entry&
                    (maxX-minX+2:maxX-minX+1+haloWidth(1)%value) &
                    = initValsImplicit(i)%values(maxX + 1:maxX+haloWidth(1)%value)
                end if

            end if
            
        end if
    end do

    !Initialize derived variables
    do i = 1,numDerivedVars
        initValsDerived(i)%name=keyVariables//"."//keyDerivedVars//"."//derivedVarNames(i)%string//"."//keyInitVals
        allocate(initValsDerived(i)%values(0))
    end do

    call jsonCont%load(initValsDerived)

    do i = 1,numDerivedVars
        if (size(initValsDerived(i)%values) > 0) then 

            if (derivedVars%isVarDist(i)) then 

                if (assertions .or. assertionLvl >= 0) call assert(size(initValsDerived(i)%values) == numX*numH*numV,&
                "Initial value vector "//initValsDerived(i)%name//" does not conform to grid")
                varCont%variables(varCont%getVarIndex(derivedVarNames(i)%string))%entry(1:(maxX-minX+1)*numH*numV) = &
                initValsDerived(i)%values((minX-1)*numH*numV+1:maxX*numH*numV)

                if (haloWidth(1)%value > 0) then
                    if (minX > 1) varCont%variables(varCont%getVarIndex(derivedVarNames(i)%string))%entry&
                    (1-haloWidth(1)%value*numH*numV:0) &
                    = initValsDerived(i)%values((minX-1-haloWidth(1)%value)*numH*numV:(minX-1)*numH*numV-1)

                    if (maxX < numX) &
                    varCont%variables(varCont%getVarIndex(derivedVarNames(i)%string))%entry&
                    ((maxX-minX+1)*numH*numV + 1:(maxX-minX+1+haloWidth(1)%value)*numH*numV) &
                    = initValsDerived(i)%values(maxX*numH*numV + 1:(maxX+haloWidth(1)%value)*numH*numV)
                end if

            else if (derivedVars%isVarScalar(i)) then 

                if (assertions .or. assertionLvl >= 0) call assert(size(initValsDerived(i)%values) == 1,&
                "Initial value vector "//initValsDerived(i)%name//" should be size 1")

                varCont%variables(varCont%getVarIndex(derivedVarNames(i)%string))%entry = &
                initValsDerived(i)%values

            else

                if (assertions .or. assertionLvl >= 0) call assert(size(initValsDerived(i)%values) == numX,&
                "Initial value vector "//initValsDerived(i)%name//" does not conform to grid")
                varCont%variables(varCont%getVarIndex(derivedVarNames(i)%string))%entry(1:maxX-minX+1) = &
                initValsDerived(i)%values(minX:maxX)

                if (haloWidth(1)%value > 0) then
                    if (minX > 1) varCont%variables(varCont%getVarIndex(derivedVarNames(i)%string))%entry(1-haloWidth(1)%value:0) &
                    = initValsDerived(i)%values(minX-haloWidth(1)%value:minX-1)

                    if (maxX < numX) &
                    varCont%variables(varCont%getVarIndex(derivedVarNames(i)%string))%entry&
                    (maxX-minX+2:maxX-minX+1+haloWidth(1)%value) &
                    = initValsDerived(i)%values(maxX + 1:maxX+haloWidth(1)%value)
                end if

            end if
            
        end if
    end do

    call jsonCont%output(derivationNames)
    call jsonCont%output(derivationReqVarNames)
    call jsonCont%output(haloWidth)

    call jsonCont%output(initValsImplicit)
    call jsonCont%output(initValsDerived)


end subroutine initVarContFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStandardTextbook(textbookObj,gridObj,geometryObj,partObj,&
                                      vSpaceObj,normObj,speciesListObj,varList,jsonCont,mpiCont)
    !! Initialize standard textbook object based on grid, geometry, partition, species, and normalization data as well as a JSON file. 
    !! Assumes that the velocity grid is normalized to electron thermal speed.
    !! The textbook includes (passed variables are in appearance order unless otherwise specified):
    !!
    !! "flowSpeedFromFlux" = flux/density (normalized to n_0 * u_0)
    !! "sonicSpeed" for each ion species  (negative IDs) = sqrt((poly_e * T_e + poly_i * T_i)/m_i) where i is the index of the ion species (normalized to u_0)
    !! "tempFromEnergy" for each species specified in config file as requiring temperature derivation = 2/3 * energy/density - mass * flux**2/(3*density**2)
    !! "leftElectronGamma" = sheath heat transmission coeffiecient for electrons at the left boundary using T_e and T_i
    !! "rightElectronGamma" = sheath heat transmission coefficient for electrons at the right boundary using T_e and T_i
    !! "densityMoment" = zeroth moment of f_0
    !! "energyMoment" = second order of f_0
    !! "cclDragCoeff" = Chang-Cooper-Langdon drag coefficient - assumes passed variable is the single harmonic f_0
    !! "cclDiffusionCoeff" = Chang-Cooper-Langdon diffusion coefficient - assumes passed variables are the single harmonic f_0 and cclWeight
    !! "cclWeight" = Chang-Cooper-Langdon interpolation weight - assumes passed variables are the single harmonic cclDragCoeff and cclDiffusionCoeff
    !! if l=1 is resolved:
    !! "fluxMoment" = first moment of f_1/3 (normalized to n_0*u_0)
    !! "heatFluxMoment" = third moment of f_1/3 (normalized to m_e * n_0 * v_th**3 / 2)
    !! if l=2 is resolved:
    !! "viscosityTensorxxMoment" = second moment of 2*f_2/15 (normalized to n_0 * e * T_0)
    !! Interpolation for staggered grids:   
    !! "gridToDual" = interpolates one variable from regular to staggered(dual) grid
    !! "dualToGrid" = interpolates one variable from staggered(dual) to regular grid
    !! "distributionInterp" = interpolates one distribution function (even harmonics to dual, odd to regular grid)
    !! Other useful derivations:
    !! "gradDeriv" = calculates gradient of variable on regular grid, interpolating on cell faces and extrapolating at boundaries
    !! "logLei" = calculates electron-ion Coulomb log taking in electron temperature and density (derivations added for each ion species)
    !! "logLee" = calculates electron self collision Coulomb log taking in electron temperature and density 
    !! "logLii" = calculates ion-ion collision frequency for each ion collision combination taking in the two species temperatures and 
    !!            densities (used as"logLiis_S" where s and S are the first and second ion species name) 
    !! "maxwellianDistribution" = calculates distribution function with Maxwellian l=0 harmonic and all other harmonics 0. Takes in temperature and density
    !! Also automatically includes all defined custom derivations

    type(Textbook)          ,intent(inout) :: textbookObj 
    type(Grid)              ,intent(in)    :: gridObj    
    type(Partition)         ,intent(in)    :: partObj
    type(Geometry)          ,intent(in)    :: geometryObj
    type(VSpace)            ,intent(in)    :: vSpaceObj
    class(Normalization)    ,intent(in)    :: normObj 
    type(SpeciesList)       ,intent(in)    :: speciesListObj  
    type(VariableList)      ,intent(in)    :: varList
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedIntegerArray) ,dimension(1) :: tempDerivIDs 
    type(NamedInteger)      ,dimension(1) :: sheathGammaIonSpeciesID
    type(NamedReal)      ,dimension(2) :: polytropicCoeffs

    integer(ik) ,allocatable ,dimension(:) :: allIDs, negativeIDs ,weightedTempIDs

    type(SimpleDerivation)                                :: flowSpeedDeriv ,tempDerivFirstTerm
    type(SimpleDerivation)   ,allocatable ,dimension(:)   :: weightedTemp ,tempDerivSecondTerm
    type(MomentDerivation)   ,allocatable ,dimension(:)   :: momentDerivs
    type(AdditiveDerivation) ,allocatable ,dimension(:)   :: tempDeriv ,sonicSpeedDeriv
    type(ElectronSheathGammaDerivation)                   :: leftGammaDeriv ,rightGammaDeriv
    type(InterpolationDerivation)                         :: gridToDualDeriv, dualToGridDeriv ,distInterp
    type(CentralDiffDerivation)                           :: gradDeriv 
    type(CoulombLogDerivation) ,allocatable ,dimension(:) :: logLeiDeriv ,logLiiDeriv
    type(CoulombLogDerivation)                            :: logLeeDeriv
    type(MaxwellianDerivation)                            :: maxwellianDeriv

    type(CCLDiffDerivation) :: cclDiffDeriv
    type(CCLDragDerivation) :: cclDragDeriv
    type(CCLWeightDerivation) :: cclWeightDeriv

    class(Species) ,allocatable :: tempSpecies ,tempSpecies2

    integer(ik) :: i ,j ,k ,numMomentDerivs ,locNumX ,minX ,maxX ,numLiiDerivs

    real(rk) :: tempNorm, speedNorm ,velSpaceNorm ,speciesTKineticRatio, densNorm
    
    logical ,allocatable ,dimension(:) :: selfIonLambdaAdded

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. textbookObj%isDefined(),"Attempted to initialize standard textbook in already defined texbook object")
        call assert(gridObj%isDefined(),"Undefined grid object passed to initStandardTextbook")
        call assert(partObj%isDefined(),"Undefined partition object passed to initStandardTextbook")
        call assert(geometryObj%isDefined(),"Undefined geometry object passed to initStandardTextbook")
        call assert(vSpaceObj%isDefined(),"Undefined velocity space object passed to initStandardTextbook")
        call assert(normObj%isDefined(),"Undefined normalization object passed to initStandardTextbook")
        call assert(speciesListObj%isDefined(),"Undefined species list passed to initStandardTextbook")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to initStandardTextbook")
    end if

    minX = partObj%getMinXAtInd(mpiCont%getWorldRank()+1) 
    maxX = partObj%getMaxXAtInd(mpiCont%getWorldRank()+1)
    locNumX =  maxX - minX + 1
    tempDerivIDs(1)%name = keyStandardTextbook//"."//keyTempDerivIDs
    allocate(tempDerivIDs(1)%values(0))

    polytropicCoeffs(1) = NamedReal(keyStandardTextbook//"."//keyElectronPolytropicC,real(1,kind=rk))
    polytropicCoeffs(2) = NamedReal(keyStandardTextbook//"."//keyIonPolytropicC,real(1,kind=rk))

    sheathGammaIonSpeciesID(1) = NamedInteger(keyStandardTextbook//"."//keySheathGammaIonSpeciesID,-1)

    call jsonCont%load(tempDerivIDs)
    call jsonCont%load(sheathGammaIonSpeciesID)
    call jsonCont%load(polytropicCoeffs)
    call jsonCont%output(tempDerivIDs)
    call jsonCont%output(polytropicCoeffs)
    call jsonCont%output(sheathGammaIonSpeciesID)

    call flowSpeedDeriv%init(real([1,-1],kind=rk))
    call tempDerivFirstTerm%init(real([1,-1],kind=rk),multConst=real(2.0d0/3.0d0,kind=rk))

    tempNorm = normObj%getNormalizationValue(keyTempEVNorm)
    speedNorm = normObj%getNormalizationValue(keySpeedNorm)
    densNorm = normObj%getNormalizationValue(keyDensNorm)
    velSpaceNorm = normObj%getNormalizationValue(keyVelGridNorm)

    allIDs = speciesListObj%getSpeciesIDs()
    negativeIDs = pack(allIDs,allIDs < 0)
    numMomentDerivs = 2 !Density and energy if only l=0 is resolved
    if (gridObj%getMaxL() > 0) numMomentDerivs = numMomentDerivs + 2 !Flow speed and heat flux
    if (gridObj%getMaxL() > 1) numMomentDerivs = 5 !Viscosity tensor xx

    call textbookObj%init()
    call textbookObj%addDerivation(flowSpeedDeriv,"flowSpeedFromFlux")
    allocate(weightedTemp(2*size(negativeIDs)))
    allocate(sonicSpeedDeriv(size(negativeIDs)))
    ! polytropicCoeff_e*T_e/m_i and polytropicCoeff_i*T_i/m_i derivations
    do i = 1,size(negativeIDs)
        if (allocated(tempSpecies)) deallocate(tempSpecies)
        allocate(tempSpecies,source = speciesListObj%getSpeciesFromID(negativeIDs(i)))
        speciesTKineticRatio = elCharge * tempNorm/(tempSpecies%getMass()*speedNorm**2)
        call weightedTemp(i)%init(real([1],kind=rk),multConst=polytropicCoeffs(1)%value*speciesTKineticRatio)
        call weightedTemp(i+size(negativeIDs))%init(real([1],kind=rk),multConst=polytropicCoeffs(2)%value*speciesTKineticRatio)
        call sonicSpeedDeriv(i)%init(2,2,resultPower=real(0.5d00,kind=rk))
        call sonicSpeedDeriv(i)%addDerivation(weightedTemp(i),[1])
        call sonicSpeedDeriv(i)%addDerivation(weightedTemp(i+size(negativeIDs)),[2])
        call textbookObj%addDerivation(sonicSpeedDeriv(i),"sonicSpeed"//tempSpecies%getName())
    end do

    allocate(tempDerivSecondTerm(size(tempDerivIDs(1)%values)))
    allocate(tempDeriv(size(tempDerivIDs(1)%values)))


    ! Allocate - m * u_0 ^2/(2*e*T_0) * u^2 term where u is the derived flow speed
    do i = 1,size(tempDerivIDs(1)%values)
        if (allocated(tempSpecies)) deallocate(tempSpecies)
        allocate(tempSpecies,source = speciesListObj%getSpeciesFromID(tempDerivIDs(1)%values(i)))
        speciesTKineticRatio = tempSpecies%getMass()*speedNorm**2/(3 * elCharge * tempNorm)
        call tempDerivSecondTerm(i)%init(real([2,-2],kind=rk),multConst=-speciesTKineticRatio)
        call tempDeriv(i)%init(2,3)
        call tempDeriv(i)%addDerivation(tempDerivFirstTerm,[1,2])
        call tempDeriv(i)%addDerivation(tempDerivSecondTerm(i),[3,2])
        call textbookObj%addDerivation(tempDeriv(i),"tempFromEnergy"//tempSpecies%getName())

    end do

    ! Moment derivations
    allocatE(momentDerivs(numMomentDerivs))
    call momentDerivs(1)%init(0,gridObj%getH(0,0,.false.),vSpaceObj) !density derivation
    call momentDerivs(2)%init(2,gridObj%getH(0,0,.false.),vSpaceObj) !energy derivation (assuming that the velocity grid is normalized to electrons)
    call textbookObj%addDerivation(momentDerivs(1),"densityMoment")
    call textbookObj%addDerivation(momentDerivs(2),"energyMoment")
    if (gridObj%getMaxL() > 0) then 
        call momentDerivs(3)%init(1,gridObj%getH(1,0,.false.),vSpaceObj,multConst=velSpaceNorm/(3*speedNorm)) !flux moment
        call momentDerivs(4)%init(3,gridObj%getH(1,0,.false.),vSpaceObj,multConst=real(1.0d0/3.0d0,kind=rk)) !heat flux moment (assuming v grid normalized to electrons)
        call textbookObj%addDerivation(momentDerivs(3),"fluxMoment")
        call textbookObj%addDerivation(momentDerivs(4),"heatFluxMoment")
    end if

    if (gridObj%getMaxL() > 1) then 
        call momentDerivs(5)%init(2,gridObj%getH(2,0,.false.),vSpaceObj,multConst=real(4.0d0/15.0d0,kind=rk)) ! xx component of viscosity tensor, normalized to n_0 * e*T_0  and assuming v grid normalized to electrons)
        call textbookObj%addDerivation(momentDerivs(5),"viscosityTensorxxMoment")
    end if

    ! Sheath heat transmission coeff derivations
    if (allocated(tempSpecies)) deallocate(tempSpecies)
    if (any(speciesListObj%getSpeciesIDs()==sheathGammaIonSpeciesID(1)%value)) then
        allocate(tempSpecies,source = speciesListObj%getSpeciesFromID(sheathGammaIonSpeciesID(1)%value))
        call leftGammaDeriv%init(elMass/tempSpecies%getMass(),1)
        call rightGammaDeriv%init(elMass/tempSpecies%getMass(),locNumX)
        call textbookObj%addDerivation(leftGammaDeriv,"leftElectronGamma")
        call textbookObj%addDerivation(rightGammaDeriv,"rightElectronGamma")
    else
        call printMessage("Warning: sheathGammaIonSpeciesID not found in species list, &
        &no electron gamma derivations added to standard textbook")
    end if

    if (gridObj%getNumX()>1) then

        call gridToDualDeriv%init(geometryObj,gridObj,minX,maxX)
        call dualToGridDeriv%init(geometryObj,gridObj,minX,maxX,inverseInterp=.true.)
        call distInterp%init(geometryObj,gridObj,minX,maxX,distInterp=.true.)

        call textbookObj%addDerivation(gridToDualDeriv,"gridToDual")
        call textbookObj%addDerivation(dualToGridDeriv,"dualToGrid")
        call textbookObj%addDerivation(distInterp,"distributionInterp")

        call gradDeriv%init(geometryObj,partObj,mpiCont%getWorldRank())
        call textbookObj%addDerivation(gradDeriv,"gradDeriv")
    end if

    call logLeeDeriv%init(real(1,kind=rk),locNumX,densNorm,tempNorm,electronLog=.true.)
    call textbookObj%addDerivation(logLeeDeriv,"logLee")

    allocate(logLeiDeriv(size(negativeIDs)))
    allocate(selfIonLambdaAdded(size(negativeIDs)))
    numLiiDerivs = size(negativeIDs)**2
    allocate(logLiiDeriv(numLiiDerivs))
    selfIonLambdaAdded = .false.
    k=1
    do i = 1, size(negativeIDs)
        if (allocated(tempSpecies)) deallocate(tempSpecies)
        allocate(tempSpecies,source = speciesListObj%getSpeciesFromID(negativeIDs(i)))
        call logLeiDeriv(i)%init(tempSpecies%getCharge(),locNumX,densNorm,tempNorm)
        call textbookObj%addDerivation(logLeiDeriv(i),"logLei"//tempSpecies%getName())
        
        do j = 1,size(negativeIDs)
            if (i /= j .or. .not. selfIonLambdaAdded(j)) then
                if (allocated(tempSpecies2)) deallocate(tempSpecies2)
                allocate(tempSpecies2,source = speciesListObj%getSpeciesFromID(negativeIDs(j)))
                call logLiiDeriv(k)%init(tempSpecies%getCharge(),locNumX,densNorm,tempNorm,ionZ2=tempSpecies2%getCharge(),&
                                         ionLog = .true. ,ionMassRatio=tempSpecies%getMass()/tempSpecies2%getMass())
                call textbookObj%addDerivation(logLiiDeriv(k),"logLii"//tempSpecies%getName()//"_"//tempSpecies2%getName())
                k = k + 1
                if (i == j) selfIonLambdaAdded(j) = .true.
            end if
        end do
    end do

    call maxwellianDeriv%init(vSpaceObj)
    call textbookObj%addDerivation(maxwellianDeriv,"maxwellianDistribution")
    
    call cclDiffDeriv%init(vSpaceObj)
    call cclDragDeriv%init(vSpaceObj)
    call cclWeightDeriv%init(vSpaceObj)

    call textbookObj%addDerivation(cclDiffDeriv,"cclDiffusionCoeff")
    call textbookObj%addDerivation(cclDragDeriv,"cclDragCoeff")
    call textbookObj%addDerivation(cclWeightDeriv,"cclWeight")

    call addCustomDerivationsToTextbook(textbookObj,gridObj,geometryObj,partObj,&
                                            vSpaceObj,normObj,speciesListObj,varList,jsonCont,mpiCont)

end subroutine initStandardTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStandardIntegrator(integratorObj,varCont,indexingObj,jsonCont,mpiCont)
    !! Initialize standard integrator based on data from a JSON file. The standard integrator is a composit integrator containing
    !! either RK, BDE, or CVODE integrators. 

    type(CompositeIntegrator) ,intent(inout) :: integratorObj 
    type(VariableContainer)   ,intent(in)    :: varCont
    type(Indexing)            ,intent(in)    :: indexingObj !! Indexing object used in BDE integrator initialization
    type(JSONController)      ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)       ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(IntegratorCallStep)   ,allocatable ,dimension(:) :: integratorSteps
    type(SimpleTimestepController)                        :: dtController
    type(ExplicitRKIntegrator) ,allocatable  :: integratorRK 
    type(PicardBDEIntegrator)  ,allocatable  :: integratorBDE
    type(CVODEIntegrator)      ,allocatable  :: integratorCVODE
    type(CVODEOptions)         ,allocatable  :: optionsCVODE

    type(NamedStringArray)                  ,dimension(1) :: integratorTags, stepTags ,modelTags
    type(NamedString)          ,allocatable ,dimension(:) :: integratorType
    type(NamedString)          ,allocatable ,dimension(:) :: stringParams
    type(NamedStringArray)     ,allocatable ,dimension(:) :: stringArrayParams
    type(NamedInteger)         ,allocatable ,dimension(:) :: integerParams
    type(NamedIntegerArray)    ,allocatable ,dimension(:) :: integerArrayParams
    type(NamedReal)            ,allocatable ,dimension(:) :: realParams 
    type(NamedRealArray)       ,allocatable ,dimension(:) :: realArrayParams 
    type(NamedLogical)         ,allocatable ,dimension(:) :: logicalParams

    type(IntArray)             ,allocatable ,dimension(:) :: groupIndices 
    type(LogicalArray)         ,allocatable ,dimension(:) :: updatesOnInternalIterations

    type(CommunicationData) ,allocatable :: commData

    real(rk) :: initialTime ,initialTimestep

    integer(ik) ,allocatable ,dimension(:) :: convIndices ,modelIndices

    integer(ik) :: i ,j ,k ,numImplicitGroups ,numGeneralGroups,integratorIndex

    logical :: found ,nonTrivialUpdate ,nonTrivialModelDataUpdate ,commNeeded

    logical ,allocatable ,dimension(:) :: updatesOnInternalIterationsModelData

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. integratorObj%isDefined(),"Attempted to initialize standard integrator in already defined object")
        call assert(varCont%isDefined(),"Undefined variable container passed to initStandardIntegrator")
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to initStandardIntegrator")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to initStandardIntegrator")
    end if

    integratorTags(1)%name = keyIntegrator//"."//keyIntegratorTags
    allocate(integratorTags(1)%values(0))
    stepTags(1)%name = keyIntegrator//"."//keyStepTags
    allocate(stepTags(1)%values(0))

    modelTags(1)%name = keyModels//"."//keyTags
    allocate(modelTags(1)%values(0))

    call jsonCont%load(integratorTags)
    call jsonCont%load(stepTags)
    call jsonCont%load(modelTags)

    if (assertions .or. assertionLvl >= 0) then 
        call assert(size(integratorTags(1)%values) > 0,"No integrator tags found by initStandardIntegrator")
        call assert(size(stepTags(1)%values) > 0,"No step tags found by initStandardIntegrator")
        call assert(size(modelTags(1)%values) > 0,"No model tags found by initStandardIntegrator")

        do i =1, size(integratorTags(1)%values)
            do j = 1, i - 1
                call assert(.not. integratorTags(1)%values(j)%string == integratorTags(1)%values(i)%string,&
                "Duplicate integrator tag detected in initStandardIntegrator")
            end do
        end do

        do i =1, size(modelTags(1)%values)
            do j = 1, i - 1
                call assert(.not. modelTags(1)%values(j)%string == modelTags(1)%values(i)%string,&
                "Duplicate model tag detected in initStandardIntegrator")
            end do
        end do

    end if

    allocate(integratorType(size(integratorTags(1)%values)))

    do i = 1, size(integratorTags(1)%values)
        integratorType(i) = NamedString(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyType,"")
    end do
    call jsonCont%load(integratorType)

    initialTime = 0
    if (varCont%isVarNameRegistered("time")) initialTime = varCont%variables(varCont%getVarIndex("time"))%entry(1)

    allocate(realParams(1))
    realParams(1) = NamedReal(keyIntegrator//"."//keyInitialTimestep,real(0.1d0,kind=rk))

    call jsonCont%load(realParams)
    call jsonCont%output(realParams)
    initialTimestep = realParams(1)%value

    allocate(logicalParams(1))
    logicalParams(1) = NamedLogical(keyIntegrator//"."//keyTimestepController//"."//keyActive,.false.)

    call jsonCont%load(logicalParams)
    call jsonCont%output(logicalParams)

    allocate(integerParams(2))
    integerParams(1) = NamedInteger(keyIntegrator//"."//keyNumImplicitGroups,1)
    integerParams(2) = NamedInteger(keyIntegrator//"."//keyNumGenGroups,1)
    call jsonCont%load(integerParams)
    call jsonCont%output(integerParams)

    numImplicitGroups = integerParams(1)%value 
    numGeneralGroups = integerParams(2)%value 

    ! Determine whether a timestep controller is used and initialize the integrator object
    if (logicalParams(1)%value) then 
        deallocate(logicalParams)
        allocate(logicalParams(2))
        logicalParams(1) = NamedLogical(keyIntegrator//"."//keyTimestepController//"."//keyUseMaxVal,.false.)
        logicalParams(2) = NamedLogical(keyIntegrator//"."//keyTimestepController//"."//keyRescaleTimestep,.false.)
        call jsonCont%load(logicalParams)
        call jsonCont%output(logicalParams)

        allocate(stringArrayParams(1))
        stringArrayParams(1)%name = keyIntegrator//"."//keyTimestepController//"."//keyReqVarNames
        allocate(stringArrayParams(1)%values(0))
        call jsonCont%load(stringArrayParams)

        if (assertions .or. assertionLvl >= 0) call assert(size(stringArrayParams(1)%values)>0,&
        "No required variables found for timestep controller in&
        & initStandardIntegrator")

        call jsonCont%output(stringArrayParams)

        allocate(realArrayParams(1))
        realArrayParams(1)%name = keyIntegrator//"."//keyTimestepController//"."//keyReqVarPowers
        allocate(realArrayParams(1)%values(size(stringArrayParams(1)%values)))
        realArrayParams(1)%values = real(1,kind=rk)

        call jsonCont%load(realArrayParams)
        call jsonCont%output(realArrayParams)

        realParams(1) = NamedReal(keyIntegrator//"."//keyTimestepController//"."//keyMultConst,real(1,kind=rk))

        call jsonCont%load(realParams)
        call jsonCont%output(realParams)

        call dtController%init(mpiCont,varCont,stringArrayParams(1)%values,realArrayParams(1)%values&
                              ,multConst=realParams(1)%value,useMaxVal=logicalParams(1)%value&
                              ,rescaleTimestep=logicalParams(2)%value)

        call integratorObj%init(initialTime,initialTimestep,size(integratorType),dtController) 
    else
        call integratorObj%init(initialTime,initialTimestep,size(integratorType)) 
    end if

    ! Add integrators
    do i = 1, size(integratorType)
        select case (integratorType(i)%value)
        case ("RK")
            if (allocated(integerParams)) deallocate(integerParams)
            allocate(integerParams(1))

            integerParams(1) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyOrder,1)

            call jsonCont%load(integerParams)
            call jsonCont%output(integerParams)

            allocate(integratorRK)
            call integratorRK%init(order=integerParams(1)%value)

            call integratorObj%addIntegrator(integratorRK)

            deallocate(integratorRK)

        case ("BDE")

            if (allocated(integerParams)) deallocate(integerParams)
            allocate(integerParams(7))

            integerParams(1) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyMaxNonlinIters,100)
            integerParams(2) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyAssociatedPETScGroup,1)
            integerParams(3) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string&
                                            //"."//keyInternalStepControl//"."//keyStartingNumSteps,1)
            integerParams(4) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string&
                                            //"."//keyInternalStepControl//"."//keyStepMultiplier,2)
            integerParams(5) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string&
                                            //"."//keyInternalStepControl//"."//keyStepDecrament,1)
            integerParams(6) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string&
                                            //"."//keyInternalStepControl//"."//keyMinNumNonlinInters,5)

            integerParams(7) = NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string&
                                            //"."//keyInternalStepControl//"."//keyMaxBDERestarts,3)

            call jsonCont%load(integerParams)
            call jsonCont%output(integerParams)

            if (allocated(logicalParams)) deallocate(logicalParams)
            allocate(logicalParams(2))

            logicalParams(1) = NamedLogical(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyUse2Norm,.false.)
            logicalParams(2) = NamedLogical(keyIntegrator//"."//integratorTags(1)%values(i)%string&
                                            //"."//keyInternalStepControl//"."//keyActive,.false.)

            call jsonCont%load(logicalParams)
            call jsonCont%output(logicalParams)

            if (allocated(realParams)) deallocate(realParams)
            allocate(realParams(3))

            realParams(1) = NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyNonlinTol,&
                                      real(1.0d-12,kind=rk))
            realParams(2) = NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyAbsTol,&
                                      real(1,kind=rk))
            realParams(3) = NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyRelaxationWeight,&
                                      real(1,kind=rk))

            call jsonCont%load(realParams)
            call jsonCont%output(realParams)

            if (allocated(stringArrayParams)) deallocate(stringArrayParams)
            allocate(stringArrayParams(1))

            stringArrayParams(1)%name = keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyConvergenceVars
            allocate(stringArrayParams(1)%values(0))

            call jsonCont%load(stringArrayParams)
            call jsonCont%output(stringArrayParams)

            allocate(integratorBDE)

            if (size(stringArrayParams(1)%values) > 0) then 

                if (allocated(convIndices)) deallocate(convIndices)
                allocate(convIndices(size(stringArrayParams(1)%values)))

                do j = 1,size(stringArrayParams(1)%values)
                    convIndices(j) = varCont%getVarIndex(stringArrayParams(1)%values(j)%string)
                end do


                if (logicalParams(2)%value) then 
                    call integratorBDE%init(indexingObj,mpiCont%getWorldRank(),nonlinTol=realParams(1)%value&
                                            ,absTol=realParams(2)%value,maxIters=integerParams(1)%value&
                                            ,convergenceIndices=convIndices,petscGroup=integerParams(2)%value&
                                            ,use2Norm=logicalParams(1)%value&
                                            ,intContOptions=InternalControllerOptions(integerParams(3)%value,&
                                                                                      integerParams(4)%value,&
                                                                                      integerParams(5)%value,&
                                                                                      integerParams(6)%value,&
                                                                                      integerParams(7)%value)&
                                            ,integratorName=integratorTags(1)%values(i)%string&
                                            ,relaxationWeight=realParams(3)%value)
                else
                    call integratorBDE%init(indexingObj,mpiCont%getWorldRank(),nonlinTol=realParams(1)%value&
                                            ,absTol=realParams(2)%value,maxIters=integerParams(1)%value&
                                            ,convergenceIndices=convIndices,petscGroup=integerParams(2)%value&
                                            ,use2Norm=logicalParams(1)%value&
                                            ,integratorName=integratorTags(1)%values(i)%string&
                                            ,relaxationWeight=realParams(3)%value)
                end if

            else

                if (logicalParams(2)%value) then 
                    call integratorBDE%init(indexingObj,mpiCont%getWorldRank(),nonlinTol=realParams(1)%value&
                                            ,absTol=realParams(2)%value,maxIters=integerParams(1)%value&
                                            ,petscGroup=integerParams(2)%value,use2Norm=logicalParams(1)%value&
                                            ,intContOptions=InternalControllerOptions(integerParams(3)%value,&
                                                                                      integerParams(4)%value,&
                                                                                      integerParams(5)%value,&
                                                                                      integerParams(6)%value)&
                                            ,integratorName=integratorTags(1)%values(i)%string&
                                            ,relaxationWeight=realParams(3)%value)
                    
                else
                    call integratorBDE%init(indexingObj,mpiCont%getWorldRank(),nonlinTol=realParams(1)%value&
                                            ,absTol=realParams(2)%value,maxIters=integerParams(1)%value&
                                            ,petscGroup=integerParams(2)%value,use2Norm=logicalParams(1)%value&
                                            ,integratorName=integratorTags(1)%values(i)%string&
                                            ,relaxationWeight=realParams(3)%value)
                end if

            end if

            call integratorObj%addIntegrator(integratorBDE)

            deallocate(integratorBDE)

        case ("CVODE")
            if (allocated(integerParams)) deallocate(integerParams)
            allocate(integerParams(3))

            integerParams(1) = &
            NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEGMRESMaxRestarts,0)
            integerParams(2) = &
            NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEMaxOrder,5)
            integerParams(3) = &
            NamedInteger(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEMaxInternalSteps,500)

            call jsonCont%load(integerParams)
            call jsonCont%output(integerParams)
            
            if (allocated(realParams)) deallocate(realParams)
            allocate(realParams(5))

            realParams(1) = &
            NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyRelTol,1d-5)

            realParams(2) = &
            NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyAbsTol,1d-10)
            realParams(3) = &
            NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEMinStep,0.d0)
            realParams(4) = &
            NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEMaxStep,0.d0)
            realParams(5) = &
            NamedReal(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEInitStep,0.d0)

            call jsonCont%load(realParams)
            call jsonCont%output(realParams)
            
            if (allocated(integerArrayParams)) deallocate(integerArrayParams)
            allocate(integerArrayParams(1))

            integerArrayParams(1) = &
            NamedIntegerArray(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyBBDPreParams,[0,0,0,0])
            
            call jsonCont%load(integerArrayParams)
            call jsonCont%output(integerArrayParams)

            if (allocated(logicalParams)) deallocate(logicalParams)
            allocate(logicalParams(2)) 

            logicalParams(1) = &
            NamedLogical(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEAM,.false.)
            logicalParams(2) = &
            NamedLogical(keyIntegrator//"."//integratorTags(1)%values(i)%string//"."//keyCVODEStabDet,.false.)
            allocate(optionsCVODE)


            call jsonCont%load(logicalParams)
            call jsonCont%output(logicalParams)

            optionsCVODE%maxRestarts = integerParams(1)%value 
            optionsCVODE%reltol = realParams(1)%value 
            optionsCVODE%abstol = realParams(2)%value 
            optionsCVODE%bbdmudq = integerArrayParams(1)%values(1)
            optionsCVODE%bbdmldq = integerArrayParams(1)%values(2)
            optionsCVODE%bbdmukeep= integerArrayParams(1)%values(3)
            optionsCVODE%bbdmlkeep= integerArrayParams(1)%values(4)
            optionsCVODE%maxOrder = integerParams(2)%value
            optionsCVODE%maxInternalSteps = integerParams(3)%value
            optionsCVODE%minTimestep = realParams(3)%value 
            optionsCVODE%maxTimestep = realParams(4)%value 
            optionsCVODE%startingTimestep = realParams(5)%value 
            optionsCVODE%amMethod = logicalParams(1)%value
            optionsCVODE%stabLimitDet = logicalParams(2)%value

            allocate(integratorCVODE)
            call integratorCVODE%init(mpiCont,optionsCVODE,&
                integratorName=integratorTags(1)%values(i)%string)

            call integratorObj%addIntegrator(integratorCVODE)

            deallocate(integratorCVODE)
            deallocate(optionsCVODE)
        case default 
            error stop "Undefined integrator type detected in initStandardIntegrator"
        end select
    end do

    !Initialize and add integration steps 

    allocate(integratorSteps(size(stepTags(1)%values)))

    do i = 1, size(integratorSteps)
        
        !Prepare integrator step parameter storage
        if (allocated(stringArrayParams)) deallocate(stringArrayParams)
        allocate(stringArrayParams(4))

        if (allocated(stringParams)) deallocate(stringParams)
        allocate(stringParams(1))

        stringParams(1) = NamedString(keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyIntegratorTagSingle,"")

        call jsonCont%load(stringParams)

        found = .false. 
        integratorIndex = 0
        do j = 1, size(integratorTags(1)%values)
            if (integratorTags(1)%values(j)%string == stringParams(1)%value) then 
                found = .true.
                integratorIndex = j 
                exit 
            end if
        end do

        if (assertions .or. assertionLvl >= 0) call assert(found,&
        "An integrator tag supplied to an integration step in initStandardIntegrator is invalid")

        integratorSteps(i)%integratorIndex = integratorIndex
        if (allocated(realParams)) deallocate(realParams)
        allocate(realParams(1))

        realParams(1) = NamedReal(keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyGlobalStepFraction,real(1,kind=rk))

        call jsonCont%load(realParams)
        call jsonCont%output(realParams)
        
        integratorSteps(i)%globalStepFraction = realParams(1)%value 

        stringArrayParams(1)%name = keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyEvolvedModels
        allocate(stringArrayParams(1)%values(0))

        stringArrayParams(2)%name = keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyCommData//"."//keyVarsToBroadcast
        allocate(stringArrayParams(2)%values(0))

        stringArrayParams(3)%name = keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyCommData//"."//keyHaloVars
        allocate(stringArrayParams(3)%values(0))

        stringArrayParams(4)%name = keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyCommData//"."//keyScalarsToBroadcast
        allocate(stringArrayParams(4)%values(0))
        if (allocated(integerArrayParams)) deallocate(integerArrayParams)
        allocate(integerArrayParams(1))

        integerArrayParams(1)%name = keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyCommData//"."//keyScalarsRoots
        allocate(integerArrayParams(1)%values(0))

        call jsonCont%load(stringArrayParams)
        call jsonCont%load(integerArrayParams)

        if (size(integerArrayParams(1)%values) == 0 .and. size(stringArrayParams(4)%values) > 0) then 
            integerArrayParams(1)%values = [(0,i=1,size(stringArrayParams(4)%values))]
        end if
        if (size(stringArrayParams(4)%values) > 0) then 
            call assert(size(integerArrayParams(1)%values)==size(stringArrayParams(4)%values),integerArrayParams(1)%name//&
            " size must conform to size of "//stringArrayParams(4)%name)
        end if
        call jsonCont%output(stringArrayParams)

        commNeeded = .false. 
        if (size(stringArrayParams(2)%values) + size(stringArrayParams(3)%values) > 0) then
            commNeeded = .true. 
            if (allocated(commData)) deallocate(commData)
            allocate(CommData)
            commData%varsToBroadcast = stringArrayParams(2)%values
            commData%haloExchangeVars = stringArrayParams(3)%values
            commData%scalarsToBroadcast = stringArrayParams(4)%values
            commData%scalarRoots = integerArrayParams(1)%values 
            
            integratorSteps(i)%commData = commData
        end if
        integratorSteps(i)%communicationNeeded = commNeeded

        if (assertions .or. assertionLvl >= 0) call assert(size(stringArrayParams(1)%values) > 0,&
        "No evolved models found for an integration step &
        &in initStandardIntegrator")

        if (allocated(modelIndices)) deallocate(modelIndices)
        allocate(modelIndices(size(stringArrayParams(1)%values)))

        do j = 1, size(stringArrayParams(1)%values)
            found = .false. 
            do k = 1,size(modelTags(1)%values)
                if (modelTags(1)%values(k)%string == stringArrayParams(1)%values(j)%string) then 
                    found = .true. 
                    modelIndices(j) = k
                    exit
                end if
            end do

            if (assertions .or. assertionLvl >= 0) &
            call assert(found,"No model found for evolved model tag in an integration step in initStandardIntegrator")
        end do

        integratorSteps(i)%modelIndices = modelIndices

        if (allocated(groupIndices)) deallocate(groupIndices)
        allocate(groupIndices(size(modelIndices)))

        if(allocated(updatesOnInternalIterations)) deallocate(updatesOnInternalIterations)
        allocate(updatesOnInternalIterations(size(modelIndices)))

        if (allocated(integerArrayParams)) deallocate(integerArrayParams)
        allocate(integerArrayParams(2*size(stringArrayParams(1)%values)))

        do j = 1, size(stringArrayParams(1)%values)
            integerArrayParams(2*(j-1)+1)%name = keyIntegrator//"."//stepTags(1)%values(i)%string//"."&
                                                //stringArrayParams(1)%values(j)%string//"."//keyGroupIndices
            integerArrayParams(2*(j-1)+1)%values = [1]

            integerArrayParams(2*(j-1)+2)%name = keyIntegrator//"."//stepTags(1)%values(i)%string//"."&
                                                //stringArrayParams(1)%values(j)%string//"."//keyInternalUpdateGroups
            integerArrayParams(2*(j-1)+2)%values = [1]
        end do

        call jsonCont%load(integerArrayParams)
        call jsonCont%output(integerArrayParams)

        if (assertions .or. assertionLvl >= 0) then 
            do j = 1, 2 * size(stringArrayParams(1)%values)
                call assert(all(integerArrayParams(j)%values > 0),&
                "Non-positive group index detected in integration step in initStandardIntegrator")
                call assert(all(integerArrayParams(j)%values <= numImplicitGroups+numGeneralGroups),&
                "Out of bounds group index detected in integration step in initStandardIntegrator")
            end do
        end if

        nonTrivialUpdate = .false. 
        do j = 1, size(stringArrayParams(1)%values)
            groupIndices(j)%entry = integerArrayParams(2*(j-1)+1)%values 
            allocate(updatesOnInternalIterations(j)%entry(size(groupIndices(j)%entry)))
            do k = 1,size(groupIndices(j)%entry)
                updatesOnInternalIterations(j)%entry(k) = any(integerArrayParams(2*(j-1)+2)%values == groupIndices(j)%entry(k))
            end do

            if (any(updatesOnInternalIterations(j)%entry)) nonTrivialUpdate = .true.
        end do

        integratorSteps(i)%groupIndices = groupIndices

        integratorSteps(i)%nonTrivialUpdate = nonTrivialUpdate
        if (nonTrivialUpdate) integratorSteps(i)%updatesOnInternalIterations = updatesOnInternalIterations

        if (allocated(logicalParams)) deallocate(logicalParams)
        allocate(logicalParams(2 + size(stringArrayParams(1)%values)))

        logicalParams(1) = NamedLogical(keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyAllowTimeEvolution,.false.)
        logicalParams(2) = NamedLogical(keyIntegrator//"."//stepTags(1)%values(i)%string//"."//keyUseInitialInput,.false.)

        if (allocated(updatesOnInternalIterationsModelData)) deallocate(updatesOnInternalIterationsModelData)
        allocate(updatesOnInternalIterationsModelData(size(stringArrayParams(1)%values)))
        do j = 1,size(stringArrayParams(1)%values)
            logicalParams(2+j) = NamedLogical(keyIntegrator//"."//stepTags(1)%values(i)%string//"."&
                                            //stringArrayParams(1)%values(j)%string//"."//keyInternalUpdateModelData,.true.)
        end do

        call jsonCont%load(logicalParams)
        call jsonCont%output(logicalParams)

        integratorSteps(i)%allowTimeEvolution = logicalParams(1)%value
        integratorSteps(i)%useInitialInput = logicalParams(2)%value


        do j = 1, size(stringArrayParams(1)%values)
            updatesOnInternalIterationsModelData(j) =  logicalParams(2+j)%value
        end do

        nonTrivialModelDataUpdate = any(updatesOnInternalIterationsModelData) 
        integratorSteps(i)%nonTrivialModelDataUpdate = nonTrivialModelDataUpdate
        if (nonTrivialModelDataUpdate) &
        integratorSteps(i)%updatesOnInternalIterationsModelData = updatesOnInternalIterationsModelData

        call integratorObj%addIntegrationStage(integratorSteps(i))
    end do

end subroutine initStandardIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStandardSignals(signalCollectionObj)
    !! Initialize signal collection with commonly used signals and their names. The signals are:
    !! 1. constant signal
    !! 2. hat signal
    !! 3. cut sine signal

    type(SignalCollection) ,intent(inout) :: signalCollectionObj 

    type(ConstSignal)   :: cSignal
    type(HatSignal)     :: hSignal
    type(CutSineSignal) :: csSignal

    call signalCollectionObj%init()

    call signalCollectionObj%addSignal(cSignal,keyConstSignal)
    call signalCollectionObj%addSignal(hSignal,keyHatSignal)
    call signalCollectionObj%addSignal(csSignal,keyCutSineSignal)

    end subroutine initStandardSignals
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setDistHarmonic(dist,h,xArr,vArr)
    !! Set the h-th harmonic of distribution function variable dist to xArr*vArr at corresponding x, and v coordinates. 
    !! vArr should be either rank 1 or 2, and if it is rank 2 it should be size(numV,size(xArr))

    real(rk) ,dimension(:)  ,intent(inout) :: dist 
    integer(ik)             ,intent(in)    :: h
    real(rk) ,dimension(:)  ,intent(in)    :: xArr 
    real(rk) ,dimension(..) ,intent(in)    :: vArr

    real(rk) ,dimension(:)   ,allocatable :: vVec 
    real(rk) ,dimension(:,:) ,allocatable :: vMat

    integer(ik) :: locNumX ,numH ,numV ,i ,indOffset

    select rank(vArr)
    rank (1)
        vVec = vArr 
        numV = size(vVec)
    rank (2)
        vMat = vArr 
        numV = size(vMat,1)
    rank default 
        error stop "Unsupported rank vArr passed to setDistHarmonic"
    end select 

    locNumX = size(xArr)
    numH = size(dist)/(locNumX*numV)

    if (assertions .or. assertionLvl >= 0) then 
        call assert(numH >= h,"Harmonic index passed to setDistHarmonic out of bounds - upper")
        call assert(h >= 1,"Harmonic index passed to setDistHarmonic out of bounds - lower")

        if (allocated(vVec)) call assert(size(vVec) == numV,"vArr must have size numV if rank 1")
        if (allocated(vMat)) then 
            call assert(size(vMat,1) == numV,"If rank of vArr is 2, the first dimension must have size numV")
            call assert(size(vMat,2) == locNumX,"If rank of vArr is 2, the second dimension must conform to&
                        & size of xArr")
        end if
    end if

    do i = 1,locNumX
        indOffset = (i-1)*numV*numH + (h-1) * numV 
        if (allocated(vMat)) then 
            dist(indOffset + 1:indOffset+numV) = xArr(i)*vMat(:,i)
        else
            dist(indOffset + 1:indOffset+numV) = xArr(i)*vVec
        end if
    end do
end subroutine setDistHarmonic
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule initialization_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
