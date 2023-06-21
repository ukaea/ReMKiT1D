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
submodule (general_mat_term_class) general_mat_term_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the general matrix term

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initGeneralTerm(this,gridObj,partitionObj,indexingObj,procRank,&
    evolvedVar,implicitVar,varCont,sTemplate,vData,normConst,&
    coordProfile,tData,evalTermGroup,mbData,copyTermName) 
    !! General matrix term initialization routine 

    class(GeneralMatrixTerm)                    ,intent(inout)  :: this
    type(Grid)                                  ,intent(in)     :: gridObj !! Grid object used to initialize stencil
    type(Partition)                             ,intent(in)     :: partitionObj !! Parition object used to determine local number of DoF
    type(Indexing)                              ,intent(in)     :: indexingObj !! Indexing object used to initialize matrix row data
    integer(ik)                                 ,intent(in)     :: procRank !! Current processor rank
    character(*)                                ,intent(in)     :: evolvedVar !! Name of evolved variable
    character(*)                                ,intent(in)     :: implicitVar !! Name of implicit variable
    type(VariableContainer)                     ,intent(in)     :: varCont !! Reference variable container
    type(StencilTemplate)                       ,intent(inout)  :: sTemplate !! StencilTemplate used 
    type(VarData)       ,optional               ,intent(in)     :: vData !! Required variable names and powers
    real(rk)            ,optional               ,intent(in)     :: normConst  !! Normalization constant 
    type(CoordProfiles) ,optional               ,intent(in)     :: coordProfile !! Multiplicative coordinate profiles (globally indexed)
    type(TimeSignalData) ,optional              ,intent(in)     :: tData !! Explicit time dependence 
    integer(ik)       ,optional                 ,intent(in)     :: evalTermGroup !! Optional term group of host model to evaluate and use as additional row variable
    class(ModelboundData) ,optional             ,intent(in)     :: mbData !! Optional modelbound data object used when initializing fixed stancil generators which might require it
    character(*) ,optional                      ,intent(in)     :: copyTermName !! Name of the term whose values should be copied and multiplied with this term's stencil

    integer(ik) :: minX, maxX ,i

    type(StringArray) ,allocatable ,dimension(:) :: reqVarsRow ,reqVarsCol
    integer(ik)     ,allocatable ,dimension(:,:) :: rowCoordsUsed

    logical :: nonTrivialRowFun ,nonTrivialColFun


    if (assertions) then 
        call assert(gridObj%isDefined(),"Undefined grid object passed to general matrix term constructor")
        call assert(partitionObj%isDefined(),"Undefined partition object passed to general matrix term constructor")
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to general matrix term constructor")
        call assert(varCont%isDefined(),"Undefined variable container object passed to general matrix term constructor")
    end if 

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    this%locNumX = maxX - minX + 1
    this%numV = gridObj%getNumV()
    this%numH = gridObj%getNumH()

    this%kineticRow = .false. 
    this%kineticCol = .false.

    nonTrivialRowFun = .false.
    nonTrivialColFun = .false.

    if (assertions) then 
        if (present(vData)) then
            if (allocated(vData%rowVars) .and. allocated(vData%rowVarPowers)) &
            call assert(size(vData%rowVars) == size(vData%rowVarPowers),&
            "If row var powers are passed to general matrix term constructor they must conform to rowVars array") 

            if (allocated(vData%colVars) .and. allocated(vData%colVarPowers)) &
            call assert(size(vData%colVars) == size(vData%colVarPowers),&
            "If col var powers are passed to general matrix term constructor they must conform to colVars array") 

            if (allocated(vData%modelboundRowVars) .and. allocated(vData%modelboundRowVarPowers)) &
            call assert(size(vData%modelboundRowVars) == size(vData%modelboundRowVarPowers),&
            "If modelbound row var powers are passed to general matrix term constructor&
            & they must conform to modelboundRowVars array") 

            if (allocated(vData%modelboundColVars) .and. allocated(vData%modelboundColVarPowers)) &
            call assert(size(vData%modelboundColVars) == size(vData%modelboundColVarPowers),&
            "If modelbound col var powers are passed to general matrix term constructor&
            & they must conform to modelboundColVars array") 
        end if

        if (present(coordProfile)) then
            if (allocated(coordProfile%xProfile)) &
            call assert(size(coordProfile%xProfile) == gridObj%getNumX(),&
            "xProfile passed to general matrix term must have same size as x-grid")

            if (allocated(coordProfile%hProfile)) &
            call assert(size(coordProfile%hProfile) == gridObj%getNumH(),&
            "hProfile passed to general matrix term must conform to number of harmonics")

            if (allocated(coordProfile%vProfile)) &
            call assert(size(coordProfile%vProfile) == gridObj%getNumV(),&
            "vProfile passed to general matrix term must have same size as v-grid")

        end if

        if (present(tData)) then 

            call assert(allocated(tData%tSignal),&
            "Signal object in time signal data passed to general matrix term not allocated")
            call assert(allocated(tData%tParams),&
            "Signal params in time signal data passed to general matrix term not allocated")

        end if

        if (allocated(sTemplate%stencilGen)) &
        call assert(sTemplate%stencilGen%isDefined(),&
        "Undefined stencil value generator passed to general matrix term constructor")

        if (allocated(sTemplate%overridingStencilCoords) .or. allocated(sTemplate%overridingStencils)) then 
            call assert(allocated(sTemplate%overridingStencilCoords) .eqv. allocated(sTemplate%overridingStencils),&
            "If overridingStencilCoords are&
            & passed to general matrix term constructor so must be overridingStencils in the stencil template, and vice versa")
            call assert(size(sTemplate%overridingStencils) == size(sTemplate%overridingStencilCoords,2),"overridingStencil and &
            &overridingStencilCoords passed to general matrix term constructor must be of same size")
        end if

        call assert(allocated(sTemplate%rowCoords),"rowCoords not allocated in sTemplate passed to general matrix term constructor")

        call assert(any(size(sTemplate%rowCoords,1) == [1,3]),"rowCoords in sTemplate passed to general matrix term constructor&
        & must be (1,:) or (3,:)")

    end if  

    this%kineticRow = size(sTemplate%rowCoords,1) == 3
    this%kineticCol = varCont%isVarDist(varCont%getVarIndex(implicitVar))

    this%isActive = partitionObj%getMinHAtInd(procRank+1) == 1
    if (this%kineticRow) this%isActive = .true. 

    if (this%isActive) then 
        rowCoordsUsed = partitionObj%filterCoords(procRank+1,sTemplate%rowCoords)

    else
        allocate(rowCoordsUsed(1,0))
    end if

    this%isActive = size(rowCoordsUsed,2) > 0

    allocate(reqVarsRow(0))
    allocate(reqVarsCol(0))

    if (present(vData)) then 

        this%vData = vData
        if (allocated(this%vData%rowVars)) then 
            if (size(this%vData%rowVars) > 0) then
                reqVarsRow = [reqVarsRow,this%vData%rowVars]
                nonTrivialRowFun = .true.
                if (.not. allocated(this%vData%rowVarPowers)) then 
                    allocate(this%vData%rowVarPowers(size(reqVarsRow)))
                    this%vData%rowVarPowers = real(1,kind=rk)
                end if
            end if
        end if

        if (allocated(this%vData%colVars)) then 
            if (size(this%vData%colVars) > 0) then
                reqVarsCol = [reqVarsCol,this%vData%colVars]
                nonTrivialColFun = .true.
                if (.not. allocated(this%vData%colVarPowers)) then 
                    allocate(this%vData%colVarPowers(size(reqVarsCol)))
                    this%vData%colVarPowers = real(1,kind=rk)
                end if
            end if
        end if

        if (allocated(this%vData%modelboundRowVars)) then 
            if (size(this%vData%modelboundRowVars) > 0) then
                nonTrivialRowFun = .true.
                allocate(this%modelboundRowVarBuffer,source=varCont%variables(varCont%getVarIndex(evolvedVar))%entry)
                if (.not. allocated(this%vData%modelboundRowVarPowers)) then 
                    allocate(this%vData%modelboundRowVarPowers(size(this%vData%modelboundRowVars)))
                    this%vData%modelboundRowVarPowers = real(1,kind=rk)
                end if
            end if
        end if

        if (allocated(this%vData%modelboundColVars)) then 
            if (size(this%vData%modelboundColVars) > 0) then 
                nonTrivialColFun = .true.
                allocate(this%modelboundColVarBuffer,source=varCont%variables(varCont%getVarIndex(implicitVar))%entry)
                if (.not. allocated(this%vData%modelboundColVarPowers)) then 
                    allocate(this%vData%modelboundColVarPowers(size(this%vData%modelboundColVars)))
                    this%vData%modelboundColVarPowers = real(1,kind=rk)
                end if
            end if
        end if

    end if

    allocate(this%reqRowVarIsDist(size(reqVarsRow)))
    allocate(this%reqColVarIsDist(size(reqVarsCol)))

    do i = 1, size(reqVarsRow)
        this%reqRowVarIsDist(i) = varCont%isVarDist(varCont%getVarIndex(reqVarsRow(i)%string))
    end do

    do i = 1, size(reqVarsCol)
        this%reqColVarIsDist(i) = varCont%isVarDist(varCont%getVarIndex(reqVarsCol(i)%string))
    end do

    !Check to make sure required variables conform
    if (assertions) then 
        if (.not. this%kineticRow) call assert(.not. any(this%reqRowVarIsDist),&
        "Required row variable in general matrix term is a distribution when the evolved variable isn't")

        if (.not. this%kineticCol) call assert(.not. any(this%reqColVarIsDist),&
        "Required column variable in general matrix term is a distribution when the implicit variable isn't")
    end if

    call this%setEvolvedAndImplicitVar(evolvedVar,implicitVar,varCont)
    call this%setReqVars(reqVarsRow,reqVarsCol,varCont)

    call this%initRowData(rowCoordsUsed,stencilMapping,indexingObj) 

    allocate(this%coordProfile%xProfile(this%locNumX))
    this%coordProfile%xProfile = real(1,kind=rk)
    allocate(this%coordProfile%hProfile(gridObj%getNumH()))
    this%coordProfile%hProfile = real(1,kind=rk)
    allocate(this%coordProfile%vProfile(this%numV))
    this%coordProfile%vProfile = real(1,kind=rk)
    if (present(coordProfile)) then 
        nonTrivialRowFun = .true.
        if (allocated(coordProfile%xProfile)) this%coordProfile%xProfile = coordProfile%xProfile(minX:maxX)
        if (allocated(coordProfile%hProfile)) this%coordProfile%hProfile = coordProfile%hProfile
        if (allocated(coordProfile%vProfile)) this%coordProfile%vProfile = coordProfile%vProfile

    end if
    if (present(tData)) then
        this%tdata = tData
        nonTrivialRowFun = .true.
    end if
    this%timeSignalMult = real(1,kind=rk)
    call this%makeDefined()

    if (this%isActive) then 
        if (present(normConst)) call this%setNormalizationConst(normConst)

        if (allocated(sTemplate%stencilGen)) then 
            if (sTemplate%fixedStencil) then
                if (present(mbData)) then
                    call this%setMultConst(sTemplate%stencilGen%calculate(varCont,mbData=mbData))
                else
                    call this%setMultConst(sTemplate%stencilGen%calculate(varCont))
                end if
            else
                allocate(this%stencilGen,source=sTemplate%stencilGen)
            end if

        end if
        this%stencilVals = this%getMultConst()

    end if

    if (present(evalTermGroup)) then 
        nonTrivialRowFun = .true.
        allocate(this%targetTermGroup,source=evalTermGroup)
        allocate(this%evalBuffer,source=varCont%variables(varCont%getVarIndex(evolvedVar))%entry)
        this%evalBuffer = real(1,kind=rk)
    end if

    call this%setNonTrivialColFun(nonTrivialColFun)
    call this%setNonTrivialRowFun(nonTrivialRowFun)

    if (present(copyTermName)) this%copyTermName = copyTermName

    contains 

    pure function stencilMapping(inputArray) result(output)

            integer(ik)    ,dimension(:)        ,intent(in) :: inputArray 
            type(intArray) ,allocatable ,dimension(:)       :: output

            integer(ik) :: j
            logical :: overridenStencil 

            overridenStencil = .false.
            if (allocated(sTemplate%overridingStencilCoords)) then

                do j = 1,size(sTemplate%overridingStencilCoords,2)
                    if (all(inputArray == sTemplate%overridingStencilCoords(:,j))) then 
                        output = sTemplate%overridingStencils(j)%mapCoords(gridObj,inputArray)
                        overridenStencil = .true.
                        exit
                    end if
                end do
            end if
            if (.not. overridenStencil) output = sTemplate%defaultStencil%mapCoords(gridObj,inputArray)

    end function stencilMapping

end subroutine initGeneralTerm
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateGeneralTerm(this,varCont,indexingData,modelData,hostModel) 
    !! Update general matrix term if using any modelbound data, updating stencil values, evaluating other terms, or if there is any explicit time dependence 

    class(GeneralMatrixTerm)        ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container - used to update stencil values and signal
    type(MatrixTermIndexingData)    ,intent(in)     :: indexingData
    class(ModelboundData) ,optional ,intent(in)     :: modelData !! Model data used to retrieve modelbound variable
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - used to evaluate term group

    real(rk)        ,allocatable ,dimension(:)      :: modelboundDataVals
    integer(ik)                                     :: i, j ,k ,mbHalo ,lboundBuffer ,mbDataDim ,bufferHalo ,haloDiff ,mbLBound
    real(rk) :: timeVal

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update undefined general matrix term object")
        call assert(varCont%isDefined(),"Attempted to update general matrix term object using undefined variable container")

        if (allocated(this%modelboundColVarBuffer) .or. allocated(this%modelboundRowVarBuffer))&
        call assert(present(modelData),&
        "Attempted to updated general matrix term object with modelbound dependencies without passing modelData")

        if (allocated(this%targetTermGroup))&
        call assert(present(hostModel),"Attempted to update general matrix term with group evaluation without passing hostModel")
    end if
    if (this%isActive) then 

        if (allocated(this%modelboundRowVarBuffer)) then
            this%modelboundRowVarBuffer = real(1,kind=rk)
            !Determine buffer halo size
            lboundBuffer = lbound(this%modelboundRowVarBuffer,1)
            bufferHalo = 1 - lboundBuffer
            if (this%kineticRow) bufferHalo = bufferHalo/(this%numH*this%numV)
            do i = 1, size(this%vData%modelboundRowVars)

                !NOTE:  gfortran select rank bug workaround
                if (.not. allocated(modelboundDataVals)) allocate(modelboundDataVals(0))

                call modelData%copyData(this%vData%modelboundRowVars(i)%string,modelboundDataVals)
                mbDataDim = modelData%getDataDim(this%vData%modelboundRowVars(i)%string)
                mbLBound = lbound(modelboundDataVals,1)
                select case (mbDataDim)
                case (0)
                    this%modelboundRowVarBuffer = this%modelboundRowVarBuffer * &
                    modelboundDataVals(1)**this%vData%modelboundRowVarPowers(i)
                case (1) 
                    mbHalo = (size(modelboundDataVals) - this%locNumX)/2
                    haloDiff = bufferHalo - mbHalo
                    if (this%kineticRow) then 

                        do j = 1, this%locNumX
                            this%modelboundRowVarBuffer((j-1)*this%numH*this%numV+1:j*this%numH*this%numV) = &
                            this%modelboundRowVarBuffer((j-1)*this%numH*this%numV+1:j*this%numH*this%numV) &
                            * modelboundDataVals(j-1+mbLBound+mbHalo)**this%vData%modelboundRowVarPowers(i)
                        end do

                    else
                        this%modelboundRowVarBuffer(lboundBuffer+haloDiff:lboundBuffer&
                        +haloDiff+size(modelboundDataVals)-1) = &
                        this%modelboundRowVarBuffer(lboundBuffer+haloDiff:lboundBuffer&
                        +haloDiff+size(modelboundDataVals)-1) * modelboundDataVals**this%vData%modelboundRowVarPowers(i)
                    end if
                case (2) 
                    mbHalo = (size(modelboundDataVals) - this%locNumX*this%numV)/(2*this%numV)
                    if (.not. this%kineticRow) &
                    error stop "Single harmonic modelbound data detected in general matrix term when evolved variable is fluid"
                    do j = 1, this%locNumX
                        do k = 1,this%numH
                            this%modelboundRowVarBuffer((j-1)*this%numH*this%numV+(k-1)*this%numV + 1 &
                                                       :(j-1)*this%numH*this%numV+k*this%numV) = &
                            this%modelboundRowVarBuffer((j-1)*this%numH*this%numV+(k-1)*this%numV + 1 &
                                                       :(j-1)*this%numH*this%numV+k*this%numV) &
                            * modelboundDataVals(mbLBound+(j+mbHalo-1)*this%numV:mbLBound+(j+mbHalo)*this%numV-1)&
                            **this%vData%modelboundRowVarPowers(i)
                        end do
                    end do
                case (3)
                    mbHalo = (size(modelboundDataVals) - this%locNumX*this%numH*this%numV)/2
                    haloDiff = bufferHalo*this%numH*this%numV - mbHalo
                    if (.not. this%kineticRow) &
                    error stop "Distribution modelbound data detected in general matrix term when evolved variable is fluid"

                    this%modelboundRowVarBuffer(lboundBuffer+haloDiff:lboundBuffer&
                        +haloDiff+size(modelboundDataVals)-1) = &
                        this%modelboundRowVarBuffer(lboundBuffer+haloDiff:lboundBuffer&
                        +haloDiff+size(modelboundDataVals)-1) * modelboundDataVals**this%vData%modelboundRowVarPowers(i)

                case default
                    error stop "Unsupported dimensionality of rank 1 modelbound data detected in general matrix term"
                end select
            end do
        end if

        if (allocated(this%modelboundColVarBuffer)) then
            this%modelboundColVarBuffer = real(1,kind=rk)

            do i = 1, size(this%vData%modelboundColVars)
                call modelData%copyData(this%vData%modelboundColVars(i)%string,modelboundDataVals)
                mbDataDim = modelData%getDataDim(this%vData%modelboundRowVars(i)%string)

                if (size(modelboundDataVals) /= size(this%modelboundColVarBuffer)) &
                    error stop "modelbound data in general matrix column buffer is not of expected size"
                
                this%modelboundColVarBuffer = &
                this%modelboundColVarBuffer * modelboundDataVals**this%vData%modelboundColVarPowers(i)
                
            end do
        end if

        if (allocated(this%tData)) then 
            timeVal = varCont%variables(varCont%getVarIndex("time"))%entry(1)
            this%timeSignalMult = this%tData%tSignal%calculate(timeVal,this%tData%tPeriod,this%tData%tParams)
        end if

        if (allocated(this%targetTermGroup)) then 
            if (this%targetTermGroup >0) then
                select type (hostModel)
                type is (Model)
                    this%evalBuffer = hostModel%evaluateTermGroup(this%targetTermGroup,varCont)
                class default 
                    error stop "Unsupported ModelSurrogate passed to updateGeneralTerm"
                end select
            end if
        end if

        if (allocated(this%copyTermName)) then

            if (assertions) then
                 call assert(present(hostModel),&
                 "hostModel must be present if copyTermName is given to GeneralMatrixTerm")

                select type (hostModel)
                type is (Model)
                    call assert(hostModel%isTermNameRegistered(this%copyTermName),&
                    "copyTermName "//this%copyTermName//"  not registered in host model passed to general term update routine")
                class default 
                    error stop "Unsupported ModelSurrogate passed to updateGeneralTerm"
                end select
            end if

            select type (hostModel)
            type is (Model)
                this%copyTermVals = hostModel%getImplicitTermRowData(hostModel%getImplicitTermIndex(this%copyTermName))
            class default 
                error stop "Unsupported ModelSurrogate passed to updateGeneralTerm"
            end select

        end if

        if (allocated(this%stencilGen)) then 
            if (present(modelData)) then 
                if (present(hostModel)) then
                    call this%stencilGen%calculateInPlace(varCont,this%stencilVals,modelData,hostModel)
                else
                    call this%stencilGen%calculateInPlace(varCont,this%stencilVals,modelData)
                end if
            else
                if (present(hostModel)) then
                    call this%stencilGen%calculateInPlace(varCont,this%stencilVals,hostModel=hostModel)
                else
                    call this%stencilGen%calculateInPlace(varCont,this%stencilVals)
                end if
            end if
        end if


        if (allocated(this%stencilGen) .or. allocated(this%copyTermName)) then

            if (.not. allocated(this%copyTermName)) then
                call this%setMultConst(this%stencilVals)
            else 
                call this%setMultConst(this%stencilVals * this%copyTermVals)

            end if

        end if

    end if


end subroutine updateGeneralTerm
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine generalRow(this,varCont,rowVals,indexingData) 
!! Multiplies all required variables and raises them to their corresponding power,using modelbound row variable buffer if allocated, and adding 
!! spatial and time profiles

    class(GeneralMatrixTerm)      ,intent(inout)   :: this
    type(VariableContainer)       ,intent(in)      :: varCont
    real(rk) ,dimension(:)        ,intent(inout)   :: rowVals
    type(MatrixTermIndexingData)  ,intent(in)      :: indexingData

    integer(ik)                                    :: i 

    if (.not. allocated(this%rowCoordsX)) then 

        allocate(this%rowCoordsX(size(rowVals)))

        do i = 1,size(rowVals)
            this%rowCoordsX(i) = indexingData%rowDataCoordsLocal(i)%rowCoords(1)
        end do

        if (this%kineticRow) then 
            allocate(this%rowCoordsH(size(rowVals)))
            allocate(this%rowCoordsV(size(rowVals)))


            do i = 1,size(rowVals)
                this%rowCoordsH(i) = indexingData%rowDataCoordsLocal(i)%rowCoords(2)
                this%rowCoordsV(i) = indexingData%rowDataCoordsLocal(i)%rowCoords(3)
            end do
        end if

    end if
    rowVals = this%coordProfile%xProfile(this%rowCoordsX) * this%timeSignalMult

    if (this%kineticRow) rowVals = rowVals * this%coordProfile%hProfile(this%rowCoordsH) * &
                                    this%coordProfile%vProfile(this%rowCoordsV)

    do i = 1, size(indexingData%rowReqVarIndices)
        if (varCont%isVarScalar(indexingData%rowReqVarIndices(i))) then 
            rowVals = rowVals * varCont%variables(indexingData%rowReqVarIndices(i))%entry(1) &
            ** this%vData%rowVarPowers(i)
        else if (this%kineticRow .and. .not. this%reqRowVarIsDist(i)) then 
            rowVals = rowVals * varCont%variables(indexingData%rowReqVarIndices(i))%entry&
            (this%rowCoordsX) &
            ** this%vData%rowVarPowers(i)
        else
            rowVals = rowVals * varCont%variables(indexingData%rowReqVarIndices(i))%entry(indexingData%localRowIndices) &
            ** this%vData%rowVarPowers(i)
        end if
    end do

    ! Assumes that the evaluated group variable and modelbound buffer conform to row variable
    if (allocated(this%modelboundRowVarBuffer)) rowVals = rowVals * this%modelboundRowVarBuffer(indexingData%localRowIndices)

    if (allocated(this%targetTermGroup)) rowVals = rowVals *this%evalBuffer(indexingData%localRowIndices)

end subroutine generalRow 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine generalCol(this,varCont,colVals,indexingData) 
    !! Multiplies all required variables and raises them to their corresponding power,using modelbound variable buffer if allocated

    class(GeneralMatrixTerm)      ,intent(inout) :: this
    type(VariableContainer)       ,intent(in)    :: varCont
    type(RealArray) ,dimension(:) ,intent(inout) :: colVals
    type(MatrixTermIndexingData)  ,intent(in)    :: indexingData

    integer(ik)                                 :: i ,j ,k

    do j = 1,size(colVals)
        colVals(j)%entry = real(1,kind=rk)
    end do 
    
    do i = 1, size(indexingData%colReqVarIndices)
        if (this%kineticCol .and. .not. this%reqColVarIsDist(i)) then 
            do j = 1,size(colVals)
                do k = 1,size(colVals(j)%entry) 
                    colVals(j)%entry(k) = colVals(j)%entry(k) * varCont%variables(indexingData%colReqVarIndices(i))%entry&
                    (indexingData%rowDataCoordsLocal(j)%colCoords(1,k)) &
                    ** this%vData%colVarPowers(i)
                end do
            end do
        else
            do j = 1,size(colVals)
                colVals(j)%entry = colVals(j)%entry &
                * varCont%variables(indexingData%colReqVarIndices(i))%entry(indexingData%localColIndices(j)%entry) &
                ** this%vData%colVarPowers(i)
            end do
        end if
    end do

    if (allocated(this%modelboundColVarBuffer)) then 
        do j = 1,size(colVals)
            colVals(j)%entry = colVals(j)%entry * this%modelboundColVarBuffer(indexingData%localColIndices(j)%entry)
        end do
    end if
end subroutine generalCol 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule general_mat_term_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
