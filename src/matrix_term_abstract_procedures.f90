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
submodule (matrix_term_abstract_class) matrix_term_abstract_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the abstract matrix term class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine unityRow(this,varCont,rowVals,indexingData) 
    !! Default row function. In general should use the passed variable container, row buffer, and matrix indexing data
    !! to provide values which is only a function of the row 

    class(MatrixTerm)             ,intent(inout)   :: this
    type(VariableContainer)       ,intent(in)      :: varCont
    real(rk) ,dimension(:)        ,intent(inout)   :: rowVals
    type(MatrixTermIndexingData)  ,intent(in)      :: indexingData

    if (assertions) then
        call assertPure(this%isDefined(),"Called unityRow with undefined matrix term object")
        call assertPure(varCont%isDefined(),"Called unityRow by passing undefined variable container")
    end if

    rowVals = real(1,kind=rk)

end subroutine unityRow 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine unityCol(this,varCont,colVals,indexingData) 
    !! Default column function. In general should provide values for each column 

    class(MatrixTerm)             ,intent(inout) :: this
    type(VariableContainer)       ,intent(in)    :: varCont
    type(RealArray) ,dimension(:) ,intent(inout) :: colVals
    type(MatrixTermIndexingData)  ,intent(in)    :: indexingData

    integer(ik) :: i 

    if (assertions) then
        call assertPure(this%isDefined(),"Called unityCol with undefined matrix term object")
        call assertPure(varCont%isDefined(),"Called unityCol by passing undefined variable container")
    end if

    do i = 1,size(colVals)
        colVals(i)%entry = real(1,kind=rk)
    end do

end subroutine unityCol  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine calculateValues(this,varCont) 
    !! Default matrix entry calculation routine - loops over each row and applies calculateRows and colFun to variables in provided container

    class(MatrixTerm)       ,intent(inout)  :: this
    type(VariableContainer) ,intent(in)     :: varCont

    integer(ik) :: i

    if (assertions) then
        call assertPure(this%isDefined(),"Attempted to calculate values of undefined matrix term object")
        call assertPure(varCont%isDefined(),&
        "Attempted to calculate values of matrix term object by passing undefined variable container")
    end if

    if (.not. this%fixedMatrixCalculated) then

        if (this%hasNonTrivialRowFun) then
            if (this%hasNonTrivialColFun) then 

                call this%calculateRows(varCont,this%rowBuffer,this%indexingData)
                call this%calculateCols(varCont,this%colBuffer,this%indexingData)
                do i = 1,size(this%rowData%values)
                    this%rowData%values(i)%entry = this%normalizationConst * this%multConst(i)%entry &
                                                * this%rowBuffer(i) &
                                                * this%colBuffer(i)%entry
                end do

            else

                call this%calculateRows(varCont,this%rowBuffer,this%indexingData)
                do i = 1,size(this%rowData%values)
                    this%rowData%values(i)%entry = this%normalizationConst * this%multConst(i)%entry &
                                                * this%rowBuffer(i) 
                end do

            end if

        else

            if (this%hasNonTrivialColFun) then

                call this%calculateCols(varCont,this%colBuffer,this%indexingData)
                do i = 1,size(this%rowData%values)
                    this%rowData%values(i)%entry = this%normalizationConst * this%multConst(i)%entry &
                                                * this%colBuffer(i)%entry
                end do

            else

                do i = 1,size(this%rowData%values)
                    this%rowData%values(i)%entry = this%normalizationConst * this%multConst(i)%entry 

                end do

            end if

        end if

    end if

    if (this%fixedMatrix) this%fixedMatrixCalculated = .true.

end subroutine calculateValues
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNormalizationConst(this,norm) 
    !! Setter for normalizationConst

    class(MatrixTerm)       ,intent(inout)  :: this
    real(rk)                ,intent(in)     :: norm

    this%normalizationConst = norm

end subroutine setNormalizationConst
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNormalizationConst (this) result(norm)
    !! Getter for normalizationConst

    class(MatrixTerm)  ,intent(in) :: this
    real(rk)                       :: norm

    if (assertions) call assertPure(this%isDefined(),"Requested normalization constant from undefined matrix term object")
    norm = this%normalizationConst

end function getNormalizationConst
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setReqVars(this,rowReqVars,colReqVars,varCont) 
    !! Set variable names required by the row and column functions and find their indices in variable container

    class(MatrixTerm)               ,intent(inout)  :: this
    type(StringArray) ,dimension(:) ,intent(in)     :: rowReqVars
    type(StringArray) ,dimension(:) ,intent(in)     :: colReqVars
    type(VariableContainer)         ,intent(in)     :: varCont

    integer(ik) :: i

    if (assertions) then 
        call assertPure(varCont%isDefined(),"Called setReqVars for matrix term by passing undefined variable container")
    end if

    allocate(this%indexingData%rowReqVarIndices(size(rowReqVars)))
    allocate(this%indexingData%colReqVarIndices(size(colReqVars)))

    do i = 1, size(rowReqVars)
        this%indexingData%rowReqVarIndices(i) = varCont%getVarIndex(rowReqVars(i)%string)
    end do

    do i = 1, size(colReqVars)
        this%indexingData%colReqVarIndices(i) = varCont%getVarIndex(colReqVars(i)%string)
    end do
end subroutine setReqVars
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setEvolvedAndImplicitVar(this,rowVarName,colVarName,varCont) 
    !! Set evolved (row) and implicit (column) variable names, and check if evolved variable is stationary using varCont

    class(MatrixTerm)         ,intent(inout)  :: this
    character(*)              ,intent(in)     :: rowVarName
    character(*)              ,intent(in)     :: colVarName
    type(VariableContainer)   ,intent(in)     :: varCont

    this%indexingData%rowVarName = rowVarName
    this%indexingData%colVarName = colVarName 

    this%stationaryEvolvedVar = varCont%isStationary(rowVarName)

end subroutine setEvolvedAndImplicitVar
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initRowData(this,rowCoords,rowToColMapping,indexingObj) 
    !! Initialize row data objects based on a set of evolved global coordinates, the row/col var name, and a function that returns 
    !! column coordinates given a row coordinate input. Requires a reference Indexing object.

    class(MatrixTerm)            ,intent(inout)  :: this
    integer(ik)  ,dimension(:,:) ,intent(in)     :: rowCoords
    procedure(coordMapping)                      :: rowToColMapping
    type(Indexing)               ,intent(in)     :: indexingObj

    integer(ik) :: i ,j

    type(IntArray) ,allocatable ,dimension(:) :: mappedColIndices ,mappedColIndicesLocal
    integer(ik) ,allocatable ,dimension(:,:)  :: allCombinationsMappedInd 

    if (assertions) then 
        call assert(indexingObj%isDefined(),"initRowData for matrix term called with undefined indexing object")
        call assert((size(rowCoords,1) == 1) .or. (size(rowCoords,1) == 3),"Row coordinate entries passed to initRowData must&
        & be of dimension (1,:) or (3,:) depending on whether the row variable is a distribution")
    end if

    call this%rowData%init()
    allocate(this%indexingData%rowDataCoordsGlobal(size(rowCoords,2)))
    allocate(this%indexingData%rowDataCoordsLocal(size(rowCoords,2)))
    allocate(this%multConst(size(rowCoords,2)))
    allocate(this%colBuffer(size(rowCoords,2)))
    allocate(this%rowBuffer(size(rowCoords,2)))
    allocate(this%indexingData%localRowIndices(size(rowCoords,2)))
    allocate(this%indexingData%localColIndices(size(rowCoords,2)))

    this%rowBuffer = real(1,kind=rk)
    if (size(rowCoords,1) == 1) then
        do i = 1,size(rowCoords,2)
            mappedColIndices = rowToColMapping(rowCoords(:,i))
            mappedColIndicesLocal = mappedColIndices
            allCombinationsMappedInd = allCombinations(mappedColIndices)
            call this%rowData%addRow(indexingObj%findIndex(this%indexingData%rowVarName,rowCoords(1,i)),&
                                     indexingObj%mapToGlobalIndices(this%indexingData%colVarName,mappedColIndices))

            this%indexingData%rowDataCoordsGlobal(i)%rowCoords = rowCoords(:,i)
            this%indexingData%localRowIndices(i) = indexingObj%findLocalXIndex(rowCoords(1,i))

            this%indexingData%rowDataCoordsLocal(i)%rowCoords = rowCoords(:,i) 
            this%indexingData%rowDataCoordsLocal(i)%rowCoords(1) = indexingObj%findLocalXIndex(rowCoords(1,i))
            this%indexingData%rowDataCoordsGlobal(i)%colCoords = allCombinationsMappedInd

            do j = 1, size(mappedColIndicesLocal(1)%entry)
                mappedColIndicesLocal(1)%entry(j) = indexingObj%findLocalXIndex(mappedColIndicesLocal(1)%entry(j),&
                                                                                locXInd=rowCoords(1,i))
            end do

            this%indexingData%rowDataCoordsLocal(i)%colCoords = allCombinations(mappedColIndicesLocal)

            allocate(this%multConst(i)%entry(size(this%rowData%values(i)%entry)))
            allocate(this%colBuffer(i)%entry(size(this%rowData%values(i)%entry)))
            this%multConst(i)%entry = real(1.0d00,kind=rk)
            this%colBuffer(i)%entry = real(1.0d00,kind=rk)

            allocate(this%indexingData%localColIndices(i)%entry(size(this%rowData%values(i)%entry)))

            if (size(mappedColIndices) == 1) then
                do j = 1, size(this%indexingData%localColIndices(i)%entry)
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findLocalXIndex(allCombinationsMappedInd(1,j),&
                                                                                                locXInd=rowCoords(1,i))
                end do
            else 
                do j = 1, size(this%indexingData%localColIndices(i)%entry)
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findDistIndex(allCombinationsMappedInd(1,j),&
                                                                                 allCombinationsMappedInd(2,j),&
                                                                                 allCombinationsMappedInd(3,j),&
                                                                                 .true.,&
                                                                                 locXInd=rowCoords(1,i))
                end do

            end if

        end do
    else 

        do i = 1,size(rowCoords,2)
            mappedColIndices = rowToColMapping(rowCoords(:,i))
            mappedColIndicesLocal = mappedColIndices
            allCombinationsMappedInd = allCombinations(mappedColIndices)
            call this%rowData%addRow(indexingObj%findIndex(this%indexingData%rowVarName,rowCoords(1,i),&
                                                                           rowCoords(2,i),&
                                                                           rowCoords(3,i)),&
                                     indexingObj%mapToGlobalIndices(this%indexingData%colVarName,mappedColIndices))

            this%indexingData%rowDataCoordsGlobal(i)%rowCoords = rowCoords(:,i)
            this%indexingData%localRowIndices(i) = indexingObj%findDistIndex(rowCoords(1,i),&
                                                                rowCoords(2,i),&
                                                                rowCoords(3,i),&
                                                                .true.)
            this%indexingData%rowDataCoordsGlobal(i)%colCoords = allCombinationsMappedInd

            this%indexingData%rowDataCoordsLocal(i)%rowCoords = rowCoords(:,i) 
            this%indexingData%rowDataCoordsLocal(i)%rowCoords(1) = indexingObj%findLocalXIndex(rowCoords(1,i))

            do j = 1, size(mappedColIndicesLocal(1)%entry)
                mappedColIndicesLocal(1)%entry(j) = indexingObj%findLocalXIndex(mappedColIndicesLocal(1)%entry(j),&
                                                                                locXInd=rowCoords(1,i))
            end do

            this%indexingData%rowDataCoordsLocal(i)%colCoords = allCombinations(mappedColIndicesLocal)

            allocate(this%multConst(i)%entry(size(this%rowData%values(i)%entry)))
            this%multConst(i)%entry = real(1.0d00,kind=rk)

            allocate(this%indexingData%localColIndices(i)%entry(size(this%rowData%values(i)%entry)))

            if (size(mappedColIndices) == 1) then
                do j = 1, size(this%indexingData%localColIndices(i)%entry)
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findLocalXIndex(allCombinationsMappedInd(1,j),&
                                                                                                locXInd=rowCoords(1,i))
                end do
            else 
                do j = 1, size(this%indexingData%localColIndices(i)%entry)
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findDistIndex(allCombinationsMappedInd(1,j),&
                                                                                allCombinationsMappedInd(2,j),&
                                                                                allCombinationsMappedInd(3,j),&
                                                                                .true.,&
                                                                                locXInd=rowCoords(1,i))
                end do

            end if
        end do

    end if

end subroutine initRowData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addRowDataToPreallocationData(this,petscPreallocData) 
    !! Add this term's row data to a PETSc preallocation object

    class(MatrixTerm)            ,intent(inout)  :: this
    type(PETScPreallocationData) ,intent(inout)  :: petscPreallocData

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to add undefined matrix term row data to PETSc preallocation data")
        call assert(petscPreallocData%isDefined(),"Attempted to add matrix term row data to undefined PETSc preallocation data")
    end if

    call petscPreallocData%addRowDataToPattern(this%rowData)

end subroutine addRowDataToPreallocationData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addRowDataPatternToController(this,petscCont) 
    !! Add this term's row data to a petsc preallocation object in PETSc controller

    class(MatrixTerm)     ,intent(in)     :: this
    type(PETScController) ,intent(inout)  :: petscCont

    call petscCont%addRowDataToPreallocation(this%rowData)

end subroutine addRowDataPatternToController
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addRowValuesToPETScMatrix(this,petscCont,mult,petscGroup) 
    !! Add this term's row values to a petsc matrix object in PETSc controller, multiplied by mult. If evolved variable is
    !! stationary mult is ignored and set to -1

    class(MatrixTerm)     ,intent(in)     :: this
    type(PETScController) ,intent(inout)  :: petscCont
    real(rk)              ,intent(in)     :: mult
    integer(ik) ,optional ,intent(in)     :: petscGroup

    real(rk)                              :: usedMult
    integer(ik)                           :: usedGroup

    if (assertions) call assert(this%isDefined(),"Attempted to add undefined matrix term values to PETSc preallocation data")

    usedMult = mult 
    if (this%stationaryEvolvedVar) usedMult = real(-1.0d0,kind=rk)

    usedGroup = 1
    if (present(petscGroup)) usedGroup = petscGroup
    call petscCont%addRowValuesToMatrix(this%rowData,usedMult,usedGroup) 

end subroutine addRowValuesToPETScMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
module function evaluateMatTerm (this,varCont) result(res)
    !! Get explicit value for the term by crudely dotting each sparse row with the locally stored implicit variable vector in varCont

    class(MatrixTerm)                    ,intent(in) :: this
    type(VariableContainer)              ,intent(in) :: varCont
    real(rk) ,allocatable ,dimension(:)              :: res

    integer(ik) :: i ,varIndexRow ,varIndexCol

    if (assertions) then
        call assertPure(this%isDefined(),"Attempted to calculate explicit version of undefined matrix term object")
        call assertPure(varCont%isDefined(),&
        "Attempted to calculate explicit version of matrix term object by passing undefined variable container")
    end if

    varIndexRow = varCont%getVarIndex(this%indexingData%rowVarName)
    varIndexCol = varCont%getVarIndex(this%indexingData%colVarName) 

    allocate(res,mold=varCont%variables(varIndexRow)%entry)
    res = 0

    do i = 1, size(this%rowData%values)
        res(this%indexingData%localRowIndices(i)) = dot_product(this%rowData%values(i)%entry,&
                                                varCont%variables(varIndexCol)%entry(this%indexingData%localColIndices(i)%entry))
    end do

end function evaluateMatTerm
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setMultConst(this,multConst) 
    !! Setter for multConst

    class(MatrixTerm)              ,intent(inout)  :: this
    type(RealArray) ,dimension(:)  ,intent(in)     :: multConst

    this%multConst = multConst

end subroutine setMultConst
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMultConst (this) result(multConst)
    !! Getter for multConst

    class(MatrixTerm)                          ,intent(in) :: this
    type(RealArray) ,allocatable ,dimension(:)             :: multConst

    if (assertions) call assertPure(this%isDefined(),"Called getMultConst on undefined matrix term")

    multConst = this%multConst

end function getMultConst
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRowVarName(this) result(name)
    !! Get name of the evolved variable of this term

    class(MatrixTerm)                    ,intent(in) :: this
    character(:) ,allocatable                        :: name

    if (assertions) call assertPure(this%isDefined(),"Called getRowVarName on undefined matrix term")

    name = this%indexingData%rowVarName 

end function getRowVarName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setFixedMatrix(this,isFixed) 
    !! Set whether the matrix is fixed -> only calculated once

    class(MatrixTerm)              ,intent(inout)  :: this
    logical                        ,intent(in)     :: isFixed

    this%fixedMatrix = isFixed

end subroutine setFixedMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNonTrivialRowFun(this,nontriv) 
    !! Set whether the matrix has a non-trivial row function

    class(MatrixTerm)              ,intent(inout)  :: this
    logical                        ,intent(in)     :: nontriv

    this%hasNonTrivialRowFun = nontriv

end subroutine setNonTrivialRowFun
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNonTrivialColFun(this,nontriv) 
    !! Set whether the matrix has a non-trivial col function

    class(MatrixTerm)              ,intent(inout)  :: this
    logical                        ,intent(in)     :: nontriv

    this%hasNonTrivialColFun = nontriv

end subroutine setNonTrivialColFun
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getIndexingData(this) result(indData)
    !! Getter for indexingData

    class(MatrixTerm)                    ,intent(in) :: this
    type(MatrixTermIndexingData)                     :: indData

    if (assertions) call assertPure(this%isDefined(),"Called getIndexingData on undefined matrix term")

    indData = this%indexingData 
    
end function getIndexingData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRowData(this) result(rowData)
    !! Getter for rowData

    class(MatrixTerm)           ,intent(in) :: this
    type(SparseRowData)                     :: rowData
    
    if (assertions) call assertPure(this%isDefined(),"Called getRowData on undefined matrix term")

    rowData = this%rowData 
    
end function getRowData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine matTermUpdate(this,varCont,modelData,hostModel) 
    !! Default matrix term update, call matrixTermUpdate 

    class(MatrixTerm)               ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont
    class(ModelboundData) ,optional ,intent(in)     :: modelData
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel

    if (present(modelData)) then 
        if (present(hostModel)) then 
            call this%matrixTermUpdate(varCont,this%indexingData,modelData=modelData,hostModel=hostModel)
        else
            call this%matrixTermUpdate(varCont,this%indexingData,modelData=modelData)
        end if
    else
        if (present(hostModel)) then 
            call this%matrixTermUpdate(varCont,this%indexingData,hostModel=hostModel)
        else
            call this%matrixTermUpdate(varCont,this%indexingData)
        end if
    end if

end subroutine matTermUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine matrixNoUpdate(this,varCont,indexingData,modelData,hostModel) 
    !! Default MatrixTerm updateMatTerm function - does nothing

    class(MatrixTerm)               ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont
    type(MatrixTermIndexingData)    ,intent(in)     :: indexingData
    class(ModelboundData) ,optional ,intent(in)     :: modelData
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel

end subroutine matrixNoUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule matrix_term_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
