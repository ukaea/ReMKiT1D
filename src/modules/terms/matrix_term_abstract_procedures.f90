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

    integer(ik) :: i ,k ,nRows
    real(rk)    :: norm ,rVal

    if (assertions) then
        call assertPure(this%isDefined(),"Attempted to calculate values of undefined matrix term object")
        call assertPure(varCont%isDefined(),&
        "Attempted to calculate values of matrix term object by passing undefined variable container")
    end if

    if (.not. this%fixedMatrixCalculated) then
        nRows = size(this%flatRowPtr) - 1
        norm = this%normalizationConst

        associate(ptr => this%flatRowPtr, fVals => this%flatValues, fMult => this%flatMultConst)
            if (this%hasNonTrivialRowFun) then
                if (this%hasNonTrivialColFun) then 

                    call this%calculateRows(varCont,this%rowBuffer,this%indexingData)
                    call this%calculateCols(varCont,this%colBuffer,this%indexingData)
                    do i = 1, nRows
                        rVal = norm * this%rowBuffer(i)
                        fVals(ptr(i):ptr(i+1)-1) = rVal * fMult(ptr(i):ptr(i+1)-1) * this%colBuffer(i)%entry
                    end do

                else

                    call this%calculateRows(varCont,this%rowBuffer,this%indexingData)
                    if (this%isOnePerRow) then
                        ! When every row has exactly one nonzero, ptr(i) == i so values map 1-to-1 with rowBuffer.
                        do i = 1, nRows
                            fVals(i) = norm * this%rowBuffer(i) * fMult(i)
                        end do
                    else
                        do i = 1, nRows
                            rVal = norm * this%rowBuffer(i)
                            do k = ptr(i), ptr(i+1) - 1
                                fVals(k) = rVal * fMult(k)
                            end do
                        end do
                    end if

                end if

            else

                if (this%hasNonTrivialColFun) then

                    call this%calculateCols(varCont,this%colBuffer,this%indexingData)
                    do i = 1, nRows
                        fVals(ptr(i):ptr(i+1)-1) = norm * fMult(ptr(i):ptr(i+1)-1) * this%colBuffer(i)%entry
                    end do

                else

                    fVals = norm * fMult

                end if

            end if
        end associate
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

    integer(ik) :: i ,j ,numRows ,numCols ,nnz

    type(IntArray) ,allocatable ,dimension(:) :: mappedColIndices ,mappedColIndicesLocal ,colVectors
    integer(ik) ,allocatable ,dimension(:)    :: rowIndices
    integer(ik) ,allocatable ,dimension(:,:)  :: allCombinationsMappedInd 

    if (assertions) then 
        call assert(indexingObj%isDefined(),"initRowData for matrix term called with undefined indexing object")
        call assert((size(rowCoords,1) == 1) .or. (size(rowCoords,1) == 3),"Row coordinate entries passed to initRowData must&
        & be of dimension (1,:) or (3,:) depending on whether the row variable is a distribution")
    end if

    numRows = size(rowCoords,2)
    allocate(rowIndices(numRows))
    allocate(colVectors(numRows))
    allocate(this%indexingData%rowDataCoordsGlobal(numRows))
    allocate(this%indexingData%rowDataCoordsLocal(numRows))
    allocate(this%colBuffer(numRows))
    allocate(this%rowBuffer(numRows))
    allocate(this%indexingData%localRowIndices(numRows))
    allocate(this%indexingData%localColIndices(numRows))

    this%rowBuffer = real(1,kind=rk)
    if (size(rowCoords,1) == 1) then
        do i = 1,numRows
            mappedColIndices = rowToColMapping(rowCoords(:,i))
            mappedColIndicesLocal = mappedColIndices
            allCombinationsMappedInd = allCombinations(mappedColIndices)
            rowIndices(i) = indexingObj%findIndex(this%indexingData%rowVarName,rowCoords(1,i))
            colVectors(i)%entry = indexingObj%mapToGlobalIndices(this%indexingData%colVarName,mappedColIndices)
            numCols = size(colVectors(i)%entry)

            if (assertions) then
                call assert(.not. any(rowIndices(1:i-1) == rowIndices(i)),&
                "Attempted to add row vector with same index as one already in sparse row data object")
            end if

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

            allocate(this%colBuffer(i)%entry(numCols))
            this%colBuffer(i)%entry = real(1.0d00,kind=rk)

            allocate(this%indexingData%localColIndices(i)%entry(numCols))

            if (size(mappedColIndices) == 1) then
                do j = 1, numCols
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findLocalXIndex(allCombinationsMappedInd(1,j),&
                                                                                                locXInd=rowCoords(1,i))
                end do
            else 
                do j = 1, numCols
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findDistIndex(allCombinationsMappedInd(1,j),&
                                                                                 allCombinationsMappedInd(2,j),&
                                                                                 allCombinationsMappedInd(3,j),&
                                                                                 .true.,&
                                                                                 locXInd=rowCoords(1,i))
                end do

            end if

        end do
    else 

        do i = 1,numRows
            mappedColIndices = rowToColMapping(rowCoords(:,i))
            mappedColIndicesLocal = mappedColIndices
            allCombinationsMappedInd = allCombinations(mappedColIndices)
            rowIndices(i) = indexingObj%findIndex(this%indexingData%rowVarName,rowCoords(1,i),&
                                                                   rowCoords(2,i),&
                                                                   rowCoords(3,i))
            colVectors(i)%entry = indexingObj%mapToGlobalIndices(this%indexingData%colVarName,mappedColIndices)
            numCols = size(colVectors(i)%entry)

            if (assertions) then
                call assert(.not. any(rowIndices(1:i-1) == rowIndices(i)),&
                "Attempted to add row vector with same index as one already in sparse row data object")
            end if

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

            allocate(this%colBuffer(i)%entry(numCols))
            this%colBuffer(i)%entry = real(1.0d00,kind=rk)

            allocate(this%indexingData%localColIndices(i)%entry(numCols))

            if (size(mappedColIndices) == 1) then
                do j = 1, numCols
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findLocalXIndex(allCombinationsMappedInd(1,j),&
                                                                                                locXInd=rowCoords(1,i))
                end do
            else 
                do j = 1, numCols
                    this%indexingData%localColIndices(i)%entry(j) = indexingObj%findDistIndex(allCombinationsMappedInd(1,j),&
                                                                                allCombinationsMappedInd(2,j),&
                                                                                allCombinationsMappedInd(3,j),&
                                                                                .true.,&
                                                                                locXInd=rowCoords(1,i))
                end do

            end if
        end do

    end if

    call this%rowData%init(rowIndices,colVectors)

    ! Build private contiguous CSR storage (authoritative for values and multConst once initRowData has run).
    if (allocated(this%multConst)) deallocate(this%multConst)
    if (allocated(this%flatRowPtr)) deallocate(this%flatRowPtr)
    if (allocated(this%flatValues)) deallocate(this%flatValues)
    if (allocated(this%flatBufferIndices)) deallocate(this%flatBufferIndices)
    if (allocated(this%flatMultConst)) deallocate(this%flatMultConst)

    allocate(this%flatRowPtr(numRows+1))
    this%flatRowPtr(1) = 1
    do i = 1, numRows
        this%flatRowPtr(i+1) = this%flatRowPtr(i) + size(colVectors(i)%entry)
    end do
    nnz = this%flatRowPtr(numRows+1) - 1

    allocate(this%flatValues(nnz))
    this%flatValues = real(0.0d0,kind=rk)
    allocate(this%flatMultConst(nnz))
    this%flatMultConst = real(1.0d0,kind=rk)
    this%isOnePerRow = (numRows > 0 .and. all(this%flatRowPtr(2:numRows+1) - this%flatRowPtr(1:numRows) == 1))

    ! Pre-flatten local (x, v) coordinates for 3D distribution terms so pure accumulateRowMoment avoids jagged lookups.
    this%momentDim = 0
    if (allocated(this%momentRowX)) deallocate(this%momentRowX)
    if (allocated(this%momentRowV)) deallocate(this%momentRowV)
    if (allocated(this%momentColV)) deallocate(this%momentColV)
    if (size(rowCoords,1) == 3 .and. numRows > 0) then
        allocate(this%momentRowX(numRows), this%momentRowV(numRows))
        this%momentDim = size(this%indexingData%rowDataCoordsLocal(1)%colCoords,1)
        if (this%momentDim == 3) allocate(this%momentColV(nnz))
        do i = 1, numRows
            this%momentRowX(i) = this%indexingData%rowDataCoordsLocal(i)%rowCoords(1)
            this%momentRowV(i) = this%indexingData%rowDataCoordsLocal(i)%rowCoords(3)
            if (this%momentDim == 3) then
                this%momentColV(this%flatRowPtr(i):this%flatRowPtr(i+1)-1) = &
                    this%indexingData%rowDataCoordsLocal(i)%colCoords(3,:)
            end if
        end do
    end if

end subroutine initRowData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addRowDataToPreallocationData(this,petscPreallocData) 
    !! Add this term's row data to a PETSc preallocation object and cache flattened value-buffer indices

    class(MatrixTerm)            ,intent(inout)  :: this
    type(PETScPreallocationData) ,intent(inout)  :: petscPreallocData

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to add undefined matrix term row data to PETSc preallocation data")
        call assert(petscPreallocData%isDefined(),"Attempted to add matrix term row data to undefined PETSc preallocation data")
    end if

    call petscPreallocData%addRowDataToPattern(this%rowData)
    call syncFlatBufferIndices(this)

end subroutine addRowDataToPreallocationData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addRowDataPatternToController(this,petscCont) 
    !! Add this term's row data to a petsc preallocation object in PETSc controller and cache flattened value-buffer indices

    class(MatrixTerm)     ,intent(inout)  :: this
    type(PETScController) ,intent(inout)  :: petscCont

    call petscCont%addRowDataToPreallocation(this%rowData)
    call syncFlatBufferIndices(this)

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

    if (assertions) then
        call assert(this%isDefined(),"Attempted to add undefined matrix term values to PETSc matrix")
        call assert(allocated(this%flatBufferIndices),&
        "Attempted to add matrix term values to PETSc matrix before preallocation registration")
        call assert(allocated(this%flatValues),&
        "Attempted to add matrix term values to PETSc matrix before initRowData")
    end if

    usedMult = mult 
    if (this%stationaryEvolvedVar) usedMult = real(-1.0d0,kind=rk)

    usedGroup = 1
    if (present(petscGroup)) usedGroup = petscGroup
    call petscCont%addFlatValuesToMatrix(this%flatBufferIndices,this%flatValues,usedMult,usedGroup)

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

    do i = 1, size(this%flatRowPtr) - 1
        res(this%indexingData%localRowIndices(i)) = dot_product( &
            this%flatValues(this%flatRowPtr(i):this%flatRowPtr(i+1)-1), &
            varCont%variables(varIndexCol)%entry(this%indexingData%localColIndices(i)%entry))
    end do

end function evaluateMatTerm
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setMultConst(this,multConst) 
    !! Setter for multConst. After initRowData, writes directly into contiguous flatMultConst; before initRowData,
    !! stores jagged multConst directly.

    class(MatrixTerm)              ,intent(inout)  :: this
    type(RealArray) ,dimension(:)  ,intent(in)     :: multConst

    integer(ik) :: i ,nRows ,p0 ,p1

    if (allocated(this%flatRowPtr)) then
        nRows = size(this%flatRowPtr) - 1
        if (assertions) call assertPure(size(multConst) == nRows,&
            "Row count passed to MatrixTerm%setMultConst does not match initialized row data")
        do i = 1, nRows
            p0 = this%flatRowPtr(i)
            p1 = this%flatRowPtr(i+1) - 1
            if (assertions) call assertPure(allocated(multConst(i)%entry),&
                "Unallocated row entry passed to MatrixTerm%setMultConst")
            ! Empty sparse rows (p1 < p0, e.g. boundary-clipped harmonic stencils) have no matrix entries to populate;
            ! stencil generators may still emit unclipped template weights for them, so only non-empty rows must match length.
            if (p1 >= p0) then
                if (assertions) call assertPure(size(multConst(i)%entry) == p1 - p0 + 1,&
                    "Row entry length passed to MatrixTerm%setMultConst does not match initialized row data")
                this%flatMultConst(p0:p1) = multConst(i)%entry
            end if
        end do
    else
        call copyRealArraysInPlace(this%multConst,multConst)
    end if

end subroutine setMultConst
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMultConst (this) result(multConst)
    !! Getter for multConst

    class(MatrixTerm)                          ,intent(in) :: this
    type(RealArray) ,allocatable ,dimension(:)             :: multConst

    integer(ik) :: i ,nRows

    if (assertions) call assertPure(this%isDefined(),"Called getMultConst on undefined matrix term")

    if (allocated(this%flatRowPtr)) then
        nRows = size(this%flatRowPtr) - 1
        allocate(multConst(nRows))
        do i = 1, nRows
            multConst(i)%entry = this%flatMultConst(this%flatRowPtr(i):this%flatRowPtr(i+1)-1)
        end do
    else
        multConst = this%multConst
    end if

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
    !! Getter for rowData (syncs authoritative flatValues into the returned SparseRowData copy)

    class(MatrixTerm)           ,intent(in) :: this
    type(SparseRowData)                     :: rowData

    integer(ik) :: i
    
    if (assertions) call assertPure(this%isDefined(),"Called getRowData on undefined matrix term")

    rowData = this%rowData
    if (allocated(this%flatRowPtr)) then
        do i = 1, size(rowData%rowIndex)
            rowData%values(i)%entry = this%flatValues(this%flatRowPtr(i):this%flatRowPtr(i+1)-1)
        end do
    end if
    
end function getRowData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine accumulateRowMoment(this,vVec,res)
    !! Accumulate velocity moment of this term's flatValues into res without copying SparseRowData.
    !! Reuses pre-allocated res(:)%entry in place on the hot path, and allocates missing entries on demand.

    class(MatrixTerm)              ,intent(in)    :: this
    real(rk) ,dimension(:)         ,intent(in)    :: vVec
    type(RealArray) ,dimension(:)  ,intent(inout) :: res

    integer(ik) :: i ,k ,nRows ,rowX ,colIdx
    real(rk)    :: vWeight
    logical     :: needsAlloc

    if (assertions) call assertPure(this%isDefined(),"Called accumulateRowMoment on undefined matrix term")

    nRows = size(this%flatRowPtr) - 1
    if (nRows == 0) then
        do i = 1, size(res)
            if (.not. allocated(res(i)%entry)) allocate(res(i)%entry(0))
            res(i)%entry(:) = 0.0d0
        end do
        return
    end if

    associate(ptr => this%flatRowPtr, fVals => this%flatValues, &
              rX => this%momentRowX, rV => this%momentRowV)
        if (this%momentDim == 3) then
            do i = 1, size(res)
                if (.not. allocated(res(i)%entry)) allocate(res(i)%entry(size(vVec)))
                res(i)%entry(:) = 0.0d0
            end do
            associate(cV => this%momentColV)
                do i = 1, nRows
                    rowX = rX(i)
                    vWeight = vVec(rV(i))
                    do k = ptr(i), ptr(i+1) - 1
                        colIdx = cV(k)
                        res(rowX)%entry(colIdx) = res(rowX)%entry(colIdx) + vWeight * fVals(k)
                    end do
                end do
            end associate
        else
            needsAlloc = .false.
            do i = 1, size(res)
                if (.not. allocated(res(i)%entry)) then
                    needsAlloc = .true.
                    exit
                end if
            end do
            if (needsAlloc) then
                do i = 1, nRows
                    rowX = rX(i)
                    if (.not. allocated(res(rowX)%entry)) allocate(res(rowX)%entry(ptr(i+1) - ptr(i)))
                end do
                do i = 1, size(res)
                    if (.not. allocated(res(i)%entry)) allocate(res(i)%entry(0))
                end do
            end if
            do i = 1, size(res)
                res(i)%entry(:) = 0.0d0
            end do
            do i = 1, nRows
                rowX = rX(i)
                vWeight = vVec(rV(i))
                res(rowX)%entry(:) = res(rowX)%entry(:) + vWeight * fVals(ptr(i):ptr(i+1)-1)
            end do
        end if
    end associate

end subroutine accumulateRowMoment
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine multiplyRowValuesInPlace(this,lhs,res)
    !! Compute element-wise product lhs * this%flatValues in place without copying SparseRowData

    class(MatrixTerm)                          ,intent(in)    :: this
    type(RealArray) ,dimension(:)              ,intent(in)    :: lhs
    type(RealArray) ,allocatable ,dimension(:) ,intent(inout) :: res

    integer(ik) :: i ,nRows ,nCols ,p0 ,p1

    if (assertions) then
        call assertPure(this%isDefined(),"Called multiplyRowValuesInPlace on undefined matrix term")
        call assertPure(size(lhs) == size(this%flatRowPtr) - 1,&
            "Row sizes in element-wise sparse row data/real array multiplication do not conform")
    end if

    nRows = size(this%flatRowPtr) - 1
    if (.not. allocated(res)) then
        allocate(res(nRows))
    else if (size(res) /= nRows) then
        deallocate(res)
        allocate(res(nRows))
    end if

    do i = 1, nRows
        p0 = this%flatRowPtr(i)
        p1 = this%flatRowPtr(i+1) - 1
        nCols = p1 - p0 + 1
        if (assertions) call assertPure(size(lhs(i)%entry) == nCols,&
            "Column vectors do not conform in element-wise sparse row data/real array multiplication")
        if (.not. allocated(res(i)%entry)) then
            allocate(res(i)%entry(nCols))
        else if (size(res(i)%entry) /= nCols) then
            deallocate(res(i)%entry)
            allocate(res(i)%entry(nCols))
        end if
        res(i)%entry(:) = lhs(i)%entry(:) * this%flatValues(p0:p1)
    end do

end subroutine multiplyRowValuesInPlace
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
pure subroutine syncFlatBufferIndices(this)
    !! Flatten rowData%cooIndices into private flatBufferIndices aligned 1-to-1 with flatRowPtr/flatValues

    class(MatrixTerm) ,intent(inout) :: this

    integer(ik) :: i ,nRows ,nnz ,p0 ,p1

    if (.not. allocated(this%flatRowPtr) .or. .not. allocated(this%rowData%cooIndices)) return

    nRows = size(this%flatRowPtr) - 1
    nnz = this%flatRowPtr(nRows+1) - 1

    if (allocated(this%flatBufferIndices)) deallocate(this%flatBufferIndices)
    allocate(this%flatBufferIndices(nnz))

    do i = 1, nRows
        p0 = this%flatRowPtr(i)
        p1 = this%flatRowPtr(i+1) - 1
        if (p1 >= p0) this%flatBufferIndices(p0:p1) = this%rowData%cooIndices(i)%entry
    end do

end subroutine syncFlatBufferIndices
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule matrix_term_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
