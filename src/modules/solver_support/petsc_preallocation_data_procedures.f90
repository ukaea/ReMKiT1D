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
submodule (petsc_preallocation_data_class) petsc_preallocation_procedures
    !! author: Stefan Mijin
    !! 
    !! Contains module procedures associated with the petsc preallocation data class


implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initPetscPreallocationData(this,indexingObj,procRank) 
    !! Petsc preallocation data object initialization routine

    class(PETScPreallocationData) ,intent(inout)  :: this
    type(Indexing)                ,intent(in)     :: indexingObj !! Indexing object to retrieve DoF data from
    integer(ik)                   ,intent(in)     :: procRank !! Current process rank

    integer(ik) ,allocatable ,dimension(:)        :: procDoFs

    if (assertions) call assertPure(indexingObj%isDefined(),"Undefined indexing object passed to petsc preallocation data &
    &constructor")

    call this%makeDefined()

    procDoFs = indexingObj%getProcDoF()

    this%localNumDoFs = procDoFs(procRank+1)
    this%localDoFOffset = sum(procDoFs(1:procRank))
    this%totalNumDoFs = sum(procDoFs)

    allocate(this%totalCOOPattern)
    call this%totalCOOPattern%init(this%localNumDoFs,sum(procDoFs),rowOffset=this%localDoFOffset)

    allocate(this%rowNumNonzeros(this%localNumDoFs))
    allocate(this%rowNumNonzerosDiag(this%localNumDoFs))
    allocate(this%rowNumNonzerosOffDiag(this%localNumDoFs))

    this%rowNumNonzeros = 0
    this%rowNumNonzerosDiag = 0
    this%rowNumNonzerosOffDiag = 0

    this%assembled = .false.

end subroutine initPetscPreallocationData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addRowDataToPattern(this,rowData) 
    !! Add sparse row data structure to total sparsity pattern 

    class(PETScPreallocationData)      ,intent(inout)  :: this
    type(SparseRowData)                ,intent(inout)  :: rowData

    integer(ik) :: i ,j

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add row data structure to undefined petsc preallocation data object")
        call assertPure(rowData%isDefined(),"Attempted to add undefined row data structure to petsc preallocation data object")
    end if

    if (allocated(rowData%cooIndices)) deallocate(rowData%cooIndices)
    allocate(rowData%cooIndices(size(rowData%rowIndex)))

    do i = 1,size(rowData%rowIndex)
        allocate(rowData%cooIndices(i)%entry(size(rowData%columnVector(i)%entry)))
        do j = 1,size(rowData%columnVector(i)%entry)
            call this%totalCOOPattern%addEntry(rowData%rowIndex(i),rowData%columnVector(i)%entry(j),&
                                               rowData%cooIndices(i)%entry(j))
        end do 
    end do

end subroutine addRowDataToPattern
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine assembleData(this) 
    !! Assemble preallocation data

    class(PETScPreallocationData)      ,intent(inout)  :: this

    integer(ik) :: i ,firstLocRow, lastLocRow ,currentRow, currentCol

    if (assertions) call assertPure(this%isDefined(),"Attempted to assemble undefined petsc preallocation data")

    firstLocRow = this%localDoFOffset + 1
    lastLocRow = this%localDoFOffset + this%localNumDoFs

    do i = 1, this%totalCOOPattern%numNonZeros
        currentRow = this%totalCOOPattern%rowIndex(i)
        currentCol = this%totalCOOPattern%colIndex(i)
        if ((currentCol >= firstLocRow) .and. (currentCol <= lastLocRow)) then 
            this%rowNumNonzerosDiag(currentRow-this%localDoFOffset) = this%rowNumNonzerosDiag(currentRow-this%localDoFOffset) + 1
        else 
            this%rowNumNonzerosOffDiag(currentRow-this%localDoFOffset) &
            = this%rowNumNonzerosOffDiag(currentRow-this%localDoFOffset) + 1
        end if
    end do

    this%rowNumNonzeros = this%rowNumNonzerosDiag + this%rowNumNonzerosOffDiag

    this%assembled = .true.

end subroutine assembleData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isAssembled (this) result(assembled)
    !! Return assembly status of preallocation data

    class(PETScPreallocationData) ,intent(in) :: this
    logical                                   :: assembled

    if (assertions) call assertPure(this%isDefined(),"Attempted to get assembly status from undefined petsc preallocation data&
    & object")

    assembled = this%assembled

end function isAssembled
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumNonzerosDiag (this) result(numNZDiag)
    !! Return number of nonzeros in the diagonal part of each row

    class(PETScPreallocationData)          ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: numNZDiag

    if (assertions) call assertPure(this%isDefined(),"Attempted to get number of diagonal block nonzeros from undefined petsc&
    & preallocation data object")

    numNZDiag = this%rowNumNonzerosDiag

end function getNumNonzerosDiag
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumNonzerosOffDiag (this) result(numNZODiag)
    !! Return number of nonzeros in the off-diagonal part of each row

    class(PETScPreallocationData)          ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: numNZODiag

    if (assertions) call assertPure(this%isDefined(),"Attempted to get number of off-diagonal block nonzeros from undefined petsc&
    & preallocation data object")

    numNZODiag = this%rowNumNonzerosOffDiag

end function getNumNonzerosOffDiag
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLocalNumDoFs (this) result(locNumDoFs)
    !! Return local number of DoFs in unknown vector

    class(PETScPreallocationData) ,intent(in) :: this
    integer(ik)                               :: locNumDoFs

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to get number of local DoFs from undefined petsc preallocation object")

    locNumDoFs = this%localNumDoFs

end function getLocalNumDoFs
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTotalNumDoFs (this) result(totalNumDoFs)
    !! Return total number of DoFs in unknown vector

    class(PETScPreallocationData) ,intent(in) :: this
    integer(ik)                               :: totalNumDoFs

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to get number of total DoFs from undefined petsc preallocation object")

    totalNumDoFs = this%totalNumDoFs

end function getTotalNumDoFs
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLocalDoFOffset (this) result(offset)
    !! Return offset of the start of this processors local vector in the global vector

    class(PETScPreallocationData) ,intent(in) :: this
    integer(ik)                               :: offset

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to get local DoF offset from undefined petsc preallocation object")

    offset = this%localDoFOffset

end function getLocalDoFOffset
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine deallocatePattern(this)
    !!  Deallocates the COO matrix pattern

    class(PETScPreallocationData)      ,intent(inout)  :: this

    if (allocated(this%totalCOOPattern)) deallocate(this%totalCOOPattern)

end subroutine deallocatePattern
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumNonzerosCOO (this) result(numNZ)
    !! Return total number of deduplicated COO entries on this processor

    class(PETScPreallocationData) ,intent(in) :: this
    integer(ik)                               :: numNZ

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to get number of COO nonzeros from undefined petsc preallocation object")

    numNZ = this%totalCOOPattern%numNonZeros

end function getNumNonzerosCOO
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCOORowIndices0 (this) result(rows0)
    !! Return 0-based global row indices of deduplicated COO pattern

    class(PETScPreallocationData)          ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: rows0

    integer(ik) :: nnz

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to get COO row indices from undefined petsc preallocation object")

    nnz = this%totalCOOPattern%numNonZeros
    allocate(rows0(max(1,nnz)))
    rows0 = 0
    if (nnz > 0) rows0(1:nnz) = this%totalCOOPattern%rowIndex(1:nnz) - 1

end function getCOORowIndices0
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCOOColIndices0 (this) result(cols0)
    !! Return 0-based global column indices of deduplicated COO pattern

    class(PETScPreallocationData)          ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: cols0

    integer(ik) :: nnz

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to get COO col indices from undefined petsc preallocation object")

    nnz = this%totalCOOPattern%numNonZeros
    allocate(cols0(max(1,nnz)))
    cols0 = 0
    if (nnz > 0) cols0(1:nnz) = this%totalCOOPattern%colIndex(1:nnz) - 1

end function getCOOColIndices0
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule petsc_preallocation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
