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
submodule (coo_sparsity_pattern_class) coo_sparsity_pattern_procedures
    !! author: Stefan Mijin 
    !! 
    !!  Contains module procedures associated with the basic coo sparsity pattern class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initCooSparsityPattern(this,numRows,numCols,bufferSize,rowOffset) 
    !! COO sparsity pattern initialization routine

    class(COOSparsityPattern)           ,intent(inout)  :: this
    integer(ik)                         ,intent(in)     :: numRows !! Matrix number of rows
    integer(ik)                         ,intent(in)     :: numCols !! Matrix number of columns
    integer(ik) ,optional               ,intent(in)     :: bufferSize !! Preallocated rolling buffer size
    integer(ik) ,optional               ,intent(in)     :: rowOffset !! Row offset for distributed arrays
    
    integer(ik) :: i

    call this%makeDefined()

    this%numNonZeros = 0
    this%numRows = numRows
    this%numCols = numCols
    this%bufferSize = 10000000
    if (present(bufferSize)) this%bufferSize = bufferSize
    allocate(this%rowIndex(this%bufferSize))
    allocate(this%colIndex(this%bufferSize))

    allocate(this%addedRow(this%numRows))
    this%addedRow = .false.
    allocate(this%addedCol(this%numCols))
    this%addedCol = .false.

    allocate(this%firstInstRow(this%numRows))
    this%firstInstRow = 0
    allocate(this%firstInstCol(this%numCols))
    this%firstInstCol = 0

    this%rowOffset = 0 
    if (present(rowOffset)) this%rowOffset = rowOffset

    allocatE(this%rowGlobalIndices(this%numRows))
    do i = 1,this%numRows
        allocate(this%rowGlobalIndices(i)%entry(0))
    end do
end subroutine initCooSparsityPattern
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function hasIndices (this,row,column) result(found)
    !! Return true if all elements with indices (row,column) are present in sparsity pattern

    class(COOSparsityPattern)   ,intent(in) :: this
    integer(ik)   ,dimension(:) ,intent(in) :: row !! Row indices
    integer(ik)   ,dimension(:) ,intent(in) :: column !! Column indices
    logical                                 :: found

    logical      ,allocatable ,dimension(:) :: foundVec

    integer(ik) :: i ,j


    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to check if indices present in undefined sparsity pattern")
        call assertPure(size(row) == size(column), "Row and column index vectors passed to hasIndices must be same size")
    end if

    allocate(foundVec(size(row)))
    foundVec = .false.

    do i = 1,size(row)

        if (this%addedRow(row(i)-this%rowOffset) .and. this%addedCol(column(i))) then 

            do j = 1,size(this%rowGlobalIndices(row(i)-this%rowOffset)%entry)
                if (this%colIndex(this%rowGlobalIndices(row(i)-this%rowOffset)%entry(j)) == column(i)) then
                    foundVec(i) = .true.
                    exit
                end if
            end do

        end if
        
    end do 

    found = all(foundVec)

end function hasIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function findLocationOfIndices (this,row,column) result(indices)
    !! Return indices with given values of (row,column)

    class(COOSparsityPattern)                ,intent(in) :: this
    integer(ik)                ,dimension(:) ,intent(in) :: row !! Row index
    integer(ik)                ,dimension(:) ,intent(in) :: column !! Column index
    integer(ik)   ,allocatable ,dimension(:)             :: indices

    integer(ik) ,allocatable ,dimension(:) :: indicesLoc
    integer(ik) :: i


    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to get location of indices from undefined sparsity pattern")
        call assertPure(this%hasIndices(row,column),"Attempted to get location of indices from sparsity pattern when those&
                                                    & indices are not present")
    end if

    allocate(indices(size(row)))
    indices = 0 
    do i = 1,size(row)

        indicesLoc = findIndices(this%colIndex(1:this%numNonZeros) == column(i) &
                                 .and. this%rowIndex(1:this%numNonZeros) == row(i))

        indices(i) = indicesLoc(1)
    end do

end function findLocationOfIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addEntry(this,row,column) 
    !! Add entry if not already present

    class(COOSparsityPattern)   ,intent(inout)  :: this
    integer(ik)                 ,intent(in)     :: row !! Row index
    integer(ik)                 ,intent(in)     :: column !! Column Index

    integer(ik) ,allocatable ,dimension(:) :: tempRow ,tempCol

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add entry to undefined sparsity pattern")
    end if

    if (.not. this%hasIndices([row],[column])) then 
        if (this%numNonZeros < size(this%rowIndex)) then
            this%numNonZeros = this%numNonZeros + 1
            this%rowIndex(this%numNonZeros) = row
            this%colIndex(this%numNonZeros) = column
        else
            tempCol = this%colIndex
            tempRow = this%rowIndex 

            deallocate(this%rowIndex)
            deallocate(this%colIndex)

            allocate(this%rowIndex(this%bufferSize+size(tempRow)))
            allocate(this%colIndex(this%bufferSize+size(tempRow)))

            this%rowIndex(1:size(tempRow)) = tempRow
            this%colIndex(1:size(tempCol)) = tempCol

            this%numNonZeros = this%numNonZeros + 1
            this%rowIndex(this%numNonZeros) = row
            this%colIndex(this%numNonZeros) = column

        end if
        
        if (.not. this%addedRow(row-this%rowOffset)) this%firstInstRow(row-this%rowOffset) = this%numNonZeros
        if (.not. this%addedCol(column)) this%firstInstCol(column) = this%numNonZeros

        this%addedRow(row-this%rowOffset) = .true.
        this%addedCol(column) = .true.

        this%rowGlobalIndices(row-this%rowOffset)%entry = [this%rowGlobalIndices(row-this%rowOffset)%entry,this%numNonZeros]
    end if

end subroutine addEntry
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addPattern(this,pattern,rowOffset,colOffset) 
    !! Add entries from existing pattern that are not already present

    class(COOSparsityPattern)      ,intent(inout)  :: this
    type(COOSparsityPattern)       ,intent(in)     :: pattern !! Pattern to add to this pattern
    integer(ik)      ,dimension(:) ,intent(in)     :: rowOffset !! Array by which to offset pattern row indices before adding
    integer(ik)      ,dimension(:) ,intent(in)     :: colOffset !! Array by which to offset pattern column indices before adding

    integer(ik) :: i

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add pattern to undefined sparsity pattern")
        call assertPure(pattern%isDefined(),"Attempted to add undefined pattern to sparsity pattern")
        call assertPure(size(colOffset) == size(rowOffset),"colOffset and rowOffset passed to addPattern must be of same size")
        call assertPure(size(colOffset) == pattern%numNonZeros,"Colum/row offsets passed to addPattern must have the same&
                        & number of elements as the added pattern")

    end if

    do i = 1,pattern%numNonZeros
        call this%addEntry(pattern%rowIndex(i)+rowOffset(i),pattern%colIndex(i)+colOffset(i))
    end do 

end subroutine addPattern
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule coo_sparsity_pattern_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
