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
submodule (sparse_row_data_class) sparse_row_data_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the sparse row data class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initSparseRowData(this,rowIndices,colVectors) 
    !! Sparse row data Object initialization routine. If indices are supplied initializes values to 0.

    class(SparseRowData)                         ,intent(inout)  :: this
    integer(ik)          ,dimension(:) ,optional ,intent(in)     :: rowIndices !! Row indices of each sparse row vector
    type(IntArray)       ,dimension(:) ,optional ,intent(in)     :: colVectors !! Arrays of column indices for each row

    integer(ik) :: i

    if (assertions) then 
        call assertPure(present(rowIndices) .eqv. present(colVectors),"If sparse row data constructor is called with initial &
        &row/column data both must be included")

        if (present(rowIndices)) call assertPure(size(rowIndices) == size(colVectors),"Size of row vector passed to sparse row &
        &data constructor must be the same as the size of the colVectors array")
    end if

    call this%makeDefined()

    if (present(rowIndices)) then 
        this%rowIndex = rowIndices
        this%columnVector = colVectors
        allocate(this%values(size(rowIndices)))
        do i = 1,size(rowIndices)
            allocate(this%values(i)%entry(size(colVectors(i)%entry)))
            this%values(i)%entry = 0
        end do
    else 
        allocate(this%rowIndex(0))
        allocate(this%columnVector(0))
        allocate(this%values(0))
    end if
end subroutine initSparseRowData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addRow(this,rowIndex,columnIndices,values) 
    !!  Add sparse row vector to data Object and initialize its values to 0, unless value vector is supplied

    class(SparseRowData)               ,intent(inout)  :: this
    integer(ik)                        ,intent(in)     :: rowIndex !! Row index of added vector
    integer(ik)          ,dimension(:) ,intent(in)     :: columnIndices !! Column indices of nonzeros in row
    real(rk) ,optional   ,dimension(:) ,intent(in)     :: values !! Values of nonzeros in row

    integer(ik) :: i

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add row vector to unitialized sparse row data object")

        call assertPure(.not. any(this%rowIndex == rowIndex),"Attempted to add row vector with same index as one already in &
        &sparse row data object")

        if (present(values)) call assertPure(size(values) == size(columnIndices),"Values passed to addRow routine of sparse row&
        & data object do not conform to passed column indices")
    end if

    this%rowIndex = [this%rowIndex,rowIndex]
    this%columnVector =  [this%columnVector,intArray(columnIndices)]
    if (present(values)) then 
        this%values = [this%values,realArray(values)]
    else 
        this%values = [this%values,realArray([(real(0.0,kind=rk),i=1,size(columnIndices))])]
    end if

end subroutine addRow   
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function multElementwise(mat1,mat2) result(res)
    !! Sparse row data Object initialization routine. If indices are supplied initializes values to 0.

    class(SparseRowData)  ,intent(in)  :: mat1
    class(SparseRowData)  ,intent(in)  :: mat2
    type(SparseRowData)                :: res

    integer(ik) :: i
    
    if (assertions) then 
        call assertPure(size(mat1%rowIndex)==size(mat2%rowIndex),&
        "Row sizes in element-wise sparse row data multiplication do not conform")

        do i = 1,size(mat1%rowIndex)
            call assertPure(size(mat1%columnVector(i)%entry) == size(mat2%columnVector(i)%entry),&
            "Column vectors do not conform in element-wise sparse row data multiplication")
        end do
    end if

    res = mat1

    do i = 1,size(res%rowIndex)
        res%values(i)%entry = res%values(i)%entry * mat2%values(i)%entry
    end do
end function multElementwise  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function multElementwiseRArray(mat1,mat2) result(res)

    type(RealArray) ,dimension(:)  ,intent(in)  :: mat1
    class(SparseRowData)           ,intent(in)  :: mat2
    type(RealArray) ,allocatable ,dimension(:)  :: res

    integer(ik) :: i
    
    if (assertions) then 
        call assertPure(size(mat1)==size(mat2%rowIndex),&
        "Row sizes in element-wise sparse row data/real array multiplication do not conform")

        do i = 1,size(mat1)
            call assertPure(size(mat1(i)%entry) == size(mat2%columnVector(i)%entry),&
            "Column vectors do not conform in element-wise sparse row data/real array multiplication")
        end do
    end if

    res = mat1

    do i = 1,size(res)
        res(i)%entry = res(i)%entry * mat2%values(i)%entry
    end do

end function multElementwiseRArray      
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule sparse_row_data_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
