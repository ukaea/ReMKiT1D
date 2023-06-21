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
module sparse_row_data_class
    !! author: Stefan Mijin 
    !! 
    !! Houses row data class in the form of distinct sparse vectors

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_types               ,only: RealArray ,IntArray

    implicit none
    private

    type ,public ,extends(Object) :: SparseRowData
        !! Sparse row matrix representation

        integer(ik)     ,allocatable ,dimension(:) ,public :: rowIndex !! Row indices of each sparse row vector
        type(IntArray)  ,allocatable ,dimension(:) ,public :: columnVector !! Arrays of column indices for each row
        type(RealArray) ,allocatable ,dimension(:) ,public :: values !! Values associated with each sparse row vector

        contains

        procedure ,public :: addRow 

        procedure ,public :: init => initSparseRowData

    end type SparseRowData

    public :: operator(*)
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initSparseRowData(this,rowIndices,colVectors) 
            !! Sparse row data Object initialization routine. If indices are supplied initializes values to 0.

            class(SparseRowData)                         ,intent(inout)  :: this
            integer(ik)          ,dimension(:) ,optional ,intent(in)     :: rowIndices !! Row indices of each sparse row vector
            type(IntArray)       ,dimension(:) ,optional ,intent(in)     :: colVectors !! Arrays of column indices for each row

        end subroutine initSparseRowData
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine addRow(this,rowIndex,columnIndices,values) 
            !!  Add sparse row vector to data Object and initialize its values to 0, unless value vector is supplied

            class(SparseRowData)               ,intent(inout)  :: this
            integer(ik)                        ,intent(in)     :: rowIndex !! Row index of added vector
            integer(ik)          ,dimension(:) ,intent(in)     :: columnIndices !! Column indices of nonzeros in row
            real(rk) ,optional   ,dimension(:) ,intent(in)     :: values !! Values of nonzeros in row

        end subroutine addRow          
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    interface operator (*)
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function multElementwise(mat1,mat2) result(res)

            class(SparseRowData)  ,intent(in)  :: mat1
            class(SparseRowData)  ,intent(in)  :: mat2
            type(SparseRowData)                :: res

        end function multElementwise      
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function multElementwiseRArray(mat1,mat2) result(res)

            type(RealArray) ,dimension(:)  ,intent(in)  :: mat1
            class(SparseRowData)           ,intent(in)  :: mat2
            type(RealArray) ,allocatable ,dimension(:)  :: res

        end function multElementwiseRArray      
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module sparse_row_data_class
!-----------------------------------------------------------------------------------------------------------------------------------
 