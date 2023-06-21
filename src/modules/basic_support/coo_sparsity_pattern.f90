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
module coo_sparsity_pattern_class
    !! author: Stefan Mijin 
    !! 
    !! Houses base coordinate list sparsity pattern

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_functions           ,only: findIndices
    use support_types

    implicit none
    private

    type ,public ,extends(Object) :: COOSparsityPattern
        !! Coordinate list format sparsity pattern

        integer(ik) ,public :: numNonZeros !! Total number of non-zero entries
        integer(ik) ,public :: numRows !! Number of matrix rows
        integer(ik) ,public :: numCols !! Number of matrix columns
 
        integer(ik) ,allocatable ,dimension(:) ,public :: rowIndex !! Array of non-zero row indices
        integer(ik) ,allocatable ,dimension(:) ,public :: colIndex !! Array of non-zero column indices

        integer(ik) ,public :: bufferSize
        integer(ik) ,public :: numAddedElements

        logical ,allocatable ,dimension(:) ,private :: addedRow
        logical ,allocatable ,dimension(:) ,private :: addedCol

        integer(ik) ,allocatable ,dimension(:) ,private :: firstInstRow
        integer(ik) ,allocatable ,dimension(:) ,private :: firstInstCol 

        type(IntArray) ,allocatable ,dimension(:) ,private :: rowGlobalIndices !! Element indices of each row

        integer(ik) ,private :: rowOffset

        contains

        procedure ,public :: hasIndices
        procedure ,public :: findLocationOfIndices
        procedure ,public :: addEntry
        procedure ,public :: addPattern

        procedure ,public :: init => initCooSparsityPattern

    end type COOSparsityPattern
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initCooSparsityPattern(this,numRows,numCols,bufferSize,rowOffset) 
            !! COO sparsity pattern initialization routine

            class(COOSparsityPattern)           ,intent(inout)  :: this
            integer(ik)                         ,intent(in)     :: numRows !! Matrix number of rows
            integer(ik)                         ,intent(in)     :: numCols !! Matrix number of columns
            integer(ik) ,optional               ,intent(in)     :: bufferSize !! Preallocated rolling buffer size
            integer(ik) ,optional               ,intent(in)     :: rowOffset !! Row offset for distributed arrays

        end subroutine initCooSparsityPattern
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function hasIndices (this,row,column) result(found)
            !! Return true if all elements with indices (row,column) are present in sparsity pattern

            class(COOSparsityPattern)   ,intent(in) :: this
            integer(ik)   ,dimension(:) ,intent(in) :: row !! Row indices
            integer(ik)   ,dimension(:) ,intent(in) :: column !! Column indices
            logical                                 :: found
 
        end function hasIndices
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function findLocationOfIndices (this,row,column) result(indices)
            !! Return indices with given values of (row,column)

            class(COOSparsityPattern)                ,intent(in) :: this
            integer(ik)                ,dimension(:) ,intent(in) :: row !! Row index
            integer(ik)                ,dimension(:) ,intent(in) :: column !! Column index
            integer(ik)   ,allocatable ,dimension(:)             :: indices
 
        end function findLocationOfIndices
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine addEntry(this,row,column) 
            !! Add entry if not already present

            class(COOSparsityPattern)   ,intent(inout)  :: this
            integer(ik)                 ,intent(in)     :: row !! Row index
            integer(ik)                 ,intent(in)     :: column !! Column Index

        end subroutine addEntry
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine addPattern(this,pattern,rowOffset,colOffset) 
            !! Add entries from existing pattern that are not already present

            class(COOSparsityPattern)      ,intent(inout)  :: this
            type(COOSparsityPattern)       ,intent(in)     :: pattern !! Pattern to add to this pattern
            integer(ik)      ,dimension(:) ,intent(in)     :: rowOffset !! Array by which to offset pattern row indices before adding
            integer(ik)      ,dimension(:) ,intent(in)     :: colOffset !! Array by which to offset pattern column indices before adding

        end subroutine addPattern
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module coo_sparsity_pattern_class
!-----------------------------------------------------------------------------------------------------------------------------------
 