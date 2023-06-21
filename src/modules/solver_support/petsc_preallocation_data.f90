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
module petsc_preallocation_data_class
    !! Houses data container class used in PETSc matrix preallocation

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_types               ,only: RealArray ,IntArray
    use indexing_class              ,only: Indexing 
    use coo_sparsity_pattern_class  ,only: COOSparsityPattern
    use sparse_row_data_class       ,only: SparseRowData

    implicit none
    private

    type ,public ,extends(Object) :: PETScPreallocationData
        !! Object used to properly preallocate PETSc matrices

        type(COOSparsityPattern) ,allocatable ,private :: totalCOOPattern !! Total (across matrix terms) sparsity pattern
        integer(ik)                           ,private :: localDoFOffset !! Local process offset of degrees of freedom
        integer(ik)                           ,private :: localNumDoFs !! Local number of degrees of freedom
        integer(ik)                           ,private :: totalNumDoFs !! Total number of degrees of freedom

        integer(ik) ,allocatable ,dimension(:) ,private :: rowNumNonzeros !! Number of nonzeros in each matrix row
        integer(ik) ,allocatable ,dimension(:) ,private :: rowNumNonzerosDiag !! Number of nonzeros in each row belonging to the local diagonal block
        integer(ik) ,allocatable ,dimension(:) ,private :: rowNumNonzerosOffDiag !! Number of nonzeros in each row not in the local diagonal block

        logical ,private :: assembled !! True if preallocation data is assembled and ready for processing by PETScController routines

        contains

        procedure ,public :: addRowDataToPattern
        procedure ,public :: assembleData 
        procedure ,public :: isAssembled

        procedure ,public :: getNumNonzerosDiag
        procedure ,public :: getNumNonzerosOffDiag

        procedure ,public :: getLocalNumDoFs
        procedure ,public :: getTotalNumDoFs
        procedure ,public :: getLocalDoFOffset

        procedure ,public :: deallocatePattern

        procedure ,public :: init => initPetscPreallocationData

    end type PETScPreallocationData
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initPetscPreallocationData(this,indexingObj,procRank) 
            !! Petsc preallocation data object initialization routine

            class(PETScPreallocationData) ,intent(inout)  :: this
            type(Indexing)                ,intent(in)     :: indexingObj !! Indexing object to retrieve DoF data from
            integer(ik)                   ,intent(in)     :: procRank !! Current process rank

        end subroutine initPetscPreallocationData
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine addRowDataToPattern(this,rowData) 
            !! Add sparse row data structure to total sparsity pattern 

            class(PETScPreallocationData)      ,intent(inout)  :: this
            type(SparseRowData)                ,intent(in)     :: rowData 

        end subroutine addRowDataToPattern
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine assembleData(this) 
            !! Assemble preallocation data

            class(PETScPreallocationData)      ,intent(inout)  :: this

        end subroutine assembleData
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isAssembled (this) result(assembled)
            !! Return assembly status of preallocation data

            class(PETScPreallocationData) ,intent(in) :: this
            logical                                   :: assembled
 
        end function isAssembled
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumNonzerosDiag (this) result(numNZDiag)
            !! Return number of nonzeros in the diagonal part of each row

            class(PETScPreallocationData)          ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)             :: numNZDiag
 
        end function getNumNonzerosDiag
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumNonzerosOffDiag (this) result(numNZODiag)
            !! Return number of nonzeros in the off-diagonal part of each row

            class(PETScPreallocationData)          ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)             :: numNZODiag
 
        end function getNumNonzerosOffDiag
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getLocalNumDoFs (this) result(locNumDoFs)
            !! Return local number of DoFs in unknown vector

            class(PETScPreallocationData) ,intent(in) :: this
            integer(ik)                               :: locNumDoFs
 
        end function getLocalNumDoFs
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getTotalNumDoFs (this) result(totalNumDoFs)
            !! Return total number of DoFs in unknown vector

            class(PETScPreallocationData) ,intent(in) :: this
            integer(ik)                               :: totalNumDoFs
 
        end function getTotalNumDoFs
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getLocalDoFOffset (this) result(offset)
            !! Return offset of the start of this processors local vector in the global vector

            class(PETScPreallocationData) ,intent(in) :: this
            integer(ik)                               :: offset
 
        end function getLocalDoFOffset
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine deallocatePattern(this)
            !!  Deallocates the COO matrix pattern

            class(PETScPreallocationData)      ,intent(inout)  :: this

        end subroutine deallocatePattern
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module petsc_preallocation_data_class
!-----------------------------------------------------------------------------------------------------------------------------------
 