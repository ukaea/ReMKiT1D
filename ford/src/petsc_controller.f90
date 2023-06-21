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
module petsc_controller_class
    !! Houses class in charge of interfacing with PETSc

#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"

    use data_kinds                       ,only: rk, ik
    use runtime_constants                ,only: debugging, assertions
    use god_objects                      ,only: Object
    use assertion_utility                ,only: assert, assertIdentical, assertPure
    use petsc_preallocation_data_class   ,only: PetscPreallocationData
    use mpi_controller_class             ,only: MPIController
    use sparse_row_data_class            ,only: SparseRowData
    use indexing_class                   ,only: Indexing
    use petscksp

    implicit none
    private

    type ,public :: SolverOptions
        !! Contains data used by the PETSc ksp solver to set various options (refer to PETSc manual)
        character(:) ,allocatable                 :: kspSolverType   !! String representing the ksp solver type
        real(rk)                                  :: solverToleranceRel = real(1.0d-17,kind=rk) !! Relative solver tolerance
        real(rk)                                  :: solverToleranceAbs = real(1.0d-20,kind=rk) !! Absolute solver tolerance
        real(rk)                                  :: solverToleranceDiv = real(1.0d07,kind=rk) !! Solver divergence tolerance
        integer(ik)                               :: maxSolverIters  = 10000 !! Maximum number of solver iterations
        character(:) ,allocatable                 :: hyprePC   !! String representing the hypre PC type
        character(:) ,allocatable                 :: petscOptions  !! String representing command line style options

    end type

    type ,public :: PETSsObjs

        Vec                          ,private :: rhs !! RHS PETSc vector object in linear solve
        Vec                          ,private :: sol !! LHS PETSc vector object in linear solve
        Mat                          ,private :: petscMat !! PETSc matrix object used in linear solve
        PC                           ,private :: preconditioner !! Precondition object used in solve 
        KSP                          ,private :: solver !! KSP solver object

    end type 

    type ,public ,extends(Object) :: PETScController
        !! Object responsible for interfacing with PETSc

        type(PetscPreallocationData) ,private :: preallocData !! Data used for preallocating PETSc objects
        type(SolverOptions)          ,private :: options !! Solver options for the KSP solver
        logical                      ,private :: objCreated  !! True if PETSc object components have been created 
        
        type(PETSsObjs) ,allocatable ,dimension(:) ,private :: objs !! Groups of PETSc objects

        KSPConvergedReason                 :: lastConvergedReason !! Converged reason for last attempted solve
        PetscInt                           :: lastNumIterations !! Number of iterations for last attemtped solve

        contains

        procedure ,public :: addRowDataToPreallocation
        procedure ,public :: createPETScObjs 
        procedure ,public :: addRowValuesToMatrix 
        procedure ,public :: linearSolve

        procedure ,public :: objectsCreated

        procedure ,public :: getLastConvergedReason
        procedure ,public :: getLastNumIterations

        procedure ,public :: init => initPETScController

        procedure ,public :: finalize


    end type PETScController
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initPETScController(this,indexingObj,mpiCont,solOptions,numObjs) 
        !! PETSc controller initialization routine

        class(PETScController) ,intent(inout)  :: this
        type(Indexing)         ,intent(in)     :: indexingObj !! Indexing object used to initialize preallocation data
        type(MPIController)    ,intent(in)     :: mpiCont !! Reference MPI controller
        type(SolverOptions)    ,intent(in)     :: solOptions !! KSP solver options
        integer(ik)  ,optional ,intent(in)     :: numObjs !! Number of PETSc object groups

    end subroutine initPETScController
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addRowDataToPreallocation(this,rowData) 
        !! Add sparse row data structure to total preallocation sparsity pattern

        class(PETScController)      ,intent(inout)  :: this
        type(SparseRowData)         ,intent(in)     :: rowData 

    end subroutine addRowDataToPreallocation
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addRowValuesToMatrix(this,rowData,multConst,objGroup) 
        !! Add row data to matrix, optionally multiplied by a real constant

        class(PETScController)      ,intent(inout)  :: this
        type(SparseRowData)         ,intent(in)     :: rowData 
        real(rk)      ,optional     ,intent(in)     :: multConst
        integer(ik)   ,optional     ,intent(in)     :: objGroup

    end subroutine addRowValuesToMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine createPETScObjs(this) 
        !! Assemble the preallocation object, create PETSc objects, and deallocate the preallocation object

        class(PETScController)      ,intent(inout)  :: this

    end subroutine createPETScObjs
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function objectsCreated(this) result(created)
        !! Return true if PETSc objects have been created

        class(PETScController)      ,intent(in)  :: this
        logical                                  :: created

    end function objectsCreated
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine linearSolve(this,knownVec,unknownVec,objGroup) 
        !! Solve the equation petscMat * unknownVec = knownVec where the vectors are local to each processor - assembles and zeros out the
        !! matrix before/after the solve

        class(PETScController)      ,intent(inout)  :: this
        real(rk)      ,dimension(:) ,intent(in)     :: knownVec 
        real(rk)      ,dimension(:) ,intent(out)    :: unknownVec
        integer(ik)   ,optional     ,intent(in)     :: objGroup

    end subroutine linearSolve
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getLastConvergedReason(this) result(conv)
        !! Getter for lastConvergedReason

        class(PETScController)      ,intent(in)  :: this
        integer(ik)                              :: conv

    end function getLastConvergedReason
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getLastNumIterations(this) result(numIters)
        !! Getter for lastNumIterations

        class(PETScController)      ,intent(in)  :: this
        integer(ik)                              :: numIters

    end function getLastNumIterations
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine finalize(this) 
        !! Finalizes PETSc and destroys all PETSc objects

        class(PETScController)      ,intent(inout)  :: this

    end subroutine finalize
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module petsc_controller_class
!-----------------------------------------------------------------------------------------------------------------------------------
 