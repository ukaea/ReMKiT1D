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
submodule (petsc_controller_class) petsc_controller_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the petsc controller class

#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initPETScController(this,indexingObj,mpiCont,solOptions,numObjs) 
    !! PETSc controller initialization routine

    class(PETScController) ,intent(inout)  :: this
    type(Indexing)         ,intent(in)     :: indexingObj !! Indexing object used to initialize preallocation data
    type(MPIController)    ,intent(in)     :: mpiCont !! Reference MPI controller
    type(SolverOptions)    ,intent(in)     :: solOptions !! KSP solver options
    integer(ik)  ,optional ,intent(in)     :: numObjs !! Number of PETSc object groups

    integer(ik)                            :: ierr 

    if (assertions) then 
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to PETSc controller constructor")
        call assert(indexingObj%isDefined(),"Undefined MPI controller passed to PETSc controller constructor")
    end if

    call this%makeDefined()
    call this%preallocData%init(indexingObj,mpiCont%getWorldRank())
    this%options = solOptions
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    if (present(numObjs)) then 
        allocate(this%objs(numObjs))
    else
        allocate(this%objs(1))
    end if
    this%objCreated = .false. 

end subroutine initPETScController
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addRowDataToPreallocation(this,rowData) 
    !! Add sparse row data structure to total preallocation sparsity pattern

    class(PETScController)      ,intent(inout)  :: this
    type(SparseRowData)         ,intent(inout)  :: rowData

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add row data structure to PETSc preallocation data using &
        &undefined PETSc controller object")
        call assertPure(rowData%isDefined(),"Attempted to add undefined row data structure to PETSc preallocation data object")
    end if

    call this%preallocData%addRowDataToPattern(rowData)

end subroutine addRowDataToPreallocation
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addRowValuesToMatrix(this,rowData,multConst,objGroup) 
    !! Add row data to matrix, optionally multiplied by a real constant

    class(PETScController)      ,intent(inout)  :: this
    type(SparseRowData)         ,intent(in)     :: rowData 
    real(rk)      ,optional     ,intent(in)     :: multConst
    integer(ik)   ,optional     ,intent(in)     :: objGroup

    real(rk) :: mult

    integer(ik) :: i ,j ,usedGroup

    if (assertions .or. assertionLvl >= 0) then
        call assert(this%isDefined(),"Attempted to add row values to matrix of undefined PETSc controller object")
        call assert(this%objCreated,&
        "Attempted to add row values to matrix of PETSc controller before PETSc objects were generated")
        call assert(allocated(rowData%cooIndices),&
        "Attempted to add row values to PETSc matrix from SparseRowData not registered during preallocation")
        call assert(size(rowData%cooIndices) == size(rowData%rowIndex),&
        "SparseRowData cooIndices size mismatch - row structure modified after preallocation registration")
    end if

    mult = real(1.0d0,kind=rk)
    usedGroup = 1 
    if (present(objGroup)) usedGroup = objGroup
    if (present(multConst)) mult = multConst 

    do i = 1,size(rowData%rowIndex)
        do j = 1,size(rowData%cooIndices(i)%entry)
            this%objs(usedGroup)%cooValues(rowData%cooIndices(i)%entry(j)) = &
                this%objs(usedGroup)%cooValues(rowData%cooIndices(i)%entry(j)) + mult * rowData%values(i)%entry(j)
        end do
    end do

end subroutine addRowValuesToMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine createPETScObjs(this) 
    !! Assemble the preallocation object, create PETSc objects, and deallocate the preallocation object

    class(PETScController)      ,intent(inout)  :: this

    integer(ik) :: ierr ,i ,nnzCOO
    PetscCount  :: ncoo
    PetscInt ,allocatable ,dimension(:) :: cooRows0 ,cooCols0
    KSPType     :: kspType

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to create PETSc objects in undefined PETSc controller")
        call assert(.not. this%objCreated,"Attempted to create PETSc objects when they are already created")
    end if

    call this%preallocData%assembleData()
    nnzCOO = this%preallocData%getNumNonzerosCOO()
    ncoo = int(nnzCOO,kind=kind(ncoo))
    cooRows0 = this%preallocData%getCOORowIndices0()
    cooCols0 = this%preallocData%getCOOColIndices0()

    do i = 1, size(this%objs)

        call PETScOptionsInsertString(PETSC_NULL_OPTIONS,this%options%petscOptions,ierr)
        CHKERRQ(ierr)

        call MatCreate(PETSc_Comm_World,this%objs(i)%petscMat,ierr)
        CHKERRQ(ierr)
        call MatSetSizes(this%objs(i)%petscMat,this%preallocData%getLocalNumDoFs(),this%preallocData%getLocalNumDoFs(),&
                                        this%preallocData%getTotalNumDoFs(),this%preallocData%getTotalNumDoFs(),ierr)
        CHKERRQ(ierr)
        call MatSetType(this%objs(i)%petscMat,MATAIJ,ierr)
        CHKERRQ(ierr)
        call MatSetFromOptions(this%objs(i)%petscMat,ierr)
        CHKERRQ(ierr)
        call MatSetOption(this%objs(i)%petscMat,MAT_NO_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
        CHKERRQ(ierr)
        call MatSetPreallocationCOO(this%objs(i)%petscMat,ncoo,cooRows0,cooCols0,ierr)
        CHKERRQ(ierr)

        allocate(this%objs(i)%cooValues(max(1,nnzCOO)))
        this%objs(i)%cooValues = 0.0d0

        call VecCreateMPI(PETSc_Comm_World,this%preallocData%getLocalNumDoFs(),&
                                        this%preallocData%getTotalNumDoFs(),&
                                        this%objs(i)%sol,ierr)
        CHKERRQ(ierr)

        call VecDuplicate(this%objs(i)%sol,this%objs(i)%rhs,ierr)
        CHKERRQ(ierr)

        call KSPCreate(PETSc_Comm_World,this%objs(i)%solver,ierr)
        CHKERRQ(ierr)

        call KSPSetType(this%objs(i)%solver,this%options%kspSolverType,ierr)
        CHKERRQ(ierr)

        call KSPSetInitialGuessNonzero(this%objs(i)%solver,PETSC_TRUE,ierr)
        CHKERRQ(ierr)

        call KSPSetFromOptions(this%objs(i)%solver,ierr)
        CHKERRQ(ierr)

        call KSPGetType(this%objs(i)%solver,kspType,ierr)
        CHKERRQ(ierr)

        if (trim(kspType) == KSPPREONLY) then
            call KSPSetInitialGuessNonzero(this%objs(i)%solver,PETSC_FALSE,ierr)
            CHKERRQ(ierr)
        end if

        call KSPGetPC(this%objs(i)%solver,this%objs(i)%preconditioner,ierr)
        CHKERRQ(ierr)

        if (len(this%options%hyprePC) > 0) then
            call PCSetType(this%objs(i)%preconditioner,PCHYPRE,ierr)
            CHKERRQ(ierr)

            call PCHYPRESetType(this%objs(i)%preconditioner,this%options%hyprePC,ierr)
            CHKERRQ(ierr)

        else
            call PCSetFromOptions(this%objs(i)%preconditioner,ierr)
            CHKERRQ(ierr)
        end if
    end do

    this%objCreated = .true.

    call this%preallocData%deallocatePattern()

end subroutine createPETScObjs
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function objectsCreated(this) result(created)
    !! Return true if PETSc objects have been created

    class(PETScController)      ,intent(in)  :: this
    logical                                  :: created

    if (assertions) call assertPure(this%isDefined(),"Attempted to get object creation status form undefined PETSc controller")

    created = this%objCreated

end function objectsCreated
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine linearSolve(this,knownVec,unknownVec,objGroup,reusePC)
    !! Solve the equation petscMat * unknownVec = knownVec where the vectors are local to each processor, using unknownVec as
    !! the initial guess - assembles and zeros out the matrix before/after the solve

    class(PETScController)      ,intent(inout)  :: this
    real(rk)      ,dimension(:) ,intent(in)     :: knownVec 
    real(rk)      ,dimension(:) ,intent(inout)  :: unknownVec
    integer(ik)   ,optional     ,intent(in)     :: objGroup
    logical       ,optional     ,intent(in)     :: reusePC !! Default false; overrides -ksp_reuse_preconditioner

    integer(ik) :: ierr ,usedGroup ,attempt
    PetscBool :: reuseThisSolve
    KSPType :: solverType

    PetscScalar ,dimension(:) ,pointer :: vecPointer ,solutionPointer

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to solve linear system using undefined PETSc controller")
        call assert(this%objCreated,"Attempted to solve linear system without first creating PETSc objects")
    end if

    usedGroup = 1
    if (present(objGroup)) usedGroup = objGroup

    reuseThisSolve = PETSC_FALSE
    if (present(reusePC)) reuseThisSolve = reusePC
    if (reuseThisSolve) then
        call KSPGetType(this%objs(usedGroup)%solver,solverType,ierr)
        CHKERRQ(ierr)
        if (trim(solverType) == KSPPREONLY) &
            error stop "PC reuse requires an iterative KSP: preonly cannot correct stale factors for a changed operator"
    end if

    call MatSetValuesCOO(this%objs(usedGroup)%petscMat,this%objs(usedGroup)%cooValues,INSERT_VALUES,ierr)
    CHKERRQ(ierr)
    this%objs(usedGroup)%cooValues = 0.0d0

    call VecGetArrayF90(this%objs(usedGroup)%rhs,vecPointer,ierr)
    CHKERRQ(ierr)
    vecPointer = knownVec
    call VecRestoreArrayF90(this%objs(usedGroup)%rhs,vecPointer,ierr)
    CHKERRQ(ierr)

    call KSPSetTolerances(this%objs(usedGroup)%solver,this%options%solverToleranceRel,&
                                     this%options%solverToleranceAbs,&
                                     this%options%solverToleranceDiv,&
                                     this%options%maxSolverIters,ierr)
    CHKERRQ(ierr)


    call KSPSetOperators(this%objs(usedGroup)%solver,this%objs(usedGroup)%petscMat,this%objs(usedGroup)%petscMat,ierr)
    CHKERRQ(ierr)

    do attempt = 1,2
        ! Set explicitly on every solve so reuse cannot leak to another caller of this group.
        call KSPSetReusePreconditioner(this%objs(usedGroup)%solver,reuseThisSolve,ierr)
        CHKERRQ(ierr)

        ! Keep the initial guess unchanged until the final result is extracted.
        call VecGetArrayF90(this%objs(usedGroup)%sol,vecPointer,ierr)
        CHKERRQ(ierr)
        vecPointer = unknownVec
        call VecRestoreArrayF90(this%objs(usedGroup)%sol,vecPointer,ierr)
        CHKERRQ(ierr)

        call KSPSolve(this%objs(usedGroup)%solver,this%objs(usedGroup)%rhs,this%objs(usedGroup)%sol,ierr)
        CHKERRQ(ierr)

        call KSPGetConvergedReason(this%objs(usedGroup)%solver,this%lastConvergedReason,ierr)
        CHKERRQ(ierr)
        call KSPGetIterationNumber(this%objs(usedGroup)%solver,this%lastNumIterations,ierr)
        CHKERRQ(ierr)

        if (this%lastConvergedReason >= 0 .or. .not. reuseThisSolve) exit
        ! Retry the installed matrix once, without reassembly or consuming the cleared COO buffer.
        reuseThisSolve = PETSC_FALSE
        call PetscPrintf(PETSC_COMM_WORLD,"Reused PC solve diverged; refreshing PC and retrying once."//achar(10),ierr)
        CHKERRQ(ierr)
    end do

    call VecGetArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
    CHKERRQ(ierr)

    unknownVec = real(solutionPointer,kind=rk)

    call VecRestoreArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
    CHKERRQ(ierr)

end subroutine linearSolve
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLastConvergedReason(this) result(conv)
    !! Getter for lastConvergedReason

    class(PETScController)      ,intent(in)  :: this
    integer(ik)                              :: conv

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to get converged reason from undefined PETSc controller")
        call assertPure(this%objCreated,"Attempted to get converged reason from PETSc controller before any PETSc objects created")
    end if

    conv = int(this%lastConvergedReason,kind=ik)

end function getLastConvergedReason
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLastNumIterations(this) result(numIters)
    !! Getter for lastNumIterations

    class(PETScController)      ,intent(in)  :: this
    integer(ik)                              :: numIters

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to get number of iterations from undefined PETSc controller")
        call assertPure(this%objCreated,&
        "Attempted to get number of iterations from PETSc controller before any PETSc objects created")
    end if

    numIters = int(this%lastNumIterations,kind=ik)

end function getLastNumIterations
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine finalize(this) 
    !! Finalizes PETSc and destroys all PETSc objects

    class(PETScController)      ,intent(inout)  :: this

    integer(ik) :: ierr ,i 

    if (this%isDefined()) then 

        do i = 1,size(this%objs)
            call VecDestroy(this%objs(i)%sol,ierr)
            CHKERRQ(ierr)
            call VecDestroy(this%objs(i)%rhs,ierr)
            CHKERRQ(ierr)
            call MatDestroy(this%objs(i)%petscMat,ierr)
            CHKERRQ(ierr)
            call KSPDestroy(this%objs(i)%solver,ierr)
            CHKERRQ(ierr)
        end do
        call PetscFinalize(ierr)
        CHKERRQ(ierr)
    end if

end subroutine finalize
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule petsc_controller_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
