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

    use ,intrinsic :: ieee_arithmetic ,only: ieee_is_finite

implicit none

    ! Submodule-scope context shared with the PETSc Fortran convergence callback relRedKSPConverged.
    ! The PETSc Fortran KSPSetConvergenceTest binding only passes an opaque integer context, so linearSolve populates
    ! relRedCtx immediately before KSPSolve (solves are strictly sequential per MPI rank).
    type :: RelRedConvContext
        PetscReal :: relRedTarget = 0.0d0 !! Requested residual reduction factor relative to initial residual ||r_0||
        PetscReal :: exactRtol = 0.0d0    !! Configured exact relative tolerance (this%options%solverToleranceRel)
        PetscReal :: absTol = 0.0d0       !! Configured absolute tolerance floor (this%options%solverToleranceAbs)
        PetscReal :: divTol = 0.0d0       !! Configured divergence factor (this%options%solverToleranceDiv)
        PetscReal :: initRnorm = 0.0d0    !! Initial KSP residual norm ||r_0|| recorded at callback iteration n == 0
        PetscReal :: rtolFloor = 0.0d0    !! Exact relative floor exactRtol * ||b_norm|| computed at n == 0
        PetscReal :: divFloor = 0.0d0     !! Divergence threshold divTol * max(||b_norm||, ||r_0||) computed at n == 0
        PetscInt  :: maxIts = 0           !! Configured KSP iteration cap (this%options%maxSolverIters)
        logical   :: hitExact = .false.   !! Set true by relRedKSPConverged iff absTol or rtolFloor was met
        Vec       :: workVec              !! Scratch PETSc Vec (this%objs(usedGroup)%workVec) for P^-1 * b at n == 0
    end type RelRedConvContext

    type(RelRedConvContext) :: relRedCtx

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
module subroutine addFlatValuesToMatrix(this,bufferIndices,values,multConst,objGroup)
    !! Scatter contiguous matrix values into the flat COO value buffer at 1-based bufferIndices,
    !! optionally scaled by a real constant

    class(PETScController)      ,intent(inout)  :: this
    integer(ik)   ,dimension(:) ,intent(in)     :: bufferIndices !! 1-based indices into the group COO value buffer
    real(rk)      ,dimension(:) ,intent(in)     :: values !! Contiguous nonzero values aligned with bufferIndices
    real(rk)      ,optional     ,intent(in)     :: multConst !! Optional multiplicative scaling factor (default 1.0)
    integer(ik)   ,optional     ,intent(in)     :: objGroup !! Optional PETSc object group index (default 1)

    real(rk)    :: mult
    integer(ik) :: i ,usedGroup

    if (assertions .or. assertionLvl >= 0) then
        call assert(this%isDefined(),"Attempted to add flat values to matrix of undefined PETSc controller object")
        call assert(this%objCreated,&
        "Attempted to add flat values to matrix of PETSc controller before PETSc objects were generated")
        call assert(size(bufferIndices) == size(values),&
        "Flat value buffer index and value array sizes do not match in addFlatValuesToMatrix")
    end if

    mult = real(1.0d0,kind=rk)
    usedGroup = 1
    if (present(objGroup)) usedGroup = objGroup
    if (present(multConst)) mult = multConst

    associate(cooVals => this%objs(usedGroup)%cooValues)
        do i = 1, size(bufferIndices)
            cooVals(bufferIndices(i)) = cooVals(bufferIndices(i)) + mult * values(i)
        end do
    end associate

end subroutine addFlatValuesToMatrix
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

        call VecDuplicate(this%objs(i)%sol,this%objs(i)%workVec,ierr)
        CHKERRQ(ierr)
        this%objs(i)%usingRelRedConv = .false.

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
module subroutine linearSolve(this,knownVec,unknownVec,objGroup,reusePC,rtol,usedExactTol)
    !! Solve the equation petscMat * unknownVec = knownVec where the vectors are local to each processor, using unknownVec as
    !! the initial guess - assembles and zeros out the matrix before/after the solve.
    !! When optional rtol > solverToleranceRel is passed for an iterative KSP, installs relRedKSPConverged to stop once the
    !! residual drops by rtol relative to ||r_0|| (floored at solverToleranceRel*||b_norm|| and solverToleranceAbs), and
    !! reports via usedExactTol whether that exact floor was satisfied.

    class(PETScController)      ,intent(inout)  :: this
    real(rk)      ,dimension(:) ,intent(in)     :: knownVec 
    real(rk)      ,dimension(:) ,intent(inout)  :: unknownVec
    integer(ik)   ,optional     ,intent(in)     :: objGroup
    logical       ,optional     ,intent(in)     :: reusePC !! Default false; overrides -ksp_reuse_preconditioner
    real(rk)      ,optional     ,intent(in)     :: rtol !! Optional relative residual reduction vs ||r_0||
    logical       ,optional     ,intent(out)    :: usedExactTol !! True iff solverToleranceRel/solverToleranceAbs was met

    external :: KSPConvergedDefault ,KSPConvergedDefaultDestroy
    integer(ik) :: ierr ,usedGroup ,attempt
    logical     :: useRelRed ,exactTol
    PetscBool :: reuseThisSolve
    KSPType :: solverType
    PetscFortranAddr :: defCtx ,convDummy

    PetscScalar ,dimension(:) ,pointer :: vecPointer ,solutionPointer

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to solve linear system using undefined PETSc controller")
        call assert(this%objCreated,"Attempted to solve linear system without first creating PETSc objects")
    end if

    usedGroup = 1
    if (present(objGroup)) usedGroup = objGroup

    ! Only enable relative-reduction early stopping when the requested rtol is strictly looser than solverToleranceRel.
    useRelRed = .false.
    if (present(rtol)) useRelRed = (rtol > this%options%solverToleranceRel)

    reuseThisSolve = PETSC_FALSE
    if (present(reusePC)) reuseThisSolve = reusePC
    if (reuseThisSolve .or. useRelRed) then
        ! Direct solvers (KSPPREONLY, e.g. LU/MUMPS) always solve directly in one application: never reuse a stale
        ! factorization across matrix updates and never install an iterative relative-reduction convergence test.
        call KSPGetType(this%objs(usedGroup)%solver,solverType,ierr)
        CHKERRQ(ierr)
        if (trim(solverType) == KSPPREONLY) then
            reuseThisSolve = PETSC_FALSE
            useRelRed = .false.
        end if
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

    if (useRelRed) then
        ! Populate submodule-scope callback context for this solve before entering KSPSolve.
        relRedCtx%relRedTarget = real(rtol,kind=kind(relRedCtx%relRedTarget))
        relRedCtx%exactRtol = real(this%options%solverToleranceRel,kind=kind(relRedCtx%exactRtol))
        relRedCtx%absTol = real(this%options%solverToleranceAbs,kind=kind(relRedCtx%absTol))
        relRedCtx%divTol = real(this%options%solverToleranceDiv,kind=kind(relRedCtx%divTol))
        relRedCtx%maxIts = int(this%options%maxSolverIters,kind=kind(relRedCtx%maxIts))
        relRedCtx%workVec = this%objs(usedGroup)%workVec
        relRedCtx%hitExact = .false.
        if (.not. this%objs(usedGroup)%usingRelRedConv) then
            ! Install custom convergence test once per KSP until an exact solve restores KSPConvergedDefault.
            convDummy = 0
            call KSPSetConvergenceTest(this%objs(usedGroup)%solver,relRedKSPConverged,convDummy,PETSC_NULL_FUNCTION,ierr)
            CHKERRQ(ierr)
            this%objs(usedGroup)%usingRelRedConv = .true.
        end if
    else if (this%objs(usedGroup)%usingRelRedConv) then
        ! Restore the PETSc default convergence test and attach KSPConvergedDefaultDestroy to own defCtx cleanup.
        call KSPConvergedDefaultCreate(defCtx,ierr)
        CHKERRQ(ierr)
        call KSPSetConvergenceTest(this%objs(usedGroup)%solver,KSPConvergedDefault,defCtx,KSPConvergedDefaultDestroy,ierr)
        CHKERRQ(ierr)
        this%objs(usedGroup)%usingRelRedConv = .false.
    end if

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

    ! Exact tolerance is satisfied either when inexact mode was not active (exact KSP / KSPPREONLY) or when the
    ! relative-reduction callback terminated on absTol or rtolFloor (relRedCtx%hitExact == .true.).
    exactTol = (.not. useRelRed) .or. relRedCtx%hitExact
    if (present(usedExactTol)) usedExactTol = exactTol

    call VecGetArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
    CHKERRQ(ierr)

    unknownVec = real(solutionPointer,kind=rk)

    call VecRestoreArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
    CHKERRQ(ierr)

end subroutine linearSolve
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine relRedKSPConverged(ksp,n,rnorm,flag,dummy,ierr)
    !! Custom PETSc KSP convergence callback implementing relative-reduction early stopping with an exact tolerance floor.
    !!
    !! Stopping criterion at Krylov iteration n with PETSc-supplied residual norm rnorm:
    !!   1. Exact ATOL:   rnorm <= relRedCtx%absTol                             -> KSP_CONVERGED_ATOL, hitExact = .true.
    !!   2. Exact RTOL:   rnorm <= relRedCtx%exactRtol * ||b_norm||             -> KSP_CONVERGED_RTOL, hitExact = .true.
    !!   3. Inexact stop: n > 0 and rnorm <= relRedCtx%relRedTarget * ||r_0||   -> KSP_CONVERGED_RTOL, hitExact = .false.
    !!   4. Divergence:   NaN/Inf rnorm, or n > 0 and rnorm >= divFloor         -> diverged
    !!   5. Max iters:    n >= relRedCtx%maxIts                                 -> KSP_DIVERGED_ITS
    !!
    !! Here ||b_norm|| matches the active PETSc KSPNormType (||P^-1 b||_2 for KSP_NORM_PRECONDITIONED, else ||b||_2).
    !! On a fresh-PC retry inside linearSolve, PETSc invokes n == 0 again with the rebuilt PC, recomputing all floors.

    KSP                ,intent(in)  :: ksp   !! Active PETSc KSP solver handle
    PetscInt           ,intent(in)  :: n     !! Current Krylov iteration count (0 before first Krylov step)
    PetscReal          ,intent(in)  :: rnorm !! Current residual norm reported by KSP under its active norm type
    KSPConvergedReason ,intent(out) :: flag  !! Convergence/divergence status returned to PETSc
    PetscFortranAddr   ,intent(in)  :: dummy !! Unused Fortran callback context (required by PETSc signature)
    PetscErrorCode     ,intent(out) :: ierr  !! PETSc error code

    KSPNormType :: normType
    PC          :: pc
    Vec         :: rhsVec
    PetscReal   :: bnorm

    ! Zero ierr while referencing dummy to satisfy -Wunused-dummy-argument.
    ierr = int(0_8 * dummy, kind=kind(ierr))
    if (.not. ieee_is_finite(real(rnorm,kind=rk))) then
        flag = KSP_DIVERGED_NANORINF
        return
    end if

    if (n == 0) then
        ! Initialize per-solve reference norms using the active PC and RHS matching the PETSc residual norm type.
        relRedCtx%hitExact = .false.
        call KSPGetNormType(ksp,normType,ierr)
        CHKERRQ(ierr)
        call KSPGetRhs(ksp,rhsVec,ierr)
        CHKERRQ(ierr)
        if (normType == KSP_NORM_PRECONDITIONED) then
            call KSPGetPC(ksp,pc,ierr)
            CHKERRQ(ierr)
            call PCApply(pc,rhsVec,relRedCtx%workVec,ierr)
            CHKERRQ(ierr)
            call VecNorm(relRedCtx%workVec,NORM_2,bnorm,ierr)
            CHKERRQ(ierr)
        else
            call VecNorm(rhsVec,NORM_2,bnorm,ierr)
            CHKERRQ(ierr)
        end if
        relRedCtx%rtolFloor = relRedCtx%exactRtol * bnorm
        relRedCtx%divFloor = relRedCtx%divTol * max(bnorm, rnorm)
        relRedCtx%initRnorm = rnorm
    end if

    ! Test exact thresholds first (ATOL before RTOL to match PETSc KSPConvergedDefault) so a solve that naturally reaches
    ! solverToleranceAbs or solverToleranceRel (including at n == 0 when the initial guess is converged) sets hitExact = .true.
    if (rnorm <= relRedCtx%absTol) then
        flag = KSP_CONVERGED_ATOL
        relRedCtx%hitExact = .true.
    else if (rnorm <= relRedCtx%rtolFloor) then
        flag = KSP_CONVERGED_RTOL
        relRedCtx%hitExact = .true.
    else if (n > 0 .and. rnorm <= relRedCtx%relRedTarget * relRedCtx%initRnorm) then
        flag = KSP_CONVERGED_RTOL
        relRedCtx%hitExact = .false.
    else if (n > 0 .and. rnorm >= relRedCtx%divFloor) then
        flag = KSP_DIVERGED_DTOL
    else if (n >= relRedCtx%maxIts) then
        flag = KSP_DIVERGED_ITS
    else
        flag = KSP_CONVERGED_ITERATING
    end if

end subroutine relRedKSPConverged
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
            call VecDestroy(this%objs(i)%workVec,ierr)
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
