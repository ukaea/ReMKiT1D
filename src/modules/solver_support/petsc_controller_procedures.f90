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

#include "petscversion.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"

implicit none

    interface
        function kspGetConvergenceTestC(ksp,converge,cctx,destroy) result(ierr) bind(C,name="KSPGetConvergenceTest")
            import :: c_ptr ,c_funptr
            type(c_ptr)    ,value       :: ksp
            type(c_funptr) ,intent(out) :: converge
            type(c_ptr)    ,intent(out) :: cctx
            type(c_funptr) ,intent(out) :: destroy
            PetscErrorCode              :: ierr
        end function kspGetConvergenceTestC

        function kspGetAndClearConvergenceTestC(ksp,converge,cctx,destroy) result(ierr) &
            bind(C,name="KSPGetAndClearConvergenceTest")
            import :: c_ptr ,c_funptr
            type(c_ptr)    ,value       :: ksp
            type(c_funptr) ,intent(out) :: converge
            type(c_ptr)    ,intent(out) :: cctx
            type(c_funptr) ,intent(out) :: destroy
            PetscErrorCode              :: ierr
        end function kspGetAndClearConvergenceTestC

        function kspConvergedDefaultC(ksp,n,rnorm,reason,ctx) result(ierr) bind(C,name="KSPConvergedDefault")
            import :: c_ptr
            type(c_ptr)        ,value       :: ksp
            PetscInt           ,value       :: n
            PetscReal          ,value       :: rnorm
            KSPConvergedReason ,intent(out) :: reason
            type(c_ptr)        ,value       :: ctx
            PetscErrorCode                  :: ierr
        end function kspConvergedDefaultC

        function kspConvergedDefaultDestroyC(ctx) result(ierr) bind(C,name="KSPConvergedDefaultDestroy")
            import :: c_ptr
            type(c_ptr)    ,value :: ctx
            PetscErrorCode        :: ierr
        end function kspConvergedDefaultDestroyC
    end interface

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

    PetscErrorCode                         :: ierr 

    if (assertions) then 
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to PETSc controller constructor")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to PETSc controller constructor")
    end if

    call this%makeDefined()
    call this%preallocData%init(indexingObj,mpiCont%getWorldRank())
    this%options = solOptions
#if PETSC_VERSION_GE(3,19,0)
    call PetscInitialize(ierr)
#else
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#endif
    CHKERRQ(ierr)
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

    PetscErrorCode :: ierr
    integer(ik)    :: i ,nnzCOO
    PetscCount  :: ncoo
    PetscInt ,allocatable ,dimension(:) :: cooRows0 ,cooCols0
    KSPType     :: kspType
    KSPNormType :: normType
    type(c_ptr)    :: kspPtr ,defCtxPtr
    type(c_funptr) :: convFn ,destFn ,defConvFn

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

        this%objs(i)%isPreonly = (trim(kspType) == KSPPREONLY)
        if (this%objs(i)%isPreonly) then
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

        this%objs(i)%allowRelRedConv = .false.
        nullify(this%objs(i)%relRedCtx)
        this%objs(i)%relRedCtxAddr = 0
        if (.not. this%objs(i)%isPreonly) then
            call KSPGetNormType(this%objs(i)%solver,normType,ierr)
            CHKERRQ(ierr)
            if (normType /= KSP_NORM_NONE) then
                kspPtr = transfer(this%objs(i)%solver%v,c_null_ptr)
                ierr = kspGetConvergenceTestC(kspPtr,convFn,defCtxPtr,destFn)
                CHKERRQ(ierr)
                defConvFn = c_funloc(kspConvergedDefaultC)
                if (c_associated(convFn,defConvFn) .and. c_associated(defCtxPtr)) then
                    ! Wrap the configured KSPConvergedDefault context once so all PETSc convergence
                    ! options (norm type, PC side, zero-RHS fallback, maxits behavior) remain PETSc-owned.
                    ierr = kspGetAndClearConvergenceTestC(kspPtr,convFn,defCtxPtr,destFn)
                    CHKERRQ(ierr)
                    allocate(this%objs(i)%relRedCtx)
                    this%objs(i)%relRedCtx%defaultCtxAddr = transfer(defCtxPtr,this%objs(i)%relRedCtx%defaultCtxAddr)
                    this%objs(i)%relRedCtxAddr = transfer(c_loc(this%objs(i)%relRedCtx),this%objs(i)%relRedCtxAddr)
                    call KSPSetConvergenceTest(this%objs(i)%solver,relRedKSPConverged,&
                                               this%objs(i)%relRedCtxAddr,PETSC_NULL_FUNCTION,ierr)
                    CHKERRQ(ierr)
                    this%objs(i)%allowRelRedConv = .true.
                end if
            end if
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
    !! When optional rtol > solverToleranceRel is passed for a compatible iterative KSP (norm-based KSPConvergedDefault),
    !! enables relative-reduction early stopping at rtol * ||r_0|| while delegating exact tolerance and divergence checks to
    !! PETSc KSPConvergedDefault, and reports via usedExactTol whether KSPConvergedDefault itself accepted the solve.

    class(PETScController)      ,intent(inout)  :: this
    real(rk)      ,dimension(:) ,intent(in)     :: knownVec 
    real(rk)      ,dimension(:) ,intent(inout)  :: unknownVec
    integer(ik)   ,optional     ,intent(in)     :: objGroup
    logical       ,optional     ,intent(in)     :: reusePC !! Default false; overrides -ksp_reuse_preconditioner
    real(rk)      ,optional     ,intent(in)     :: rtol !! Optional relative residual reduction vs ||r_0||
    logical       ,optional     ,intent(out)    :: usedExactTol !! True if the configured KSP policy or exact floor accepted the solve

    PetscErrorCode :: ierr
    integer(ik)    :: usedGroup ,attempt
    logical     :: useRelRed ,exactTol
    PetscBool :: reuseThisSolve

    PetscScalar ,dimension(:) ,pointer :: vecPointer ,solutionPointer

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to solve linear system using undefined PETSc controller")
        call assert(this%objCreated,"Attempted to solve linear system without first creating PETSc objects")
    end if

    usedGroup = 1
    if (present(objGroup)) usedGroup = objGroup

    ! Only enable relative-reduction early stopping when the KSP uses wrapped KSPConvergedDefault
    ! and the requested rtol is strictly looser than solverToleranceRel.
    useRelRed = .false.
    if (present(rtol)) useRelRed = this%objs(usedGroup)%allowRelRedConv .and. (rtol > this%options%solverToleranceRel)

    reuseThisSolve = PETSC_FALSE
    if (present(reusePC)) reuseThisSolve = reusePC
    ! Direct solvers (KSPPREONLY, e.g. LU/MUMPS) always solve directly in one application: never reuse a stale factorization.
    if (this%objs(usedGroup)%isPreonly) reuseThisSolve = PETSC_FALSE

    call MatSetValuesCOO(this%objs(usedGroup)%petscMat,this%objs(usedGroup)%cooValues,INSERT_VALUES,ierr)
    CHKERRQ(ierr)
    this%objs(usedGroup)%cooValues = 0.0d0

#if PETSC_VERSION_GE(3,23,0)
    call VecGetArray(this%objs(usedGroup)%rhs,vecPointer,ierr)
#else
    call VecGetArrayF90(this%objs(usedGroup)%rhs,vecPointer,ierr)
#endif
    CHKERRQ(ierr)
    vecPointer = knownVec
#if PETSC_VERSION_GE(3,23,0)
    call VecRestoreArray(this%objs(usedGroup)%rhs,vecPointer,ierr)
#else
    call VecRestoreArrayF90(this%objs(usedGroup)%rhs,vecPointer,ierr)
#endif
    CHKERRQ(ierr)

    call KSPSetTolerances(this%objs(usedGroup)%solver,this%options%solverToleranceRel,&
                                     this%options%solverToleranceAbs,&
                                     this%options%solverToleranceDiv,&
                                     this%options%maxSolverIters,ierr)
    CHKERRQ(ierr)

    if (this%objs(usedGroup)%allowRelRedConv) then
        if (useRelRed) then
            this%objs(usedGroup)%relRedCtx%relRedTarget = real(rtol,kind=kind(this%objs(usedGroup)%relRedCtx%relRedTarget))
        else
            this%objs(usedGroup)%relRedCtx%relRedTarget = 0.0d0
        end if
        this%objs(usedGroup)%relRedCtx%hitExact = .false.
    end if

    call KSPSetOperators(this%objs(usedGroup)%solver,this%objs(usedGroup)%petscMat,this%objs(usedGroup)%petscMat,ierr)
    CHKERRQ(ierr)

    do attempt = 1,2
        ! Set explicitly on every solve so reuse cannot leak to another caller of this group.
        call KSPSetReusePreconditioner(this%objs(usedGroup)%solver,reuseThisSolve,ierr)
        CHKERRQ(ierr)

        ! Keep the initial guess unchanged until the final result is extracted.
#if PETSC_VERSION_GE(3,23,0)
        call VecGetArray(this%objs(usedGroup)%sol,vecPointer,ierr)
#else
        call VecGetArrayF90(this%objs(usedGroup)%sol,vecPointer,ierr)
#endif
        CHKERRQ(ierr)
        vecPointer = unknownVec
#if PETSC_VERSION_GE(3,23,0)
        call VecRestoreArray(this%objs(usedGroup)%sol,vecPointer,ierr)
#else
        call VecRestoreArrayF90(this%objs(usedGroup)%sol,vecPointer,ierr)
#endif
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

    ! Exact tolerance is satisfied only when the solve converged and either inexact mode was not active or
    ! KSPConvergedDefault itself declared convergence (relRedCtx%hitExact == .true.).
    if (this%lastConvergedReason <= 0) then
        exactTol = .false.
    else if (useRelRed) then
        exactTol = this%objs(usedGroup)%relRedCtx%hitExact
    else
        exactTol = .true.
    end if
    if (present(usedExactTol)) usedExactTol = exactTol

#if PETSC_VERSION_GE(3,23,0)
    call VecGetArrayRead(this%objs(usedGroup)%sol,solutionPointer,ierr)
#else
    call VecGetArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
#endif
    CHKERRQ(ierr)

    unknownVec = real(solutionPointer,kind=rk)

#if PETSC_VERSION_GE(3,23,0)
    call VecRestoreArrayRead(this%objs(usedGroup)%sol,solutionPointer,ierr)
#else
    call VecRestoreArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
#endif
    CHKERRQ(ierr)

end subroutine linearSolve
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine relRedKSPConverged(ksp,n,rnorm,flag,ctxAddr,ierr)
    !! Thin wrapper around PETSc KSPConvergedDefault adding optional relative-reduction early stopping.
    !! Delegates exact ATOL/RTOL, norm type, PC side, zero-RHS, NaN/Inf, divergence, and maxits handling to
    !! KSPConvergedDefault first; only if still iterating at n > 0 with relRedTarget > 0 checks
    !! rnorm <= relRedTarget * ||r_0||.

    KSP                ,intent(in)  :: ksp     !! Active PETSc KSP solver handle
    PetscInt           ,intent(in)  :: n       !! Current Krylov iteration count (0 before first Krylov step)
    PetscReal          ,intent(in)  :: rnorm   !! Current residual norm reported by KSP under its active norm type
    KSPConvergedReason ,intent(out) :: flag    !! Convergence/divergence status returned to PETSc
    PetscFortranAddr   ,intent(in)  :: ctxAddr !! Opaque C address of the owning group RelRedConvContext
    PetscErrorCode     ,intent(out) :: ierr    !! PETSc error code

    type(c_ptr)                      :: cptr
    type(RelRedConvContext) ,pointer :: ctx

    cptr = transfer(ctxAddr,cptr)
    call c_f_pointer(cptr,ctx)

    if (n == 0) then
        ctx%hitExact = .false.
        ctx%initRnorm = rnorm
    end if

    ierr = kspConvergedDefaultC(transfer(ksp%v,c_null_ptr),n,rnorm,flag,transfer(ctx%defaultCtxAddr,c_null_ptr))
    CHKERRQ(ierr)

    if (flag /= KSP_CONVERGED_ITERATING) then
        ctx%hitExact = (flag > 0)
        return
    end if

    if (n > 0 .and. ctx%relRedTarget > 0.0d0 .and. rnorm <= ctx%relRedTarget * ctx%initRnorm) then
        flag = KSP_CONVERGED_RTOL
        ctx%hitExact = .false.
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

    PetscErrorCode :: ierr
    integer(ik)    :: i 

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
            if (associated(this%objs(i)%relRedCtx)) then
                if (this%objs(i)%relRedCtx%defaultCtxAddr /= 0) then
                    ierr = kspConvergedDefaultDestroyC(transfer(this%objs(i)%relRedCtx%defaultCtxAddr,c_null_ptr))
                    CHKERRQ(ierr)
                    this%objs(i)%relRedCtx%defaultCtxAddr = 0
                end if
                deallocate(this%objs(i)%relRedCtx)
                nullify(this%objs(i)%relRedCtx)
            end if
            this%objs(i)%relRedCtxAddr = 0
        end do
        call PetscFinalize(ierr)
        CHKERRQ(ierr)
    end if

end subroutine finalize
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule petsc_controller_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
