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
    type(SparseRowData)         ,intent(in)     :: rowData 

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

    integer(ik) :: i ,ierr ,usedGroup

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to add row values to matrix of undefined PETSc controller object")
        call assert(this%objCreated,&
        "Attempted to add row values to matrix of PETSc controller before PETSc objects were generated")
    end if

    mult = real(1.0d0,kind=rk)
    usedGroup = 1 
    if (present(objGroup)) usedGroup = objGroup
    if (present(multConst)) mult = multConst 

    do i = 1,size(rowData%rowIndex)
        call MatSetValues(this%objs(usedGroup)%petscMat,1, [rowData%rowIndex(i) - 1], &
                                        size(rowData%columnVector(i)%entry),rowData%columnVector(i)%entry - 1,&
                                        mult * rowData%values(i)%entry,ADD_VALUES,ierr)
        CHKERRQ(ierr)
    end do  

end subroutine addRowValuesToMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine createPETScObjs(this) 
    !! Assemble the preallocation object, create PETSc objects, and deallocate the preallocation object

    class(PETScController)      ,intent(inout)  :: this

    integer(ik) :: ierr ,i

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to create PETSc objects in undefined PETSc controller")
        call assert(.not. this%objCreated,"Attempted to create PETSc objects when they are already created")
    end if

    call this%preallocData%assembleData()

    do i = 1, size(this%objs)

        call PETScOptionsInsertString(PETSC_NULL_OPTIONS,this%options%petscOptions,ierr)
        CHKERRQ(ierr)

        call MatCreateAIJ(PETSc_Comm_World,this%preallocData%getLocalNumDoFs(),this%preallocData%getLocalNumDoFs(),&
                                        this%preallocData%getTotalNumDoFs(),this%preallocData%getTotalNumDoFs(),&
                                        1,this%preallocData%getNumNonzerosDiag(),&
                                        1,this%preallocData%getNumNonzerosOffDiag(),&
                                        this%objs(i)%petscMat,ierr)
        CHKERRQ(ierr) 

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
module subroutine linearSolve(this,knownVec,unknownVec,objGroup) 
    !! Solve the equation petscMat * unknownVec = knownVec where the vectors are local to each processor - assembles and zeros out the
    !! matrix before/after the solve

    class(PETScController)      ,intent(inout)  :: this
    real(rk)      ,dimension(:) ,intent(in)     :: knownVec 
    real(rk)      ,dimension(:) ,intent(out)    :: unknownVec
    integer(ik)   ,optional     ,intent(in)     :: objGroup

    integer(ik) :: i ,ierr ,offset ,usedGroup

    PetscScalar ,dimension(:) ,pointer :: solutionPointer

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to solve linear system using undefined PETSc controller")
        call assert(this%objCreated,"Attempted to solve linear system without first creating PETSc objects")
    end if

    usedGroup = 1
    if (present(objGroup)) usedGroup = objGroup

    call MatAssemblyBegin(this%objs(usedGroup)%petscMat,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRQ(ierr)
    call MatAssemblyEnd(this%objs(usedGroup)%petscMat,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRQ(ierr)
    

    offset = this%preallocData%getLocalDoFOffset()
    call VecSetValues(this%objs(usedGroup)%rhs,this%preallocData%getLocalNumDoFs(),&
                     [(i + offset - 1,i=1,this%preallocData%getLocalNumDoFs())],&
                              knownVec,INSERT_VALUES,ierr)
    CHKERRQ(ierr)

    call VecAssemblyBegin(this%objs(usedGroup)%rhs,ierr)
    CHKERRQ(ierr)
    
    call VecAssemblyEnd(this%objs(usedGroup)%rhs,ierr)
    CHKERRQ(ierr)

    call VecSetValues(this%objs(usedGroup)%sol,this%preallocData%getLocalNumDoFs(),&
                     [(i + offset - 1,i=1,this%preallocData%getLocalNumDoFs())],&
                              knownVec,INSERT_VALUES,ierr)
    CHKERRQ(ierr)

    call VecAssemblyBegin(this%objs(usedGroup)%sol,ierr)
    CHKERRQ(ierr)
    
    call VecAssemblyEnd(this%objs(usedGroup)%sol,ierr)
    CHKERRQ(ierr)

    call KSPSetTolerances(this%objs(usedGroup)%solver,this%options%solverToleranceRel,&
                                     this%options%solverToleranceAbs,&
                                     this%options%solverToleranceDiv,&
                                     this%options%maxSolverIters,ierr)
    CHKERRQ(ierr)

    call KSPSetInitialGuessNonzero(this%objs(usedGroup)%solver,PETSC_TRUE,ierr)
    CHKERRQ(ierr)

    call KSPSetOperators(this%objs(usedGroup)%solver,this%objs(usedGroup)%petscMat,this%objs(usedGroup)%petscMat,ierr)
    CHKERRQ(ierr)

    call KSPSolve(this%objs(usedGroup)%solver,this%objs(usedGroup)%rhs,this%objs(usedGroup)%sol,ierr)
    CHKERRQ(ierr)

    call KSPGetConvergedReason(this%objs(usedGroup)%solver,this%lastConvergedReason,ierr)
    CHKERRQ(ierr)

    call KSPGetIterationNumber(this%objs(usedGroup)%solver,this%lastNumIterations,ierr)
    CHKERRQ(ierr)

    call VecGetArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
    CHKERRQ(ierr)

    unknownVec = real(solutionPointer,kind=rk)

    call VecRestoreArrayReadF90(this%objs(usedGroup)%sol,solutionPointer,ierr)
    CHKERRQ(ierr)

    call MatZeroEntries(this%objs(usedGroup)%petscMat,ierr) 
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
