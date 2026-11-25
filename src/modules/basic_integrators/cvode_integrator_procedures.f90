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
submodule (cvode_integrator_class) cvode_integrator_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the CVODE integrator class

    implicit none

    type :: CVODECallbackContext
        class(CVODEIntegrator)  ,pointer               :: integrator => null()
        type(Modeller)          ,pointer               :: modeller => null()
        type(IntArray)          ,pointer ,dimension(:) :: termGroups => null()
        integer(ik)             ,pointer ,dimension(:) :: modelIndices => null()
        type(IntArray)          ,pointer ,dimension(:) :: varIndices => null()
        type(LogicalArray)      ,pointer ,dimension(:) :: updateRules => null()
        logical                 ,pointer ,dimension(:) :: modelDataUpdateRules => null()
        type(CommunicationData) ,pointer               :: commData => null()
        integer(ik)                                    :: timeVarIndex = 0
        logical                                        :: commNeeded = .false.
        logical                                        :: nonTrivialUpdate = .false.
        logical                                        :: nonTrivialModelDataUpdate = .false.
    end type CVODECallbackContext

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine checkCVODEError(ierr, callName)
        !! Abort with diagnostic message if a SUNDIALS/CVODE routine returned a non-zero error code

        integer(c_int) ,intent(in) :: ierr
        character(*)   ,intent(in) :: callName
        character(len=32)          :: codeStr

        if (ierr /= 0) then
            write(codeStr,'(I0)') ierr
            call assert(.false., "SUNDIALS/CVODE call "//callName//" failed with error code "//trim(codeStr))
        end if

    end subroutine checkCVODEError
!-----------------------------------------------------------------------------------------------------------------------------------
    integer(c_int) function cvodeRhsCallback(t, sunvecY, sunvecYdot, uData) result(retval) bind(C)
        !! Module-level bind(C) RHS callback for CVODE receiving Fortran state via FCVodeSetUserData

        real(c_double) ,value :: t
        type(N_Vector)        :: sunvecY
        type(N_Vector)        :: sunvecYdot
        type(c_ptr)    ,value :: uData

        type(CVODECallbackContext) ,pointer :: ctx
        integer(ik)                         :: j ,k

        call c_f_pointer(uData, ctx)

        ctx%integrator%rhsVec => FN_VGetArrayPointer(sunvecYdot)
        ctx%integrator%yVec => FN_VGetArrayPointer(sunvecY)

        ctx%integrator%copyBufferVals = ctx%integrator%yVec

        call ctx%integrator%bufferY%copyNamedVarsFromVec(ctx%integrator%copyBufferVals, ctx%integrator%evolvedVars)

        if (ctx%integrator%bufferY%isVarNameRegistered("time")) then
            ctx%integrator%bufferY%variables(ctx%timeVarIndex)%entry(1) = t
        end if

        call ctx%modeller%callManipulator(0, ctx%integrator%bufferY, ctx%integrator%bufferY)
        if (ctx%commNeeded) then
            call ctx%modeller%safeCommAndDeriv(ctx%commData, ctx%integrator%bufferY, derivPriority=0)
        else
            call ctx%integrator%bufferY%calculateDerivedVars(derivPriority=0)
        end if

        call ctx%integrator%bufferRHS%zeroVars(ctx%integrator%evolvedVars)

        do j = 1, size(ctx%modelIndices)
            if (ctx%nonTrivialModelDataUpdate) then
                if (ctx%modelDataUpdateRules(j)) &
                    call ctx%modeller%updateModelData(ctx%modelIndices(j), ctx%integrator%bufferY, updatePriority=0)
            end if
            do k = 1, size(ctx%termGroups(j)%entry)
                if (ctx%nonTrivialUpdate) then
                    if (ctx%updateRules(j)%entry(k)) &
                        call ctx%modeller%updateModelTermGroup(ctx%modelIndices(j), ctx%termGroups(j)%entry(k), &
                                                               ctx%integrator%bufferY)
                end if
                call ctx%modeller%calculateMatGroupValsInModel(ctx%modelIndices(j), &
                                                               ctx%termGroups(j)%entry(k), &
                                                               ctx%integrator%bufferY)
                ctx%integrator%bufferRHS%variables(ctx%varIndices(j)%entry(k))%entry = &
                    ctx%integrator%bufferRHS%variables(ctx%varIndices(j)%entry(k))%entry + &
                    ctx%modeller%evaluateModelTermGroup(ctx%modelIndices(j), ctx%termGroups(j)%entry(k), &
                                                        ctx%integrator%bufferY)
            end do
        end do

        call ctx%integrator%bufferRHS%copyNamedVarsToVec(ctx%integrator%copyBufferVals, ctx%integrator%evolvedVars)

        ctx%integrator%rhsVec = ctx%integrator%copyBufferVals

        retval = 0

    end function cvodeRhsCallback
!-----------------------------------------------------------------------------------------------------------------------------------
    integer(c_int) function cvodeRhsGFunCallback(nnlocal, t, sunvecY, sunvecG, uData) result(retval) bind(C)
        !! Module-level bind(C) local RHS callback for CVODE BBD preconditioner

        integer(c_int64_t) ,value :: nnlocal
        real(c_double)     ,value :: t
        type(N_Vector)            :: sunvecY
        type(N_Vector)            :: sunvecG
        type(c_ptr)        ,value :: uData

        if (nnlocal < 0) then
            retval = -1
            return
        end if

        retval = cvodeRhsCallback(t, sunvecY, sunvecG, uData)

    end function cvodeRhsGFunCallback
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrateCVODE(this,manipulatedModeller,outputVars,inputVars)
        !! Implementation of abstract manipulate routine for the SUNDIALS CVODE integrator

        class(CVODEIntegrator)                ,intent(inout) :: this
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

        if (assertions) then
            call assert(this%isDefined(),"Integration routine called by undefined CVODE integrator")
            call assert(manipulatedModeller%isDefined(),"Undefined Modeller object passed to CVODE integration routine")
            call assert(outputVars%isDefined(),"outputVars passed to integrateCVODE not defined")
            call assert(inputVars%isDefined(),"inputVars passed to integrateCVODE not defined")
        end if

        select type (manipulatedModeller)
        type is (Modeller)
            call runCVODEIntegration(this, manipulatedModeller, outputVars, inputVars)
        class default
            error stop "Unsupported surrogate passed to CVODE integrator"
        end select

    end subroutine integrateCVODE
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine runCVODEIntegration(this,manipulatedModeller,outputVars,inputVars)
        !! Execute CVODE integration with target-qualified arguments so callback context is standard-conforming

        class(CVODEIntegrator)   ,intent(inout) ,target :: this
        type(Modeller)           ,intent(inout) ,target :: manipulatedModeller
        class(VariableContainer) ,intent(inout)         :: outputVars
        class(VariableContainer) ,intent(in)            :: inputVars

        type(IntArray) ,allocatable ,dimension(:) ,target :: termGroups
        integer(ik)    ,allocatable ,dimension(:) ,target :: modelIndices

        integer(ik)                         :: numSteps ,i ,j ,k ,mpiErr
        real(rk) ,allocatable ,dimension(:) :: dt
        real(rk)                            :: fullTimestep ,controllerTimestep ,currentTimeVal

        type(IntArray)     ,allocatable ,dimension(:) ,target :: varIndices
        type(LogicalArray) ,allocatable ,dimension(:) ,target :: updateRules
        integer(ik)                                           :: timeVarIndex
        logical            ,allocatable ,dimension(:) ,target :: modelDataUpdateRules
        integer(ik)        ,allocatable ,dimension(:)         :: indices

        logical                           :: commNeeded
        logical                           :: nonTrivialUpdate
        logical                           :: nonTrivialModelDataUpdate
        logical                           :: timeEvolving
        type(CommunicationData) ,target   :: commData
        type(CVODECallbackContext) ,target :: cbCtx

        integer(c_int)   :: ierr
        integer(c_int)   :: maxres ,lastOrder(1)
        integer(c_int64_t) :: nlocal ,neq ,mu ,ml ,mudq ,mldq ,numCVODESteps(1)
        integer(c_int)   :: iPretype0 !! Preconditioner type (left)
        integer(c_int)   :: iGStype   !! Gram-Schmidt orthogonalization type (classical)
        real(c_double)   :: t(1) ,rtol ,atol

        iPretype0 = 1
        iGStype = 1
        timeVarIndex = 0

        termGroups = this%getTermGroups()
        modelIndices = this%getModelIndices()
        nonTrivialUpdate = this%hasNonTrivialUpdate()
        if (nonTrivialUpdate) updateRules = this%getUpdateRules()
        nonTrivialModelDataUpdate = this%hasNonTrivialModelDataUpdate()
        if (nonTrivialModelDataUpdate) modelDataUpdateRules = this%getModelDataUpdateRules()

        fullTimestep = this%getTimestep()
        timeEvolving = this%isTimeEvolving()
        ! Calculate the controller timestep, and if it is smaller than the full timestep divide the integration into numSteps > 1
        if (this%hasTimestepController()) then
            controllerTimestep = this%getTimestepFromController(inputVars)
            numSteps = ceiling(fullTimestep/controllerTimestep)
            allocate(dt(numSteps))
            if (numSteps > 1) then
                dt(1:numSteps-1) = controllerTimestep
                dt(numSteps) = fullTimestep - sum(dt(1:numSteps-1))
            else
                dt = fullTimestep
            end if
        else
            numSteps = 1
            allocate(dt(1))
            dt = fullTimestep
        end if

        ! Check if communication needed
        commNeeded = this%isCommunicationNeeded()
        if (commNeeded) commData = this%getCommunicationData()

        currentTimeVal = 0.0d0
        if (inputVars%isVarNameRegistered("time")) then
            timeVarIndex = inputVars%getVarIndex("time")
            currentTimeVal = inputVars%variables(timeVarIndex)%entry(1)
        end if

        allocate(varIndices(size(modelIndices)))
        do i = 1, size(modelIndices)
            allocate(varIndices(i)%entry(size(termGroups(i)%entry)))
            do j = 1, size(termGroups(i)%entry)
                call assert(.not. manipulatedModeller%isModelTermGroupMixed(termGroups(i)%entry(j),modelIndices(i)),&
                    "CVODE integrator does not support mixed term groups")
                varIndices(i)%entry(j) = inputVars%getVarIndex(manipulatedModeller&
                                                              %getEvolvedVarInTermGroup(termGroups(i)%entry(j),modelIndices(i)))
                call assert(.not.inputVars%isStationary(manipulatedModeller&
                %getEvolvedVarInTermGroup(termGroups(i)%entry(j),modelIndices(i))),&
                    "CVODE integrator detected stationary variable among evolved variables - this is unsupported")
            end do
        end do

        mu = int(this%options%bbdmukeep, kind=c_int64_t)
        ml = int(this%options%bbdmlkeep, kind=c_int64_t)
        mudq = int(this%options%bbdmudq, kind=c_int64_t)
        mldq = int(this%options%bbdmldq, kind=c_int64_t)

        if (.not. allocated(this%evolvedVars)) then
            allocate(indices(0))
            do i = 1, size(modelIndices)
                do j = 1, size(termGroups(i)%entry)
                    indices = [indices, varIndices(i)%entry(j)]
                end do
            end do

            indices = removeDupeInts(indices)
            allocate(this%evolvedVars(size(indices)))
            do i = 1, size(indices)
                this%evolvedVars(i)%string = inputVars%getVarName(indices(i))
            end do

            allocate(this%bufferRHS, source=inputVars)
            allocate(this%bufferY, source=inputVars)
            this%bufferY%variables = inputVars%variables
            this%bufferRHS%variables = inputVars%variables

            call this%bufferY%copyNamedVarsToVec(this%copyBufferVals, this%evolvedVars)

            nlocal = int(size(this%copyBufferVals), kind=c_int64_t)
            call MPI_Allreduce(nlocal, neq, 1, MPI_INTEGER8, MPI_SUM, this%mpiComm, mpiErr)
            call assert(mpiErr == 0, "MPI_Allreduce failed while computing global CVODE equation count")

            ! Lazily initialize SUNDIALS context on first integration so shallow copies before integration own no handles
            ierr = FSUNContext_Create(this%mpiComm%MPI_VAL, this%sunctx)
            call checkCVODEError(ierr, "FSUNContext_Create")

            this%sunVecY => FN_VNew_Parallel(this%mpiComm%MPI_VAL, nlocal, neq, this%sunctx)
            this%sunVecYDot => FN_VNew_Parallel(this%mpiComm%MPI_VAL, nlocal, neq, this%sunctx)
            call assert(associated(this%sunVecY) .and. associated(this%sunVecYDot), &
                        "FN_VNew_Parallel returned null N_Vector")

            this%rhsVec(1:nlocal) => FN_VGetArrayPointer(this%sunVecYDot)
            this%yVec(1:nlocal) => FN_VGetArrayPointer(this%sunVecY)

            if (this%options%amMethod) then
                this%cvode = FCVodeCreate(CV_ADAMS, this%sunctx)
            else
                this%cvode = FCVodeCreate(CV_BDF, this%sunctx)
            end if
            call assert(c_associated(this%cvode), "FCVodeCreate returned null CVODE handle")

            ierr = FCVodeInit(this%cvode, c_funloc(cvodeRhsCallback), currentTimeVal, this%sunVecY)
            call checkCVODEError(ierr, "FCVodeInit")

            ierr = FCVodeSetMaxOrd(this%cvode, this%options%maxOrder)
            call checkCVODEError(ierr, "FCVodeSetMaxOrd")

            if (this%options%stabLimitDet) then
                ierr = FCVodeSetStabLimDet(this%cvode, SUNTRUE)
                call checkCVODEError(ierr, "FCVodeSetStabLimDet")
            end if

            ierr = FCVodeSetMaxNumSteps(this%cvode, int(this%options%maxInternalSteps, kind=c_long))
            call checkCVODEError(ierr, "FCVodeSetMaxNumSteps")

            if (this%options%minTimestep > 0.0d0) then
                ierr = FCVodeSetMinStep(this%cvode, this%options%minTimestep)
                call checkCVODEError(ierr, "FCVodeSetMinStep")
            end if
            if (this%options%maxTimestep > 0.0d0) then
                ierr = FCVodeSetMaxStep(this%cvode, this%options%maxTimestep)
                call checkCVODEError(ierr, "FCVodeSetMaxStep")
            end if
            if (this%options%startingTimestep > 0.0d0) then
                ierr = FCVodeSetInitStep(this%cvode, this%options%startingTimestep)
                call checkCVODEError(ierr, "FCVodeSetInitStep")
            end if

            this%sunls => FSUNLinSol_SPGMR(this%sunVecY, iPretype0, 0, this%sunctx)
            call assert(associated(this%sunls), "FSUNLinSol_SPGMR returned null linear solver")
            this%sunmat => null()

            ierr = FCVodeSetLinearSolver(this%cvode, this%sunls, this%sunmat)
            call checkCVODEError(ierr, "FCVodeSetLinearSolver")

            rtol = this%options%reltol
            atol = this%options%abstol

            ierr = FSUNLinSol_SPGMRSetGSType(this%sunls, iGStype)
            call checkCVODEError(ierr, "FSUNLinSol_SPGMRSetGSType")

            maxres = this%options%maxRestarts

            ierr = FSUNLinSol_SPGMRSetMaxRestarts(this%sunls, maxres)
            call checkCVODEError(ierr, "FSUNLinSol_SPGMRSetMaxRestarts")

            ierr = FCVodeSStolerances(this%cvode, rtol, atol)
            call checkCVODEError(ierr, "FCVodeSStolerances")

            ierr = FCVBBDPrecInit(this%cvode, nlocal, mudq, mldq, mu, ml, 0.0d0, &
                                  c_funloc(cvodeRhsGFunCallback), c_null_ptr)
            call checkCVODEError(ierr, "FCVBBDPrecInit")
        else
            this%bufferY%variables = inputVars%variables
        end if

        ! Bind callback context for this integration step
        cbCtx%integrator => this
        cbCtx%modeller => manipulatedModeller
        cbCtx%termGroups => termGroups
        cbCtx%modelIndices => modelIndices
        cbCtx%varIndices => varIndices
        if (nonTrivialUpdate) cbCtx%updateRules => updateRules
        if (nonTrivialModelDataUpdate) cbCtx%modelDataUpdateRules => modelDataUpdateRules
        if (commNeeded) cbCtx%commData => commData
        cbCtx%timeVarIndex = timeVarIndex
        cbCtx%commNeeded = commNeeded
        cbCtx%nonTrivialUpdate = nonTrivialUpdate
        cbCtx%nonTrivialModelDataUpdate = nonTrivialModelDataUpdate

        ierr = FCVodeSetUserData(this%cvode, c_loc(cbCtx))
        call checkCVODEError(ierr, "FCVodeSetUserData")

        ! Perform start-of-timestep updates
        do j = 1, size(modelIndices)
            if (nonTrivialModelDataUpdate) then
                if (.not. modelDataUpdateRules(j)) &
                    call manipulatedModeller%updateModelData(modelIndices(j), this%bufferY, updatePriority=0)
            end if
            do k = 1, size(termGroups(j)%entry)
                if (nonTrivialUpdate) then
                    if (.not. updateRules(j)%entry(k)) &
                        call manipulatedModeller%updateModelTermGroup(modelIndices(j), termGroups(j)%entry(k), this%bufferY)
                end if
            end do
        end do

        call manipulatedModeller%callManipulator(2, this%bufferY, this%bufferY)
        outputVars%variables = this%bufferY%variables
        if (commNeeded) then
            call manipulatedModeller%safeCommAndDeriv(commData, outputVars)
        else
            ! A little bit overkill - could specialize to update only those variables that could have changed
            call outputVars%calculateDerivedVars()
        end if

        call this%bufferY%copyNamedVarsToVec(this%copyBufferVals, this%evolvedVars)
        this%yVec => FN_VGetArrayPointer(this%sunvecY)
        this%yVec = this%copyBufferVals
        ierr = FCVodeReInit(this%cvode, currentTimeVal, this%sunVecY)
        call checkCVODEError(ierr, "FCVodeReInit")
        ierr = FCVBBDPrecReInit(this%cvode, mudq, mldq, 0.0d0)
        call checkCVODEError(ierr, "FCVBBDPrecReInit")

        t(1) = currentTimeVal
        do i = 1, numSteps
            ierr = FCVode(this%cvode, t(1)+dt(i), this%sunVecY, t, CV_NORMAL)
            call checkCVODEError(ierr, "FCVode")

            ierr = FCVodeGetNumSteps(this%cvode, numCVODESteps)
            call checkCVODEError(ierr, "FCVodeGetNumSteps")
            call printNamedValue(this%integratorName//": CVODE number of steps", int(numCVODESteps(1), kind=ik))
            ierr = FCVodeGetLastOrder(this%cvode, lastOrder)
            call checkCVODEError(ierr, "FCVodeGetLastOrder")
            call printNamedValue(this%integratorName//": Last internal step order", int(lastOrder(1), kind=ik))
            if (numSteps > 1) then
                this%yVec => FN_VGetArrayPointer(this%sunvecY)
                this%copyBufferVals = this%yVec

                call printMessage(this%integratorName//": Consolidating internal CVODE steps")
                call this%bufferY%copyNamedVarsFromVec(this%copyBufferVals, this%evolvedVars)

                if (this%bufferY%isVarNameRegistered("time")) &
                    this%bufferY%variables(timeVarIndex)%entry(1) = t(1)
                call manipulatedModeller%callManipulator(1, this%bufferY, this%bufferY)
                if (commNeeded) then
                    call manipulatedModeller%safeCommAndDeriv(commData, this%bufferY)
                else
                    call this%bufferY%calculateDerivedVars()
                end if
                call this%bufferY%copyNamedVarsToVec(this%copyBufferVals, this%evolvedVars)

                this%yVec = this%copyBufferVals
            end if
        end do

        ! Detach stack-allocated callback context before returning
        ierr = FCVodeSetUserData(this%cvode, c_null_ptr)
        call checkCVODEError(ierr, "FCVodeSetUserData(null)")

        this%yVec => FN_VGetArrayPointer(this%sunvecY)
        this%copyBufferVals = this%yVec

        call this%bufferY%copyNamedVarsFromVec(this%copyBufferVals, this%evolvedVars)

        if (this%bufferY%isVarNameRegistered("time")) then
            this%bufferY%variables(timeVarIndex)%entry(1) = currentTimeVal
            if (timeEvolving) this%bufferY%variables(timeVarIndex)%entry(1) = t(1)
        end if
        outputVars%variables = this%bufferY%variables
        call manipulatedModeller%callManipulator(2, outputVars, outputVars)
        if (commNeeded) then
            call manipulatedModeller%safeCommAndDeriv(commData, outputVars)
        else
            call outputVars%calculateDerivedVars()
        end if

    end subroutine runCVODEIntegration
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCVODEIntegrator(this,mpiCont,options,modelList,termGroups,&
            evolvesTimeVar,dtController,initialTimestep,integratorName)

        class(CVODEIntegrator)                    ,intent(inout) :: this
        type(CVODEOptions)                        ,intent(in)    :: options !! Options object containing CVODE options
        type(MPIController)                       ,intent(in)    :: mpiCont !! MPI Controller used to initialise the sundials context
        integer(ik)    ,optional    ,dimension(:) ,intent(in)    :: modelList !! List of models this integrator will be responsible for
        type(IntArray) ,optional    ,dimension(:) ,intent(in)    :: termGroups !! Term groups this integrator is responsible for - should conform with modelList
        logical        ,optional                  ,intent(in)    :: evolvesTimeVar !! Set to true if this integrator is allowed to change the "time" varible (if present in passed variable container)
        class(TimestepController) ,optional       ,intent(in)    :: dtController !! User-supplied timestep controller object
        real(rk)                  ,optional       ,intent(in)    :: initialTimestep !! Default timestep
        character(*)                    ,optional ,intent(in)    :: integratorName !! Name of integrator used in printing

        this%mpiComm = mpiCont%getComm()
        this%numProcs = mpiCont%getWorldSize()

        if (assertions .or. assertionLvl >= 0) then

            if (present(termGroups)) call assert(present(modelList),&
                "Term groups object passed to CVODE integrator constructor without model list")
            if (present(modelList) .and. present(termGroups)) call assert(size(modelList)==size(termGroups),"If both model list and&
            & term groups arguments are passed to initCVODEIntegrator both must have the same length")
            if (present(dtController)) call assert(dtController%isDefined(),"Undefined timestep controller passed to CVODE&
            & integrator constructor")
        end if

        if (present(evolvesTimeVar)) call this%setTimeEvolving(evolvesTimeVar)
        if (present(modelList)) call this%setModelIndices(modelList)
        if (present(termGroups)) call this%setTermGroups(termGroups)
        if (present(initialTimestep)) call this%setTimestep(initialTimestep)
        if (present(dtController)) call this%setTimestepController(dtController)
        call this%setNonTrivialModelDataUpdate(.false.)
        call this%setNonTrivialUpdate(.false.)

        this%integratorName = "IntegratorCVODE"
        if (present(integratorName)) this%integratorName = integratorName
        this%options = options
        call this%makeDefined()

    end subroutine initCVODEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    impure elemental module subroutine finalizeCVODEIntegrator(this)
        !! Safely release owned SUNDIALS/CVODE resources in reverse creation order

        type(CVODEIntegrator) ,intent(inout) :: this

        integer(c_int) :: ierr

        if (c_associated(this%cvode)) then
            call FCVodeFree(this%cvode)
            this%cvode = c_null_ptr
        end if
        if (associated(this%sunls)) then
            ierr = FSUNLinSolFree(this%sunls)
            nullify(this%sunls)
        end if
        if (associated(this%sunVecYDot)) then
            call FN_VDestroy(this%sunVecYDot)
            nullify(this%sunVecYDot)
        end if
        if (associated(this%sunVecY)) then
            call FN_VDestroy(this%sunVecY)
            nullify(this%sunVecY)
        end if
        nullify(this%sunmat)
        nullify(this%rhsVec)
        nullify(this%yVec)
        if (c_associated(this%sunctx)) then
            ierr = FSUNContext_Free(this%sunctx)
            this%sunctx = c_null_ptr
        end if

    end subroutine finalizeCVODEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule cvode_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
