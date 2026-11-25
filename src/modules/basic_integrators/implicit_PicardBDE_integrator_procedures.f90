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
submodule (implicit_PicardBDE_integrator_class) implicit_PicardBDE_integrator_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the implicit BDE integrator class

    use ,intrinsic :: ieee_arithmetic ,only: ieee_is_finite

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initBDEIntegrator(this,indexingObj,procRank,nonlinTol,absTol,maxIters,convergenceIndices,&
    modelList,termGroups,evolvesTimeVar,dtController,initialTimestep&
    ,use2Norm,petscGroup,intContOptions,integratorName,relaxationWeight,reusePreconditioner,andersonDepth)
    !! BDE integrator constructor 

    class(PicardBDEIntegrator)                ,intent(inout) :: this
    type(Indexing)                            ,intent(in)    :: indexingObj !! Indexing object to be used in initializing the implicit vectors 
    integer(ik)                               ,intent(in)    :: procRank !! Rank of current processors
    real(rk)       ,optional                  ,intent(in)    :: nonlinTol !! Picard iteration relative tolerance
    real(rk)       ,optional                  ,intent(in)    :: absTol !! Picard iteration absolute tolerance in epsilon units
    integer(ik)    ,optional                  ,intent(in)    :: maxIters !! Maximum number of Picard iterations
    integer(ik)    ,optional    ,dimension(:) ,intent(in)    :: convergenceIndices !! Variable indices in VariableContainer used to determine convergence of Picard iterations
    integer(ik)    ,optional    ,dimension(:) ,intent(in)    :: modelList !! List of models this integrator is responsible for
    type(IntArray) ,optional    ,dimension(:) ,intent(in)    :: termGroups !! Term groups this integrator is responsible for - should conform with modelList
    logical        ,optional                  ,intent(in)    :: evolvesTimeVar !! Set to true if this integrator is allowed to change the "time" varible (if present in passed variable container)
    class(TimestepController) ,optional       ,intent(in)    :: dtController !! User-supplied timestep controller object
    real(rk)                  ,optional       ,intent(in)    :: initialTimestep !! Default timestep
    logical        ,optional                  ,intent(in)    :: use2Norm !! Set to true if this integrator uses the 2-norm instead of the local inf-norm for checking individual variable convergence 
    integer(ik)    ,optional                  ,intent(in)    :: petscGroup !! PETSc obj group this solver should interact with (defaults to 1)
    type(InternalControllerOptions) ,optional ,intent(in)    :: intContOptions !! InternalControlOptions (if present turns on internal control)
    character(*)                    ,optional ,intent(in)    :: integratorName !! Name of integrator used in printing
    real(rk)                  ,optional       ,intent(in)    :: relaxationWeight !! relaxationWeight
    logical                   ,optional       ,intent(in)    :: reusePreconditioner !! Reuse within each attempted substep
    integer(ik)               ,optional       ,intent(in)    :: andersonDepth !! Type-II Anderson acceleration history depth m >= 0 (default 0 = disabled; allocates no AA workspace when 0)
    
    integer(ik) ,allocatable ,dimension(:) :: procDoFs
    integer(ik)                            :: nLoc

    if (assertions .or. assertionLvl >= 0) then
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to BDE integrator constructor")
        if (present(termGroups)) call assert(present(modelList),"Term groups object passed to BDE integrator constructor without&
                                                               & model list")
        if (present(modelList) .and. present(termGroups)) call assert(size(modelList)==size(termGroups),"If both model list and &
        & term groups arguments are passed to initBDEIntegrator both must have the same length")
        if (present(dtController)) call assert(dtController%isDefined(),"Undefined timestep controller passed to BDE&
        & integrator constructor")
        if (present(andersonDepth)) call assert(andersonDepth >= 0,"andersonDepth passed to initBDEIntegrator must be >= 0")
    end if

    procDoFs = indexingObj%getProcDoF()
    nLoc = procDoFs(procRank+1)

    allocate(this%implicitVectorNew(nLoc))
    allocate(this%implicitVectorOld(nLoc))
    allocate(this%implicitVectorInit(nLoc))

    this%implicitVectorNew = 0
    this%implicitVectorOld = 0
    this%implicitVectorInit = 0

    this%timesCalled = 0 
    this%totNumIters = 0

    !Default values
    this%nonlinTol = real(1.0d-12,kind=rk)
    this%maxIterations = 100
    this%absTol = real(1,kind=rk)
    this%reusePreconditioner = .false.
    if (present(reusePreconditioner)) this%reusePreconditioner = reusePreconditioner
    this%use2Norm = .false.
    if (present(use2Norm)) this%use2Norm = use2Norm
    this%associatedPETScObjGroup = 1
    if (present(petscGroup)) this%associatedPETScObjGroup = petscGroup

    if (present(nonlinTol)) this%nonlinTol = nonlinTol
    if (present(absTol)) this%absTol = absTol
    if (present(maxIters)) this%maxIterations = maxIters
    if (present(convergenceIndices)) this%convergenceTestVars = convergenceIndices

    if (present(evolvesTimeVar)) call this%setTimeEvolving(evolvesTimeVar) 
    if (present(modelList)) call this%setModelIndices(modelList)
    if (present(termGroups)) call this%setTermGroups(termGroups)
    if (present(initialTimestep)) call this%setTimestep(initialTimestep)
    if (present(dtController)) call this%setTimestepController(dtController)
    this%internalStepControl = .false.
    if (present(intContOptions)) then 
        this%internalControlOpts = intContOptions
        this%internalStepControl = .true.
    end if

    this%integratorName = "IntegratorBDE"
    if (present(integratorName)) this%integratorName = integratorName
    call this%setNonTrivialModelDataUpdate(.false.)
    call this%setNonTrivialUpdate(.false.)

    this%relaxationWeight = 1.0d0 
    if (present(relaxationWeight)) this%relaxationWeight = relaxationWeight

    this%andersonDepth = 0
    if (present(andersonDepth)) this%andersonDepth = andersonDepth
    if (this%andersonDepth > 0) then
        allocate(this%aaDeltaX(nLoc,this%andersonDepth))
        allocate(this%aaDeltaF(nLoc,this%andersonDepth))
        allocate(this%aaQ(nLoc,this%andersonDepth))
        allocate(this%aaPrevX(nLoc))
        allocate(this%aaPrevF(nLoc))
        allocate(this%aaWorkF(nLoc))
        this%aaDeltaX = 0.0d0
        this%aaDeltaF = 0.0d0
        this%aaQ = 0.0d0
        this%aaPrevX = 0.0d0
        this%aaPrevF = 0.0d0
        this%aaWorkF = 0.0d0
    end if

    call this%makeDefined()

end subroutine initBDEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine integrateBDE(this,manipulatedModeller,outputVars,inputVars) 
    !! Integration routine for BDE integrator - implementation of abstract manipulate routine

    class(PicardBDEIntegrator)            ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
    class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
    class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    integer(ik) :: timeVarIndex, nonlinIter
    real(rk)    :: dt, elapsedTime
    real(rk)    :: fullTimestep, controllerTimestep 

    logical                 :: solveSuccess
    logical                 :: timeEvolving
    logical                 :: firstStep

    real (rk) :: startingTime
    character(len = 140)                      :: tmpstring

    if (assertions) then 
        call assert(this%isDefined(),"Integration routine called by undefined BDE integrator")
        call assert(manipulatedModeller%isDefined(),"Undefined modeller object passed to BDE integration routine")
        call assert(outputVars%isDefined(),"outputVars passed to integrateBDE not defined")
        call assert(inputVars%isDefined(),"inputVars passed to integrateBDE not defined")
    end if

    fullTimestep = this%getTimestep()

    timeEvolving = this%isTimeEvolving()
    if (inputVars%isVarNameRegistered("time")) then 
        timeVarIndex = inputVars%getVarIndex("time")
        startingTime = inputVars%variables(timeVarIndex)%entry(1)
    end if

    elapsedTime = 0
    dt = fullTimestep
    this%internalControlOpts%restartCount = 0
    firstStep = .true.

    if (.not. this%internalControlOpts%lazyEval) then
        do 
            solveSuccess = .true.
            if (this%hasTimestepController() .and. .not. this%internalStepControl) then 
                if (firstStep) then 
                    controllerTimestep = this%getTimestepFromController(inputVars)
                else
                    controllerTimestep = this%getTimestepFromController(this%buffer)
                end if
                dt = min(dt,controllerTimestep)
            end if
            if (this%internalControlOpts%restartCount > this%internalControlOpts%maxRestarts) & 
                 error stop "Max BDE restarts reached"
            call printNamedValue(this%integratorName//": Attempting step with length",dt)

            if (firstStep) then 
                call tryIntegrate(this,manipulatedModeller,inputVars,dt,nonlinIter,solveSuccess)
            else
                call tryIntegrate(this,manipulatedModeller,this%buffer,dt,nonlinIter,solveSuccess)
            end if

            if (nonlinIter >=  this%maxIterations) then 
                call printMessage("WARNING: "//this%integratorName//" reached maximum number of iterations")

                if (this%internalStepControl) solveSuccess = .false.
            end if

            if (.not. solveSuccess) then  

                    call printMessage("Applying internal step control and restarting BDE integration")
                    dt = dt / this%internalControlOpts%stepMultiplier
                    this%internalControlOpts%restartCount = this%internalControlOpts%restartCount + 1
            end if

            if (solveSuccess) then
                elapsedTime = elapsedTime + dt
                dt = min(dt,fullTimestep - elapsedTime)
                this%totNumIters = this%totNumIters + nonlinIter
                this%timesCalled = this%timesCalled + 1
                this%internalControlOpts%restartCount = 0
                call printNamedValue(this%integratorName//": last number of iterations",nonlinIter)
                call printNamedValue(this%integratorName//": total number of iterations",this%totNumIters)
                call printNamedValue(this%integratorName//": average number of iterations"&
                                    ,real(this%totNumIters,kind=rk)/this%timesCalled)
                 
                write(tmpstring,'(A,ES11.4,A,ES11.4)') 'BDE substep successful. Elapsed time: ',elapsedTime, &
                                                        '. Requested timestep: ', fullTimestep
                call printMessage(trim(tmpstring))
                tmpstring=''

                if (firstStep .and. nonlinIter == 1 .and. dt > 0 .and. this%internalControlOpts%allowLazyEval) then 

                    this%internalControlOpts%lazyEval = .true.
                    exit 
                end if
                firstStep = .false.

            end if

            if (this%internalStepControl .and. solveSuccess) then
                if (nonlinIter < this%internalControlOpts%minNonlinIters) then
                    call printMessage(this%integratorName//&
                    ": Minimum nonlinear iterations reached - attempting to reduce number of BDE substeps")
                    dt = min(fullTimestep - elapsedTime,1.5*dt)
                end if
            end if


            if (solveSuccess .and. fullTimestep - elapsedTime < 100*epsilon(fullTimestep)) then
                exit
            end if 
        end do
   
        call copyRealArraysInPlace(outputVars%variables,this%buffer%variables)
        if (inputVars%isVarNameRegistered("time") .and. (.not. timeEvolving)) &
            outputVars%variables(timeVarIndex)%entry = startingTime

        select type (manipulatedModeller)
        type is (Modeller)
            call manipulatedModeller%callManipulator(2,outputVars,outputVars)
        class default
            error stop "Unsupported surrogate passed to BDE integrator"
        end select

    else 
                    call printMessage(this%integratorName//&
                    ": Lazy evaluation triggered, evolving only time")
        call copyRealArraysInPlace(outputVars%variables,inputVars%variables)
        if (inputVars%isVarNameRegistered("time") .and. timeEvolving) &
            outputVars%variables(timeVarIndex)%entry = startingTime + fullTimestep
    end if

end subroutine integrateBDE
!-----------------------------------------------------------------------------------------------------------------------------------
function checkConvergence(oldVars,newVars,indicesToCheck,nonlinTol,absTol,use2Norm,convergenceCounter) result(conv)
    !! Checks whether all variables determined by indicesToCheck have converged based on a given nonlinear tolerance

    type(RealArray) ,dimension(:) ,intent(in) :: oldVars !! Previous variable values
    type(RealArray) ,dimension(:) ,intent(in) :: newVars !! New variable values
    integer(ik)     ,dimension(:) ,intent(in) :: indicesToCheck !! Indices of oldVars/newVars to check for convergence
    real(rk)                      ,intent(in) :: nonlinTol !! Relative convergence tolerances for each variable 
    real(rk)                      ,intent(in) :: absTol !! Absolute tolerance in epsilon units for each variable
    logical                       ,intent(in) :: use2Norm !! True if this should use 2-norm instead of local inf-norm
    integer(ik)     ,dimension(:) ,intent(inout) :: convergenceCounter !! Incremented if the variable with a given index has not yet converged

    logical :: conv ,varConverged

    integer(ik) :: i ,haloDataChunkSize ,ind ,nonHaloLen
    real(rk) :: relError ,absError ,oldNorm ,absTolScaled

    conv = .true.
    absTolScaled = epsilon(absError)*absTol
    do i = 1, size(indicesToCheck)
        ind = indicesToCheck(i)
        haloDataChunkSize = 1 - lbound(newVars(ind)%entry,1)
        nonHaloLen = size(newVars(ind)%entry) - 2*haloDataChunkSize

        if (use2Norm) then

            absError = norm2(oldVars(ind)%entry(1:nonHaloLen)-newVars(ind)%entry(1:nonHaloLen))
            oldNorm = norm2(oldVars(ind)%entry(1:nonHaloLen))

            if (oldNorm > absTolScaled) then
                relError = absError / oldNorm
            else
                relError = absError
            end if
        else

            absError = maxval(abs(oldVars(ind)%entry(1:nonHaloLen)-newVars(ind)%entry(1:nonHaloLen)))

            if (all(abs(oldVars(ind)%entry(1:nonHaloLen)) > absTolScaled)) then 
                relError = maxval(abs(oldVars(ind)%entry(1:nonHaloLen)-newVars(ind)%entry(1:nonHaloLen))&
                    /abs(oldVars(ind)%entry(1:nonHaloLen)))
            else
                relError = absError
            end if
            
        end if
        varConverged = relError < nonlinTol .or. absError < absTolScaled
        if (.not. varConverged) then
            conv = .false.
            convergenceCounter(i) = convergenceCounter(i) + 1
        end if
    end do

end function
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNonlinTol(this,nonlinTol)
    !! Setter for nonlinTol value

    class(PicardBDEIntegrator)     ,intent(inout)  :: this
    real(rk)                       ,intent(in)     :: nonlinTol !! Picard iteration tolerance

    this%nonlinTol = nonlinTol

end subroutine setNonlinTol
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNonlinTol(this) result(nonlinTol)
    !! Getter for nonlinTol value

    class(picardBDEIntegrator) ,intent(in) :: this
    real(rk)                               :: nonlinTol

    if (assertions) call assertPure(this%isDefined(),"getNonlinTol called on undefined BDE integrator")

    nonlinTol = this%nonlinTol

end function getNonlinTol
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setMaxIterations(this,maxIter)
    !! Setter for maxIterations value

    class(picardBDEIntegrator)     ,intent(inout)  :: this
    integer(ik)                    ,intent(in)     :: maxIter

    this%maxIterations = maxIter

end subroutine setMaxIterations
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxIterations(this) result(maxIter)
    !! Geter for maxIterations value

    class(picardBDEIntegrator) ,intent(in) :: this
    integer(ik)                            :: maxIter

    if (assertions) call assertPure(this%isDefined(),"getMaxIterations called on undefined BDE integrator")

    maxIter = this%maxIterations

end function getMaxIterations
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setConvergenceIndices(this,convVars)
    !! Setter for convergenceTestVars array

    class(picardBDEIntegrator)        ,intent(inout)  :: this
    integer(ik) ,dimension(:)         ,intent(in)     :: convVars ! Variable indices in VariableContainer used to determine convergence of Picard iterations

    this%convergenceTestVars = convVars

end subroutine setConvergenceIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getConvergenceIndices(this) result(convVars)
    !! Getter for convergenceTestVars array

    class(picardBDEIntegrator)    ,intent(in) :: this
    integer(ik) ,dimension(:) ,allocatable    :: convVars

    if (assertions) call assertPure(this%isDefined(),"getConvergenceIndices called on undefined BDE integrator")

    convVars = this%convergenceTestVars

end function getConvergenceIndices
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine tryIntegrate(this,manipulatedModeller,inputVars,dt,nonlinIter,solveSuccess) 

    type(PicardBDEIntegrator)             ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
    class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    real(rk)                            ,intent(in)      :: dt
    logical                             ,intent(inout)   :: solveSuccess
    integer(ik)                         ,intent(out)     :: nonlinIter

    type(IntArray) ,allocatable ,dimension(:)  :: termGroups 
    integer(ik)    ,allocatable ,dimension(:)  :: modelIndices 

    type(LogicalArray) ,allocatable ,dimension(:) :: updateRules
    logical            ,allocatable ,dimension(:) :: modelDataUpdateRules
    integer(ik)                                   :: timeVarIndex ,j ,k ,convReason

    logical                 :: commNeeded
    logical                 :: nonTrivialUpdate
    logical                 :: nonTrivialModelDataUpdate
    logical                 :: nonTrivialConvergenceCheck
    logical                 :: tolReached

    type(CommunicationData) :: commData

    
    integer(ik)    ,allocatable ,dimension(:)  :: convergenceCounter
    integer(ik)    ,allocatable ,dimension(:)  :: aaOrder
    integer(ik)                                :: aaHistCount ,aaNextSlot
    logical                                    :: aaHasPrev ,aaApplied
    logical                                    :: forceExactTol !! When true, bypass inexact rtol and solve to exact solverToleranceRel
    logical                                    :: solveWasExact !! True iff the latest linear solve satisfied solverToleranceRel/Abs
    real(rk) ,parameter                        :: inexactKSPRelRed = real(1.0d-2,kind=rk) !! Per-Picard-iteration KSP residual reduction target

    termGroups = this%getTermGroups()
    modelIndices = this%getModelIndices()
    nonTrivialUpdate = this%hasNonTrivialUpdate()
    if (nonTrivialUpdate) updateRules = this%getUpdateRules()
    nonTrivialModelDataUpdate = this%hasNonTrivialModelDataUpdate()
    if (nonTrivialModelDataUpdate) modelDataUpdateRules = this%getModelDataUpdateRules()
    nonTrivialConvergenceCheck = allocated(this%convergenceTestVars) 

    if (nonTrivialConvergenceCheck) allocate(convergenceCounter(size(this%convergenceTestVars)))
    if (this%andersonDepth > 0) allocate(aaOrder(this%andersonDepth))

    !Check if communication needed
    commNeeded = this%isCommunicationNeeded()
    if (commNeeded) commData = this%getCommunicationData()

    select type (manipulatedModeller)
    type is (Modeller)
        if (.not. allocated(this%buffer)) then
            allocate(this%buffer,source=inputVars)
        else
            call copyRealArraysInPlace(this%buffer%variables,inputVars%variables)
        end if
        
            if (commNeeded) then 
                call manipulatedModeller%safeCommAndDeriv(commData,this%buffer,derivPriority=0)
            else
                !Calculate only highest priority derived variables in internal iterations
                call this%buffer%calculateDerivedVars(derivPriority=0)
            end if

        if (inputVars%isVarNameRegistered("time")) then 
            timeVarIndex = inputVars%getVarIndex("time")
            this%buffer%variables(timeVarIndex)%entry(1) = this%buffer%variables(timeVarIndex)%entry(1) + dt
        end if
        tolReached = .false.
        if (nonTrivialConvergenceCheck) convergenceCounter = 1
        ! Start every attempted substep/restart with an empty AA history so stale secants never cross timestep changes.
        aaHasPrev = .false.
        aaHistCount = 0
        aaNextSlot = 0
        aaApplied = .false.
        if (this%andersonDepth > 0) aaOrder = 0
        ! Stationary solves (dt == 0) and single-iteration linear steps always run at exact solverToleranceRel.
        forceExactTol = (dt == 0.0d0 .or. this%maxIterations <= 1)

        call this%buffer%copyImplicitVarsToVec(this%implicitVectorInit,ignoreStationary=.true.)
        nonlinIter = 0
        do while (nonlinIter < this%maxIterations)
            call this%buffer%copyImplicitVarsToVec(this%implicitVectorOld)
            this%implicitVectorNew = this%implicitVectorOld

            do j = 1,size(modelIndices)
                if (nonlinIter == 0) then 
                    call manipulatedModeller%updateModelData(modelIndices(j),this%buffer)
                else if (nonTrivialModelDataUpdate) then 
                    if (modelDataUpdateRules(j)) &
                    call manipulatedModeller%updateModelData(modelIndices(j),this%buffer,updatePriority=0)
                end if
                do k = 1,size(termGroups(j)%entry)
                    if (nonTrivialUpdate) then 
                        if (updateRules(j)%entry(k))&
                        call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                    else if (nonlinIter == 0) then 
                        call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                    end if
                    call manipulatedModeller%calculateMatGroupValsInModel(modelIndices(j),&
                                                                            termGroups(j)%entry(k),this%buffer)
                    call manipulatedModeller%addModelMatGroupToPETSc(modelIndices(j),termGroups(j)%entry(k),-dt,&
                                                                    this%associatedPETScObjGroup)
                end do
            end do

            ! Refresh PC on the first Picard solve (nonlinIter == 0) of every attempted substep/split, and never reuse for dt == 0.
            ! Intermediate Picard solves use relative reduction inexactKSPRelRed against ||r_0||; the final allowed iteration or
            ! any verification pass after forceExactTol is triggered runs at exact solverToleranceRel.
            if (forceExactTol .or. nonlinIter == this%maxIterations - 1) then
                call manipulatedModeller%linearSolvePETSc(this%implicitVectorInit,this%implicitVectorNew,.true.,&
                                                        convReason,this%associatedPETScObjGroup,&
                                                        reusePC=this%reusePreconditioner .and. nonlinIter > 0 .and. dt /= 0.0d0)
                solveWasExact = .true.
            else
                call manipulatedModeller%linearSolvePETSc(this%implicitVectorInit,this%implicitVectorNew,.true.,&
                                                        convReason,this%associatedPETScObjGroup,&
                                                        reusePC=this%reusePreconditioner .and. nonlinIter > 0 .and. dt /= 0.0d0,&
                                                        rtol=inexactKSPRelRed,usedExactTol=solveWasExact)
            end if

            ! Check if PETSc diverged 
            if (convReason < -1) then 
                call printNamedValue(this%integratorName//": PETScConvergedReason:",convReason)
                solveSuccess = .false.
                if (.not. this%internalStepControl) then
                    error stop "PETSc convergence reason negative and internal step control disabled"
                end if
                return
            end if
            if (nonTrivialConvergenceCheck) then
                if (.not. allocated(this%oldBufferVals)) then
                    allocate(this%oldBufferVals(size(this%buffer%variables)))
                else if (size(this%oldBufferVals) /= size(this%buffer%variables)) then
                    deallocate(this%oldBufferVals)
                    allocate(this%oldBufferVals(size(this%buffer%variables)))
                end if
                do j = 1, size(this%convergenceTestVars)
                    k = this%convergenceTestVars(j)
                    call copyRealArrayInPlace(this%oldBufferVals(k)%entry,this%buffer%variables(k)%entry)
                end do
            end if

            aaApplied = .false.
            if (this%andersonDepth > 0) then
                call andersonExtrapolate(this,manipulatedModeller,this%implicitVectorOld,this%implicitVectorNew,&
                                         aaHasPrev,aaHistCount,aaNextSlot,aaOrder,aaApplied)
            else
                this%implicitVectorNew = this%relaxationWeight*this%implicitVectorNew &
                                        + (1.0d0-this%relaxationWeight)*this%implicitVectorOld
            end if

            call applyIterateAndCheckConvergence(this,manipulatedModeller,commNeeded,commData,&
                                                 nonTrivialConvergenceCheck,this%implicitVectorNew,&
                                                 convergenceCounter,tolReached)

            ! Guard against false AA convergence: because convergence tests step change ||x_{k+1} - x_k||,
            ! an extrapolated step can stall (x_{k+1}^{AA} ~= x_k) while the true residual f_k = g(x_k) - x_k is still large.
            ! Re-verify the unextrapolated Picard map g(x_k) = x_k + f_k. Reusing convergenceCounter is safe because
            ! the AA candidate already passed on every rank (so no counter was incremented on the first check).
            if (tolReached .and. aaApplied) then
                this%implicitVectorNew = this%implicitVectorOld + this%aaWorkF
                call applyIterateAndCheckConvergence(this,manipulatedModeller,commNeeded,commData,&
                                                     nonTrivialConvergenceCheck,this%implicitVectorNew,&
                                                     convergenceCounter,tolReached)
                if (.not. tolReached) then
                    ! Discard stalled AA subspace; when relaxationWeight /= 1, restore damped iterate x_k + w*f_k
                    ! (and its buffer/derived state) rather than advancing from the unrelaxed g(x_k) used above.
                    aaHistCount = 0
                    aaHasPrev = .false.
                    if (this%relaxationWeight /= 1.0d0) then
                        this%implicitVectorNew = this%implicitVectorOld + this%relaxationWeight * this%aaWorkF
                        call applyImplicitIterate(this,manipulatedModeller,commNeeded,commData,this%implicitVectorNew)
                    end if
                end if
            end if

            ! Guarantee the accepted iterate was verified at exact solverToleranceRel: if nonlinear convergence fired on an
            ! inexact linear solve, reject termination and force exact KSP tolerance on the next iteration(s) so conservation
            ! and fixed-point accuracy match full-tolerance solves.
            if (tolReached .and. .not. solveWasExact) then
                tolReached = .false.
                forceExactTol = .true.
            end if

            nonlinIter = nonlinIter + 1
            
            if (tolReached) exit
        end do

        if (nonTrivialConvergenceCheck) then
            do j = 1, size(convergenceCounter)
                if (convergenceCounter(j) == nonlinIter) &
                    call printMessage(this%integratorName//": convergence bottleneck: "&
                                // inputVars%getVarName(this%convergenceTestVars(j)),.true.)
            end do
        end if

        ! Call manipulator with priority 1
        call manipulatedModeller%callManipulator(1,this%buffer,this%buffer)

        ! Calculate all variables
        if (commNeeded) then 
            call manipulatedModeller%safeCommAndDeriv(commData,this%buffer)
        else
            call this%buffer%calculateDerivedVars()
        end if


    class default
        error stop "Unsupported surrogate passed to BDE integrator"
    end select

end subroutine tryIntegrate
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine applyImplicitIterate(this,manipulatedModeller,commNeeded,commData,implicitVec)
    !! Unpack an implicit iterate vector into this%buffer and run priority-0 manipulator, halo communication, and derivations.
    !! Shared by candidate evaluation, AA raw-map verification, and damped fallback so all paths leave this%buffer consistent.

    type(PicardBDEIntegrator) ,intent(inout) :: this !! Integrator whose buffer is updated in place
    type(Modeller)            ,intent(inout) :: manipulatedModeller !! Modeller providing manipulator and communication callbacks
    logical                   ,intent(in)    :: commNeeded !! True if inter-rank communication is required during Picard iterations
    type(CommunicationData)   ,intent(in)    :: commData !! Communication pattern metadata for safeCommAndDeriv
    real(rk) ,dimension(:)    ,intent(in)    :: implicitVec !! Local rank slice of the candidate implicit state vector (size nLoc)

    call this%buffer%extractImplicitVars(implicitVec)
    call manipulatedModeller%callManipulator(0,this%buffer,this%buffer)
    if (commNeeded) then
        call manipulatedModeller%safeCommAndDeriv(commData,this%buffer,derivPriority=0)
    else
        !Calculate only highest priority derived variables in internal iterations
        call this%buffer%calculateDerivedVars(derivPriority=0)
    end if

end subroutine applyImplicitIterate
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine applyIterateAndCheckConvergence(this,manipulatedModeller,commNeeded,commData,nonTrivialConvergenceCheck,&
                                           implicitVec,convergenceCounter,tolReached)
    !! Apply an implicit iterate to this%buffer via applyImplicitIterate and evaluate global nonlinear convergence across all ranks.
    !! Compares the updated state against this%oldBufferVals (or this%implicitVectorOld when no convergenceTestVars are set).

    type(PicardBDEIntegrator)              ,intent(inout) :: this !! Integrator owning buffer, oldBufferVals, implicitVectorOld, and tolerances
    type(Modeller)                         ,intent(inout) :: manipulatedModeller !! Modeller used for postprocessing and global all-true reduction
    logical                                ,intent(in)    :: commNeeded !! True if inter-rank communication is required during Picard iterations
    type(CommunicationData)                ,intent(in)    :: commData !! Communication pattern metadata for safeCommAndDeriv
    logical                                ,intent(in)    :: nonTrivialConvergenceCheck !! True to check named convergenceTestVars; false to check relative L2 norm of implicitVec
    real(rk) ,dimension(:)                 ,intent(in)    :: implicitVec !! Local rank slice of the candidate implicit state vector (size nLoc)
    integer(ik) ,allocatable ,dimension(:) ,intent(inout) :: convergenceCounter !! Per-variable failure counters for bottleneck reporting (incremented only on variables failing tol)
    logical                                ,intent(out)   :: tolReached !! True iff every MPI rank satisfies nonlinear convergence tolerance

    logical :: locConverged

    call applyImplicitIterate(this,manipulatedModeller,commNeeded,commData,implicitVec)

    if (nonTrivialConvergenceCheck) then
        locConverged = checkConvergence(this%oldBufferVals,this%buffer%variables,this%convergenceTestVars,&
                                        this%nonlinTol,this%absTol,this%use2Norm,convergenceCounter)
    else
        locConverged = (norm2(this%implicitVectorOld-implicitVec)/norm2(this%implicitVectorOld)) &
                    < this%nonlinTol
    end if
    tolReached = manipulatedModeller%isTrueEverywhere(locConverged)

end subroutine applyIterateAndCheckConvergence
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine andersonExtrapolate(this,manipulatedModeller,xOld,gAndOut,aaHasPrev,aaHistCount,aaNextSlot,aaOrder,aaApplied)
    !! Type-II Anderson Acceleration (AA) in unconstrained difference form using distributed CGS2 QR.
    !!
    !! Given current iterate x_k (= xOld) and unrelaxed Picard map g(x_k) (= gAndOut on entry), this routine:
    !!   1. Forms the raw Picard residual f_k = g(x_k) - x_k (stored in this%aaWorkF).
    !!   2. Appends the newest step difference pair (Delta x, Delta f) = (x_k - x_{k-1}, f_k - f_{k-1}) to the
    !!      circular buffers aaDeltaX/aaDeltaF (up to m = this%andersonDepth columns, tracked newest-to-oldest).
    !!   3. Builds a distributed CGS2 QR factorization DeltaF = Q * R over the nActive linearly independent
    !!      columns of DeltaF = [Delta f_{slot(1)}, ..., Delta f_{slot(nActive)}] (newest first, dropping
    !!      tiny or nearly collinear columns).
    !!   4. Solves min_{gamma} ||f_k - DeltaF * gamma||_2 (i.e. R * gamma = Q^T * f_k) for the least-squares
    !!      coefficient vector gamma(1:nActive).
    !!   5. Forms the accelerated iterate with relaxation weight w = this%relaxationWeight:
    !!        x_{k+1} = (x_k + w * f_k) - sum_{p=1}^{nActive} gamma(p) * (Delta x_{slot(p)} + w * Delta f_{slot(p)}).
    !!
    !! Because the unconstrained formulation works on step differences, each gamma(p) adds -gamma(p) and +gamma(p)
    !! to a pair of consecutive iterates in the underlying affine combination of previous damped Picard steps
    !! (x_j + w * f_j). The total weight therefore always sums to 1, even when collinear columns are dropped.
    !!
    !! On exit: gAndOut holds x_{k+1} (or damped Picard x_k + w * f_k if history is empty, all columns are dropped,
    !! or gamma fails conditioning guards), and this%aaWorkF retains raw f_k = g(x_k) - x_k for convergence
    !! verification and fallback in tryIntegrate.

    type(PicardBDEIntegrator) ,intent(inout) :: this !! Integrator holding relaxationWeight, andersonDepth, and AA workspace arrays
    type(Modeller)            ,intent(inout) :: manipulatedModeller !! Modeller providing global vector sum reductions (getGlobalSumVec)
    real(rk) ,dimension(:)    ,intent(in)    :: xOld !! Current nonlinear iterate x_k on this rank (size nLoc)
    real(rk) ,dimension(:)    ,intent(inout) :: gAndOut !! On entry: raw Picard solve g(x_k); on exit: proposed next iterate x_{k+1}
    logical                   ,intent(inout) :: aaHasPrev !! True once aaPrevX/aaPrevF hold a valid previous (x_{k-1}, f_{k-1}) pair; reset on fallback
    integer(ik)               ,intent(inout) :: aaHistCount !! Number of valid difference columns stored in aaDeltaX/aaDeltaF (0 .. m)
    integer(ik)               ,intent(inout) :: aaNextSlot !! Circular buffer write pointer (1 .. m) for the newest difference column
    integer(ik) ,dimension(:) ,intent(inout) :: aaOrder !! Slot indices ordered newest-to-oldest (aaOrder(1) = newest); size m
    logical                   ,intent(out)   :: aaApplied !! True iff an AA correction was actually applied (false on damped Picard fallback)

    integer(ik) :: m ,j ,p ,pass ,slot ,nActive
    real(rk)    :: normV ,normRes ,normFloor ,rhsVal ,maxGammaSum ,locNorm(1) ,globNorm(1)
    logical     :: solveOk
    real(rk) ,allocatable ,dimension(:)   :: v ,locDots ,globDots ,rCol ,rhs ,gamma
    real(rk) ,allocatable ,dimension(:,:) :: R
    integer(ik) ,allocatable ,dimension(:) :: activeSlots

    m = this%andersonDepth
    this%aaWorkF = gAndOut - xOld

    ! Overwrite the oldest circular slot without shifting nLoc-sized vectors; keep aaOrder sorted newest-to-oldest.
    if (aaHasPrev) then
        aaNextSlot = mod(aaNextSlot,m) + 1
        this%aaDeltaX(:,aaNextSlot) = xOld - this%aaPrevX
        this%aaDeltaF(:,aaNextSlot) = this%aaWorkF - this%aaPrevF
        do j = min(aaHistCount + 1,m), 2, -1
            aaOrder(j) = aaOrder(j - 1)
        end do
        aaOrder(1) = aaNextSlot
        aaHistCount = min(aaHistCount + 1,m)
    end if

    this%aaPrevX = xOld
    this%aaPrevF = this%aaWorkF
    aaHasPrev = .true.

    ! Pre-populate damped Picard x_k + w*f_k so any early exit or guard failure falls back cleanly.
    gAndOut = xOld + this%relaxationWeight * this%aaWorkF
    aaApplied = .false.
    if (aaHistCount <= 0) return

    allocate(v(size(xOld)))
    allocate(locDots(aaHistCount),globDots(aaHistCount),rCol(aaHistCount),rhs(aaHistCount),gamma(aaHistCount))
    allocate(R(aaHistCount,aaHistCount),activeSlots(aaHistCount))

    nActive = 0
    R = 0.0d0
    ! normFloor guards against zero/near-underflow norms; maxGammaSum caps L1 extrapolation amplification.
    normFloor = sqrt(tiny(1.0d0))
    maxGammaSum = 1.0d3

    ! Build orthonormal basis Q for DeltaF newest-to-oldest so fresher secant directions take priority
    ! when older columns are linearly dependent.
    do j = 1, aaHistCount
        slot = aaOrder(j)
        v = this%aaDeltaF(:,slot)
        locNorm(1) = dot_product(v,v)
        call manipulatedModeller%getGlobalSumVec(locNorm,globNorm)
        if (.not. ieee_is_finite(globNorm(1))) cycle
        normV = sqrt(max(globNorm(1),0.0d0))
        if (normV <= normFloor) cycle

        ! Two-pass Classical Gram-Schmidt (CGS2): restores near-machine orthogonality while batching
        ! all nActive inner products into 2 MPI all-reduces per column instead of nActive sequential reductions.
        if (nActive > 0) then
            rCol(1:nActive) = 0.0d0
            do pass = 1, 2
                do p = 1, nActive
                    locDots(p) = dot_product(this%aaQ(:,p),v)
                end do
                call manipulatedModeller%getGlobalSumVec(locDots(1:nActive),globDots(1:nActive))
                rCol(1:nActive) = rCol(1:nActive) + globDots(1:nActive)
                do p = 1, nActive
                    v = v - globDots(p) * this%aaQ(:,p)
                end do
            end do
        end if

        locNorm(1) = dot_product(v,v)
        call manipulatedModeller%getGlobalSumVec(locNorm,globNorm)
        if (.not. ieee_is_finite(globNorm(1))) cycle
        normRes = sqrt(max(globNorm(1),0.0d0))

        ! Drop nearly collinear columns (||v_perp|| <= 1e-10 * ||v||) to keep R well-conditioned.
        if (normRes > 1.0d-10 * normV .and. normRes > normFloor) then
            nActive = nActive + 1
            activeSlots(nActive) = slot
            this%aaQ(:,nActive) = v / normRes
            if (nActive > 1) R(1:nActive-1,nActive) = rCol(1:nActive-1)
            R(nActive,nActive) = normRes
        end if
    end do

    if (nActive <= 0) return

    ! Project raw residual f_k onto Q; R and rhs are small (<= m) and replicated identically on every rank.
    do p = 1, nActive
        locDots(p) = dot_product(this%aaQ(:,p),this%aaWorkF)
    end do
    call manipulatedModeller%getGlobalSumVec(locDots(1:nActive),rhs(1:nActive))

    ! Guarded back-substitution for R * gamma = rhs: test |rhsVal| <= maxGammaSum * R(p,p) before dividing
    ! to prevent division blow-up or overflow FPE traps on small diagonals.
    solveOk = .true.
    gamma(1:nActive) = 0.0d0
    do p = nActive, 1, -1
        if (p < nActive) then
            rhsVal = rhs(p) - dot_product(R(p,p+1:nActive),gamma(p+1:nActive))
        else
            rhsVal = rhs(p)
        end if
        if (.not. ieee_is_finite(rhsVal) .or. abs(rhsVal) > maxGammaSum * R(p,p)) then
            solveOk = .false.
            exit
        end if
        gamma(p) = rhsVal / R(p,p)
    end do

    if (solveOk .and. all(ieee_is_finite(gamma(1:nActive))) .and. sum(abs(gamma(1:nActive))) <= maxGammaSum) then
        do p = 1, nActive
            slot = activeSlots(p)
            gAndOut = gAndOut - gamma(p) * (this%aaDeltaX(:,slot) + this%relaxationWeight * this%aaDeltaF(:,slot))
        end do
        aaApplied = .true.
    else
        ! gAndOut already holds damped Picard x_k + w*f_k; clear ill-conditioned history.
        aaHistCount = 0
        aaHasPrev = .false.
    end if

end subroutine andersonExtrapolate
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule implicit_PicardBDE_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
