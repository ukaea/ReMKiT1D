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

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initBDEIntegrator(this,indexingObj,procRank,nonlinTol,absTol,maxIters,convergenceIndices,&
    modelList,termGroups,evolvesTimeVar,dtController,initialTimestep,use2Norm,petscGroup,intContOptions,integratorName)
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
    
    integer(ik) ,allocatable ,dimension(:) :: procDoFs

    if (assertions) then
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to BDE integrator constructor")
        if (present(termGroups)) call assert(present(modelList),"Term groups object passed to BDE integrator constructor without&
                                                               & model list")
        if (present(modelList) .and. present(termGroups)) call assert(size(modelList)==size(termGroups),"If both model list and &
        & term groups arguments are passed to initBDEIntegrator both must have the same length")
        if (present(dtController)) call assert(dtController%isDefined(),"Undefined timestep controller passed to BDE&
        & integrator constructor")
    end if

    procDoFs = indexingObj%getProcDoF()

    allocate(this%implicitVectorNew(procDoFs(procRank+1)))
    allocate(this%implicitVectorOld(procDoFs(procRank+1)))

    this%implicitVectorNew = 0
    this%implicitVectorOld = 0

    this%timesCalled = 0 
    this%totNumIters = 0

    !Default values
    this%nonlinTol = real(1.0d-12,kind=rk)
    this%maxIterations = 100
    this%absTol = real(1,kind=rk)
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
    call this%makeDefined()

end subroutine initBDEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine integrateBDE(this,manipulatedModeller,outputVars,inputVars) 
    !! Integration routine for BDE integrator - implementation of abstract manipulate routine

    class(PicardBDEIntegrator)            ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
    class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
    class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    integer(ik)                         :: numSteps 
    real(rk) ,allocatable ,dimension(:) :: dt
    real(rk)                            :: fullTimestep, controllerTimestep 

    logical                 :: solveSuccess


    if (assertions) then 
        call assert(this%isDefined(),"Integration routine called by undefined BDE integrator")
        call assert(manipulatedModeller%isDefined(),"Undefined modeller object passed to BDE integration routine")
        call assert(outputVars%isDefined(),"outputVars passed to integrateBDE not defined")
        call assert(inputVars%isDefined(),"inputVars passed to integrateBDE not defined")
    end if

    fullTimestep = this%getTimestep()

    do 
        solveSuccess = .true.
        !Calculate the controller timestep, and if it is smaller than the full timestep divide the integration into numSteps > 1
        if (this%hasTimestepController() .and. .not. this%internalStepControl) then 
            controllerTimestep = this%getTimestepFromController(inputVars)
            numSteps = ceiling(fullTimestep/controllerTimestep)
            allocate(dt(numSteps))
            if (numSteps > 1) then 
                dt(1:numSteps-1) = controllerTimestep
                dt(numSteps) = fullTimestep - sum(dt(1:numSteps-1))
            else
                dt = fullTimestep
            end if
        else if (this%internalStepControl) then
            numSteps = this%internalControlOpts%currentNumSubsteps
            if (allocated(dt)) deallocate(dt)
            allocate(dt(numSteps))
            dt = fullTimestep/numSteps
        else
            numSteps = 1
            allocate(dt(1))
            dt = fullTimestep
        end if
        if (this%internalControlOpts%restartCount > this%internalControlOpts%hardMaxRestarts) &
            error stop "Max BDE restarts reached - hard maximum"
        if (this%internalControlOpts%restartCount > this%internalControlOpts%maxRestarts .and. &
            this%internalControlOpts%stepsSinceLastConsolidation > 1) error stop "Max BDE restarts reached"
        call tryIntegrate(this,manipulatedModeller,outputVars,inputVars,numSteps,dt,solveSuccess)
        if (solveSuccess) then
            this%internalControlOpts%restartCount = 0
            this%internalControlOpts%stepsSinceLastConsolidation = this%internalControlOpts%stepsSinceLastConsolidation + 1

            if (this%internalControlOpts%stepsSinceLastConsolidation >= this%internalControlOpts%consolidationInterval) then 
                this%internalControlOpts%currentNumSubsteps = 1
                this%internalControlOpts%stepsSinceLastConsolidation = 0
                call printMessage(this%integratorName//": Consolidating internal BDE steps")
            end if
            exit
        end if 
    end do
   
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

    logical :: conv 
    logical ,allocatable ,dimension(:) :: varConverged

    integer(ik) :: i ,haloDataChunkSize ,ind ,nonHaloLen
    real(rk) :: relError ,absError

    allocate(varConverged(size(indicesToCheck)))
    varConverged = .false. 
    do i = 1, size(indicesToCheck)
        ind = indicesToCheck(i)
        haloDataChunkSize = 1 - lbound(newVars(ind)%entry,1)
        nonHaloLen = size(newVars(ind)%entry) - haloDataChunkSize

        if (use2Norm) then

            absError = norm2(oldVars(ind)%entry(1:nonHaloLen)-newVars(ind)%entry(1:nonHaloLen))

            if (norm2(oldVars(ind)%entry(1:nonHaloLen))>epsilon(absError)*absTol) then
                relError = norm2(oldVars(ind)%entry(1:nonHaloLen)-newVars(ind)%entry(1:nonHaloLen))&
                        /norm2(oldVars(ind)%entry(1:nonHaloLen))
            else
                relError = absError
            end if
        else

            absError = maxval(abs((oldVars(ind)%entry(1:nonHaloLen)-newVars(ind)%entry(1:nonHaloLen))))

            if (all(abs(oldVars(ind)%entry(1:nonHaloLen)) > epsilon(absError)*absTol)) then 
                relError = maxval(abs((oldVars(ind)%entry(1:nonHaloLen)-newVars(ind)%entry(1:nonHaloLen)))&
                    /abs(oldVars(ind)%entry(1:nonHaloLen)))
                else
                    relError = absError
            end if
            
        end if
        varConverged(i) = relError < nonlinTol .or. absError < epsilon(absError)*absTol
        if (.not. varConverged(i)) convergenceCounter(i) = convergenceCounter(i) + 1
    end do

    conv = all(varConverged)

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
subroutine tryIntegrate(this,manipulatedModeller,outputVars,inputVars,numSteps,dt,solveSuccess) 

    type(PicardBDEIntegrator)            ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
    class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
    class(VariableContainer)  ,value      ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    integer(ik)                         ,intent(in)      :: numSteps 
    real(rk)              ,dimension(:) ,intent(in)      :: dt
    logical                             ,intent(inout)   :: solveSuccess

    type(IntArray) ,allocatable ,dimension(:)  :: termGroups 
    integer(ik)    ,allocatable ,dimension(:)  :: modelIndices 

    type(LogicalArray) ,allocatable ,dimension(:) :: updateRules
    logical            ,allocatable ,dimension(:) :: modelDataUpdateRules
    integer(ik)                                   :: timeVarIndex ,i ,j ,k , nonlinIter ,convReason

    logical                 :: commNeeded
    logical                 :: nonTrivialUpdate
    logical                 :: nonTrivialModelDataUpdate
    logical                 :: nonTrivialConvergenceCheck
    logical                 :: tolReached ,locConverged
    logical                 :: timeEvolving

    type(CommunicationData) :: commData

    type(RealArray) ,allocatable ,dimension(:) :: oldBufferVals
    real(rk)        ,allocatable ,dimension(:) :: implicitVectorInit

    real (rk) :: startingTime 
    
    integer(ik)    ,allocatable ,dimension(:)  :: convergenceCounter 

    termGroups = this%getTermGroups()
    modelIndices = this%getModelIndices()
    nonTrivialUpdate = this%hasNonTrivialUpdate()
    if (nonTrivialUpdate) updateRules = this%getUpdateRules()
    nonTrivialModelDataUpdate = this%hasNonTrivialModelDataUpdate()
    if (nonTrivialModelDataUpdate) modelDataUpdateRules = this%getModelDataUpdateRules()
    nonTrivialConvergenceCheck = allocated(this%convergenceTestVars) 

    if (nonTrivialConvergenceCheck) allocate(convergenceCounter(size(this%convergenceTestVars)))

    !Check if communication needed
    commNeeded = this%isCommunicationNeeded()
    if (commNeeded) commData = this%getCommunicationData()

    timeEvolving = this%isTimeEvolving()
    select type (manipulatedModeller)
    type is (Modeller)
        if (.not. allocated(this%buffer))allocate(this%buffer,source=inputVars)
        allocate(oldBufferVals,source=inputVars%variables)
        this%buffer%variables = inputVars%variables
        
        if (inputVars%isVarNameRegistered("time")) then 
            timeVarIndex = inputVars%getVarIndex("time")
            startingTime = this%buffer%variables(timeVarIndex)%entry(1)
        end if

        allocate(implicitVectorInit,source=this%implicitVectorOld)
        do i = 1, numSteps
            if (inputVars%isVarNameRegistered("time")) &
            this%buffer%variables(timeVarIndex)%entry(1) = this%buffer%variables(timeVarIndex)%entry(1) + dt(i)
            tolReached = .false.
            if (nonTrivialConvergenceCheck) convergenceCounter = 1

            call this%buffer%copyImplicitVarsToVec(implicitVectorInit,ignoreStationary=.true.)
            do nonlinIter = 1, this%maxIterations
                call this%buffer%copyImplicitVarsToVec(this%implicitVectorOld)
                this%implicitVectorNew = this%implicitVectorOld

                do j = 1,size(modelIndices)
                    if (nonlinIter == 1) then 
                        call manipulatedModeller%updateModelData(modelIndices(j),this%buffer)
                    else if (nonTrivialModelDataUpdate) then 
                        if (modelDataUpdateRules(j)) &
                        call manipulatedModeller%updateModelData(modelIndices(j),this%buffer,updatePriority=0)
                    end if
                    do k = 1,size(termGroups(j)%entry)
                        if (nonTrivialUpdate) then 
                            if (updateRules(j)%entry(k))&
                            call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                        else if (nonlinIter == 1) then 
                            call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                        end if
                        call manipulatedModeller%calculateMatGroupValsInModel(modelIndices(j),&
                                                                                termGroups(j)%entry(k),this%buffer)
                        call manipulatedModeller%addModelMatGroupToPETSc(modelIndices(j),termGroups(j)%entry(k),-dt(i),&
                                                                        this%associatedPETScObjGroup)
    
                    end do
                end do

                call manipulatedModeller%linearSolvePETSc(implicitVectorInit,this%implicitVectorNew,.true.,&
                                                        convReason,this%associatedPETScObjGroup)

                ! Check if PETSc diverged 
                if (convReason < -1) then 
                    if (this%internalStepControl) then

                        call printNamedValue(this%integratorName//": PETScConvergedReason:",convReason)
                        call printMessage("Applying internal step control and restarting BDE integration")
                        this%internalControlOpts%currentNumSubsteps = this%internalControlOpts%currentNumSubsteps &
                                                                    * this%internalControlOpts%stepMultiplier
                        solveSuccess = .false.
                        this%internalControlOpts%restartCount = this%internalControlOpts%restartCount + 1
                        return
                    else
                        call printNamedValue(this%integratorName//": PETScConvergedReason:",convReason)
                        error stop "PETSc convergence reason negative and internal step control disabled"
                    end if
                end if
                oldBufferVals = this%buffer%variables
                call this%buffer%extractImplicitVars(this%implicitVectorNew)
                call manipulatedModeller%callManipulator(0,this%buffer,this%buffer)
                if (commNeeded) then 
                    call manipulatedModeller%safeCommAndDeriv(commData,this%buffer,derivPriority=0)
                else
                    !Calculate only highest priority derived variables in internal iterations
                    call this%buffer%calculateDerivedVars(derivPriority=0)
                end if

                if (nonTrivialConvergenceCheck) then
                    locConverged = &
                    checkConvergence(oldBufferVals,this%buffer%variables,this%convergenceTestVars,&
                                    this%nonlinTol,this%absTol,this%use2Norm,convergenceCounter)
                else
                    locConverged = (norm2(this%implicitVectorOld-this%implicitVectorNew)/norm2(this%implicitVectorOld)) &
                                < this%nonlinTol
                end if
                tolReached =  manipulatedModeller%isTrueEverywhere(locConverged)

                if (tolReached) exit
            end do
            if (nonlinIter > this%maxIterations) then 
                call printMessage("WARNING: "//this%integratorName//" reached maximum number of iterations")

                if (this%internalStepControl) then
                    call printMessage("Applying internal step control and restarting BDE integration")
                    this%internalControlOpts%currentNumSubsteps = this%internalControlOpts%currentNumSubsteps &
                                                                * this%internalControlOpts%stepMultiplier
                    solveSuccess = .false.
                    this%internalControlOpts%restartCount = this%internalControlOpts%restartCount + 1
                    return
                end if
            end if

            if (nonTrivialConvergenceCheck) then
                do j = 1, size(convergenceCounter)
                    if (convergenceCounter(j) == nonlinIter) &
                        call printMessage(this%integratorName//": convergence bottleneck: "&
                                    // inputVars%getVarName(this%convergenceTestVars(j)),.true.)
                end do
            end if

            if (solveSuccess) then
                this%totNumIters = this%totNumIters + nonlinIter
                this%timesCalled = this%timesCalled + 1
                call printNamedValue(this%integratorName//": last number of iterations",nonlinIter)
                call printNamedValue(this%integratorName//": total number of iterations",this%totNumIters)
                call printNamedValue(this%integratorName//": average number of iterations"&
                                    ,real(this%totNumIters,kind=rk)/this%timesCalled)
            end if

            if (this%internalStepControl .and. solveSuccess) then
                if (nonlinIter < this%internalControlOpts%minNonlinIters &
                    .and. this%internalControlOpts%currentNumSubsteps > 1) then
                    call printMessage(this%integratorName//&
                    ": Minumum nonlinear iterations reached - attempting to reduce number of BDE substeps")
                    this%internalControlOpts%currentNumSubsteps = max(this%internalControlOpts%currentNumSubsteps &
                                                                - this%internalControlOpts%stepDecrament,1)
                end if
            end if
            ! Call manipulator with priority 1
            call manipulatedModeller%callManipulator(1,this%buffer,this%buffer)

            ! Calculate all variables
            if (commNeeded) then 
                call manipulatedModeller%safeCommAndDeriv(commData,this%buffer)
            else
                call this%buffer%calculateDerivedVars()
            end if

        end do

        outputVars%variables = this%buffer%variables

        if (inputVars%isVarNameRegistered("time") .and. (.not. timeEvolving)) &
            outputVars%variables(timeVarIndex)%entry = startingTime
        call manipulatedModeller%callManipulator(2,outputVars,outputVars)
    class default
        error stop "Unsupported surrogate passed to BDE integrator"
    end select

end subroutine tryIntegrate
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule implicit_PicardBDE_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
