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
    modelList,termGroups,evolvesTimeVar,dtController,initialTimestep&
    ,use2Norm,petscGroup,intContOptions,integratorName,relaxationWeight)
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
    
    integer(ik) ,allocatable ,dimension(:) :: procDoFs

    if (assertions .or. assertionLvl >= 0) then
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
    allocate(this%implicitVectorInit(procDoFs(procRank+1)))

    this%implicitVectorNew = 0
    this%implicitVectorOld = 0
    this%implicitVectorInit = 0

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

    this%relaxationWeight = 1.0d0 
    if (present(relaxationWeight)) this%relaxationWeight = relaxationWeight
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
   
        outputVars%variables = this%buffer%variables
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
        outputVars%variables = inputVars%variables
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

    logical :: conv 
    logical ,allocatable ,dimension(:) :: varConverged

    integer(ik) :: i ,haloDataChunkSize ,ind ,nonHaloLen
    real(rk) :: relError ,absError

    allocate(varConverged(size(indicesToCheck)))
    varConverged = .false. 
    do i = 1, size(indicesToCheck)
        ind = indicesToCheck(i)
        haloDataChunkSize = 1 - lbound(newVars(ind)%entry,1)
        nonHaloLen = size(newVars(ind)%entry) - 2*haloDataChunkSize

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
    logical                 :: tolReached ,locConverged

    type(CommunicationData) :: commData

    
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

    select type (manipulatedModeller)
    type is (Modeller)
        if (.not. allocated(this%buffer)) allocate(this%buffer,source=inputVars)
        if (.not. allocated(this%oldBufferVals)) allocate(this%oldBufferVals,source=inputVars%variables)
        this%buffer%variables = inputVars%variables
        this%oldBufferVals = inputVars%variables
        
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

            call manipulatedModeller%linearSolvePETSc(this%implicitVectorInit,this%implicitVectorNew,.true.,&
                                                    convReason,this%associatedPETScObjGroup)

            ! Check if PETSc diverged 
            if (convReason < -1) then 
                call printNamedValue(this%integratorName//": PETScConvergedReason:",convReason)
                solveSuccess = .false.
                if (.not. this%internalStepControl) then
                    error stop "PETSc convergence reason negative and internal step control disabled"
                end if
                return
            end if
            this%oldBufferVals = this%buffer%variables

            this%implicitVectorNew = this%relaxationWeight*this%implicitVectorNew &
                                    + (1.0d0-this%relaxationWeight)*this%implicitVectorOld
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
                checkConvergence(this%oldBufferVals,this%buffer%variables,this%convergenceTestVars,&
                                this%nonlinTol,this%absTol,this%use2Norm,convergenceCounter)
            else
                locConverged = (norm2(this%implicitVectorOld-this%implicitVectorNew)/norm2(this%implicitVectorOld)) &
                            < this%nonlinTol
            end if
            tolReached =  manipulatedModeller%isTrueEverywhere(locConverged)

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
end submodule implicit_PicardBDE_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
