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

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrateCVODE(this,manipulatedModeller,outputVars,inputVars) 
        !! Implementation of abstract manipulate routine for the case of Runge-Kutta integrator

        class(CVODEIntegrator)                ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

        type(IntArray) ,allocatable ,dimension(:)  :: termGroups 
        integer(ik)    ,allocatable ,dimension(:)  :: modelIndices 

        integer(ik)                         :: numSteps ,i ,j ,k ,l 
        real(rk) ,allocatable ,dimension(:) :: dt
        real(rk)                            :: fullTimestep, controllerTimestep, currentTimeVal

        type(IntArray)     ,allocatable ,dimension(:) :: varIndices
        type(LogicalArray) ,allocatable ,dimension(:) :: updateRules
        type(LogicalArray) ,allocatable ,dimension(:) :: mixedGroup
        integer(ik)                                   :: timeVarIndex
        logical            ,allocatable ,dimension(:) :: modelDataUpdateRules
        integer(ik) ,allocatable ,dimension(:) :: indices 

        logical                 :: commNeeded
        logical                 :: nonTrivialUpdate
        logical                 :: nonTrivialModelDataUpdate
        logical                 :: timeEvolving
        type(CommunicationData) :: commData

        integer(c_int) :: ierr
        integer(c_int) :: maxres 
        integer(c_long) :: nlocal, neq, mu, ml, mudq, mldq, numCVODESteps(1)
        integer(c_int) :: iPretype0 = 1 !! Preconditioner type (left)
        integer(c_int) :: iGStype = 1 !! Gram-Schmidt orthoganlization type (classical)
        real(c_double) :: t(1) ,rtol ,atol  

        termGroups = this%getTermGroups()
        modelIndices = this%getModelIndices()
        nonTrivialUpdate = this%hasNonTrivialUpdate()
        if (nonTrivialUpdate) updateRules = this%getUpdateRules()
        nonTrivialModelDataUpdate = this%hasNonTrivialModelDataUpdate()
        if (nonTrivialModelDataUpdate) modelDataUpdateRules = this%getModelDataUpdateRules()
    
        if (assertions) then 
            call assert(this%isDefined(),"Integration routine called by undefined CVODE integrator")
            call assert(manipulatedModeller%isDefined(),"Undefined Modeller object passed to CVODE integration routine")
            call assert(outputVars%isDefined(),"outputVars passed to integrateCVODE not defined")
            call assert(inputVars%isDefined(),"inputVars passed to integrateCVODE not defined")
        end if

        fullTimestep = this%getTimestep()
        timeEvolving = this%isTimeEvolving()
        !Calculate the controller timestep, and if it is smaller than the full timestep divide the integration into numSteps > 1
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

        !Check if communication needed
        commNeeded = this%isCommunicationNeeded()
        if (commNeeded) commData = this%getCommunicationData()

        currentTimeVal = 0
        if (inputVars%isVarNameRegistered("time")) then
            timeVarIndex = inputVars%getVarIndex("time")
            currentTimeVal = inputVars%variables(timeVarIndex)%entry(1) 
        end if
        select type (manipulatedModeller)
        type is (Modeller)
            
            allocate(varIndices(size(modelIndices)))
            do i = 1, size(modelIndices)
                allocate(varIndices(i)%entry(size(termGroups(i)%entry)))
                do j = 1, size(termGroups(i)%entry)
                    call assert (.not. manipulatedModeller%isModelTermGroupMixed(termGroups(i)%entry(j),modelIndices(i)),&
                        "CVODE integrator does not support mixed term groups")
                    varIndices(i)%entry(j) = inputVars%getVarIndex(manipulatedModeller&
                                                                  %getEvolvedVarInTermGroup(termGroups(i)%entry(j),modelIndices(i)))
                    call assert(.not.inputVars%isStationary(manipulatedModeller&
                    %getEvolvedVarInTermGroup(termGroups(i)%entry(j),modelIndices(i))),&
                        "CVODE integrator detected stationary variable among evolved variables - this is unsupported")
                    
                end do
            end do

            mu = this%options%bbdmukeep
            ml = this%options%bbdmlkeep
            mudq = this%options%bbdmudq
            mldq = this%options%bbdmldq
            if (.not. allocated(this%evolvedVars)) then 
                allocate(indices(0))
                do i = 1,size(modelIndices) 
                    do j = 1,size(termGroups(i)%entry)
                       indices = [indices,varIndices(i)%entry(j)]

                   end do
               end do

               indices = removeDupeInts(indices) 
               allocate(this%evolvedVars(size(indices)))
               do i = 1,size(indices) 
                   this%evolvedVars(i)%string = inputVars%getVarName(indices(i))

               end do


               allocate(this%bufferRHS,source = inputVars) 
               allocate(this%bufferY,source = inputVars)
               this%bufferY%variables = inputVars%variables
               this%bufferRHS%variables = inputVars%variables
               
               call this%bufferY%copyNamedVarsToVec(this%copyBufferVals,this%evolvedVars)

               nlocal = size(this%copyBufferVals)
               neq = nlocal * this%numProcs
               this%sunVecY => &
               FN_VNew_Parallel(this%mpiComm%MPI_VAL,nlocal,neq,this%sunctx)
               this%sunVecYDot => &
               FN_VNew_Parallel(this%mpiComm%MPI_VAL,nlocal,neq,this%sunctx)
               
               this%rhsVec(1:nlocal) => FN_VGetArrayPointer(this%sunVecYDot)
               this%yVec(1:nlocal) => FN_VGetArrayPointer(this%sunVecY)

               this%cvode = FCVodeCreate(CV_BDF, this%sunctx)
               
               ierr = FCVodeInit(this%cvode,c_funloc(rhs), currentTimeVal, this%sunVecY)

               this%sunls => FSUNLinSol_SPGMR(this%sunVecY,iPretype0,0,this%sunctx)
               this%sunmat => null()

               ierr = FCVodeSetLinearSolver(this%cvode,this%sunls,this%sunmat)

               rtol = this%options%reltol
               atol = this%options%abstol

               ierr = FSUNLinSol_SPGMRSetGSType(this%sunls,iGStype)

               maxres = this%options%maxRestarts

               ierr = FSUNLinSol_SPGMRSetMaxRestarts(this%sunls,maxres)

               ierr = FCVodeSStolerances(this%cvode,rtol,atol)

               ierr = FCVBBDPrecInit(this%cvode, nlocal, mudq, mldq, mu, ml, 0.d0, &
                                        c_funloc(rhsGFun), c_null_ptr)
            else 
                this%bufferY%variables = inputVars%variables
            end if

            
            !Perform start-of-timestep updates
            do j = 1,size(modelIndices)
                if (nonTrivialModelDataUpdate) then 
                    if (.not. modelDataUpdateRules(j)) &
                    call manipulatedModeller%updateModelData(modelIndices(j),this%bufferY,updatePriority=0)
                end if
                do k = 1,size(termGroups(j)%entry)
                    if (nonTrivialUpdate) then 
                        if (.not. updateRules(j)%entry(k))&
                            call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%bufferY)
                    end if
                end do
            end do

            call manipulatedModeller%callManipulator(2,this%bufferY,this%bufferY)
            outputVars%variables = this%bufferY%variables
            if (commNeeded) then 
                call manipulatedModeller%safeCommAndDeriv(commData,outputVars)
            else
                !A little bit overkill - could specialize to update only those variables that could have changed
                call outputVars%calculateDerivedVars()
            end if 

            call this%bufferY%copyNamedVarsToVec(this%copyBufferVals,this%evolvedVars)        
            this%yVec => FN_VGetArrayPointer(this%sunvecY)
            this%yVec = this%copyBufferVals
            ierr = FCVodeReInit(this%cvode, currentTimeVal, this%sunVecY)
            ierr = FCVBBDPrecReInit(this%cvode,mudq,mldq,0.0d0)

            t(1) = currentTimeVal
            do i=1,numSteps
                ierr = FCVode(this%cvode,t(1)+dt(i),this%sunVecY,t,CV_NORMAL)
                
                ierr = FCVodeGetNumSteps(this%cvode,numCVODESteps)
                call printNamedValue(this%integratorName//": CVODE number of steps",int(numCVODESteps(1),kind=ik))
                if (numSteps > 1) then
                    this%yVec => FN_VGetArrayPointer(this%sunvecY)
                    this%copyBufferVals = this%yVec

                    call printMessage(this%integratorName//": Consolidating internal BDE steps")
                    call this%bufferY%copyNamedVarsFromVec(this%copyBufferVals,this%evolvedVars)
                    
                    if (this%bufferY%isVarNameRegistered("time")) &
                        this%bufferY%variables(timeVarIndex)%entry(1) = t(1) 
                    call manipulatedModeller%callManipulator(1,this%bufferY,this%bufferY)
                    if (commNeeded) then 
                        call manipulatedModeller%safeCommAndDeriv(commData,this%bufferY)
                    else
                        call this%bufferY%calculateDerivedVars()
                    end if 
                    call this%bufferY%copyNamedVarsToVec(this%copyBufferVals,this%evolvedVars)        
                    
                    this%yVec = this%copyBufferVals
                end if
            end do

            this%yVec => FN_VGetArrayPointer(this%sunvecY)
            this%copyBufferVals = this%yVec

            call this%bufferY%copyNamedVarsFromVec(this%copyBufferVals,this%evolvedVars)


            if (this%bufferY%isVarNameRegistered("time")) then
                    this%bufferY%variables(timeVarIndex)%entry(1) = currentTimeVal
                    if (timeEvolving) this%bufferY%variables(timevarIndex)%entry(1) = t(1)
            end if
            outputVars%variables = this%bufferY%variables
            call manipulatedModeller%callManipulator(2,outputVars,outputVars)
            if (commNeeded) then 
                call manipulatedModeller%safeCommAndDeriv(commData,outputVars)
            else
                call outputVars%calculateDerivedVars()
            end if 
            
            
        class default
            error stop "Unsupported surrogate passed to CVODE integrator"
        end select
        
    contains
    integer(c_int) function rhs(t, sunvecY, sunvecYdot, uData) result(retval) bind(C) 
      real(c_double), value :: t            ! current time
      type(N_Vector)        :: sunvecY     ! solution N_Vector
      type(N_Vector)        :: sunvecYdot  ! rhs N_Vector
      type(c_ptr)           :: uData    ! user-defined data
        
      select type (manipulatedModeller)
      class is (Modeller)
        this%rhsVec => FN_VGetArrayPointer(sunvecYdot)
        this%yVec => FN_VGetArrayPointer(sunvecY)

        this%copyBufferVals = this%yVec

        call this%bufferY%copyNamedVarsFromVec(this%copyBufferVals,this%evolvedVars)

        
        if (this%bufferY%isVarNameRegistered("time")) then
            this%bufferY%variables(timeVarIndex)%entry(1) = t
        end if

        call manipulatedModeller%callManipulator(0,this%bufferY,this%bufferY) 
        if (commNeeded) then 
            call manipulatedModeller%safeCommAndDeriv(commData,this%bufferY,derivPriority=0)
        else
            call this%bufferY%calculateDerivedVars(derivPriority=0)
        end if 

        call this%bufferRHS%zeroVars(this%evolvedVars)

        do j = 1,size(modelIndices)
            if (nonTrivialModelDataUpdate) then 
                if (modelDataUpdateRules(j)) &
                call manipulatedModeller%updateModelData(modelIndices(j),this%bufferY,updatePriority=0)
            end if
            do k = 1,size(termGroups(j)%entry)
                if (nonTrivialUpdate) then 
                    if (updateRules(j)%entry(k))&
                    call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%bufferY)
                end if
                call manipulatedModeller%calculateMatGroupValsInModel(modelIndices(j)&
                                                                     ,termGroups(j)%entry(k)&
                                                                     ,this%bufferY)
                this%bufferRHS%variables(varIndices(j)%entry(k))%entry = &
                this%bufferRHS%variables(varIndices(j)%entry(k))%entry + & 
                manipulatedModeller%evaluateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%bufferY)
            end do
        end do
        
        call this%bufferRHS%copyNamedVarsToVec(this%copyBufferVals,this%evolvedVars)        
        
        this%rhsVec = this%copyBufferVals

      class default
        error stop "Unsupported surrogate passed to CVODE integrator"
      end select

      retval = 0

    end function rhs

    integer(c_int) function rhsGFun(nnlocal, t, sunvecY, sunvecG, uData) &
      result(retval) bind(C)

      real(c_double), value :: t            ! current time
      integer(c_long)       :: nnlocal      ! local space
      type(N_Vector)        :: sunvecY     ! solution N_Vector
      type(N_Vector)        :: sunvecG     ! output g N_Vector
      type(c_ptr)           :: uData    ! user-defined data

      integer :: ierrLoc

      ierrLoc = rhs(t, sunvecY, sunvecG, uData)

      retval = 0              ! Return with success
    end function rhsGFun
    
    end subroutine integrateCVODE
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

        integer(c_int) :: ierr 

        this%mpiComm = mpiCont%getComm()
        this%numProcs = mpiCont%getWorldSize()

        ierr = FSUNContext_Create(this%mpiComm%MPI_VAL,this%sunctx) 

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

        type(CVODEIntegrator) ,intent(inout) :: this

        !This should be necessary... but it causes a segfault... no segfault when removed... no issues so far... 
       ! call FCVodeFree(this%cvode)
       ! call FN_VDestroy(this%sunVecY)

    end subroutine finalizeCVODEIntegrator 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule cvode_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
