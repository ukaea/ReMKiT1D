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
submodule (explicit_rk_integrator_class) explicit_rk_integrator_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the explicit RK integrator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initRKIntegrator(this,modelList,termGroups,order,tableau,evolvesTimeVar,dtController,initialTimestep)
    !! Constructs an RK integrator with initial evaluated model indices and correspoding term groups if those
    !! are present. If neither order or tableau are provided defaults to forward Euler time-stepping. Possible order values are 1,2,3,4
    !! where 2 is the midpoint method, 3 is SSPRK3, and 4 is the standard RK4 method. 

    implicit none 

    class(ExplicitRKIntegrator)               ,intent(inout) :: this
    integer(ik)    ,optional    ,dimension(:) ,intent(in)    :: modelList !! List of models this integrator will be responsible for
    type(IntArray) ,optional    ,dimension(:) ,intent(in)    :: termGroups !! Term groups this integrator is responsible for - should conform with modelList
    integer(ik)    ,optional                  ,intent(in)    :: order !! Runge-Kutta order
    type(BTableau) ,optional                  ,intent(in)    :: tableau !! User-defined Butcher tableau
    logical        ,optional                  ,intent(in)    :: evolvesTimeVar !! Set to true if this integrator is allowed to change the "time" varible (if present in passed variable container)
    class(TimestepController) ,optional       ,intent(in)    :: dtController !! User-supplied timestep controller object
    real(rk)                  ,optional       ,intent(in)    :: initialTimestep !! Default timestep

    integer(ik)                                         :: i,j

    if (assertions) then
        
        if (present(termGroups)) call assert(present(modelList),"Term groups object passed to RK integrator constructor without&
                                                               & model list")
        if (present(modelList) .and. present(termGroups)) call assert(size(modelList)==size(termGroups),"If both model list and &
        & term groups arguments are passed to initRKIntegrator both must have the same length")
        if (present(order)) call assert(.not. present(tableau),"Both tableau and order arguments cannot be present in RK integrator&
                                       & constructor")
        if (present(order)) call assert(any([1,2,3,4] == order) ,"Supported default orders for explicit RK integrator & 
                                            &are 1, 2, 3, and 4 - for non-default orders and tableaus supply user tableau")
        if (present(tableau)) then 
            call assertIdentical([size(tableau%b),size(tableau%c),size(tableau%a)+1],"Tableau data passed to RK integrator &
                                &constructor does not conform in size")
            do i = 2, size(tableau%b)   
                call assert(size(tableau%a(i-1)%entry)==i-1,"Tableau A matrix not lower triangular")
            end do
        end if
        if (present(dtController)) call assert(dtController%isDefined(),"Undefined timestep controller passed to explicit RK&
        & integrator constructor")
    end if

    if (present(order)) then
        select case (order)
        case(1)

            allocate(this%tableau%a(0))
            this%tableau%b = real([1],kind=rk)
            this%tableau%c = real([0],kind=rk)

        case(2)

            allocate(this%tableau%a(1))
            this%tableau%a(1)%entry = real([0.5d0],kind=rk)
            this%tableau%b = real([0,1],kind=rk)
            this%tableau%c = real([0d00,0.5d0],kind=rk)
        case(3)
            allocate(this%tableau%a(2))
            this%tableau%a(1)%entry = real([1d00],kind=rk)
            this%tableau%a(2)%entry = real([0.25d0,0.25d0],kind=rk)

            this%tableau%b = real([1.0d0/6.0d0,1.0d0/6.0d0,2.0d0/3.0d0],kind=rk)
            this%tableau%c = real([0d00,1d00,0.5d0],kind=rk)
        case(4)

            allocate(this%tableau%a(3))
            this%tableau%a(1)%entry = real([0.5d0],kind=rk)
            this%tableau%a(2)%entry = real([0d00,0.5d0],kind=rk)
            this%tableau%a(3)%entry = real([0,0,1],kind=rk)

            this%tableau%b = real([1.0d0/6.0d0,1.0d0/3.0d0,1.0d0/3.0d0,1.0d0/6.0d0],kind=rk)
            this%tableau%c = real([0d00,0.5d0,0.5d0,1d00],kind=rk)
        case default 
            error stop "Unsupported order passed to RK integrator constructor"
        end select
    end if

    if (present(tableau)) this%tableau = tableau 
    ! Default to Forward Euler
    if (.not. present(order) .and. .not. present(tableau)) then 
        allocate(this%tableau%a(0))
        this%tableau%b = real([1],kind=rk)
        this%tableau%c = real([0],kind=rk)
    end if

    !Find tableau non-zeros 
    allocate(this%tableauANonzeros(size(this%tableau%a)))

    do i = 1, size(this%tableau%a)
        allocate(this%tableauANonzeros(i)%entry(0))

        do j = 1, size(this%tableau%a(i)%entry)
            if (abs(this%tableau%a(i)%entry(j)) > 100*epsilon(this%tableau%a(i)%entry(j))) &
            this%tableauANonzeros(i)%entry = [this%tableauANonzeros(i)%entry,j]
        end do  

    end do

    if (present(evolvesTimeVar)) call this%setTimeEvolving(evolvesTimeVar) 
    if (present(modelList)) call this%setModelIndices(modelList)
    if (present(termGroups)) call this%setTermGroups(termGroups)
    if (present(initialTimestep)) call this%setTimestep(initialTimestep)
    if (present(dtController)) call this%setTimestepController(dtController)
    call this%setNonTrivialModelDataUpdate(.false.)
    call this%setNonTrivialUpdate(.false.)

    call this%makeDefined()

end subroutine initRKIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine integrateRK(this,manipulatedModeller,outputVars,inputVars) 
    !! Implementation of abstract manipulate routine for the case of Runge-Kutta integrator

    class(ExplicitRKIntegrator)           ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
    class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
    class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine


    type(IntArray) ,allocatable ,dimension(:)  :: termGroups 
    integer(ik)    ,allocatable ,dimension(:)  :: modelIndices 

    integer(ik)                         :: RKOrder ,numSteps ,i ,j ,k ,RKStage ,l ,nonzeroIndex
    real(rk) ,allocatable ,dimension(:) :: dt
    real(rk)                            :: fullTimestep, controllerTimestep

    type(IntArray)     ,allocatable ,dimension(:) :: varIndices
    type(LogicalArray) ,allocatable ,dimension(:) :: updateRules
    type(LogicalArray) ,allocatable ,dimension(:) :: mixedGroup
    integer(ik)                                   :: timeVarIndex
    logical            ,allocatable ,dimension(:) :: modelDataUpdateRules


    logical                 :: commNeeded
    logical                 :: nonTrivialUpdate
    logical                 :: nonTrivialModelDataUpdate
    logical                 :: timeEvolving
    type(CommunicationData) :: commData

    termGroups = this%getTermGroups()
    modelIndices = this%getModelIndices()
    nonTrivialUpdate = this%hasNonTrivialUpdate()
    if (nonTrivialUpdate) updateRules = this%getUpdateRules()
    nonTrivialModelDataUpdate = this%hasNonTrivialModelDataUpdate()
    if (nonTrivialModelDataUpdate) modelDataUpdateRules = this%getModelDataUpdateRules()
    
    if (assertions) then 
        call assert(this%isDefined(),"Integration routine called by undefined explicit RK integrator")
        call assert(manipulatedModeller%isDefined(),"Undefined Modeller object passed to explicit RK integration routine")
        call assert(outputVars%isDefined(),"outputVars passed to integrateRK not defined")
        call assert(inputVars%isDefined(),"inputVars passed to integrateRK not defined")
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

    select type (manipulatedModeller)
    type is (Modeller)
        
        RKOrder = size(this%tableau%b) 

        if (.not. allocated(this%buffer)) then 
            allocate(this%buffer,source=inputVars)
            allocate(this%intermediateValues(size(this%buffer%variables),RKOrder+1))
            do i = 1, RKOrder
                do j = 1, size(this%buffer%variables)
                    allocate(this%intermediateValues(j,i)%entry(size(this%buffer%variables(j)%entry)))
                end do
            end do
        end if

        allocate(mixedGroup(size(modelIndices)))
        allocate(varIndices(size(modelIndices)))
        do i = 1, size(modelIndices)
            allocate(varIndices(i)%entry(size(termGroups(i)%entry)))
            allocate(mixedGroup(i)%entry(size(termGroups(i)%entry)))
            do j = 1, size(termGroups(i)%entry)
                mixedGroup(i)%entry(j) = manipulatedModeller%isModelTermGroupMixed(termGroups(i)%entry(j),modelIndices(i))
                if (.not. mixedGroup(i)%entry(j)) then
                    varIndices(i)%entry(j) = inputVars%getVarIndex(manipulatedModeller&
                                                              %getEvolvedVarInTermGroup(termGroups(i)%entry(j),modelIndices(i)))
                    if (assertions) call assert(.not.inputVars%isStationary(manipulatedModeller&
                    %getEvolvedVarInTermGroup(termGroups(i)%entry(j),modelIndices(i))),&
                    "Explicit RK integrator detected stationary variable among evolved variables - this is unsupported")
                end if
            end do
        end do

        if (inputVars%isVarNameRegistered("time")) timeVarIndex = inputVars%getVarIndex("time")
        
        this%intermediateValues(:,RKOrder+1) = inputVars%variables
        do i = 1,numSteps 
            this%buffer%variables = this%intermediateValues(:,RKOrder+1)
            do j = 1,size(modelIndices)
                call manipulatedModeller%updateModelData(modelIndices(j),this%buffer)
                do k = 1,size(termGroups(j)%entry)
                    if (.not. mixedGroup(j)%entry(k)) then
                        call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                        call manipulatedModeller%calculateMatGroupValsInModel(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                        this%intermediateValues(varIndices(j)%entry(k),1)%entry = &
                        manipulatedModeller%evaluateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                    end if
                end do
            end do

            do RKStage = 2, RKOrder

                this%buffer%variables = this%intermediateValues(:,RKOrder+1)
                do l = 1,size(this%tableauANonzeros(RKStage-1)%entry)
                    nonzeroIndex = this%tableauANonzeros(RKStage-1)%entry(l)
                    do j = 1,size(modelIndices)
                        do k = 1,size(termGroups(j)%entry)
                            if (.not. mixedGroup(j)%entry(k)) then
                                    this%buffer%variables(varIndices(j)%entry(k))%entry = &
                                    this%buffer%variables(varIndices(j)%entry(k))%entry + &
                                    this%intermediateValues(varIndices(j)%entry(k),nonzeroIndex)%entry &
                                    * dt(i) * this%tableau%a(RKStage-1)%entry(nonzeroIndex)
                            end if
                        end do
                    end do
                end do
                
                if (this%buffer%isVarNameRegistered("time")) &
                    this%buffer%variables(timeVarIndex)%entry(1) = &
                    this%buffer%variables(timeVarIndex)%entry(1) + dt(i) * this%tableau%c(RKStage)

                call manipulatedModeller%callManipulator(0,this%buffer,this%buffer)
                if (commNeeded) then 
                    call manipulatedModeller%safeCommAndDeriv(commData,this%buffer,derivPriority=0)
                else
                    call this%buffer%calculateDerivedVars(derivPriority=0)
                end if 
                do j = 1,size(modelIndices)
                    if (nonTrivialModelDataUpdate) then 
                        if (modelDataUpdateRules(j)) &
                        call manipulatedModeller%updateModelData(modelIndices(j),this%buffer,updatePriority=0)
                    end if
                    do k = 1,size(termGroups(j)%entry)
                        if (.not. mixedGroup(j)%entry(k)) then
                            if (nonTrivialUpdate) then 
                                if (updateRules(j)%entry(k))&
                                call manipulatedModeller%updateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                            end if
                            call manipulatedModeller%calculateMatGroupValsInModel(modelIndices(j)&
                                                                                 ,termGroups(j)%entry(k)&
                                                                                 ,this%buffer)
                            this%intermediateValues(varIndices(j)%entry(k),RKStage)%entry = &
                            manipulatedModeller%evaluateModelTermGroup(modelIndices(j),termGroups(j)%entry(k),this%buffer)
                        end if
                    end do
                end do
               
            end do

            do RKStage = 1,RKOrder
                do j = 1,size(modelIndices)
                    do k = 1,size(termGroups(j)%entry)
                        if (.not. mixedGroup(j)%entry(k)) then
                            this%intermediateValues(varIndices(j)%entry(k),RKOrder+1)%entry = &
                            this%intermediateValues(varIndices(j)%entry(k),RKOrder+1)%entry + &
                            this%intermediateValues(varIndices(j)%entry(k),RKStage)%entry * dt(i) * this%tableau%b(RKStage)
                        end if
                    end do
                end do
            end do

            !Evolve time if allowed to and time is explicitly present
            if (timeEvolving .and. this%buffer%isVarNameRegistered("time")) &
                this%intermediateValues(timeVarIndex,RKOrder+1)%entry(1) = &
                this%intermediateValues(timeVarIndex,RKOrder+1)%entry(1) + dt(i)
            call manipulatedModeller%callManipulator(1,this%buffer,this%buffer)
            if (commNeeded) then 
                call manipulatedModeller%safeCommAndDeriv(commData,this%buffer)
            else
                call this%buffer%calculateDerivedVars()
            end if 
        end do

        outputVars%variables = this%intermediateValues(:,RKOrder+1)
        call manipulatedModeller%callManipulator(2,this%buffer,this%buffer)
        if (commNeeded) then 
            call manipulatedModeller%safeCommAndDeriv(commData,outputVars)
        else
            !A little bit overkill - could specialize to update only those variables that could have changed
            call outputVars%calculateDerivedVars()
        end if 
    class default
        error stop "Unsupported surrogate passed to explicit RK integrator"
    end select

end subroutine integrateRK
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule explicit_rk_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
