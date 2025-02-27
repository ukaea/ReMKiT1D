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
submodule (composite_integrator_class) composite_integrator_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains procedures associated with the composite integrator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initCompositeIntegrator(this,initialTime,initialTimestep,numIntegrators,dtController) 
    !! Composite integraotr initialization routine

    class(CompositeIntegrator)          ,intent(inout)  :: this
    real(rk)                            ,intent(in)     :: initialTime !! Time value before first step
    real(rk)                            ,intent(in)     :: initialTimestep !! Default timestep
    integer(ik)                         ,intent(in)     :: numIntegrators !! Number of integrators expected
    class(TimestepController) ,optional ,intent(in)     :: dtController !! Optional timestep controller

    if (assertions .and. present(dtController)) call assertPure(dtController%isDefined(),"Undefined timestep controller passed to&
    &composite integrator constructor")
    
    this%globalTimestep = initialTimestep
    this%initialTimestep = initialTimestep
    this%currentTime = initialTime 

    allocate(this%integrators(numIntegrators))

    this%numIntegratorsAdded = 0
    this%allIntegratorsAdded = .false.

    allocate(this%integrationStage(0))

    if (present(dtController)) call this%setTimestepController(dtController)
    call this%setRequestedTimestep(initialTimestep)

    call this%makeDefined()

end subroutine initCompositeIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCurrentTime(this) result(time)
    !! Getter for currentTime

    class(CompositeIntegrator)  ,intent(in) :: this
    real(rk)                                :: time

    if (assertions) call assertPure(this%isDefined(),"Attempted to get current time from undefined composite integrator")

    time = this%currentTime

end function getCurrentTime
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addIntegrator(this,integ) 
    !! Add Integrator object to composite integrator

    class(CompositeIntegrator)          ,intent(inout)  :: this
    class(Integrator)                   ,intent(in)     :: integ

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add integrator component to undefined composite integrator")
        call assertPure(integ%isDefined(),"Attempted to add undefined integrator component to composite integrator")
        call assertPure(.not. this%allIntegratorsAdded,"Attempted to add integrator component to composite integrator with no free&
        & integrator slots")
    end if

    this%numIntegratorsAdded = this%numIntegratorsAdded + 1

    allocate(this%integrators(this%numIntegratorsAdded)%entry,source=integ)

    if (this%numIntegratorsAdded == size(this%integrators)) this%allIntegratorsAdded = .true.

end subroutine addIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addIntegrationStage(this,integStage) 
    !! Add integration stage to composite integrator

    class(CompositeIntegrator)          ,intent(inout)  :: this
    type(IntegratorCallStep)            ,intent(in)     :: integStage

    if (assertions) call assertPure(this%isDefined(),"Attempted to add integration stage to undefined composite integrator")

    this%integrationStage = [this%integrationStage,integStage]

end subroutine addIntegrationStage
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setTimestepController(this,controller) 
    !! Setter for dtController

    class(CompositeIntegrator)        ,intent(inout)  :: this
    class(TimestepController)         ,intent(in)     :: controller

    if(allocated(this%dtController)) deallocate(this%dtController)
    allocate(this%dtController,source=controller)

end subroutine setTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setRequestedTimestep(this,timestep) 
    !! Setter for requestedTimestep

    class(CompositeIntegrator)        ,intent(inout)  :: this
    real(rk)                          ,intent(in)     :: timestep 

    this%requestedTimestep = timestep

end subroutine setRequestedTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine resetRequestedTimestep(this) 
    !! Resets the requested timestep to the initial timestep

    class(CompositeIntegrator)        ,intent(inout)  :: this

    this%requestedTimestep = this%initialTimestep

end subroutine resetRequestedTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine integrateAll(this,manipulatedModeller,outputVars,inputVars) 
    !! Call all integrators based on the integration stages and global timestep. The global timestep is updated at the start if there is
    !! an allocated timestep controller.

    class(CompositeIntegrator)            ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller object used in callback
    class(VariableContainer)              ,intent(inout) :: outputVars !! Container for integration output
    class(VariableContainer)              ,intent(in)    :: inputVars !! Integration input variables

    integer(ik) :: i ,integratorIndex ,numStages

    if (assertions) then 

        call assertPure(this%isDefined(),"Attempted to integrate modeller using undefined composite integrator")
        call assertPure(manipulatedModeller%isDefined(),"Attempted to integrate undefined modeller using composite integrator")
        call assertPure(outputVars%isDefined(),"Undefined outputVars passed to integrate all routine of composite integrator")
        call assertPure(inputVars%isDefined(),"Undefined inputVars passed to integrate all routine of composite integrator")
        call assertPure(size(this%integrationStage) >= 1,"At least one integration stage must be added to a composite integrator")
        call assertPure(this%allIntegratorsAdded,"Attempted to integrate modeller when not all integrators have been added to &
        &composite integrator object")

    end if

    if (size(this%integrationStage) > 1) then
        if (.not. allocated(this%stepBuffer)) allocate(this%stepBuffer,source=inputVars)
    end if

    this%globalTimestep = this%initialTimestep
    if (allocated(this%dtController)) this%globalTimestep = this%dtController%evaluateTimestep(inputVars,this%globalTimestep)
    this%globalTimestep = min(this%globalTimestep,this%requestedTimestep)

    numStages = size(this%integrationStage)

    if (inputVars%isVarNameRegistered("time")) this%currentTime = inputVars%variables(outputVars%getVarIndex("time"))%entry(1)
    do i = 1,numStages
        integratorIndex = this%integrationStage(i)%integratorIndex

        call this%integrators(integratorIndex)%entry%setTimeEvolving(this%integrationStage(i)%allowTimeEvolution)
        call this%integrators(integratorIndex)%entry%setTimestep(this%globalTimestep*this%integrationStage(i)%globalStepFraction)
        call this%integrators(integratorIndex)%entry%setTermGroups(this%integrationStage(i)%groupIndices)
        call this%integrators(integratorIndex)%entry%setModelIndices(this%integrationStage(i)%modelIndices)
        call this%integrators(integratorIndex)%entry%setCommunicationNeeded(this%integrationStage(i)%communicationNeeded)
        if (this%integrationStage(i)%communicationNeeded) &
        call this%integrators(integratorIndex)%entry%setCommunicationData(this%integrationStage(i)%commData)
        call this%integrators(integratorIndex)%entry%setNonTrivialUpdate(this%integrationStage(i)%nonTrivialUpdate)
        if (this%integrationStage(i)%nonTrivialUpdate) then
            call this%integrators(integratorIndex)%entry%setUpdateRules(this%integrationStage(i)%updatesOnInternalIterations)
        end if
        call this%integrators(integratorIndex)&
        %entry%setNonTrivialModelDataUpdate(this%integrationStage(i)%nonTrivialModelDataUpdate)
        if (this%integrationStage(i)%nonTrivialModelDataUpdate) then
            call this%integrators(integratorIndex)&
            %entry%setModelDataUpdateRules(this%integrationStage(i)%updatesOnInternalIterationsModelData)
        end if
       
        if (this%integrationStage(i)%useInitialInput .or. (i == 1)) then
            if (i == numStages) then 
                call this%integrators(integratorIndex)%entry%affect(manipulatedModeller,outputVars,inputVars)
            else
                call this%integrators(integratorIndex)%entry%affect(manipulatedModeller,this%stepBuffer,inputVars)
            end if
        else
            if (i == numStages) then 
                call this%integrators(integratorIndex)%entry%affect(manipulatedModeller,outputVars,this%stepBuffer)
            else
                call this%integrators(integratorIndex)%entry%affect(manipulatedModeller,this%stepBuffer,this%stepBuffer)
            end if
        end if
        
    end do

    this%currentTime = this%currentTime + this%globalTimestep 
    if (outputVars%isVarNameRegistered("time")) outputVars%variables(outputVars%getVarIndex("time"))%entry(1) = this%currentTime

end subroutine integrateAll
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule composite_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
