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
module composite_integrator_class
    !! Houses container/controller for multiple integration routines

    use data_kinds                           ,only: rk, ik
    use runtime_constants                    ,only: debugging, assertions
    use god_objects                          ,only: object
    use assertion_utility                    ,only: assert, assertIdentical, assertPure
    use variable_container_class             ,only: VariableContainer
    use integrator_abstract_class            ,only: Integrator
    use timestep_controller_abstract_class   ,only: TimestepController
    use modeller_surrogate_class             ,only: ModellerSurrogate
    use manipulator_abstract_class           ,only: Manipulator
    use support_types                        ,only: IntArray, LogicalArray
    use mpi_controller_class                 ,only: CommunicationData

    implicit none
    private

    type ,public :: IntegratorContainer
        !! Container allowing for heterogeneous Integrator arrays
        class(Integrator) ,allocatable ,public :: entry
    end type

    type ,public :: IntegratorCallStep
        !! Object containing data pertaining to single integration step

        integer(ik)                                     ,public :: integratorIndex !! Index of integrator active during this step 
        type(IntArray) ,allocatable ,dimension(:)       ,public :: groupIndices !! Term groups evolved during step 
        integer(ik)    ,allocatable ,dimension(:)       ,public :: modelIndices !! Model indices of evolved groups
        real(rk)                                        ,public :: globalStepFraction !! Fraction of global timestep allocated to this integration stage
        logical                                         ,public :: allowTimeEvolution =.false. !! Set to true to allow this step to change "time" variable if it is present
        logical                                         ,public :: useInitialInput = .false.  !! Set to true to make this timestep use data from the start of the composite integration
        logical                                         ,public :: communicationNeeded = .false. !! Set to true to make integrator perform communication 
        logical                                         ,public :: nonTrivialUpdate = .false. !! Set to true to update terms during internal integrator iterations
        logical                                         ,public :: nonTrivialModelDataUpdate = .false. !! Set to true to update model data during internal integrator iterations
        type(CommunicationData)                         ,public :: commData !! Data used for communication if applicable
        type(LogicalArray) ,allocatable ,dimension(:)   ,public :: updatesOnInternalIterations !! True for those terms to be updated during internal integrator iterations
        logical       ,allocatable ,dimension(:)        ,public :: updatesOnInternalIterationsModelData !! True for those models whose data should be updated during internal integrator iterations
        
    end type 

    type ,public ,extends(Manipulator) :: CompositeIntegrator
        !! Composite integrator object allowing for different integration stages

        type(IntegratorContainer) ,allocatable ,dimension(:) ,private :: integrators !! Integrators contained in this composite integrator

        integer(ik)                                          ,private :: numIntegratorsAdded !! Counter keeping track of how many integrators have been added

        real(rk)                                             ,private :: currentTime !! Current time value (should agree with "time" variable if present)
        real(rk)                                             ,private :: globalTimestep !! Global timestep value
        real(rk)                                             ,private :: initialTimestep !! Initial timestep value
        real(rk)                                             ,private :: requestedTimestep !! Externally requested timestep, the global timestep will be set to the minimum between this value and the initialTimestep

        type(VariableContainer) ,allocatable                 ,private :: stepBuffer !! Variable container buffer used between steps

        type(IntegratorCallStep) ,allocatable ,dimension(:)  ,private :: integrationStage !! Integration stages 

        logical                                              ,private :: allIntegratorsAdded !! True when all integrators have been allocated

        class(TimestepController) ,allocatable               ,private :: dtController !! Optional timestep controller object at composite integrator level

        contains

        procedure ,public :: affect => integrateAll 
        procedure ,public :: getCurrentTime 
        procedure ,public :: addIntegrator
        procedure ,public :: addIntegrationStage 
        procedure ,public :: setTimestepController
        procedure ,public :: setRequestedTimestep


        procedure ,public :: init => initCompositeIntegrator

    end type CompositeIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initCompositeIntegrator(this,initialTime,initialTimestep,numIntegrators,dtController) 
        !! Composite integraotr initialization routine

        class(CompositeIntegrator)          ,intent(inout)  :: this
        real(rk)                            ,intent(in)     :: initialTime !! Time value before first step
        real(rk)                            ,intent(in)     :: initialTimestep !! Default timestep
        integer(ik)                         ,intent(in)     :: numIntegrators !! Number of integrators expected
        class(TimestepController) ,optional ,intent(in)     :: dtController !! Optional timestep controller

    end subroutine initCompositeIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCurrentTime(this) result(time)
        !! Getter for currentTime

        class(CompositeIntegrator)  ,intent(in) :: this
        real(rk)                                :: time

    end function getCurrentTime
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addIntegrator(this,integ) 
        !! Add Integrator object to composite integrator

        class(CompositeIntegrator)          ,intent(inout)  :: this
        class(Integrator)                   ,intent(in)     :: integ

    end subroutine addIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addIntegrationStage(this,integStage) 
        !! Add integration stage to composite integrator

        class(CompositeIntegrator)          ,intent(inout)  :: this
        type(IntegratorCallStep)            ,intent(in)     :: integStage

    end subroutine addIntegrationStage
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setTimestepController(this,controller) 
        !! Setter for dtController

        class(CompositeIntegrator)        ,intent(inout)  :: this
        class(TimestepController)         ,intent(in)     :: controller

    end subroutine setTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setRequestedTimestep(this,timestep) 
        !! Setter for requestedTimestep

        class(CompositeIntegrator)        ,intent(inout)  :: this
        real(rk)                          ,intent(in)     :: timestep 

    end subroutine setRequestedTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrateAll(this,manipulatedModeller,outputVars,inputVars) 
        !! Call all integrators based on the integration stages and global timestep. The global timestep is updated at the start if there is
        !! an allocated timestep controller.

        class(CompositeIntegrator)            ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller object used in callback
        class(VariableContainer)              ,intent(inout) :: outputVars !! Container for integration output
        class(VariableContainer)              ,intent(in)    :: inputVars !! Integration input variables

    end subroutine integrateAll
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module composite_integrator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
