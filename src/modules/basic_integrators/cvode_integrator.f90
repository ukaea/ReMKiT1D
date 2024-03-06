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
module cvode_integrator_class
    !! author: Stefan Mijin
    !! 
    !! Contains an Sundials CVODE integrator class 

    use data_kinds                         ,only: rk ,ik
    use god_objects                        ,only: Object
    use runtime_constants                  ,only: debugging, assertions, assertionLvl
    use assertion_utility                  ,only: assert, assertIdentical, assertPure
    use modeller_surrogate_class           ,only: ModellerSurrogate
    use variable_container_class           ,only: VariableContainer
    use modeller_class                     ,only: Modeller 
    use integrator_abstract_class          ,only: Integrator
    use support_types                      ,only: IntArray ,RealArray ,LogicalArray
    use mpi_controller_class               ,only: CommunicationData,MPIController
    use timestep_controller_abstract_class ,only: TimestepController
    use, intrinsic :: iso_c_binding

    use fcvode_mod                ! Access CVode
    implicit none
    private

    type ,public ,extends(Integrator) :: CVODEIntegrator

        type(c_ptr), private :: sunctx !! Sundial context 
        real(rk) ,allocatable ,dimension(:) ,private :: evolvedVars

        contains

        procedure ,public :: affect => integrateCVODE

        procedure ,public :: init => initCVODEIntegrator

    end type CVODEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrateCVODE(this,manipulatedModeller,outputVars,inputVars) 
        !! Implementation of abstract manipulate routine for the case of Runge-Kutta integrator

        class(CVODEIntegrator)                ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    end subroutine integrateCVODE
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCVODEIntegrator(this,mpiCont,modelList,termGroups,evolvesTimeVar,dtController,initialTimestep)

        class(CVODEIntegrator)                    ,intent(inout) :: this
        type(MPIController)                       ,intent(in)    :: mpiCont !! 
        integer(ik)    ,optional    ,dimension(:) ,intent(in)    :: modelList !! List of models this integrator will be responsible for
        type(IntArray) ,optional    ,dimension(:) ,intent(in)    :: termGroups !! Term groups this integrator is responsible for - should conform with modelList
        logical        ,optional                  ,intent(in)    :: evolvesTimeVar !! Set to true if this integrator is allowed to change the "time" varible (if present in passed variable container)
        class(TimestepController) ,optional       ,intent(in)    :: dtController !! User-supplied timestep controller object
        real(rk)                  ,optional       ,intent(in)    :: initialTimestep !! Default timestep

    end subroutine initCVODEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module cvode_integrator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
