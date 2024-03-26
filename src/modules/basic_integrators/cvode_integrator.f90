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
    use support_types                      ,only: IntArray ,RealArray ,LogicalArray, StringArray
    use mpi_controller_class               ,only: CommunicationData,MPIController
    use timestep_controller_abstract_class ,only: TimestepController
    use support_functions
    use status_printing                    
    use, intrinsic :: iso_c_binding
    use mpi_f08

    use fcvode_mod                ! Access CVode
    use fsundials_core_mod
    use fnvector_parallel_mod 
    use fsunlinsol_spgmr_mod
    implicit none
    private

    type ,public :: CVODEOptions 

        integer(ik) :: bbdmukeep !! BBD preconditioner size of kept upper diagonal band 
        integer(ik) :: bbdmlkeep !! BBD preconditioner size of kept lower diagonal band 
        integer(ik) :: bbdmldq   !! BBD preconditioner size of lower differentiation stencil band
        integer(ik) :: bbdmudq    !! BBD preconditioner size of upper differentitation stencil band

        real(rk) :: reltol !! Linear solver relative tolerance 
        real(rk) :: abstol !! Linear solver absolute tolerance 

        integer(ik) :: maxRestarts !! Maximum number of GMRES restarts

    end type CVODEOptions

    type ,public ,extends(Integrator) :: CVODEIntegrator

        type(c_ptr)                                  ,private :: sunctx !! Sundials context 
        type(SUNLinearSolver) ,pointer               ,private :: sunls !! Sundials linear solver pointer 
        type(N_Vector)        ,pointer               ,private :: sunVecY, sunVecYDot !! Solution and RHS sundials vector objects
        type(c_ptr)                                  ,private :: cvode !! CVODE solver context
        type(SUNMatrix)       ,pointer               ,private :: sunmat !! Sundials matrix (will be nulled)
        real(c_double)        ,pointer ,dimension(:) ,private :: rhsVec,yVec !! Pointers used for accessing the underlying N_Vectors 

        type(VariableContainer) ,allocatable               ,private :: bufferRHS ,bufferY !! ReMKiT1D variable containers corresponding to the RHS and solutions - enable derivation and communication calls
        real(rk)                ,allocatable ,dimension(:) ,private :: copyBufferVals !! Buffer for moving data between the C pointers and the ReMKiT1D variable containers 
        type(StringArray)       ,allocatable ,dimension(:) ,private :: evolvedVars !! List of variables evolved by this solver

        type(MPI_Comm) ,private :: mpiComm !! MPI communicator object needed by some CVODE routines
        integer(ik)    ,private :: numProcs !! Number of processes needed for initialising the number of variables/equations in CVODE 

        type(CVODEOptions) ,private :: options
        
        character(:) ,allocatable                    ,private :: integratorName !! Integrator named used in printing 

        contains

        procedure ,public :: affect => integrateCVODE

        procedure ,public :: init => initCVODEIntegrator

        final :: finalizeCVODEIntegrator

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

    end subroutine initCVODEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    impure elemental module subroutine finalizeCVODEIntegrator(this) 

        type(CVODEIntegrator) ,intent(inout) :: this

    end subroutine finalizeCVODEIntegrator 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module cvode_integrator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
