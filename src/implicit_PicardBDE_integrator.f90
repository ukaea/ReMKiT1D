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
module implicit_PicardBDE_integrator_class
    !! author: Stefan Mijin 
    !! 
    !! Contains an implicit Picard Backwards Euler integrator class 

    use data_kinds                         ,only: rk ,ik
    use god_objects                        ,only: Object
    use runtime_constants                  ,only: debugging, assertions
    use assertion_utility                  ,only: assert, assertIdentical, assertPure
    use modeller_surrogate_class           ,only: ModellerSurrogate
    use variable_container_class           ,only: VariableContainer
    use modeller_class                     ,only: Modeller 
    use integrator_abstract_class          ,only: Integrator
    use support_types                      ,only: IntArray ,RealArray ,LogicalArray
    use mpi_controller_class               ,only: CommunicationData
    use timestep_controller_abstract_class ,only: TimestepController
    use indexing_class                     ,only: Indexing
    use status_printing                    

    implicit none
    private

    type ,public :: InternalControllerOptions 

        integer(ik) :: currentNumSubsteps = 1
        integer(ik) :: stepMultiplier = 2
        integer(ik) :: stepDecrament = 1
        
        integer(ik) :: minNonlinIters = 5

        integer(ik) :: maxRestarts = 3
        integer(ik) :: restartCount = 0

    end type InternalControllerOptions

    type ,public ,extends(Integrator) :: PicardBDEIntegrator
        !! Implicit Backwards Euler integrator with Picard iterations for handling non-linearities

        real(rk) ,allocatable ,dimension(:)          ,private :: implicitVectorOld !! Buffer for implicit vector from previous Picard iteration
        real(rk) ,allocatable ,dimension(:)          ,private :: implicitVectorNew !! Buffer for implicit vector from current Picard iteration
        integer(ik) ,allocatable ,dimension(:)       ,private :: convergenceTestVars !! Variable indices in VariableContainer used to determine convergence of Picard iterations
        real(rk)                                     ,private :: nonlinTol !! Picard iteration relative convergence tolerance
        real(rk)                                     ,private :: absTol !! Picard iteration absolute convergence tolerance in epsilon units
        integer(ik)                                  ,private :: maxIterations !! Maximum allowed number of Picard iterations
        type(VariableContainer) ,allocatable         ,private :: buffer !! VariableContainer buffer for passing to Modeller routines

        integer(ik)                                  ,private :: timesCalled !! Tracker for number of times called
        integer(ik)                                  ,private :: totNumIters !! Tracker for total number of iterations

        logical                                      ,private :: use2Norm !! True if the norm to be used is the 2norm instead of local inf-norm

        integer(ik)                                  ,private :: associatedPETScObjGroup !! PETSc obj group this solver should interact with (defaults to 1)

        logical                                      ,private :: internalStepControl !! True if using internal step control - will attempt restarts on failed solves and will ignore timestep controller

        type(InternalControllerOptions)              ,private :: internalControlOpts 

        character(:) ,allocatable                    ,private :: integratorName !! Integrator named used in printing 
        contains

        procedure ,public :: affect => integrateBDE

        procedure ,public :: init => initBDEIntegrator

        procedure ,public :: setNonlinTol
        procedure ,public :: getNonlinTol 

        procedure ,public :: setMaxIterations
        procedure ,public :: getMaxIterations

        procedure ,public :: setConvergenceIndices
        procedure ,public :: getConvergenceIndices

    end type PicardBDEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrateBDE(this,manipulatedModeller,outputVars,inputVars) 
        !! Integration routine for BDE integrator - implementation of abstract manipulate routine

        class(PicardBDEIntegrator)            ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    end subroutine integrateBDE
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

    end subroutine initBDEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNonlinTol(this,nonlinTol)
        !! Setter for nonlinTol value

        class(PicardBDEIntegrator)     ,intent(inout)  :: this
        real(rk)                       ,intent(in)     :: nonlinTol !! Picard iteration tolerance

    end subroutine setNonlinTol
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNonlinTol(this) result(nonlinTol)
        !! Getter for nonlinTol value

        class(PicardBDEIntegrator) ,intent(in) :: this
        real(rk)                               :: nonlinTol

    end function getNonlinTol
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setMaxIterations(this,maxIter)
        !! Setter for maxIterations value

        class(PicardBDEIntegrator)     ,intent(inout)  :: this
        integer(ik)                    ,intent(in)     :: maxIter !! Maximum number of Picard iterations

    end subroutine setMaxIterations
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getMaxIterations(this) result(maxIter)
        !! Geter for maxIterations value

        class(PicardBDEIntegrator) ,intent(in) :: this
        integer(ik)                            :: maxIter

    end function getMaxIterations
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setConvergenceIndices(this,convVars)
        !! Setter for convergenceTestVars array

        class(PicardBDEIntegrator)        ,intent(inout)  :: this
        integer(ik) ,dimension(:)         ,intent(in)     :: convVars ! Variable indices in VariableContainer used to determine convergence of Picard iterations

    end subroutine setConvergenceIndices
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getConvergenceIndices(this) result(convVars)
        !! Getter for convergenceTestVars array

        class(PicardBDEIntegrator)    ,intent(in) :: this
        integer(ik) ,dimension(:) ,allocatable    :: convVars

    end function getConvergenceIndices
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module implicit_PicardBDE_integrator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 