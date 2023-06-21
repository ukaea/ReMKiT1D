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
module integrator_abstract_class
    !! author: Stefan Mijin 
    !! 
    !! ouses abstract Integrator object, a ,anipulator with a timestep size and term group component and corresponding getters/setters

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use manipulator_abstract_class            ,only: Manipulator
    use variable_container_class              ,only: VariableContainer
    use timestep_controller_abstract_class    ,only: TimestepController
    use support_types                         ,only: IntArray ,LogicalArray
    use mpi_controller_class                  ,only: CommunicationData

    implicit none
    private

    type ,public ,extends(Manipulator), abstract :: Integrator
        !! Abstract Integrator object

        real(rk)                                  ,private :: dt !! The timestep this Integrator needs to take
        type(IntArray) ,allocatable ,dimension(:) ,private :: termGroups !! Term groups for each model 
                                                                         !! this Integrator is responsible for
        integer(ik)    ,allocatable ,dimension(:) ,private :: modelIndices !! Model indices this Integrator is responsible for

        logical ,private :: evolvesTimeVar !! If true and a variable named "time" exists in passed variable container the Integrator
                                           !! will be allowed to evolve it

        class(TimestepController) ,allocatable ,private :: dtController !! Optional timestep controller that computes the timestep
                                                                        !! based on input variables

        logical :: communicationNeeded !! True if this Integrator requires MPI communication during evolution

        type(CommunicationData) :: commData !! Communication data object for this Integrator

        type(LogicalArray) ,allocatable ,dimension(:) ,private :: updateOnInternalIteration !! Array that conforms to termGroups 
        !! and is true when given group should be updated on each internal iteration, otherwise they are only to be updated 
        !! every timestep

        logical :: nonTrivialUpdate !! True if this integrator should perform non-trivial updates of terms

        logical ,allocatable ,dimension(:) ,private :: updateModelDataOnInternalIteration !! Array that conforms with modelIndices and
        !! is true when the given model has model data that should be updated on every internal iteration

        logical :: nonTrivialModelDataUpdate !! True if this integrator should perform non-trivial updates of model data

        contains

        procedure ,public :: setTimestep
        procedure ,public :: getTimestep 

        procedure ,public :: setTermGroups
        procedure ,public :: getTermGroups

        procedure ,public :: setModelIndices
        procedure ,public :: getModelIndices

        procedure ,public :: setTimeEvolving
        procedure ,public :: isTimeEvolving

        procedure ,public :: setTimestepController
        procedure ,public :: hasTimestepController
        procedure ,public :: getTimestepFromController

        procedure ,public :: isCommunicationNeeded
        procedure ,public :: setCommunicationNeeded 

        procedure ,public :: getCommunicationData
        procedure ,public :: setCommunicationData

        procedure ,public :: getUpdateRules
        procedure ,public :: setUpdateRules

        procedure ,public :: getModelDataUpdateRules
        procedure ,public :: setModelDataUpdateRules

        procedure ,public :: hasNonTrivialUpdate
        procedure ,public :: setNonTrivialUpdate

        procedure ,public :: hasNonTrivialModelDataUpdate
        procedure ,public :: setNonTrivialModelDataUpdate

    end type Integrator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setTimestep(this,timestep)
        !! Setter for dt

        class(Integrator)              ,intent(inout)  :: this
        real(rk)                       ,intent(in)     :: timestep

    end subroutine setTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTimestep(this) result(timestep)
        !! Getter for dt

        class(Integrator)    ,intent(in) :: this
        real(rk)                         :: timestep

    end function getTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setTermGroups(this,groups)
        !! Setter for termGroups

        class(Integrator)               ,intent(inout)  :: this
        type(IntArray)   ,dimension(:)  ,intent(in)     :: groups

    end subroutine setTermGroups
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTermGroups(this) result(groups)
        !! Getter for termGroups

        class(Integrator)                     ,   intent(in) :: this
        type(IntArray) ,allocatable ,dimension(:)            :: groups

    end function getTermGroups
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setModelIndices(this,indices)
        !! Setter for modelIndices

        class(Integrator)            ,intent(inout)  :: this
        integer(ik)   ,dimension(:)  ,intent(in)     :: indices

    end subroutine setModelIndices
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getModelIndices(this) result(indices)
        !! Getter for modelIndices

        class(Integrator)                     ,intent(in) :: this
        integer(ik) ,allocatable ,dimension(:)            :: indices

    end function getModelIndices
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setTimeEvolving(this,evo)
        !! Setter for evolvesTimeVar

        class(Integrator)          ,intent(inout)  :: this
        logical                    ,intent(in)     :: evo

    end subroutine setTimeEvolving
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isTimeEvolving(this) result(evo)
        !! Check if this Integrator is allowed to evolve a time variable if present

        class(Integrator) ,intent(in) :: this
        logical                       :: evo

    end function isTimeEvolving
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setTimestepController(this,controller)
        !! Setter for dtController

        class(Integrator)          ,intent(inout)  :: this
        class(TimestepController)  ,intent(in)     :: controller

    end subroutine setTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function hasTimestepController(this) result(hasController)
        !! Check if this Integrator has an allocated timestep controller

        class(Integrator) ,intent(in) :: this
        logical                       :: hasController

    end function hasTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
    module function getTimestepFromController(this,inputVars) result(timestep)
        !! Get the individual timestep size if Integrator has a timestep controller

        class(Integrator)                     ,intent(inout) :: this 
        class(VariableContainer)              ,intent(in)    :: inputVars
        real(rk)                                             :: timestep
        
    end function getTimestepFromController
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setCommunicationNeeded(this,commNeeded)
        !! Setter for communicationNeeded

        class(Integrator)          ,intent(inout)  :: this
        logical                    ,intent(in)     :: commNeeded

    end subroutine setCommunicationNeeded
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isCommunicationNeeded(this) result(commNeeded)
        !! Check whether this Integrator requires MPI communication

        class(Integrator) ,intent(in) :: this
        logical                       :: commNeeded

    end function isCommunicationNeeded
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setCommunicationData(this,commData)
        !! Setter for commData

        class(Integrator)            ,intent(inout)  :: this
        type(CommunicationData)      ,intent(in)     :: commData

    end subroutine setCommunicationData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCommunicationData(this) result(commData)
        !! Getter for commData

        class(Integrator)                     ,intent(in) :: this
        type(CommunicationData) ,allocatable              :: commData

    end function getCommunicationData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setUpdateRules(this,updateRules)
        !! Setter for updateOnInternalIteration

        class(Integrator)                ,intent(inout)  :: this
        type(LogicalArray) ,dimension(:) ,intent(in)     :: updateRules

    end subroutine setUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getUpdateRules(this) result(updateRules)
        !! Getter for updateOnInternalIteration

        class(Integrator)                            ,intent(in) :: this
        type(LogicalArray) ,dimension(:) ,allocatable            :: updateRules

    end function getUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNonTrivialUpdate(this,nonTrivialUpdate)
        !! Setter for nonTrivialUpdate

        class(Integrator)          ,intent(inout)  :: this
        logical                    ,intent(in)     :: nonTrivialUpdate

    end subroutine setNonTrivialUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function hasNonTrivialUpdate(this) result(nonTrivialUpdate)
        !! Getter for nonTrivialUpdate

        class(Integrator) ,intent(in) :: this
        logical                       :: nonTrivialUpdate

    end function hasNonTrivialUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setModelDataUpdateRules(this,updateRules)
        !! Setter for updateModelDataOnInternalIteration

        class(Integrator)     ,intent(inout)  :: this
        logical ,dimension(:) ,intent(in)     :: updateRules

    end subroutine setModelDataUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getModelDataUpdateRules(this) result(updateRules)
        !! Getter for updateModelDataOnInternalIteration

        class(Integrator)                 ,intent(in) :: this
        logical  ,dimension(:) ,allocatable           :: updateRules

    end function getModelDataUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNonTrivialModelDataUpdate(this,nonTrivialUpdate)
        !! Setter for nonTrivialModelDataUpdate

        class(Integrator)          ,intent(inout)  :: this
        logical                    ,intent(in)     :: nonTrivialUpdate

    end subroutine setNonTrivialModelDataUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function hasNonTrivialModelDataUpdate(this) result(nonTrivialUpdate)
        !! Getter for nonTrivialModelDataUpdate
    
        class(Integrator) ,intent(in) :: this
        logical                       :: nonTrivialUpdate

    end function hasNonTrivialModelDataUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module integrator_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 