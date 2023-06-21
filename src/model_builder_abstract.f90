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
module model_builder_abstract_class
    !! author: Stefan Mijin
    !!
    !! Houses abstract Model builder object, used to construct and add Model objects to Modeller objects

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use support_types                         
    use mpi_controller_class                  ,only: MPIController
    use json_controller_class                 ,only: JSONController
    use modeller_class                        ,only: Modeller

    implicit none
    private

    type ,public ,extends(object), abstract :: ModelBuilder
        !! Abstract ModelBuilder object, housing basic functionality for the construction of Model objects and passing those objects
        !! to their final Modeller destination 

        type(NamedScalarContainer)          ,private :: scalarParams !! Named scalar parameters for reading/writing JSON config files and use in Model construction
        type(NamedArrayContainer)           ,private :: arrayParams !! Named array parameters for reading/writing JSON config files and use in Model construction

        integer(ik)                         ,private :: numImplicitGroups !! Number of implicit Term groups contained in the Model to be built
        integer(ik)                         ,private :: numGeneralGroups !! Number of general Term groups contained in th Model to be built

        contains

        procedure ,public :: loadParams
        procedure ,public :: outputUsedParams 

        procedure ,public :: getScalarParams
        procedure ,public :: getArrayParams
        
        procedure ,public :: setScalarParams
        procedure ,public :: setArrayParams

        procedure(addModel) ,deferred :: addModelToModeller

    end type ModelBuilder
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        subroutine addModel(this,modellerObj) 
            !! Abstract routine for adding Models to a Modeller object after building them
            import :: ModelBuilder ,Modeller

            class(ModelBuilder)          ,intent(inout) :: this 
            class(Modeller)              ,intent(inout) :: modellerObj
            !! Modeller object that will house the constructed Model after it is built by the builder

        end subroutine addModel
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadParams(this,jsonCont,mpiCont) 
        !! Load parameters from "./config.json"

        class(ModelBuilder)   ,intent(inout)  :: this
        type(JSONController)  ,intent(inout)  :: jsonCont
        !! JSONController object responsible for reading the config file
        type(MPIController)   ,intent(inout)  :: mpiCont
        !! MPIController object to be used with JSON IO
    
    end subroutine loadParams
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine outputUsedParams(this,jsonCont,mpiCont) 
        !! Output used parameters to "./used_config.json"

        class(ModelBuilder)   ,intent(inout)  :: this
        type(JSONController)  ,intent(inout)  :: jsonCont
        !! JSONController object responsible for writing the config file
        type(MPIController)   ,intent(inout)  :: mpiCont
        !! MPIController object to be used with JSON IO
    
    end subroutine outputUsedParams
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getScalarParams(this) result(params)
        !! Getter for ModelBuilder scalarParams

        class(ModelBuilder)   ,intent(in)  :: this
        type(NamedScalarContainer)         :: params
    
    end function getScalarParams
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getArrayParams(this) result(params)
        !! Getter for ModelBuilder arrayParams

        class(ModelBuilder)  ,intent(in)  :: this
        type(NamedArrayContainer)         :: params
    
    end function getArrayParams
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setScalarParams(this,params) 
        !! Setter for ModelBuilder scalarParams

        class(ModelBuilder)        ,intent(inout)  :: this
        type(NamedScalarContainer) ,intent(in)     :: params
    
    end subroutine setScalarParams
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setArrayParams(this,params) 
        !! Setter for ModelBuilder scalarParams

        class(ModelBuilder)        ,intent(inout)  :: this
        type(NamedArrayContainer)  ,intent(in)     :: params
    
    end subroutine setArrayParams
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module model_builder_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 