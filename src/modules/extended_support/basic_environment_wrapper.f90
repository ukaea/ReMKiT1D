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
module basic_environment_wrapper
    !! author: Stefan Mijin 
    !!
    !! Houses wrapper containing non-modeller objects for streamlined initialization and passing

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use grid_class                  ,only: Grid
    use geometry_class              ,only: Geometry 
    use mpi_controller_class        ,only: MPIController 
    use v_space_class               ,only: VSpace 
    use partition_class             ,only: Partition 
    use indexing_class              ,only: Indexing 
    use json_controller_class       ,only: JSONController
    use petsc_controller_class      ,only: PETScController
    use variable_list_class         ,only: VariableList
    use variable_container_class    ,only: VariableContainer
    use HDF5_controller_class       ,only: HDF5Controller
    use textbook_class              ,only: Textbook
    use species_list_class          ,only: SpeciesList
    use signal_collection_class     ,only: SignalCollection
    use status_printing 
    use initialization_support
    use support_types

    implicit none
    private

    type ,public ,extends(Object) :: EnvironmentWrapper
        !! EnvironmentWrapper object storing publicly accessible grid, variable, and controller data

        type(Grid)                         ,public :: gridObj 
        type(Geometry)                     ,public :: geometryObj 
        type(Partition)                    ,public :: partitionObj 
        type(Indexing)                     ,public :: indexingObj 
        type(VSpace)                       ,public :: vSpaceObj

        type(MPIController)                ,public :: mpiCont 
        type(JSONController)  ,pointer     ,public :: jsonCont => null()
        type(HDF5Controller)               ,public :: hdf5Cont
        type(PETScController) ,allocatable ,public :: petscCont
        type(Textbook)                     ,public :: textbookObj
        type(SpeciesList)                  ,public :: speciesListObj
        type(SignalCollection)             ,public :: signalCollectionObj

        type(VariableContainer)            ,public :: externalVars
        type(VariableList)                 ,public :: allVars

        contains

        procedure ,public :: init => initEnvironmentFromJSON
        procedure ,public :: finishInit 

        final :: finalizeEnv

    end type EnvironmentWrapper
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initEnvironmentFromJSON(this,textbookObj,normObj)
        !! Initialize environment using data from config.json and textbook object for initializing external variables. Also initializes
        !! the status_printing module. If the textbook is not supplied the envirnment is not considered defined and finishInit must be 
        !! called before the environment is passed around.

        class(EnvironmentWrapper)  ,intent(inout) :: this
        type(Textbook) ,optional   ,intent(in)    :: textbookObj !! Textbook object used to retrieve derivation rules for variable container 
        class(Normalization) ,optional ,intent(in):: normObj !! Normalization object used to retrieve grid normalization. If not present grid is assumed normalized

    end subroutine initEnvironmentFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine finishInit(this,textbookObj)
        !! Finalize the initialization of the environment by setting the textbook object and performing initializations that require it
        
        class(EnvironmentWrapper)  ,intent(inout) :: this
        type(Textbook)             ,intent(in)    :: textbookObj !! Textbook object used to retrieve derivation rules for variable container 

    end subroutine finishInit
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeEnv(this) 

        type(EnvironmentWrapper) ,intent(inout) :: this

    end subroutine finalizeEnv 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module basic_environment_wrapper
!-----------------------------------------------------------------------------------------------------------------------------------
 