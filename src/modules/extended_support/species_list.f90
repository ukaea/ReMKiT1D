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
module species_list_class
    !! author: Stefan Mijin
    !! 
    !! Houses species list class containing species initialized from config.json and accessed either via their ID or name

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions, assertionLvl
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use species_class                         ,only: Species
    use mpi_controller_class                  ,only: MPIController
    use json_controller_class                 ,only: JSONController
    use key_names
    use support_types                
    

    implicit none
    private

    type ,public ,extends(Object) :: SpeciesList
        !! Species list object providing centralized access and initialization of species objects

        type(Species)          ,allocatable ,dimension(:) ,private :: speciesEntries
        integer(ik)            ,allocatable ,dimension(:) ,private :: speciesIDs
        type(StringArray)      ,allocatable ,dimension(:) ,private :: speciesNames

        contains

        procedure ,public :: getSpeciesFromID
        procedure ,public :: getSpeciesFromName

        procedure ,public :: getSpeciesVarFromID
        procedure ,public :: getSpeciesVarFromName

        procedure ,public :: getSpeciesIDs 

        procedure ,public :: init => initSpeciesList

    end type SpeciesList
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initSpeciesList(this,jsonCont,mpiCont)
        !! SpeciesList initialization from config.json.

        class(SpeciesList)      ,intent(inout) :: this
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
        
    end subroutine initSpeciesList 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getSpeciesFromID (this,id) result(spec)
        !! Getter for species based on id

        class(SpeciesList)   ,intent(in) :: this
        integer(ik)          ,intent(in) :: id 
        type(Species) ,allocatable       :: spec

    end function getSpeciesFromID
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getSpeciesIDs (this) result(ids)
        !! Getter for speciesIDs

        class(SpeciesList)                      ,intent(in) :: this
        integer(ik)  ,allocatable, dimension(:)             :: ids

    end function getSpeciesIDs
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getSpeciesFromName (this,name) result(spec)
        !! Getter for species based on name

        class(SpeciesList)   ,intent(in) :: this
        character(*)         ,intent(in) :: name 
        type(Species) ,allocatable      :: spec

    end function getSpeciesFromName
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getSpeciesVarFromName (this,name,ind) result(var)
        !! Return associated var with index ind of species with given name

        class(SpeciesList)   ,intent(in) :: this
        character(*)         ,intent(in) :: name
        integer(ik)          ,intent(in) :: ind 
        character(:) ,allocatable        :: var

    end function getSpeciesVarFromName
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getSpeciesVarFromID (this,id,ind) result(var)
        !! Return associated var with index ind of species with given id

        class(SpeciesList)   ,intent(in) :: this
        integer(ik)          ,intent(in) :: id 
        integer(ik)          ,intent(in) :: ind 
        character(:) ,allocatable        :: var

    end function getSpeciesVarFromID
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module species_list_class
!-----------------------------------------------------------------------------------------------------------------------------------
 