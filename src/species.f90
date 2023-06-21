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
module species_class
    !! author: Stefan Mijin
    !! 
    !! Houses base species class containing particle properties and associated variable names

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use mpi_controller_class                  ,only: MPIController
    use json_controller_class                 ,only: JSONController
    use key_names
    use support_types                         
    use physical_constants
    

    implicit none
    private

    type ,public ,extends(Object) :: Species
        !! Base species class containing particle properties and associated variable names

        integer(ik)                                  ,private :: id 
        character(:)      ,allocatable               ,private :: name 
        real(rk)                                     ,private :: charge !! in elementary charge units
        real(rk)                                     ,private :: mass !! in kg
        type(StringArray) ,allocatable ,dimension(:) ,private :: associatedVars

        contains

        procedure ,public :: getID
        procedure ,public :: getName
        procedure ,public :: getCharge 
        procedure ,public :: getMass 
        procedure ,public :: getAssociatedVars

        procedure ,public :: initFromJSON
        procedure ,public :: init => initSpecies

    end type Species
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initFromJSON(this,name,jsonCont,mpiCont)
        !! Initialize species from JSON file. Here for extensibility 

        class(Species)          ,intent(inout) :: this
        character(*)            ,intent(in)    :: name 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
        
    end subroutine initFromJSON  
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initSpecies(this,id,name,Z,A,vars)
        !! Species initialization. Z is the charge and A the atomic mass (in amus). A/Z are ignored if id = 0, which is reserved for electrons.  
        !! vars is a StringArray containing names of associated variables for easy access 

        class(Species)                  ,intent(inout) :: this
        integer(ik)                     ,intent(in)    :: id 
        character(*)                    ,intent(in)    :: name 
        real(rk)                        ,intent(in)    :: Z 
        real(rk)                        ,intent(in)    :: A 
        type(StringArray) ,dimension(:) ,intent(in)    :: vars 
        
    end subroutine initSpecies  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getID (this) result(id)
        !! Getter for id

        class(Species)       ,intent(in) :: this
        integer(ik)                      :: id

    end function getID
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getMass (this) result(mass)
        !! Getter for mass

        class(Species)    ,intent(in) :: this
        real(rk)                      :: mass

    end function getMass
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCharge (this) result(charge)
        !! Getter for charge

        class(Species)    ,intent(in) :: this
        real(rk)                      :: charge

    end function getCharge
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getName (this) result(name)
        !! Getter for name

        class(Species)    ,intent(in) :: this
        character(:) ,allocatable     :: name

    end function getName
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getAssociatedVars (this) result(associatedVars)
        !! Getter for associatedVars

        class(Species)                      ,intent(in) :: this
        type(StringArray) ,allocatable ,dimension(:)    :: associatedVars

    end function getAssociatedVars
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module species_class
!-----------------------------------------------------------------------------------------------------------------------------------
 