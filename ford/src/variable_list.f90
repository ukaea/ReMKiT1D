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
module variable_list_class
    !! Houses list of variables containing their names and whether they're a distribution function

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_types               ,only: StringArray

    implicit none
    private

    type ,public ,extends(Object) :: VariableList
        !! Contains a list of variable names and records whether they are distribution functions

        type(StringArray) ,allocatable ,dimension(:) ,private :: names !! Names of variables
        logical           ,allocatable ,dimension(:) ,private :: distf !! True for any variable which is a distribution function
        logical           ,allocatable ,dimension(:) ,private :: isSingleHarmonic !! True for any distribution variable which represent only a single harmonic
        logical           ,allocatable ,dimension(:) ,private :: isScalar !! True for any scalar variables
        logical           ,allocatable ,dimension(:) ,private :: isOnDualGrid !! True for any variable which lives on the dual grid
        logical           ,allocatable ,dimension(:) ,private :: isStationary !! True for any variable which has d/dt = 0
        integer(ik)       ,allocatable ,dimension(:) ,private :: priority !! Integer priority for uses such as variable derivation. Defaults to 0. 

        contains

        procedure ,public :: getVarNames
        procedure ,public :: getVarName
        procedure ,public :: getNumVars
        procedure ,public :: getVarPriority
        procedure ,public :: addVar

        procedure ,public :: isVarDist
        procedure ,public :: isVarSingleHarmonic
        procedure ,public :: isVarScalar
        procedure ,public :: isVarNameRegistered
        procedure ,public :: isVarOnDualGrid
        procedure ,public :: isVarStationary
        procedure ,public :: getVarIndex 

        procedure ,public :: combineWith

        procedure ,public :: init => initVarList

    end type VariableList
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initVarList(this) 
            !! Variable list initialization routine

            class(VariableList)           ,intent(inout)  :: this

        end subroutine initVarList
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumVars(this) result(numVars)
            !! Returns number of variables in list

            class(VariableList)  ,intent(in) :: this
            integer(ik)                      :: numVars
 
        end function getNumVars
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarNames(this) result(names)
            !! Getter of names

            class(VariableList)                          ,intent(in) :: this
            type(StringArray) ,allocatable ,dimension(:)             :: names
 
        end function getVarNames
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarName(this,ind) result(name)
            !! Return variable name at index ind

            class(VariableList)                  ,intent(in) :: this
            integer(ik)                          ,intent(in) :: ind
            character(:) ,allocatable                        :: name
 
        end function getVarName
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarPriority(this,ind) result(priority)
            !! Return priority of variable with given index 
        
            class(VariableList)  ,intent(in) :: this
            integer(ik)          ,intent(in) :: ind
            integer(ik)                      :: priority
        
        end function getVarPriority
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine addVar(this,name,isDist,isSingleHarmonic,isScalar,isOnDualGrid,isStationary,priority) 
            !! Add variable with given name to list; isDist determines if variable is a distribution function and is .false. by default;
            !! isSingleHarmonic sets whether a given distribution is only a single harmonic - i.e. a function of just x and v
            !! isScalar tags variable as a scalar (will be stored as a dimension 1 array)
            !! isOnDualGrid marks variable as living on the dual/staggered grid (or having staggered harmonics if it's a distribution)
            !! isStationary marks variable as having d/dt=0
            !! priority is an integer governing operations such as variable derivation

            class(VariableList) ,intent(inout)  :: this
            character(*)        ,intent(in)     :: name
            logical ,optional   ,intent(in)     :: isDist
            logical ,optional   ,intent(in)     :: isSingleHarmonic
            logical ,optional   ,intent(in)     :: isScalar
            logical ,optional   ,intent(in)     :: isOnDualGrid
            logical ,optional   ,intent(in)     :: isStationary
            integer(ik) ,optional ,intent(in)   :: priority

        end subroutine addVar
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarDist(this,ind) result(distf)
            !! Check whether variable with given index is a full distribution function

            class(VariableList) ,intent(in)  :: this
            integer(ik)         ,intent(in)  :: ind
            logical                          :: distf

        end function isVarDist
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarSingleHarmonic(this,ind) result(singleH)
            !! Check whether variable with given index is a single harmonic

            class(VariableList) ,intent(in)  :: this
            integer(ik)         ,intent(in)  :: ind
            logical                          :: singleH

        end function isVarSingleHarmonic
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarScalar(this,ind) result(scal)
            !! Check whether variable with given index is a scalar

            class(VariableList) ,intent(in)  :: this
            integer(ik)         ,intent(in)  :: ind
            logical                          :: scal

        end function isVarScalar
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarOnDualGrid(this,ind) result(dual)
            !! Check whether variable with given index is a on dual grid

            class(VariableList) ,intent(in)  :: this
            integer(ik)         ,intent(in)  :: ind
            logical                          :: dual

        end function isVarOnDualGrid
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarStationary(this,ind) result(stat)
            !! Check whether variable with given index is stationary

            class(VariableList) ,intent(in)  :: this
            integer(ik)         ,intent(in)  :: ind
            logical                          :: stat

        end function isVarStationary
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarNameRegistered(this,name) result(reg)
            !! Check whether variable with given name is registered

            class(VariableList) ,intent(in)  :: this
            character(*)        ,intent(in)  :: name
            logical                          :: reg

        end function isVarNameRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarIndex(this,name) result(ind)
            !! Get index of variable with given name

            class(VariableList)  ,intent(in) :: this
            character(*)         ,intent(in) :: name
            integer(ik)                      :: ind
 
        end function getVarIndex
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function combineWith(this,other) result(res)
            !! Combine two lists into one
        
            class(VariableList)  ,intent(in) :: this
            type(VariableList)   ,intent(in) :: other
            type(VariableList)               :: res
        
        end function combineWith
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module variable_list_class
!-----------------------------------------------------------------------------------------------------------------------------------
 