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
module variable_container_class
    !! author: Stefan Mijin 
    !!
    !! Houses all variables, both implicit and derived, as well as relevant support types 

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_types               ,only: StringArray ,RealArray ,IntArray
    use variable_list_class         ,only: VariableList
    use partition_class             ,only: Partition
    use indexing_class              ,only: Indexing 
    use derivation_abstract_class   ,only: Derivation
    use support_functions           ,only: findIndices

    implicit none
    private

    type ,public ,extends(Object) :: CalculationRule 
        !! Object used to calculate derived variables using a derivation object and set of required variable names

        class(Derivation) ,allocatable               ,public :: derivationMethod !! Derivation object used to calculate derived variables
        type(StringArray) ,allocatable ,dimension(:) ,public :: requiredVarNames !! Names of variables required for the derivation 

        contains 

        procedure ,public :: init => initCalculationRule

    end type CalculationRule

    type ,public ,extends(Object) :: VariableContainer
        !! Container for variable data, allows for handling conversion between local data storage and flattened vectors used for PETSc
        !! NOTE: Does not support single harmonic variables

        type(RealArray)       ,allocatable ,dimension(:) ,public  :: variables !! Value arrays for each stored variable in implicit,derived order
        type(VariableList)                               ,private :: implicitVarList !! List of implicit variables in container
        type(VariableList)                               ,private :: derivedVarList !! List of derived variables in container
        type(CalculationRule) ,allocatable ,dimension(:) ,private :: derivationRules !! List of derivation rules - conforms in size with derived list
        type(IntArray)        ,allocatable ,dimension(:) ,private :: requiredDerivationIndices !! Indices of variables required by each derivation rule
        type(IntArray)        ,allocatable ,dimension(:) ,private :: implicitVarIndices !! Indices in local implicit vector corresponding to each variabls implicitToLocIndex
        type(IntArray)        ,allocatable ,dimension(:) ,private :: implicitToLocIndex !! Indices in each variable corresponding to their indices in the implicit vector
        integer(ik)           ,allocatable ,dimension(:) ,private :: derivationDepth !! Numbers of layers of derived variables on which each derivation depends (e.g. 0 if all required variables are implicit, 1 if all required variables have derivation depth 0, etc.) - -1 for implicit variables

        integer(ik) ,private :: maxDerivPriority !! Highest priority value (lowest priority) among derived variables, used in calculating derived quantities
        integer(ik) ,allocatable ,dimension(:) ,private :: varLens !! Lengths of variables not including the halos 
        contains

        procedure ,public :: getVarIndex
        procedure ,public :: isVarNameRegistered
        procedure ,public :: calculateDerivedVars
        procedure ,public :: extractImplicitVars
        procedure ,public :: copyImplicitVarsToVec
        procedure ,public :: isVarDist
        procedure ,public :: isVarScalar
        procedure ,public :: getVarName
        procedure ,public :: getAllVarNames
        procedure ,public :: getImplicitVarNames
        procedure ,public :: isVarImplicit
        procedure ,public :: isVarOnDualGrid
        procedure ,public :: getVarDepth
        procedure ,public :: getMaxDepth
        procedure ,public :: copyNamedVarsToVec
        procedure ,public :: copyNamedVarsFromVec
        procedure ,public :: getVarLens 

        procedure ,public :: isStationary

        procedure ,public :: init => initVarContainer

    end type VariableContainer
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initCalculationRule(this,deriv,names) 
            !! Calculation rule object initialization routine

            class(CalculationRule)                    ,intent(inout)  :: this
            class(Derivation) ,optional               ,intent(in)     :: deriv !! Derivation component
            type(StringArray) ,optional ,dimension(:) ,intent(in)     :: names !! Required variable names

        end subroutine initCalculationRule
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initVarContainer(this,&
                                                implicitVars,&
                                                derivedVars,&
                                                derivationRules,&
                                                indexingObj,&
                                                partitionObj,&
                                                xHaloWidth,&
                                                procRank) 
        !! Variable container initialization routine
        
            class(VariableContainer)            ,intent(inout)  :: this
            type(VariableList)                  ,intent(in)     :: implicitVars !! List of implicit variables
            type(VariableList)                  ,intent(in)     :: derivedVars !! List of derived variables
            type(CalculationRule) ,dimension(:) ,intent(in)     :: derivationRules  !! Derivation rules corresponding to derived variables
            type(Indexing)                      ,intent(in)     :: indexingObj !! Reference indexing object used to map between variables and implicit vectors
            type(Partition)                     ,intent(in)     :: partitionObj !! Reference partition object 
            integer(ik)                         ,intent(in)     :: xHaloWidth !! Width of halo in x-direction 
            integer(ik)                         ,intent(in)     :: procRank !! Current processor rank

        end subroutine initVarContainer
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarIndex(this,name) result(ind)
            !! Get index of variable with given name

            class(VariableContainer)  ,intent(in) :: this
            character(*)              ,intent(in) :: name
            integer(ik)                           :: ind
 
        end function getVarIndex
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine calculateDerivedVars(this,derivPriority,derivDepth)
        !! Calculate derived variables from implicit variables using derivation rules. If derivPriority is supplied only those variables with
        !! priority <= derivPriority will be derived. If derivDepth is present only variables at that derivation depth are calculated,
        !! otherwise all variables are calculated

            class(VariableContainer) ,intent(inout)  :: this
            integer(ik) ,optional    ,intent(in)     :: derivPriority
            integer(ik) ,optional    ,intent(in)     :: derivDepth

        end subroutine calculateDerivedVars
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine extractImplicitVars(this,vec)
            !! Extract implicit variables from a locally linearly indexed real vector

            class(VariableContainer)    ,intent(inout)  :: this
            real(rk) ,dimension(:)      ,intent(in)     :: vec

        end subroutine extractImplicitVars
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine copyImplicitVarsToVec(this,vec,ignoreStationary)
            !! Copy variables into a locally linearly indexed real vector. If ignoreStationary is true copies any stationary variables as 0

            class(VariableContainer)    ,intent(inout)  :: this
            real(rk) ,dimension(:)      ,intent(inout)  :: vec
            logical ,optional           ,intent(in)     :: ignoreStationary

        end subroutine copyImplicitVarsToVec
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarDist(this,ind) result(distf)
            !! Check whether variable with given index is a distribution function

            class(VariableContainer) ,intent(in)  :: this
            integer(ik)              ,intent(in)  :: ind
            logical                               :: distf

        end function isVarDist
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarScalar(this,ind) result(scal)
            !! Check whether variable with given index is a scalar 

            class(VariableContainer) ,intent(in)  :: this
            integer(ik)              ,intent(in)  :: ind
            logical                               :: scal

        end function isVarScalar
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarNameRegistered(this,name) result(reg)
            !! Check whether variable with given name is registered in either the implicit or derived list

            class(VariableContainer) ,intent(in)  :: this
            character(*)             ,intent(in)  :: name
            logical                               :: reg

        end function isVarNameRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarImplicit(this,name) result(imp)
            !! Check whether variable with given name is registered in either the implicit list

            class(VariableContainer) ,intent(in)  :: this
            character(*)             ,intent(in)  :: name
            logical                               :: imp

        end function isVarImplicit
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isVarOnDualGrid(this,name) result(dual)
            !! Check whether variable with given name is on dual grid

            class(VariableContainer) ,intent(in)  :: this
            character(*)             ,intent(in)  :: name
            logical                               :: dual

        end function isVarOnDualGrid
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarName(this,ind) result(name)
            !! Return variable name at index ind

            class(VariableContainer)             ,intent(in) :: this
            integer(ik)                          ,intent(in) :: ind
            character(:) ,allocatable                        :: name
 
        end function getVarName
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getAllVarNames(this) result(names)
            !! Return all variable names in this container

            class(VariableContainer)             ,intent(in) :: this
            type(StringArray) ,allocatable ,dimension(:)     :: names
 
        end function getAllVarNames
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getImplicitVarNames(this) result(names)
            !! Return all implicit variable names in this container

            class(VariableContainer)             ,intent(in) :: this
            type(StringArray) ,allocatable ,dimension(:)     :: names
 
        end function getImplicitVarNames
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarDepth(this,name) result(depth)
            !! Get depth of variable with given name

            class(VariableContainer)  ,intent(in) :: this
            character(*)              ,intent(in) :: name
            integer(ik)                           :: depth
 
        end function getVarDepth
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMaxDepth(this) result(depth)
            !! Get max of derivationDepth

            class(VariableContainer)  ,intent(in) :: this
            integer(ik)                           :: depth
 
        end function getMaxDepth
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isStationary(this,name) result(stationary)
            !! Check whether variable with given name is stationary

            class(VariableContainer) ,intent(in)  :: this
            character(*)             ,intent(in)  :: name
            logical                               :: stationary

        end function isStationary
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine copyNamedVarsToVec(this,vec,names)
            !! Copy variables into locally indexed vector by name 
            
            class(VariableContainer)                 ,intent(inout)  :: this
            real(rk) ,allocatable ,dimension(:)      ,intent(inout)  :: vec
            type(StringArray) ,dimension(:)          ,intent(in) :: names 

        end subroutine copyNamedVarsToVec
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine copyNamedVarsFromVec(this,vec,names)
            !! Copy variables from locally indexed vector by name 
            
            class(VariableContainer)                 ,intent(inout)  :: this
            real(rk)  ,dimension(:)                  ,intent(in)  :: vec
            type(StringArray) ,dimension(:)          ,intent(in) :: names 

        end subroutine copyNamedVarsFromVec
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVarLens(this,names) result(lens)
            !! Get lengths of variable data vectors (not including halos) based on a list of names
            
            class(VariableContainer)                 ,intent(in)  :: this
            type(StringArray) ,dimension(:)          ,intent(in)  :: names 
            integer(ik) ,allocatable ,dimension(:)                :: lens 

        end function getVarLens
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module variable_container_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
