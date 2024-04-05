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
module derivation_explicit_term_class
    !! author: Stefan Mijin 
    !!
    !! Houses term calculated as modelbound data variable * derivation result

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions, assertionLvl
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray ,IntArray ,StringArray
    use variable_container_class               ,only: VariableContainer
    use partition_class                        ,only: Partition
    use grid_class                             ,only: Grid
    use indexing_class                         ,only: Indexing
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate
    use model_class                            ,only: Model
    use support_functions                      ,only: removeDupeInts ,allCombinations
    use explicit_term_abstract_class           ,only: ExplicitTerm
    use derivation_abstract_class              ,only: Derivation

    implicit none
    private

    type ,public ,extends(ExplicitTerm) :: DerivationTerm
        !! Term of the form modelbound variable * derivation result, where the modelbound variable is optional 

        real(rk) ,allocatable ,dimension(:) ,private :: resultBuffer !! Buffer holding the result

        logical                             ,private :: isActive 

        logical                             ,private :: kineticRow !! True if row variable is a distribution

        character(:) ,allocatable           ,private :: mbVarName !! Name of the mb variable to multiply de derivation result with

        class(Derivation), allocatable      ,private :: derivObj !! The main derivaiton object 

        integer(ik)       ,dimension(:) ,allocatable      ,private :: derivReqIndices !! Required indices for the derivation object

        integer(ik) ,private :: locNumX !! Local number of spatial cells - used when handling modelbound data 
        integer(ik) ,private :: numH !! Local number of harmonics - used when handling modelbound data 
        integer(ik) ,private :: numV !! Local number of velocity cells - used when handling modelbound data 

        contains

        procedure ,public :: update => updateDerivationTerm

        procedure ,public :: outerFun => derivationOuterFun
        procedure ,public :: init => initDerivationTerm

    end type DerivationTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDerivationTerm(this,gridObj,partitionObj,procRank,evolvedVar,varCont,derivObj,derivIndices,mbVarName)
        !! Explicit derivation term initialization routine 
    
        class(DerivationTerm)                       ,intent(inout)  :: this
        type(Grid)                                  ,intent(in)     :: gridObj
        type(Partition)                             ,intent(in)     :: partitionObj !! Parition object used to determine local number of DoF
        integer(ik)                                 ,intent(in)     :: procRank !! Current processor rank
        character(*)                                ,intent(in)     :: evolvedVar !! Name of evolved variable
        type(VariableContainer)                     ,intent(in)     :: varCont !! Reference variable container
        class(Derivation)                           ,intent(in)     :: derivObj !! Derivation object used by the term 
        integer(ik)              ,dimension(:)      ,intent(in)     :: derivIndices !! Required variable indices for the derivation object
        character(*)          ,optional             ,intent(in)     :: mbVarName !! Optional modelbound variable with which to multiply the derivation result 

    end subroutine initDerivationTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateDerivationTerm(this,varCont,modelData,hostModel) 
        !! Update function, does not call update on operator, assuming it is never set. Only updates the modelbound data buffer if
        !! required, and calculates the result buffer

        class(DerivationTerm)           ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update 
        class(ModelboundData) ,optional ,intent(in)     :: modelData !! Reference model data - unused
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused

    end subroutine updateDerivationTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    module function derivationOuterFun(this,varCont) result(res)
        !! Outer function simply returning the stored result buffer

        class(DerivationTerm)           ,intent(in)   :: this
        type(VariableContainer)         ,intent(in)   :: varCont

        real(rk) ,allocatable           ,dimension(:) :: res  

    end function derivationOuterFun
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module derivation_explicit_term_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
