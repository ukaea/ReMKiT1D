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
module explicit_term_abstract_class
    !! author: Stefan Mijin 
    !! 
    !! Houses abstract explicit term interfaces

    use data_kinds                     ,only: rk, ik
    use runtime_constants              ,only: debugging, assertions
    use god_objects                    ,only: Object
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use support_types                  ,only: RealArray ,IntArray ,StringArray
    use variable_list_class            ,only: VariableList
    use variable_container_class       ,only: VariableContainer
    use term_abstract_class            ,only: Term
    use operator_abstract_class        ,only: Operator
    use modelbound_data_abstract_class ,only: ModelboundData
    use model_surrogate_class          ,only: ModelSurrogate

    implicit none
    private

    type ,public ,extends(term) ,abstract :: ExplicitTerm
        !! Abstract term providing basic interfaces for build terms used in explicit timestepping

        character(:)        ,allocatable               ,private :: evolvedVarName !! Name of evolved variable 
        class(Operator)     ,allocatable               ,private :: termOperator !! Optional operator used by this term

        contains

        procedure ,public  :: setOperator
        procedure ,public  :: setEvolvedVar
        procedure ,public  :: evaluate => evaluateExpTerm

        procedure ,private :: outerFun => unityFun
        procedure ,private :: innerFun  => unityFun

        procedure ,public  :: update => simpleUpdate

        procedure ,public  :: getVarName => getEvolvedVarName

    end type ExplicitTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine simpleUpdate(this,varCont,modelData,hostModel) 
        !! Default update function - just calls update on Operator if allocated

        class(ExplicitTerm)             ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update 
        class(ModelboundData) ,optional ,intent(in)     :: modelData !! Reference model data - unused
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused

    end subroutine simpleUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    module function unityFun(this,varCont) result(res)
    !! Default multiplicative function - returns a real vector of ones based on evolvedVarName. In general should use passed
    !! variable container to construct a vector conforming to the evolved variable

        class(ExplicitTerm)             ,intent(in)   :: this
        type(VariableContainer)         ,intent(in)   :: varCont

        real(rk) ,allocatable           ,dimension(:) :: res  

    end function unityFun
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setOperator(this,op) 
        !! Setter for termOperator

        class(ExplicitTerm)       ,intent(inout)  :: this
        class(Operator)           ,intent(in)     :: op

    end subroutine setOperator
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setEvolvedVar(this,evolvedVarName) 
        !! Set evolved variable name

        class(ExplicitTerm)       ,intent(inout)  :: this
        character(*)              ,intent(in)     :: evolvedVarName

    end subroutine setEvolvedVar
!-----------------------------------------------------------------------------------------------------------------------------------
    module function evaluateExpTerm (this,varCont) result(res)
        !! Evaluates the term as outerFun * [* operatorTerm%actOn(innerFun)] depending on whether the Operator is allocated

        class(ExplicitTerm)                  ,intent(in) :: this
        type(VariableContainer)              ,intent(in) :: varCont
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateExpTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getEvolvedVarName(this) result(name)
        !! Get name of the evolved variable of this term

        class(ExplicitTerm)                  ,intent(in) :: this
        character(:) ,allocatable                        :: name

    end function getEvolvedVarName
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module explicit_term_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 