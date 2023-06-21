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
    
    type ,public :: ExplicitTermIndexingData
        !! Indexing data used by the default abstract explicit term

        integer(ik)         ,allocatable ,dimension(:) :: reqVarIndices !! Indices of variables required by the multiplicative function
        character(:)        ,allocatable               :: evolvedVarName !! Name of evolved variable
        character(:)        ,allocatable               :: operatedVarName !! Name of variable to which the term's operator should be applied
        
    end type

    type ,public ,extends(term) ,abstract :: ExplicitTerm
        !! Abstract term providing basic interfaces for build terms used in explicit timestepping

        real(rk)            ,allocatable ,dimension(:) ,private :: multConst !! Multiplicative constant with size conforming to the variable being evolved
        real(rk)                                       ,private :: normalizationConst = real(1.0d0,kind=rk) !! Normalization constant of this term 

        type(ExplicitTermIndexingData)                 ,private :: indexingData !! Indexing data used by this term 
        
        class(Operator)     ,allocatable               ,private :: termOperator !! Optional operator used by this term

        contains

        procedure ,public  :: setOperator
        procedure ,public  :: setNormalizationConst
        procedure ,public  :: getNormalizationConst
        procedure ,public  :: setReqVars
        procedure ,public  :: setEvolvedAndOperatedVar
        procedure ,public  :: evaluate => evaluateExpTerm

        procedure ,public  :: getMultConst
        procedure ,public  :: setMultConst

        procedure ,private :: multFun => unityFun
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
    pure module function unityFun(this,varCont,indexingData) result(res)
        !! Default multiplicative function - returns a real vector of ones based on evolvedVarName. In general should use passed
        !! variable container and indexing data to construct a vector conforming to the evolved variable

        class(ExplicitTerm)             ,intent(in)   :: this
        type(VariableContainer)         ,intent(in)   :: varCont
        type(ExplicitTermIndexingData)  ,intent(in)   :: indexingData

        real(rk) ,allocatable           ,dimension(:) :: res  

    end function unityFun
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setOperator(this,op) 
        !! Setter for termOperator

        class(ExplicitTerm)       ,intent(inout)  :: this
        class(Operator)           ,intent(in)     :: op

    end subroutine setOperator
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNormalizationConst(this,norm) 
        !! Setter for normalizationConst

        class(ExplicitTerm)       ,intent(inout)  :: this
        real(rk)                  ,intent(in)     :: norm

    end subroutine setNormalizationConst
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNormalizationConst (this) result(norm)
        !! Getter for normalizationConst

        class(ExplicitTerm)  ,intent(in) :: this
        real(rk)                         :: norm

    end function getNormalizationConst
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setReqVars(this,reqVars,varCont) 
        !! Set variable names required by the multiplicative function and find their indices in variable container

        class(ExplicitTerm)               ,intent(inout)  :: this
        type(StringArray) ,dimension(:)   ,intent(in)     :: reqVars
        type(VariableContainer)           ,intent(in)     :: varCont

    end subroutine setReqVars
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setEvolvedAndOperatedVar(this,evolvedVarName,operatedVarName) 
        !! Set evolved and operated (if applicable) variable names

        class(ExplicitTerm)       ,intent(inout)  :: this
        character(*)              ,intent(in)     :: evolvedVarName
        character(*)   ,optional  ,intent(in)     :: operatedVarName

    end subroutine setEvolvedAndOperatedVar
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function evaluateExpTerm (this,varCont) result(res)
        !! Evaluates the term as multConst * multFun [* operatorTerm%actOn(<operatedVar>)] depending on whether the Operator is allocated

        class(ExplicitTerm)                  ,intent(in) :: this
        type(VariableContainer)              ,intent(in) :: varCont
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateExpTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setMultConst(this,multConst) 
        !! Setter for multConst

        class(ExplicitTerm)           ,intent(inout)  :: this
        real(rk)        ,dimension(:) ,intent(in)     :: multConst

    end subroutine setMultConst
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getMultConst (this) result(multConst)
        !! Getter for multConst

        class(ExplicitTerm)                          ,intent(in) :: this
        real(rk)          ,allocatable ,dimension(:)             :: multConst

    end function getMultConst
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
 