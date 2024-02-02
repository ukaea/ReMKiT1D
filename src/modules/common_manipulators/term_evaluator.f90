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
module term_evaluator_class
    !! author: Stefan Mijin
    !!
    !! Houses TermEvaluator class, a Manipulator that evaluates terms specified my modelIndex,termName tuples
    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions, assertionLvl
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use modeller_surrogate_class              ,only: ModellerSurrogate
    use variable_container_class              ,only: VariableContainer
    use manipulator_abstract_class            ,only: Manipulator
    use modeller_class                        ,only: Modeller
    use support_types

    implicit none
    private

    type ,public ,extends(Manipulator) :: TermEvaluator
        !! Manipulator that evaluates a specific model,term pairs and writes the total result into variable with given index

        integer(ik)                                  ,private :: resultVarIndex 
        integer(ik)       ,allocatable ,dimension(:) ,private :: evaluatedModelIndex 
        type(StringArray) ,allocatable ,dimension(:) ,private :: evaluatedTermName

        contains

        procedure ,public :: init => initEvaluator
        procedure ,public :: affect => evaluate

    end type TermEvaluator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initEvaluator(this,resultVarIndex,modelIndices,termNames) 
        !! TermEvaluator initialization routine

        class(TermEvaluator)            ,intent(inout)  :: this
        integer(ik)                     ,intent(in)     :: resultVarIndex !! Index of variable to write the result in
        integer(ik)       ,dimension(:) ,intent(in)     :: modelIndices   !! Indices of models whose named term should be evaluated
        type(StringArray) ,dimension(:) ,intent(in)     :: termNames     !! Name of evaluated term corresponding to each model 

    end subroutine initEvaluator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine evaluate(this,manipulatedModeller,outputVars,inputVars) 
        !! Implementation of abstract manipulate routine for the evaluator

        class(TermEvaluator)                  ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during manipulation
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the manipulation output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the manipulation routine

    end subroutine evaluate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module term_evaluator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 