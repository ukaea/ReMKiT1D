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
module group_evaluator_class
    !! author: Stefan Mijin
    !!
    !! Houses GroupEvaluator class, a Manipulator that evaluates specific model term groups

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use modeller_surrogate_class              ,only: ModellerSurrogate
    use variable_container_class              ,only: VariableContainer
    use manipulator_abstract_class            ,only: Manipulator
    use modeller_class                        ,only: Modeller

    implicit none
    private

    type ,public ,extends(Manipulator) :: GroupEvaluator
        !! Manipulator that evaluates a specific model term group and writes the result into variable with given index

        integer(ik) ,private :: resultVarIndex 
        integer(ik) ,private :: evaluatedModelIndex 
        integer(ik) ,private :: evaluatedTermGroup

        contains

        procedure ,public :: init => initEvaluator
        procedure ,public :: affect => evaluate

    end type GroupEvaluator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initEvaluator(this,resultVarIndex,modelIndex,termGroup) 
        !! GroupEvaluator initialization routine

        class(GroupEvaluator)    ,intent(inout)  :: this
        integer(ik)              ,intent(in)     :: resultVarIndex !! Index of variable to write the result in
        integer(ik)              ,intent(in)     :: modelIndex     !! Index of model whose term group should be evaluated
        integer(ik)              ,intent(in)     :: termGroup     !! Term group to evaluate in model

    end subroutine initEvaluator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine evaluate(this,manipulatedModeller,outputVars,inputVars) 
        !! Implementation of abstract manipulate routine for the evaluator

        class(GroupEvaluator)                 ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during manipulation
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the manipulation output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the manipulation routine

    end subroutine evaluate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module group_evaluator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 