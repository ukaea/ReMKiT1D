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
submodule (group_evaluator_class) group_evaluator_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the GroupEvaluator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initEvaluator(this,resultVarIndex,modelIndex,termGroup) 
    !! GroupEvaluator initialization routine

    class(GroupEvaluator)    ,intent(inout)  :: this
    integer(ik)              ,intent(in)     :: resultVarIndex !! Index of variable to write the result in
    integer(ik)              ,intent(in)     :: modelIndex     !! Index of model whose term group should be evaluated
    integer(ik)              ,intent(in)     :: termGroup     !! Term group to evaluate in model

    this%resultVarIndex = resultVarIndex 
    this%evaluatedModelIndex = modelIndex 
    this%evaluatedTermGroup = termGroup

    call this%makeDefined()

end subroutine initEvaluator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine evaluate(this,manipulatedModeller,outputVars,inputVars) 
    !! Implementation of abstract manipulate routine for the evaluator

    class(GroupEvaluator)                      ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during manipulation
    class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the manipulation output 
    class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the manipulation routine

    if (assertions) then 
        call assert(this%isDefined(),"evaluate routine called by undefined GroupEvaluator")
        call assert(manipulatedModeller%isDefined(),"Undefined Modeller object passed to GroupEvaluator evaluate routine")
        call assert(outputVars%isDefined(),"outputVars passed to evaluate not defined")
        call assert(inputVars%isDefined(),"inputVars passed to evaluate not defined")
    end if

    select type (manipulatedModeller)
    type is (Modeller)
        outputVars%variables(this%resultVarIndex)%entry = &
        manipulatedModeller%evaluateModelTermGroup(this%evaluatedModelIndex,this%evaluatedTermGroup,inputVars)
    class default
        error stop "Unsupported surrogate passed to GroupEvaluator"
    end select

end subroutine evaluate
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule group_evaluator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
