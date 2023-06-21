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
submodule (term_evaluator_class) term_evaluator_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the TermEvaluator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initEvaluator(this,resultVarIndex,modelIndices,termNames) 
    !! TermEvaluator initialization routine

    class(TermEvaluator)            ,intent(inout)  :: this
    integer(ik)                     ,intent(in)     :: resultVarIndex !! Index of variable to write the result in
    integer(ik)       ,dimension(:) ,intent(in)     :: modelIndices   !! Indices of models whose named term should be evaluated
    type(StringArray) ,dimension(:) ,intent(in)     :: termNames     !! Name of evaluated term corresponding to each model 

    if (assertions) then 
        call assertPure(size(termNames)==size(modelIndices)&
        ,"termNames and modelIndices passed to TermEvaluator constructor must have the same size")
        call assertPure(size(modelIndices)>0,"modelIndices and termNames passed to TermEvaluator constructor must have size > 0")
    end if
    this%resultVarIndex = resultVarIndex 
    this%evaluatedModelIndex = modelIndices 
    this%evaluatedTermName = termNames

    call this%makeDefined()

end subroutine initEvaluator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine evaluate(this,manipulatedModeller,outputVars,inputVars) 
    !! Implementation of abstract manipulate routine for the evaluator

    class(TermEvaluator)                  ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during manipulation
    class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the manipulation output 
    class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the manipulation routine

    integer(ik) :: i

    if (assertions) then 
        call assert(this%isDefined(),"evaluate routine called by undefined TermEvaluator")
        call assert(manipulatedModeller%isDefined(),"Undefined Modeller object passed to TermEvaluator evaluate routine")
        call assert(outputVars%isDefined(),"outputVars passed to evaluate not defined")
        call assert(inputVars%isDefined(),"inputVars passed to evaluate not defined")
    end if

    select type (manipulatedModeller)
    type is (Modeller)
        outputVars%variables(this%resultVarIndex)%entry = &
        manipulatedModeller%evaluateModelTermByName(this%evaluatedModelIndex(1),this%evaluatedTermName(1)%string,inputVars)

        do i = 2,size(this%evaluatedModelIndex)
            outputVars%variables(this%resultVarIndex)%entry = outputVars%variables(this%resultVarIndex)%entry + &
            manipulatedModeller%evaluateModelTermByName(this%evaluatedModelIndex(i),this%evaluatedTermName(i)%string,inputVars)
        end do
    class default
        error stop "Unsupported surrogate passed to TermEvaluator"
    end select

end subroutine evaluate
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule term_evaluator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
