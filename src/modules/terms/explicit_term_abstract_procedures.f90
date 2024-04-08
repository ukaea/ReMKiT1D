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
submodule (explicit_term_abstract_class) explicit_term_abstract_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the abstract explicit term class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine simpleUpdate(this,varCont,modelData,hostModel) 
    !! Default update function - just calls update on Operator if allocated

    class(ExplicitTerm)             ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update 
    class(ModelboundData) ,optional ,intent(in)     :: modelData !! Reference model data - unused
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused

    if (assertions) then
        call assert(this%isDefined(),"Attempted to update undefined explicit term object")
        call assert(varCont%isDefined(),"Attempted to update explicit term object by passing undefined variable container")
    end if

    if (allocated(this%termOperator)) call this%termOperator%update(varCont)

end subroutine simpleUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
module function unityFun(this,varCont) result(res)
    !! Default multiplicative function - returns a real vector of ones based on evolvedVarName. In general should use passed
    !! variable container to construct a vector conforming to the evolved variable

    class(ExplicitTerm)             ,intent(in)   :: this
    type(VariableContainer)         ,intent(in)   :: varCont

    real(rk) ,allocatable           ,dimension(:) :: res  

    if (assertions) then
        call assertPure(this%isDefined(),"Called unityFun with undefined explicit term object")
        call assertPure(varCont%isDefined(),"Called unityFun by passing undefined variable container")
    end if

    allocate(res,mold=varCont%variables(varCont%getVarIndex(this%evolvedVarName))%entry)
    res = real(1.0d00,kind=rk)

end function unityFun
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setOperator(this,op) 
    !! Setter for termOperator

    class(ExplicitTerm)       ,intent(inout)  :: this
    class(Operator)           ,intent(in)     :: op

    if (allocated(this%termOperator)) deallocate(this%termOperator)
    allocate(this%termOperator,source=op)

end subroutine setOperator
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setEvolvedVar(this,evolvedVarName)
    !! Set evolved variable names

    class(ExplicitTerm)       ,intent(inout)  :: this
    character(*)              ,intent(in)     :: evolvedVarName

    this%evolvedVarName = evolvedVarName

end subroutine setEvolvedVar
!-----------------------------------------------------------------------------------------------------------------------------------
module function evaluateExpTerm (this,varCont) result(res)
    !! Evaluates the term as multConst * multFun [* operatorTerm%actOn(<operatedVar>)] depending on whether the Operator is allocated

    class(ExplicitTerm)                  ,intent(in) :: this
    type(VariableContainer)              ,intent(in) :: varCont
    real(rk) ,allocatable ,dimension(:)              :: res

    if (assertions) then
        call assertPure(this%isDefined(),"Attempted to evaluate undefined explicit term")
        call assertPure(varCont%isDefined(),&
        "Attempted to calculate explicit term by passing undefined variable container")
    end if
        res = this%outerFun(varCont)
    if (allocated(this%termOperator)) & 
        res = res * this%termOperator%actOn(this%innerFun(varCont))

end function evaluateExpTerm
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getEvolvedVarName(this) result(name)
    !! Get name of the evolved variable of this term

    class(ExplicitTerm)                  ,intent(in) :: this
    character(:) ,allocatable                        :: name

    if (assertions) call assertPure(this%isDefined(),"Attempted to get evolved variable name from undefined explicit term")

    name = this%evolvedVarName

end function getEvolvedVarName
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule explicit_term_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
