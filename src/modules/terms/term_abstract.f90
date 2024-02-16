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
module term_abstract_class
    !! author: Stefan Mijin 
    !!
    !! Houses abstract term class determining basic term interfaces

    use data_kinds                     ,only: rk
    use god_objects                    ,only: object
    use variable_container_class       ,only: VariableContainer
    use runtime_constants              ,only: debugging, assertions
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use modelbound_data_abstract_class ,only: ModelboundData
    use model_surrogate_class          ,only: ModelSurrogate

    implicit none
    private

    type ,public ,extends(Object) ,abstract :: Term
        !! General abstract term object

        contains

        procedure                 ,public   :: update => noUpdate
        procedure(evaluateTerm)   ,deferred :: evaluate
        procedure(getEvolvedVar)  ,deferred :: getVarName

    end type Term
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    function evaluateTerm(this,varCont) result(res)
        !! Get explicit value for the term 

        import :: Term, VariableContainer ,rk

        class(Term)                          ,intent(in) :: this
        type(VariableContainer)              ,intent(in) :: varCont
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    pure function getEvolvedVar(this) result(name)
        !! Get name of the evolved variable of this term

        import :: Term

        class(Term)                          ,intent(in) :: this
        character(:) ,allocatable                        :: name

    end function getEvolvedVar
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine noUpdate(this,varCont,modelData,hostModel) 
            !! Default Term update function - does nothing

            class(Term)                     ,intent(inout)  :: this
            type(VariableContainer)         ,intent(in)     :: varCont
            class(ModelboundData) ,optional ,intent(in)     :: modelData
            class(ModelSurrogate) ,optional ,intent(in)     :: hostModel

        end subroutine noUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module term_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 