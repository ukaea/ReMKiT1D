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
module test_explicit_term_class

    use data_kinds                     ,only: rk, ik
    use runtime_constants              ,only: debugging, assertions
    use god_objects                    ,only: Object
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use explicit_term_abstract_class   ,only: ExplicitTerm
    use variable_container_class       ,only: VariableContainer


    implicit none
    private

    type ,public ,extends(ExplicitTerm) :: TestExplicitTerm

        contains

        procedure ,public :: init => initTestExplicitTerm
        procedure ,public :: outerFun => testFun

    end type TestExplicitTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initTestExplicitTerm(this,evolvedVarName,varCont) 

        class(TestExplicitTerm)       ,intent(inout)  :: this
        character(*)                  ,intent(in)     :: evolvedVarName
        type(VariableContainer)       ,intent(in)     :: varCont

        end subroutine initTestExplicitTerm
!-----------------------------------------------------------------------------------------------------------------------------------
module function testFun(this,varCont) result(res)

    class(TestExplicitTerm)             ,intent(in)   :: this
    type(VariableContainer)         ,intent(in)   :: varCont

    real(rk) ,allocatable           ,dimension(:) :: res  

end function testFun
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module test_explicit_term_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
