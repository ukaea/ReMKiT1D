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
submodule (test_explicit_term_class) test_explicit_term_procedures

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initTestExplicitTerm(this,evolvedVarName,varCont) 

    class(TestExplicitTerm)       ,intent(inout)  :: this
    character(*)                  ,intent(in)     :: evolvedVarName
    type(VariableContainer)       ,intent(in)     :: varCont

    real(rk) ,allocatable ,dimension(:) :: multConst

    call this%makeDefined()

    allocate(multConst,mold=varCont%variables(varCont%getVarIndex(evolvedVarName))%entry)

    multConst = 1.00d00

    call this%setMultConst(multConst)

end subroutine initTestExplicitTerm
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule test_explicit_term_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
