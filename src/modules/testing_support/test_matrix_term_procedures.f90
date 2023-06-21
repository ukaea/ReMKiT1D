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
submodule (test_matrix_term_class) test_matrix_term_procedures

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine testRow(this,varCont,rowVals,indexingData) 

    class(TestMatrixTerm)         ,intent(inout)   :: this
    type(VariableContainer)       ,intent(in)      :: varCont
    real(rk) ,dimension(:)        ,intent(inout)   :: rowVals
    type(MatrixTermIndexingData)  ,intent(in)      :: indexingData

    rowVals = varCont%variables(indexingData%rowReqVarIndices(1))%entry(indexingData%localRowIndices) *&
     varCont%variables(indexingData%rowReqVarIndices(2))%entry(indexingData%localRowIndices)

end subroutine testRow 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine testCol(this,varCont,colVals,indexingData) 

    class(TestMatrixTerm)         ,intent(inout) :: this
    type(VariableContainer)       ,intent(in)    :: varCont
    type(RealArray) ,dimension(:) ,intent(inout) :: colVals
    type(MatrixTermIndexingData)  ,intent(in)    :: indexingData

    integer(ik) :: i

    do i = 1,size(colVals)
        colVals(i)%entry = varCont%variables(indexingData%colReqVarIndices(1))%entry(indexingData%localColIndices(i)%entry)
    end do

end subroutine testCol 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule test_matrix_term_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
