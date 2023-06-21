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
module test_matrix_term_class

    use data_kinds                     ,only: rk, ik
    use runtime_constants              ,only: debugging, assertions
    use god_objects                    ,only: Object
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use variable_container_class       ,only: VariableContainer
    use basic_interfaces               ,only: coordMapping
    use indexing_class                 ,only: Indexing
    use matrix_term_abstract_class     ,only: MatrixTerm ,MatrixTermIndexingData
    use support_types

    implicit none
    private

    type ,public ,extends(MatrixTerm) :: TestMatrixTerm

        contains

        procedure ,private :: calculateRows => testRow
        procedure ,private :: calculateCols => testCol

    end type TestMatrixTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine testRow(this,varCont,rowVals,indexingData) 

        class(TestMatrixTerm)         ,intent(inout)   :: this
        type(VariableContainer)       ,intent(in)      :: varCont
        real(rk) ,dimension(:)        ,intent(inout)   :: rowVals
        type(MatrixTermIndexingData)  ,intent(in)      :: indexingData

    end subroutine testRow 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine testCol(this,varCont,colVals,indexingData) 
        !! Default column function. In general should provide values for each column 

        class(TestMatrixTerm)         ,intent(inout) :: this
        type(VariableContainer)       ,intent(in)    :: varCont
        type(RealArray) ,dimension(:) ,intent(inout) :: colVals
        type(MatrixTermIndexingData)  ,intent(in)    :: indexingData
        
    end subroutine testCol  
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module test_matrix_term_class
!-----------------------------------------------------------------------------------------------------------------------------------
 