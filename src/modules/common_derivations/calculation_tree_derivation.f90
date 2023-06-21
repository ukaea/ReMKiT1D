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
module calculation_tree_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation class with an attached calculation tree object

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation
    use calculation_tree_class      ,only: CalculationTree, FlatTree

    implicit none
    private

    type ,public ,extends(Derivation) :: CalculationTreeDerivation
        !! Derivation object that has a calculation tree component, used to evaluate the derivation result independently of the passed 
        !! indices. The derivation is constructed using a FlatTree, which needs to be unpacked before the derivation can be used. 
        !! Note: This derivation should only be used if the passed object is the global VariableContainers variable array.

        type(CalculationTree) ,allocatable ,private :: tree 
        type(FlatTree)                     ,private :: flattenedTree 

        integer(ik) ,pointer ,private :: testPointer => null() !! Used to check whether the tree needs to be reconstructed. Should become unassociated if this derivation object goes through a copy/allocation which is pointer unsafe.

        contains

        procedure ,public :: init => initCalculationTreeDeriv

        procedure ,public :: calculate => calculateTree

    end type CalculationTreeDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCalculationTreeDeriv(this,flattenedTree)
        !! Initialize calculation tree derivation using flattened tree

        class(CalculationTreeDerivation)   ,intent(inout) :: this
        type(FlatTree)                     ,intent(in)    :: flattenedTree

    end subroutine initCalculationTreeDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateTree(this,inputArray,indices) result(output)

        class(CalculationTreeDerivation)            ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateTree
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module calculation_tree_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 