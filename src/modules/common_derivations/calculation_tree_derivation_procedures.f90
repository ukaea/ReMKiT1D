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
submodule (calculation_tree_derivation_class) calculation_tree_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the calculation tree derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCalculationTreeDeriv(this,flattenedTree)
    !! Initialize calculation tree derivation using flattened tree

    class(CalculationTreeDerivation)   ,intent(inout) :: this
    type(FlatTree)                     ,intent(in)    :: flattenedTree

    this%flattenedTree = flattenedTree

    call this%makeDefined()

end subroutine initCalculationTreeDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateTree(this,inputArray,indices) result(output)

    class(CalculationTreeDerivation)    ,intent(inout) :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    if (assertions) call assert(this%isDefined(),"calculateTree called on undefined CalculationTreeDerivation")

    if (.not. associated(this%testPointer)) then

        if (allocated(this%tree)) deallocate(this%tree)

        allocate(this%tree)
        call this%tree%initFromFlatTree(this%flattenedTree)
        allocate(this%testPointer)
        this%testPointer = 1

    end if

    output = this%tree%evaluate(inputArray)

end function calculateTree
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule calculation_tree_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
