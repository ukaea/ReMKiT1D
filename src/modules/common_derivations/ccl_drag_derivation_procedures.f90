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
submodule (ccl_drag_derivation_class) ccl_drag_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the CCL drag coefficient derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCCLDrag(this,vSpaceObj)
    !! Initialize Chang-Cooper-Langdon drag coefficient derivation

    class(CCLDragDerivation)     ,intent(inout) :: this
    type(VSpace)                 ,intent(in)    :: vSpaceObj

    if (assertions) call assert(vSpaceObj%isDefined(),"Undefined VSpace object passed to initCCLDrag")

    this%numV = vSpaceObj%getNumV()
    this%v2dv = vSpaceObj%getVGrid()**2 * vSpaceObj%getVCellWidths()
    call this%makeDefined()
    
end subroutine initCCLDrag  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateCCLDrag(this,inputArray,indices) result(output)

    class(CCLDragDerivation)           ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i  ,j ,inferredNumX ,lboundInput ,indOffset

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateCCLDrag called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateCCLDrag must be 1")
        call assertPure(all(indices>0),"indices passed to calculateCCLDrag out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateCCLDrag out of bounds - upper")

    end if

    inferredNumX = size(inputArray(indices(1))%entry)/this%numV

    lboundInput = lbound(inputArray(indices(1))%entry,1)

    allocate(output(inferredNumX*this%numV))
    output = 0
    do i = 1, inferredNumX
        indOffset = (i-1)*this%numV + lboundInput - 1

        do j = 1,this%numV-1
            output((i-1)*this%numV+j) = 4*pi*dot_product(inputArray(1)%entry(indOffset+1:indOffset+j),this%v2dv(1:j))
        end do
    end do

end function calculateCCLDrag
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule ccl_drag_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
