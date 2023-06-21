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
submodule (ccl_diff_derivation_class) ccl_diff_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the CCL diff coefficient derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCCLDiff(this,vSpaceObj)
    !! Initialize Chang-Cooper-Langdon diffusion coefficient derivation

    class(CCLDiffDerivation)     ,intent(inout) :: this
    type(VSpace)                 ,intent(in)    :: vSpaceObj

    real(rk) ,allocatable ,dimension(:) :: vGrid, dv

    if (assertions) call assert(vSpaceObj%isDefined(),"Undefined VSpace object passed to initCCLDiff")

    this%numV = vSpaceObj%getNumV()

    vGrid = vSpaceObj%getVGrid()
    dv = vSpaceObj%getVCellWidths()
    this%v2dv = vGrid**2 * dv 
    allocate(this%vdvPlus(this%numV))
    allocate(this%vStar(this%numV))
    this%vdvPlus = 0
    this%vStar = 0

    this%vStar(1:this%numV-1) = (vGrid(2:this%numV) + vGrid(1:this%numV-1))/2

    this%vdvPlus(1:this%numV-1) = (vGrid(2:this%numV) - vGrid(1:this%numV-1))*(vGrid(2:this%numV) + vGrid(1:this%numV-1))/2
    
    call this%makeDefined()
    
end subroutine initCCLDiff  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateCCLDiff(this,inputArray,indices) result(output)

    class(CCLDiffDerivation)           ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i  ,j ,k,inferredNumX ,lboundF ,indOffsetF, lboundDelta ,indOffsetDelta

    real(rk) ,allocatable ,dimension(:) :: fInterp ,dProduct

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateCCLDiff called on undefined derivation object")
        call assertPure(size(indices) == 2,"Number of indices passed to calculateCCLDiff must be 2")
        call assertPure(all(indices>0),"indices passed to calculateCCLDiff out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateCCLDiff out of bounds - upper")

    end if

    inferredNumX = size(inputArray(indices(1))%entry)/this%numV

    lboundF = lbound(inputArray(indices(1))%entry,1)
    lboundDelta = lbound(inputArray(indices(2))%entry,1)

    allocate(output(inferredNumX*this%numV))
    output = 0
    allocate(fInterp(this%numV))
    allocate(dProduct(this%numV-1))
    
    do i = 1, inferredNumX
        indOffsetF = (i-1)*this%numV + lboundF - 1
        indOffsetDelta = (i-1)*this%numV + lboundDelta - 1
        fInterp = 0
        fInterp(1:this%numV-1) = (real(1,kind=rk)-inputArray(indices(2))%entry(indOffsetDelta+1:indOffsetDelta+this%numV-1))&
                                 *inputArray(indices(1))%entry(indOffsetF+1:indOffsetF+this%numV-1)&
                                 + inputArray(indices(2))%entry(indOffsetDelta+1:indOffsetDelta+this%numV-1)&
                                 *inputArray(indices(1))%entry(indOffsetF+2:indOffsetF+this%numV)
        do j = 1,this%numV-1
            dProduct(j) = dot_product(fInterp(j:this%numV-1),this%vdvPlus(j:this%numV-1))*this%v2dv(j)
        end do

        do j = 1,this%numV-1

            output((i-1)*this%numV+j) = sum(dProduct(1:j))
            output((i-1)*this%numV+j) = output((i-1)*this%numV+j) *4*pi/this%vStar(j)
        end do
    end do

end function calculateCCLDiff
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule ccl_diff_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
