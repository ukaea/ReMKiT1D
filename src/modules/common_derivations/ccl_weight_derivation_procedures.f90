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
submodule (ccl_weight_derivation_class) ccl_weight_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the CCL interpolation weights derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCCLWeights(this,vSpaceObj)
    !! Initialize Chang-Cooper-Langdon interpolation weight derivation

    class(CCLWeightDerivation)     ,intent(inout) :: this
    type(VSpace)                   ,intent(in)    :: vSpaceObj

    real(rk) ,allocatable ,dimension(:) :: vGrid

    if (assertions) call assert(vSpaceObj%isDefined(),"Undefined VSpace object passed to initCCLWeights")

    this%numV = vSpaceObj%getNumV()
    allocate(this%dvPlus(this%numV))
    this%dvPlus = 0
    vGrid = vSpaceObj%getVGrid()
    this%dvPlus(1:this%numV-1) = vGrid(2:this%numV) -  vGrid(1:this%numV-1)
    call this%makeDefined()
    
end subroutine initCCLWeights  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateCCLWeights(this,inputArray,indices) result(output)

    class(CCLWeightDerivation)         ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,inferredNumX ,indOffsetC ,lboundC ,indOffsetD ,lboundD

    real(rk) ,allocatable ,dimension(:) :: wBuffer

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateCCLWeights called on undefined derivation object")
        call assertPure(size(indices) == 2,"Number of indices passed to calculateCCLWeights must be 2")
        call assertPure(all(indices>0),"indices passed to calculateCCLWeights out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateCCLWeights out of bounds - upper")

    end if

    inferredNumX = size(inputArray(indices(1))%entry)/this%numV

    allocate(wBuffer(inferredNumX*this%numV))

    lBoundC = lbound(inputArray(indices(1))%entry,1)
    lBoundD = lbound(inputArray(indices(2))%entry,1)

    do i = 1, inferredNumX
        indOffsetC = (i-1)*this%numV + lboundC - 1
        indOffsetD = (i-1)*this%numV + lboundD - 1

        wBuffer((i-1)*this%numV +1 : i*this%numV) = this%dvPlus * &
                                                    inputArray(indices(1))%entry(indOffsetC+1:indOffsetC+this%numV)&
                                                    /inputArray(indices(2))%entry(indOffsetD+1:indOffsetD+this%numV)

        wBuffer(i*this%numV) = ieee_value(wBuffer(i*this%numV),ieee_positive_inf)
    end do

    output =real(1,kind=rk) -( real(1,kind=rk)/wBuffer - real(1,kind=rk)/(exp(wBuffer)-real(1,kind=rk)))

end function calculateCCLWeights
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule ccl_weight_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
