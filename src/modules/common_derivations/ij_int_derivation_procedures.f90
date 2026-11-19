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
submodule (ij_int_derivation_class) ij_int_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the Shkarofsky I/J integral derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initIJInt(this,vSpaceObj,ind,isJInt)
    !! Initialize Shkarofsky I/J integral derivation

    class(IJIntDerivation)     ,intent(inout) :: this
    type(VSpace)               ,intent(in)    :: vSpaceObj
    integer(ik)                ,intent(in)    :: ind !! Index of integral
    logical ,optional          ,intent(in)    :: isJInt !! If true the upper-triangular J integral is calculated instead of the I integral. Defaults to false.

    real(rk) ,allocatable ,dimension(:) :: vGrid ,vWidths
    if (assertions) call assert(vSpaceObj%isDefined(),"Undefined VSpace object passed to initIJInt")

    this%numV = vSpaceObj%getNumV()
    
    this%isJInt = .false.
    if (present(isJInt)) this%isJInt = isJInt

    vGrid = vSpaceObj%getVGrid()
    vWidths = vSpaceObj%getVCellWidths()
    this%weights = 4.0d0*pi*(vGrid**(2+ind))*vWidths
    this%invVPow = 1.0d0/(vGrid**ind)

    call this%makeDefined()
    
end subroutine initIJInt  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateIJInt(this,inputArray,indices) result(output)
    !! Evaluate the Shkarofsky I (forward prefix sum) or J (backward suffix sum) integral over each velocity slice of inputArray(indices(1))

    class(IJIntDerivation)             ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,j ,inferredHarmonics ,indOffset ,outOffset ,lboundInput ,inIdx
    real(rk)    :: runningSum ,term

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateIJInt called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateIJInt must be 1")
        call assertPure(all(indices>0),"indices passed to calculateIJInt out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateIJInt out of bounds - upper")

    end if

    inIdx = indices(1)
    inferredHarmonics = size(inputArray(inIdx)%entry)/this%numV

    allocate(output(inferredHarmonics*this%numV))
    lBoundInput = lbound(inputArray(inIdx)%entry,1)

    if (this%isJInt) then
        do i = 1, inferredHarmonics
            outOffset = (i-1)*this%numV
            indOffset = outOffset + lboundInput - 1
            runningSum = 0.0d0
            do j = this%numV, 1, -1
                term = this%weights(j) * inputArray(inIdx)%entry(indOffset+j)
                output(outOffset+j) = (runningSum + 0.5d0*term) * this%invVPow(j)
                runningSum = runningSum + term
            end do
        end do
    else
        do i = 1, inferredHarmonics
            outOffset = (i-1)*this%numV
            indOffset = outOffset + lboundInput - 1
            runningSum = 0.0d0
            do j = 1, this%numV
                term = this%weights(j) * inputArray(inIdx)%entry(indOffset+j)
                output(outOffset+j) = (runningSum + 0.5d0*term) * this%invVPow(j)
                runningSum = runningSum + term
            end do
        end do
    end if

end function calculateIJInt
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule ij_int_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
