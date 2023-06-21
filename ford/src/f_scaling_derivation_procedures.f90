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
submodule (f_scaling_derivation_class) f_scaling_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the distribution scaling extrapolation derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initFScaling(this,partitionObj,procRank,numV,leftBoundary,staggeredVars,extrapolateToBoundary)
    !! Initialize distribution scaling extrapolation derivation object

    class(FScalingDerivation)     ,intent(inout)  :: this
    type(Partition)               ,intent(in)     :: partitionObj !! Partition object used to determine local processor grid chunk
    integer(ik)                   ,intent(in)     :: procRank !! Current processor rank
    integer(ik)                   ,intent(in)     :: numV !! Number of expected velocity grid points
    logical  ,optional            ,intent(in)     :: leftBoundary !! True if extrapolating to left boundary. Defaults to false.
    logical  ,optional            ,intent(in)     :: staggeredVars !! True if the distribution has staggered harmonics. Defaults to false.
    logical  ,optional            ,intent(in)     :: extrapolateToBoundary !! True if the extrapolation should be performed to the cell boundary and not the last cell centre before the boundary.

    integer(ik) :: minX ,maxX, inferredGridSize

    if (assertions) call assert(partitionObj%isDefined(),"Partition passed to initFScaling not defined")

    this%leftBoundary = .false. 
    if (present(leftBoundary)) this%leftBoundary = leftBoundary

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    inferredGridSize = maxval(partitionObj%getMaxX())

    this%hasBoundary = (minX == 1 .and. this%leftBoundary) .or. (maxX == inferredGridSize .and. .not. this%leftBoundary)

    this%numV = numV 
    this%numH = maxval(partitionObj%getMaxH())
    if (this%hasBoundary) then 

        this%staggeredVars = .false. 
        if (present(staggeredVars)) this%staggeredVars = staggeredVars

        this%extrapolateToBoundary = .false. 
        if (present(extrapolateToBoundary)) this%extrapolateToBoundary = extrapolateToBoundary

        this%exterpCoords = [maxX-minX+1,maxX-minX]

        if (this%leftBoundary) this%exterpCoords = [1,1]

    end if

    call this%makeDefined()

end subroutine initFScaling  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateFScaling(this,inputArray,indices) result(output)

    class(FScalingDerivation)           ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:)  ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:,:)               :: output

    integer(ik) :: expectedNumIndices ,i ,offset

    real(rk) ,dimension(2) :: scalingFactor 

    if (assertions) call assertPure(this%isDefined(),"calculateFScaling called on undefined derivation object")

    expectedNumIndices = 2 
    if (this%staggeredVars) expectedNumIndices = expectedNumIndices + 1
    if (this%extrapolateToBoundary) expectedNumIndices = expectedNumIndices + 1

    if (assertions) then 

        call assertPure(size(indices) == expectedNumIndices,"Unexpected number of indices passed to calculateFScaling")
        call assertPure(all(indices>0),"indices passed to calculateFScaling out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateFScaling out of bounds - upper")

    end if

    allocate(output(this%numH,this%numV))

    output = 0 

    if (this%hasBoundary) then 
        scalingFactor = this%getScalingFactors(inputArray,indices)
        
        if (this%staggeredVars) then 

            do i = 1,this%numH
                if (mod(i-1,2)==1) then
                    offset = (this%exterpCoords(2)-1)*this%numH*this%numV + (i-1)*this%numV
                    output(i,:) = scalingFactor(2)*inputArray(indices(1))%entry(offset+1:offset+this%numV)
                else
                    offset = (this%exterpCoords(1)-1)*this%numH*this%numV + (i-1)*this%numV
                    output(i,:) = scalingFactor(1)*inputArray(indices(1))%entry(offset+1:offset+this%numV)
                end if

            end do
        else 

            do i = 1,this%numH
                offset = (this%exterpCoords(1)-1)*this%numH*this%numV + (i-1)*this%numV
                output(i,:) = scalingFactor(1)*inputArray(indices(1))%entry(offset+1:offset+this%numV)
            end do
        end if
    end if

end function calculateFScaling
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getScalingFactors(this,inputArray,indices) result(output)

    class(FScalingDerivation)            ,intent(in)    :: this 
    type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:)  ,intent(in)    :: indices           
    real(rk)              ,dimension(2)                 :: output

    integer(ik) :: expectedNumIndices 

    if (assertions) call assertPure(this%isDefined(),"getScalingFactors called on undefined derivation object")

    expectedNumIndices = 2 
    if (this%staggeredVars) expectedNumIndices = expectedNumIndices + 1
    if (this%extrapolateToBoundary) expectedNumIndices = expectedNumIndices + 1

    if (assertions) then 

        call assertPure(size(indices) == expectedNumIndices,"Unexpected number of indices passed to getScalingFactors")
        call assertPure(all(indices>0),"indices passed to getScalingFactors out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to getScalingFactors out of bounds - upper")

    end if
    output = real(1,kind=rk)
    if(this%extrapolateToBoundary) &
        output = inputArray(indices(expectedNumIndices))%entry(1)/inputArray(indices(2))%entry(this%exterpCoords(1))
        
     if (this%staggeredVars) then 
            output(2) =  inputArray(indices(2))%entry(this%exterpCoords(1))&
                                /inputArray(indices(3))%entry(this%exterpCoords(2))

            if(this%extrapolateToBoundary) &
            output(2) = inputArray(indices(4))%entry(1)/inputArray(indices(3))%entry(this%exterpCoords(2))
     end if
end function getScalingFactors
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule f_scaling_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
