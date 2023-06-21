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
submodule (bounded_ext_derivation_class) bounded_ext_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the bounded extrapolation derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initBExt(this,partitionObj,procRank,extrapolationObj,fixedLowerBound,fixedUpperBound,&
                        expectLowerBoundVar,expectUpperBoundVar)
    !! Initialize bounded extrapolation derivation object

    class(BoundedExtDerivation)   ,intent(inout)  :: this
    type(Partition)               ,intent(in)     :: partitionObj !! Partition object used to determine local processor grid chunk
    integer(ik)                   ,intent(in)     :: procRank !! Current processor rank
    class(Extrapolation)          ,intent(in)     :: extrapolationObj !! Extrapolation object used to calculate the extrapolated values
    real(rk) ,optional            ,intent(in)     :: fixedLowerBound !! Optional fixed lower bound value
    real(rk) ,optional            ,intent(in)     :: fixedUpperBound !! Optional fixed upper bound value
    logical  ,optional            ,intent(in)     :: expectLowerBoundVar !! True if lower bound variable is expected. Defaults to false. 
    logical  ,optional            ,intent(in)     :: expectUpperBoundVar !! True if lower bound variable is expected. Defaults to false. 

    logical :: sVars

    integer(ik) :: minX ,maxX, inferredGridSize

    if (assertions) call assert(partitionObj%isDefined(),"Partition passed to initBExt not defined")

    if (present(fixedLowerBound)) this%fixedLowerBound = fixedLowerBound
    if (present(fixedUpperBound)) this%fixedUpperBound = fixedUpperBound

    if (allocated(this%fixedLowerBound) .and. (allocated(this%fixedUpperBound))) &
    call assert(abs(this%fixedUpperBound)>abs(this%fixedLowerBound),&
    "Fixed upper bound passed to initBExt greater than passed lower bound")

    this%expectLowerBoundVar = .false. 
    this%expectUpperBoundVar = .false. 

    if (present(expectLowerBoundVar)) this%expectLowerBoundVar = expectLowerBoundVar
    if (present(expectUpperBoundVar)) this%expectUpperBoundVar = expectUpperBoundVar

    allocate(this%extrapolationObj,source=extrapolationObj)

    call this%makeDefined()

end subroutine initBExt  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateBExt(this,inputArray,indices) result(output)

    class(BoundedExtDerivation)         ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:)  ,intent(in)    :: indices
    real(rk) ,allocatable ,dimension(:)                 :: output

    integer(ik) :: expectedNumIndices

    real(rk) :: lowerBound ,upperBound

    logical :: applyLBound ,applyUBound

    if (assertions) call assertPure(this%isDefined(),"calculateBExt called on undefined derivation object")

    expectedNumIndices = 1 
    if (this%expectLowerBoundVar) expectedNumIndices = expectedNumIndices + 1
    if (this%expectUpperBoundVar) expectedNumIndices = expectedNumIndices + 1

    if (assertions) then 

        call assertPure(size(indices) == expectedNumIndices,"Unexpected number of indices passed to calculateBExt")
        call assertPure(all(indices>0),"indices passed to calculateBExt out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateBExt out of bounds - upper")

    end if

    allocate(output(1))
    output(1) = 0
    if (this%extrapolationObj%hasBoundary()) then
        !Perform initial linear extrapolation
        output(1) = this%extrapolationObj%extrapolate(inputArray(indices(1))%entry)
        !Check for bounds 
        applyLBound = allocated(this%fixedLowerBound) .or. this%expectLowerBoundVar
        applyUBound = allocated(this%fixedUpperBound) .or. this%expectUpperBoundVar

        if (applyLBound) then 
            if (this%expectLowerBoundVar) then 
                if (size(inputArray(indices(2))%entry) == 1) then 
                    lowerBound = inputArray(indices(2))%entry(1)
                else
                    lowerBound = this%extrapolationObj%extrapolate(inputArray(indices(2))%entry)
                end if
            else
                lowerBound = this%fixedLowerBound
            end if
        end if

        if (applyUBound) then 
            if (this%expectUpperBoundVar) then 
                if (size(inputArray(indices(expectedNumIndices))%entry) == 1) then 
                    upperBound = inputArray(indices(expectedNumIndices))%entry(1)
                else
                    upperBound = this%extrapolationObj%extrapolate(inputArray(indices(expectedNumIndices))%entry)
                end if
            else
                upperBound = this%fixedUpperBound
            end if
        end if

        if (applyLBound .and. applyUBound) &
        call assertPure(abs(lowerBound) < abs(upperBound),"Lower bound not smaller than upper bound in calculateBExt")
        
        if (this%extrapolationObj%isLeftBoundary()) then 
            if (applyLBound) output(1) = min(output(1),-abs(lowerBound))
            if (applyUBound) output(1) = max(output(1),-abs(upperBound))

        else
            if (applyLBound) output(1) = max(output(1),abs(lowerBound))
            if (applyUBound) output(1) = min(output(1),abs(upperBound))
        end if
    end if
    
end function calculateBExt
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule bounded_ext_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
