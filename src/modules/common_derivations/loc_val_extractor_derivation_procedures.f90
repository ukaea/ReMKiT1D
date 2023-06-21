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
submodule (loc_val_extractor_derivation_class) loc_val_extractor_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the local value extractor derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initLocValDeriv(this,partObj,numProc,targetX)
    !! Initialize harmonic extractor derivation

    class(LocValExtractorDerivation)    ,intent(inout) :: this
    type(Partition)                     ,intent(in)    :: partObj
    integer(ik)                         ,intent(in)    :: numProc
    integer(ik)                         ,intent(in)    :: targetX !! x location to extract from - global indexing

    if (assertions) call assert(partObj%isDefined(),"Undefined VSpace object passed to initLocValDeriv")

    this%isActive = (partObj%getMinXAtInd(numProc+1) <= targetX) .and. (partObj%getMaxXAtInd(numProc+1) >= targetX) 

    this%targetX = targetX - partObj%getMinXAtInd(numProc+1) + 1

    this%locNumX = partObj%getMaxXAtInd(numProc+1) - partObj%getMinXAtInd(numProc+1) + 1

    call this%makeDefined()

end subroutine initLocValDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateLocVal(this,inputArray,indices) result(output)

    class(LocValExtractorDerivation)   ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: lboundInput, inferredHaloWidth

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateLocVal called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateLocVal must be 1")
        call assertPure(all(indices>0),"indices passed to calculateLocVal out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateLocVal out of bounds - upper")

    end if

    allocate(output(1))
    output = 0
    lboundInput = lbound(inputArray(indices(1))%entry,1)
    inferredHaloWidth = (size(inputArray(indices(1))%entry) - this%locNumX)/2
    if (this%isActive) output = inputArray(indices(1))%entry(lboundInput+inferredHaloWidth+this%targetX-1)

end function calculateLocVal
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule loc_val_extractor_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
