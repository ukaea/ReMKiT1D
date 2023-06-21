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
submodule (harmonic_extractor_derivation_class) harmonic_extractor_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the harmonic extractor derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initHEDerivation(this,vSpaceObj,targetH)
    !! Initialize harmonic extractor derivation

    class(HExtractorDerivation)    ,intent(inout) :: this
    type(VSpace)                   ,intent(in)    :: vSpaceObj
    integer(ik)                    ,intent(in)    :: targetH !! Harmonic to extract

    if (assertions) call assert(vSpaceObj%isDefined(),"Undefined VSpace object passed to initHEDerivation")

    this%numH = vSpaceObj%getNumH()
    this%numV = vSpaceObj%getNumV()

    this%targetH = targetH

    call this%makeDefined()

end subroutine initHEDerivation  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateHE(this,inputArray,indices) result(output)

    class(HExtractorDerivation)        ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,inferredNumX ,indOffset ,lboundInput

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateHE called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateHE must be 1")
        call assertPure(all(indices>0),"indices passed to calculateHE out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateHE out of bounds - upper")

    end if

    inferredNumX = size(inputArray(indices(1))%entry)/(this%numH*this%numV)
    lBoundInput = lbound(inputArray(indices(1))%entry,1)
    allocate(output(inferredNumX*this%numV))

    do i = 1, inferredNumX 
        indOffset = (i-1)*this%numV*this%numH + (this%targetH-1)*this%numV + lboundInput - 1
        output((i-1)*this%numV+1:i*this%numV) = inputArray(indices(1))%entry(indOffset+1:indOffset+this%numV)
    end do

end function calculateHE
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule harmonic_extractor_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
