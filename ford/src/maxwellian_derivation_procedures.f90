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
submodule (maxwellian_derivation_class) maxwellian_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the Maxwellian derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMaxwellianDeriv(this,vSpaceObj)
    !! Initialize Maxwellian derivation

    class(MaxwellianDerivation)      ,intent(inout) :: this
    type(VSpace)                     ,intent(in)    :: vSpaceObj

    if (assertions) call assert(vSpaceObj%isDefined(),"Undefined VSpace passed to iniMaxwellianDeriv")

    this%vGridCopy = vSpaceObj%getVGrid()
    this%numH = vSpaceObj%getNumH()

    call this%makeDefined()

end subroutine initMaxwellianDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateMaxwellian(this,inputArray,indices) result(output)

    class(MaxwellianDerivation)        ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,inferredNumX ,numV ,lBoundInput
    real(rk) :: n, T

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateMaxwellian called on undefined derivation object")
        call assertPure(size(indices) == 2,"Number of indices passed to calculateMaxwellian must be 2")
        call assertPure(all(indices>0),"indices passed to calculateMaxwellian out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateMaxwellian out of bounds - upper")

    end if
    inferredNumX = size(inputArray(indices(1))%entry)
    numV = size(this%vGridCopy)
    allocate(output(inferredNumX*this%numH*numV))
    output = 0
    lBoundInput = lbound(inputArray(indices(1))%entry,1)
    do i = 1, inferredNumX
        T = inputArray(indices(1))%entry(lBoundInput-1+i)
        n = inputArray(indices(2))%entry(lBoundInput-1+i)
        output((i-1)*numV*this%numH+1:(i-1)*numV*this%numH+numV) = normMaxwellian(n,T,this%vGridCopy)
    end do

end function calculateMaxwellian
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule maxwellian_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
