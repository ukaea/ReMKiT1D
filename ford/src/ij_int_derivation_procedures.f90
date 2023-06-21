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
    !! Initialize Shkarofsky I/J integral derivation derivation

    class(IJIntDerivation)     ,intent(inout) :: this
    type(VSpace)               ,intent(in)    :: vSpaceObj
    integer(ik)                ,intent(in)    :: ind !! Index of indegral
    logical ,optional          ,intent(in)    :: isJInt !! If true the lower triangular J integral is calculated instead of the I integral. Defaults to false.

    logical :: jInt 
    if (assertions) call assert(vSpaceObj%isDefined(),"Undefined VSpace object passed to initIJInt")

    this%numV = vSpaceObj%getNumV()
    
    jInt = .false. 
    if (present(isJInt)) jInt = isJInt 

    if (jInt) then 
        this%ijMatBuffer = vSpaceObj%getShkarofskyJMat(ind)
    else
        this%ijMatBuffer = vSpaceObj%getShkarofskyIMat(ind)
    end if
    call this%makeDefined()
    
end subroutine initIJInt  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateIJInt(this,inputArray,indices) result(output)

    class(IJIntDerivation)             ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,j ,inferredHarmonics ,indOffset ,lboundInput

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateIJInt called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateIJInt must be 1")
        call assertPure(all(indices>0),"indices passed to calculateIJInt out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateIJInt out of bounds - upper")

    end if

    inferredHarmonics = size(inputArray(indices(1))%entry)/this%numV

    allocate(output(inferredHarmonics*this%numV))
    output = 0
    lBoundInput = lbound(inputArray(indices(1))%entry,1)

    do i = 1, inferredHarmonics
        indOffset = (i-1)*this%numV + lboundInput -1
        do j = 1,size(this%ijMatBuffer%values)
            output((i-1)*this%numV+this%ijMatBuffer%rowIndex(j)) = dot_product(this%ijMatBuffer%values(j)%entry &
                                                            ,inputArray(1)%entry(indOffset+this%ijMatBuffer%columnVector(j)%entry))
        end do
    end do

end function calculateIJInt
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule ij_int_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
