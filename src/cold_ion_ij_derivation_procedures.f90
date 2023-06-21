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
submodule (cold_ion_ij_int_derivation_class) cold_ion_ij_int_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the cold ion Shkarofsky I/J integral derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initColdIJInt(this,gridObj,ind,isJInt)
    !! Initialize col dion Shkarofsky I/J integral derivation derivation

    class(ColdIonIJIntDerivation)     ,intent(inout) :: this
    type(Grid)                        ,intent(in)    :: gridObj
    integer(ik)                       ,intent(in)    :: ind !! Index of indegral
    logical ,optional                 ,intent(in)    :: isJInt !! If true the lower triangular J integral is calculated instead of the I integral. Defaults to false.


    if (assertions) call assert(gridObj%isDefined(),"Undefined VSpace object passed to initColdIJInt")
    
    this%isJInt = .false. 
    if (present(isJInt)) this%isJInt = isJInt 

    this%ind = ind 
    this%lGridCopy = gridObj%getLGrid()
    this%vGridCopy = gridObj%getVGrid()

    call this%makeDefined()
    
end subroutine initColdIJInt  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateColdIJInt(this,inputArray,indices) result(output)

    class(ColdIonIJIntDerivation)      ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,j ,k ,inferredNumX ,lboundInput ,offset ,lNum

    logical ,allocatable ,dimension(:) :: heavisideStep


    if (assertions) then 

        call assertPure(this%isDefined(),"calculateColdIJInt called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateColdIJInt must be 1")
        call assertPure(all(indices>0),"indices passed to calculateColdIJInt out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateColdIJInt out of bounds - upper")

    end if

    inferredNumX = size(inputArray(indices(1))%entry)

    allocate(output(inferredNumX*size(this%vGridCopy)*size(this%lGridCopy)))
    output = 0
    lBoundInput = lbound(inputArray(indices(1))%entry,1)

    allocate(heavisideStep(size(this%vGridCopy)))

    do i = 1, inferredNumX
        if (this%isJInt) then 
            heavisideStep = this%vGridCopy < abs(inputArray(indices(1))%entry(i+lboundInput-1))
        else
            heavisideStep = this%vGridCopy > abs(inputArray(indices(1))%entry(i+lboundInput-1))

        end if
        do j = 1,size(this%lGridCopy)
            
            offset = (i-1)*size(this%lGridCopy)*size(this%vGridCopy) + (j-1)*size(this%vGridCopy)
            lNum = this%lGridCopy(j)
            do k = 1,size(this%vGridCopy)
                if (heavisideStep(k)) &
                output(offset+k) = (2*lNum+1)*sign(real(1,kind=rk),inputArray(indices(1))%entry(i+lboundInput-1))**lNum &
                                 * abs(inputArray(indices(1))%entry(i+lboundInput-1)) ** this%ind / this%vGridCopy(k)**this%ind 
            end do
        end do
    end do

end function calculateColdIJInt
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule cold_ion_ij_int_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
