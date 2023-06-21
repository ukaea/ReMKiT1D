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
submodule (range_filter_derivation_class) range_filter_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the range filter derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initRangeFilterDeriv(this,derivObj,controlIndices,controlRanges,derivIndices)
        !! Initialize range filter derivation object

        class(RangeFilterDerivation)        ,intent(inout) :: this
        class(Derivation)                   ,intent(in)    :: derivObj
        integer(ik) ,dimension(:)           ,intent(in)    :: controlIndices
        type(RealArray) ,dimension(:)       ,intent(in)    :: controlRanges
        integer(ik) ,optional ,dimension(:) ,intent(in)    :: derivIndices

        integer(ik) :: i

        if (assertions) then 
            call assert(derivObj%isDefined(),"derivObj passed to initRangeFilterDeriv not defined")
            call assert(size(controlIndices) == size(controlRanges),&
            "controlRanges and controlIndices passed to initRangeFilterDeriv must conform")
            
            do i = 1,size(controlRanges)
                call assert(size(controlRanges(i)%entry)==2,"all controlRanges passed to initRangeFilterDeriv must be size 2")
                call assert(controlRanges(i)%entry(1)<controlRanges(i)%entry(2),"all controlRanges must be ordered lower to upper")
            end do
        end if

        allocate(this%derivObj,source=derivObj)

        this%controlIndices=controlIndices
        this%controlRanges=controlRanges
        
        if (present(derivIndices)) this%derivIndices=derivIndices

        call this%makeDefined()

    end subroutine initRangeFilterDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateRangeFilter(this,inputArray,indices) result(output)

        class(RangeFilterDerivation)       ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

        logical ,allocatable ,dimension(:) :: mask

        integer ,allocatable ,dimension(:) :: zeroIndices

        integer(ik) :: i
        
        if (assertions) then 

            call assertPure(this%isDefined(),"calculateRangeFilter called on undefined derivation object")
            
            call assertPure(all(indices>0),"indices passed to calculateRangeFilter out of bounds - lower")
            call assertPure(all(indices<=size(inputArray)),"indices passed to calculateRangeFilter out of bounds - upper")
    
        end if

        if (allocated(this%derivIndices)) then
            output = this%derivObj%calculate(inputArray,indices(this%derivIndices))
        else
            output = this%derivObj%calculate(inputArray,indices)
        end if

        allocate(mask(size(output)))

        mask = .true.

        do i = 1,size(this%controlRanges)
            mask = mask .and. inputArray(indices(this%controlIndices(i)))%entry > this%controlRanges(i)%entry(1) .and. &
                   inputArray(indices(this%controlIndices(i)))%entry < this%controlRanges(i)%entry(2)
        end do

        zeroIndices = findIndices(.not. mask)

        output(zeroIndices) = 0 

    end function calculateRangeFilter
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule range_filter_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
