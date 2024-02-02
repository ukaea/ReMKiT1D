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
submodule (lin_interpnd_derivation_class) lin_interpnd_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the multidimensional linear interpolation derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initInterpDeriv(this,interpObj,data)
    !! Initialize n-linear interpolation derivation object

    class(NDInterpolationDerivation)   ,intent(inout) :: this
    type(InterpolationND)              ,intent(in)    :: interpObj
    type(FlatNDData)                   ,intent(in)    :: data

    if (assertions) then
         call assert(interpObj%isDefined(),"Undefined InterpolationND object passed to derivation constructor")
         call assert(data%isDefined(),"Undefined FlatNDData object passed to derivation constructor")
    end if 

    this%interpObj = interpObj
    this%data = data
    call this%makeDefined()

end subroutine initInterpDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateInterp(this,inputArray,indices) result(output)

    class(NDInterpolationDerivation)    ,intent(inout) :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,varLength

    real(rk) ,dimension(:,:) ,allocatable :: interpolationPoints 

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateInterp called on undefined NDInterpolationDerivation object")
        call assertPure(size(indices) == size(this%data%getDims()),&
                        "The number of indices passed to calculateInterp on &
                        &NDInteprolationDerivation must be the same as the interpolation data dimensionality")

        do i = 1, size(indices)
            call assertPure(size(inputArray(indices(i))%entry) == size(inputArray(indices(1))%entry),"All indices passed to &
            &calculateInterp on NDInterpolationDerivation must correspond to variables with the same length")
        end do
    end if

    varLength = size(inputArray(indices(1))%entry)
    allocate(interpolationPoints(size(indices),varLength))

    do i = 1, size(indices)
        interpolationPoints(i,:) = inputArray(indices(i))%entry
    end do

    call this%interpObj%updateInterpolationPoints(interpolationPoints)
    output = this%interpObj%interpolate(this%data)

end function calculateInterp
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule lin_interpnd_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
