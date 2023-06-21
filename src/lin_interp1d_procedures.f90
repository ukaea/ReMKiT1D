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
submodule (lin_interp1d_class) lin_interp1d_procedures
!! author: Stefan Mijin 
!! 
!!  Contains module procedures associated with the 1d linear interpolation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initInterpolation(this,gridPoints,interpolationPoints) 
    !! Initialization routine for 1D interpolation object

    class(Interpolation1D)           ,intent(inout)  :: this
    real(rk) ,dimension(:)           ,intent(in)     :: gridPoints 
    real(rk) ,dimension(:)           ,intent(in)     :: interpolationPoints

    integer(ik) :: i 

    integer(ik) ,dimension(2) :: pair

    if (assertions) then 
        call assertPure(all(interpolationPoints >= gridPoints(1)),"Attempted to construct interpolation object when one or more&
        & interpolation points fall outside supplied grid (lower bound)")
        call assertPure(all(interpolationPoints <= gridPoints(size(gridPoints))),"Attempted to construct interpolation object &
        &when one or more interpolation points fall outside supplied grid (upper bound)")
    end if

    this%interpPoints = interpolationPoints 
    allocate(this%interpWeights(size(this%interpPoints)))
    allocate(this%firstDataIndex(size(this%interpPoints)))

    do i = 1, size(this%interpPoints)
        pair = findNearestPointsInArray(gridPoints,interpolationPoints(i))
        this%firstDataIndex(i) = pair(1)
        this%interpWeights(i) = (interpolationPoints(i) - gridPoints(pair(1)))/(gridPoints(pair(2)) - gridPoints(pair(1)))
    end do

    call this%makeDefined()

end subroutine initInterpolation
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFirstDataIndices (this) result(inds)
    !! Getter for firstDataIndex

    class(Interpolation1D)                ,intent(in) :: this
    integer(ik)    ,allocatable ,dimension(:)         :: inds

    if (assertions) call assertPure(this%isDefined(),"getFirstDataIndices called for undefined interpolation object")

    inds = this%firstDataIndex

end function getFirstDataIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpWeights (this) result(weights)
    !! Getter for interpWeights

    class(Interpolation1D)                ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)             :: weights

    if (assertions) call assertPure(this%isDefined(),"getInterpWeights called for undefined interpolation object")

    weights = this%interpWeights 

end function getInterpWeights
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpPoints (this) result(points)
    !! Getter for interpPoints

    class(Interpolation1D)                ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)             :: points

    if (assertions) call assertPure(this%isDefined(),"getInterpPoints called for undefined interpolation object")

    points = this%interpPoints 

end function getInterpPoints
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function interpolate (this,targetArray) result(interpVals)
    !! Interpolate input array defined on the same grid as the interpolation object at points 

    class(Interpolation1D)                ,intent(in) :: this
    real(rk) ,dimension(:)                ,intent(in) :: targetArray 
    real(rk)   ,allocatable ,dimension(:)             :: interpVals

    if (assertions) call assertPure(this%isDefined(),"interpolate called for undefined interpolation object")

    interpVals = targetArray(this%firstDataIndex) * (real(1,kind=rk) - this%interpWeights) &
                + this%interpWeights*targetArray(this%firstDataIndex + 1)
end function interpolate
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule lin_interp1d_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
