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
    !! Initialization routine for 1D Interpolation1D object

    class(Interpolation1D)           ,intent(inout)  :: this
    real(rk) ,dimension(:)           ,intent(in)     :: gridPoints 
    real(rk) ,dimension(:) ,optional ,intent(in)     :: interpolationPoints

    if (assertions) then 
        if (present(interpolationPoints)) then 
        call assertPure(all(interpolationPoints >= gridPoints(1)),"Attempted to construct Interpolation1D object when one or more&
        & interpolation points fall outside supplied grid (lower bound)")
        call assertPure(all(interpolationPoints <= gridPoints(size(gridPoints))),"Attempted to construct Interpolation1D object &
        &when one or more interpolation points fall outside supplied grid (upper bound)")
        end if
    end if

    this%gridBuffer = gridPoints

    if (present(interpolationPoints)) call this%updateInterpolationPoints(interpolationPoints)

    call this%makeDefined()

end subroutine initInterpolation
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine updateInterpolationPoints(this,interpolationPoints) 
    !! Updates the interpolation points and weightss

    class(Interpolation1D)           ,intent(inout)  :: this
    real(rk) ,dimension(:)           ,intent(in)     :: interpolationPoints

    integer(ik) :: i 

    integer(ik) ,dimension(2) :: pair

    this%interpPoints = interpolationPoints 
    if (.not. allocated(this%interpWeights)) then 
        allocate(this%interpWeights(size(this%interpPoints)))
        allocate(this%firstDataIndex(size(this%interpPoints)))
    else if (size(this%interpWeights) .ne. size(this%interpPoints)) then
        deallocate(this%interpWeights)
        deallocate(this%firstDataIndex)

        allocate(this%interpWeights(size(this%interpPoints)))
        allocate(this%firstDataIndex(size(this%interpPoints)))

    end if

    do i = 1, size(this%interpPoints)
        pair = findNearestPointsInArray(this%gridBuffer,interpolationPoints(i))
        this%firstDataIndex(i) = pair(1)
        this%interpWeights(i) = (interpolationPoints(i) - this%gridBuffer(pair(1)))&
                                /(this%gridBuffer(pair(2)) - this%gridBuffer(pair(1)))
    end do

end subroutine updateInterpolationPoints
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFirstDataIndices (this) result(inds)
    !! Getter for firstDataIndex

    class(Interpolation1D)                ,intent(in) :: this
    integer(ik)    ,allocatable ,dimension(:)         :: inds

    if (assertions) then 
        call assertPure(this%isDefined(),"getFirstDataIndices called for undefined Interpolation1D object")
        call assertPure(allocated(this%firstDataIndex),&
                        "getFirstDataIndices called before interpolation points defined for Interpolation1D object")
    end if
    inds = this%firstDataIndex

end function getFirstDataIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpWeights (this) result(weights)
    !! Getter for interpWeights

    class(Interpolation1D)                ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)             :: weights

    if (assertions) then 
        call assertPure(this%isDefined(),"getInterpWeights called for undefined Interpolation1D object")
        call assertPure(allocated(this%interpWeights),&
                        "getInterpWeights called before interpolation points defined for Interpolation1D object")
    end if

    weights = this%interpWeights 

end function getInterpWeights
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpPoints (this) result(points)
    !! Getter for interpPoints

    class(Interpolation1D)                ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)             :: points

    if (assertions) then 
        call assertPure(this%isDefined(),"getInterpPoints called for undefined Interpolation1D object")
        call assertPure(allocated(this%InterpPoints),&
                        "getInterpPoints called before interpolation points defined for Interpolation1D object")
    end if

    points = this%interpPoints 

end function getInterpPoints
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function interpolate (this,targetArray) result(interpVals)
    !! Interpolate input array defined on the same grid as the Interpolation1D object at points 

    class(Interpolation1D)                ,intent(in) :: this
    real(rk) ,dimension(:)                ,intent(in) :: targetArray 
    real(rk)   ,allocatable ,dimension(:)             :: interpVals

    if (assertions) then 
        call assertPure(this%isDefined(),"interpolate called for undefined Interpolation1D object")
        call assertPure(allocated(this%InterpPoints),&
                        "interpolate called before interpolation points defined for Interpolation1D object")
    end if

    interpVals = targetArray(this%firstDataIndex) * (real(1,kind=rk) - this%interpWeights) &
                + this%interpWeights*targetArray(this%firstDataIndex + 1)
end function interpolate
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule lin_interp1d_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
