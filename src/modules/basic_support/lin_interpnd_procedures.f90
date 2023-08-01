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
submodule (lin_interpnd_class) lin_interpnd_procedures
!! author: Stefan Mijin 
!! 
!!  Contains module procedures associated with the nd linear interpolation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initInterpolation(this,interpolationObjs) 
    !! Initialization routine for ND interpolation object

    class(InterpolationND)           ,intent(inout)  :: this
    type(Interpolation1D) ,dimension(:) ,intent(in)  :: interpolationObjs

    integer(ik) :: i,j,k

    this%interpObjs = interpolationObjs

    allocate(this%hyperCube(2**size(this%interpObjs)))
    !Generate hypercube coordinates using binary representations of 1 to 2**Nd
    do i = 1, size(this%hyperCube)
        allocate(this%hyperCube(i)%entry(size(this%interpObjs)))
        k=i
        do j = 1, size(this%interpObjs)
            this%hyperCube(i)%entry(j) = mod(k,2)
            k = k/2
        end do
    end do

    call this%makeDefined()

end subroutine initInterpolation
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine updateInterpolationPoints(this,interpolationPoints) 
    !! Update the interpolation points and weights of all individual 1D interpolation objects and the weights/data indices used

    class(InterpolationND)           ,intent(inout)  :: this
    real(rk) ,dimension(:,:)         ,intent(in)     :: interpolationPoints

    integer(ik) :: i ,j ,k

    real(rk) ,dimension(:,:) ,allocatable :: interpWeights

    if (assertions) then 
        call assertPure(this%isDefined(),"updateInterpolationPoints called on undefined InterpolationND object")
        call assertPure(size(interpolationPoints,1) == size(this%interpObjs),&
                        "interpolationPoints passed to updateInterpolationPointsRoutine of InterpolationND&
                        & objects do not conform to expected dimensionality")
    end if

    allocate(interpWeights,mold=interpolationPoints)
    if (.not. allocated(this%firstDataIndices)) allocate(this%firstDataIndices,mold=int(interpolationPoints,ik))
    if (.not. allocated(this%weights)) allocate(this%weights(2**size(this%interpObjs),size(interpolationPoints,2)))

    do i = 1,size(this%interpObjs)
        call this%interpObjs(i)%updateInterpolationPoints(interpolationPoints(i,:))
        this%firstDataIndices(i,:) = this%interpObjs(i)%getFirstDataIndices()
        interpWeights(i,:) = this%interpObjs(i)%getInterpWeights()
    end do

    this%weights = 1
    do j = 1, size(this%weights,2)
        do i = 1,size(this%weights,1)
            do k = 1, size(this%interpObjs)
                this%weights(i,j) = this%weights(i,j) * (1-interpWeights(k,j))**(1-this%hyperCube(i)%entry(k))&
                                                        *interpWeights(k,j)**this%hyperCube(i)%entry(k)
            end do
        end do
    end do

end subroutine updateInterpolationPoints
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFirstDataIndicesForDim (this,dim) result(inds)
    !! Getter for firstDataIndex of component 1D interpolation object with index dim

    class(InterpolationND)                ,intent(in) :: this
    integer(ik)                           ,intent(in) :: dim
    integer(ik)   ,allocatable ,dimension(:)          :: inds

    if (assertions) then 
        call assertPure(this%isDefined(),"getFirstDataIndicesForDim called on undefined InterpolationND object")
        call assertPure(dim <= size(this%interpObjs),"dim passed to getFirstDataIndicesForDim out of range - upper")
        call assertPure(dim > 0,"dim passed to getFirstDataIndicesForDim out of range - lower")

    end if

    inds = this%interpObjs(dim)%getFirstDataIndices()

end function getFirstDataIndicesForDim
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpWeightsForDim (this,dim) result(weights)
    !! Getter for interpWeights of component 1D interpolation object with index dim


    class(InterpolationND)                ,intent(in) :: this
    integer(ik)                           ,intent(in) :: dim
    real(rk)   ,allocatable ,dimension(:)             :: weights

    if (assertions) then 
        call assertPure(this%isDefined(),"getFirstDataIndicesForDim called on undefined InterpolationND object")
        call assertPure(dim <= size(this%interpObjs),"dim passed to getFirstDataIndicesForDim out of range - upper")
        call assertPure(dim > 0,"dim passed to getFirstDataIndicesForDim out of range - lower")

    end if

    weights = this%interpObjs(dim)%getInterpWeights()

end function getInterpWeightsForDim
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpPointsForDim (this,dim) result(points)
    !! Getter for interpolationPoints of component 1D interpolation object with index dim

    class(InterpolationND)                ,intent(in) :: this
    integer(ik)                           ,intent(in) :: dim
    real(rk)   ,allocatable ,dimension(:)             :: points

    if (assertions) then 
        call assertPure(this%isDefined(),"getInterpPointsForDim called on undefined InterpolationND object")
        call assertPure(dim <= size(this%interpObjs),"dim passed to getInterpPointsForDim out of range - upper")
        call assertPure(dim > 0,"dim passed to getInterpPointsForDim out of range - lower")

    end if

    points = this%interpObjs(dim)%getInterpPoints()

end function getInterpPointsForDim
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function interpolate (this,targetArray) result(interpVals)
    !! Interpolate FlatNDData onto interpolation points using weighted averages of vertices of the containing hypercube

    class(InterpolationND)                ,intent(in) :: this
    type(FlatNDData)                      ,intent(in) :: targetArray 
    real(rk) ,allocatable ,dimension(:)               :: interpVals

    integer(ik) :: i ,j

    allocate(interpVals(size(this%weights,2)))
    interpVals = 0
    do i = 1, size(interpVals)
        do j = 1,size(this%hyperCube)
            interpVals(i) = interpVals(i) &
                          + this%weights(j,i) * targetArray%getValue(this%firstDataIndices(:,i)+this%hyperCube(j)%entry)
        end do
    end do  
end function interpolate
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule lin_interpnd_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
