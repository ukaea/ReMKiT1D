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
module lin_interpnd_class
    !! author: Stefan Mijin 
    !! 
    !! Houses multi-linear interpolation class on n-dimensional data

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use lin_interp1d_class          ,only: Interpolation1D
    use flat_nd_data_class               ,only: FlatNDData
    use support_types               ,only: IntArray

    implicit none
    private

    type ,public ,extends(Object) :: InterpolationND
        !! Linear interpolation object housing information necessary for linearly interpolating N-dimensional data using weighted averages. If interpolation points are outside of the grid will return 0 for interpolated value. NOTE: Does not support changing the number of interpolation points after the first interpolation call.

        type(Interpolation1D) ,dimension(:) ,allocatable   ,private :: interpObjs !! Individual 1D interpolation objects. They must correspond to the axes of the data this object should be interpolating
        integer(ik)           ,dimension(:,:) ,allocatable ,private :: firstDataIndices !! An array of shape (Nd,:) where Nd is the dimensionality of the data being interpolated containing the coordinates of the origin of each hyper-cube of points used in interpolation

        type(IntArray)        ,dimension(:)   ,allocatable ,private :: hyperCube !! Indices of individual hyper-cube points (e.g. [0,1,1,1,0] etc.)

        real(rk)              ,dimension(:,:)  ,allocatable ,private :: weights !! Hypercube vertex weigths (2^Nd,:) for each interpolation point
        contains

        procedure ,public :: getFirstDataIndicesForDim
        procedure ,public :: getInterpWeightsForDim
        procedure ,public :: getInterpPointsForDim

        procedure ,public :: updateInterpolationPoints

        procedure ,public :: interpolate

        procedure ,public :: init => initInterpolation

    end type InterpolationND
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initInterpolation(this,interpolationObjs) 
            !! Initialization routine for ND interpolation object

            class(InterpolationND)           ,intent(inout)  :: this
            type(Interpolation1D) ,dimension(:) ,intent(in)  :: interpolationObjs

        end subroutine initInterpolation
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine updateInterpolationPoints(this,interpolationPoints) 
            !! Update the interpolation points and weights of all individual 1D interpolation objects and the weights/data indices used

            class(InterpolationND)           ,intent(inout)  :: this
            real(rk) ,dimension(:,:)         ,intent(in)     :: interpolationPoints

        end subroutine updateInterpolationPoints
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getFirstDataIndicesForDim (this,dim) result(inds)
            !! Getter for firstDataIndex of component 1D interpolation object with index dim
 
            class(InterpolationND)                ,intent(in) :: this
            integer(ik)                           ,intent(in) :: dim
            integer(ik)   ,allocatable ,dimension(:)          :: inds
 
     end function getFirstDataIndicesForDim
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getInterpWeightsForDim (this,dim) result(weights)
            !! Getter for interpWeights of component 1D interpolation object with index dim
 

            class(InterpolationND)                ,intent(in) :: this
            integer(ik)                           ,intent(in) :: dim
            real(rk)   ,allocatable ,dimension(:)             :: weights

        end function getInterpWeightsForDim
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getInterpPointsForDim (this,dim) result(points)
            !! Getter for interpolationPoints of component 1D interpolation object with index dim

            class(InterpolationND)                ,intent(in) :: this
            integer(ik)                           ,intent(in) :: dim
            real(rk)   ,allocatable ,dimension(:)             :: points

        end function getInterpPointsFordim
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function interpolate (this,targetArray) result(interpVals)
            !! Interpolate FlatNDData onto interpolation points using weighted averages of vertices of the containing hypercube

            class(InterpolationND)                ,intent(in) :: this
            type(FlatNDData)                      ,intent(in) :: targetArray 
            real(rk) ,allocatable ,dimension(:)               :: interpVals

        end function interpolate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module lin_interpnd_class
!-----------------------------------------------------------------------------------------------------------------------------------
 