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
module lin_interp1D_class
    !! author: Stefan Mijin 
    !! 
    !! Houses linear interpolation class

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_functions           ,only: findNearestPointsInArray

    implicit none
    private

    type ,public ,extends(Object) :: Interpolation1D
        !! Linear interpolation object housing information necessary for linearly interpolating 1D data

        integer(ik) ,allocatable ,dimension(:) ,private :: firstDataIndex !! Array containing the first interpolation index (the second is 1 above the first)
        real(rk)    ,allocatable ,dimension(:) ,private :: interpWeights !! Array containing interpolation weights for each interpolation point
        real(rk)    ,allocatable ,dimension(:) ,private :: interpPoints !! Points at which this interpolation object provides values

        contains

        procedure ,public :: getFirstDataIndices
        procedure ,public :: getInterpWeights
        procedure ,public :: getInterpPoints

        procedure ,public :: interpolate

        procedure ,public :: init => initInterpolation

    end type Interpolation1D
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initInterpolation(this,gridPoints,interpolationPoints) 
            !! Initialization routine for 1D interpolation object

            class(Interpolation1D)           ,intent(inout)  :: this
            real(rk) ,dimension(:)           ,intent(in)     :: gridPoints 
            real(rk) ,dimension(:)           ,intent(in)     :: interpolationPoints

        end subroutine initInterpolation
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getFirstDataIndices (this) result(inds)
            !! Getter for firstDataIndex
 
            class(Interpolation1D)                ,intent(in) :: this
            integer(ik)   ,allocatable ,dimension(:)          :: inds
 
     end function getFirstDataIndices
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getInterpWeights (this) result(weights)
            !! Getter for interpWeights

            class(Interpolation1D)                ,intent(in) :: this
            real(rk)   ,allocatable ,dimension(:)             :: weights

        end function getInterpWeights
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getInterpPoints (this) result(points)
            !! Getter for interpPoints

            class(Interpolation1D)                ,intent(in) :: this
            real(rk)   ,allocatable ,dimension(:)             :: points

        end function getInterpPoints
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function interpolate (this,targetArray) result(interpVals)
            !! Interpolate input array defined on the same grid as the interpolation object at points 

            class(Interpolation1D)                ,intent(in) :: this
            real(rk) ,dimension(:)                ,intent(in) :: targetArray 
            real(rk)   ,allocatable ,dimension(:)             :: interpVals

        end function interpolate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module lin_interp1D_class
!-----------------------------------------------------------------------------------------------------------------------------------
 