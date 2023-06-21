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
module bounded_ext_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that extrapolates a quantity, optionally applying a lower and/or upper bound to the extrapolated value

    use data_kinds                   ,only: rk ,ik
    use runtime_constants            ,only: debugging, assertions
    use assertion_utility            ,only: assert, assertIdentical, assertPure
    use god_objects                  ,only: Object
    use support_types                ,only: RealArray 
    use derivation_abstract_class    ,only: Derivation
    use partition_class              ,only: Partition
    use extrapolation_abstract_class ,only: Extrapolation

    implicit none
    private

    type ,public ,extends(Derivation) :: BoundedExtDerivation
        !! Calculates extrapolated value of a fluid quantity at a boundary cell edge, optionally applying a lower and/or upper bound to the extrapolated value.
        !! Accepts scalar variables for lower/upper bounds. Expects 1-3 variables, The first expected variable name is the interpolated variable. 
        !! The second is the upper bound if no lower bound variable is expected, otherwise it is the lower bound. The third is the upper bound (if expected).
        !! Bounds are all expected to have positive values, and this will be ensured by taking their absolute values.
        !! If applied to the left boundary, will apply -abs(lowerBound) as upper bound and -abs(upperBound) as lower bound, i.e. the bounds are reflected around 0.

        real(rk) ,allocatable ,private :: fixedLowerBound !! Optional fixed lower bound value. Ignored if lower bound var expected.
        real(rk) ,allocatable ,private :: fixedUpperBound !! Optional fixed upper bound value. Ignored if upper bound var expected.

        logical ,private :: expectLowerBoundVar !! True if lower bound variable is expected
        logical ,private :: expectUpperBoundVar !! True if upper bound variable is expected 

        class(Extrapolation) ,allocatable :: extrapolationObj !! Extrapolation object used to calculate the extrapolated values

        contains

        procedure ,public :: init => initBExt

        procedure ,public :: calculate => calculateBExt

    end type BoundedExtDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initBExt(this,partitionObj,procRank,extrapolationObj,fixedLowerBound,fixedUpperBound,&
        expectLowerBoundVar,expectUpperBoundVar)
        !! Initialize bounded extrapolation derivation object

        class(BoundedExtDerivation)   ,intent(inout)  :: this
        type(Partition)               ,intent(in)     :: partitionObj !! Partition object used to determine local processor grid chunk
        integer(ik)                   ,intent(in)     :: procRank !! Current processor rank
        class(Extrapolation)          ,intent(in)     :: extrapolationObj !! Extrapolation object used to calculate the extrapolated values
        real(rk) ,optional            ,intent(in)     :: fixedLowerBound !! Optional fixed lower bound value
        real(rk) ,optional            ,intent(in)     :: fixedUpperBound !! Optional fixed upper bound value
        logical  ,optional            ,intent(in)     :: expectLowerBoundVar !! True if lower bound variable is expected. Defaults to false. 
        logical  ,optional            ,intent(in)     :: expectUpperBoundVar !! True if lower bound variable is expected. Defaults to false. 

    end subroutine initBExt  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateBExt(this,inputArray,indices) result(output)

        class(BoundedExtDerivation)         ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:)  ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                 :: output

    end function calculateBExt
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module bounded_ext_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 