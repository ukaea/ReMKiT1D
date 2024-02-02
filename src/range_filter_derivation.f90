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
module range_filter_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses composite derivation class where a derivations result is zeroed out wherever a set of passed variables is not within their
    !! defined ranges 

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray ,IntArray 
    use derivation_abstract_class   ,only: Derivation
    use support_functions           ,only: findIndices

    implicit none
    private

    type ,public ,extends(Derivation) :: RangeFilterDerivation
        !! Composite derivation class containing a single calculation rule, which is overriden with zeros wherever any one in a set of control
        !! variables is outside their defined range. If derivIndices isn't passed will pass all indices to the derivation.

        class(Derivation) ,allocatable ,private :: derivObj !! Derivation providing default values

        integer(ik) ,allocatable ,dimension(:) ,private :: derivIndices !! Subset of index vector passed to  derivation
        integer(ik) ,allocatable ,dimension(:) ,private :: controlIndices !! Subset of index vector used to determine which elements should be zeroed out

        type(RealArray) ,allocatable ,dimension(:) ,private :: controlRanges !! Array of size 2 vectors representing the lower and upper bound for each control value

        contains

        procedure ,public :: init => initRangeFilterDeriv

        procedure ,public :: calculate => calculateRangeFilter

    end type RangeFilterDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initRangeFilterDeriv(this,derivObj,controlIndices,controlRanges,derivIndices)
        !! Initialize range filter derivation object

        class(RangeFilterDerivation)        ,intent(inout) :: this
        class(Derivation)                   ,intent(in)    :: derivObj
        integer(ik) ,dimension(:)           ,intent(in)    :: controlIndices
        type(RealArray) ,dimension(:)       ,intent(in)    :: controlRanges
        integer(ik) ,optional ,dimension(:) ,intent(in)    :: derivIndices

    end subroutine initRangeFilterDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateRangeFilter(this,inputArray,indices) result(output)

        class(RangeFilterDerivation)       ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateRangeFilter
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module range_filter_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 