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
module lin_interpnd_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that linearly interpolates on multidimensional data based on some input variable values. 

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation
    use lin_interpnd_class          ,only: InterpolationND
    use flat_nd_data_class          ,only: FlatNDData


    implicit none
    private

    type ,public ,extends(Derivation) :: NDInterpolationDerivation
        !! Interpolates n-dimensional data linearly based on n-input variables (must all be the same length)

        type(InterpolationND) ,allocatable :: interpObj !! n-dimensional interpolation object respnsible for interpolation 
        type(FlatNDData)      ,allocatable :: data !! Data to interpolate on

        contains

        procedure ,public :: init => initInterpDeriv

        procedure ,public :: calculate => calculateInterp

    end type NDInterpolationDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initInterpDeriv(this,interpObj,data)
        !! Initialize n-linear interpolation derivation object

        class(NDInterpolationDerivation)   ,intent(inout) :: this
        type(InterpolationND)              ,intent(in)    :: interpObj
        type(FlatNDData)                   ,intent(in)    :: data

    end subroutine initInterpDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateInterp(this,inputArray,indices) result(output)

        class(NDInterpolationDerivation)    ,intent(inout) :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateInterp
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module lin_interpnd_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 