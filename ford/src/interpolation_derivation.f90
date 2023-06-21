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
module interpolation_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that interpolates one quantity from the original to the dual grid or vice versa. 

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation
    use geometry_class              ,only: Geometry
    use grid_class                  ,only: Grid

    implicit none
    private

    type ,public ,extends(Derivation) :: InterpolationDerivation
        !! Derivation that interpolates a variable from the original to the dual grid, or vise versa. Distributions are interpolated
        !! by assuming that even l harmonics live on the original grid and odd harmonics on the dual grid. 

        real(rk) ,allocatable ,dimension(:) ,private :: linInterp !! Linear interpolation coefficients for inner points
        real(rk) ,allocatable ,dimension(:) ,private :: linInterpDual !! Linear interpolation coefficients for inner points on dual grid
        real(rk)                            ,private :: linExterp !! Linear extrapolation coefficients for the right boundary point
        real(rk)                            ,private :: linExterpRDual !! Linear extrapolation coefficient for the right boundary point during inverse interpolation on non-periodic grid
        real(rk)                            ,private :: linExterpLDual !! Linear extrapolation coefficient for the left boundary point during inverse interpolation on non-periodic grid


        logical ,private :: periodicGrid
        logical ,private :: inverseInterp !! True if interpolating from dual to original grid. Defaults to false
        logical ,private :: distInterp !! True if interpolating a distribution function

        logical ,private :: containsLeftBoundary !! True if the local x-grid contains the left boundary
        logical ,private :: containsRightBoundary !! True if the local x-grid contains the right boundary

        integer(ik) ,private :: locNumX !! Size of local grid chunk (without any halos)
        integer(ik) ,private :: numV
        integer(ik) ,private :: numH 

        logical(ik) ,allocatable ,dimension(:) :: oddHarmonic !! True for odd l harmonics


        contains

        procedure ,public :: init => initInterpDeriv

        procedure ,public :: calculate => calculateInterp

    end type InterpolationDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initInterpDeriv(this,geometryObj,gridObj,minX,maxX,inverseInterp,distInterp)
        !! Initialize interpolation derivation object

        class(InterpolationDerivation)   ,intent(inout) :: this
        type(Geometry)                   ,intent(in)    :: geometryObj
        type(Grid)                       ,intent(in)    :: gridObj
        integer(ik)                      ,intent(in)    :: minX 
        integer(ik)                      ,intent(in)    :: maxX
        logical ,optional                ,intent(in)    :: inverseInterp
        logical ,optional                ,intent(in)    :: distInterp

    end subroutine initInterpDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateInterp(this,inputArray,indices) result(output)

        class(InterpolationDerivation)     ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateInterp
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module interpolation_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 