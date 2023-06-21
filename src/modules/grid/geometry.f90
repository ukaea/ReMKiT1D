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
module geometry_class
    !! author: Stefan Mijin 
    !!
    !! Houses geometry object responsible for storing cell widths and Jacobians

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure

    implicit none
    private

    type ,public ,extends(Object) :: Geometry
        !! Geometry object storing spatial cell widths, cell face and cell centre jacobians, as well as linear interpolation coefficients

        real(rk) ,allocatable ,dimension(:) ,private :: cellWidths     !! x cell widths
        real(rk) ,allocatable ,dimension(:) ,private :: jacobianLeft   !! Jacobian of left cell faces
        real(rk) ,allocatable ,dimension(:) ,private :: jacobianRight  !! Jacobian of right cell faces
        real(rk) ,allocatable ,dimension(:) ,private :: jacobianCentre !! Jacobian at cell centre (avg of faces)

        real(rk) ,allocatable ,dimension(:) ,private :: linInterp      !! Linear interpolation coefficient for inner faces

        logical ,private :: periodicGrid

        contains

        procedure ,public :: getCellWidths
        procedure ,public :: getJacobianLeft
        procedure ,public :: getJacobianRight
        procedure ,public :: getJacobianCentre
        procedure ,public :: isPeriodic

        procedure ,public :: getLinInterp

        procedure ,public :: init => initGeometry

    end type Geometry
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initGeometry(this,cellWidths,jLeft,jRight,periodicGrid) 
        !! Geometry initialization routine

        class(Geometry)           ,intent(inout)  :: this
        real(rk)    ,dimension(:) ,intent(in)     :: cellWidths !! Spatial cell widths
        real(rk)    ,dimension(:) ,intent(in)     :: jLeft !! Left face jacobian values
        real(rk)    ,dimension(:) ,intent(in)     :: jRight !! Right face jacobian values
        logical     ,optional     ,intent(in)     :: periodicGrid !! Set to true if the grid is periodic

    end subroutine initGeometry
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCellWidths (this,dualGrid,extendedBoundaryCells) result(dx)
        !! Getter for cellWidths, if dualGrid is true returns dual grid values based off of cellWidths. If extendedBoundaryCells is
        !! true will extend dual grid cells at boundaries to the boundaries themselves if not periodic (defaults to true)

        class(Geometry)                       ,intent(in) :: this
        logical ,optional                     ,intent(in) :: dualGrid
        logical ,optional                     ,intent(in) :: extendedBoundaryCells
        real(rk)   ,allocatable ,dimension(:)             :: dx

    end function getCellWidths
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getJacobianLeft (this,dualGrid) result(jLeft)
        !! Getter for jacobianLeft, if dualGrid is true returns dual grid values 

        class(Geometry)                       ,intent(in) :: this
        logical ,optional                     ,intent(in) :: dualGrid
        real(rk)   ,allocatable ,dimension(:)             :: jLeft

    end function getJacobianLeft
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getJacobianRight (this,dualGrid,extendedBoundaryCells) result(jRight)
        !! Getter for jacobianRight, if dualGrid is true returns dual grid values 

        class(Geometry)                       ,intent(in) :: this
        logical ,optional                     ,intent(in) :: dualGrid
        logical ,optional                     ,intent(in) :: extendedBoundaryCells
        real(rk)   ,allocatable ,dimension(:)             :: jRight

    end function getJacobianRight
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getJacobianCentre (this,dualGrid,extendedBoundaryCells) result(jCentre)
        !! Getter for jacobianCentre, if dualGrid is true returns dual grid values 

        class(Geometry)                       ,intent(in) :: this
        logical ,optional                     ,intent(in) :: dualGrid
        logical ,optional                     ,intent(in) :: extendedBoundaryCells
        real(rk)   ,allocatable ,dimension(:)             :: jCentre

    end function getJacobianCentre
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getLinInterp (this,dualGrid) result(linInterp)
        !! Getter for linInterp, if dualGrid is true returns dual grid values (0.5 everywhere)

        class(Geometry)                       ,intent(in) :: this
        logical ,optional                     ,intent(in) :: dualGrid
        real(rk)   ,allocatable ,dimension(:)             :: linInterp

    end function getLinInterp
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isPeriodic (this) result(periodic)
        !! Getter for periodicGrid

        class(Geometry) ,intent(in) :: this
        logical                     :: periodic

    end function isPeriodic
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module geometry_class
!-----------------------------------------------------------------------------------------------------------------------------------
 