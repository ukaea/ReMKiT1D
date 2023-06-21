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
module stencil1d_class
    !! author: Stefan Mijin
    !!
    !! Houses 1D stencil component class for matrix terms

    use data_kinds                     ,only: rk, ik
    use runtime_constants              ,only: debugging, assertions
    use god_objects                    ,only: Object
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use support_types                  ,only: IntArray
    use grid_class                     ,only: Grid
    use support_functions              ,only: withinBounds
    implicit none
    private
   

    type ,public ,extends(Object) :: Stencil1D
        !! 1D stencil component for global stencil construction

        integer(ik) ,allocatable ,dimension(:) ,private :: rawStencil 

        type(IntArray) ,allocatable ,dimension(:) ,private :: fixedStencil !! Optional fixed vStencil 

        contains

        procedure ,public :: getMask
        procedure ,public :: getFixedStencil
        procedure ,public :: getStencilDims
        procedure ,public :: mapCoords
        procedure ,public :: isStencilFixed
        procedure ,public :: init => initStencil

    end type Stencil1D
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initStencil(this,rawStencil,fixedStencil) 
        !! Stencil1D object initialization routine

        class(Stencil1D)                    ,intent(inout)  :: this
        integer(ik) ,optional ,dimension(:) ,intent(in)     :: rawStencil   !! Optional raw stencil, defaults to [0] - a diagonal stencil.
        type(IntArray) ,optional ,dimension(:) ,intent(in)  :: fixedStencil !! Optional fixed stencil 

    end subroutine initStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function mapCoords(this,inputCoord,dimSize,periodic) result(output)
    !! Stencil1D coordinate mapping routine

    class(Stencil1D)       ,intent(in) :: this
    integer(ik)            ,intent(in) :: inputCoord  !! Input coordinate value
    integer(ik)            ,intent(in) :: dimSize !! Size of dimension in which mapping is done
    logical      ,optional ,intent(in) :: periodic !! True if dimension is periodic
    integer(ik) ,allocatable ,dimension(:)       :: output 

    end function mapCoords
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getMask(this,coord,dimSize,periodic) result(res)
        !! Get logical mask for included stencil points for given coordinate and dimension size. If periodic, the stencil dimension is treated as being periodic with periodicity equal to dimSize.

        class(Stencil1D)                    ,intent(in) :: this
        integer(ik)                         ,intent(in) :: coord  
        integer(ik)                         ,intent(in) :: dimSize
        logical      ,optional              ,intent(in) :: periodic
        logical ,allocatable ,dimension(:)              :: res 

    end function getMask
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getStencilDims(this) result(dim)
        !! Returns size of stencil

        class(Stencil1D)       ,intent(in) :: this
        integer(ik)                        :: dim 

    end function getStencilDims
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isStencilFixed(this) result(stencilIsFixed)
        !! Check if stencil is fixed

        class(Stencil1D) ,intent(in) :: this
        logical                       :: stencilIsFixed 

    end function isStencilFixed
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getFixedStencil(this) result(fixedStencil)
        !! Return values of fixed stencil

        class(Stencil1D)                          ,intent(in) :: this
        type(IntArray) ,allocatable ,dimension(:)             :: fixedStencil 

    end function getFixedStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module stencil1d_class
!-----------------------------------------------------------------------------------------------------------------------------------
 