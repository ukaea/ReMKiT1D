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
module stencil_class
    !! author: Stefan Mijin
    !!
    !! Houses stencil class which provides stencil-based coordinate mapping for matrix terms

    use data_kinds                     ,only: rk, ik
    use runtime_constants              ,only: debugging, assertions
    use god_objects                    ,only: Object
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use support_types                  ,only: IntArray
    use grid_class                     ,only: Grid
    use support_functions              ,only: withinBounds
    use stencil1d_class                ,only: Stencil1D
    implicit none
    private
   

    type ,public ,extends(Object) :: Stencil
        !! Stencil object for construction of matrix patterns by mapping row coordinates to sets of column coordinates 

        type(Stencil1D) ,private :: xStencil !! Stencil in x-direction
        type(Stencil1D) ,allocatable ,private :: vStencil !! Stencil in v-direction (used as absolute coordinates if mapping from just x-coordinate)
        type(Stencil1D) ,allocatable ,private :: hStencil !! Stencil in total harmonics (used as absolute coordinates if mapping from just x-coordinate)

        logical                                ,private :: xPeriodic !! True if the x grid is periodic

        contains

        procedure ,public :: mapCoords 

        procedure ,public :: init => initStencil

    end type Stencil
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initStencil(this,xStencil,vStencil,hStencil,mapToDist,xPeriodic,vStencilFixed) 
        !! Stencil object initialization routine

        class(stencil)                      ,intent(inout)  :: this
        integer(ik) ,optional ,dimension(:) ,intent(in)     :: xStencil !! Stencil in x-direction
        integer(ik) ,optional ,dimension(:) ,intent(in)     :: vStencil !! Stencil in v-direction 
        integer(ik) ,optional ,dimension(:) ,intent(in)     :: hStencil !! l harmonic stencil
        logical     ,optional               ,intent(in)     :: mapToDist !! Set to true if stencil maps to harmonic/velocity space
        logical     ,optional               ,intent(in)     :: xPeriodic !! Set to true if x-grid is periodic
        type(IntArray) ,optional ,dimension(:) ,intent(in)  :: vStencilFixed !! Optional fixed vStencil 

    end subroutine initStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function mapCoords(this,gridObj,inputArray) result(output)
        !! Stencil coordinate mapping routine

        class(stencil)                      ,intent(in) :: this
        type(Grid)                          ,intent(in) :: gridObj !! Grid used to construct the mapping
        integer(ik)    ,dimension(:)        ,intent(in) :: inputArray  !! Input array of coordinates (size 1 - [x], or size 3 [x,h,v])
        type(IntArray) ,allocatable ,dimension(:)       :: output 

    end function mapCoords
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module stencil_class
!-----------------------------------------------------------------------------------------------------------------------------------
 