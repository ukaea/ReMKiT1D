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
module grid_class
    !! author: Stefan Mijin 
    !!
    !!  Houses Grid object responsible for storing spatial, velocity, and harmonic grid vertices


    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_functions           ,only: findIndices

    implicit none
    private

    type ,public ,extends(Object) :: Grid
        !! Grid object storing cell positions of x and v grid, as well as Legendre harmonic data of the grid

        real(rk) ,allocatable ,dimension(:) ,private :: xGrid !! Positions of x-grid cell centres
        real(rk) ,allocatable ,dimension(:) ,private :: vGrid !! Positions of v-grid cell centres
 
        integer(ik) ,allocatable ,dimension(:) ,private :: lGrid !! l-numbers of each resolved harmonic
        integer(ik) ,allocatable ,dimension(:) ,private :: mGrid !! m-numbers of each resolved harmonic
        logical     ,allocatable ,dimension(:) ,private :: imaginaryHarmonic !! True if resolved harmonic is imaginary

        integer(ik) ,private :: maxL !! Highest l-number resolved on grid
        integer(ik) ,private :: maxM !! Highest m-number resolved on grid

        contains

        procedure ,public :: getXGrid
        procedure ,public :: getVGrid

        procedure ,public :: getMaxL
        procedure ,public :: getMaxM
        procedure ,public :: getNumX
        procedure ,public :: getNumV
        procedure ,public :: getNumH

        procedure ,public :: getLGrid
        procedure ,public :: getMGrid
        procedure ,public :: getHarmonicIm
        procedure ,public :: getH
        procedure ,public :: getL 
        procedure ,public :: getM
        procedure ,public :: isImaginary

        procedure ,public :: init => initGrid

    end type Grid
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initGrid(this,x,v,maxL,maxM) 
            !! Grid initialization routine

            class(Grid)           ,intent(inout)  :: this
            real(rk)    ,dimension(:) ,intent(in) :: x !! Positions of x-grid cell centres
            real(rk)    ,dimension(:) ,intent(in) :: v !! Positions of v-grid cell centres
            integer(ik)               ,intent(in) :: maxL !! Highest resolved l-harmonic
            integer(ik)               ,intent(in) :: maxM !! Highest resolved m-harmonic

        end subroutine initGrid
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getXGrid (this) result(x)
            !! Getter for xGrid

            class(Grid)                           ,intent(in) :: this
            real(rk)   ,allocatable ,dimension(:)             :: x
 
        end function getXGrid
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getVGrid (this) result(v)
            !! Getter for vGrid

            class(Grid)                           ,intent(in) :: this
            real(rk)   ,allocatable ,dimension(:)             :: v
 
        end function getVGrid
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMaxL (this) result(maxL)
            !! Getter for maxL

            class(Grid) ,intent(in) :: this
            integer(ik)             :: maxL
 
        end function getMaxL
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMaxM (this) result(maxM)
            !! Getter for maxM

            class(Grid) ,intent(in) :: this
            integer(ik)             :: maxM
 
        end function getMaxM
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumH (this) result(numH)
            !! Return total number of resolved harmonic on grid

            class(Grid) ,intent(in) :: this
            integer(ik)             :: numH
 
        end function getNumH
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumX (this) result(numX)
            !! Return number of x points on grid

            class(Grid) ,intent(in) :: this
            integer(ik)             :: numX
 
        end function getNumX
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumV (this) result(numV)
            !! Return number of v points on grid

            class(Grid) ,intent(in) :: this
            integer(ik)             :: numV
 
        end function getNumV
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getLGrid (this) result(l)
            !! Getter for lGrid

            class(Grid)                           ,intent(in) :: this
            integer(ik)   ,allocatable ,dimension(:)          :: l
 
        end function getLGrid
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMGrid (this) result(m)
            !! Getter for mGrid

            class(Grid)                           ,intent(in) :: this
            integer(ik)   ,allocatable ,dimension(:)          :: m
 
        end function getMGrid
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getHarmonicIm (this) result(im)
            !! Getter for imaginaryHarmonic

            class(Grid)                       ,intent(in) :: this
            logical   ,allocatable ,dimension(:)          :: im
 
        end function getHarmonicIm
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getH (this,l,m,im) result(h)
            !! Return index of harmonic l,m, (if im=true returns the imaginary component)

            class(Grid)  ,intent(in) :: this
            integer(ik)  ,intent(in) :: l,m
            logical      ,intent(in) :: im
            integer(ik)              :: h
            
        end function getH
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getL (this,ind) result(l)
            !! Return l number for given harmonic index

            class(Grid)  ,intent(in) :: this
            integer(ik)  ,intent(in) :: ind
            integer(ik)              :: l
            
        end function getL
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getM (this,ind) result(m)
            !! Return m number for given harmonic index

            class(Grid)  ,intent(in) :: this
            integer(ik)  ,intent(in) :: ind
            integer(ik)              :: m
            
        end function getM
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function isImaginary (this,ind) result(im)
            !! Return true if harmonic with given index is imaginary

            class(Grid)  ,intent(in) :: this
            integer(ik)  ,intent(in) :: ind
            logical                  :: im
            
        end function isImaginary
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module grid_class
!-----------------------------------------------------------------------------------------------------------------------------------
 