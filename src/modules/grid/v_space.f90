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
module v_space_class
    !! author: Stefan Mijin 
    !! 
    !!  Houses object responsible for moment-taking and velocity space interpolation

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use grid_class                  ,only: Grid
    use physics_functions           
    use support_functions           ,only: findNearestPointsInArray
    use sparse_row_data_class       ,only: SparseRowData
    use physical_constants

    implicit none
    private

    type ,public ,extends(Object) :: VSpace
        !! Object responsible for taking velocity moments of distributions and other velocity grid manipulation

        real(rk) ,allocatable ,dimension(:) ,private :: vGrid    !! Copy of velocity grid
        real(rk) ,allocatable ,dimension(:) ,private :: vWidths  !! Velocity grid cell widths

        real(rk) ,allocatable ,dimension(:) ,private :: linInterp !! Linear interpolation coefficient for cell faces in velocity space 

        integer(ik) :: numV !! Number of velocity cells 
        integer(ik) :: numH !! Number of harmonics

        contains

        procedure ,public :: getVCellWidths
        procedure ,public :: getVGrid 
        procedure ,public :: getNumH
        procedure ,public :: getNumV

        procedure ,public :: getVLinInterp

        procedure ,public :: calculateMoment

        procedure ,public :: getNearestPoints
        procedure ,public :: getContainingVCell

        procedure ,public :: getShkarofskyIMat
        procedure ,public :: getShkarofskyJMat

        procedure ,public :: init => initVSpace

    end type VSpace
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initVSpace(this,gridObj) 
        !! VSpace initialization routine

        class(VSpace)           ,intent(inout)  :: this
        type(Grid)              ,intent(in)     :: gridObj !! Grid object used to initialize VSpace

    end subroutine initVSpace
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNumH (this) result(numH)
        !! Return total number of resolved harmonic on grid

        class(VSpace) ,intent(in) :: this
        integer(ik)               :: numH

    end function getNumH
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNumV (this) result(numV)
        !! Return number of v points on grid

        class(VSpace) ,intent(in) :: this
        integer(ik)             :: numV

    end function getNumV
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getVGrid (this) result(v)
       !! Getter for vGrid

        class(VSpace)                         ,intent(in) :: this
        real(rk)   ,allocatable ,dimension(:)             :: v

    end function getVGrid
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getVCellWidths (this) result(dv)
        !! Getter for vWidths

        class(VSpace)                       ,intent(in) :: this
        real(rk)   ,allocatable ,dimension(:)           :: dv

    end function getVCellWidths
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getVLinInterp (this) result(linInterp)
        !! Getter for linInterp

        class(VSpace)                       ,intent(in) :: this
        real(rk)   ,allocatable ,dimension(:)           :: linInterp

    end function getVLinInterp
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNearestPoints (this,v) result(pair)
        !! Return two nearest points for given velocity value v. If the first point is 0, the velocity is below the lowest velocity in the 
        !! grid, and if the second point is 0 the velocity is above the greatest v in the grid.

        class(VSpace)             ,intent(in) :: this
        real(rk)                  ,intent(in) :: v
        integer(ik)   ,dimension(2)           :: pair

    end function getNearestPoints
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getContainingVCell (this,v) result(ind)
        !! Return index of cell which containes the given velocity value v. If the returned index is 0, the point is outside of the grid. 

        class(VSpace)             ,intent(in) :: this
        real(rk)                  ,intent(in) :: v
        integer(ik)                           :: ind

    end function getContainingVCell
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function calculateMoment (this,f,h,mOrder,g,gDependsOnX) result(res)
        !! Calculate moment of the h-th harmonic of passed local distribution function f, optionally multiplied by g(v). If gDependsOnX 
        !! is .true. assumes that g is given as a strided array with stride numV, and it is used to allocate the result - 
        !! this allows g to not have a halo while f does. The moment is of order mOrder. 

        class(VSpace)                       ,intent(in) :: this
        real(rk)    ,dimension(:)           ,intent(in) :: f
        integer(ik)                         ,intent(in) :: h, mOrder
        real(rk)    ,optional ,dimension(:) ,intent(in) :: g
        logical     ,optional               ,intent(in) :: gDependsOnX

        real(rk)   ,allocatable ,dimension(:)           :: res

    end function calculateMoment
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getShkarofskyIMat (this,index) result(res)
        !! Return sparse row data format for the lower triangular Shkarofsky I_index integral.

        class(VSpace)                       ,intent(in) :: this
        integer(ik)                         ,intent(in) :: index
        type(SparseRowData)                             :: res

    end function getShkarofskyIMat
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getShkarofskyJMat (this,index) result(res)
        !! Return sparse row data format for the upper triangular Shkarofsky J_index integral.

        class(VSpace)                       ,intent(in) :: this
        integer(ik)                         ,intent(in) :: index
        type(SparseRowData)                             :: res

    end function getShkarofskyJMat
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module v_space_class
!-----------------------------------------------------------------------------------------------------------------------------------
 