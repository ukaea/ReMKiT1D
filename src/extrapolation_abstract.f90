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
module extrapolation_abstract_class
    !! author: Stefan Mijin
    !!
    !! Houses abstract Extrapolation object

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use partition_class                       ,only: Partition
    use grid_class                            ,only: Grid

    implicit none
    private

    type ,public ,extends(Object), abstract :: Extrapolation
        !! Abstract Extrapolation object, responsible for extrapolating fluid quantites onto x grid boundaries

        logical ,private :: leftBoundary = .false.
        logical ,private :: staggeredVars = .false. 
        logical ,private :: active = .false.

        integer(ik) ,private :: haloWidth = 0

        contains

        procedure ,public :: setOnBoundary

        procedure ,public :: setLeftBoundary 
        procedure ,public :: setStaggeredVars 
        procedure ,public :: setHaloWidth

        procedure ,public :: hasBoundary
        procedure ,public :: isLeftBoundary
        procedure ,public :: usesStaggeredVars
        procedure ,public :: getHaloWidth

        procedure(extrap) ,deferred :: extrapolate

    end type Extrapolation
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function extrap(this,input) result(res)
            !! Abstract routine for extrapolation

            import :: Extrapolation ,rk

            class(Extrapolation)     ,intent(in) :: this 
            real(rk) ,dimension(:)   ,intent(in) :: input
            real(rk)                             :: res
 
        end function extrap
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setOnBoundary(this,partObj,gridObj,numProc) 
        !! Check if this processor has the boundary corresponding with this extrapolation and set active

        class(Extrapolation) ,intent(inout) :: this
        type(Partition)      ,intent(in)    :: partObj
        type(Grid)           ,intent(in)    :: gridObj
        integer(ik)          ,intent(in)    :: numProc
    
    end subroutine setOnBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setLeftBoundary(this,leftBoundary) 
        !! Setter for leftBoundary

        class(Extrapolation)  ,intent(inout) :: this
        logical               ,intent(in)    :: leftBoundary
    
    end subroutine setLeftBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setStaggeredVars(this,stagVars) 
        !! Setter for staggeredVars

        class(Extrapolation)  ,intent(inout) :: this
        logical               ,intent(in)    :: stagVars
    
    end subroutine setStaggeredVars
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setHaloWidth(this,haloWidth) 
        !! Setter for haloWidth

        class(Extrapolation)  ,intent(inout) :: this
        integer(ik)           ,intent(in)    :: haloWidth
    
    end subroutine setHaloWidth
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function hasBoundary(this) result(active)
        !! Getter for active

        class(Extrapolation)   ,intent(in)  :: this
        logical                             :: active
    
    end function hasBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isLeftBoundary(this) result(leftBoundary)
        !! Getter for leftBoundary

        class(Extrapolation)   ,intent(in)  :: this
        logical                             :: leftBoundary
    
    end function isLeftBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function usesStaggeredVars(this) result(staggeredVars)
        !! Getter for staggeredVars

        class(Extrapolation)   ,intent(in)  :: this
        logical                             :: staggeredVars
    
    end function usesStaggeredVars
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getHaloWidth(this) result(haloWidth)
        !! Getter for haloWidth

        class(Extrapolation)   ,intent(in)  :: this
        integer(ik)                         :: haloWidth
    
    end function getHaloWidth
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module extrapolation_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 