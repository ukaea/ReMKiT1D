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
module lin_extrapolation_class
    !! author: Stefan Mijin
    !!
    !! Houses linear Extrapolation object

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use partition_class                       ,only: Partition
    use grid_class                            ,only: Grid
    use geometry_class                        ,only: Geometry
    use extrapolation_abstract_class          ,only: Extrapolation

    implicit none
    private

    type ,public ,extends(Extrapolation) :: LinExtrapolation
        !! Linear extrapolation object

        real(rk) ,private :: linExterp !! Linear extrapolation coefficient used

        integer(ik) ,dimension(2) :: exterpCoords !! Coordinates used for extrapolation. Should be the indices of the data closest and second closest to the boundary, respectively.

        contains

        procedure ,public :: init => initLinExtrap

        procedure ,public :: extrapolate => extrapolateLin

    end type LinExtrapolation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initLinExtrap(this,partObj,gridObj,numProc,haloWidth,geometryObj,leftBoundary,staggeredVars) 
        !! Initialization routine for LinExtrapolation object
    
        class(LinExtrapolation) ,intent(inout) :: this
        type(Partition)         ,intent(in)    :: partObj
        type(Grid)              ,intent(in)    :: gridObj
        integer(ik)             ,intent(in)    :: numProc
        integer(ik)             ,intent(in)    :: haloWidth
        type(Geometry)          ,intent(in)    :: geometryObj
        logical                 ,intent(in)    :: leftBoundary
        logical                 ,intent(in)    :: staggeredVars
    
    end subroutine initLinExtrap
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function extrapolateLin(this,input) result(res)
        !! Linear extrapolation

        class(LinExtrapolation)  ,intent(in)    :: this 
        real(rk) ,dimension(:)   ,intent(in)    :: input
        real(rk)                                :: res

    end function extrapolateLin
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module lin_extrapolation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 