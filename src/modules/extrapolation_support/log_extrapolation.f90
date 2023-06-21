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
module log_extrapolation_class
    !! author: Stefan Mijin
    !!
    !! Houses logarithmic Extrapolation object

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

    type ,public ,extends(Extrapolation) :: LogExtrapolation
        !! Linear extrapolation object

        real(rk) ,private :: logExterp !! Logarithmic extrapolation coefficient used
        real(rk) ,private :: linInterp !! Linear interpolation coefficient used when interpolationg

        integer(ik) ,dimension(2) :: exterpCoords !! Coordinates used for extrapolation. Should be the indices of the data closest and second closest to the boundary, respectively.

        logical  ,private :: interpolate !! Interpolates the second closest value to the boundary using the extrapolation coordinates instead of using the value
                                         !! Associated with the second coordinate
        contains

        procedure ,public :: init => initLogExtrap

        procedure ,public :: extrapolate => extrapolateLog

    end type LogExtrapolation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initLogExtrap(this,partObj,gridObj,numProc,haloWidth,geometryObj,leftBoundary,staggeredVars,interpolate) 
        !! Initialization routine for LogExtrapolation object
    
        class(LogExtrapolation) ,intent(inout) :: this
        type(Partition)         ,intent(in)    :: partObj
        type(Grid)              ,intent(in)    :: gridObj
        integer(ik)             ,intent(in)    :: numProc
        integer(ik)             ,intent(in)    :: haloWidth
        type(Geometry)          ,intent(in)    :: geometryObj
        logical                 ,intent(in)    :: leftBoundary
        logical                 ,intent(in)    :: staggeredVars
        logical                 ,intent(in)    :: interpolate
    
    end subroutine initLogExtrap
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function extrapolateLog(this,input) result(res)
        !! Logarithmic extrapolation

        class(LogExtrapolation)  ,intent(in)    :: this 
        real(rk) ,dimension(:)   ,intent(in)    :: input
        real(rk)                                :: res

    end function extrapolateLog
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module log_extrapolation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 