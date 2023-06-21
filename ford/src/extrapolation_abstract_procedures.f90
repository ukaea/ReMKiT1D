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
submodule (extrapolation_abstract_class) extrapolation_abstract_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the abstract extrapolation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setOnBoundary(this,partObj,gridObj,numProc) 
    !! Check if this processor has the boundary corresponding with this extrapolation and set active

    class(Extrapolation) ,intent(inout) :: this
    type(Partition)      ,intent(in)    :: partObj
    type(Grid)           ,intent(in)    :: gridObj
    integer(ik)          ,intent(in)    :: numProc

    if (assertions) then 
        call assert(partObj%isDefined(),"Undefined partition object passed to setOnBoundary")
        call assert(gridObj%isDefined(),"Undefined grid object passed to setOnBoundary")
    end if

    if (this%isLeftBoundary()) then 
        this%active = partObj%getMinXAtInd(numProc+1) == 1
    else
        this%active = partObj%getMaxXAtInd(numProc+1) == gridObj%getNumX()
    end if

end subroutine setOnBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setLeftBoundary(this,leftBoundary) 
    !! Setter for leftBoundary

    class(Extrapolation)  ,intent(inout) :: this
    logical               ,intent(in)    :: leftBoundary

    this%leftBoundary = leftBoundary 

end subroutine setLeftBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setStaggeredVars(this,stagVars) 
    !! Setter for staggeredVars

    class(Extrapolation)  ,intent(inout) :: this
    logical               ,intent(in)    :: stagVars

    this%staggeredVars = stagVars

end subroutine setStaggeredVars
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setHaloWidth(this,haloWidth) 
    !! Setter for haloWidth

    class(Extrapolation)  ,intent(inout) :: this
    integer(ik)           ,intent(in)    :: haloWidth

    this%haloWidth = haloWidth

end subroutine setHaloWidth
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function hasBoundary(this) result(active)
    !! Getter for active

    class(Extrapolation)   ,intent(in)  :: this
    logical                             :: active

    active = this%active

end function hasBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isLeftBoundary(this) result(leftBoundary)
    !! Getter for leftBoundary

    class(Extrapolation)   ,intent(in)  :: this
    logical                             :: leftBoundary

    leftBoundary = this%leftBoundary

end function isLeftBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function usesStaggeredVars(this) result(staggeredVars)
    !! Getter for staggeredVars

    class(Extrapolation)   ,intent(in)  :: this
    logical                             :: staggeredVars

    staggeredVars = this%staggeredVars

end function usesStaggeredVars
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getHaloWidth(this) result(haloWidth)
    !! Getter for haloWidth

    class(Extrapolation)   ,intent(in)  :: this
    integer(ik)                         :: haloWidth

    haloWidth = this%haloWidth

end function getHaloWidth
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule extrapolation_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
