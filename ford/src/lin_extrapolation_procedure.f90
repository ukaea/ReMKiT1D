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
submodule (lin_extrapolation_class) lin_extrapolation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the LinExtrapolation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
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

    real(rk) ,allocatable ,dimension(:) :: dx,dxReg,linInterp
    real(rk)                            :: lInterp
    integer(ik) :: minX ,maxX

    if (assertions) then 
        call assert(partObj%isDefined(),"Undefined partition object passed to initLinExtrap")
        call assert(gridObj%isDefined(),"Undefined grid object passed to initLinExtrap")
        call assert(geometryObj%isDefined(),"Undefined geometry object passed to initLinExtrap")
    end if

    call this%setStaggeredVars(staggeredVars)
    call this%setLeftBoundary(leftBoundary)
    call this%setOnBoundary(partObj,gridObj,numProc)
    call this%setHaloWidth(haloWidth)

    minX = partObj%getMinXAtInd(numProc+1)
    maxX = partObj%getMaxXAtInd(numProc+1)

    if (this%hasBoundary()) then
        dx = geometryObj%getCellWidths(dualGrid=staggeredVars)
        dxReg = geometryObj%getCellWidths()
        allocate(linInterp,source=geometryObj%getLinInterp())

        lInterp = linInterp(gridObj%getNumX()-1)
        if (leftBoundary) lInterp = linInterp(1)
        this%linExterp = real(1,kind=rk) - lInterp
        if (leftBoundary) this%linExterp = lInterp

        if (staggeredVars) then 
            this%linExterp = dxReg(gridObj%getNumX()-1)/dxReg(gridObj%getNumX()-2)
            if (leftBoundary) this%linExterp = dxReg(1)/dxReg(2)
        end if

        this%exterpCoords = [maxX-minX+1,maxX-minX]
        if (staggeredVars) this%exterpCoords = this%exterpCoords - 1

        if (leftBoundary) this%exterpCoords = [1,2]
    end if

    call this%makeDefined()

end subroutine initLinExtrap
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function extrapolateLin(this,input) result(res)
    !! Linear extrapolation

    class(LinExtrapolation)  ,intent(in)    :: this 
    real(rk) ,dimension(:)   ,intent(in)    :: input
    real(rk)                                :: res

    integer(ik) :: lBoundData, offset

    res = 0

    if (this%hasBoundary()) then 

        lBoundData = lbound(input,1)
        offset = lBoundData - 1 + this%getHaloWidth()

        res = (real(1,kind=rk)+this%linExterp) * input(this%exterpCoords(1)+offset) &
              - this%linExterp * input(this%exterpCoords(2)+offset)

    end if

end function extrapolateLin
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule lin_extrapolation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
