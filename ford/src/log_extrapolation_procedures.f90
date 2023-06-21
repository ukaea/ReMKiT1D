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
submodule (log_extrapolation_class) log_extrapolation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the LogExtrapolation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
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

    real(rk) ,allocatable ,dimension(:) :: dx,linInterp
    real(rk)                            :: dx1,dx2 ,a,b
    integer(ik) :: minX ,maxX

    if (assertions) then 
        call assert(partObj%isDefined(),"Undefined partition object passed to initLogExtrap")
        call assert(gridObj%isDefined(),"Undefined grid object passed to initLogExtrap")
        call assert(geometryObj%isDefined(),"Undefined geometry object passed to initLogExtrap")
    end if

    call this%setStaggeredVars(staggeredVars)
    call this%setLeftBoundary(leftBoundary)
    call this%setOnBoundary(partObj,gridObj,numProc)
    call this%setHaloWidth(haloWidth)

    minX = partObj%getMinXAtInd(numProc+1)
    maxX = partObj%getMaxXAtInd(numProc+1)

    if (this%hasBoundary()) then

        this%interpolate = interpolate

        dx = geometryObj%getCellWidths()
        allocate(linInterp,source=geometryObj%getLinInterp(dualGrid=staggeredVars))

        this%linInterp = real(1,kind=rk) - linInterp(gridObj%getNumX()-1)
        dx1 = dx(gridObj%getNumX())
        dx2 = dx(gridObj%getNumX()-1)
        if (leftBoundary) then 
            this%linInterp = linInterp(1)
            dx1 = dx(1)
            dx2 = dx(2)
        end if

        this%exterpCoords = [maxX-minX+1,maxX-minX]
        if (staggeredVars) this%exterpCoords = this%exterpCoords - 1

        if (leftBoundary) this%exterpCoords = [1,2]

        b = dx1/2
        if (staggeredVars) b = dx1

        if (staggeredVars) then 
            a = dx2 
            if (interpolate) a = dx2/2
        else
            a = dx1/2 + dx2/2
            if (interpolate) a = dx1/2
        end if

        this%logExterp = b/a + real(1,kind=rk)
    end if

    call this%makeDefined()

end subroutine initLogExtrap
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function extrapolateLog(this,input) result(res)
    !! Logarithmic extrapolation

    class(LogExtrapolation)  ,intent(in)    :: this 
    real(rk) ,dimension(:)   ,intent(in)    :: input
    real(rk)                                :: res

    integer(ik) :: lBoundData, offset

    real(rk) :: val1 ,val2

    res = 0

    if (this%hasBoundary()) then 

        lBoundData = lbound(input,1)
        offset = lBoundData - 1 + this%getHaloWidth()

        val1 = input(this%exterpCoords(1)+offset)
        val2 = input(this%exterpCoords(2)+offset)

        if (this%interpolate) val2 = (real(1,kind=rk) - this%linInterp)*val1 + this%linInterp*val2

        res = val1**this%logExterp * val2**(real(1,kind=rk)-this%logExterp)

    end if

end function extrapolateLog
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule log_extrapolation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
