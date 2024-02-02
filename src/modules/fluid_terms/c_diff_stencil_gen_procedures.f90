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
submodule (c_diff_stencil_gen_class) c_diff_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the central difference stencil value generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCValGen(this,partitionObj,procRank,innerJ,outerJ,xPeriodic,staggeredGridMode) 
    !! Central differentiation stencil value generator initialization routine

    class(CDiffStencilValGenerator)           ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    real(rk) ,dimension(:)                    ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star (should conform to x-grid)
    real(rk) ,dimension(:)                    ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star (should conform to x-grid)
    logical ,optional                         ,intent(in)     :: xPeriodic !! Used to determine if outer processors should have their stencils trimmed.  Defaults to .false.
    logical ,optional                         ,intent(in)     :: staggeredGridMode !! If true will remove last row/column if not periodic

    integer(ik) :: inferredGridSize

    if (assertions .or. assertionLvl >= 0) call assert(partitionObj%isDefined(),&
    "Undefined partition object passed to central diff val generator constructor")

    inferredGridSize = maxval(partitionObj%getMaxX())

    if (assertions .or. assertionLvl >= 0) then 
        call assert(size(innerJ) == inferredGridSize,"innerJ passed to initCValGen does not conform to inferred grid size")
        call assert(size(outerJ) == inferredGridSize,"outerJ passed to initCValGen does not conform to inferred grid size")
    end if

    this%periodicGrid = .false. 

    if (present(xPeriodic)) this%periodicGrid = xPeriodic

    this%minX = partitionObj%getMinXAtInd(procRank+1)
    this%maxX = partitionObj%getMaxXAtInd(procRank+1)

    this%innerJ = innerJ 
    this%outerJ = outerJ

    this%staggeredGridMode = .false. 
    if (present(staggeredGridMode)) this%staggeredGridMode = staggeredGridMode

    call this%makeDefined()

end subroutine initCValGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcCVals(this,varCont,res,mbData,hostModel)
    !! Calculate central diff stencil values in place (does not depend on varCont)

    class(CDiffStencilValGenerator)             ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    logical :: containsLeftBoundary ,containsRightBoundary

    integer(ik) :: i ,localSize

    if (assertions) then 
        call assert(this%isDefined(),"calcCVals called from undefined stencil value generator")
        call assert(varCont%isDefined(),"Undefined variable container passed to calcCVals")
    end if

    containsLeftBoundary = this%minX == 1
    containsRightBoundary = this%maxX == size(this%innerJ)

    localSize = this%maxX - this%minX + 1
    if (.not. this%periodicGrid .and. this%staggeredGridMode .and. containsRightBoundary) localSize = localSize - 1

    if (allocated(res)) then 
        if (assertions) &
        call assert(size(res) == localSize,"res passed to calcCVals does not conform to local x-grid size")
    else
        allocate(res(localSize))
    end if
    
    if (containsLeftBoundary) then 
        if (this%periodicGrid) then 
            res(1)%entry = this%outerJ(this%minX)*[-this%innerJ(size(this%innerJ)),this%innerJ(this%minX+1)]
        else
            res(1)%entry = this%outerJ(this%minX)*[this%innerJ(this%minX+1)]
        end if
    else 
        res(1)%entry = this%outerJ(this%minX)*[-this%innerJ(this%minX-1),this%innerJ(this%minX+1)]
    end if

    do i = 1,this%maxX - this%minX - 1
        res(i+1)%entry = this%outerJ(this%minX+i)*[-this%innerJ(this%minX-1+i),this%innerJ(this%minX+i+1)]
    end do

    if (containsRightBoundary) then 
        if (this%periodicGrid) then 
            res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX-1),this%innerJ(1)]
        else if (this%staggeredGridMode) then
            res(size(res))%entry = this%outerJ(this%maxX-1)*[-this%innerJ(this%maxX-2),real(0,kind=rk)] !Ignore rightmost dual grid point
        else
            res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX-1)]
        end if
    else
        res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX-1),this%innerJ(this%maxX+1)]
    end if

end subroutine calcCVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule c_diff_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
