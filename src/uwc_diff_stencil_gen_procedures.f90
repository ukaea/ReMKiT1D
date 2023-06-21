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
submodule (uwc_diff_stencil_gen_class) uwc_diff_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the central/upwinding difference stencil value generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initUWCDiffStencil(this,partitionObj,procRank,innerJ,outerJ,linInterp,xPeriodic,interpVarIndex,upwindingMode,&
    staggeredGridMode) 
!! Central/upwind differentiation stencil value generator initialization routine

class(UWCDiffStencilValGenerator)         ,intent(inout)  :: this
type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
real(rk) ,dimension(:)                    ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star (should conform to x-grid)
real(rk) ,dimension(:)                    ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star (should conform to x-grid)
real(rk) ,dimension(:)                    ,intent(in)     :: linInterp !! Linear interpolation coefficients to right cell faces (should be size(xGrid)+1)
logical ,optional                         ,intent(in)     :: xPeriodic !! Used to determine if outer processors should have their stencils trimmed. Defaults to .false.  
integer(ik) ,optional                     ,intent(in)     :: interpVarIndex !! Optional interpolated variable index 
integer(ik) ,optional                     ,intent(in)     :: upwindingMode !! 0 - no upwinding, 1 - upwinding with stencil containing the flux jacobian, 2- upwinding with stencil not including flux jacobian. Defaults to 0.
logical ,optional                         ,intent(in)     :: staggeredGridMode !! Removes last row if true

    integer(ik) :: inferredGridSize

    if (assertions) call assert(partitionObj%isDefined(),&
    "Undefined partition object passed to central diff val generator constructor")

    inferredGridSize = maxval(partitionObj%getMaxX())

    if (assertions) then 
        call assert(size(innerJ) == inferredGridSize,"innerJ passed to initUWCDiffStencil does not conform to inferred grid size")
        call assert(size(outerJ) == inferredGridSize,"outerJ passed to initUWCDiffStencil does not conform to inferred grid size")
        call assert(size(linInterp) == inferredGridSize+1,&
        "linInterp passed to initUWCDiffStencil does not conform to inferred grid size")

        if (present(upwindingMode)) then 
            call assert (any(upwindingMode == [0,1,2]),"Incompatible upwinding mode passed to initUWCDiffStencil")
            if (any(upwindingMode == [1,2])) call assert(present(interpVarIndex),"If upwindingMode is 1 or 2 in initUWCDiffStencil&
            & the interpolated variable index must be passed as well")
        end if
    end if

    this%periodicGrid = .false. 

    if (present(xPeriodic)) this%periodicGrid = xPeriodic

    this%minX = partitionObj%getMinXAtInd(procRank+1)
    this%maxX = partitionObj%getMaxXAtInd(procRank+1)

    this%innerJ = innerJ 
    this%outerJ = outerJ
    allocate(this%linInterp(0:this%maxX-this%minX+1))
    this%linInterp = linInterp(this%minX:this%maxX+1)

    allocate(this%linInterpEffective(0:this%maxX-this%minX+1))
    this%linInterpEffective = linInterp(this%minX:this%maxX+1)

    if (present(interpVarIndex)) allocate(this%interpVarIndex,source=interpVarIndex)

    allocate(this%interpVarBuffer(0:this%maxX-this%minX+1))
        this%interpVarBuffer = real(1,kind=rk)

    this%upwindingMode = 0

    if (present(upwindingMode)) this%upwindingMode = upwindingMode

    this%staggeredGridMode = .false. 
    if (present(staggeredGridMode)) this%staggeredGridMode = staggeredGridMode

    call this%makeDefined()

end subroutine initUWCDiffStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcUWCDiffVals(this,varCont,res,mbData,hostModel)
    !! Calculate central/upwind diff stencil values in place 

    class(UWCDiffStencilValGenerator)           ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    logical :: containsLeftBoundary ,containsRightBoundary

    integer(ik) :: i ,localSize

    if (assertions) then 
        call assert(this%isDefined(),"calcUWCDiffVals called from undefined stencil value generator")
        call assert(varCont%isDefined(),"Undefined variable container passed to calcUWCDiffVals")
    end if

    containsLeftBoundary = this%minX == 1
    containsRightBoundary = this%maxX == size(this%innerJ)

    localSize = this%maxX - this%minX + 1
    if (containsRightBoundary .and. (.not. this%periodicGrid) .and. this%staggeredGridMode) localSize = localSize - 1
    if (allocated(res)) then 
        if (assertions) &
        call assert(size(res) == localSize,"res passed to calcUWCDiffVals does not conform to local x-grid size")
    else
        allocate(res(localSize))
    end if
    if (allocated(this%interpVarIndex)) then 
        this%interpVarBuffer(0:size(res)) = &
        varCont%variables(this%interpVarIndex)%entry(0:size(res))*(1-this%linInterp(0:size(res))) &
        + varCont%variables(this%interpVarIndex)%entry(1:size(res)+1)*this%linInterp(0:size(res))
        if (this%upwindingMode>0) then 
            this%linInterpEffective = real(1,kind=rk)

            do i = 0,size(res)
                if (this%interpVarBuffer(i) > 0) this%linInterpEffective(i) = 0
            end do

            if (this%upwindingMode == 2) this%interpVarBuffer = real(1,kind=rk) !TODO: This is dangerous, should be fixed to handle sign change!
        end if
    end if
    if (containsLeftBoundary) then 
        if (this%periodicGrid) then 
            res(1)%entry = &
            this%outerJ(this%minX)*&
            [-this%innerJ(size(this%innerJ))*(real(1,kind=rk)-this%linInterpEffective(0))*this%interpVarBuffer(0),& ! left cell
            this%innerJ(this%minX)*(real(1,kind=rk)-this%linInterpEffective(1))*this%interpVarBuffer(1)& !central cell
            -this%innerJ(size(this%innerJ))*this%linInterpEffective(0)*this%interpVarBuffer(0),&
            this%innerJ(this%minX)*this%linInterpEffective(1)*this%interpVarBuffer(1)] !right cell
        else
            res(1)%entry = &
            this%outerJ(this%minX)*&
            [this%innerJ(this%minX)*(real(1,kind=rk)-this%linInterpEffective(1))*this%interpVarBuffer(1),& !central cell
            this%innerJ(this%minX)*this%linInterpEffective(1)*this%interpVarBuffer(1)] !right cell
        end if
    else 
        res(1)%entry = &
            this%outerJ(this%minX)*&
            [-this%innerJ(this%minX-1)*(real(1,kind=rk)-this%linInterpEffective(0))*this%interpVarBuffer(0),& ! left cell
            this%innerJ(this%minX)*(real(1,kind=rk)-this%linInterpEffective(1))*this%interpVarBuffer(1)& !central cell
            -this%innerJ(this%minX-1)*this%linInterpEffective(0)*this%interpVarBuffer(0),&
            this%innerJ(this%minX)*this%linInterpEffective(1)*this%interpVarBuffer(1)] !right cell
    end if

    do i = 1,localSize-2
        res(i+1)%entry = &
            this%outerJ(this%minX + i)*&
            [-this%innerJ(this%minX + i -1)*(real(1,kind=rk)-this%linInterpEffective(i))*this%interpVarBuffer(i),& ! left cell
            this%innerJ(this%minX + i) *(real(1,kind=rk)-this%linInterpEffective(i+1))*this%interpVarBuffer(i+1)& !central cell
            -this%innerJ(this%minX + i -1)*this%linInterpEffective(i)*this%interpVarBuffer(i),&
            this%innerJ(this%minX + i )*this%linInterpEffective(i+1)*this%interpVarBuffer(i+1)] !right cell
    end do

    if (this%periodicGrid .or. .not. containsRightBoundary) then 
        res(size(res))%entry = &
        this%outerJ(this%minX + localSize -1)*&
        [-this%innerJ(this%minX + localSize -2)*(real(1,kind=rk)-this%linInterpEffective(size(res)-1))&
        *this%interpVarBuffer(size(res)-1),& ! left cell
        this%innerJ(this%minX + localSize -1) *(real(1,kind=rk)-this%linInterpEffective(size(res)))*this%interpVarBuffer(size(res))& !central cell
        -this%innerJ(this%minX + localSize -2)*this%linInterpEffective(size(res) -1)*this%interpVarBuffer(size(res)-1),&
        this%innerJ(this%minX + localSize -1)*this%linInterpEffective(size(res))*this%interpVarBuffer(size(res))] !right cell
    else
        res(size(res))%entry = &
        this%outerJ(this%minX + localSize -1)*&
        [-this%innerJ(this%minX + localSize -2)*(real(1,kind=rk)-this%linInterpEffective(size(res)-1))&
        *this%interpVarBuffer(size(res)-1),& ! left cell
        -this%innerJ(this%minX + localSize -2)*this%linInterpEffective(size(res) -1)*this%interpVarBuffer(size(res)-1)] ! central cell
    end if
end subroutine calcUWCDiffVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule uwc_diff_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
