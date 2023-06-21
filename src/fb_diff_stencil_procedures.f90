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
submodule (fb_diff_stencil_gen_class) fb_diff_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the forward/backward difference stencil value generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initFBValGen(this,partitionObj,procRank,innerJ,outerJ,xPeriodic,backwardsDiff,staggeredGridMode) 
    !! Forward/backward differentiation stencil value generator initialization routine

    class(FBDiffStencilValGenerator)          ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    real(rk) ,dimension(:)                    ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star (should conform to x-grid)
    real(rk) ,dimension(:)                    ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star (should conform to x-grid)
    logical ,optional                         ,intent(in)     :: xPeriodic !! Used to determine if outer processors should have their stencils trimmed. 
                                                                           !! Defaults to .false.  
    logical ,optional                         ,intent(in)     :: backwardsDiff !! If true the difference is backwards. Defaults to .false.
    logical ,optional                         ,intent(in)     :: staggeredGridMode !! If true will zero out the rightmost contribution to backwards differencing or remove last row for forward differencing

    integer(ik) :: inferredGridSize

    if (assertions) call assert(partitionObj%isDefined(),"Undefined partition object passed to fb diff val generator constructor")

    inferredGridSize = maxval(partitionObj%getMaxX())

    if (assertions) then 
        call assert(size(innerJ) == inferredGridSize,"innerJ passed to initFBValGen does not conform to inferred grid size")
        call assert(size(outerJ) == inferredGridSize,"outerJ passed to initFBValGen does not conform to inferred grid size")
    end if

    this%periodicGrid = .false. 
    this%backwardsDiff = .false. 

    if (present(xPeriodic)) this%periodicGrid = xPeriodic
    if (present(backwardsDiff)) this%backwardsDiff = backwardsDiff 

    this%minX = partitionObj%getMinXAtInd(procRank+1)
    this%maxX = partitionObj%getMaxXAtInd(procRank+1)

    this%innerJ = innerJ 
    this%outerJ = outerJ

    this%staggeredGridMode = .false. 
    if (present(staggeredGridMode)) this%staggeredGridMode = staggeredGridMode
    
    call this%makeDefined()

end subroutine initFBValGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcFBVals(this,varCont,res,mbData,hostModel)
    !! Calculate forward/backwards diff stencil values in place (does not depend on varCont)

    class(FBDiffStencilValGenerator)            ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    logical :: containsLeftBoundary ,containsRightBoundary

    integer(ik) :: i ,localSize

    if (assertions) then 
        call assert(this%isDefined(),"calcFBVals called from undefined stencil value generator")
        call assert(varCont%isDefined(),"Undefined variable container passed to calcFBVals")
    end if

    containsLeftBoundary = this%minX == 1
    containsRightBoundary = this%maxX == size(this%innerJ)

    localSize = this%maxX - this%minX + 1
    if (.not. this%periodicGrid .and. this%staggeredGridMode) then 
        if (.not. this%backwardsDiff .and. containsRightBoundary) localSize = localSize - 1
    end if

    if (allocated(res)) then 
        if (assertions) &
        call assert(size(res) == localSize,"res passed to calcFBVals does not conform to local x-grid size")
    else
        allocate(res(localSize))
    end if
        
    if (this%backwardsDiff) then 
        if (containsLeftBoundary) then 
            if (this%periodicGrid) then 
                res(1)%entry = this%outerJ(this%minX)*[-this%innerJ(size(this%innerJ)),this%innerJ(this%minX)]
            else
                res(1)%entry = this%outerJ(this%minX)*[this%innerJ(this%minX)]
            end if
        else 
            res(1)%entry = this%outerJ(this%minX)*[-this%innerJ(this%minX-1),this%innerJ(this%minX)]
        end if

        do i = 1,this%maxX - this%minX - 1
            res(i+1)%entry = this%outerJ(this%minX+i)*[-this%innerJ(this%minX-1+i),this%innerJ(this%minX+i)]
        end do

        if (this%staggeredGridMode .and. (.not. this%periodicGrid ).and. containsRightBoundary) then 
            res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX-1),real(0,kind=rk)]
        else
            res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX-1),this%innerJ(this%maxX)]
        end if

    else

        do i = 1,this%maxX - this%minX
            res(i)%entry = this%outerJ(this%minX+i-1)*[-this%innerJ(this%minX-1+i),this%innerJ(this%minX+i)]
            
        end do

        if (containsRightBoundary) then 
            if (this%periodicGrid) then 
                res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX),this%innerJ(1)]
            else if (.not. this%staggeredGridMode) then
                res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX)]
            end if
        else 
            res(size(res))%entry = this%outerJ(this%maxX)*[-this%innerJ(this%maxX),this%innerJ(this%maxX+1)]
        end if
    end if

end subroutine calcFBVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule fb_diff_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
