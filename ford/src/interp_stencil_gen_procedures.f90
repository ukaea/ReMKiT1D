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
submodule (interp_stencil_gen_class) interp_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the interpolation stencil value generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initInterpValGen(this,partitionObj,geometryObj,procRank,staggeredGridMode) 
    !! Central differentiation stencil value generator initialization routine

    class(InterpStencilGenerator)             ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    type(Geometry)                            ,intent(in)     :: geometryObj !! Geometry object used to get interpolation/extrapolation coefficients
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    logical ,optional                         ,intent(in)     :: staggeredGridMode !!  If true will interpolate from staggered to regular grid.

    real(rk) ,allocatable ,dimension(:) :: dx ,lInterp
    integer(ik) :: minX, maxX

    if (assertions) call assert(partitionObj%isDefined(),"Undefined partition object passed to interp val generator constructor")

    if (assertions) call assert(geometryObj%isDefined(),"Undefined geometry object passed to initInteprValGen")

    this%periodicGrid = geometryObj%isPeriodic()
    dx = geometryObj%getCellWidths()
    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)

    this%locNumX = maxX - minX + 1
    this%containsLeftBoundary = minX == 1
    this%containsRightBoundary = maxX == size(dx)

    this%staggeredGridMode = .false. 
    if (present(staggeredGridMode)) this%staggeredGridMode = staggeredGridMode

    allocate(lInterp(0:size(dx)))
    lInterp = geometryObj%getLinInterp(dualGrid=this%staggeredGridMode)

    this%linInterp = lInterp(minX:maxX)
    
    this%linExterpL = dx(1)/(2*dx(2))
    this%linExterpR = dx(size(dx))/(2*dx(size(dx)-1))

    call this%makeDefined()

end subroutine initInterpValGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcInterpVals(this,varCont,res,mbData,hostModel)
    !! Calculate interpolation stencil values in place (does not depend on varCont)

    class(InterpStencilGenerator)               ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    integer(ik) :: i ,localSize

    if (assertions) then 
        call assert(this%isDefined(),"calcInterpVals called from undefined stencil value generator")
        call assert(varCont%isDefined(),"Undefined variable container passed to calcInterpVals")
    end if

    localSize = this%locNumX
    if (.not. this%periodicGrid .and. .not. this%staggeredGridMode .and. this%containsRightBoundary) localSize = localSize - 1

    if (allocated(res)) then 
        if (assertions) &
        call assert(size(res) == localSize,"res passed to calcInterpVals does not conform to local x-grid size")
    else
        allocate(res(localSize))
    end if

    if (this%staggeredGridMode) then 

        if (this%containsLeftBoundary .and. .not. this%periodicGrid) then 
            res(1)%entry =[real(1,kind=rk) + this%linExterpL,-this%linExterpL]
        else
            res(1)%entry =[real(1,kind=rk)-this%linInterp(1),this%linInterp(1)]
        end if

        do i = 2,localSize - 1
            res(i)%entry = [real(1,kind=rk)-this%linInterp(i),this%linInterp(i)]
        end do

        if (this%containsRightBoundary .and. .not. this%periodicGrid) then 
            res(localSize)%entry =[-this%linExterpR,real(1,kind=rk) + this%linExterpR]
        else
            res(localSize)%entry =[real(1,kind=rk)-this%linInterp(localSize),this%linInterp(localSize)]
        end if

    else 
        do i = 1,localSize
            res(i)%entry = [real(1,kind=rk)-this%linInterp(i),this%linInterp(i)]
        end do
    end if

end subroutine calcInterpVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule interp_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
