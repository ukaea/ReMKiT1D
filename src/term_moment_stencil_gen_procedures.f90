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
submodule (term_moment_stencil_gen_class) term_moment_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the term moment stencil generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initTermMomentGen(this,partitionObj,vspaceObj,procRank,momentOrder,termName,removeLastCell) 
    !! Term moment stencil value generator initialization routine

    class(TermMomentStencilGenerator)         ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get interpolation object
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    integer(ik)                               ,intent(in)     :: momentOrder !! Order of moment to be taken
    character(*)                              ,intent(in)     :: termName !! Name of term in host model whose moment should be taken
    logical ,optional                         ,intent(in)     :: removeLastCell !! Set to true if the row variable is staggered and the grid is not periodic. Defaults to false.

    integer(ik) :: minX ,maxX

    logical :: rLastCell

    if (assertions) then 
        call assert(partitionObj%isDefined(),"Undefined partition passed to term moment stencil value generator constructor")
        call assert(vspaceObj%isDefined(),"Undefined VSpace passed to term moment stencil value generator constructor")
    end if

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    this%locNumX = maxX - minX + 1

    this%termName = termName 

    rLastCell = .false. 
    if (present(removeLastCell)) rLastCell = removeLastCell

    if (rLastCell .and. maxX == maxval(partitionObj%getMaxX())) this%locNumX = this%locNumX - 1

    this%vVec = 4*pi*vspaceObj%getVGrid()**(2+momentOrder)*vspaceObj%getVCellWidths()

    this%numV = vspaceObj%getNumV()

    call this%makeDefined()

end subroutine initTermMomentGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcTermMomentVals(this,varCont,res,mbData,hostModel)
    !! Calculate term moment stencil values in place. Requires hostModel.

    class(TermMomentStencilGenerator)           ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    integer(ik) :: i ,dim

    if (assertions) then 
        call assert(this%isDefined(),"calcTermMomentVals called from undefined stencil generator")
        call assert(present(hostModel),"No hostModel passed to calcTermMomentVals when expected")
    end if

    select type(hostModel)
    type is (Model)
        if (.not. allocated(this%indexingDataBuffer)) then 
            allocate(this%indexingDataBuffer)
            this%indexingDataBuffer = hostModel%getImplicitTermIndexingData(hostModel%getImplicitTermIndex(this%termName))
        end if

        this%matBuffer = hostModel%getImplicitTermRowData(hostModel%getImplicitTermIndex(this%termName))
    class default 
        error stop "Unexpected ModelSurrogate detected in calcTermMomentVals"
    end select

    dim = size(this%indexingDataBuffer%rowDataCoordsLocal(1)%colCoords,1)

    if (allocated(res)) then
        if (assertions) call assert(size(res) == this%locNumX,"res passed to calcTermMomentVals has unexpected size")
    else
        allocate(res(this%locNumX))
        do i = 1,size(this%indexingDataBuffer%rowDataCoordsLocal)
            if (.not. allocated(res(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(1))%entry)) &
            allocate(res(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(1))%entry(&
                     size(this%indexingDataBuffer%rowDataCoordsLocal(i)%colCoords(dim,:))))
        end do
    end if
    do i = 1,size(res)
        res(i)%entry = 0
    end do
    if (dim == 3) then
        do i = 1,size(this%indexingDataBuffer%rowDataCoordsLocal)
            res(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(1))%entry(&
                this%indexingDataBuffer%rowDataCoordsLocal(i)%colCoords(dim,:)) = &
                res(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(1))%entry(&
                this%indexingDataBuffer%rowDataCoordsLocal(i)%colCoords(dim,:)) &
                + this%vVec(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(3)) * this%matBuffer%values(i)%entry
        end do
    else
        do i = 1,size(this%indexingDataBuffer%rowDataCoordsLocal)
            res(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(1))%entry = &
            res(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(1))%entry &
                + this%vVec(this%indexingDataBuffer%rowDataCoordsLocal(i)%rowCoords(3)) * this%matBuffer%values(i)%entry
        end do
    end if
end subroutine calcTermMomentVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule term_moment_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
