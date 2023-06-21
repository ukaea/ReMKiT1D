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
submodule (stencil_generator1d_class) stencil_generator1d_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the StencilGenerator1D class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initGenerator(this,stencilObj,columnVecs,periodicDim,coordInterval) 
    !! 1D fixed stencil value generator initialization routine

    class(StencilGenerator1D)             ,intent(inout)  :: this
    type(Stencil1D)                       ,intent(in)     :: stencilObj 
    type(RealArray) ,dimension(:)         ,intent(in)     :: columnVecs 
    logical ,optional                     ,intent(in)     :: periodicDim 
    integer(ik) ,dimension(2) ,optional   ,intent(in)     :: coordInterval

    integer(ik) :: i ,j

    if (assertions) then 
        call assert(stencilObj%isDefined(),"Undefined stencil object passed to initGenerator routine of StencilGenerator1D")
        call assert(size(columnVecs) == stencilObj%getStencilDims(),&
        "columnVecs passed to StencilGenerator1D constructor do not conform to stencil size")

        do i =1,size(columnVecs)
            call assert(size(columnVecs(i)%entry)==size(columnVecs(1)%entry),&
            "All columnVecs passed to StencilGenerator1D constructor must have the same size")
        end do
    end if

    this%columnVectors = columnVecs

    this%periodicDim = .false.
    if (present(periodicDim)) this%periodicDim = periodicDim

    this%coordInterval = [1,size(columnVecs(1)%entry)]
    if (present(coordInterval)) then 
        call assert(all(coordInterval>0 .and. coordInterval<=size(columnVecs(1)%entry)),&
        "coordInterval coordinates passed to StencilGenerator1D constructor out of range")
        call assert(coordInterval(2)>=coordInterval(1),"coordInterval passed to StencilGenerator1D not a valid interval")
        this%coordInterval = coordInterval
    end if

    allocate(this%presentColumns(this%coordInterval(2)-this%coordInterval(1)+1))

    do i = 1,size(this%presentColumns)
        this%presentColumns(i)%entry = pack([(j,j=1,stencilObj%getStencilDims())],&
        mask=stencilObj%getMask(i-1+this%coordInterval(1),size(columnVecs(1)%entry),this%periodicDim))
    end do

    call this%makeDefined()
end subroutine initGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcVals(this,varCont,res,mbData,hostModel)
    !! Calculate fixed 1D stencil values in place (does not depend on varCont,mbData, or hostModel)

    class(StencilGenerator1D)                   ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    integer(ik) :: i,j

    if (assertions) then 
        call assert(this%isDefined(),"calcVals called from undefined StencilGenerator1D")
    end if

    if (allocated(res)) then 
        if (assertions) &
        call assert(size(res) == size(this%presentColumns),"res passed to calcVals on StencilGenerator1D not of expected size")
    else
        allocate(res(size(this%presentColumns)))
    end if

    do i = 1, size(this%presentColumns)
        if (allocated(res(i)%entry)) then
            if (assertions) &
            call assert(size(res) == size(this%presentColumns(i)%entry),&
            "res passed to calcVals on StencilGenerator1D not of expected size")
        else
            allocate(res(i)%entry(size(this%presentColumns(i)%entry)))
        end if
        do j = 1,size(this%presentColumns(i)%entry)
            res(i)%entry(j) = this%columnVectors(this%presentColumns(i)%entry(j))%entry(i-1+this%coordInterval(1))
        end do
    end do

end subroutine calcVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule stencil_generator1d_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
