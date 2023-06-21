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
submodule (fluid_gen1d_class) fluid_gen1d_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the FluidStencilGen1D class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initGenerator(this,stencilObj,columnVecs,varContColVarNames,mbColVarNames,periodicDim,coordInterval) 
    !! 1D fluid variable stencil value generator initialization routine

    class(FluidStencilGen1D)              ,intent(inout)  :: this
    type(Stencil1D)                       ,intent(in)     :: stencilObj 
    type(RealArray) ,dimension(:)         ,intent(in)     :: columnVecs 
    type(StringArray) ,dimension(:)       ,intent(in)     :: varContColVarNames
    type(StringArray) ,dimension(:)       ,intent(in)     :: mbColVarNames
    logical ,optional                     ,intent(in)     :: periodicDim 
    integer(ik) ,dimension(2) ,optional   ,intent(in)     :: coordInterval

    call this%fixedStencilGen%init(stencilObj,columnVecs,periodicDim,coordInterval)

    if (assertions) then
        call assert(size(varContColVarNames)==stencilObj%getStencilDims(),&
        "varContColVarNames passed to FluidStencilGen1D constructor must conform to stencil size")
        call assert(size(mbColVarNames)==stencilObj%getStencilDims(),&
        "mbColVarNames passed to FluidStencilGen1D constructor must conform to stencil size")
    end if

    this%mbColVarNames = mbColVarNames
    this%varContColVarNames = varContColVarNames

    call this%makeDefined()

end subroutine initGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcVals(this,varCont,res,mbData,hostModel)
    !! Calculate variable fluid 1D stencil values in place (does not depend on hostModel)

    class(FluidStencilGen1D)                    ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    integer(ik) :: i,j ,localSize ,mbDataDim ,mbLBound ,mbHalo

    real(rk)        ,allocatable ,dimension(:)      :: modelboundDataVals

    if (.not. allocated(this%fixedStencilVals)) &
    call this%fixedStencilGen%calculateInPlace(varCont,this%fixedStencilVals,mbData,hostModel)

    if (.not. allocated(this%varContVarIndices)) then 
        allocate(this%varContVarIndices(size(this%varContColVarNames)))
        this%varContVarIndices = 0
        do i = 1, size(this%varContVarIndices)
            if (allocated(this%varContColVarNames(i)%string)) &
            this%varContVarIndices(i) = varCont%getVarIndex(this%varContColVarNames(i)%string)
        end do
    end if

    localSize = this%fixedStencilGen%coordInterval(2) - this%fixedStencilGen%coordInterval(1) + 1

    if (.not. allocated(this%columnValueBuffers)) allocate(this%columnValueBuffers(size(this%varContColVarNames)))

    ! Fill out variable buffers
    do i = 1, size(this%varContColVarNames)
        if (.not. allocated(this%columnValueBuffers(i)%entry)) allocate(this%columnValueBuffers(i)%entry(localSize))
        if (this%varContVarIndices(i)>0) then
            this%columnValueBuffers(i)%entry = varCont%variables(this%varContVarIndices(i))%entry(1:localSize)
        else
            this%columnValueBuffers(i)%entry = real(1,kind=rk)
        end if

        if (allocated(this%mbColVarNames(i)%string)) then
            if (assertions) call assert(present(mbData),&
            "No modelbound data present in fluid gen 1d calc call when mb col vars present")

            ! Gfortran select rank workaround
            if (.not. allocated(modelboundDataVals)) allocate(modelboundDataVals(0))

            call mbData%copyData(this%mbColVarNames(i)%string,modelboundDataVals)
            mbDataDim = mbData%getDataDim(this%mbColVarNames(i)%string)

            if (assertions) call assert (mbDataDim==1,&
            "Only fluid variables are support as modelbound data column variables in FluidGen1D")

            mbLBound = lbound(modelboundDataVals,1)
            mbHalo = (size(modelboundDataVals) - size(this%columnValueBuffers(i)%entry))/2

            this%columnValueBuffers(i)%entry = this%columnValueBuffers(i)%entry &
                                               * modelboundDataVals(mbLBound+mbHalo:mbLBound+localSize+mbHalo-1)

        end if
    end do

    res = this%fixedStencilVals

    do i = 1, size(res)
        do j = 1,size(res(i)%entry)
            res(i)%entry(j) = this%columnValueBuffers(this%fixedStencilGen%presentColumns(i)%entry(j))%entry(i) *  res(i)%entry(j)
        end do
    end do  

end subroutine calcVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule fluid_gen1d_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
