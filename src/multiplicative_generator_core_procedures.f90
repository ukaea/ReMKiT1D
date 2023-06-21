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
submodule (multiplicative_generator_core_class) multiplicative_generator_core_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the multiplicative generator core class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCore(this,gridObj,partitionObj,procRank,rowCoords,fluidCol,fixedVVals,vValsDependOnX) 
    !! Multiplicative generator core initialization routine

    class(MultiplicativeGeneratorCore)        ,intent(inout)  :: this
    type(Grid)                                ,intent(in)     :: gridObj !! Grid object used to determine local number of rows
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    integer(ik) ,dimension(:,:)               ,intent(in)     :: rowCoords !! Global row coordinate values
    logical ,optional                         ,intent(in)     :: fluidCol !! True if column variable is fluid. Defaults to false.
    real(rk) ,optional ,dimension(:)          ,intent(in)     :: fixedVVals !! Fixed velocity stencil values when the row variable is fluid and the column variable is kinetic
    logical ,optional                         ,intent(in)     :: vValsDependOnX !! True if v stencil has spatial dependence. Defaults to false.

    if (assertions) then 
        call assert(gridObj%isDefined(),"Undefined grid object passed to multiplicative generator core constructor")
        call assert(partitionObj%isDefined(),"Undefined partition object passed to multiplicative generator core constructor")
    end if

    this%localRowCoords = partitionObj%filterCoords(procRank+1,rowCoords,normalize=.true.)
    this%expXSize = size(removeDupeInts(this%localRowCoords(1,:)))
    if (size(rowCoords,1) > 1) this%expVSize = size(removeDupeInts(this%localRowCoords(3,:)))

    this%fluidCol = .false. 
    if (present(fluidCol)) this%fluidCol = fluidCol

    if (.not. this%fluidCol) then 

        if (present(fixedVVals) .and. size(rowCoords,1) == 1) this%fixedVVals = fixedVVals 
    end if
    this%vValsDependOnX = .false.
    if (present(vValsDependOnX)) this%vValsDependOnX = vValsDependOnX

    call this%makeDefined() 

end subroutine initCore
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcCoreVals(this,res,xVals,hVals,vVals)
    !! Calculate multiplicative core values. All passed stencil values should be indexed starting at the first evolved point. If there are
    !! any gaps in the row values those must also be included for indexing to work.

    class(MultiplicativeGeneratorCore)          ,intent(inout) :: this
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    type(RealArray)           ,dimension(:)     ,intent(in)    :: xVals !! x stencil values 
    type(RealArray) ,optional ,dimension(:)     ,intent(in)    :: hVals !! h stencil values
    type(RealArray) ,optional ,dimension(:)     ,intent(in)    :: vVals !! v stencil values

    integer(ik) :: i 

    if (assertions) then 
        call assert(this%isDefined(),"calcCoreVals called on undefined generator core")
        if (.not. this%fluidCol) then 
            call assert(present(hVals),"hVals must be present in calcCoreVals if the column variable is kinetic")
            call assert(present(vVals),"vVals must be present in calcCoreVals if the column variable is kinetic")
        end if
    end if

    if (allocated(res)) then
        call assert(size(res) == size(this%localRowCoords,2),&
        "res passed to calcCoreVals does not conform to expected number of rows")
    else
        allocate(res(size(this%localRowCoords,2)))
    end if

    if (this%fluidCol) then 
        do i = 1,size(res)
            res(i)%entry = xVals(this%localRowCoords(1,i))%entry
        end do
    else
        if (size(this%localRowCoords,1) == 1) then
            if (assertions) call assert(allocated(this%fixedVVals),&
            "calcCoreVals called for fluid row and kinetic column variable with no fixedVVals allocated in the core")
            !Note: this assumes that there is only one harmonic in the harmonic stencil
            do i = 1,size(res)
                res(i)%entry = flatTensorProduct(xVals(this%localRowCoords(1,i))%entry,&
                                                            [real(1,kind=rk)],&
                                                            this%fixedVVals)
            end do
        else

            if (this%vValsDependOnX) then 
                do i = 1,size(res) 
                    res(i)%entry = flatTensorProduct(xVals(this%localRowCoords(1,i))%entry,&
                                                    hVals(this%localRowCoords(2,i))%entry,&
                                               vVals((this%localRowCoords(1,i)-1)*this%expVSize + this%localRowCoords(3,i))%entry)
                end do
            else
                do i = 1,size(res) 
                    res(i)%entry = flatTensorProduct(xVals(this%localRowCoords(1,i))%entry,&
                                                    hVals(this%localRowCoords(2,i))%entry,&
                                                    vVals(this%localRowCoords(3,i))%entry)
                end do
            end if
        end if
    end if

end subroutine calcCoreVals
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule multiplicative_generator_core_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
