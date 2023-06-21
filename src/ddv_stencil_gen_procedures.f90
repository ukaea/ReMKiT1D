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
submodule (ddv_stencil_gen_class) ddv_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the ddv stencil generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDDVValGen(this,partitionObj,vspaceObj,procRank,activeXCoords,fixedC,fixedInterp,mbC,mbInterp,cfAtZero) 
    !! d/dv stencil value generator initialization routine

    class(DDVStencilGenerator)                ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get interpolation object
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    integer(ik) ,dimension(:)                 ,intent(in)     :: activeXCoords !! List of active global x coords
    real(rk) ,optional   ,dimension(:)        ,intent(in)     :: fixedC !! Fixed C values (size numV). Defaults to ones.
    real(rk) ,optional   ,dimension(:)        ,intent(in)     :: fixedInterp !! Fixed vel interpolation values (size numV). Defaults to values from VSpace.
    character(*) ,optional                    ,intent(in)     :: mbC !! Optional modelbound value for C (single harmonic variable). Overrides fixed values.
    character(*) ,optional                    ,intent(in)     :: mbInterp !! Optional modelbound value for interpolation weight (single harmonic variable). Overrides fixed values.
    real(rk) ,optional ,dimension(2)          ,intent(in)     :: cfAtZero !! Optional extrapolation of C*f at zero in the form A1*f(v1)+A2*f(v2) where all A's are fixed. Defaults to 0.

    integer(ik) :: i ,minX ,maxX

    if (assertions) then 
        call assert(partitionObj%isDefined(),"Undefined partition passed to ddv stencil value generator constructor")
        call assert(vspaceObj%isDefined(),"Undefined VSpace passed to ddv stencil value generator constructor")
    end if

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    this%locNumX = maxX - minX + 1

    this%usedXCoords = pack(activeXCoords,activeXCoords >= minX .and. activeXCoords <= maxX) - minX + 1

    this%numV = vspaceObj%getNumV()

    this%bufferC = [(real(1,kind=rk),i=1,this%numV)]
    this%bufferInterp = vspaceObj%getVLinInterp()

    if (present(fixedC)) then
        if (assertions) call assert(size(fixedC)==this%numV,&
        "fixedC passed to ddv stencil value generator constructor not of size numV")

        this%bufferC = fixedC
    end if

    if (present(fixedInterp)) then
        if (assertions) call assert(size(fixedInterp)==this%numV,&
        "fixedInterp passed to ddv stencil value generator constructor not of size numV")

        this%bufferInterp = fixedInterp
    end if

    if (present(mbC)) then
        deallocate(this%bufferC) 
        this%modelboundC = mbC 
    end if 

    if (present(mbInterp)) then 
        deallocate(this%bufferInterp)
        this%modelboundInterp = mbInterp
    end if

    this%cfAtZero = 0
    if (present(cfAtZero)) this%cfAtZero = cfAtZero

    this%vGridWidths = vspaceObj%getVCellWidths()

    call this%makeDefined()

end subroutine initDDVValGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcDDVVals(this,varCont,res,mbData,hostModel)
    !! Calculate d/dv stencil values in place (does not depend on varCont)

    class(DDVStencilGenerator)                  ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    integer(ik) :: i ,j ,k ,offsetC ,offsetInterp
    if (assertions) then 
        call assert(this%isDefined(),"calcDDVVals called from undefined stencil generator")
        call assert(varCont%isDefined(),"Undefined varCont passed to calcDDVVals")
        if (allocated(this%modelboundC) .or. allocated(this%modelboundInterp)) &
        call assert(present(mbData),"No mbData passed to calcDDVVals when expected")
    end if

    if (allocated(res)) then
        if (assertions) call assert(size(res) == size(this%usedXCoords)*this%numV,"res passed to calcDDVVals has unexpected size")
    else
        allocate(res(size(this%usedXCoords)*this%numV))
    end if

    if (allocated(this%modelboundC)) then 
        call mbData%copyData(this%modelboundC,this%bufferC)
        offsetC = (size(this%bufferC) - this%numV*this%locNumX)/2 + lbound(this%bufferC,1) - 1
        offsetC = offsetC/this%numV 
    end if
    if (allocated(this%modelboundInterp)) then 
        call mbData%copyData(this%modelboundInterp,this%bufferInterp)
        offsetInterp = (size(this%bufferInterp) - this%numV*this%locNumX)/2 + lbound(this%bufferInterp,1) - 1

        offsetInterp = offsetInterp/this%numV 
    end if
    k = 1

    if (allocated(this%modelboundC)) then 
        if (allocated(this%modelboundInterp)) then 
            do i = 1,size(this%usedXCoords)

                    !At first v cell
                    res(k)%entry = [this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+1)*(real(1,kind=rk)&
                                   -this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+1))-this%cfAtZero(1),&
                                    this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+1)&
                                    *this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+1)-this%cfAtZero(2)]

                    res(k)%entry = res(k)%entry/this%vGridWidths(1) 

                    k = k + 1 

                    do j = 2,this%numV-1
                        res(k)%entry = [-this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j-1)&
                        *(real(1,kind=rk)-this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+j-1)),& 
                                        -this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j-1)&
                                        *this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+j-1)&
                                        +this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j)&
                                        *(real(1,kind=rk)-this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+j)),&
                                        this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j)&
                                        *this%bufferInterp((i+offsetInterp-1)*this%numV+j)]

                        res(k)%entry = res(k)%entry/this%vGridWidths(j) 

                        k = k + 1
                    end do
                    ! Last v cell
                    res(k)%entry = [-this%bufferC((this%usedXCoords(i)+offsetC)*this%numV-1)*(real(1,kind=rk)&
                                    -this%bufferInterp((this%usedXCoords(i)+offsetInterp)*this%numV-1)),& 
                                    -this%bufferC((this%usedXCoords(i)+offsetC)*this%numV-1)&
                                    *this%bufferInterp((this%usedXCoords(i)+offsetInterp)*this%numV-1)]
                    
                    res(k)%entry = res(k)%entry/this%vGridWidths(this%numV) 
                    
                    k = k + 1

            end do
        else

            do i = 1,size(this%usedXCoords)

                !At first v cell
                res(k)%entry = [this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+1)*(real(1,kind=rk)&
                               -this%bufferInterp(1))-this%cfAtZero(1),&
                                this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+1)&
                                *this%bufferInterp(1)-this%cfAtZero(2)]
                
                res(k)%entry = res(k)%entry/this%vGridWidths(1) 

                k = k + 1 
                do j = 2,this%numV-1
                    res(k)%entry = [-this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j-1)&
                    *(real(1,kind=rk)-this%bufferInterp(j-1)),& 
                                    -this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j-1)*this%bufferInterp(j-1)&
                                    +this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j)&
                                    *(real(1,kind=rk)-this%bufferInterp(j)),&
                                    this%bufferC((this%usedXCoords(i)+offsetC-1)*this%numV+j)*this%bufferInterp(j)]

                    res(k)%entry = res(k)%entry/this%vGridWidths(j) 

                    k = k + 1
                end do
                ! Last v cell
                res(k)%entry = [-this%bufferC((this%usedXCoords(i)+offsetC)*this%numV-1)*(real(1,kind=rk)&
                                -this%bufferInterp(this%numV-1)),& 
                                -this%bufferC((this%usedXCoords(i)+offsetC)*this%numV-1)*this%bufferInterp(this%numV-1)]

                res(k)%entry = res(k)%entry/this%vGridWidths(this%numV) 
                
                k = k + 1

            end do
    
        end if
    else
        if (allocated(this%modelboundInterp)) then 
            do i = 1,size(this%usedXCoords)

                !At first v cell
                res(k)%entry = [this%bufferC(1)*(real(1,kind=rk)&
                               -this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+1))-this%cfAtZero(1),&
                                this%bufferC(1)&
                                *this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+1)-this%cfAtZero(2)]

                res(k)%entry = res(k)%entry/this%vGridWidths(1) 
                
                k = k + 1 
                do j = 2,this%numV-1
                    res(k)%entry = [-this%bufferC(j-1)&
                    *(real(1,kind=rk)-this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+j-1)),& 
                                    -this%bufferC(j-1)*this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+j-1)&
                                    +this%bufferC(j)&
                                    *(real(1,kind=rk)-this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+j)),&
                                    this%bufferC(j)*this%bufferInterp((this%usedXCoords(i)+offsetInterp-1)*this%numV+j)]

                    res(k)%entry = res(k)%entry/this%vGridWidths(j) 

                    k = k + 1
                end do
                ! Last v cell
                res(k)%entry = [-this%bufferC(this%numV-1)*(real(1,kind=rk)&
                                -this%bufferInterp((this%usedXCoords(i)+offsetInterp)*this%numV-1)),& 
                                -this%bufferC(this%numV-1)*this%bufferInterp((i+offsetInterp)*this%numV-1)]

                res(k)%entry = res(k)%entry/this%vGridWidths(this%numV) 
                
                k = k + 1

            end do
        else
            do i = 1,size(this%usedXCoords)

                !At first v cell
                res(k)%entry = [this%bufferC(1)*(real(1,kind=rk)&
                               -this%bufferInterp(1))-this%cfAtZero(1),&
                                this%bufferC(1)&
                                *this%bufferInterp(1)-this%cfAtZero(2)]
                res(k)%entry = res(k)%entry/this%vGridWidths(1) 
                
                k = k + 1 
                do j = 2,this%numV-1
                    res(k)%entry = [-this%bufferC(j-1)&
                    *(real(1,kind=rk)-this%bufferInterp(j-1)),& 
                                    -this%bufferC(j-1)*this%bufferInterp(j-1)&
                                    +this%bufferC(j)&
                                    *(real(1,kind=rk)-this%bufferInterp(j)),&
                                    this%bufferC(j)*this%bufferInterp(j)]

                    res(k)%entry = res(k)%entry/this%vGridWidths(j) 
                    
                    k = k + 1
                end do
                ! Last v cell
                res(k)%entry = [-this%bufferC(this%numV-1)*(real(1,kind=rk)&
                                -this%bufferInterp(this%numV-1)),& 
                                -this%bufferC(this%numV-1)*this%bufferInterp(this%numV-1)]
                
                res(k)%entry = res(k)%entry/this%vGridWidths(this%numV) 
                
                k = k + 1

            end do
    
        end if
    end if

end subroutine calcDDVVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule ddv_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
