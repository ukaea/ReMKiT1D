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
submodule (v_diffusion_stencil_gen_class) v_diffusion_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the d(Ad/dv)/dv stencil generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVDiffValGen(this,partitionObj,vspaceObj,procRank,activeXCoords,fixedA,mbA,adfAtZero) 
    !! d(Ad/dv)/dv stencil value generator initialization routine

    class(VDiffStencilGen)                    ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get interpolation object
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    integer(ik) ,dimension(:)                 ,intent(in)     :: activeXCoords !! List of active global x coords
    real(rk) ,optional   ,dimension(:)        ,intent(in)     :: fixedA !! Fixed A values (size numV). Defaults to ones.
    character(*) ,optional                    ,intent(in)     :: mbA !! Optional modelbound value for A (single harmonic variable). Overrides fixed values.
    real(rk) ,optional ,dimension(2)          ,intent(in)     :: adfAtZero !! Optional extrapolation of A*df/dv at zero in the form A1*f(v1)+A2*f(v2) where all A's are fixed. 

    integer(ik) :: i ,minX ,maxX

    real(rk) ,allocatable ,dimension(:) :: vGridCopy

    if (assertions) then 
        call assert(partitionObj%isDefined(),"Undefined partition passed to velocity diffusion stencil value generator constructor")
        call assert(vspaceObj%isDefined(),"Undefined VSpace passed to velocity diffusion stencil value generator constructor")
    end if

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    this%locNumX = maxX - minX + 1

    this%usedXCoords = pack(activeXCoords,activeXCoords >= minX .and. activeXCoords <= maxX) - minX + 1

    this%numV = vspaceObj%getNumV()

    this%bufferA = [(real(1,kind=rk),i=1,this%numV)]

    if (present(fixedA)) then
        if (assertions) call assert(size(fixedA)==this%numV,&
        "fixedA passed to velocity diffusion stencil value generator constructor not of size numV")

        this%bufferA = fixedA
    end if

    if (present(mbA)) then
        deallocate(this%bufferA) 
        this%modelboundA = mbA 
    end if 

    this%adfAtZero = 0
    if (present(adfAtZero)) this%adfAtZero = adfAtZero

    this%vGridWidths = vspaceObj%getVCellWidths()
    vGridCopy = vspaceObj%getVGrid()
    allocate(this%dvPlus(this%numV))
    this%dvPlus = 0
    this%dvPlus(1:this%numV-1) = vGridCopy(2:this%numV) -  vGridCopy(1:this%numV-1)

    call this%makeDefined()

end subroutine initVDiffValGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcVDiffVals(this,varCont,res,mbData,hostModel)
    !! Calculate d(Ad/dv)/dv stencil values in place (does not depend on varCont)

    class(VDiffStencilGen)                      ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    integer(ik) :: i ,j ,k ,offsetA
    if (assertions) then 
        call assert(this%isDefined(),"calcVDiffVals called from undefined stencil generator")
        call assert(varCont%isDefined(),"Undefined varCont passed to calcVDiffVals")
        if (allocated(this%modelboundA)) &
        call assert(present(mbData),"No mbData passed to calcVDiffVals when expected")
    end if

    if (allocated(res)) then
        if (assertions) call assert(size(res) == size(this%usedXCoords)*this%numV,"res passed to calcVDiffVals has unexpected size")
    else
        allocate(res(size(this%usedXCoords)*this%numV))
    end if

    if (allocated(this%modelboundA)) then 
        call mbData%copyData(this%modelboundA,this%bufferA)
        offsetA = (size(this%bufferA) - this%numV*this%locNumX)/2 + lbound(this%bufferA,1) - 1
        offsetA = offsetA/this%numV 
    end if
    
    k = 1

    if (allocated(this%modelboundA)) then 
        
            do i = 1,size(this%usedXCoords)

                !At first v cell
                res(k)%entry = [-this%bufferA((this%usedXCoords(i)+offsetA-1)*this%numV+1)/this%dvPlus(1)-this%adfAtZero(1),&
                                this%bufferA((this%usedXCoords(i)+offsetA-1)*this%numV+1)&
                                /this%dvPlus(1)-this%adfAtZero(2)]
                
                res(k)%entry = res(k)%entry/this%vGridWidths(1) 

                k = k + 1 
                do j = 2,this%numV-1
                    res(k)%entry = [this%bufferA((this%usedXCoords(i)+offsetA-1)*this%numV+j-1)/this%dvPlus(j-1),& 
                                    -this%bufferA((this%usedXCoords(i)+offsetA-1)*this%numV+j-1)/this%dvPlus(j-1)&
                                    -this%bufferA((this%usedXCoords(i)+offsetA-1)*this%numV+j)/this%dvPlus(j),&
                                    this%bufferA((this%usedXCoords(i)+offsetA-1)*this%numV+j)/this%dvPlus(j)]

                    res(k)%entry = res(k)%entry/this%vGridWidths(j) 

                    k = k + 1
                end do
                ! Last v cell
                res(k)%entry = [this%bufferA((this%usedXCoords(i)+offsetA)*this%numV-1)/this%dvPlus(this%numV-1),& 
                                -this%bufferA((this%usedXCoords(i)+offsetA)*this%numV-1)/this%dvPlus(this%numV-1)]

                res(k)%entry = res(k)%entry/this%vGridWidths(this%numV) 
                
                k = k + 1

            end do

    else
        
            do i = 1,size(this%usedXCoords)

                !At first v cell
                res(k)%entry = [-this%bufferA(1)/this%dvPlus(1)-this%adfAtZero(1),&
                                this%bufferA(1)/this%dvPlus(1)-this%adfAtZero(2)]
                res(k)%entry = res(k)%entry/this%vGridWidths(1) 
                
                k = k + 1 
                do j = 2,this%numV-1
                    res(k)%entry = [this%bufferA(j-1)/this%dvPlus(j-1),& 
                                    -this%bufferA(j-1)/this%dvPlus(j-1)&
                                    -this%bufferA(j)/this%dvPlus(j),&
                                    this%bufferA(j)/this%dvPlus(j)]

                    res(k)%entry = res(k)%entry/this%vGridWidths(j) 
                    
                    k = k + 1
                end do
                ! Last v cell
                res(k)%entry = [this%bufferA(this%numV-1)/this%dvPlus(this%numV-1),& 
                                -this%bufferA(this%numV-1)/this%dvPlus(this%numV-1)]
                
                res(k)%entry = res(k)%entry/this%vGridWidths(this%numV) 
                
                k = k + 1

            end do
    
    end if

end subroutine calcVDiffVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule v_diffusion_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
