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
submodule (diffusion_stencil_gen_class) diffusion_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the diffusion stencil value generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDiffusionStencilGen(this,partitionObj,procRank,innerJ,outerJ,linInterp&
    ,xPeriodic,diffCoeffDerivRule,xHaloWidth,doNotInterpolateD) 
    !! Diffusion stencil value generator initialization routine

    class(DiffusionStencilValGenerator)       ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    real(rk) ,dimension(:)                    ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star (should conform to x-grid)
    real(rk) ,dimension(:)                    ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star (should conform to x-grid)
    real(rk) ,dimension(:)                    ,intent(in)     :: linInterp !! Linear interpolation coefficients to right cell faces (should be size(xGrid)+1)
    logical ,optional                         ,intent(in)     :: xPeriodic !! Used to determine if outer processors should have their stencils trimmed. Defaults to .false.  
    type(CalculationRule) ,optional           ,intent(in)     :: diffCoeffDerivRule !! Rule for deriving the optional diffusion coefficient
    integer(ik) ,optional                     ,intent(in)     :: xHaloWidth !! Halo width in the x direction. Defaults to 1 and must always be >0
    logical ,optional                         ,intent(in)     :: doNotInterpolateD !! Assume the diffusion coefficient is already calcuated at cell boundaries. Defaults to False.

    integer(ik) :: inferredGridSize ,usedHaloWidth

    if (assertions .or. assertionLvl >= 0) call assert(partitionObj%isDefined(),&
    "Undefined partition object passed to central diff val generator constructor")

    inferredGridSize = maxval(partitionObj%getMaxX())

    if (assertions .or. assertionLvl >= 0) then 
        call assert(size(innerJ) == inferredGridSize,&
        "innerJ passed to initDiffusionStencilGen does not conform to inferred grid size")
        call assert(size(outerJ) == inferredGridSize,&
        "outerJ passed to initDiffusionStencilGen does not conform to inferred grid size")
        call assert(size(linInterp) == inferredGridSize+1,&
        "linInterp passed to initDiffusionStencilGen does not conform to inferred grid size")

        if (present(diffCoeffDerivRule)) call assert(allocated(diffCoeffDerivRule%derivationMethod),&
        "diffCoeffDerivRule passed to initDiffusionStencilGen must have an allocated derivationMethod")

        if (present(xHaloWidth)) call assert(xHaloWidth > 0,"xHaloWidth passed to initDiffusionStencilGen must be greater than 0")
    end if

    this%periodicGrid = .false. 

    if (present(xPeriodic)) this%periodicGrid = xPeriodic

    usedHaloWidth = 1 

    if (present(xHaloWidth)) usedHaloWidth = xHaloWidth

    this%minX = partitionObj%getMinXAtInd(procRank+1)
    this%maxX = partitionObj%getMaxXAtInd(procRank+1)

    this%innerJ = innerJ 
    this%outerJ = outerJ
    allocate(this%linInterp(0:this%maxX-this%minX+1))
    this%linInterp = linInterp(this%minX:this%maxX+1)

    if (present(diffCoeffDerivRule)) allocate(this%diffCoeffDerivRule,source=diffCoeffDerivRule)

    allocate(this%diffCoeffBuffer(1-usedHaloWidth:this%maxX-this%minX+1+usedHaloWidth))
    allocate(this%diffCoeffBufferInterp(0:this%maxX-this%minX+1))
    this%diffCoeffBuffer = real(1,kind=rk)
    this%diffCoeffBufferInterp = real(1,kind=rk)

    this%doNotInterpolateD = .false.
    if (present(doNotInterpolateD)) this%doNotInterpolateD = doNotInterpolateD

    call this%makeDefined()

end subroutine initDiffusionStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcDiffusionStencilGen(this,varCont,res,mbData,hostModel)
    !! Calculate diffusion stencil values in place 

    class(DiffusionStencilValGenerator)         ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    logical :: containsLeftBoundary ,containsRightBoundary

    integer(ik) :: i

    if (assertions) then 
        call assert(this%isDefined(),"calcDiffusionStencilGen called from undefined stencil value generator")
        call assert(varCont%isDefined(),"Undefined variable container passed to calcDiffusionStencilGen")
    end if

    containsLeftBoundary = this%minX == 1
    containsRightBoundary = this%maxX == size(this%innerJ)

    if (allocated(res)) then 
        if (assertions) &
        call assert(size(res) == this%maxX - this%minX + 1,&
        "res passed to calcDiffusionStencilGen does not conform to local x-grid size")
    else
        allocate(res(this%maxX - this%minX + 1))
    end if
    if (allocated(this%diffCoeffDerivRule)) then 

        if (.not. allocated(this%diffCoeffDerivIndices)) then 
            allocate(this%diffCoeffDerivIndices(size(this%diffCoeffDerivRule%requiredVarNames)))
            do i = 1, size(this%diffCoeffDerivRule%requiredVarNames)
                this%diffCoeffDerivIndices(i) = varCont%getVarIndex(this%diffCoeffDerivRule%requiredVarNames(i)%string)
            end do
        end if

        this%diffCoeffBuffer = this%diffCoeffDerivRule%derivationMethod%calculate(varCont%variables,this%diffCoeffDerivIndices)

        if (this%doNotInterpolateD) then

            this%diffCoeffBufferInterp = this%diffCoeffBuffer 

        else
            this%diffCoeffBufferInterp = &
            this%diffCoeffBuffer(0:size(res))*(1-this%linInterp) &
            + this%diffCoeffBuffer(1:size(res)+1)*this%linInterp
        end if

    end if
    if (containsLeftBoundary) then 
        if (this%periodicGrid) then 
            res(1)%entry = &
            this%outerJ(this%minX)*&
            [this%innerJ(size(this%innerJ))*this%diffCoeffBufferInterp(0),& ! left cell
            -this%innerJ(this%minX)*this%diffCoeffBufferInterp(1)& !central cell
            -this%innerJ(size(this%innerJ))*this%diffCoeffBufferInterp(0),&
            this%innerJ(this%minX)*this%diffCoeffBufferInterp(1)] !right cell
        else
            res(1)%entry = &
            this%outerJ(this%minX)*&
            [-this%innerJ(this%minX)*this%diffCoeffBufferInterp(1),& !central cell
            this%innerJ(this%minX)*this%diffCoeffBufferInterp(1)] !right cell
        end if
    else 
        res(1)%entry = &
            this%outerJ(this%minX)*&
            [this%innerJ(this%minX-1)*this%diffCoeffBufferInterp(0),& ! left cell
            -this%innerJ(this%minX)*this%diffCoeffBufferInterp(1)& !central cell
            -this%innerJ(this%minX-1)*this%diffCoeffBufferInterp(0),&
            this%innerJ(this%minX)*this%diffCoeffBufferInterp(1)] !right cell
    end if

    do i = 1,this%maxX - this%minX -1
        res(i+1)%entry = &
            this%outerJ(this%minX + i)*&
            [this%innerJ(this%minX + i -1)*this%diffCoeffBufferInterp(i),& ! left cell
            -this%innerJ(this%minX + i) *this%diffCoeffBufferInterp(i+1)& !central cell
            -this%innerJ(this%minX + i -1)*this%diffCoeffBufferInterp(i),&
            this%innerJ(this%minX + i )*this%diffCoeffBufferInterp(i+1)] !right cell
    end do

    if (this%periodicGrid .or. .not. containsRightBoundary) then 
        res(size(res))%entry = &
        this%outerJ(this%maxX)*&
        [this%innerJ(this%maxX-1)*this%diffCoeffBufferInterp(size(res)-1),& ! left cell
        -this%innerJ(this%maxX) *this%diffCoeffBufferInterp(size(res))& !central cell
        -this%innerJ(this%maxX -1)*this%diffCoeffBufferInterp(size(res)-1),&
        this%innerJ(this%maxX )*this%diffCoeffBufferInterp(size(res))] !right cell
    else
        res(size(res))%entry = &
        this%outerJ(this%maxX)*&
        [this%innerJ(this%maxX-1)*this%diffCoeffBufferInterp(size(res)-1),& ! left cell
        -this%innerJ(this%maxX -1)*this%diffCoeffBufferInterp(size(res)-1)] ! central cell
    end if

end subroutine calcDiffusionStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule diffusion_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
