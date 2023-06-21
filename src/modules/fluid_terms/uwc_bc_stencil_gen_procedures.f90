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
submodule (uwc_bc_stencil_gen_class) uwc_bc_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the flux-like boundary stencil value generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initUWCBCStencil(this,partitionObj,procRank,innerJ,outerJ,linInterp,extrapolate,&
    interpVarIndex,leftBoundary,lowerBoundVarIndex,fixedLowerBound,linExterp) 
    !! Flux-like boundary condition stencil value generator initialization routine

    class(UWCBCStencilValGenerator)           ,intent(inout)  :: this
    type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
    integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
    real(rk)                                  ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star at boundary face
    real(rk)                                  ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star at boundary cell
    real(rk)                                  ,intent(in)     :: linInterp !! Linear interpolation coefficient at face before boundary
    logical ,optional                         ,intent(in)     :: extrapolate !! Extrapolate implicit/column variables. Defaults to .false.
    logical ,optional                         ,intent(in)     :: leftBoundary !! True if the stencil values are for the left boundary condition
    integer(ik) ,optional                     ,intent(in)     :: interpVarIndex !! Variable to optionally interpolate as part of stencil (flux jacobian)
    integer(ik) ,optional                     ,intent(in)     :: lowerBoundVarIndex !! Variable index of lower jacobian bound variable
    real(rk)    ,optional                     ,intent(in)     :: fixedLowerBound !! A constant lower bound value (ignored if lowerBoundVarIndex is allocated)
    real(rk)    ,optional                     ,intent(in)     :: linExterp !! Linear extrapolation coefficient for column variables. Computed from linInterp by default

    if (assertions) call assert(partitionObj%isDefined(),&
    "Undefined partition object passed to flux-like boundary val generator constructor")

    if (assertions) then 

        if (present(lowerBoundVarIndex)) call assert(present(interpVarIndex),"If lowerBoundVarIndex passed to initUWCBCStencil&
            & the interpolated/flux jacobian variable index must be passed as well")

        if (present(fixedLowerBound)) call assert(present(interpVarIndex),"If fixedLowerBound passed to initUWCBCStencil&
            & the interpolated/flux jacobian variable index must be passed as well")
    end if

    this%minX = partitionObj%getMinXAtInd(procRank+1)
    this%maxX = partitionObj%getMaxXAtInd(procRank+1)
    this%inferredGridSize = maxval(partitionObj%getMaxX())

    this%innerJ = innerJ 
    this%outerJ = outerJ
    this%linInterp = linInterp

    if (present(interpVarIndex)) allocate(this%interpVarIndex,source=interpVarIndex) 

    if (present(lowerBoundVarIndex)) allocate(this%lowerBoundVarIndex,source=lowerBoundVarIndex)

    this%extrapolate = .false.

    if (present(extrapolate)) this%extrapolate = extrapolate

    this%leftBoundary = .false. 

    if (present(leftBoundary)) this%leftBoundary = leftBoundary

    this%linExterp = this%linInterp
    if (.not. this%leftBoundary) this%linExterp = real(1,kind=rk) - this%linInterp

    if (present(linExterp)) this%linExterp = linExterp

    if (present(fixedLowerBound)) allocate(this%fixedLowerBound,source=fixedLowerBound)

    call this%makeDefined()

end subroutine initUWCBCStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcUWCBCVals(this,varCont,res,mbData,hostModel)
    !! Calculate flux-like boundary stencil values in place 

    class(UWCBCStencilValGenerator)             ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    logical :: containsLeftBoundary ,containsRightBoundary

    real(rk) :: fluxJacAtBoundary ,lowerBoundAtBoundary

    if (assertions) then 
        call assert(this%isDefined(),"calcUWCBCVals called from undefined stencil value generator")
        call assert(varCont%isDefined(),"Undefined variable container passed to calcUWCBCVals")
    end if

    containsLeftBoundary = this%minX == 1
    containsRightBoundary = this%maxX == this%inferredGridSize

    fluxJacAtBoundary = real(1,kind=rk)

    if (containsLeftBoundary .and. this%leftBoundary) then 

        if (allocated(res)) deallocate(res)
        allocate(res(1))

        if (allocated(this%interpVarIndex)) &
        fluxJacAtBoundary = varCont%variables(this%interpVarIndex)%entry(1) * (real(1,kind=rk)+this%linInterp) &
                          - this%linInterp * varCont%variables(this%interpVarIndex)%entry(2)

        if (allocated(this%lowerBoundVarIndex)) then 
            lowerBoundAtBoundary = varCont%variables(this%lowerBoundVarIndex)%entry(1) * (real(1,kind=rk)+this%linInterp) &
            - this%linInterp * varCont%variables(this%lowerBoundVarIndex)%entry(2)

            fluxJacAtBoundary = min(fluxJacAtBoundary,-lowerBoundAtBoundary) !min because of normal sign

        else if (allocated(this%fixedLowerBound)) then 
            fluxJacAtBoundary = min(fluxJacAtBoundary,-this%fixedLowerBound) !min because of normal sign
        end if

        if (this%extrapolate) then 
            res(1)%entry = -this%outerJ*this%innerJ*fluxJacAtBoundary*[real(1,kind=rk)+this%linExterp,-this%linExterp]
        else
            res(1)%entry = -this%outerJ*this%innerJ*fluxJacAtBoundary*real([1,0],kind=rk)
        end if
    else if (containsRightBoundary .and. .not. this%leftBoundary) then 

        if (allocated(res)) deallocate(res)
        allocate(res(1))

        if (allocated(this%interpVarIndex)) &
        fluxJacAtBoundary = varCont%variables(this%interpVarIndex)%entry(this%maxX-this%minX+1) &
                          * (real(2,kind=rk)-this%linInterp) &
                          - (real(1,kind=rk) - this%linInterp) * varCont%variables(this%interpVarIndex)%entry(this%maxX-this%minX)

        if (allocated(this%lowerBoundVarIndex)) then 
            lowerBoundAtBoundary = varCont%variables(this%lowerBoundVarIndex)%entry(this%maxX-this%minX+1) &
            * (real(2,kind=rk)-this%linInterp) &
            - (real(1,kind=rk) - this%linInterp) * varCont%variables(this%lowerBoundVarIndex)%entry(this%maxX-this%minX)

            fluxJacAtBoundary = max(fluxJacAtBoundary,lowerBoundAtBoundary)
        else if (allocated(this%fixedLowerBound)) then 
            fluxJacAtBoundary = max(fluxJacAtBoundary,this%fixedLowerBound) 
        end if

        if (this%extrapolate) then 
            res(1)%entry = &
            this%outerJ*this%innerJ*fluxJacAtBoundary*[real(1,kind=rk)+this%linExterp,-this%linExterp]
        else
            res(1)%entry = this%outerJ*this%innerJ*fluxJacAtBoundary*real([1,0],kind=rk)
        end if

    else

        if (allocated(res)) deallocate(res)
        allocate(res(0))

    end if

end subroutine calcUWCBCVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule uwc_bc_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
