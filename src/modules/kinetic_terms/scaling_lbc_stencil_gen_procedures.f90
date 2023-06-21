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
submodule (scaling_lbc_stencil_gen_class) scaling_lbc_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the scaling logical boundary condition stencil generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initScalingLBCGen(this,vspaceObj,isActive,fDeriv,fDerivReqIndices,colL,dx,includedDecompHarmonics) 
    !! Scaling logical boundary condition value generator initialization routine

    class(ScalingLBCStencilGen)               ,intent(inout)  :: this
    type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get grid data 
    logical                                   ,intent(in)     :: isActive !! Set to true if the process containing this has a boundary where this is to be applied 
    type(FScalingDerivation)                  ,intent(in)     :: fDeriv !! Derivation used to get scaling factors 
    integer(ik) ,dimension(:)                 ,intent(in)     :: fDerivReqIndices !! Indices required for scaling extrapolation derivation
    integer(ik)                               ,intent(in)     :: colL !! Column l harmonic
    real(rk)                                  ,intent(in)     :: dx !! Boundary cell width
    integer(ik) ,optional ,dimension(:)       ,intent(in)     :: includedDecompHarmonics  !! Harmonics to be included in the decomposition for this generator (useful for right boundary on staggered grid). Defaults to all harmonics.

    integer(ik) :: i 

    if (assertions) then 
        call assert(vspaceObj%isDefined(),"Undefined VSpace object passed to initScalingLBCGen")
        call assert(fDeriv%isDefined(),"Undefined FScalingDerivation object passed to initScalingLBCGen")
    end if

    this%isActive = isActive 

    this%fDeriv = fDeriv
    this%fDerivReqVarIndices = fDerivReqIndices

    this%colL = colL 

    this%vGrid = vspaceObj%getVGrid()
    this%vGridWidths = vspaceObj%getVCellWidths()

    allocate(this%bufferPl(size(this%vGrid),vspaceObj%getNumH()))

    if (present(includedDecompHarmonics)) then 
        this%includedDecompHarmonics = includedDecompHarmonics
    else
        this%includedDecompHarmonics = [(i,i=1,vspaceObj%getNumH())]
    end if

    this%dx = dx
    call this%makeDefined()

end subroutine initScalingLBCGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcScalingLBCVals(this,varCont,res,mbData,hostModel)
    !! Calculate scaling logical boundary condition values in place 

    class(ScalingLBCStencilGen)                  ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    integer(ik) :: i ,j ,coCell 

    real(rk) ,dimension(2) :: interpCoords, interpWidths ,interpCoeffs ,scalingFactors

    if (assertions) then 
        call assert(this%isDefined(),"calcScalingLBCVals called from undefined stencil generator")
        call assert(varCont%isDefined(),"Undefined VariableContainer passed to calcScalingLBCVals")
        call assert(present(mbData),"mbData must be present when calling calcScalingLBCVals")
        call assert(mbData%isDefined(),"Undefined mbData passed to calcScalingLBCVals")
    end if

    if (this%isActive) then

        if (.not. allocated(res)) then 
            allocate(res(size(this%vGrid)))
            do i = 1,size(res)-1
                allocate(res(i)%entry(2*size(this%includedDecompHarmonics)))
            end do
            allocate(res(size(res))%entry(size(this%includedDecompHarmonics)))
        else
            call assert(size(res)==size(this%vGrid),"res passed to calcScalingLBCVals not of expected size")
            do i = 1,size(res)-1
                call assert(size(res(i)%entry)==2*size(this%includedDecompHarmonics),&
                "res passed to calcScalingLBCVals not of expected size")
            end do
            call assert(size(res(size(res))%entry)==size(this%includedDecompHarmonics),&
                "res passed to calcScalingLBCVals not of expected size")
        end if

        select type (mbData)
        type is (ModelboundLBCData)

            this%bufferPl = mbData%getPll(this%colL)
            interpCoeffs = mbData%getInterpCoeffs()
            interpCoords = mbData%getInterpCoords()
            interpWidths = mbData%getInterpWidths()
            coCell = mbData%getCoCell()
            scalingFactors = this%fDeriv%getScalingFactors(varCont%variables,this%fDerivReqVarIndices)/this%dx

            !Cells before interpolation
            do i = 1,coCell - 2
                res(i)%entry = 0
                do j = 1,size(this%includedDecompHarmonics)
                    res(i)%entry(2*j-1) = scalingFactors(1+mod(this%includedDecompHarmonics(j)-1,2))&
                                     *this%vGrid(i)*this%bufferPl(i,this%includedDecompHarmonics(j))
                end do
            end do

            !Cells around interpolation
            do j = 1,size(this%includedDecompHarmonics)
                res(coCell-1)%entry(2*j-1) = scalingFactors(1+mod(this%includedDecompHarmonics(j)-1,2))&
                                         *interpCoords(1)**2*interpWidths(1)*interpCoeffs(1)&
                                         *this%bufferPl(coCell-1,this%includedDecompHarmonics(j))&
                                         /(this%vGrid(coCell-1)*this%vGridWidths(coCell-1))
                                         
                res(coCell-1)%entry(2*j) = scalingFactors(1+mod(this%includedDecompHarmonics(j)-1,2))&
                                         *interpCoords(1)**2*interpWidths(1)&
                                         *(real(1,kind=rk) - interpCoeffs(1))&
                                         *this%bufferPl(coCell-1,this%includedDecompHarmonics(j))&
                                         /(this%vGrid(coCell-1)*this%vGridWidths(coCell-1))
            end do

            do j = 1,size(this%includedDecompHarmonics)

                res(coCell)%entry(2*j-1) = scalingFactors(1+mod(this%includedDecompHarmonics(j)-1,2))&
                                            *interpCoords(2)**2*interpWidths(2)&
                                            *interpCoeffs(2)*this%bufferPl(coCell,this%includedDecompHarmonics(j))&
                                            /(this%vGrid(coCell)*this%vGridWidths(coCell))

               res(coCell)%entry(2*j) = scalingFactors(1+mod(this%includedDecompHarmonics(j)-1,2))&
                                    *interpCoords(2)**2*interpWidths(2)*(real(1,kind=rk)-interpCoeffs(2))&
                                    *this%bufferPl(coCell,this%includedDecompHarmonics(j))&
                                    /(this%vGrid(coCell)*this%vGridWidths(coCell))

            end do

            ! Cells after interpolation 

            do i = coCell+1,size(res)-1
                res(i)%entry = 0
                do j = 1,size(this%includedDecompHarmonics)
                    res(i)%entry(2*j-1) = scalingFactors(1+mod(this%includedDecompHarmonics(j)-1,2))&
                                      *this%vGrid(i)*this%bufferPl(i,this%includedDecompHarmonics(j))
                end do
            end do

            do j = 1,size(this%includedDecompHarmonics)
                res(size(res))%entry(j) = scalingFactors(1+mod(this%includedDecompHarmonics(j)-1,2))&
                                  *this%vGrid(size(res))*this%bufferPl(size(res),this%includedDecompHarmonics(j))
            end do

        class default 

            error stop "unsupported modelbound data detected in calcScalingLBCVals"

        end select

    else
        if (.not. allocated(res)) then 
            allocate(res(0))
        else
            call assert(size(res)==0,"res passed to calcScalingLBCVals not of expected size")
        end if
    end if

end subroutine calcScalingLBCVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule scaling_lbc_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
