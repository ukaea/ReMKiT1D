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
submodule (fixed_boltzmann_stencil_gen_class) fixed_boltzmann_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the fixed mapping Boltzmann stencil generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initBoltzGen(this,vspaceObj,transitionIndex,fixedWIndex,lNum,absorptionTerm,dbTerm) 
    !! Boltzmann emission/absorption stencil value generator initialization routine

    class(FixedBoltzmannStencilGen)           ,intent(inout)  :: this
    type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get grid
    integer(ik)                               ,intent(in)     :: transitionIndex !! Index of transition whose cross-section is required
    integer(ik)                               ,intent(in)     :: fixedWIndex !! Fixed inelastic mapping index associated with the transition
    integer(ik) ,optional                     ,intent(in)     :: lNum !! Harmonic index of cross-section required. Used only if absorption term
    logical     ,optional                     ,intent(in)     :: absorptionTerm !! True if this is an absorption term stencil generator. Defaults to false.
    logical     ,optional                     ,intent(in)     :: dbTerm !! True if this is a detailed balance term stencil generator. Defaults to false.

    if (assertions) call assert(vspaceObj%isDefined(),"Undefined VSpace object passed to initBoltzGen")

    this%vGridCopy = vspaceObj%getVGrid()

    this%transitionIndex = transitionIndex 
    this%fixedWIndex = fixedWIndex
    this%lNum = 0

    if (present(lNum)) this%lNum = lNum

    this%absorptionTerm = .false.

    if (present(absorptionTerm)) this%absorptionTerm = absorptionTerm

    this%dbTerm = .false. 

    if (present(dbTerm)) this%dbTerm = dbTerm

    call this%makeDefined()

end subroutine initBoltzGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcBoltzVals(this,varCont,res,mbData,hostModel)
    !! Calculate Boltzmann emission/absorption stencil values in place (does not depend on varCont)

    class(FixedBoltzmannStencilGen)             ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    real(rk) ,allocatable ,dimension(:) :: csVals ,fixedEmissionVec

    integer(ik) :: i ,j ,inferredXSize

    if (assertions) then 
        call assert(this%isDefined(),"calcBoltzVals called from undefined stencil generator")
        call assert(varCont%isDefined(),"Undefined VariableContainer passed to calcBoltzVals")
        call assert(present(mbData),"mbData must be present when calling calcBoltzVals")
    end if

    select type (mbData)

    type is (ModelboundCRMData) 

        if (this%absorptionTerm) then 
            csVals =  mbData%getTransitionCrossSection(this%lNum+1,this%transitionIndex)
        else
            csVals =  mbData%getTransitionCrossSection(1,this%transitionIndex) 
        end if

        if (allocated(res)) then
            call assert(size(res)==size(csVals),"res passed to calcBoltzVals not of expected size")
        else
            allocate(res(size(csVals)))
        end if
        inferredXSize = size(csVals) / size(this%vGridCopy)
        if (this%absorptionTerm) then 

            if (.not. allocated(this%bufferW)) this%bufferW = mbData%getFixedW(this%fixedWIndex)

            do i = 1,size(res)
                if (.not. allocated(res(i)%entry)) allocate(res(i)%entry(0))
                res(i)%entry = 0
            end do

            do j = 1, inferredXSize
                do i = 1,size(this%bufferW%rowIndex)
                    res((j-1)*size(this%vGridCopy) + this%bufferW%rowIndex(i))%entry = this%bufferW%values(i)%entry &
                    * csVals((j-1)*size(this%vGridCopy)+ this%bufferW%columnVector(i)%entry)&
                    * this%vGridCopy(this%bufferW%columnVector(i)%entry)
                end do
            end do
        else

            if (.not. this%dbTerm) then
                fixedEmissionVec = mbData%getFixedEmissionVector(this%fixedWIndex)
            else
                fixedEmissionVec = [(real(1,kind=rk),i=1,size(this%vGridCopy))]
            end if
            do j = 1, inferredXSize
                do i = 1,size(this%vGridCopy)
                    res((j-1) * size(this%vGridCopy) + i)%entry =&
                    [-csVals((j-1) * size(this%vGridCopy) +i)*this%vGridCopy(i)*fixedEmissionVec(i)]
                end do
            end do  
        end if
    class default 

        error stop "calcBoltzVals expects ModelboundCRMData for mbData"

    end select 

end subroutine calcBoltzVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule fixed_boltzmann_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
