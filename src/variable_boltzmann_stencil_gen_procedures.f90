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
submodule (variable_boltzmann_stencil_gen_class) variable_boltzmann_stencil_gen_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the variable mapping Boltzmann stencil generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVarBoltzGen(this,vspaceObj,transitionIndex,lNum,absorptionTerm,superelasticTerm) 
    !! Boltzmann emission/absorption stencil value generator initialization routine with variable cross-sections and weights

    class(VariableBoltzmannStencilGen)        ,intent(inout)  :: this
    type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get grid
    integer(ik)                               ,intent(in)     :: transitionIndex !! Index of transition whose cross-section is required
    integer(ik) ,optional                     ,intent(in)     :: lNum !! Harmonic index of cross-section required. Used only if absorption term
    logical     ,optional                     ,intent(in)     :: absorptionTerm !! True if this is an absorption term stencil generator. Defaults to false.
    logical     ,optional                     ,intent(in)     :: superelasticTerm !! True if this is for a superelastic term (negative transition cost). Defaults to false.

    if (assertions) call assert(vspaceObj%isDefined(),"Undefined VSpace object passed to initVarBoltzGen")

    this%vGridCopy = vspaceObj%getVGrid()

    this%transitionIndex = transitionIndex 
    this%lNum = 0

    if (present(lNum)) this%lNum = lNum

    this%absorptionTerm = .false.

    if (present(absorptionTerm)) this%absorptionTerm = absorptionTerm

    this%superelasticTerm = .false.

    if (present(superelasticTerm)) this%superelasticTerm = superelasticTerm

    call this%makeDefined()

end subroutine initVarBoltzGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcVarBoltzVals(this,varCont,res,mbData,hostModel)
    !! Calculate Boltzmann emission/absorption stencil values in place (does not depend on varCont)

    class(VariableBoltzmannStencilGen)          ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    real(rk) ,allocatable ,dimension(:) :: csVals ,emissionVec ,energyCost

    integer(ik) :: i ,j ,inferredXSize

    if (assertions) then 
        call assert(this%isDefined(),"calcVarBoltzVals called from undefined stencil generator")
        call assert(varCont%isDefined(),"Undefined VariableContainer passed to calcVarBoltzVals")
        call assert(present(mbData),"mbData must be present when calling calcVarBoltzVals")
    end if

    select type (mbData)

    type is (ModelboundCRMData) 

        if (this%absorptionTerm) then 
            csVals =  mbData%getTransitionCrossSection(this%lNum+1,this%transitionIndex)
        else
            csVals =  mbData%getTransitionCrossSection(1,this%transitionIndex) 
        end if

        if (allocated(res)) then
            call assert(size(res)==size(csVals),"res passed to calcVarBoltzVals not of expected size")
        else
            allocate(res(size(csVals)))
        end if
        inferredXSize = size(csVals) / size(this%vGridCopy)

        energyCost = mbData%getTransitionEnergy(this%transitionIndex)
        if (assertions) then 
            call assert(size(energyCost) == inferredXSize,"energyCost vector in calcVarBoltzVals does not conform with&
                                                            & inferred spatial size from cross-section values")
            call assert(all(energyCost<0) .eqv. this%superelasticTerm,"energyCost vector in calcVarBoltzVals must have values with&
                                                                    & sign according to whether the term is superelastic or not")
        end if
        if (this%absorptionTerm) then 

            do j = 1, inferredXSize

                if (.not. allocated(this%bufferW)) then 
                    allocate(this%bufferW)
                    call this%bufferW%init(rowIndices=[(i,i=1,size(this%vGridCopy))],&
                                           colVectors=triangularIntArray(size(this%vGridCopy),this%superelasticTerm))
                    
                end if

                call mbData%interpolateW(energyCost(j),this%bufferW)

                do i = 1,size(this%bufferW%rowIndex)
                    res((j-1)*size(this%vGridCopy) + this%bufferW%rowIndex(i))%entry = this%bufferW%values(i)%entry &
                    * csVals((j-1)*size(this%vGridCopy)+ this%bufferW%columnVector(i)%entry)&
                    * this%vGridCopy(this%bufferW%columnVector(i)%entry)
                end do
            end do
        else

            do j = 1, inferredXSize
                emissionVec = mbData%getInterpolatedEmissionVector(energyCost(j))
                do i = 1,size(this%vGridCopy)
                    res((j-1) * size(this%vGridCopy) + i)%entry =&
                    [-csVals((j-1) * size(this%vGridCopy) +i)*this%vGridCopy(i)*emissionVec(i)]
                end do
            end do  
        end if
    class default 

        error stop "calcVarBoltzVals expects ModelboundCRMData for mbData"

    end select 

end subroutine calcVarBoltzVals
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule variable_boltzmann_stencil_gen_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
