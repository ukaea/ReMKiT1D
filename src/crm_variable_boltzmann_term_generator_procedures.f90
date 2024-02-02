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
submodule (crm_variable_boltzmann_term_generator_class) crm_variable_boltzmann_term_generator_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the CRMBVarBoltzTermGenerator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMVarBoltzTermGenerator(this,envObj,normObj,crmData,distributionName,includedTransitionIndices,&
    evolvedHarmonic,generatorTag,absorptionTerms,superelasticTerms,&
    associatedVariableIndex) 
    !! Constructor routine for CRM Variable Boltzmann term generator 

    class(CRMVarBoltzTermGenerator)     ,intent(inout) :: this 
    type(EnvironmentWrapper) ,target    ,intent(in)    :: envObj !! Environment wrapper used to get species data, partition info, etc.
    class(Normalization)                ,intent(in)    :: normObj !! Normalization object used to calculate the normalization constant for the generated terms
    type(ModelboundCRMData)             ,intent(in)    :: crmData !! CRM data object used to get transition data
    character(*)                        ,intent(in)    :: distributionName !! Evolved and implicit distribution name
    integer(ik)           ,dimension(:) ,intent(in)    :: includedTransitionIndices !! List of transitions to be included 
                                                !! in calculations by this generator.
    integer(ik)                         ,intent(in)    :: evolvedHarmonic !! Index of harmonic for which this generator should be generating terms
    character(*)                        ,intent(in)    :: generatorTag
    logical     ,optional               ,intent(in)    :: absorptionTerms !! True if this is an absorption term stencil generator. Defaults to false.
    logical     ,optional               ,intent(in)    :: superelasticTerms !! True if this is a superelastic term stencil generator. Defaults to false.
    integer(ik) ,optional               ,intent(in)    :: associatedVariableIndex !! Density index in associated variable array. Defaults to 1.

    integer(ik) ,allocatable ,dimension(:,:) :: densityDat
    
    integer(ik) :: i ,j ,usedAVarInd

    type(VarData) ,allocatable :: tempVData 

    if (assertions .or. assertionLvl >= 0) then 

        call assert(envObj%isDefined(),"Undefined environment wrapper passed to CRMVarBoltzTermGenerator constructor")
        call assert(crmData%isDefined(),"Undefined modelbound CRM data passed to CRMVarBoltzTermGenerator constructor")

    end if

    this%envPointer => envObj

    this%distributionName = distributionName

    this%evolvedHarmonic = evolvedHarmonic

    this%absorptionTerms = .false. 
    
    if (present(absorptionTerms)) this%absorptionTerms = absorptionTerms

    this%superelasticTerms = .false. 

    if (present(superelasticTerms)) this%superelasticTerms = superelasticTerms

    this%transitionIndices = includedTransitionIndices

    allocate(this%vData(size(includedTransitionIndices)))

    usedAVarInd = 1
    if (present(associatedVariableIndex)) usedAVarInd = associatedVariableIndex
    do i = 1,size(includedTransitionIndices)

        if (allocated(tempVData)) deallocate(tempVData)
        allocate(tempVData)
        densityDat = crmData%getRequiredDensityData(includedTransitionIndices(i))

        allocate(tempVData%rowVars(size(densityDat,1)))
        do j = 1,size(tempVData%rowVars)
            !Assumes density is associated as the first variable
            tempVData%rowVars(j)%string = envObj%speciesListObj%getSpeciesVarFromID(densityDat(j,1),usedAVarInd) 
        end do

        tempVData%rowVarPowers = real(densityDat(:,2),kind=rk)

        this%vData(i) = tempVData

    end do 

    this%normConst = normObj%getNormalizationValue(keyDensNorm) * & 
                     normObj%getNormalizationValue(keyVelGridNorm) * & 
                     normObj%getNormalizationValue(keyTimeNorm) * &
                     normObj%getNormalizationValue(keyCrossSectionNorm)
    
    call this%setGeneratorPrefix(generatorTag)
    call this%makeDefined()

end subroutine initCRMVarBoltzTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine generateVarBoltzTerms(this,mbData) 
    !! Generates and allocates parent implicit variable Boltzmann terms 

    class(CRMVarBoltzTermGenerator) ,intent(inout) :: this 
    class(ModelboundData) ,optional ,intent(in)    :: mbData

    type(MatTermContainer) ,allocatable, dimension(:) :: matTerms 
    type(TermContainer)    ,allocatable, dimension(:) :: genTerms

    type(GeneralMatrixTerm) ,allocatable :: termBuffer

    type(StencilTemplate) ,allocatable :: templateObj

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(this%isDefined(),"generate called from undefined CRMVarBoltzTermGenerator")
        call assert(present(mbData) ,"generate on CRMVarBoltzTermGenerator requires passed modelbound data")
    end if
    allocate(genTerms(0))
    allocate(matTerms(size(this%vData)))

    do i = 1, size(this%vData)

        if (allocated(templateObj)) deallocate(templateObj)
        allocate(templateObj)

        select type (mbData)
        type is (ModelboundCRMData)
            call initVariableBoltzmannStencilDirect(templateObj,this%envPointer,this%distributionName,this%distributionName,mbData,&
                 this%evolvedHarmonic,this%transitionIndices(i),this%absorptionTerms,this%superelasticTerms) 
        class default 
            error stop "generate on CRMVarBoltzTermGenerator expects ModelboundCRMData type"
        end select

        allocate(termBuffer)
        call termBuffer%init(this%envPointer%gridObj,this%envPointer%partitionObj,this%envPointer%indexingObj,&
                            this%envPointer%mpiCont%getWorldRank(),&
                            this%distributionName,this%distributionName,this%envPointer%externalVars,templateObj,&
                            vData=this%vData(i),normConst=this%normConst,mbData=mbData)

        call move_alloc(termBuffer,matTerms(i)%entry)
    end do

    call this%setImplicitTerms(matTerms)
    call this%setGeneralTerms(genTerms)

end subroutine generateVarBoltzTerms
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeCRMVarBoltzTermGenerator(this) 
    !! Deallocate pointer component

    type(CRMVarBoltzTermGenerator)                    ,intent(inout) :: this

    nullify(this%envPointer)

end subroutine finalizeCRMVarBoltzTermGenerator 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule crm_variable_boltzmann_term_generator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
