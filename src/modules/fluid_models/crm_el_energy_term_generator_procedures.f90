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
submodule (crm_el_energy_term_generator_class) crm_dens_term_generator_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the CRMElEnergyTermGenerator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMElEnergyTermGenerator(this,envObj,crmData,generatorTag,electronEnergyVarName,includedTransitionIndices) 
    !! Constructor routine for CRM energy source/sink term generator 

    class(CRMElEnergyTermGenerator)     ,intent(inout) :: this 
    type(EnvironmentWrapper) ,target    ,intent(in)    :: envObj !! Environment wrapper used to get species data, partition info, etc.
    type(ModelboundCRMData)             ,intent(in)    :: crmData !! CRM data object used to get transition data
    character(*)                        ,intent(in)    :: electronEnergyVarName !! Name of electron energy variable (used as evolved var)
    character(*)                        ,intent(in)    :: generatorTag
    integer(ik) ,optional ,dimension(:) ,intent(in)    :: includedTransitionIndices !! Optional list of transitions to be included 
                                                                                    !! in source calculations by this generator.
                                                                                    !! Defaults to all transitions that include electrons in ingoingStates

    integer(ik) ,allocatable ,dimension(:) :: usedTransitionIndices ,inStates

    integer(ik) ,allocatable ,dimension(:,:) :: densityDat
    
    integer(ik) :: i ,j ,l ,numTerms ,implicitSpeciesID

    type(VarData) ,allocatable :: tempVData 

    character(len=80) :: transIndexBuffer

    if (assertions) then 

        call assert(envObj%isDefined(),"Undefined environment wrapper passed to CRMElEnergyTermGenerator constructor")
        call assert(crmData%isDefined(),"Undefined modelbound CRM data passed to CRMElEnergyTermGenerator constructor")
        call assert(envObj%externalVars%isVarNameRegistered(electronEnergyVarName),"electronEnergyVarName not registered in &
        &environment wrapper passed to CRMElEnergyTermGenerator constructor")
    end if

    this%envPointer => envObj

    this%electronEnergyVarName = electronEnergyVarName
    allocate(usedTransitionIndices(0))
    if (present(includedTransitionIndices)) then 
        do i = 1,size(includedTransitionIndices)
            if (any(crmData%getTransitionIngoingStates(includedTransitionIndices(i)) == crmData%getElState())) &
            usedTransitionIndices = [usedTransitionIndices,includedTransitionIndices(i)]
        end do
    else

        do i = 1,crmData%getNumTransitions()
            if (any(crmData%getTransitionIngoingStates(i) == crmData%getElState())) &
            usedTransitionIndices = [usedTransitionIndices,i]
        end do
    end if

    numTerms = size(usedTransitionIndices)

    allocate(this%vData(numTerms))
    allocate(this%implicitVars(numTerms))

    do j = 1,size(usedTransitionIndices)

        if (allocated(tempVData)) deallocate(tempVData)
        allocate(tempVData)
        densityDat = crmData%getRequiredDensityData(usedTransitionIndices(j),removeLastState=.true.)

        allocate(tempVData%rowVars(size(densityDat,1)))
        do l = 1,size(tempVData%rowVars)
            !Assumes density is associated as the first variable
            tempVData%rowVars(l)%string = envObj%speciesListObj%getSpeciesVarFromID(densityDat(l,1),1) 
        end do

        tempVData%rowVarPowers = real(densityDat(:,2),kind=rk)
        write(transIndexBuffer,'(I0)') usedTransitionIndices(j)
        tempVData%modelboundRowVars = [StringArray("rate2index"//trim(transIndexBuffer))]

        inStates = crmData%getTransitionIngoingStates(usedTransitionIndices(j))

        implicitSpeciesID = inStates(size(inStates))

        this%vData(j) = tempVData
        this%implicitVars(j)%string = envObj%speciesListObj%getSpeciesVarFromID(implicitSpeciesID,1)
               
    end do 
    
    this%numX = envObj%gridObj%getNumX()

    call this%setGeneratorPrefix(generatorTag)
    call this%makeDefined()
    
end subroutine initCRMElEnergyTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine generateElEnergySourceTerms(this,mbData) 
    !! Generates and allocates parent implicit energy source terms 

    class(CRMElEnergyTermGenerator)  ,intent(inout) :: this 
    class(ModelboundData) ,optional  ,intent(in)    :: mbData

    type(MatTermContainer) ,allocatable, dimension(:) :: matTerms 
    type(TermContainer)    ,allocatable, dimension(:) :: genTerms

    type(GeneralMatrixTerm) ,allocatable :: termBuffer
    type(StencilTemplate) , allocatable :: templateObj

    integer(ik) :: i

    if (assertions) call assert(this%isDefined(),"generate called from undefined CRMElEnergyTermGenerator")

    allocate(genTerms(0))
    allocate(matTerms(size(this%vData)))

    do i = 1, size(this%vData)

        if (allocated(templateObj)) deallocate(templateObj)
        allocate(templateObj)
        call initDiagonalStencilTemplateDirect(templateObj,this%envPointer,this%electronEnergyVarName,this%implicitVars(i)%string)
        allocate(termBuffer)
        call termBuffer%init(this%envPointer%gridObj,this%envPointer%partitionObj,this%envPointer%indexingObj,&
                            this%envPointer%mpiCont%getWorldRank(),&
                            this%electronEnergyVarName ,this%implicitVars(i)%string,this%envPointer%externalVars,templateObj,&
                            vData=this%vData(i),normConst=real(-1,kind=rk)) ! -1 normConst since inelastic collisions have positive transition energies

        call move_alloc(termBuffer,matTerms(i)%entry)
    end do

    call this%setImplicitTerms(matTerms)
    call this%setGeneralTerms(genTerms)

end subroutine generateElEnergySourceTerms
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeCRMElEnergyGenerator(this) 
!! Deallocate pointer component

    type(CRMElEnergyTermGenerator)                    ,intent(inout) :: this

    nullify(this%envPointer)

end subroutine finalizeCRMElEnergyGenerator 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule crm_dens_term_generator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
