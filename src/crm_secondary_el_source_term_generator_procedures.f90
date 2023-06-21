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
submodule (crm_secondary_el_source_term_generator_class) crm_secondary_el_source_term_generator_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the CRMSecElTermGenerator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMSecElTermGenerator(this,envObj,crmData,distributionName,generatorTag,includedTransitionIndices) 
    !! Constructor routine for CRM secondary electron source/sink term generator 

    class(CRMSecElTermGenerator)        ,intent(inout) :: this 
    type(EnvironmentWrapper) ,target    ,intent(in)    :: envObj !! Environment wrapper used to get species data, partition info, etc.
    type(ModelboundCRMData)             ,intent(in)    :: crmData !! CRM data object used to get transition data
    character(*)                        ,intent(in)    :: distributionName !! Evolved distribution name
    character(*)                        ,intent(in)    :: generatorTag
    integer(ik) ,optional ,dimension(:) ,intent(in)    :: includedTransitionIndices !! Optional list of transitions to be included 
                                                                                    !! in source calculations by this generator.
                                                                                    !! Defaults to all transitions

    integer(ik) ,allocatable ,dimension(:) ::  usedTransitionIndices ,inStates

    integer(ik) ,allocatable ,dimension(:,:) :: popChangeMat ,densityDat
    
    integer(ik) :: i,l,k ,numTerms ,implicitSpeciesID

    type(VarData) ,allocatable :: tempVData 

    character(len=80) :: transIndexBuffer

    if (assertions) then 

        call assert(envObj%isDefined(),"Undefined environment wrapper passed to CRMDensTermGenerator constructor")
        call assert(crmData%isDefined(),"Undefined modelbound CRM data passed to CRMDensTermGenerator constructor")

    end if

    this%envPointer => envObj

    if (present(includedTransitionIndices)) then 
        usedTransitionIndices = includedTransitionIndices
    else
        usedTransitionIndices = [(i,i=1,crmData%getNumTransitions())]
    end if
    popChangeMat = crmData%getPopulationChangeMatrix([0],usedTransitionIndices)

    numTerms = count(popChangeMat /= 0)

    allocate(this%popChange(numTerms))
    allocate(this%vData(numTerms))
    allocate(this%implicitVars(numTerms))

    k=1
    do i = 1,size(usedTransitionIndices)

        if (allocated(tempVData)) deallocate(tempVData)
        allocate(tempVData)
        densityDat = crmData%getRequiredDensityData(usedTransitionIndices(i),removeLastState=.true.)

        allocate(tempVData%rowVars(size(densityDat,1)))
        do l = 1,size(tempVData%rowVars)
            !Assumes density is associated as the first variable
            tempVData%rowVars(l)%string = envObj%speciesListObj%getSpeciesVarFromID(densityDat(l,1),1) 
        end do

        tempVData%rowVarPowers = real(densityDat(:,2),kind=rk)
        write(transIndexBuffer,'(I0)') usedTransitionIndices(i)
        tempVData%modelboundRowVars = [StringArray("rate0index"//trim(transIndexBuffer))]

        inStates = crmData%getTransitionIngoingStates(usedTransitionIndices(i))

        implicitSpeciesID = inStates(size(inStates))

        if (popChangeMat(1,i) /=0 ) then

            this%popChange(k) = popChangeMat(1,i)
            this%vData(k) = tempVData
            this%implicitVars(k)%string = envObj%speciesListObj%getSpeciesVarFromID(implicitSpeciesID,1)
            k=k+1
        end if

    end do 

    this%distributionName = distributionName

    call this%setGeneratorPrefix(generatorTag)
    call this%makeDefined()
    
end subroutine initCRMSecElTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine generateSecElSourceTerms(this,mbData) 
    !! Generates and allocates parent implicit secondary electron source/sink terms 

    class(CRMSecElTermGenerator)  ,intent(inout) :: this 
    class(ModelboundData) ,optional ,intent(in) :: mbData

    type(MatTermContainer) ,allocatable, dimension(:) :: matTerms 
    type(TermContainer)    ,allocatable, dimension(:) :: genTerms

    type(GeneralMatrixTerm) ,allocatable :: termBuffer

    type(StencilTemplate) ,allocatable :: templateObj

    type(CoordProfiles) :: cProfs

    integer(ik) :: i

    integer(ik) ,dimension(0) :: dummyXCells

    if (assertions) call assert(this%isDefined(),"generate called from undefined CRMDensTermGenerator")

    allocate(genTerms(0))
    allocate(matTerms(size(this%vData)))

    cProfs%vProfile = 1/(4 * pi * this%envPointer%vSpaceObj%getVGrid()**2 * this%envPointer%vSpaceObj%getVCellWidths())

    do i = 1, size(this%vData)

        if (allocated(templateObj)) deallocate(templateObj)
        allocate(templateObj)

        call initKinDiagonalStencilTemplateDirect(templateObj,this%envPointer,this%distributionName,&
                                                  this%implicitVars(i)%string,dummyXCells,[1],[1])
        allocate(termBuffer)
        call termBuffer%init(this%envPointer%gridObj,this%envPointer%partitionObj,this%envPointer%indexingObj,&
                            this%envPointer%mpiCont%getWorldRank(),&
                            this%distributionName ,this%implicitVars(i)%string,this%envPointer%externalVars,templateObj,&
                            vData=this%vData(i),normConst=real(this%popChange(i),kind=rk),coordProfile=cProfs,mbData=mbData)

        call move_alloc(termBuffer,matTerms(i)%entry)
    end do

    call this%setImplicitTerms(matTerms)
    call this%setGeneralTerms(genTerms)

end subroutine generateSecElSourceTerms
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeCRMSecElGenerator(this) 
!! Deallocate pointer component

    type(CRMSecElTermGenerator)                    ,intent(inout) :: this
    
    nullify(this%envPointer)

end subroutine finalizeCRMSecElGenerator 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule crm_secondary_el_source_term_generator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
