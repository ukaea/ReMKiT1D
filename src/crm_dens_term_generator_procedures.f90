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
submodule (crm_dens_term_generator_class) crm_dens_term_generator_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the CRMDensTermGenerator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMDensTermGenerator(this,envObj,crmData,generatorTag,evolvedSpeciesIDs,includedTransitionIndices) 
    !! Constructor routine for CRM particle source term generator 

    class(CRMDensTermGenerator)         ,intent(inout) :: this 
    type(EnvironmentWrapper) ,target    ,intent(in)    :: envObj !! Environment wrapper used to get species data, partition info, etc.
    type(ModelboundCRMData)             ,intent(in)    :: crmData !! CRM data object used to get transition data
    character(*)                        ,intent(in)    :: generatorTag 
    integer(ik) ,optional ,dimension(:) ,intent(in)    :: evolvedSpeciesIDs !! Optional list of species IDs for which this 
                                                                            !! generator will build source. Defaults to all species.
    integer(ik) ,optional ,dimension(:) ,intent(in)    :: includedTransitionIndices !! Optional list of transitions to be included 
                                                                                    !! in source calculations by this generator.
                                                                                    !! Defaults to all transitions

    integer(ik) ,allocatable ,dimension(:) :: usedSpeciesIDs, usedTransitionIndices ,inStates

    integer(ik) ,allocatable ,dimension(:,:) :: popChangeMat ,densityDat
    
    integer(ik) :: i ,j ,k ,l ,numTerms ,implicitSpeciesID

    type(VarData) ,allocatable :: tempVData 

    character(len=80) :: transIndexBuffer

    if (assertions) then 

        call assert(envObj%isDefined(),"Undefined environment wrapper passed to CRMDensTermGenerator constructor")
        call assert(crmData%isDefined(),"Undefined modelbound CRM data passed to CRMDensTermGenerator constructor")

    end if

    this%envPointer => envObj

    if (present(evolvedSpeciesIDs)) then 
        usedSpeciesIDs = evolvedSpeciesIDs
    else
        usedSpeciesIDs = envObj%speciesListObj%getSpeciesIDs()
    end if

    if (present(includedTransitionIndices)) then 
        usedTransitionIndices = includedTransitionIndices
    else
        usedTransitionIndices = [(i,i=1,crmData%getNumTransitions())]
    end if

    popChangeMat = crmData%getPopulationChangeMatrix(usedSpeciesIDs,usedTransitionIndices)

    numTerms = count(popChangeMat /= 0)

    allocate(this%popChange(numTerms))
    allocate(this%vData(numTerms))
    allocate(this%evolvedVars(numTerms))
    allocate(this%implicitVars(numTerms))

    k = 1 

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
        tempVData%modelboundRowVars = [StringArray("rate0index"//trim(transIndexBuffer))]

        inStates = crmData%getTransitionIngoingStates(usedTransitionIndices(j))

        implicitSpeciesID = inStates(size(inStates))

        do i = 1,size(usedSpeciesIDs)

            if (popChangeMat(i,j) /=0 ) then

                this%popChange(k) = popChangeMat(i,j)
                this%vData(k) = tempVData
                this%evolvedVars(k)%string = envObj%speciesListObj%getSpeciesVarFromID(usedSpeciesIDs(i),1)
                this%implicitVars(k)%string = envObj%speciesListObj%getSpeciesVarFromID(implicitSpeciesID,1)
                k = k + 1
            end if

        end do
    end do 
    
    call this%setGeneratorPrefix(generatorTag)
    this%numX = envObj%gridObj%getNumX()
    call this%makeDefined()
    
end subroutine initCRMDensTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine generatePartSourceTerms(this,mbData) 
    !! Generates and allocates parent implicit particle source terms 

    class(CRMDensTermGenerator)  ,intent(inout) :: this 
    class(ModelboundData) ,optional ,intent(in) :: mbData

    type(MatTermContainer) ,allocatable, dimension(:) :: matTerms 
    type(TermContainer)    ,allocatable, dimension(:) :: genTerms

    type(GeneralMatrixTerm) ,allocatable :: termBuffer

    type(StencilTemplate) ,allocatable :: templateObj

    integer(ik) :: i, j
    character(len=30) :: intToStrBuffer
    character(len=30) :: buffer

    if (assertions) call assert(this%isDefined(),"generate called from undefined CRMDensTermGenerator")

    allocate(genTerms(0))
    allocate(matTerms(size(this%evolvedVars))) 


    do i = 1, size(this%evolvedVars)
        intToStrBuffer = ""
        write(intToStrBuffer,'(I0)') i
        call printMessage("Generating CRM density term "//this%getGeneratorPrefix()//"_implicit_"//trim(intToStrBuffer))

        call printMessage("  Term evolved var: "//this%evolvedVars(i)%string)
        call printMessage("  Term implicit var: "//this%implicitVars(i)%string)
        call printMessage("  Term modelbound row var: "//this%vData(i)%modelboundRowVars(1)%string)

        buffer=""
        do j = 1,size(this%vData(i)%rowVars)
            buffer = trim(buffer)//" "//this%vData(i)%rowVars(j)%string
        end do

        call printMessage("  Term row vars: "//buffer)
        buffer=""
        do j = 1,size(this%vData(i)%rowVars)
            intToStrBuffer = ""
            write(intToStrBuffer,'(I0)') int(this%vData(i)%rowVarPowers(j),kind=ik)
            buffer = trim(buffer)//" "//trim(intToStrBuffer)
        end do
        call printMessage("  Term row powers: "//buffer)
        if (allocated(templateObj)) deallocate(templateObj)
        allocate(templateObj)
        call initDiagonalStencilTemplateDirect(templateObj,this%envPointer,this%evolvedVars(i)%string,this%implicitVars(i)%string)
        allocate(termBuffer)
        call termBuffer%init(this%envPointer%gridObj,this%envPointer%partitionObj,this%envPointer%indexingObj,&
                            this%envPointer%mpiCont%getWorldRank(),&
                            this%evolvedVars(i)%string ,this%implicitVars(i)%string,this%envPointer%externalVars,templateObj,&
                            vData=this%vData(i),normConst=real(this%popChange(i),kind=rk))

        call move_alloc(termBuffer,matTerms(i)%entry)
    end do

    call this%setImplicitTerms(matTerms)
    call this%setGeneralTerms(genTerms)

end subroutine generatePartSourceTerms
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeCRMDensGenerator(this) 
    !! Deallocate pointer component

    type(CRMDensTermGenerator)                    ,intent(inout) :: this

    nullify(this%envPointer)

end subroutine finalizeCRMDensGenerator 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule crm_dens_term_generator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
