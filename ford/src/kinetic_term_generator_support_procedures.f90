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
submodule(kinetic_term_generator_support) kinetic_term_generator_support_procedures
!! author: Stefan Mijin 
!! 
!! Contains implementation of kinetic term generator support routines 

implicit none

contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMFixedBoltzTermGeneratorFromJSON(termGenObj,normObj,modelObj,envObj,jsonPrefix,generatorTag)
    !! Initialize term generator object as a CRMBoltzTermGenerator using the JSON config file
    
    class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
    class(Normalization)              ,intent(in)    :: normObj
    type(Model)                       ,intent(in)    :: modelObj
    type(EnvironmentWrapper)          ,intent(inout) :: envObj
    character(*)                      ,intent(in)    :: jsonPrefix
    character(*)                      ,intent(in)    :: generatorTag

    class(ModelboundData) ,allocatable :: mbData
    
    type(CRMBoltzTermGenerator) ,allocatable :: termGenBuffer

    type(NamedIntegerArray) ,dimension(1) :: includedTransitionIDs ,fixedEnergyIndices

    type(NamedString) ,dimension(1) :: distributionName 

    type(NamedLogical) ,dimension(1) :: absorptionTerm, dbTerm
 
    type(NamedInteger) ,dimension(1) :: evolvedHarmonic ,associatedVarIndex

    call modelObj%copyModelData(mbData)

    if (assertions) call assert(allocated(mbData),&
    "initCRMFixedBoltzTermGeneratorFromJSON unable to detect modelbound data in modelObj")

    select type (mbData)
    class is (ModelboundCRMData)
        allocate(termGenBuffer)
        
        fixedEnergyIndices(1)%name = jsonPrefix//"."//keyFixedEnergyIndices
        allocate(fixedEnergyIndices(1)%values(0))

        includedTransitionIDs(1)%name = jsonPrefix//"."//keyIncludedTransitionInds
        allocate(includedTransitionIDs(1)%values(0))

        absorptionTerm(1) = NamedLogical(jsonPrefix//"."//keyAbsorptionTerm,.false.)
        dbTerm(1) = NamedLogical(jsonPrefix//"."//keyDBTerm,.false.)

        distributionName(1) = NamedString(jsonPrefix//"."//keyDistributionVarName,"")

        evolvedHarmonic(1) = NamedInteger(jsonPrefix//"."//keyEvolvedHarmonicSingle,0)

        associatedVarIndex(1) = NamedInteger(jsonPrefix//"."//keyAssociatedVarIndex,1)

        call envObj%jsonCont%load(fixedEnergyIndices)
        call envObj%jsonCont%load(includedTransitionIDs)
        call envObj%jsonCont%load(absorptionTerm)
        call envObj%jsonCont%load(dbTerm)
        call envObj%jsonCont%load(distributionName)
        call envObj%jsonCont%load(evolvedHarmonic)
        call envObj%jsonCont%load(associatedVarIndex)

        if (assertions) then 
            call assert(size(fixedEnergyIndices(1)%values) > 0,fixedEnergyIndices(1)%name//" must have nonzero size")
            call assert(size(fixedEnergyIndices(1)%values) > 0,includedTransitionIDs(1)%name//" must have nonzero size")
            call assert(envObj%externalVars%isVarNameRegistered(distributionName(1)%value),&
                       distributionName(1)%name//" not registered in environment object")

            call assert(evolvedHarmonic(1)%value>0,evolvedHarmonic(1)%name//" must be positive")
            call assert(evolvedHarmonic(1)%value<=envObj%gridObj%getNumH(),evolvedHarmonic(1)%name//" out of bounds")

        end if

        call termGenBuffer%init(envObj,normObj,mbData,distributionName(1)%value,includedTransitionIDs(1)%values,&
                                fixedEnergyIndices(1)%values,evolvedHarmonic(1)%value,generatorTag,&
                                absorptionTerm(1)%value,dbTerm(1)%value,associatedVarIndex(1)%value)

        call move_alloc(termGenBuffer,termGenObj)
    class default 
        error stop "initCRMFixedBoltzTermGeneratorFromJSON detected unsupported data "
    end select

end subroutine initCRMFixedBoltzTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMSecElTermGeneratorFromJSON(termGenObj,modelObj,envObj,jsonPrefix,generatorTag)
    !! Initialize term generator object as a CRMSecElTermGenerator using the JSON config file
    
    class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
    type(Model)                       ,intent(in)    :: modelObj
    type(EnvironmentWrapper)          ,intent(inout) :: envObj
    character(*)                      ,intent(in)    :: jsonPrefix
    character(*)                      ,intent(in)    :: generatorTag

    class(ModelboundData) ,allocatable :: mbData
    
    type(CRMSecElTermGenerator) ,allocatable :: termGenBuffer

    type(NamedIntegerArray) ,dimension(1) :: includedTransitionIDs
    type(NamedString)       ,dimension(1) :: distributionName

    call modelObj%copyModelData(mbData)

    if (assertions) call assert(allocated(mbData),&
    "initCRMSecElTermGeneratorFromJSON unable to detect modelbound data in modelObj")

    select type (mbData)
    class is (ModelboundCRMData)
        allocate(termGenBuffer)

        includedTransitionIDs(1)%name = jsonPrefix//"."//keyIncludedTransitionInds
        allocate(includedTransitionIDs(1)%values(0))

        distributionName(1) = NamedString(jsonPrefix//"."//keyDistributionVarName,"")

        call envObj%jsonCont%load(distributionName)
        call envObj%jsonCont%load(includedTransitionIDs)
        if (assertions) then 
            call assert(envObj%externalVars%isVarNameRegistered(distributionName(1)%value),&
                       distributionName(1)%name//" not registered in environment object")
        end if

        if (size(includedTransitionIDs(1)%values) > 0) then
            call termGenBuffer%init(envObj,mbData,distributionName(1)%value,generatorTag,&
                                    includedTransitionIndices=includedTransitionIDs(1)%values)
        else
            call termGenBuffer%init(envObj,mbData,distributionName(1)%value,generatorTag)
        end if

        call move_alloc(termGenBuffer,termGenObj)
    class default 
        error stop "initCRMSecElTermGeneratorFromJSON detected unsupported data "
    end select

end subroutine initCRMSecElTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMVarBoltzTermGeneratorFromJSON(termGenObj,normObj,modelObj,envObj,jsonPrefix,generatorTag)
    !! Initialize term generator object as a CRMVarBoltzTermGenerator using the JSON config file
    
    class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
    class(Normalization)              ,intent(in)    :: normObj
    type(Model)                       ,intent(in)    :: modelObj
    type(EnvironmentWrapper)          ,intent(inout) :: envObj
    character(*)                      ,intent(in)    :: jsonPrefix
    character(*)                      ,intent(in)    :: generatorTag

    class(ModelboundData) ,allocatable :: mbData
    
    type(CRMVarBoltzTermGenerator) ,allocatable :: termGenBuffer

    type(NamedIntegerArray) ,dimension(1) :: includedTransitionIDs 

    type(NamedString) ,dimension(1) :: distributionName 

    type(NamedLogical) ,dimension(1) :: absorptionTerm, superelasticTerm
 
    type(NamedInteger) ,dimension(1) :: evolvedHarmonic ,associatedVarIndex

    call modelObj%copyModelData(mbData)

    if (assertions) call assert(allocated(mbData),&
    "initCRMVarBoltzTermGeneratorFromJSON unable to detect modelbound data in modelObj")

    select type (mbData)
    class is (ModelboundCRMData)
        allocate(termGenBuffer)

        includedTransitionIDs(1)%name = jsonPrefix//"."//keyIncludedTransitionInds
        allocate(includedTransitionIDs(1)%values(0))

        absorptionTerm(1) = NamedLogical(jsonPrefix//"."//keyAbsorptionTerm,.false.)
        superelasticTerm(1) = NamedLogical(jsonPrefix//"."//keySuperelasticTerm,.false.)

        distributionName(1) = NamedString(jsonPrefix//"."//keyDistributionVarName,"")

        evolvedHarmonic(1) = NamedInteger(jsonPrefix//"."//keyEvolvedHarmonicSingle,0)

        associatedVarIndex(1) = NamedInteger(jsonPrefix//"."//keyAssociatedVarIndex,1)

        call envObj%jsonCont%load(includedTransitionIDs)
        call envObj%jsonCont%load(absorptionTerm)
        call envObj%jsonCont%load(superelasticTerm)
        call envObj%jsonCont%load(distributionName)
        call envObj%jsonCont%load(evolvedHarmonic)
        call envObj%jsonCont%load(associatedVarIndex)

        if (assertions) then 
            call assert(envObj%externalVars%isVarNameRegistered(distributionName(1)%value),&
                       distributionName(1)%name//" not registered in environment object")

            call assert(evolvedHarmonic(1)%value>0,evolvedHarmonic(1)%name//" must be positive")
            call assert(evolvedHarmonic(1)%value<=envObj%gridObj%getNumH(),evolvedHarmonic(1)%name//" out of bounds")

        end if

        call termGenBuffer%init(envObj,normObj,mbData,distributionName(1)%value,includedTransitionIDs(1)%values,&
                                evolvedHarmonic(1)%value,generatorTag,&
                                absorptionTerm(1)%value,superelasticTerm(1)%value,associatedVarIndex(1)%value)

        call move_alloc(termGenBuffer,termGenObj)
    class default 
        error stop "initCRMVarBoltzTermGeneratorFromJSON detected unsupported data "
    end select

end subroutine initCRMVarBoltzTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------  
!-----------------------------------------------------------------------------------------------------------------------------------  
end submodule kinetic_term_generator_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------  
