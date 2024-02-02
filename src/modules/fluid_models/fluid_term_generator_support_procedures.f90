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
submodule(fluid_term_generator_support) fluid_term_generator_support_procedures
!! author: Stefan Mijin 
!! 
!! Contains implementation of fluid term generator support routines 

implicit none

contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMDensTermGeneratorFromJSON(termGenObj,modelObj,envObj,jsonPrefix,generatorTag)
    !! Initialize term generator object as a CRMDensTermGenerator using the JSON config file
    
    class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
    type(Model)                       ,intent(in)    :: modelObj
    type(EnvironmentWrapper)          ,intent(inout) :: envObj
    character(*)                      ,intent(in)    :: jsonPrefix
    character(*)                      ,intent(in)    :: generatorTag

    class(ModelboundData) ,allocatable :: mbData
    
    type(CRMDensTermGenerator) ,allocatable :: termGenBuffer

    type(NamedIntegerArray) ,dimension(1) :: evolvedSpeciesIDs, includedTransitionIDs

    call modelObj%copyModelData(mbData)

    if (assertions .or. assertionLvl >= 0) call assert(allocated(mbData),&
    "initCRMDensTermGeneratorFromJSON unable to detect modelbound data in modelObj")

    select type (mbData)
    class is (ModelboundCRMData)
        allocate(termGenBuffer)
        
        evolvedSpeciesIDs(1)%name = jsonPrefix//"."//keyEvolvedSpeciesIDs
        allocate(evolvedSpeciesIDs(1)%values(0))

        includedTransitionIDs(1)%name = jsonPrefix//"."//keyIncludedTransitionInds
        allocate(includedTransitionIDs(1)%values(0))

        call envObj%jsonCont%load(evolvedSpeciesIDs)
        call envObj%jsonCont%load(includedTransitionIDs)

        if (size(evolvedSpeciesIDs(1)%values) > 0) then 
            if (size(includedTransitionIDs(1)%values) > 0) then
                call termGenBuffer%init(envObj,mbData,generatorTag,evolvedSpeciesIDs=evolvedSpeciesIDs(1)%values,&
                                        includedTransitionIndices=includedTransitionIDs(1)%values)
            else
                call termGenBuffer%init(envObj,mbData,generatorTag,evolvedSpeciesIDs=evolvedSpeciesIDs(1)%values)

            end if
        else
            if (size(includedTransitionIDs(1)%values) > 0) then
                call termGenBuffer%init(envObj,mbData,generatorTag,&
                                        includedTransitionIndices=includedTransitionIDs(1)%values)
            else
                call termGenBuffer%init(envObj,mbData,generatorTag)
            end if
        end if

        call move_alloc(termGenBuffer,termGenObj)
    class default 
        error stop "initCRMDensTermGeneratorFromJSON detected unsupported data "
    end select

end subroutine initCRMDensTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCRMElEnergyTermGeneratorFromJSON(termGenObj,modelObj,envObj,jsonPrefix,generatorTag)
    !! Initialize term generator object as a CRMElEnergyTermGenerator using the JSON config file
    
    class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
    type(Model)                       ,intent(in)    :: modelObj
    type(EnvironmentWrapper)          ,intent(inout) :: envObj
    character(*)                      ,intent(in)    :: jsonPrefix
    character(*)                      ,intent(in)    :: generatorTag

    class(ModelboundData) ,allocatable :: mbData
    
    type(CRMElEnergyTermGenerator) ,allocatable :: termGenBuffer

    type(NamedIntegerArray) ,dimension(1) :: includedTransitionIDs
    type(NamedString)       ,dimension(1) :: electronEnergyVarName 

    call modelObj%copyModelData(mbData)

    if (assertions .or. assertionLvl >= 0) call assert(allocated(mbData),&
    "initCRMDElEnergyTermGeneratorFromJSON unable to detect modelbound data in modelObj")

    select type (mbData)
    class is (ModelboundCRMData)
        allocate(termGenBuffer)

        includedTransitionIDs(1)%name = jsonPrefix//"."//keyIncludedTransitionInds
        allocate(includedTransitionIDs(1)%values(0))

        electronEnergyVarName(1) = NamedString(jsonPrefix//"."//keyElectronEnergyVar,"")

        call envObj%jsonCont%load(electronEnergyVarName)
        call envObj%jsonCont%load(includedTransitionIDs)

        if (size(includedTransitionIDs(1)%values) > 0) then
            call termGenBuffer%init(envObj,mbData,electronEnergyVarName(1)%value,generatorTag,&
                                    includedTransitionIndices=includedTransitionIDs(1)%values)
        else
            call termGenBuffer%init(envObj,mbData,generatorTag,electronEnergyVarName(1)%value)
        end if

        call move_alloc(termGenBuffer,termGenObj)
    class default 
        error stop "initCRMElEnergyTermGeneratorFromJSON detected unsupported data "
    end select

end subroutine initCRMElEnergyTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------  
!-----------------------------------------------------------------------------------------------------------------------------------  
end submodule fluid_term_generator_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------  
