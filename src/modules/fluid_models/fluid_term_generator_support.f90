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
module fluid_term_generator_support
    !! author: Stefan Mijin 
    !! 
    !! Contains routines for initializing various fluid term generators from JSON data

    use data_kinds                         ,only: rk, ik
    use runtime_constants                  ,only: debugging, assertions, assertionLvl
    use assertion_utility                  ,only: assert, assertIdentical, assertPure
    use term_generator_abstract_class      ,only: TermGenerator
    use crm_dens_term_generator_class      ,only: CRMDensTermGenerator
    use crm_el_energy_term_generator_class ,only: CRMElEnergyTermGenerator
    use basic_environment_wrapper          ,only: EnvironmentWrapper
    use modelbound_CRM_data_class          ,only: ModelboundCRMData
    use modelbound_data_abstract_class     ,only: ModelboundData
    use model_class                        ,only: Model
    use support_types
    use key_names

    implicit none 

    public

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCRMDensTermGeneratorFromJSON(termGenObj,modelObj,envObj,jsonPrefix,generatorTag)
        !! Initialize term generator object as a CRMDensTermGenerator using the JSON config file
        
        class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
        type(Model)                       ,intent(in)    :: modelObj
        type(EnvironmentWrapper)          ,intent(inout) :: envObj
        character(*)                      ,intent(in)    :: jsonPrefix
        character(*)                      ,intent(in)    :: generatorTag

    end subroutine initCRMDensTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCRMElEnergyTermGeneratorFromJSON(termGenObj,modelObj,envObj,jsonPrefix,generatorTag)
        !! Initialize term generator object as a CRMElEnergyTermGenerator using the JSON config file
        
        class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
        type(Model)                       ,intent(in)    :: modelObj
        type(EnvironmentWrapper)          ,intent(inout) :: envObj
        character(*)                      ,intent(in)    :: jsonPrefix
        character(*)                      ,intent(in)    :: generatorTag

    end subroutine initCRMElEnergyTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module fluid_term_generator_support
!-----------------------------------------------------------------------------------------------------------------------------------