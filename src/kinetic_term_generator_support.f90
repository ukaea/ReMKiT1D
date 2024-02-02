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
module kinetic_term_generator_support
    !! author: Stefan Mijin 
    !! 
    !! Contains routines for initializing various kinetic term generators from JSON data

    use data_kinds                                   ,only: rk, ik
    use runtime_constants                            ,only: debugging, assertions, assertionLvl
    use assertion_utility                            ,only: assert, assertIdentical, assertPure
    use term_generator_abstract_class                ,only: TermGenerator
    use crm_fixed_boltzmann_term_generator_class     ,only: CRMBoltzTermGenerator
    use crm_secondary_el_source_term_generator_class ,only: CRMSecElTermGenerator
    use crm_variable_boltzmann_term_generator_class  ,only: CRMVarBoltzTermGenerator
    use basic_environment_wrapper                    ,only: EnvironmentWrapper
    use modelbound_CRM_data_class                    ,only: ModelboundCRMData
    use modelbound_data_abstract_class               ,only: ModelboundData
    use model_class                                  ,only: Model
    use normalization_abstract_class                 ,only: Normalization
    use support_types
    use key_names

    implicit none 

    public

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCRMFixedBoltzTermGeneratorFromJSON(termGenObj,normObj,modelObj,envObj,jsonPrefix,generatorTag)
        !! Initialize term generator object as a CRMBoltzTermGenerator using the JSON config file
        
        class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
        class(Normalization)              ,intent(in)    :: normObj
        type(Model)                       ,intent(in)    :: modelObj
        type(EnvironmentWrapper)          ,intent(inout) :: envObj
        character(*)                      ,intent(in)    :: jsonPrefix
        character(*)                      ,intent(in)    :: generatorTag

    end subroutine initCRMFixedBoltzTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCRMSecElTermGeneratorFromJSON(termGenObj,modelObj,envObj,jsonPrefix,generatorTag)
        !! Initialize term generator object as a CRMSecElTermGenerator using the JSON config file
        
        class(TermGenerator) ,allocatable ,intent(inout) :: termGenObj 
        type(Model)                       ,intent(in)    :: modelObj
        type(EnvironmentWrapper)          ,intent(inout) :: envObj
        character(*)                      ,intent(in)    :: jsonPrefix
        character(*)                      ,intent(in)    :: generatorTag

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

    end subroutine initCRMVarBoltzTermGeneratorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module kinetic_term_generator_support
!-----------------------------------------------------------------------------------------------------------------------------------