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
module modelbound_data_support
    !! author: Stefan Mijin
    !!
    !! Contains support for adding modelbound data to models based on JSON data

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions, assertionLvl
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use variable_container_class               ,only: VariableContainer, CalculationRule
    use variable_list_class                    ,only: VariableList
    use basic_environment_wrapper              ,only: EnvironmentWrapper
    use partition_class                        ,only: Partition
    use geometry_class                         ,only: Geometry
    use model_class                            ,only: Model
    use modelbound_data_varlike_class          ,only: ModelboundDataVarlike
    use derivation_abstract_class              ,only: Derivation ,DerivationContainer
    use modelbound_CRM_data_class              ,only: ModelboundCRMData
    use inelastic_grid_data_class              ,only: InelasticGridData
    use transition_abstract_class              ,only: Transition
    use simple_transition_class                ,only: SimpleTransition
    use fixed_ecs_transition_class             ,only: FixedECSTransition
    use variable_ecs_transition_class          ,only: VariableECSTransition
    use db_transition_class                    ,only: DBTransition
    use derived_transition_class               ,only: DerivedTransition
    use normalization_abstract_class           ,only: Normalization
    use modelbound_lbc_data_class              ,only: ModelboundLBCData
    use mat_derivation_abstract_class          ,only: MatDerivation
    use f_scaling_derivation_class             ,only: FScalingDerivation
    use janev_crm_data_support                             
    use support_types
    use key_names

    implicit none 
    private 

    public :: addModelboundDataToModel

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addModelboundDataToModel(modelObj,modelTag,envObj,normObj)
        !! Initialize modelbound data and add to corresponding model object

        type(Model)                ,intent(inout) :: modelObj
        character(*)               ,intent(in)    :: modelTag 
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        class(Normalization)       ,intent(in)    :: normObj
        
    end subroutine addModelboundDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addCustomVarlikeMBDataToModel(modelObj,modelTag,envObj)
        !! Initialize varlike modelbound data and add to corresponding model object

        type(Model)                ,intent(inout) :: modelObj
        character(*)               ,intent(in)    :: modelTag 
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        
    end subroutine addCustomVarlikeMBDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addCustomCRMDataToModel(modelObj,modelTag,envObj,normObj)
        !! Initialize custom CRM modelbound data and add to corresponding model object

        type(Model)                ,intent(inout) :: modelObj
        character(*)               ,intent(in)    :: modelTag 
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        class(Normalization)       ,intent(in)    :: normObj
        
    end subroutine addCustomCRMDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addCustomLBCDataToModel(modelObj,modelTag,envObj,normObj)
        !! Initialize custom LBC modelbound data and add to corresponding model object

        type(Model)                ,intent(inout) :: modelObj
        character(*)               ,intent(in)    :: modelTag 
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        class(Normalization)       ,intent(in)    :: normObj
        
    end subroutine addCustomLBCDataToModel
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module modelbound_data_support
!-----------------------------------------------------------------------------------------------------------------------------------