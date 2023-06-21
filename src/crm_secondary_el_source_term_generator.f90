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
module crm_secondary_el_source_term_generator_class
    !! author: Stefan Mijin
    !!
    !! Houses TermGenerator object used to generate secondary electron source/sink terms for based on associated variables in species list and CRM
    !! modelbound data 

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use support_types                         
    use model_class                           ,only: Model, MatTermContainer ,TermContainer
    use modelbound_CRM_data_class             ,only: ModelboundCRMData
    use general_mat_term_class                ,only: GeneralMatrixTerm ,StencilTemplate ,VarData ,CoordProfiles
    use basic_environment_wrapper             ,only: EnvironmentWrapper
    use term_generator_abstract_class         ,only: TermGenerator
    use modelbound_data_abstract_class        ,only: ModelboundData
    use kinetic_stencil_templates

    implicit none
    private

    type ,public ,extends(TermGenerator) :: CRMSecElTermGenerator
        !! TermGenerator object used to generate secondary electron source/sink terms based on associated variables in species list and CRM
        !! modelbound data. 

        integer(ik) ,allocatable ,dimension(:)       ,private :: popChange !! Population change associated with each source term reaction
        type(VarData) ,allocatable ,dimension(:)     ,private :: vData !! Individual term variable data including densities and required modelbound data
        character(:) ,allocatable                    ,private :: distributionName !! Evolved distribution name
        type(StringArray) ,allocatable ,dimension(:) ,private :: implicitVars !! Implicit density names (corresponding to last inState in transitions)

        type(EnvironmentWrapper) ,pointer            ,private :: envPointer => null()

        contains

        procedure ,public :: init => initCRMSecElTermGenerator

        procedure ,public :: generate => generateSecElSourceTerms

        final             :: finalizeCRMSecElGenerator

    end type CRMSecElTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
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

    end subroutine initCRMSecElTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine generateSecElSourceTerms(this,mbData) 
        !! Generates and allocates parent implicit secondary electron source/sink terms 

        class(CRMSecElTermGenerator)  ,intent(inout) :: this 
        class(ModelboundData) ,optional ,intent(in) :: mbData

    end subroutine generateSecElSourceTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeCRMSecElGenerator(this) 
        !! Deallocate pointer component

        type(CRMSecElTermGenerator)                    ,intent(inout) :: this

    end subroutine finalizeCRMSecElGenerator 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module crm_secondary_el_source_term_generator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 