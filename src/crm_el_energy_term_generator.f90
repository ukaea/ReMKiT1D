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
module crm_el_energy_term_generator_class
    !! author: Stefan Mijin
    !!
    !! Houses TermGenerator object used to generate electron energy implicit source terms based on associated variables in species list and CRM
    !! modelbound data 

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions, assertionLvl
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use support_types                         
    use model_class                           ,only: Model, MatTermContainer ,TermContainer
    use modelbound_CRM_data_class             ,only: ModelboundCRMData
    use general_mat_term_class                ,only: GeneralMatrixTerm ,StencilTemplate ,VarData
    use basic_environment_wrapper             ,only: EnvironmentWrapper
    use term_generator_abstract_class         ,only: TermGenerator
    use modelbound_data_abstract_class        ,only: ModelboundData
    use fluid_stencil_templates

    implicit none
    private

    type ,public ,extends(TermGenerator) :: CRMElEnergyTermGenerator
        !! TermGenerator object used to generate electronj energy source and sink terms based on associated variables in species list and CRM
        !! modelbound data. Assumes rates are normalized so that the resulting source term is normalized to densNorm*eVTemperature/timeNorm. 

        type(VarData) ,allocatable ,dimension(:)     ,private :: vData !! Individual term variable data including densities and required modelbound data
        type(StringArray) ,allocatable ,dimension(:) ,private :: implicitVars !! Implicit density names (corresponding to last inState in transitions)
        character(:) ,allocatable                    ,private :: electronEnergyVarName !! Name of implicit electron energy variable

        type(EnvironmentWrapper) ,pointer            ,private :: envPointer => null()

        integer(ik) ,private :: numX !! Local copy of x-grid size

        contains

        procedure ,public :: init => initCRMElEnergyTermGenerator

        procedure ,public :: generate => generateElEnergySourceTerms

        final             :: finalizeCRMElEnergyGenerator

    end type CRMElEnergyTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
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

    end subroutine initCRMElEnergyTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine generateElEnergySourceTerms(this,mbData) 
        !! Generates and allocates parent implicit energy source terms 

        class(CRMElEnergyTermGenerator)  ,intent(inout) :: this 
        class(ModelboundData) ,optional  ,intent(in)    :: mbData

    end subroutine generateElEnergySourceTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeCRMElEnergyGenerator(this) 
        !! Deallocate pointer component

        type(CRMElEnergyTermGenerator)                    ,intent(inout) :: this

    end subroutine finalizeCRMElEnergyGenerator 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module crm_el_energy_term_generator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 