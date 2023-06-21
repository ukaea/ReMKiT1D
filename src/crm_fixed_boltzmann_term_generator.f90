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
module crm_fixed_boltzmann_term_generator_class
    !! author: Stefan Mijin
    !!
    !! Houses TermGenerator object used to generate Boltzmann terms with fixed inelastic mappings based on associated variables in species list and CRM
    !! modelbound data 

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use support_types                         
    use model_class                           ,only: Model, MatTermContainer ,TermContainer
    use modelbound_CRM_data_class             ,only: ModelboundCRMData
    use general_mat_term_class                ,only: GeneralMatrixTerm ,StencilTemplate ,VarData
    use basic_environment_wrapper             ,only: EnvironmentWrapper
    use term_generator_abstract_class         ,only: TermGenerator
    use normalization_abstract_class          ,only: Normalization
    use key_names
    use kinetic_stencil_templates

    implicit none
    private

    type ,public ,extends(TermGenerator) :: CRMBoltzTermGenerator
        !! TermGenerator object used to generate electron Boltzmann terms with fixed inelastic mappings based on associated variables in species list and CRM
        !! modelbound data. Assumes that all passed transition indices correspond to fixed inlastic mapping transitions.
        !! Assumes that the cross-section values in the transitions include any nonlinear normalization dependence (e.g. recombination etc.)

        type(VarData) ,allocatable ,dimension(:)     ,private :: vData !! Individual term variable data including densities and required modelbound data
        character(:) ,allocatable                    ,private :: distributionName !! Evolved and implicit distribution name

        integer(ik) ,allocatable ,dimension(:)       ,private :: transitionIndices !! Included transition indices
        integer(ik) ,allocatable ,dimension(:)       ,private :: fixedEnergyIndices !! Fixed energy indices associated with included transition in inelData
        logical                                      ,private :: absorptionTerms !! True if generating absorption terms. Defaults to false.
        logical                                      ,private :: dbTerms !! True if generating detailed balance terms (non-fixed stencils). Defaults to false.
        integer(ik)                                  ,private :: evolvedHarmonic !! Index of harmonic for which this generator should be generating terms
        type(EnvironmentWrapper) ,pointer            ,private :: envPointer => null()

        real(rk)                                     ,private :: normConst !! Normalization constant for the generated terms

        contains

        procedure ,public :: init => initCRMBoltzTermGenerator

        procedure ,public :: generate => generateBoltzTerms

        final             :: finalizeCRMBoltzTermGenerator

    end type CRMBoltzTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCRMBoltzTermGenerator(this,envObj,normObj,crmData,distributionName,includedTransitionIndices,&
                                                fixedEnergyIndices,evolvedHarmonic,generatorTag,absorptionTerms,dbTerms,&
                                                associatedVariableIndex) 
        !! Constructor routine for CRM Boltzmann term generator 

        class(CRMBoltzTermGenerator)        ,intent(inout) :: this 
        type(EnvironmentWrapper) ,target    ,intent(in)    :: envObj !! Environment wrapper used to get species data, partition info, etc.
        class(Normalization)                ,intent(in)    :: normObj !! Normalization object used to calculate the normalization constant for the generated terms
        type(ModelboundCRMData)             ,intent(in)    :: crmData !! CRM data object used to get transition data
        character(*)                        ,intent(in)    :: distributionName !! Evolved and implicit distribution name
        integer(ik)           ,dimension(:) ,intent(in)    :: includedTransitionIndices !! List of transitions to be included 
                                                                                        !! in calculations by this generator.
        integer(ik)           ,dimension(:) ,intent(in)    :: fixedEnergyIndices !! Fixed energy indices associated with included transitions
        integer(ik)                         ,intent(in)    :: evolvedHarmonic !! Index of harmonic for which this generator should be generating terms
        character(*)                        ,intent(in)    :: generatorTag
        logical     ,optional               ,intent(in)    :: absorptionTerms !! True if this is an absorption term stencil generator. Defaults to false.
        logical     ,optional               ,intent(in)    :: dbTerms !! True if this is a detailed balance term stencil generator. Defaults to false.
        integer(ik) ,optional               ,intent(in)    :: associatedVariableIndex !! Density index in associated variable array. Defaults to 1.

    end subroutine initCRMBoltzTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine generateBoltzTerms(this,mbData) 
        !! Generates and allocates parent implicit Boltzmann terms 
    
        class(CRMBoltzTermGenerator)    ,intent(inout) :: this 
        class(ModelboundData) ,optional ,intent(in)    :: mbData

    end subroutine generateBoltzTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeCRMBoltzTermGenerator(this) 
        !! Deallocate pointer component

        type(CRMBoltzTermGenerator)                    ,intent(inout) :: this

    end subroutine finalizeCRMBoltzTermGenerator 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module crm_fixed_boltzmann_term_generator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 