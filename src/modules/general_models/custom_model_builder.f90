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
module custom_model_builder_class
    !! author: Stefan Mijin
    !!
    !! Houses model builder for user-defined models

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions ,assertionLvl
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use support_types                         
    use mpi_controller_class                  ,only: MPIController
    use json_controller_class                 ,only: JSONController
    use modeller_class                        ,only: Modeller
    use model_builder_abstract_class          ,only: ModelBuilder
    use grid_class                            ,only: Grid 
    use partition_class                       ,only: Partition 
    use geometry_class                        ,only: Geometry 
    use indexing_class                        ,only: Indexing 
    use variable_container_class              ,only: VariableContainer
    use model_class                           ,only: Model
    use matrix_term_abstract_class            ,only: MatrixTerm
    use basic_environment_wrapper             ,only: EnvironmentWrapper
    use normalization_abstract_class          ,only: Normalization
    use general_mat_term_class                ,only: VarData ,GeneralMatrixTerm ,StencilTemplate ,CoordProfiles
    use term_generator_abstract_class         ,only: TermGenerator ,TermGeneratorContainer
    use signal_abstract_class                 ,only: TimeSignalData
    use status_printing
    use fluid_term_generator_support
    use kinetic_term_generator_support
    use fluid_stencil_templates
    use kinetic_stencil_templates
    use modelbound_data_support
    use derivation_explicit_term_class        ,only: DerivationTerm
    use term_abstract_class                   ,only: Term

    use key_names

    implicit none
    private

    type ,public ,extends(ModelBuilder) :: CustomModelBuilder
        !! Builder for custom models with user defined terms and modelbound data

        class(Model) ,allocatable ,private :: modelBuffer !! Model buffer filled during initialization

        contains

        procedure ,public :: init => initCustomBuilder

        procedure ,public :: addModelToModeller => addCustomModel
        
        procedure ,private :: addTermToModel

        procedure ,private :: applyTermGenerator

    end type CustomModelBuilder
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCustomBuilder(this,env,normObject,modelTag) 
        !! Constructs the model for this builder and sets it into defined state
    
        class(CustomModelBuilder)                 ,intent(inout)  :: this 
        class(EnvironmentWrapper)                 ,intent(inout)  :: env 
        class(Normalization)                      ,intent(in)     :: normObject !! Reference normalization object
        character(*)                              ,intent(in)     :: modelTag !! Tag of this model

        end subroutine initCustomBuilder
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine addCustomModel(this,modellerObj) 
            !! Adds the model built by the builder and resets the builder to become undefined for further use

            class(CustomModelBuilder)           ,intent(inout) :: this 
            class(Modeller)                     ,intent(inout) :: modellerObj

        end subroutine addCustomModel
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine addTermToModel(this,termJSONPrefix,termTag,env,normObject,mbData) 
            !! Adds individual term to model buffer based on JSON file data
        
            class(CustomModelBuilder)                 ,intent(inout)  :: this 
            character(*)                              ,intent(in)     :: termJSONPrefix 
            character(*)                              ,intent(in)     :: termTag 
            class(EnvironmentWrapper)                 ,intent(inout)  :: env 
            class(Normalization)                      ,intent(in)     :: normObject
            class(ModelboundData) ,optional           ,intent(in)     :: mbData

        end subroutine addTermToModel
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine applyTermGenerator(this,env,normObject,modelTag,currentNumTerms,mbData) 
            !! Checks for associated term generator in JSON file and applies to custom model
        
            class(CustomModelBuilder)                 ,intent(inout)  :: this 
            class(EnvironmentWrapper)                 ,intent(inout)  :: env 
            class(Normalization)                      ,intent(in)     :: normObject !! Reference normalization object
            character(*)                              ,intent(in)     :: modelTag !! Tag of this model
            integer(ik)                               ,intent(in)     :: currentNumTerms
            class(ModelboundData) ,optional           ,intent(in)     :: mbData
    
        end subroutine applyTermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module custom_model_builder_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
