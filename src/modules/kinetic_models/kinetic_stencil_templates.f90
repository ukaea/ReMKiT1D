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
module kinetic_stencil_templates
    !! author: Stefan Mijin
    !!
    !! Contains kinetic stencil template generation for use in custom model construction

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions, assertionLvl
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use variable_container_class               ,only: VariableContainer, CalculationRule
    use basic_environment_wrapper              ,only: EnvironmentWrapper
    use partition_class                        ,only: Partition
    use geometry_class                         ,only: Geometry
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use interp_stencil_gen_class               ,only: InterpStencilGenerator
    use uwc_diff_stencil_gen_class             ,only: UWCDiffStencilValGenerator
    use fb_diff_stencil_gen_class              ,only: FBDiffStencilValGenerator
    use uwc_bc_stencil_gen_class               ,only: UWCBCStencilValGenerator
    use diffusion_stencil_gen_class            ,only: DiffusionStencilValGenerator
    use derivation_abstract_class              ,only: Derivation
    use mat_derivation_abstract_class          ,only: MatDerivation
    use general_mat_term_class                 ,only: StencilTemplate
    use model_class                            ,only: Model
    use multiplicative_generator_core_class    ,only: MultiplicativeGeneratorCore
    use multiplicative_stencil_generator_class ,only: MultiplicativeStencilGen
    use ddv_stencil_gen_class                  ,only: DDVStencilGenerator 
    use v_diffusion_stencil_gen_class          ,only: VDiffStencilGen 
    use stencil_class                          ,only: Stencil
    use sparse_row_data_class                  ,only: SparseRowData
    use fixed_boltzmann_stencil_gen_class      ,only: FixedBoltzmannStencilGen
    use modelbound_CRM_data_class              ,only: ModelboundCRMData
    use modelbound_data_abstract_class         ,only: ModelboundData
    use scaling_lbc_stencil_gen_class          ,only: ScalingLBCStencilGen
    use f_scaling_derivation_class             ,only: FScalingDerivation
    use term_moment_stencil_gen_class          ,only: TermMomentStencilGenerator
    use variable_boltzmann_stencil_gen_class   ,only: VariableBoltzmannStencilGen
    use physical_constants
    use support_functions
    use support_types
    use key_names

    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initKineticStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
        !! Initialize kinetic stencil template based on existing model, environment object, and JSON file

        type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)         ,intent(inout) :: envObj
        character(*)                     ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)                     ,intent(in)    :: evolvedVar
        character(*)                     ,intent(in)    :: implicitVar
        class(ModelboundData) ,optional ,intent(in)     :: mbData

    end subroutine initKineticStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initKinDiagonalStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize diagonal stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initKinDiagonalStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initMomentStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize distribution moment stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initMomentStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initSpatialDiffStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize d/dx kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initSpatialDiffStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDDVStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize d/dv kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initDDVStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVelDiffusionStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize d^2/d^2v kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initVelDiffusionStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initIJStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize Shkarofsky I/J integral kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initIJStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initFixedBoltzmannStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
        !! Initialize fixed mapping Boltzmann kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        class(ModelboundCRMData)   ,intent(in)    :: mbData

    end subroutine initFixedBoltzmannStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initScalingLBCStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize scaling logical boundary condition kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initScalingLBCStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initFixedBoltzmannStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,mbData,&
        evolvedHarmonic,transitionIndex,fixedEnergyIndex,absorptionTerm,dbTerm)
        !! Initialize fixed mapping Boltzmann kinetic stencil template based on environment non-JSON inputs. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        class(ModelboundCRMData)   ,intent(in)    :: mbData
        integer(ik)                ,intent(in)    :: evolvedHarmonic
        integer(ik)                ,intent(in)    :: transitionIndex
        integer(ik)                ,intent(in)    :: fixedEnergyIndex
        logical                    ,intent(in)    :: absorptionTerm
        logical                    ,intent(in)    :: dbTerm

    end subroutine initFixedBoltzmannStencilDirect 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initKinDiagonalStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
        evolvedXCells,evolvedHarmonics,evolvedVCells)
        !! Initialize diagonal stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        integer(ik) ,dimension(:)  ,intent(in)    :: evolvedXCells
        integer(ik) ,dimension(:)  ,intent(in)    :: evolvedHarmonics
        integer(ik) ,dimension(:)  ,intent(in)    :: evolvedVCells

    end subroutine initKinDiagonalStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initMomentStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,momentHarmonic,momentOrder)
        !! Initialize distribution moment stencil template based on direct input data
    
        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        integer(ik)                ,intent(in)    :: momentHarmonic
        integer(ik)                ,intent(in)    :: momentOrder

    end subroutine initMomentStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initSpatialDiffStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,colHarmonic)
        !! Initialize d/dx kinetic stencil template based on direct input. 
    
        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        integer(ik)                ,intent(in)    :: rowHarmonic
        integer(ik)                ,intent(in)    :: colHarmonic

    end subroutine initSpatialDiffStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDDVStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,colHarmonic,&
        fixedC,fixedInterp,cfAtZero,modelboundC,modelboundInterp)
        !! Initialize d/dv kinetic stencil template based on direct input. 

        type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)         ,intent(inout) :: envObj
        character(*)                     ,intent(in)    :: evolvedVar
        character(*)                     ,intent(in)    :: implicitVar
        integer(ik)                      ,intent(in)    :: rowHarmonic
        integer(ik)                      ,intent(in)    :: colHarmonic
        real(rk) ,optional ,dimension(:) ,intent(in)    :: fixedC
        real(rk) ,optional ,dimension(:) ,intent(in)    :: fixedInterp
        real(rk) ,optional ,dimension(2) ,intent(in)    :: cfAtZero
        character(*) ,optional           ,intent(in)    :: modelboundC
        character(*) ,optional           ,intent(in)    :: modelboundInterp
    
    end subroutine initDDVStencilTemplateDirect 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVelDiffusionStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,&
                                                            colHarmonic,fixedA,adfAtZero,modelboundA)
        !! Initialize d^2/d^2v kinetic stencil template based on direct inputs. 

        type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)         ,intent(inout) :: envObj
        character(*)                     ,intent(in)    :: evolvedVar
        character(*)                     ,intent(in)    :: implicitVar
        integer(ik)                      ,intent(in)    :: rowHarmonic
        integer(ik)                      ,intent(in)    :: colHarmonic
        real(rk) ,optional ,dimension(:) ,intent(in)    :: fixedA
        real(rk) ,optional ,dimension(2) ,intent(in)    :: adfAtZero
        character(*) ,optional           ,intent(in)    :: modelboundA

    end subroutine initVelDiffusionStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initIJStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
        rowHarmonic,colHarmonic,integralIndex,jIntegral)
        !! Initialize Shkarofsky I/J integral kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        integer(ik)                ,intent(in)    :: rowHarmonic
        integer(ik)                ,intent(in)    :: colHarmonic
        integer(ik)                ,intent(in)    :: integralIndex
        logical                    ,intent(in)    :: jIntegral

    end subroutine initIJStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initScalingLBCStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,colHarmonic,&
        decompHarmonics,derivName,derivReqVars,leftBoundary)
        !! Initialize scaling logical boundary condition kinetic stencil template based on direct input. 

        type(StencilTemplate)           ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)        ,intent(inout) :: envObj
        character(*)                    ,intent(in)    :: evolvedVar
        character(*)                    ,intent(in)    :: implicitVar
        integer(ik)                     ,intent(in)    :: rowHarmonic
        integer(ik)                     ,intent(in)    :: colHarmonic
        integer(ik) ,dimension(:)       ,intent(in)    :: decompHarmonics
        character(*)                    ,intent(in)    :: derivName
        type(StringArray) ,dimension(:) ,intent(in)    :: derivReqVars
        logical ,optional               ,intent(in)    :: leftBoundary

    end subroutine initScalingLBCStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initTermMomentStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize term moment kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initTermMomentStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initTermMomentStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,momentOrder,colHarmonic,termName)
        !! Initialize term moment kinetic stencil template based on direct input.
    
        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        integer(ik)                ,intent(in)    :: momentOrder
        integer(ik)                ,intent(in)    :: colHarmonic
        character(*)               ,intent(in)    :: termName

    end subroutine initTermMomentStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVariableBoltzmannStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
        !! Initialize variable mapping Boltzmann kinetic stencil template based on environment object and JSON file. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        class(ModelboundCRMData)   ,intent(in)    :: mbData

    end subroutine initVariableBoltzmannStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVariableBoltzmannStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,mbData,&
        evolvedHarmonic,transitionIndex,absorptionTerm,superelasticTerm)
        !! Initialize variable mapping Boltzmann kinetic stencil template based on environment non-JSON inputs. 

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        class(ModelboundCRMData)   ,intent(in)    :: mbData
        integer(ik)                ,intent(in)    :: evolvedHarmonic
        integer(ik)                ,intent(in)    :: transitionIndex
        logical                    ,intent(in)    :: absorptionTerm
        logical                    ,intent(in)    :: superelasticTerm

    end subroutine initVariableBoltzmannStencilDirect 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module kinetic_stencil_templates
!-----------------------------------------------------------------------------------------------------------------------------------