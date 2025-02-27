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
module fluid_stencil_templates
    !! author: Stefan Mijin
    !!
    !! Contains stencil template generation for use in custom model construction

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
    use fluid_gen1d_class                      ,only: FluidStencilGen1D
    use stencil_generator1d_class              ,only: StencilGenerator1D
    use derivation_abstract_class              ,only: Derivation
    use general_mat_term_class                 ,only: StencilTemplate
    use stencil1d_class                        ,only: Stencil1D
    use kinetic_stencil_templates
    use support_types
    use key_names

    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initFluidStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initFluidStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDiagonalStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize diagonal stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initDiagonalStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDiagonalStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,evolvedXCells)
        !! Initialize diagonal stencil template based on direct inputs

        type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)         ,intent(inout) :: envObj
        character(*)                     ,intent(in)    :: evolvedVar
        character(*)                     ,intent(in)    :: implicitVar
        integer(ik) ,optional ,dimension(:) ,intent(in) :: evolvedXCells

    end subroutine initDiagonalStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCentralDifferenceInterpTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize central difference with face interpolation stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initCentralDifferenceInterpTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCentralDifferenceInterpTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                                ignoreJacobian,interpolatedVarName)
        !! Initialize central difference with face interpolation stencil template based on direct inputs

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        logical                    ,intent(in)    :: ignoreJacobian
        character(*) ,optional     ,intent(in)    :: interpolatedVarName

    end subroutine initCentralDifferenceInterpTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initStaggeredDifferenceTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize staggered difference stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initStaggeredDifferenceTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initStaggeredDifferenceTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian)
        !! Initialize staggered difference stencil template based on direct inputs

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        logical                    ,intent(in)    :: ignoreJacobian

    end subroutine initStaggeredDifferenceTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initUpwindingDifferenceTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize upwinding difference stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initUpwindingDifferenceTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initUpwindingDifferenceTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                            ignoreJacobian,fluxJacVar)
        !! Initialize upwinding difference with stencil template based on direct inputs
    
        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        logical                    ,intent(in)    :: ignoreJacobian
        character(*)               ,intent(in)    :: fluxJacVar

    end subroutine initUpwindingDifferenceTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initBCTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize flux-like boundary condition stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initBCTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initBCTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian,leftBoundary,&
        dontExtrapolate,noLowerBound,fixedLowerBound,fluxJacVar,lowerBoundVar)
        !! Initialize flux-like boundary condition stencil template based on direct inputs 
                
        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar
        logical                    ,intent(in)    :: ignoreJacobian
        logical                    ,intent(in)    :: leftBoundary
        logical                    ,intent(in)    :: dontExtrapolate
        logical                    ,intent(in)    :: noLowerBound
        real(rk)                   ,intent(in)    :: fixedLowerBound
        character(*) ,optional     ,intent(in)    :: fluxJacVar
        character(*) ,optional     ,intent(in)    :: lowerBoundVar

    end subroutine initBCTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDiffusionStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
        !! Initialize diffusion stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
        character(*)               ,intent(in)    :: evolvedVar
        character(*)               ,intent(in)    :: implicitVar

    end subroutine initDiffusionStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDiffusionStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                 ignoreJacobian,calcRule,doNotInterpolateD)
        !! initialize diffusion stencil template based on direct inputs

        type(stenciltemplate)           ,intent(inout) :: stencilTemplateObj
        type(environmentwrapper)        ,intent(inout) :: envObj
        character(*)                    ,intent(in)    :: evolvedVar
        character(*)                    ,intent(in)    :: implicitVar
        logical                         ,intent(in)    :: ignoreJacobian
        type(calculationrule) ,optional ,intent(in)    :: calcRule
        logical ,optional               ,intent(in)    :: doNotInterpolateD
    
    end subroutine initDiffusionStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCustomFluid1DStencil(stencilTemplateObj,envObj,jsonPrefix)
        !! Initialize custom fluid 1D stencil template based on environment object and JSON file

        type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)   ,intent(inout) :: envObj
        character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template

    end subroutine initCustomFluid1DStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCustomFluid1DStencilDirect(stencilTemplateObj,envObj&
                                                    ,xStencil,columnVecs,varContColVarNames,mbColVarNames)
        !! Initialize custom fluid 1D stencil template based on direct inputs
    
        type(StencilTemplate)           ,intent(inout) :: stencilTemplateObj
        type(EnvironmentWrapper)        ,intent(inout) :: envObj
        integer(ik) ,dimension(:)       ,intent(in)    :: xStencil 
        type(RealArray) ,dimension(:)   ,intent(in)    :: columnVecs
        type(StringArray) ,dimension(:) ,intent(in)    :: varContColVarNames 
        type(StringArray) ,dimension(:) ,intent(in)    :: mbColVarNames 
    
    end subroutine initCustomFluid1DStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module fluid_stencil_templates
!-----------------------------------------------------------------------------------------------------------------------------------
