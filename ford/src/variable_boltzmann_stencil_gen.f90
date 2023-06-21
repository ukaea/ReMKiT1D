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
module variable_boltzmann_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for v stencil of Boltzmann collision operator emission/absorption terms with variable cross-sections and energies

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray
    use variable_container_class               ,only: VariableContainer
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use v_space_class                          ,only: VSpace
    use modelbound_data_abstract_class         ,only: ModelboundData
    use modelbound_CRM_data_class              ,only: ModelboundCRMData
    use sparse_row_data_class                  ,only: SparseRowData
    use model_surrogate_class                  ,only: ModelSurrogate
    use support_functions

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: VariableBoltzmannStencilGen
        !! JaggedArrayGenerator for calculating velocity component of Boltzmann collision operator emission/absorption terms with variable inelastic mappings.
        !! Expects upper/lower triangular vStencil (depending on energy sign), otherwise vStencil = [0]. 
        !! Generates W*v*sigma(v) for absorption term or - v*sigma(v) for emission.
        !! Assumes no halo in cross-section.

        real(rk)    ,allocatable ,dimension(:) ,private :: vGridCopy !! Copy of velocity grid 

        integer(ik) ,private :: transitionIndex !! Transition index associated with this stencil generator
        integer(ik) ,private :: lNum !! Harmonic index of cross-section required, valid only if emission term

        type(SparseRowData) ,allocatable ,private :: bufferW

        logical ,private :: absorptionTerm !! True if this is for an absorption term
        logical ,private :: superelasticTerm !! True if this is for a superelastic term (negative transition cost)

        contains

        procedure ,public :: init => initVarBoltzGen

        procedure ,public :: calculateInPlace => calcVarBoltzVals

    end type VariableBoltzmannStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVarBoltzGen(this,vspaceObj,transitionIndex,lNum,absorptionTerm,superelasticTerm) 
        !! Boltzmann emission/absorption stencil value generator initialization routine with variable cross-sections and weights
    
        class(VariableBoltzmannStencilGen)        ,intent(inout)  :: this
        type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get grid
        integer(ik)                               ,intent(in)     :: transitionIndex !! Index of transition whose cross-section is required
        integer(ik) ,optional                     ,intent(in)     :: lNum !! Harmonic index of cross-section required. Used only if absorption term
        logical     ,optional                     ,intent(in)     :: absorptionTerm !! True if this is an absorption term stencil generator. Defaults to false.
        logical     ,optional                     ,intent(in)     :: superelasticTerm !! True if this is for a superelastic term (negative transition cost). Defaults to false.
        
    end subroutine initVarBoltzGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcVarBoltzVals(this,varCont,res,mbData,hostModel)
        !! Calculate Boltzmann emission/absorption stencil values in place (does not depend on varCont)

        class(VariableBoltzmannStencilGen)          ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcVarBoltzVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module variable_boltzmann_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 