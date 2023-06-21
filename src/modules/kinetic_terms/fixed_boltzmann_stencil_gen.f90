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
module fixed_boltzmann_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for v stencil of Boltzmann collision operator emission/absorption terms with fixed cross-sections

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

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: FixedBoltzmannStencilGen
        !! JaggedArrayGenerator for calculating velocity component of Boltzmann collision operator emissionabsorption terms with fixed inelastic mappings.
        !! Expects fixed vStencil derived from the weights matrix if absorption term, otherwise vStencil = [0]. 
        !! Generates W*v*sigma(v) for absorption term or - v*sigma(v) for emission.
        !! If cross-section depends on x, it assumes no halo.

        real(rk)    ,allocatable ,dimension(:) ,private :: vGridCopy !! Copy of velocity grid 

        integer(ik) ,private :: transitionIndex !! Transition index associated with this stencil generator
        integer(ik) ,private :: fixedWIndex !! Inelastic mapping index associated with this stencil generator
        integer(ik) ,private :: lNum !! Harmonic index of cross-section required, valid only if emission term

        type(SparseRowData) ,allocatable ,private :: bufferW

        logical ,private :: absorptionTerm !! True if this is for an absorption term
        logical ,private :: dbTerm !! True if this is for a detailed balance term

        contains

        procedure ,public :: init => initBoltzGen

        procedure ,public :: calculateInPlace => calcBoltzVals

    end type FixedBoltzmannStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initBoltzGen(this,vspaceObj,transitionIndex,fixedWIndex,lNum,absorptionTerm,dbTerm) 
        !! Boltzmann emission/absorption stencil value generator initialization routine
    
        class(FixedBoltzmannStencilGen)           ,intent(inout)  :: this
        type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get grid
        integer(ik)                               ,intent(in)     :: transitionIndex !! Index of transition whose cross-section is required
        integer(ik)                               ,intent(in)     :: fixedWIndex !! Fixed inelastic mapping index associated with the transition
        integer(ik) ,optional                     ,intent(in)     :: lNum !! Harmonic index of cross-section required. Used only if absorption term
        logical     ,optional                     ,intent(in)     :: absorptionTerm !! True if this is an absorption term stencil generator. Defaults to false.
        logical     ,optional                     ,intent(in)     :: dbTerm !! True if this is a detailed balance term stencil generator. Defaults to false.
        
    end subroutine initBoltzGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcBoltzVals(this,varCont,res,mbData,hostModel)
        !! Calculate Boltzmann emission/absorption stencil values in place (does not depend on varCont)

        class(FixedBoltzmannStencilGen)             ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcBoltzVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module fixed_boltzmann_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 