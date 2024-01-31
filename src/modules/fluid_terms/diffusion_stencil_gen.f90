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
module diffusion_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for diffusion in x direction

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions, assertionLvl
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray
    use variable_container_class               ,only: VariableContainer ,CalculationRule
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use partition_class                        ,only: Partition
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: DiffusionStencilValGenerator
        !! JaggedArrayGenerator for calculating a diffusion stencil, optionally using a derivation object to calculate the diffusion coefficient
        !! at faces. The default diffusion coefficient is unity.

        real(rk) ,allocatable ,dimension(:) ,private :: outerJ !! Outer/row inverse jacobian/hodge star
        real(rk) ,allocatable ,dimension(:) ,private :: linInterp !! linear interpolation coefficients to right cell faces
        real(rk) ,allocatable ,dimension(:) ,private :: innerJ !! Inner/column jacobian/hodge star (right cell faces)

        logical ,private :: periodicGrid !! True if stencil is calculated on periodic grid

        integer(ik) ,private :: minX !! Local minX (obtained from partition data)
        integer(ik) ,private :: maxX !! Local maxX

        integer(ik) ,allocatable ,dimension(:) ,private :: diffCoeffDerivIndices !! Required variable indices for the optional diffusion coefficient derivation
        type(CalculationRule) ,allocatable ,private :: diffCoeffDerivRule !! Rule for deriving the optional diffusion coefficient
        real(rk) ,allocatable ,dimension(:) ,private :: diffCoeffBuffer
        real(rk) ,allocatable ,dimension(:) ,private :: diffCoeffBufferInterp

        logical ,private :: doNotInterpolateD !! Assume the diffusion coefficient is already calculated at cell boundaries. Defaults to false. 

        contains

        procedure ,public :: init => initDiffusionStencilGen

        procedure ,public :: calculateInPlace => calcDiffusionStencilGen

    end type DiffusionStencilValGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDiffusionStencilGen(this,partitionObj,procRank,innerJ,outerJ,linInterp&
                                             ,xPeriodic,diffCoeffDerivRule,xHaloWidth,doNotInterpolateD) 
        !! Diffusion stencil value generator initialization routine
    
        class(DiffusionStencilValGenerator)       ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        real(rk) ,dimension(:)                    ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star (should conform to x-grid)
        real(rk) ,dimension(:)                    ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star (should conform to x-grid)
        real(rk) ,dimension(:)                    ,intent(in)     :: linInterp !! Linear interpolation coefficients to right cell faces (should be size(xGrid)+1)
        logical ,optional                         ,intent(in)     :: xPeriodic !! Used to determine if outer processors should have their stencils trimmed. Defaults to .false.  
        type(CalculationRule) ,optional           ,intent(in)     :: diffCoeffDerivRule !! Rule for deriving the optional diffusion coefficient
        integer(ik) ,optional                     ,intent(in)     :: xHaloWidth !! Halo width in the x direction. Defaults to 1 and must always be >0
        logical ,optional                         ,intent(in)     :: doNotInterpolateD !! Assume the diffusion coefficient is already calcuated at cell boundaries. Defaults to False.
 
    end subroutine initDiffusionStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcDiffusionStencilGen(this,varCont,res,mbData,hostModel)
        !! Calculate diffusion stencil values in place 

        class(DiffusionStencilValGenerator)         ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcDiffusionStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module diffusion_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 