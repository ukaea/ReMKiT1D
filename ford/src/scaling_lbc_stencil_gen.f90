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
module scaling_lbc_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for logical boundary condition with scaled distribution extrapolation

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray
    use variable_container_class               ,only: VariableContainer
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use partition_class                        ,only: Partition
    use v_space_class                          ,only: VSpace
    use modelbound_data_abstract_class         ,only: ModelboundData
    use modelbound_lbc_data_class              ,only: ModelboundLBCData
    use f_scaling_derivation_class             ,only: FScalingDerivation
    use model_surrogate_class                  ,only: ModelSurrogate

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: ScalingLBCStencilGen
        !! JaggedArrayGenerator for calculating logical boundary condition for electron distribution function. Contains x,h,v dependencies. 
        !! Expects xStencil = [0] or [-1], depending on whether the boundary is left or right, and whether variables are staggered. 
        !! The hStencil should cover all harmonics and assumes no m>0 harmonics. Expects vStencil = [0,1].
        !! Calculates v*f_{cl}, without taking the x width or boundary sign into account.


        type(FScalingDerivation) ,private :: fDeriv !! Derivation used to get scaling factors 
        integer(ik) ,allocatable ,dimension(:) ,private :: fDerivReqVarIndices !! Required variable indices for scaling extrapolation derivation

        real(rk) ,allocatable ,dimension(:,:) ,private :: bufferPl !! Tensor buffer

        real(rk) ,allocatable ,dimension(:) ,private :: vGrid !! Copy of v grid
        real(rk) ,allocatable ,dimension(:) ,private :: vGridWidths !! Copy of v grid widths

        integer(ik) ,private :: colL !! l harmonic index of the column harmonic
        integer(ik) ,allocatable ,dimension(:) :: includedDecompHarmonics !! Harmonics to be included in the decomposition for this generator (useful for right boundary on staggered grid)

        real(rk) ,private :: dx !! Width of boundary cell
        logical ,private :: isActive

        contains

        procedure ,public :: init => initScalingLBCGen

        procedure ,public :: calculateInPlace => calcScalingLBCVals

    end type ScalingLBCStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initScalingLBCGen(this,vspaceObj,isActive,fDeriv,fDerivReqIndices,colL,dx,includedDecompHarmonics) 
        !! Scaling logical boundary condition value generator initialization routine
    
        class(ScalingLBCStencilGen)               ,intent(inout)  :: this
        type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get grid data 
        logical                                   ,intent(in)     :: isActive !! Set to true if the process containing this has a boundary where this is to be applied 
        type(FScalingDerivation)                  ,intent(in)     :: fDeriv !! Derivation used to get scaling factors 
        integer(ik) ,dimension(:)                 ,intent(in)     :: fDerivReqIndices !! Indices required for scaling extrapolation derivation
        integer(ik)                               ,intent(in)     :: colL !! Column l harmonic
        real(rk)                                  ,intent(in)     :: dx !! Boundary cell width
        integer(ik) ,optional ,dimension(:)       ,intent(in)     :: includedDecompHarmonics  !! Harmonics to be included in the decomposition for this generator (useful for right boundary on staggered grid). Defaults to all harmonics.
    
    end subroutine initScalingLBCGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcScalingLBCVals(this,varCont,res,mbData,hostModel)
        !! Calculate scaling logical boundary condition values in place 

        class(ScalingLBCStencilGen)                  ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcScalingLBCVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module scaling_lbc_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 