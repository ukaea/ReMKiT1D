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
module v_diffusion_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for d(Ad/dv)/dv type terms 

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
    use model_surrogate_class                  ,only: ModelSurrogate

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: VDiffStencilGen
        !! JaggedArrayGenerator for calculating d(Ad/dv)/dv type stencils. Formally a d/dv term where the flux is
        !! A*df/dv where df/dv at cell edge n+1/2 is (f_{n+1}-f_{n})/(v_{n+1}-v_{n}).
        !! A is assumed to be defined on right velocity cell boundaries.
        !! Assumes that the velocity stencil is [-1,0,1]

        integer(ik) ,allocatable ,dimension(:) ,private :: usedXCoords !! List of local x coordinates for which this stencil generator creates values

        real(rk)    ,allocatable ,dimension(:) ,private :: bufferA !! Buffer for diffusion coefficient (with or without spatial dependence). Defaults to 1. 

        character(:) ,allocatable           ,private :: modelboundA !! Optional modelbound value for A (single harmonic dist)

        real(rk) ,dimension(2)              ,private :: adfAtZero !! Extrapolation of A*df/dv at zero in the form A1*f(v1)+A2*f(v2) where all A's are fixed (default = 0)

        integer(ik) ,private :: numV !! Local copy of vGrid size
        integer(ik) ,private :: locNumX !! Local x grid size
        real(rk) ,allocatable ,dimension(:) ,private :: vGridWidths !! Copy of v grid widths
        real(rk) ,allocatable ,dimension(:) ,private:: dvPlus !! v_{n+1}-v_{n} used to calculate cell edge derivatives

        contains

        procedure ,public :: init => initVDiffValGen

        procedure ,public :: calculateInPlace => calcVDiffVals

    end type VDiffStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVDiffValGen(this,partitionObj,vspaceObj,procRank,activeXCoords,fixedA,mbA,adfAtZero) 
        !! d(Ad/dv)/dv stencil value generator initialization routine
    
        class(VDiffStencilGen)                    ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get interpolation object
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        integer(ik) ,dimension(:)                 ,intent(in)     :: activeXCoords !! List of active global x coords
        real(rk) ,optional   ,dimension(:)        ,intent(in)     :: fixedA !! Fixed A values (size numV). Defaults to ones.
        character(*) ,optional                    ,intent(in)     :: mbA !! Optional modelbound value for A (single harmonic variable). Overrides fixed values.
        real(rk) ,optional ,dimension(2)          ,intent(in)     :: adfAtZero !! Optional extrapolation of A*df/dv at zero in the form A1*f(v1)+A2*f(v2) where all A's are fixed. Defaults to 0.

    end subroutine initVDiffValGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcVDiffVals(this,varCont,res,mbData,hostModel)
        !! Calculate d(Ad/dv)/dv stencil values in place (does not depend on varCont)

        class(VDiffStencilGen)                      ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcVDiffVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module v_diffusion_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 