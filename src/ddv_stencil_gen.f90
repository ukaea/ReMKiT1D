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
module ddv_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for d/dv type terms 

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
   
    type ,public ,extends(JaggedArrayGenerator) :: DDVStencilGenerator
        !! JaggedArrayGenerator for calculating d/dv type stencils for single harmonics. Automatically includes spatial dependence.
        !! Has options for custom velocity grid interpolation of the distribution function and definition of the velocity
        !! cell boundary flux in the form C*f where f is the distribution function and C is a custom value. 
        !! The interpolation weights and C are assumed to be defined on right velocity cell boundaries.
        !! Assumes that the velocity stencil is [-1,0,1]

        integer(ik) ,allocatable ,dimension(:) ,private :: usedXCoords !! List of local x coordinates for which this stencil generator creates values

        real(rk)    ,allocatable ,dimension(:) ,private :: bufferC !! Buffer for velocity cell boundary flux jacobian (with or without spatial dependence). Defaults to 1. 
        real(rk)    ,allocatable ,dimension(:) ,private :: bufferInterp !! Buffer for velocity cell boundary flux jacobian (with or without spatial dependence). 
                                                                        !! Defaults to interpolation values from vSpace. 

        character(:) ,allocatable           ,private :: modelboundC !! Optional modelbound value for C (single harmonic dist)
        character(:) ,allocatable           ,private :: modelboundInterp !! Optional modelbound value for interpolation weight (single harmonic dist)

        real(rk) ,dimension(2)              ,private :: cfAtZero !! Extrapolation of C*f at zero in the form A1*f(v1)+A2*f(v2) where all A's are fixed (default = 0)

        integer(ik) ,private :: numV !! Local copy of vGrid size
        integer(ik) ,private :: locNumX !! Local x grid size
        real(rk) ,allocatable ,dimension(:) ,private :: vGridWidths !! Copy of v grid widths

        contains

        procedure ,public :: init => initDDVValGen

        procedure ,public :: calculateInPlace => calcDDVVals

    end type DDVStencilGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDDVValGen(this,partitionObj,vspaceObj,procRank,activeXCoords,fixedC,fixedInterp,mbC,mbInterp,cfAtZero) 
        !! d/dv stencil value generator initialization routine
    
        class(DDVStencilGenerator)                ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get interpolation object
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        integer(ik) ,dimension(:)                 ,intent(in)     :: activeXCoords !! List of active global x coords
        real(rk) ,optional   ,dimension(:)        ,intent(in)     :: fixedC !! Fixed C values (size numV). Defaults to ones.
        real(rk) ,optional   ,dimension(:)        ,intent(in)     :: fixedInterp !! Fixed vel interpolation values (size numV). Defaults to values from VSpace.
        character(*) ,optional                    ,intent(in)     :: mbC !! Optional modelbound value for C (single harmonic variable). Overrides fixed values.
        character(*) ,optional                    ,intent(in)     :: mbInterp !! Optional modelbound value for interpolation weight (single harmonic variable). Overrides fixed values.
        real(rk) ,optional ,dimension(2)          ,intent(in)     :: cfAtZero !! Optional extrapolation of C*f at zero in the form A1*f(v1)+A2*f(v2) where all A's are fixed. Defaults to 0.

    end subroutine initDDVValGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcDDVVals(this,varCont,res,mbData,hostModel)
        !! Calculate d/dv stencil values in place (does not depend on varCont)

        class(DDVStencilGenerator)                  ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcDDVVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module ddv_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 