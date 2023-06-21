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
module multiplicative_stencil_generator_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for multiplicative stencil

    use data_kinds                               ,only: rk, ik
    use runtime_constants                        ,only: debugging, assertions
    use god_objects                              ,only: Object
    use assertion_utility                        ,only: assert, assertIdentical, assertPure
    use support_types                            ,only: RealArray
    use variable_container_class                 ,only: VariableContainer
    use jagged_array_generator_class             ,only: JaggedArrayGenerator
    use multiplicative_generator_core_class      ,only: MultiplicativeGeneratorCore
    use modelbound_data_abstract_class           ,only: ModelboundData
    use model_surrogate_class                    ,only: ModelSurrogate

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: MultiplicativeStencilGen
        !! Generates stencil values based on independent data for x,h,v dimensions

        type(MultiplicativeGeneratorCore) ,allocatable ,private :: core !!Calculation core generated using stencil and rowCoord data

        class(JaggedArrayGenerator) ,allocatable ,private :: xStencilGen !! Optional x stencil generator
        class(JaggedArrayGenerator) ,allocatable ,private :: hStencilGen !! Optional harmonic stencil generator
        class(JaggedArrayGenerator) ,allocatable ,private :: vStencilGen !! Optional velocity space stencil generator

        type(RealArray) ,allocatable ,dimension(:) ,private :: xVals !! Buffer for x stencil values
        type(RealArray) ,allocatable ,dimension(:) ,private :: hVals !! Buffer for h stencil values
        type(RealArray) ,allocatable ,dimension(:) ,private :: vVals !! Buffer for v stencil values

        logical ,private :: fluidCol !! True if the column variable is fluid. Defaults to false.

        contains

        procedure ,public :: init => initMultGen

        procedure ,public :: calculateInPlace => calcMultVals 

        procedure ,public :: setXGen
        procedure ,public :: setHGen
        procedure ,public :: setVGen

    end type MultiplicativeStencilGen
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initMultGen(this,coreObj,fluidCol,initXVals,initHVals,initVVals) 
        !! Multiplicative stencil value generator initialization routine
    
        class(MultiplicativeStencilGen)                ,intent(inout)  :: this
        type(MultiplicativeGeneratorCore) ,allocatable ,intent(inout)  :: coreObj !! Multiplicative core - will be deallocated after call
        logical ,optional                              ,intent(in)     :: fluidCol !! True if the column variable for this stencil is fluid. 
        type(RealArray) ,optional ,dimension(:)        ,intent(in)     :: initXVals !! Optional initial raw x stencil values. Defaults to unallocated. 
        type(RealArray) ,optional ,dimension(:)        ,intent(in)     :: initHVals !! Optional initial raw h stencil values. Defaults to unallocated. 
        type(RealArray) ,optional ,dimension(:)        ,intent(in)     :: initVVals !! Optional initial raw v stencil values. Defaults to unallocated. 

    end subroutine initMultGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcMultVals (this,varCont,res,mbData,hostModel)
        !! Calculate multiplicative stencil values in place

        class(MultiplicativeStencilGen)             ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcMultVals 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setXGen (this,gen)
        !! Setter for xStencilGen.

        class(MultiplicativeStencilGen)               ,intent(inout) :: this
        class(JaggedArrayGenerator)                   ,intent(in)    :: gen

    end subroutine setXGen 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setHGen (this,gen)
        !! Setter for hStencilGen.

        class(MultiplicativeStencilGen)               ,intent(inout) :: this
        class(JaggedArrayGenerator)                   ,intent(in)    :: gen

    end subroutine setHGen 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setVGen (this,gen)
        !! Setter for vStencilGen.

        class(MultiplicativeStencilGen)               ,intent(inout) :: this
        class(JaggedArrayGenerator)                   ,intent(in)    :: gen

    end subroutine setVGen 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module multiplicative_stencil_generator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 