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
module fluid_gen1d_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator extending StencilGenerator1D to use variable column vectors based on fluid variables

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray ,IntArray ,StringArray
    use variable_container_class               ,only: VariableContainer
    use stencil_generator1d_class              ,only: StencilGenerator1D
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate
    use stencil1d_class                        ,only: Stencil1D

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: FluidStencilGen1D
        !! JaggedArrayGenerator for extending the fixed stencil capabilities of StencilGenerator1D to allow for an additional multiplicative
        !! dependence on fluid variables in the passed variable container or modelbound data. Allows associating different variables to each column.

        type(StencilGenerator1D) ,private :: fixedStencilGen !! Fixed component generator
        type(StringArray) ,allocatable ,dimension(:) ,private :: varContColVarNames !! Names of variable container variables associated with each of the stencil columns. Max one variable ber column
        type(StringArray) ,allocatable ,dimension(:) ,private :: mbColVarNames !! Names of modelbound variables associated with each of the stencil columns. Max one variable ber column

        integer(ik) ,allocatable ,dimension(:) ,private :: varContVarIndices !! Indices of variable container variables associated with each of the stencil columns. Max one variable ber column

        type(RealArray) ,allocatable ,dimension(:) ,private :: fixedStencilVals !! Fixed stencil values calculated using parent StencilGenerator1D

        type(RealArray) ,allocatable ,dimension(:) ,private :: columnValueBuffers !! Buffer for the variable component of the column vectors

        contains

        procedure ,public :: init => initGenerator

        procedure ,public :: calculateInPlace => calcVals

    end type FluidStencilGen1D
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initGenerator(this,stencilObj,columnVecs,varContColVarNames,mbColVarNames,periodicDim,coordInterval) 
        !! 1D fluid variable stencil value generator initialization routine
    
        class(FluidStencilGen1D)              ,intent(inout)  :: this
        type(Stencil1D)                       ,intent(in)     :: stencilObj 
        type(RealArray) ,dimension(:)         ,intent(in)     :: columnVecs 
        logical ,optional                     ,intent(in)     :: periodicDim 
        integer(ik) ,dimension(2) ,optional   ,intent(in)     :: coordInterval
        type(StringArray) ,dimension(:)       ,intent(in)     :: varContColVarNames
        type(StringArray) ,dimension(:)       ,intent(in)     :: mbColVarNames

    end subroutine initGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcVals(this,varCont,res,mbData,hostModel)
        !! Calculate variable fluid 1D stencil values in place (does not depend on hostModel)

        class(FluidStencilGen1D)                   ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module fluid_gen1d_class
!-----------------------------------------------------------------------------------------------------------------------------------
 