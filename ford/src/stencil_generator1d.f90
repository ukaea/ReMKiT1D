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
module stencil_generator1d_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for fixed stencils that exist only in one dimension

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray ,IntArray
    use variable_container_class               ,only: VariableContainer
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate
    use stencil1d_class                        ,only: Stencil1D

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: StencilGenerator1D
        !! JaggedArrayGenerator for calculating a fixed stencil based on fixed column-wise value vectors and a Stencil1D object. Each column value vector is associated with its entry in the Stencil1D object.

        type(RealArray) ,allocatable ,dimension(:) ,private :: columnVectors !! Column vectors corresponding to each of the stencil points. They must all be the same length corresponding to the global size of the dimension.
        logical                                    ,private :: periodicDim !! True if the dimension is periodic (used by the stencil to determine which columns are present in edge cases)

        type(IntArray) ,allocatable ,dimension(:)  ,public :: presentColumns !! Jagged integer array with the same shape as the result, containing the corresponding stencil point index used to sample from columnVectors

        integer(ik) ,dimension(2) ,public :: coordInterval !! Local coordinate interval (bounds inclusive) - Defaults to the entirety of the dimension

        contains

        procedure ,public :: init => initGenerator

        procedure ,public :: calculateInPlace => calcVals

    end type StencilGenerator1D
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initGenerator(this,stencilObj,columnVecs,periodicDim,coordInterval) 
        !! 1D fixed stencil value generator initialization routine
    
        class(StencilGenerator1D)             ,intent(inout)  :: this
        type(Stencil1D)                       ,intent(in)     :: stencilObj 
        type(RealArray) ,dimension(:)         ,intent(in)     :: columnVecs 
        logical ,optional                     ,intent(in)     :: periodicDim 
        integer(ik) ,dimension(2) ,optional   ,intent(in)     :: coordInterval

    end subroutine initGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcVals(this,varCont,res,mbData,hostModel)
        !! Calculate fixed 1D stencil values in place (does not depend on varCont,mbData, or hostModel)

        class(StencilGenerator1D)                   ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module stencil_generator1d_class
!-----------------------------------------------------------------------------------------------------------------------------------
 