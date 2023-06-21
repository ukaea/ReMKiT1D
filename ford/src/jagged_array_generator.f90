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
module jagged_array_generator_class
    !! author: Stefan Mijin
    !!
    !! Houses abstract class defining an interface for the calculation of jagged real arrays used in matrix construction based on variable data

    use god_objects                    ,only: Object
    use support_types                  ,only: RealArray
    use variable_container_class       ,only: VariableContainer
    use modelbound_data_abstract_class ,only: ModelboundData
    use model_surrogate_class          ,only: ModelSurrogate

    implicit none
    private
   

    type ,public ,extends(Object) ,abstract :: JaggedArrayGenerator
        !! JaggedArrayGenerator object for calculating jagged real arrays based on variable data

        contains

        procedure(calculateStencilValsInPlace) ,deferred :: calculateInPlace
        procedure ,public :: calculate

    end type JaggedArrayGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine calculateStencilValsInPlace(this,varCont,res,mbData,hostModel)
        !! Calculate stencil values in place (a RealArray that can be used as a multConst for a MatrixTerm)

        import :: JaggedArrayGenerator, VariableContainer ,RealArray ,ModelboundData ,ModelSurrogate

        class(JaggedArrayGenerator)                ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calculateStencilValsInPlace
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
        module function calculate(this,varCont,mbData,hostModel) result(res)
            !! Use in place version of stencil calculation to return values

            class(JaggedArrayGenerator)                ,intent(inout) :: this
            type(VariableContainer)                     ,intent(in)    :: varCont
            type(RealArray) ,allocatable ,dimension(:)                 :: res
            class(ModelboundData) ,optional             ,intent(in)    :: mbData
            class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

        end function calculate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module jagged_array_generator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 