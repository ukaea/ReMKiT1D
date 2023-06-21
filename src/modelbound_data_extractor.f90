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
module modelbound_data_extractor_class
    !! author: Stefan Mijin
    !!
    !! Houses ModelboundDataExtractor class, a Manipulator that extract data contained in a specific model's modelbound data object

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use modeller_surrogate_class              ,only: ModellerSurrogate
    use variable_container_class              ,only: VariableContainer
    use manipulator_abstract_class            ,only: Manipulator
    use modeller_class                        ,only: Modeller

    implicit none
    private

    type ,public ,extends(Manipulator) :: ModelboundDataExtractor
        !! Manipulator that extracts a (rank 1!) modelbound data variable from a model into a variable 

        integer(ik) ,private :: resultVarIndex !! Variable index where the data should be copied into (must be compatible with the data)
        integer(ik) ,private :: modelIndex !! Index of the model whose modelbound data should be copied
        character(:) ,allocatable ,private :: modelboundDataName !! Name of the (rank 1) data requested from modelbound data

        contains

        procedure ,public :: init => initExtractor
        procedure ,public :: affect => extract

    end type ModelboundDataExtractor
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initExtractor(this,resultVarIndex,modelIndex,modelboundDataName) 
        !! ModelboundDataExtractor initialization routine

        class(ModelboundDataExtractor) ,intent(inout)  :: this
        integer(ik)                    ,intent(in)     :: resultVarIndex !! Index of variable to write the result in
        integer(ik)                    ,intent(in)     :: modelIndex     !! Index of model housing required modelbound data
        character(*)                   ,intent(in)     :: modelboundDataName !! Name of data to extract

    end subroutine initExtractor
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine extract(this,manipulatedModeller,outputVars,inputVars) 
        !! Implementation of abstract manipulate routine for the evaluator

        class(ModelboundDataExtractor)        ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during manipulation
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the manipulation output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the manipulation routine

    end subroutine extract
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module modelbound_data_extractor_class
!-----------------------------------------------------------------------------------------------------------------------------------
 