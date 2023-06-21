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
module manipulator_support
    !! author: Stefan Mijin
    !!
    !! Contains support for constructing manipulators based on JSON data

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use composite_manipulator_class            ,only: CompositeManipulator
    use group_evaluator_class                  ,only: GroupEvaluator
    use term_evaluator_class                   ,only: TermEvaluator
    use modelbound_data_extractor_class        ,only: ModelboundDataExtractor
    use normalization_abstract_class           ,only: Normalization
    use basic_environment_wrapper              ,only: EnvironmentWrapper
    use support_types
    use support_functions
    use key_names

    implicit none 
    private 

    public :: initCompositeManipulatorFromJSON

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCompositeManipulatorFromJSON(manip,envObj,normObj)
        !! Initialize modelbound data and add to corresponding model object

        type(CompositeManipulator) ,allocatable ,intent(inout) :: manip
        type(EnvironmentWrapper)                ,intent(inout) :: envObj
        class(Normalization)                    ,intent(in)    :: normObj
        
    end subroutine initCompositeManipulatorFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addGroupEvaluatorToCompositeManipulator(manip,envObj,normObj,jsonPrefix)
        !! Add GroupEvaluator type manipulator based json file

        type(CompositeManipulator)              ,intent(inout) :: manip
        type(EnvironmentWrapper)                ,intent(inout) :: envObj
        class(Normalization)                    ,intent(in)    :: normObj
        character(*)                            ,intent(in)    :: jsonPrefix
        
    end subroutine addGroupEvaluatorToCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addTermEvaluatorToCompositeManipulator(manip,envObj,normObj,jsonPrefix)
        !! Add TermEvaluator type manipulator based json file

        type(CompositeManipulator)              ,intent(inout) :: manip
        type(EnvironmentWrapper)                ,intent(inout) :: envObj
        class(Normalization)                    ,intent(in)    :: normObj
        character(*)                            ,intent(in)    :: jsonPrefix
        
    end subroutine addTermEvaluatorToCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addMBDataExtractorToCompositeManipulator(manip,envObj,normObj,jsonPrefix)
        !! Add ModelbounDataExtractor type manipulator based json file

        type(CompositeManipulator)              ,intent(inout) :: manip
        type(EnvironmentWrapper)                ,intent(inout) :: envObj
        class(Normalization)                    ,intent(in)    :: normObj
        character(*)                            ,intent(in)    :: jsonPrefix
        
    end subroutine addMBDataExtractorToCompositeManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module manipulator_support
!-----------------------------------------------------------------------------------------------------------------------------------