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
module standard_modeller_assembly
    !! Contains various support routines used in data initialization

    use data_kinds                          ,only: rk, ik
    use runtime_constants                   ,only: debugging, assertions, assertionLvl
    use assertion_utility                   ,only: assert, assertIdentical, assertPure
    use json_controller_class               ,only: JSONController 
    use mpi_controller_class                ,only: MPIController ,CommunicationData
    use normalization_abstract_class        ,only: Normalization
    use composite_integrator_class          ,only: CompositeIntegrator
    use basic_environment_wrapper           ,only: EnvironmentWrapper
    use modeller_class                      ,only: Modeller 
    use custom_model_builder_class    ,only: CustomModelBuilder
    use composite_manipulator_class         ,only: CompositeManipulator
    use status_printing
    use manipulator_support
    use initialization_support
    use support_types
    use key_names

    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initStandardModeller(modellerObj,envObj,normObj)
        !! Initialize standard modeller based on config file and normalization and environment objects

        type(Modeller)            ,intent(inout) :: modellerObj
        class(EnvironmentWrapper) ,intent(inout) :: envObj    
        class(Normalization)      ,intent(in)    :: normObj   

    end subroutine initStandardModeller
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module standard_modeller_assembly
!-----------------------------------------------------------------------------------------------------------------------------------