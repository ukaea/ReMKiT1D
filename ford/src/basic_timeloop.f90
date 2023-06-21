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
module basic_timeloop_class
    !! author: Stefan Mijin
    !! 
    !! Houses basic timeloop class responsible for the main computation-output loop in simulations

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use modeller_class                        ,only: Modeller 
    use basic_environment_wrapper             ,only: EnvironmentWrapper
    use normalization_abstract_class          ,only: Normalization
    use variable_container_class              ,only: VariableContainer
    use HDF5_controller_class                 ,only: HDF5Controller
    use status_printing
    use support_types                         
    use key_names

    implicit none
    private

    type ,public ,extends(Object) :: Timeloop
        !! Object responsible for the main integration loop

        integer(ik) ,private :: modeTimeloop !! Loop mode: 0 - fixed step number, 1 - target elapsed time (default 0)
        integer(ik) ,private :: numTimesteps !! Fixed step number for loop mode 0
        real(rk)    ,private :: targetTime !! Target time value for loop mode 1
        integer(ik) ,private :: modeSave !! Save mode: 0 - fixed save frequency, 1 - save with minimum time interval
        integer(ik) ,private :: saveSteps !! Fixed step number for save mode 0
        real(rk)    ,private :: minSaveInterval !! Minimum time interval between saves for save mode 1

        logical     ,private :: loadRestart !! True if restart files should be loaded at the start of the loop
        logical     ,private :: loadSerial !! True if initial values are loaded from a serial file
        logical     ,private :: saveRestart !! True if restart files should be saved 
        integer(ik) ,private :: restartFrequency !! Number of steps between restart saves

        character(:)      ,allocatable       ,private :: loadFilename !! Serial load filename - default "ReMKiT1DVarInput"
        
        type(VariableContainer) ,allocatable ,private :: bufferVars

        type(HDF5Controller)    ,allocatable ,private :: inputHDF5Controller !! Controller used to input variables

        contains

        procedure ,public :: init => initStandardTimeloop
        procedure ,public :: loop

    end type Timeloop
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initStandardTimeloop(this,envObj,normObj)
        !! Initialize timeloop parameters using config file and time normalization from normObj

        class(Timeloop)           ,intent(inout) :: this
        class(EnvironmentWrapper) ,intent(inout) :: envObj    
        class(Normalization)      ,intent(in)    :: normObj  

    end subroutine initStandardTimeloop  
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loop(this,envObj,modellerObj)
        !! Perform loop based on timeloop parameters and using passed modeller and environment

        class(Timeloop)           ,intent(inout) :: this
        class(EnvironmentWrapper) ,intent(inout) :: envObj    
        type(Modeller)            ,intent(inout) :: modellerObj

    end subroutine loop  
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module basic_timeloop_class
!-----------------------------------------------------------------------------------------------------------------------------------
 