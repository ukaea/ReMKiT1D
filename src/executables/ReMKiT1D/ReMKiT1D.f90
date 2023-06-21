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
program ReMKiT1D 
!! author: Stefan Mijin  
!!
!! Main project executable 

    use data_kinds                 ,only: ik, rk
    use runtime_constants          ,only: debugging, assertions
    use mpi_controller_class       ,only: MPIController
    use modeller_class             ,only: Modeller  
    use petsc_controller_class     ,only: PETScController
    use json_controller_class      ,only: JSONController
    use basic_environment_wrapper  ,only: EnvironmentWrapper
    use textbook_class             ,only: Textbook
    use basic_normalization_class  ,only: BasicNormalization
    use basic_timeloop_class       ,only: Timeloop
    use standard_modeller_assembly ,only: initStandardModeller
    use initialization_support     ,only: initStandardTextbook
    use status_printing
    implicit none 

    type(JSONController)  ,target :: normalizationJSONCont !! JSON Controller used for normalization initialization 
    type(MPIController)      :: normalizationMPICont  !! MPI Controller used only for normalization initialization
    type(Textbook)           :: textbookObj           !! Textbook object used for initializing the environment
    type(EnvironmentWrapper) :: envObj                !! Environment wrapper containing most of the main computational elements 
    type(BasicNormalization) :: normObj               !! Normalization object used in term building
    type(Modeller)           :: modellerObj           !! Main container of models/terms/data
    type(Timeloop)           :: timeloopObj           !! Object controlling the main time loop

    integer(ik) :: i 
    character(len=128) :: arg 
    character(:) ,allocatable :: alternativeConfigPath

    ! Read alternative config path from command line

    do i = 1, command_argument_count()
        call get_command_argument(i,arg)

        if (arg(:18) =='-with_config_path=') alternativeConfigPath=trim(arg(19:))

    end do

    ! Prepare normalization 

    call normalizationMPICont%init()
    if (allocated(alternativeConfigPath)) call normalizationJSONCont%setAlternativeJSONPath(alternativeConfigPath)
    call normalizationJSONCont%loadFile(normalizationMPICont)
    call normObj%init(normalizationJSONCont,normalizationMPICont)
    ! First Environment initialization pass
    envObj%jsonCont => normalizationJSONCont
    call envObj%init(normObj=normObj)
    if (debugging) call printMessage("WARNING: running code with debugging enabled. Likely compiled in debug mode.")
    if (assertions) call printMessage("WARNING: running code with assertions enabled. Likely compiled in debug mode.")

    ! Textbook initialization 

    call initStandardTextbook(textbookObj,envObj%gridObj,envObj%geometryObj,envObj%partitionObj,envObj%vSpaceObj,normObj,&
                              envObj%speciesListObj,envObj%allVars,envObj%jsonCont,envObj%mpiCont)

    ! Finish environment initialization
    call envObj%finishInit(textbookObj)
    call printMessage("Environment setup complete")

    ! Set up modeller 

    call initStandardModeller(modellerObj,envObj,normObj)
    call printMessage("Modeller setup complete")

    ! Initialize timeloop

    call timeloopObj%init(envObj,normObj)

    !Close config file and save changes 
    call envObj%jsonCont%closeFile(envObj%mpiCont,saveFile=.true.)
    ! Loop 

    call timeloopObj%loop(envObj,modellerObj)
    if (debugging) call printMessage("WARNING: Code ran with debugging enabled. Likely compiled in debug mode.")
    if (assertions) call printMessage("WARNING: Code ran with assertions enabled. Likely compiled in debug mode.")

    ! Clean up PETSc 
    call normalizationJSONCont%closeFile(normalizationMPICont)

    call envObj%petscCont%finalize()

end program 
!-----------------------------------------------------------------------------------------------------------------------------------