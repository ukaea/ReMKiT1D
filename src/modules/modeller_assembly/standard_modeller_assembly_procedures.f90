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
submodule (standard_modeller_assembly) standard_modeller_assembly_procedures
!! author: Stefan Mijin  
!!
!! Contains the implementation of standard modeller assembly routines

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStandardModeller(modellerObj,envObj,normObj)
    !! Initialize standard modeller based on config file and normalization and environment objects

    type(Modeller)            ,intent(inout) :: modellerObj
    class(EnvironmentWrapper) ,intent(inout) :: envObj    
    class(Normalization)      ,intent(in)    :: normObj  

    type(CommunicationData)                  :: commData

    type(NamedStringArray)     ,dimension(4) :: stringArrayParams
    type(NamedString)          ,dimension(1) :: modelType
    type(NamedIntegerArray)    ,dimension(1) :: intArrayParams
    type(CompositeIntegrator)                :: integratorObj
    type(CustomModelBuilder)                 :: builderCustom

    type(CompositeManipulator) ,allocatable :: manipObj

    logical :: hasComm

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(envObj%isDefined(),"Undefined environment wrapper passed to initStandardModeller")
        call assert(normObj%isDefined(),"Undefined normalization object passed to initStandardModeller")
    end if

    stringArrayParams(1)%name = keyModels//"."//keyTags
    allocate(stringArrayParams(1)%values(0))

    stringArrayParams(2)%name = keyMPI//"."//keyCommData//"."//keyVarsToBroadcast
    allocate(stringArrayParams(2)%values(0))

    stringArrayParams(3)%name = keyMPI//"."//keyCommData//"."//keyHaloVars
    allocate(stringArrayParams(3)%values(0))

    stringArrayParams(4)%name = keyMPI//"."//keyCommData//"."//keyScalarsToBroadcast
    allocate(stringArrayParams(4)%values(0))

    intArrayParams(1)%name = keyMPI//"."//keyCommData//"."//keyScalarsRoots
    allocate(intArrayParams(1)%values(0))

    call envObj%jsonCont%load(stringArrayParams)
    call envObj%jsonCont%load(intArrayParams)
    
    if (size(intArrayParams(1)%values) == 0 .and. size(stringArrayParams(4)%values) > 0) then 
        intArrayParams(1)%values = [(0,i=1,size(stringArrayParams(4)%values))]
    end if
    if (assertions .or. assertionLvl >= 0) then 
        call assert(size(stringArrayParams(1)%values) > 0,"No models detected by initStandardModeller")

        if (size(stringArrayParams(4)%values) > 0) then 
            call assert(size(intArrayParams(1)%values)==size(stringArrayParams(4)%values),intArrayParams(1)%name//" size must &
            &conform to size of "//stringArrayParams(4)%name)
        end if

    end if

    hasComm = size(stringArrayParams(2)%values) + size(stringArrayParams(3)%values) > 0

    if (hasComm) then
        commData%varsToBroadcast = stringArrayParams(2)%values
        commData%haloExchangeVars = stringArrayParams(3)%values
        commData%scalarsToBroadcast = stringArrayParams(4)%values
        commData%scalarRoots = intArrayParams(1)%values
    end if


    call printMessage("Initializing modeller")
    if (allocated(envObj%petscCont)) then 

        if (hasComm) then 
            call modellerObj%init(size(stringArrayParams(1)%values),envObj%externalVars,envObj%mpiCont,envObj%petscCont,commData)
        else
            call modellerObj%init(size(stringArrayParams(1)%values),envObj%externalVars,envObj%mpiCont,envObj%petscCont)
        end if 
        call printMessage("Calculating PETSc identity matrix")
        call modellerObj%calculateIdentityMat(envObj%indexingObj)

    else

        if (hasComm) then 
            call modellerObj%init(size(stringArrayParams(1)%values),envObj%externalVars,envObj%mpiCont,commData=commData)
        else
            call modellerObj%init(size(stringArrayParams(1)%values),envObj%externalVars,envObj%mpiCont)
        end if

    end if

    call printMessage("Initializing integrators")
    call initStandardIntegrator(integratorObj,envObj%externalVars,envObj%indexingObj,envObj%jsonCont,envObj%mpiCont)

    call modellerObj%setIntegrator(integratorObj)

    do i = 1,size(stringArrayParams(1)%values)
        modelType(1) = NamedString(keyModels//"."//stringArrayParams(1)%values(i)%string//"."//keyType,"")
        call envObj%jsonCont%load(modelType)
        select case (modelType(1)%value)
        case (keyCustomModel)
            call printMessage("Building model: "//stringArrayParams(1)%values(i)%string)
            call builderCustom%init(envObj,normObj,stringArrayParams(1)%values(i)%string)
            call builderCustom%addModelToModeller(modellerObj)
        case default 
            error stop "Unsupported model type detected by initStandardModeller"
        end select
    end do

    call initCompositeManipulatorFromJSON(manipObj,envObj,normObj)

    if (allocated(manipObj)) then 
        call modellerObj%setManipulator(manipObj)
    end if
    call printMessage("Assembling modeller")
    call modellerObj%assemble(withIdentityMat=allocated(envObj%petscCont))


end subroutine initStandardModeller
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule standard_modeller_assembly_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
