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
submodule (basic_environment_wrapper) basic_environment_wrapper_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the basic environment wrapper

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initEnvironmentFromJSON(this,textbookObj,normObj)
    !! Initialize environment using data from config.json and textbook object for initializing external variables. Also initializes
    !! the status_printing module. If the textbook is not supplied the envirnment is not considered defined and finishInit must be 
    !! called before the environment is passed around.

    class(EnvironmentWrapper)  ,intent(inout) :: this
    type(Textbook) ,optional   ,intent(in)    :: textbookObj !! Textbook object used to retrieve derivation rules for variable container 
    class(Normalization) ,optional ,intent(in):: normObj !! Normalization object used to retrieve grid normalization. If not present grid is assumed normalized

    type(NamedInteger) ,dimension(3) :: intParams
    type(NamedLogical) ,dimension(2) :: logicalParams
    type(NamedLogical) ,dimension(1) :: lengthInMeters 

    type(VariableList) :: implicitVars
    type(VariableList) :: derivedVars

    real(rk) :: gridNorm

    if (assertions) then 

        call assert(.not. this%isDefined(),"initEnvironmentFromJSON called on already defined environment")

    end if

    call this%mpiCont%init()
    if (.not. associated(this%jsonCont)) then
        allocate(this%jsonCont)
        call this%jsonCont%loadFile(this%mpiCont)
    end if
    !Make sure that the original config file is preserved in debug
    if (debugging) then
        call this%jsonCont%closeFile(this%mpiCont,saveFile=.true.,filepath="./starting_config.json")
        call this%jsonCont%loadFile(this%mpiCont)
    end if
    logicalParams(1) = NamedLogical(keyXGrid//"."//keyPeriodic,.false.)

    logicalParams(2) = NamedLogical(keyPETSc//"."//keyActive,.true.)
    
    intParams(1) = NamedInteger(keyMPI//"."//keyXHaloWidth,1)

    intParams(2) = NamedInteger(keyMPI//"."//keyNumPX,this%mpiCont%getWorldSize())
    intParams(3) = NamedInteger(keyMPI//"."//keyNumPH,1)

    lengthInMeters(1) = NamedLogical(keyXGrid//"."//keyLengthInMeters,.false.)

    call this%jsonCont%load(logicalParams)
    call this%jsonCont%load(lengthInMeters)
    call this%jsonCont%load(intParams)

    call this%jsonCont%output(intParams)
    call this%jsonCont%output(logicalParams)
    call this%jsonCont%output(lengthInMeters)
    call this%mpiCont%setUpRows(intParams(2)%value,intParams(3)%value)
    call this%mpiCont%initializeNeighbourPairs(logicalParams(1)%value)

    if (present(normObj)) then 
        gridNorm = normObj%getNormalizationValue(keyLengthNorm)
        call initGridFromJSON(this%gridObj,this%jsonCont,this%mpiCont,gridNorm)
    else
        call assert(.not. lengthInMeters(1)%value,lengthInMeters(1)%name//" set to true when no normalization &
        &object supplied to initEnvironmentFromJSON")
        call initGridFromJSON(this%gridObj,this%jsonCont,this%mpiCont)
    end if

    call initPartFromJSON(this%partitionObj,this%gridObj,this%jsonCont,this%mpiCont)

    call this%mpiCont%calculateRowDistData(this%partitionObj,intParams(1)%value,this%gridObj%getNumV())

    call initGeometryFromJSON(this%geometryObj,this%gridObj,this%jsonCont,this%mpiCont)

    call initVarListFromJSON(implicitVars,this%jsonCont,this%mpiCont)
    call initVarListFromJSON(derivedVars,this%jsonCont,this%mpiCont,.true.)

    this%allVars = implicitVars%combineWith(derivedVars)

    if (logicalParams(1)%value) then 
        call this%indexingObj%init(this%partitionObj,this%gridObj,implicitVars,intParams(1)%value)
    else 
        call this%indexingObj%init(this%partitionObj,this%gridObj,implicitVars)
    end if
    call this%vSpaceObj%init(this%gridObj)
    
    call this%speciesListObj%init(this%jsonCont,this%mpiCont)
    
    if (logicalParams(2)%value) then 
        allocate(this%petscCont)
        call initPETScContFromJSON(this%petscCont,this%indexingObj,this%jsonCont,this%mpiCont)
    end if
    call initStandardSignals(this%signalCollectionObj)

    call preparePrinter(this%mpiCont%getWorldRank())

    if (present(textbookObj)) call this%finishInit(textbookObj)

end subroutine initEnvironmentFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine finishInit(this,textbookObj)
    !! Finalize the initialization of the environment by setting the textbook object and performing initializations that require it
    
    class(EnvironmentWrapper)  ,intent(inout) :: this
    type(Textbook)             ,intent(in)    :: textbookObj !! Textbook object used to retrieve derivation rules for variable container 

    if (assertions) then 

        call assert(.not. this%isDefined(),"finishInit called on already defined environment")
        call assert(textbookObj%isDefined(),"Undefined textbook object passed to finishInit")

    end if

    this%textbookObj = textbookObj

    call initVarContFromJSON(this%externalVars,this%indexingObj,this%partitionObj,this%textbookObj,this%jsonCont,this%mpiCont)

    call initHDF5ContFromJSON(this%hdf5Cont,this%externalVars,this%jsonCont,this%mpiCont)

    call this%makeDefined()

end subroutine finishInit
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeEnv(this) 

    type(EnvironmentWrapper) ,intent(inout) :: this

    nullify(this%jsonCont)

end subroutine finalizeEnv 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule basic_environment_wrapper_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
