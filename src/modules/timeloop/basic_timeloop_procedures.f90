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
submodule (basic_timeloop_class) basic_timeloop_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the basic timeloop class 

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStandardTimeloop(this,envObj,normObj)
    !! Initialize timeloop parameters using config file and time normalization from normObj

    class(Timeloop)           ,intent(inout) :: this
    class(EnvironmentWrapper) ,intent(inout) :: envObj    
    class(Normalization)      ,intent(in)    :: normObj  

    type(NamedString)      ,dimension(2) :: modes 
    type(NamedInteger)     ,dimension(1) :: timestepNum, saveFreq ,restartFreq, initialOutputIndex 
    type(NamedReal)        ,dimension(1) :: targetTime ,minInterval
    type(NamedLogical)     ,dimension(3) :: restartSwitches
    type(NamedLogical)     ,dimension(1) :: loadFromHDF5
    type(NamedString)      ,dimension(1) :: hdf5Filename ,hdf5Filepath
    type(NamedStringArray) ,dimension(1) :: hdf5InputVars
    type(NamedRealArray)   ,dimension(1) :: outputDrivenPoints

    real(rk) :: timeNorm
    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(envObj%isDefined(),"Undefined environment wrapper passed to initStandardTimeloop")
        call assert(normObj%isDefined(),"Undefined normalization object passed to initStandardTimeloop")
    end if

    modes(1) = NamedString(keyTimeloop//"."//keyMode,keyFixedSteps)
    modes(2) = NamedString(keyTimeloop//"."//keyOutputMode,keyFixedSteps)

    timestepNum(1) = NamedInteger(keyTimeloop//"."//keyNumTimestep,1)
    initialOutputIndex(1) = NamedInteger(keyTimeloop//"."//keyRestart//"."//keyInitialOutputIndex,0)
    saveFreq(1) = NamedInteger(keyTimeloop//"."//keySaveInterval,1)
    restartFreq(1) = NamedInteger(keyTimeloop//"."//keyRestart//"."//keyFrequency,1)

    targetTime(1) = NamedReal(keyTimeloop//"."//keyTimeTarget,real(1,kind=rk))
    minInterval(1) = NamedReal(keyTimeloop//"."//keyMinSaveInterval,real(1,kind=rk))

    restartSwitches(1) = NamedLogical(keyTimeloop//"."//keyRestart//"."//keySave,.false.)
    restartSwitches(2) = NamedLogical(keyTimeloop//"."//keyRestart//"."//keyLoad,.false.)
    restartSwitches(3) = NamedLogical(keyTimeloop//"."//keyRestart//"."//keyResetTime,.false.)

    loadFromHDF5(1) = NamedLogical(keyTimeloop//"."//keyHDF5LoadInit,.false.)

    hdf5Filename(1) = NamedString(keyTimeloop//"."//keyInitPath,"ReMKiT1DVarInput")
    hdf5Filepath(1) = NamedString(keyHDF5//"."//keyFilepath,"./")

    hdf5InputVars(1)%name = keyHDF5//"."//keyInputVars
    hdf5InputVars(1)%values = envObj%externalVars%getAllVarNames()

    outputDrivenPoints(1)%name = keyTimeloop//"."//keyOutputPoints
    allocate(outputDrivenPoints(1)%values(0))

    call envObj%jsonCont%load(modes)

    timeNorm = normObj%getNormalizationValue(keyTimeNorm)

    select case(modes(1)%value)
    case (keyFixedSteps)
        this%modeTimeloop = 0
        call envObj%jsonCont%load(timestepNum)
        call envObj%jsonCont%output(timestepNum)

        this%numTimesteps = timestepNum(1)%value

    case (keyNormTimeTarget)
        this%modeTimeloop = 1
        call envObj%jsonCont%load(targetTime)
        call envObj%jsonCont%output(targetTime)

        this%targetTime = targetTime(1)%value
    case (keyRealTimeTarget)

        this%modeTimeloop = 1
        call envObj%jsonCont%load(targetTime)
        call envObj%jsonCont%output(targetTime)

        this%targetTime = targetTime(1)%value/timeNorm

    case (keyOutputDrivenMode)

        this%modeTimeloop = 2 
        call envObj%jsonCont%load(outputDrivenPoints)
        call envObj%jsonCont%output(outputDrivenPoints)

        if (assertions .or. assertionLvl >= 0) then 
            call assert(size(outputDrivenPoints(1)%values)>0,&
                "Output points must be supplied for output-driven timelooop mode")
            call assert(outputDrivenPoints(1)%values(1)>0,&
                "Output points must all be strictly positive for output-driven timeloop mode")
            do i=2,size(outputDrivenPoints(1)%values)
               call assert(outputDrivenPoints(1)%values(i)>outputDrivenPoints(1)%values(i-1),&
                   "Output points must be monotonic for output-driven timeloop mode") 
           end do
        end if

        this%outputPoints = outputDrivenPoints(1)%values

    case default 
        error stop "Unrecognized timeloop mode detected in initStandardTimeloop"
    end select

    call envObj%jsonCont%load(initialOutputIndex)
    call envObj%jsonCont%output(initialOutputIndex)
    
    this%initialOutputIndex = initialOutputIndex(1)%value

    select case(modes(2)%value)
    case (keyFixedSteps)
        this%modeSave = 0
        call envObj%jsonCont%load(saveFreq)
        call envObj%jsonCont%output(saveFreq)

        this%saveSteps = saveFreq(1)%value

    case (keyMinSaveInterval)
        this%modeSave = 1
        call envObj%jsonCont%load(minInterval)
        call envObj%jsonCont%output(minInterval)

        this%minSaveInterval = minInterval(1)%value
    case default 
        error stop "Unrecognized save mode detected in initStandardTimeloop"
    end select

    call envObj%jsonCont%load(restartSwitches)
    call envObj%jsonCont%load(loadFromHDF5)
    call envObj%jsonCont%output(restartSwitches)
    call envObj%jsonCont%output(loadFromHDF5)

    this%saveRestart = restartSwitches(1)%value
    this%loadRestart = restartSwitches(2)%value
    this%resetTimeRestart = restartSwitches(3)%value
    this%loadSerial = loadFromHDF5(1)%value 

    if (this%loadSerial) call assert(.not. this%loadRestart,&
        "Timeloop initialization found both serial and restart initialization")

    if (this%loadSerial) then

        call envObj%jsonCont%load(hdf5Filename)
        call envObj%jsonCont%load(hdf5Filepath)
        call envObj%jsonCont%load(hdf5InputVars)
        call envObj%jsonCont%output(hdf5Filename)
        call envObj%jsonCont%output(hdf5InputVars)
        allocate(this%inputHDF5Controller)
        call this%inputHDF5Controller%init(envObj%externalVars,hdf5InputVars(1)%values,hdf5Filepath(1)%value)
        this%loadFilename = hdf5Filename(1)%value 

    end if

    if (this%saveRestart) then 

        call envObj%jsonCont%load(restartFreq)
        call envObj%jsonCont%output(restartFreq)

        this%restartFrequency = restartFreq(1)%value 

    end if

    allocate(this%bufferVars,source=envObj%externalVars)

    call this%makeDefined()

end subroutine initStandardTimeloop  
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loop(this,envObj,modellerObj)
    !! Perform loop based on timeloop parameters and using passed modeller and environment

    class(Timeloop)           ,intent(inout) :: this
    class(EnvironmentWrapper) ,intent(inout) :: envObj    
    type(Modeller)            ,intent(inout) :: modellerObj

    integer(ik) :: timestepIndex ,timestepsSinceLastOutput ,timestepsSinceLastRestartDump ,outputIndex
    real(rk)    :: currentTime ,timeElapsedSinceLastOutput ,previousTime

    character(len = 140)                      :: tmpstring

    logical :: outputVals ,dumpRestart ,endOfLoopReached

    integer(ik) :: currentOutputPoint
    real(rk)    :: requestedStep


    if (assertions) then 
        call assert(this%isDefined(),"loop called on undefined Timeloop object")
        call assert(envObj%isDefined(),"Undefined environment wrapper passed to basic timeloop's loop routine")
        call assert(modellerObj%isDefined(),"Undefined modeller passed to basic timeloop's loop routine")
    end if

    call printMessage("Starting main timeloop routine")
    if (this%loadRestart) then 
        call printMessage("Loading restart files")
        call envObj%hdf5Cont%loadRestartFiles(envObj%mpiCont,this%bufferVars)

        if (this%resetTimeRestart) then 
            if (this%bufferVars%isVarNameRegistered("time")) then 
                this%bufferVars%variables(this%bufferVars%getVarIndex("time"))%entry = 0
            end if
        end if
        call modellerObj%copyVarValuesFrom(this%bufferVars)
        call printMessage("Restart files loaded")
    end if

    if (this%loadSerial) then
        call printMessage("Loading initial data from HDF5 file")
        call this%inputHDF5Controller%loadVarsSerial(envObj%mpiCont,this%bufferVars,filename=this%loadFilename)
        call modellerObj%copyVarValuesFrom(this%bufferVars)
        call printMessage("Initial values loaded from HDF5 file")
    end if

    call modellerObj%safeCommAndDeriv()
    call modellerObj%updateAllModelData()
    call modellerObj%callManipulator(4) 
    call printMessage("-------------------------------------------")
    call printMessage("Performing 0-length timestep")
    call printMessage("-------------------------------------------")
    call modellerObj%integrate(requestedTimestep=0d0)
    call printMessage("-------------------------------------------")
    call printMessage("0-length timestep successful")
    call printMessage("-------------------------------------------")
    call modellerObj%safeCommAndDeriv()
    call modellerObj%callManipulator(4) 
    call modellerObj%copyVarValuesTo(this%bufferVars)
    currentTime = 0
    timeElapsedSinceLastOutput = 0
    timestepIndex = 0 
    timestepsSinceLastOutput = 0
    timestepsSinceLastRestartDump = 0 
    outputIndex = this%initialOutputIndex
    if (this%bufferVars%isVarNameRegistered("time")) &
    currentTime = this%bufferVars%variables(this%bufferVars%getVarIndex("time"))%entry(1)

    if (outputIndex == 0) then 
        call printMessage("Outputting initial values and grid")
        call envObj%hdf5Cont%outputVarsSerial(envObj%mpiCont,this%bufferVars,IDNum=outputIndex)
        call envObj%hdf5Cont%outputGridDataSerial(envObj%mpiCont,envObj%gridObj)
        call printMessage("Initial values and grid written do hdf5 files")
    else
        call printMessage("Nonzero initial output index given - not outputting initial values") 
        outputIndex = outputIndex - 1
    end if
    call printMessage("Entering loop")

    endOfLoopReached = .false.

    currentOutputPoint = 1

    do 

        previousTime = currentTime
        timestepIndex = timestepIndex + 1
        timestepsSinceLastOutput = timestepsSinceLastOutput + 1
        timestepsSinceLastRestartDump = timestepsSinceLastRestartDump + 1

        write(tmpstring,'(A,I0)') 'Starting step ',timestepIndex
        call printMessage(trim(tmpstring))
        tmpstring=''
        
        if (this%modeTimeloop == 2) then 
            
            requestedStep = this%outputPoints(currentOutputPoint) - currentTime
            
            call modellerObj%integrate(requestedTimestep=requestedStep)
        else
            call modellerObj%integrate()
        end if
        call modellerObj%callManipulator(3) ! Call manipulators with priority 3 or lower
        call modellerObj%safeCommAndDeriv()

        currentTime = modellerObj%getCurrentTime()
        timeElapsedSinceLastOutput = timeElapsedSinceLastOutput + currentTime - previousTime

        if (this%modeTimeloop == 1) then 
            write(tmpstring,'(A,ES11.4,A,ES11.4)') 'Integration call complete. Current time (normalized): ',currentTime, &
                                                    '. Target time: ', this%targetTime
        else
            write(tmpstring,'(A,ES11.4)') 'Integration call complete. Current time (normalized): ',currentTime
        end if
        call printMessage(trim(tmpstring))
        tmpstring=''

        outputVals = .false. 
        dumpRestart = .false.

        select case (this%modeTimeloop)
        case (0)
            endOfLoopReached = timestepIndex >= this%numTimesteps
        case (1)
            endOfLoopReached = currentTime > this%targetTime
        case (2)
            endOfLoopReached = currentTime > this%outputPoints(size(this%outputPoints)) - 1d-12*currentTime
        end select

        select case (this%modeSave)
        case (0)
            outputVals = timestepsSinceLastOutput >= this%saveSteps
        case (1)
            outputVals = timeElapsedSinceLastOutput > this%minSaveInterval
        end select

        if (this%modeTimeloop == 2) then 
            !Using separate counters for genera output and this to avoid potential out-of-bounds issues
            outputVals = currentTime > this%outputPoints(currentOutputPoint) - 1d-12*currentTime
            if (outputVals) currentOutputPoint = currentOutputPoint + 1
        end if
        if (endOfLoopReached) outputVals = .true.

        if (outputVals) then 
            call modellerObj%callManipulator(4) ! Call manipulators with priority 4 or lower
            call modellerObj%safeCommAndDeriv()
            outputIndex = outputIndex + 1
            call printNamedValue("Outputting variable data to output file with index",outputIndex)
            call modellerObj%copyVarValuesTo(this%bufferVars)
            call envObj%hdf5Cont%outputVarsSerial(envObj%mpiCont,this%bufferVars,IDNum = outputIndex)
            timestepsSinceLastOutput = 0
            timeElapsedSinceLastOutput = 0
        end if

        if (this%saveRestart) then 
            dumpRestart = timestepsSinceLastRestartDump >= this%restartFrequency
            if (endOfLoopReached) dumpRestart = .true.
        end if

        if (dumpRestart) then 
            call printMessage("Dumping restart data")
            call envObj%hdf5Cont%dumpRestartFiles(envObj%mpiCont,this%bufferVars)
            timestepsSinceLastRestartDump = 0
        end if

        if (this%modeTimeloop == 0) then 
            write(tmpstring,'(A,I0,A,I0)') 'Completed step ',timestepIndex,"/",this%numTimesteps
        else
            write(tmpstring,'(A,I0)') 'Completed step ',timestepIndex
        end if
        call printMessage(trim(tmpstring))
        tmpstring=''

        if (endOfLoopReached) exit

    end do

    call printMessage("Main loop completed")

    call envObj%mpiCont%barrier()

end subroutine loop  
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule basic_timeloop_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
