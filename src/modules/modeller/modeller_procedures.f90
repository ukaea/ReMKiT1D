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
submodule (modeller_class) modeller_procedures
    !! author: Stefan Mijin 
    !!
    !! Contains module procedures associated with the modeller class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initModeller(this,numModels,initVars,mpiCont,petscCont,commData) 
    !! Modeller initialization routine

    class(Modeller)                   ,intent(inout) :: this
    integer(ik)                       ,intent(in)    :: numModels !! Number of models this modeller expects to be added
    type(VariableContainer)           ,intent(in)    :: initVars !! Initial variable container
    type(MPIController)               ,intent(in)    :: mpiCont !! Main MPIController 
    type(PETScController)   ,optional ,intent(in)    :: petscCont !! Optional PETSc controller - should be supplied if any integration/manipulation routine uses PETSc
    type(CommunicationData) ,optional ,intent(in)    :: commData !! Default MPI communication data

    if (assertions .or. assertionLvl >= 0) then
        call assert(initVars%isDefined(),"Initial variableContainer passed to modeller constructor is undefined")
        call assert(mpiCont%isDefined(),"MPI controller passed to modeller constructor is undefined")
        if (present(petscCont)) call assert(petscCont%isDEfined(),"Undefined PETSc controller passed to modeller constructor")
    end if

    this%numModelsAdded = 0 
    this%assembled = .false. 

    this%vars = initVars 
    this%MPICont = mpiCont 

    allocate(this%models(numModels))
    if (present(petscCont)) this%petscCont = petscCont 
    if (present(commData)) this%commData = commData

    call this%makeDefined()

end subroutine initModeller
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isAssembled(this) result(assembled)
    !! Return true if modeller is assembled and ready to evolve variables

    class(Modeller)  ,intent(in) :: this
    logical                      :: assembled

    if (assertions) call assertPure(this%isDefined(),"Attempted to get assembly status of undefined modeller")

    assembled = this%assembled

end function isAssembled
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine copyVarValuesFrom(this,varCont)
    !! Copy values of variables from outside variable container into this modeller's container

    class(Modeller)          ,intent(inout)  :: this
    type(VariableContainer)  ,intent(in)     :: varCont

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to copy variable values into undefined modeller")
        call assertPure(varCont%isDefined(),"Attempted to copy variable values into modeller from undefined variable container")
    end if

    this%vars%variables = varCont%variables

end subroutine copyVarValuesFrom
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine copyVarValuesTo(this,varCont)
    !! Copy values of variables to outside variable container from this modeller's container

    class(Modeller)          ,intent(inout)  :: this
    type(VariableContainer)  ,intent(inout)  :: varCont

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to copy variable values from undefined modeller")
        call assertPure(varCont%isDefined(),"Attempted to copy variable values from modeller into undefined variable container")
    end if

    varCont%variables = this%vars%variables 

end subroutine copyVarValuesTo
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setIntegrator(this,integ)
    !! Setter for integrator object

    class(Modeller)                          ,intent(inout)  :: this
    class(Manipulator)                       ,intent(in)     :: integ

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"Attempted to set integrator in undefined modeller")
        call assertPure(integ%isDefined(),"Attempted to set undefined integrator in modeller")
    end if

    allocate(this%integ,source = integ)

end subroutine setIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setManipulator(this,manip)
    !! Setter for manipulator object

    class(Modeller)                          ,intent(inout)  :: this
    class(CompositeManipulator)              ,intent(in)     :: manip

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"Attempted to set manipulator in undefined modeller")
        call assertPure(manip%isDefined(),"Attempted to set undefined manipulator in modeller")
    end if

    allocate(this%manip,source = manip)

end subroutine setManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateModelTermGroup(this,modelIndex,groupIndex,varCont)
    !! Call the update routine of model with given index for the given term group - optionally use variable container other than the
    !! one stored in the modeller

    class(Modeller)                   ,intent(inout)  :: this
    integer(ik)                       ,intent(in)     :: modelIndex
    integer(ik)                       ,intent(in)     :: groupIndex
    type(VariableContainer) ,optional ,intent(in)     :: varCont

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update model term group in undefined modeller object")
        call assert(this%isAssembled(),"Attempted to update model term group in unassembled modeller object")
        if (present(varCont)) call assert(varCont%isDefined(),"Attempted to update model term group in modeller object by passing&
        & undefined variable container")
    end if

    if (present(varCont)) then 
        call this%models(modelIndex)%entry%updateTermGroup(groupIndex,varCont)
    else
        call this%models(modelIndex)%entry%updateTermGroup(groupIndex,this%vars)
    end if

end subroutine updateModelTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function evaluateModelTermGroup(this,modelIndex,groupIndex,varCont) result(res)
    !! Call the evaluateTermGroup routine on model with given index and for given term group - optionally use variable container other
    !! than the one stored in the modeller

    class(Modeller)                      ,intent(in) :: this
    integer(ik)                          ,intent(in) :: modelIndex
    integer(ik)                          ,intent(in) :: groupIndex
    type(VariableContainer)  ,optional   ,intent(in) :: varCont
    real(rk) ,allocatable ,dimension(:)              :: res

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to evaluate model term group in undefined modeller object")
        call assertPure(this%isAssembled(),"Attempted to evaluate model term group in unassembled modeller object")
        if (present(varCont)) call assertPure(varCont%isDefined(),"Attempted to evaluate model term group in modeller object by &
        &passing undefined variable container")
    end if

    if (present(varCont)) then 
        res = this%models(modelIndex)%entry%evaluateTermGroup(groupIndex,varCont)
    else
        res = this%models(modelIndex)%entry%evaluateTermGroup(groupIndex,this%vars)
    end if

end function evaluateModelTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function evaluateModelTermByName(this,modelIndex,name,varCont) result(res)
        !! Call the evaluateTermByName routine on model with given index and for given term name - optionally use variable container other
        !! than the one stored in the modeller

        class(Modeller)                      ,intent(in) :: this
        integer(ik)                          ,intent(in) :: modelIndex
        character(*)                         ,intent(in) :: name
        type(VariableContainer)  ,optional   ,intent(in) :: varCont
        real(rk) ,allocatable ,dimension(:)              :: res

        if (assertions) then 
            call assertPure(this%isDefined(),"Attempted to evaluate model term in undefined modeller object")
            call assertPure(this%isAssembled(),"Attempted to evaluate model term in unassembled modeller object")
            if (present(varCont)) call assertPure(varCont%isDefined(),"Attempted to evaluate model term in modeller object by &
            &passing undefined variable container")
        end if
    
        if (present(varCont)) then 
            res = this%models(modelIndex)%entry%evaluateTermByName(name,varCont)
        else
            res = this%models(modelIndex)%entry%evaluateTermByName(name,this%vars)
        end if

    end function evaluateModelTermByName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine calculateMatGroupValsInModel(this,modelIndex,groupIndex,varCont)
    !! Calculate matrix value in implicit term group given by groupIndex in model given by modelIndex, and optionally using variable
    !! container other than the one stored in the modeller

    class(Modeller)                   ,intent(inout)  :: this
    integer(ik)                       ,intent(in)     :: modelIndex
    integer(ik)                       ,intent(in)     :: groupIndex
    type(VariableContainer) ,optional ,intent(in)     :: varCont

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to calculate mat terms in undefined modeller object")
        call assertPure(this%isAssembled(),"Attempted to calculate mat terms in unassembled modeller object")
        if (present(varCont)) call assertPure(varCont%isDefined(),"Attempted to calculate mat terms in modeller object by passing&
        & undefined variable container")
    end if

    if (present(varCont)) then 
        call this%models(modelIndex)%entry%calculateMatGroupValues(groupIndex,varCont)
    else
        call this%models(modelIndex)%entry%calculateMatGroupValues(groupIndex,this%vars)
    end if

end subroutine calculateMatGroupValsInModel
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addModelMatGroupToPETSc(this,modelIndex,groupIndex,mult,petscGroup)
    !! Send off matrix values of given term group and model to the PETSc controller, multiplied by mult

    class(Modeller)         ,intent(inout)    :: this
    integer(ik)             ,intent(in)       :: modelIndex
    integer(ik)             ,intent(in)       :: groupIndex
    real(rk)                ,intent(in)       :: mult
    integer(ik) ,optional   ,intent(in)       :: petscGroup

    integer(ik) :: usedGroup

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to add matrices to PETSc in undefined modeller object")
        call assert(this%isAssembled(),"Attempted to add matrices to PETSc in unassembled modeller object")
        call assert(allocated(this%petscCont),"Attempted to add matrices to PETSc controller when it is not allocated")
    end if

    usedGroup = 1 
    if (present(petscGroup)) usedGroup = petscGroup
    call this%models(modelIndex)%entry%addMatGroupValsToPETSc(groupIndex,this%petscCont,mult,usedGroup)

end subroutine addModelMatGroupToPETSc
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addModel(this,newModel)
    !! Add a Model to this modeller

    class(Modeller)                  ,intent(inout)  :: this
    class(Model) ,allocatable        ,intent(inout)  :: newModel

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"Attempted to add a model to undefined modeller")
        call assertPure(newModel%isDefined(),"Attempted to add undefined model to modeller")
        call assertPure(this%numModelsAdded < size(this%models),"Attempted to add model to modeller with no free model slots")
    end if

    this%numModelsAdded = this%numModelsAdded + 1
    call move_alloc(newModel,this%models(this%numModelsAdded)%entry)

end subroutine addModel
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine assemble(this,withIdentityMat)
    !! Assemble all of the matrix terms, preallocate PETScController objects and make sure modeller is ready for use. 
    !! If withIdentityMat is true preallocates diagonal elements for the identity matrix.

    class(Modeller)          ,intent(inout) :: this
    logical ,optional        ,intent(in)    :: withIdentityMat

    integer(ik) :: i

    logical :: addIdentityMat 

    addIdentityMat = .false. 
    if (present(withIdentityMat)) addIdentityMat = withIdentityMat

    if (assertions .or. assertionLvl >= 0) then 
        call assert(this%isDefined(),"Attempted to assemble undefined modeller")
        call assert(this%numModelsAdded == size(this%models),"Attempted to assemble modeller object before all models added")
        if (addIdentityMat) call assert(allocated(this%identityMat),"linearSolvePETSc called with addIdentityMat before identity&
        & matrix calculated for the modeller")
    end if

    call printMessage("Assembling models")

    if (allocated(this%petscCont)) then 
        do i = 1, size(this%models)
            call printNamedValue("Assembling model with index:",i)
            call this%models(i)%entry%assemble(this%petscCont)
        end do
        if (.not. this%petscCont%objectsCreated()) then 
            if (addIdentityMat) call this%petscCont%addRowDataToPreallocation(this%identityMat)
    
            call this%petscCont%createPETScObjs()
    
        end if
    else 
        do i = 1, size(this%models)
            call this%models(i)%entry%assemble()
        end do
    end if

    this%assembled = .true.

end subroutine assemble
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine performComm(this,commData,varCont,onlyDepth)
    !! Perform communications using the modeller's MPI controller - optionally uses external CommunicationData instructions or 
    !! exchanges data in external variable container

    class(Modeller)                   ,intent(inout) :: this
    type(CommunicationData) ,optional ,intent(in)    :: commData !! Optional non-default communication data 
    type(VariableContainer) ,optional ,intent(inout) :: varCont !! Optional variable container to perform communications on instead of the modeller's 
    integer(ik)             ,optional ,intent(in)    :: onlyDepth !! Only communicate variables at given derivation depth

    integer(ik) :: i
    type(CommunicationData) :: usedData
    logical :: varIsDist ,communicate

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to perform MPI communications from undefined modeller")
        if (present(varCont)) call assert(varCont%isDefined(),"Attempted to perform MPI communications by passing undefined &
        &variable container")
        call assert(present(commData) .or. allocated(this%commData),"performComm called without either a passed&
        &commData or allocated modeller commData component")
    end if
    if (allocated(this%commData)) usedData = this%commData
    if (present(commData)) usedData = commData 

    do i = 1,size(usedData%varsToBroadcast)

        communicate = .true. 

        if (present(onlyDepth)) communicate = onlyDepth == this%vars%getVarDepth(usedData%varsToBroadcast(i)%string)
        if (communicate) then
            varIsDist = this%vars%isVarDist(this%vars%getVarIndex(usedData%varsToBroadcast(i)%string))
            if (present(varCont)) then
                if (varIsDist) then 
                    call this%MPICont%exchangeDistVarInRow(varCont,usedData%varsToBroadcast(i)%string)
                else
                    call this%MPICont%broadcastVarInRow(varCont,usedData%varsToBroadcast(i)%string)
                end if
            else
                if (varIsDist) then 
                    call this%MPICont%exchangeDistVarInRow(this%vars,usedData%varsToBroadcast(i)%string)
                else
                    call this%MPICont%broadcastVarInRow(this%vars,usedData%varsToBroadcast(i)%string)
                end if
            end if
        end if
    end do

    do i = 1,size(usedData%haloExchangeVars)
        communicate = .true. 
        if (present(onlyDepth)) communicate = onlyDepth == this%vars%getVarDepth(usedData%haloExchangeVars(i)%string)
        if (communicate) then
            varIsDist = this%vars%isVarDist(this%vars%getVarIndex(usedData%haloExchangeVars(i)%string))
            if (present(varCont)) then
                call this%MPICont%exchangeVarXHalos(varCont,usedData%haloExchangeVars(i)%string,varIsDist)
            else
                call this%MPICont%exchangeVarXHalos(this%vars,usedData%haloExchangeVars(i)%string,varIsDist)
            end if
        end if
    end do

    do i = 1,size(usedData%scalarsToBroadcast)
        communicate = .true. 
        if (present(onlyDepth)) communicate = onlyDepth == this%vars%getVarDepth(usedData%scalarsToBroadcast(i)%string)
        if (communicate) then
            if (present(varCont)) then
                call this%MPICont%broadcastReal(varCont%variables(&
                                                varCont%getVarIndex(usedData%scalarsToBroadcast(i)%string))%entry,&
                                                usedData%scalarRoots(i))
            else
                call this%MPICont%broadcastReal(this%vars%variables(&
                                                this%vars%getVarIndex(usedData%scalarsToBroadcast(i)%string))%entry&
                                                ,usedData%scalarRoots(i))
            end if
        end if
    end do

end subroutine performComm
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine safeCommAndDeriv(this,commData,varCont,derivPriority)
    !! Perform communications interlaced with derivation calls of increasing derivation depth

    class(Modeller)                   ,intent(inout) :: this
    type(CommunicationData) ,optional ,intent(in)    :: commData !! Optional non-default communication data 
    type(VariableContainer) ,optional ,intent(inout) :: varCont !! Optional variable container to perform communications on instead of the modeller's 
    integer(ik)             ,optional ,intent(in)    :: derivPriority !! Derivation priority for interlaced calls

    integer(ik) :: i ,maxDepth
    type(CommunicationData) :: usedData

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to perform MPI communications from undefined modeller")
        if (present(varCont)) call assert(varCont%isDefined(),"Attempted to perform MPI communications by passing undefined &
        &variable container")
    end if

    if (.not. (allocated(this%commData) .or. present(commData))) then 
        if (present(varCont)) then 
            if (present(derivPriority)) then 
                call varCont%calculateDerivedVars(derivPriority=derivPriority)
            else
                call varCont%calculateDerivedVars()
            end if
        else
            if (present(derivPriority)) then 
                call this%vars%calculateDerivedVars(derivPriority=derivPriority)
            else
                call this%vars%calculateDerivedVars()
            end if
        end if
    else

        if (allocated(this%commData)) usedData = this%commData
        if (present(commData)) usedData = commData 

        maxDepth = this%vars%getMaxDepth()

        !Communicate only implicit variables
        if (present(varCont)) then 
            call this%performComm(usedData,varCont,onlyDepth=-1)
        else
            call this%performComm(usedData,this%vars,onlyDepth=-1)
        end if

        do i = 0, maxDepth
            if (present(varCont)) then 
                if (present(derivPriority)) then 
                    call varCont%calculateDerivedVars(derivPriority=derivPriority,derivDepth=i)
                else
                    call varCont%calculateDerivedVars(derivDepth=i)
                end if

                call this%performComm(usedData,varCont,onlyDepth=i)
            else
                if (present(derivPriority)) then 
                    call this%vars%calculateDerivedVars(derivPriority=derivPriority,derivDepth=i)
                else
                    call this%vars%calculateDerivedVars(derivDepth=i)
                end if

                call this%performComm(usedData,this%vars,onlyDepth=i)
            end if
        end do
    end if

end subroutine safeCommAndDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine callManipulator(this,priority,inVars,outVars)
    !! Call optional CompositeManipulator manipulate routine with inVars and outVars - the default for the optional VariableContainers is the modeller's VariableContainer

    class(Modeller)                   ,intent(inout) :: this
    integer(ik)                       ,intent(in)    :: priority
    type(VariableContainer) ,optional ,intent(in)    :: inVars
    type(VariableContainer) ,optional ,intent(inout) :: outVars

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to call manipulater from undefined modeller")
        if (present(inVars)) call assert(inVars%isDefined(),"inVars passed to callManipulator not defined")
        if (present(outVars)) call assert(outVars%isDefined(),"outVars passed to callManipulator not defined")
    end if

    if (allocated(this%manip)) then
        if (present(inVars)) then 
            if (present(outVars)) then 
                call this%manip%manipulate(this,outVars,inVars,priority)
            else
                call this%manip%manipulate(this,this%Vars,inVars,priority)
            end if
        else
            if (present(outVars)) then 
                call this%manip%manipulate(this,outVars,this%vars,priority)
            else
                call this%manip%manipulate(this,this%vars,this%vars,priority)
            end if
        end if
    end if

end subroutine callManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine integrate(this,inVars,outVars)
    !! Call integrator affect routine with inVars and outVars - the default for the optional VariableContainers is the modeller's VariableContainer

    class(Modeller)                   ,intent(inout) :: this
    type(VariableContainer) ,optional ,intent(in)    :: inVars
    type(VariableContainer) ,optional ,intent(inout) :: outVars

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to call integrator from undefined modeller")
        call assert(allocated(this%integ),"Attempted to call integrator from modeller when no integrator allocated")
        if (present(inVars)) call assert(inVars%isDefined(),"inVars passed to callIntegrator not defined")
        if (present(outVars)) call assert(outVars%isDefined(),"outVars passed to callIntegrator not defined")
    end if

    if (present(inVars)) then 
        if (present(outVars)) then 
            call this%integ%affect(this,outVars,inVars)
        else
            call this%integ%affect(this,this%Vars,inVars)
        end if
    else
        if (present(outVars)) then 
            call this%integ%affect(this,outVars,this%vars)
        else
            call this%integ%affect(this,this%vars,this%vars)
        end if
    end if

end subroutine integrate
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isModelTermGroupMixed(this,groupIndex,modelIndex) result(mixed)
!! Return true if group with given index is mixed in model with given index

    class(Modeller)  ,intent(in) :: this
    integer(ik)      ,intent(in) :: groupIndex
    integer(ik)      ,intent(in) :: modelIndex
    logical                      :: mixed

    if (assertions) then 
        call assertPure(this%isDefined(),"Called isModelTermGroupMixed on undefined modeller")
        call assertPure(this%isAssembled(),"Called isModelTermGroupMixed on unassembled modeller")
    end if

    mixed = this%models(modelIndex)%entry%isGroupMixed(groupIndex)

end function isModelTermGroupMixed
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getEvolvedVarInTermGroup(this,groupIndex,modelIndex) result(name)
    !! Return variable name for given group index and model index

    class(Modeller)           ,intent(in) :: this
    integer(ik)               ,intent(in) :: groupIndex
    integer(ik)               ,intent(in) :: modelIndex
    character(:) ,allocatable             :: name

    if (assertions) then 
        call assertPure(this%isDefined(),"Called getEvolvedVarInTermGroup on undefined modeller")
        call assertPure(this%isAssembled(),"Called getEvolvedVarInTermGroup on unassembled modeller")
    end if

    name = this%models(modelIndex)%entry%getGroupVarName(groupIndex)

end function getEvolvedVarInTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCurrentTime(this) result(time)
    !! Return current value of the integrator time variable if the integrator object is a composite integrator or the value of the 
    !! "time" variable if it isn't

    class(Modeller)           ,intent(in) :: this
    real(rk)                              :: time

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to get current time from undefined modeller")
        select type (integ => this%integ)
        type is (compositeIntegrator)
        !Do nothing
        class default
        !Check if there is a time variable in the variable container
        call assertPure(this%vars%isVarNameRegistered("time"),"Attempted to get current time from modeller where the integrator is&
        & not a composite integrator and there is no time variable")
        end select
    end if

    select type (integ => this%integ)
        type is (compositeIntegrator)
        time = integ%getCurrentTime()
        class default
        if (this%vars%isVarNameRegistered("time")) then 
            time = this%vars%variables(this%vars%getVarIndex("time"))%entry(1)
        else
            time = 0.0
        end if
    end select

end function getCurrentTime
!-----------------------------------------------------------------------------------------------------------------------------------
module function isTrueEverywhere (this,input) result(isTrue)
    !! Return true if input is true on every processor using the MPI controller of the modeller

    class(Modeller) ,intent(inout) :: this
    logical         ,intent(inout) :: input
    logical                        :: isTrue

    if (assertions) call assert(this%isDefined(),"isTrueEverywhere called on unidentified modeller")

    isTrue = this%MPICont%isTrueEverywhere(input)

end function isTrueEverywhere
!-----------------------------------------------------------------------------------------------------------------------------------
module function getGlobalMin (this,input) result(min)
    !! Return minimum value of real input computed on all processes using the MPI controller of the modeller

    class(Modeller)        ,intent(inout) :: this
    real(rk)               ,intent(inout) :: input
    real(rk)                              :: min

    if (assertions) call assert(this%isDefined(),"getGlobalMin called on unidentified modeller")

    min = this%MPICont%allreduceMin(input)

end function getGlobalMin
!-----------------------------------------------------------------------------------------------------------------------------------
module function getGlobalMax (this,input) result(max)
    !! Return maximum value of real input computed on all processes using the MPI controller of the modeller

    class(Modeller)        ,intent(inout) :: this
    real(rk)               ,intent(inout) :: input
    real(rk)                              :: max

    if (assertions) call assert(this%isDefined(),"getGlobalMax called on unidentified modeller")

    max = this%MPICont%allreduceMax(input)

end function getGlobalMax
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine linearSolvePETSc(this,knownVec,unknownVec,addIdentityMat,convReason,petscGroup) 
    !! Call linearSolve routine on PETScController of this modeller. Creates PETSc objects if they've not yet been created. If 
    !! addIdentityMat is true adds an identity matrix to the PETSc matrix.

    class(Modeller)             ,intent(inout)  :: this
    real(rk)      ,dimension(:) ,intent(in)     :: knownVec 
    real(rk)      ,dimension(:) ,intent(out)    :: unknownVec
    logical                     ,intent(in)     :: addIdentityMat
    integer(ik)                 ,intent(inout)  :: convReason
    integer(ik)   ,optional     ,intent(in)     :: petscGroup

    integer(ik) :: usedGroup

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to call PETSc linear solve in undefined modeller object")
        call assert(this%isAssembled(),"Attempted to call PETSc linear solve in unassembled modeller object")
        call assert(allocated(this%petscCont),"Attempted to call PETSc linear solve in modeller where PETSc&
        & controller is not allocated")
        if (addIdentityMat) call assert(allocated(this%identityMat),"linearSolvePETSc called with addIdentityMat before identity&
        & matrix calculated for the modeller")
    end if

    usedGroup = 1
    if (present(petscGroup)) usedGroup = petscGroup

    if (addIdentityMat) call this%petscCont%addRowValuesToMatrix(this%identityMat,objGroup=usedGroup)

    call this%petscCont%linearSolve(knownVec,unknownVec,usedGroup)

    convReason = this%petscCont%getLastConvergedReason()

end subroutine linearSolvePETSc
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calculateIdentityMat(this,indexingObj)
    !!  Initializes the identity matrix of this modeller based on Indexing object

    class(Modeller)          ,intent(inout)  :: this
    type(Indexing)           ,intent(in)     :: indexingObj

    integer(ik) ,allocatable ,dimension(:)   :: procDoFs
    integer(ik)                              :: dofOffset, locNumDoFs

    type(StringArray) ,allocatable ,dimension(:)     :: implicitVarNames
    integer(ik)       ,allocatable ,dimension(:)     :: allVarIndices
    type(IntArray)    ,allocatable ,dimension(:)     :: colArray
    integer(ik) :: i ,j ,varIndex

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to calculate identity matrix for modeller object before modeller is defined")
        call assert(indexingObj%isDefined(),"Undefined indexing object passed to calculateIdentityMat")
    end if

    if (.not. allocated(this%identityMat)) then 
        allocate(this%identityMat)
        procDoFs = indexingObj%getProcDoF()
        locNumDoFs = procDoFs(this%MPICont%getWorldRank()+1)
        dofOffset = sum(procDoFs(1:this%MPICont%getWorldRank()))

        allocate(colArray(locNumDoFs))
        do i = 1,locNumDoFs
            colArray(i) = IntArray([i+dofOffset])
        end do
        call this%identityMat%init([(i+dofOffset, i = 1,locNumDoFs)],colArray)

        do i = 1,locNumDoFs
            this%identityMat%values(i)%entry = real(1.d00,kind=rk)
        end do

        ! Set all values corresponding to stationary variables to 0 - NOTE: might want to just not include them above
        implicitVarNames = this%vars%getImplicitVarNames()

        do i = 1,size(implicitVarNames)
            if (this%vars%isStationary(implicitVarNames(i)%string)) then
                varIndex = this%vars%getVarIndex(implicitVarNames(i)%string)
                allVarIndices = indexingObj%getAllIndicesOfVar(varIndex,this%MPICont%getWorldRank())
                
                do j = 1, size(allVarIndices)
                    this%identityMat%values(allVarIndices(j))%entry = 0
                end do
            end if
        end do
        
    end if  
    
end subroutine calculateIdentityMat
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateModelData(this,modelIndex,varCont,updatePriority)
    !! Update modelbound data of model with given index if that data is allocated

    class(Modeller)         ,intent(inout)  :: this
    integer(ik)             ,intent(in)     :: modelIndex !! Model index for model whose data is to be update
    type(VariableContainer) ,intent(in)     :: varCont !! Variable container to be used in update
    integer(ik) ,optional   ,intent(in)     :: updatePriority !! Priority for this update call 

    if (assertions) call assert(this%isDefined(),"updateModelData called on undefined modeller")

    if (present(updatePriority)) then 
        call this%models(modelIndex)%entry%updateModelData(varCont,updatePriority=updatePriority)
    else
        call this%models(modelIndex)%entry%updateModelData(varCont)
    end if

end subroutine updateModelData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine copyDataFromModel(this,modelIndex,varName,container) 
    !! Copy modelbound variable data with given name from model with given index.

    class(Modeller)                       ,intent(in)    :: this 
    integer(ik)                           ,intent(in)    :: modelIndex
    character(*)                          ,intent(in)    :: varName !! Name of data
    real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into 

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to retrieve model data entry through undefined modeller")
    end if

    call this%models(modelIndex)%entry%copyModelboundDataEntry(varName,container)

end subroutine copyDataFromModel
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule modeller_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
