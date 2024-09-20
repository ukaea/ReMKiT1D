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
submodule (model_class) model_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the model class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initModel(this,numImplicitTerms,numGeneralTerms,numImplicitGroups,numGeneralGroups) 
    !! Model initialization routine

    class(Model)          ,intent(inout)  :: this
    integer(ik) ,optional ,intent(in)     :: numImplicitTerms !! Number of MatrixTerm objects this model expects to be added
    integer(ik) ,optional ,intent(in)     :: numGeneralTerms !! Number of general Term objects this model expects to be added
    integer(ik) ,optional ,intent(in)     :: numImplicitGroups !! Number of implicit/matrix term groups registered with this model 
    integer(ik) ,optional ,intent(in)     :: numGeneralGroups !! Number of general term groups registered with this model

    integer(ik) :: i

    call this%makeDefined()

    this%assembled = .false.
    this%numAddedMatrixTerms = 0
    this%numAddedGeneralTerms = 0

    this%setupCounter = .false. 

    if (present(numImplicitTerms)) call this%setNumImplicitTerms(numImplicitTerms)
    if (present(numGeneralTerms)) call this%setNumGeneralTerms(numGeneralTerms)
    if (present(numImplicitGroups)) call this%setNumImplicitGroups(numImplicitGroups)
    if (present(numGeneralGroups)) call this%setNumGeneralGroups(numGeneralGroups)

end subroutine initModel
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isAssembled(this) result(assembled)
    !! Check if model is assembled and ready to use

    class(Model)  ,intent(in) :: this
    logical                   :: assembled

    if (assertions) call assertPure(this%isDefined(),"Attempted to get assembly status from undefined model object")

    assembled = this%assembled

end function isAssembled
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNumImplicitTerms(this,numImplicitTerms) 
    !! Set number of implicit terms and perform allocation

    class(Model)             ,intent(inout)  :: this
    integer(ik)              ,intent(in)     :: numImplicitTerms !! Number of MatrixTerm objects this model expects to be added

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"setNumImplicitTerms called on undefined model object")
        call assertPure(.not. this%setupCounter(1),"setNumImplicitTerms called more than once")
    end if

    allocate(this%implicitTerms(numImplicitTerms))
    allocate(this%implicitTermNames(numImplicitTerms))
    allocate(this%skipPattern(numImplicitTerms))
    this%skipPattern = .false.
    this%setupCounter(1) = .true.

end subroutine setNumImplicitTerms
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNumGeneralTerms(this,numGeneralTerms) 
    !! Set number of general terms and perform allocation

    class(Model)             ,intent(inout)  :: this
    integer(ik)              ,intent(in)     :: numGeneralTerms !! Number of general Term objects this model expects to be 
    
    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"setNumGeneralTerms called on undefined model object")
        call assertPure(.not. this%setupCounter(2),"setNumGeneralTerms called more than once")
    end if

    allocate(this%generalTerms(numGeneralTerms))
    allocate(this%generalTermNames(numGeneralTerms))

    this%setupCounter(2) = .true.

end subroutine setNumGeneralTerms
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNumImplicitGroups(this,numImplicitGroups) 
    !! Set number of implicit groups and perform allocation

    class(Model)             ,intent(inout)  :: this
    integer(ik)              ,intent(in)     :: numImplicitGroups !! Number of implicit/matrix term groups registered with this model 

    integer(ik) :: i 

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"setNumImplicitGroups called on undefined model object")
        call assertPure(.not. this%setupCounter(3),"setNumImplicitGroups called more than once")
    end if

    allocate(this%implicitTermGroup(numImplicitGroups))
    allocate(this%implicitGroupMixed(numImplicitGroups))

    this%implicitGroupMixed = .true.

    do i = 1, numImplicitGroups
        allocate(this%implicitTermGroup(i)%entry(0))
    end do

    this%setupCounter(3) = .true.

end subroutine setNumImplicitGroups
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNumGeneralGroups(this,numGeneralGroups) 
    !! Set number of general groups and perform allocation

    class(Model)             ,intent(inout)  :: this
    integer(ik)              ,intent(in)     :: numGeneralGroups !! Number of general term groups registered with this model
    
    integer(ik) :: i 

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"setNumGeneralGroups called on undefined model object")
        call assertPure(.not. this%setupCounter(4),"setNumGeneralGroups called more than once")
    end if

    allocate(this%generalTermGroup(numGeneralGroups))

    allocate(this%generalGroupMixed(numGeneralGroups))

    this%generalGroupMixed = .true.

    do i = 1, numGeneralGroups
        allocate(this%generalTermGroup(i)%entry(0))
    end do

    this%setupCounter(4) = .true.

end subroutine setNumGeneralGroups
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addImplicitTerm(this,impTerm,implicitGroups,generalGroups,termName,skipPattern)
    !! Add a MatrixTerm object to the model (deallocating the source!), and specify which implicit and general groups it belongs to

    class(Model)                          ,intent(inout)  :: this
    class(MatrixTerm) ,allocatable        ,intent(inout)  :: impTerm !! MatrixTerm object to be reallocated to this model
    integer(ik)             ,dimension(:) ,intent(in)     :: implicitGroups !! Implicit groups the added term should belong to
    integer(ik)             ,dimension(:) ,intent(in)     :: generalGroups !! General groups the added term should belong to
    character(*)                          ,intent(in)     :: termName !! Name of added term for indexing purposes
    logical  ,optional                    ,intent(in)     :: skipPattern !! True if the matrix term pattern should not be added to PETSc preallocation

    integer(ik) :: i ,prevIndex ,prevGeneralTermIndex

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"Attempted to add implicit term to undefined model object")
        call assertPure(impTerm%isDefined(),"Attempted to add undefined implicit term to model object")
        call assertPure(all(this%setupCounter),"Attempted to add implicit term to model object before setting up allocations")
        call assertPure(this%numAddedMatrixTerms < size(this%implicitTerms),"Attempted to add implicit term to model when no&
        & free term slots available")
        call assertPure(size(implicitGroups) >= 1,&
        "At least one implicit group must be specifyied when adding implicit term to model")
    end if

    this%numAddedMatrixTerms = this%numAddedMatrixTerms + 1
    if (present(skipPattern)) this%skipPattern(this%numAddedMatrixTerms) = skipPattern
    call move_alloc(impTerm,this%implicitTerms(this%numAddedMatrixTerms)%entry)
    this%implicitTermNames(this%numAddedMatrixTerms)%string = termName
    do i = 1, size(implicitGroups)
        this%implicitTermGroup(implicitGroups(i))%entry = &
        [this%implicitTermGroup(implicitGroups(i))%entry,this%numAddedMatrixTerms]
        prevIndex = size(this%implicitTermGroup(implicitGroups(i))%entry) - 1
        !Determine whether the group is mixed based on previous term added to the group
        if (prevIndex > 0) then 
            if (.not. this%implicitGroupMixed(implicitGroups(i))) &
            this%implicitGroupMixed(implicitGroups(i)) = this%implicitTerms(this%numAddedMatrixTerms)%entry%getVarName() /= &
            this%implicitTerms(this%implicitTermGroup(implicitGroups(i))%entry(prevIndex))%entry%getVarName()
        else
            this%implicitGroupMixed(implicitGroups(i)) = .false.
        end if
    end do

    do i = 1, size(generalGroups)
        this%generalTermGroup(generalGroups(i))%entry = &
        [this%generalTermGroup(generalGroups(i))%entry,this%numAddedMatrixTerms]
        prevIndex = size(this%generalTermGroup(generalGroups(i))%entry) - 1
        !Determine whether the group is mixed based on previous term added to the group

        if (prevIndex > 0) then 
            prevGeneralTermIndex = this%generalTermGroup(generalGroups(i))%entry(prevIndex)
            if (.not. this%generalGroupMixed(generalGroups(i))) then 
                if (prevGeneralTermIndex > size(this%implicitTerms)) then
                    !Weird gfortran bug workaround (if the following is used to index the array directly when size(this%implicitTerms)=0
                    !a segfault occurs)
                    prevGeneralTermIndex = prevGeneralTermIndex - size(this%implicitTerms)
                    this%generalGroupMixed(generalGroups(i)) = this%implicitTerms(this%numAddedMatrixTerms)%entry%getVarName() /=&
                    this%generalTerms(prevGeneralTermIndex)%entry%getVarName()
                else
                    this%generalGroupMixed(generalGroups(i)) = this%implicitTerms(this%numAddedMatrixTerms)%entry%getVarName() /=&
                    this%implicitTerms(prevGeneralTermIndex)%entry%getVarName()
                end if
            end if
        else
            this%generalGroupMixed(generalGroups(i)) = .false.
        end if
    end do

   
    
end subroutine addImplicitTerm
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addGeneralTerm(this,genTerm,generalGroups,termName)
    !! Add a Term object to the model (deallocating the source!), and specify which general groups it belongs to

    class(Model)                   ,intent(inout)  :: this
    class(Term)       ,allocatable ,intent(inout)  :: genTerm !! General Term object to be reallocated to this model
    integer(ik)      ,dimension(:) ,intent(in)     :: generalGroups !! General groups the added term should belong to
    character(*)                          ,intent(in)     :: termName !! Name of added term for indexing purposes

    integer(ik) :: i ,prevIndex ,prevGeneralTermIndex

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"Attempted to add general term to undefined model object")
        call assertPure(genTerm%isDefined(),"Attempted to add undefined general term to model object")
        call assertPure(all(this%setupCounter),"Attempted to add general term to model object before setting up allocations")
        call assertPure(this%numAddedGeneralTerms < size(this%generalTerms),"Attempted to add general term to model when no&
        & free term slots available")
        call assertPure(size(generalGroups) >= 1,&
        "At least one general group must be specifyied when adding general term to model")
    end if

    this%numAddedGeneralTerms = this%numAddedGeneralTerms + 1
    call move_alloc(genTerm,this%generalTerms(this%numAddedGeneralTerms)%entry)
    this%generalTermNames(this%numAddedGeneralTerms)%string = termName

    do i = 1, size(generalGroups)
        this%generalTermGroup(generalGroups(i))%entry = &
        [this%generalTermGroup(generalGroups(i))%entry,this%numAddedGeneralTerms+size(this%implicitTerms)]
        prevIndex = size(this%generalTermGroup(generalGroups(i))%entry) - 1
        if (prevIndex > 0) then 
            if (.not. this%generalGroupMixed(generalGroups(i))) then 
                prevGeneralTermIndex = this%generalTermGroup(generalGroups(i))%entry(prevIndex)
                !Determine whether the group is mixed based on previous term added to the group
                if (.not. this%generalGroupMixed(generalGroups(i))) then 
                    if (prevGeneralTermIndex > size(this%implicitTerms)) then
                        !Weird gfortran bug workaround (if the following is used to index the array directly when size(this%implicitTerms)=0
                        !a segfault occurs)
                        prevGeneralTermIndex = prevGeneralTermIndex - size(this%implicitTerms)
                        this%generalGroupMixed(generalGroups(i)) = this%generalTerms(this%numAddedGeneralTerms)%entry%getVarName() &
                        /= this%generalTerms(prevGeneralTermIndex)%entry%getVarName()
                    else
                        this%generalGroupMixed(generalGroups(i)) = this%generalTerms(this%numAddedGeneralTerms)%entry%getVarName() &
                        /= this%implicitTerms(prevGeneralTermIndex)%entry%getVarName()
                    end if
                end if
            end if
        else
            this%generalGroupMixed(generalGroups(i)) = .false.
        end if
    end do

end subroutine addGeneralTerm
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateTermGroup(this,groupIndex,varCont)
    !! Update a term group - if groupIndex > size(implicitGroup) it is taken to be in the general group

    class(Model)            ,intent(inout)  :: this
    integer(ik)             ,intent(in)     :: groupIndex !! Group index to be updated 
    type(VariableContainer) ,intent(in)     :: varCont !! Variable container used to update this group

    integer(ik) :: i ,trueGroupIndex ,termIndex

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to update term group in undefined model object")
        call assertPure(this%isAssembled(),"Attempted to update term group in unassembled model object")
        call assertPure(varCont%isDefined(),&
        "Attempted to update term group in model object by passing undefined variable container")
    end if

    if (groupIndex > size(this%implicitTermGroup)) then 
        trueGroupIndex = groupIndex - size(this%implicitTermGroup)
        do i = 1, size(this%generalTermGroup(trueGroupIndex)%entry)
            if (this%generalTermGroup(trueGroupIndex)%entry(i) > size(this%implicitTerms)) then
                termIndex = this%generalTermGroup(trueGroupIndex)%entry(i) - size(this%implicitTerms)
                if (allocated(this%modelData)) then 
                    call this%generalTerms(termIndex)%entry%update(varCont,this%modelData,hostModel=this)
                else
                    call this%generalTerms(termIndex)%entry%update(varCont,hostModel=this)
                end if
            else 
                termIndex = this%generalTermGroup(trueGroupIndex)%entry(i)
                if (allocated(this%modelData)) then 
                    call this%implicitTerms(termIndex)%entry%update(varCont,this%modelData,hostModel=this)
                else
                    call this%implicitTerms(termIndex)%entry%update(varCont,hostModel=this)
                end if
            end if
        end do
    else 
        do i = 1, size(this%implicitTermGroup(groupIndex)%entry)
            termIndex = this%implicitTermGroup(groupIndex)%entry(i)
            if (allocated(this%modelData)) then 
                call this%implicitTerms(termIndex)%entry%update(varCont,this%modelData,hostModel=this)
            else
                call this%implicitTerms(termIndex)%entry%update(varCont,hostModel=this)
            end if
        end do
    end if
    
end subroutine updateTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
module function evaluateTermGroup(this,groupIndex,varCont) result(res)
    !! Evaluate a term group, returning the sum of all explicit results from the term group 
    !! - if groupIndex > size(implicitGroup) it is taken to be in the general group

    class(Model)                         ,intent(in) :: this
    integer(ik)                          ,intent(in) :: groupIndex !! Group index to evaluate
    type(VariableContainer)              ,intent(in) :: varCont !! Variable container used to evaluate this group
    real(rk) ,allocatable ,dimension(:)              :: res

    integer(ik) :: i ,trueGroupIndex ,termIndex

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to evaluate term group in undefined model object")
        call assertPure(this%isAssembled(),"Attempted to evaluate term group in unassembled model object")
        call assertPure(varCont%isDefined(),&
        "Attempted to evaluate term group in model object by passing undefined variable container")
        call assertPure(groupIndex > 0 ,"Negative group index passed to evaluateTermGroup")
        call assertPure(groupIndex <= size(this%implicitTermGroup) + size(this%generalTermGroup),&
        "Out of range group index passed to evaluateTermGroup")
        if (groupIndex > size(this%implicitTermGroup)) then 
            trueGroupIndex = groupIndex - size(this%implicitTermGroup)
            call assertPure(.not. this%generalGroupMixed(trueGroupIndex),"Attempted to evaluate mixed term group in model")
        else
            call assertPure(.not. this%implicitGroupMixed(groupIndex),"Attempted to evaluate mixed term group in model")
        end if
    end if

    if (groupIndex > size(this%implicitTermGroup)) then 
        trueGroupIndex = groupIndex - size(this%implicitTermGroup)
        do i = 1, size(this%generalTermGroup(trueGroupIndex)%entry)
            if (this%generalTermGroup(trueGroupIndex)%entry(i) > size(this%implicitTerms)) then
                termIndex = this%generalTermGroup(trueGroupIndex)%entry(i) - size(this%implicitTerms)
                if (.not.allocated(res)) then 
                    allocate(res,source=this%generalTerms(termIndex)%entry%evaluate(varCont))
                else
                    res = res + this%generalTerms(termIndex)%entry%evaluate(varCont)
                end if
            else 
                termIndex = this%generalTermGroup(trueGroupIndex)%entry(i)
                if (.not.allocated(res)) then 
                    allocate(res,source=this%implicitTerms(termIndex)%entry%evaluate(varCont))
                else
                    res = res + this%implicitTerms(termIndex)%entry%evaluate(varCont)
                end if
            end if
        end do
    else 
        do i = 1, size(this%implicitTermGroup(groupIndex)%entry)
            termIndex = this%implicitTermGroup(groupIndex)%entry(i)
            if (.not.allocated(res)) then 
                res=this%implicitTerms(termIndex)%entry%evaluate(varCont)
            else
                res = res + this%implicitTerms(termIndex)%entry%evaluate(varCont)
            end if
        end do
    end if

end function evaluateTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine calculateMatGroupValues(this,groupIndex,varCont)
    !! Calculate matrix value in implicit term group given by groupIndex

    class(Model)            ,intent(inout)  :: this
    integer(ik)             ,intent(in)     :: groupIndex !! Group index of terms whose matrices should be calculated
    type(VariableContainer) ,intent(in)     :: varCont !! Variable container used for the calculation

    integer(ik) :: i ,termIndex

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to calculate matrix values of term group in undefined model object")
        call assertPure(this%isAssembled(),"Attempted to calculate matrix values of term group in unassembled model object")
        call assertPure(varCont%isDefined(),&
        "Attempted to calculate matrix values of term group in model object by passing undefined variable container")
        call assertPure(groupIndex > 0 ,"Negative group index passed to calculateMatGroupValues")
        call assertPure(groupIndex <= size(this%implicitTermGroup)+size(this%generalTermGroup),&
        "Out of range group index passed to calculateMatGroupValues")
    end if

    if (groupIndex <= size(this%implicitTermGroup)) then
        do i = 1,size(this%implicitTermGroup(groupIndex)%entry)
            termIndex = this%implicitTermGroup(groupIndex)%entry(i)
            call this%implicitTerms(termIndex)%entry%calculateValues(varCont)
        end do
    end if
end subroutine calculateMatGroupValues
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addMatGroupValsToPETSc(this,groupIndex,petscCont,mult,petscGroup)
    !! Send off matrix values of given term group to the PETSc controller, multiplied by mult

    class(Model)            ,intent(in)    :: this
    integer(ik)             ,intent(in)    :: groupIndex !! Group index of terms whose matrices should be sent to PETSc
    type(PETScController)   ,intent(inout) :: petscCont !! PETScController object housing PETSc matrices 
    real(rk)                ,intent(in)    :: mult !! Multiplier used when adding matrices to PETSc - usually -dt 
    integer(ik) ,optional   ,intent(in)    :: petscGroup

    integer(ik) :: i ,termIndex ,usedGroup

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to send matrix values to PETSc from undefined model object")
        call assert(this%isAssembled(),"Attempted to send matrix values to PETSc from unassembled model object")
        call assert(petscCont%isDefined(),"Attempted to send matrix values to undefined PETSc controller from model object")
        call assertPure(groupIndex > 0 ,"Negative group index passed to addMatGroupValsToPETSc")
        call assertPure(groupIndex <= size(this%implicitTermGroup),&
        "Out of range group index passed to addMatGroupValsToPETSc")
    end if

    usedGroup = 1 
    if (present(petscGroup)) usedGroup = petscGroup
    do i = 1,size(this%implicitTermGroup(groupIndex)%entry)
        termIndex = this%implicitTermGroup(groupIndex)%entry(i)
        call this%implicitTerms(termIndex)%entry%addRowValuesToPETScMatrix(petscCont,mult,usedGroup) 
    end do

end subroutine addMatGroupValsToPETSc
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine assemble(this,petscCont)
    !! Assemble all of the matrix terms, preallocate PETScController objects and make sure model is ready for use

    class(Model)                     ,intent(inout) :: this
    type(PETScController)  ,optional ,intent(inout) :: petscCont !! Optional PETScController - should be present if the model has any implicitly evaluated terms

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(this%isDefined(),"Attempted to assemble undefined model object")
        call assertPure(all(this%setupCounter),"Attempted to assemble model object before setting up allocations")
        if (present(petscCont)) call assert(petscCont%isDefined(),"Attempted to assemble model object using undefined PETSc object")
        call assert(this%numAddedGeneralTerms == size(this%generalTerms),"Attempted to assemble model object before all general&
        & terms have been added")
        call assert(this%numAddedMatrixTerms == size(this%implicitTerms),"Attempted to assemble model object before all &
        &implicit terms have been added")
    end if

    if (present(petscCont)) then
        do i = 1, size(this%implicitTerms)
            if (.not. this%skipPattern(i)) &
                call this%implicitTerms(i)%entry%addRowDataPatternToController(petscCont) 
        end do
    end if

    this%assembled = .true.

end subroutine assemble
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isGroupMixed(this,groupIndex) result(mixed)
    !! Return true if group with given index is mixed

    class(Model)  ,intent(in) :: this
    integer(ik)   ,intent(in) :: groupIndex
    logical                   :: mixed

    if (assertions) then 
        call assertPure (this%isDefined(),"Called isGroupMixed from undefined model object")
        call assertPure(groupIndex > 0 ,"Negative group index passed to isGroupMixed")
        call assertPure(groupIndex <= size(this%implicitTermGroup)+size(this%generalTermGroup),&
        "Out of range group index passed to isGroupMixed")
    end if
    if (groupIndex > size(this%implicitTermGroup)) then
        mixed = this%generalGroupMixed(groupIndex-size(this%implicitTermGroup))
    else
        mixed = this%implicitGroupMixed(groupIndex)
    end if

end function isGroupMixed
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getGroupVarName(this,groupIndex) result(name)
    !! Return evolved variable name for given group index - works only for groups that are not mixed

    class(Model)              ,intent(in) :: this
    integer(ik)                ,intent(in) :: groupIndex
    character(:) ,allocatable              :: name

    integer(ik) :: trueGroupIndex ,termIndex ,trueTermIndex

    if (assertions) then 
        call assertPure (this%isDefined(),"Called getGroupVarName from undefined model object")
        call assertPure(groupIndex > 0 ,"Negative group index passed to getGroupVarName")
        call assertPure(groupIndex <= size(this%implicitTermGroup)+size(this%generalTermGroup),&
        "Out of range group index passed to getGroupVarName")

        if (groupIndex > size(this%implicitTermGroup)) then 
            trueGroupIndex = groupIndex - size(this%implicitTermGroup)
            call assertPure(.not. this%generalGroupMixed(trueGroupIndex),&
            "Attempted to get variable name for mixed term group in model")
        else
            call assertPure(.not. this%implicitGroupMixed(groupIndex),&
            "Attempted to get variable name for mixed term group in model")
        end if
        
    end if
    if (groupIndex > size(this%implicitTermGroup)) then
        trueGroupIndex = groupIndex - size(this%implicitTermGroup)
        termIndex = this%generalTermGroup(trueGroupIndex)%entry(1)
        if (termIndex > size(this%implicitTerms)) then 
            trueTermIndex = termIndex-size(this%implicitTerms)
            name = this%generalTerms(trueTermIndex)%entry%getVarName()
        else
            name = this%implicitTerms(termIndex)%entry%getVarName()
        end if
    else
        trueGroupIndex = groupIndex 
        termIndex = this%implicitTermGroup(trueGroupIndex)%entry(1)
        name = this%implicitTerms(termIndex)%entry%getVarName()
    end if

end function getGroupVarName
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateModelData(this,varCont,updatePriority)
    !! Update this model's modelbound data if allocated

    class(Model)            ,intent(inout)  :: this
    type(VariableContainer) ,intent(in)     :: varCont  !! Variable container used in the update
    integer(ik) ,optional   ,intent(in)     :: updatePriority !! Priority for this update call 

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update model data on undefined model")
        call assert(varCont%isDefined(),"Attempted to update model data using undefined variable container")
    end if

    if (allocated(this%modelData)) then 
        if (present(updatePriority)) then 
            call this%modelData%update(this,varCont,updatePriority=updatePriority)
        else
            call this%modelData%update(this,varCont)
        end if
    end if
end subroutine updateModelData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setModelData(this,modelData)
    !! Setter for modelData

    class(Model)            ,intent(inout)  :: this
    class(ModelboundData)   ,intent(in)     :: modelData

    if (assertions) call assert(modelData%isDefined(),"Attempted to set model data by passing undefined modelboundData object")

    if (allocated(this%modelData)) deallocate(this%modelData)
    allocate(this%modelData,source=modelData)

end subroutine setModelData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine copyModelData(this,modelData)
    !! Copy model data from this model into modelData

    class(Model)                        ,intent(in)     :: this
    class(ModelboundData) ,allocatable  ,intent(inout)  :: modelData

    if (allocated(this%modelData)) then 

        if (allocated(modelData)) deallocate(modelData)
        allocate(modelData,source=this%modelData)

    end if

end subroutine copyModelData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getImplicitTermRowData(this,termIndex) result(rowData)
    !! Return row data of implicit term with given term index

    class(Model)              ,intent(in) :: this
    integer(ik)               ,intent(in) :: termIndex
    type(SparseRowData)                   :: rowData

    if (assertions) then 
        call assertPure(this%isDefined(),"getImplicitTermRowData called on undefined model")
        call assertPure(this%assembled,"getImplicitTermRowData called on unassmbled model")
    end if

    rowData = this%implicitTerms(termIndex)%entry%getRowData()

end function getImplicitTermRowData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getImplicitTermIndexingData(this,termIndex) result(indexingData)
    !! Return indexing data of implicit term with given term index

    class(Model)              ,intent(in) :: this
    integer(ik)               ,intent(in) :: termIndex
    type(MatrixTermIndexingData)          :: indexingData

    if (assertions) then 
        call assertPure(this%isDefined(),"getImplicitTermIndexingData called on undefined model")
        call assertPure(this%assembled,"getImplicitTermIndexingData called on unassmbled model")
    end if

    indexingData = this%implicitTerms(termIndex)%entry%getIndexingData()

end function getImplicitTermIndexingData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine copyModelboundDataEntry(this,varName,container)
    !! Copy modelbound variable data with given name. Will throw error if no modelbound data is found

    class(Model)                          ,intent(in)    :: this 
    character(*)                          ,intent(in)    :: varName !! Name of data
    real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into 

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to retrieve model data entry from undefined model")
    end if
    call assert(allocated(this%modelData),"Attempted to retrieve model data entry from model with no modelbound data")

    call this%modelData%copyData(varName,container)

end subroutine copyModelboundDataEntry
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isTermNameRegistered(this,name) result(reg)
    !! Check whether term with given name is registered

    class(Model)        ,intent(in)  :: this
    character(*)        ,intent(in)  :: name
    logical                          :: reg

    integer(ik) :: i 

    if (assertions) call assertPure(this%isDefined(),"Attempted to get term name registration status from undefined model object")

    reg = this%isTermNameImplicit(name)

    if (.not. reg) then
        do i = 1,size(this%generalTermNames)
            if (this%generalTermNames(i)%string == name) then 
                reg = .true. 
                exit 
            end if
        end do
    end if

end function isTermNameRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isTermNameImplicit(this,name) result(reg)
    !! Check whether term with given name is implicit

    class(Model)        ,intent(in)  :: this
    character(*)        ,intent(in)  :: name
    logical                          :: reg

    integer(ik) :: i 

    if (assertions) call assertPure(this%isDefined(),"Attempted to get term name implicit status from undefined model object")

    reg = .false.
    do i = 1,size(this%implicitTermNames)
        if (this%implicitTermNames(i)%string == name) then 
            reg = .true. 
            exit 
        end if
    end do

end function isTermNameImplicit
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getImplicitTermIndex(this,name) result(ind)
    !! Get index of implicit term with given name

    class(Model)          ,intent(in) :: this
    character(*)         ,intent(in) :: name
    integer(ik)                      :: ind

    logical :: found
    integer(ik) :: i

    if (assertions) call assertPure(this%isDefined(),"getImplicitTermIndex called for undefined model")
    
    found = .false.

    do i = 1,size(this%implicitTermNames)
        if (this%implicitTermNames(i)%string == name) then 
            found = .true. 
            ind = i 
            exit 
        end if
    end do

    if (assertions) call assertPure(found,"Attempted to getImplicitTermIndex with name "//name//" not in model")

end function getImplicitTermIndex
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getGeneralTermIndex(this,name) result(ind)
    !! Get index of general term with given name

    class(Model)          ,intent(in) :: this
    character(*)         ,intent(in) :: name
    integer(ik)                      :: ind

    logical :: found
    integer(ik) :: i

    if (assertions) call assertPure(this%isDefined(),"getGeneralTermIndex called for undefined model")
    
    found = .false.

    do i = 1,size(this%generalTermNames)
        if (this%generalTermNames(i)%string == name) then 
            found = .true. 
            ind = i 
            exit 
        end if
    end do

    if (assertions) call assertPure(found,"Attempted to getGeneralTermIndex with name not in model")

end function getGeneralTermIndex
!-----------------------------------------------------------------------------------------------------------------------------------
module function evaluateTermByName(this,name,varCont) result(res)
    !! Evaluate a term by name

    class(Model)                         ,intent(in) :: this
    character(*)                         ,intent(in) :: name
    type(VariableContainer)              ,intent(in) :: varCont !! Variable container used to evaluate this term
    real(rk) ,allocatable ,dimension(:)              :: res

    integer(ik) :: termIndex

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to evaluate term in undefined model object")
        call assertPure(this%isAssembled(),"Attempted to evaluate term in unassembled model object")
        call assertPure(varCont%isDefined(),&
        "Attempted to evaluate term in model object by passing undefined variable container")
        call assertPure(this%isTermNameRegistered(name),"Attempted to evaluate term whose name is not registered in model object")
    end if

    if (this%isTermNameImplicit(name)) then 
        termIndex = this%getImplicitTermIndex(name)
        res = this%implicitTerms(termIndex)%entry%evaluate(varCont)
    else
        termIndex = this%getGeneralTermIndex(name)
        res = this%generalTerms(termIndex)%entry%evaluate(varCont)
    end if

end function evaluateTermByName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine calculateMatValsByName(this,name,varCont) 
    !! Calculate matrix values of a term by name

    class(Model)                         ,intent(inout) :: this
    character(*)                         ,intent(in)    :: name
    type(VariableContainer)              ,intent(in)    :: varCont 

    integer(ik) :: termIndex

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to calculate matrix of term in undefined model object")
        call assertPure(this%isAssembled(),"Attempted to calculate matrix of term in unassembled model object")
        call assertPure(varCont%isDefined(),&
        "Attempted to calculat matrix of term in model object by passing undefined variable container")
        call assertPure(this%isTermNameRegistered(name),&
            "Attempted to calculate matrix of term whose name is not registered in model object")
    end if

    if (this%isTermNameImplicit(name)) then 
        termIndex = this%getImplicitTermIndex(name)
        call this%implicitTerms(termIndex)%entry%calculateValues(varCont)
    end if

end subroutine calculateMatValsByName
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateTermByName(this,name,varCont)
    !! Update a term by name 

    class(Model)            ,intent(inout)  :: this
    character(*)            ,intent(in)     :: name
    type(VariableContainer) ,intent(in)     :: varCont 

    integer(ik) :: i ,trueGroupIndex ,termIndex

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to update term in undefined model object")
        call assertPure(this%isAssembled(),"Attempted to update term in unassembled model object")
        call assertPure(varCont%isDefined(),&
        "Attempted to update term in model object by passing undefined variable container")
    end if

    if (this%isTermNameImplicit(name)) then 
        termIndex = this%getImplicitTermIndex(name)
        if (allocated(this%modelData)) then 
            call this%implicitTerms(termIndex)%entry%update(varCont,this%modelData,hostModel=this)
        else
            call this%implicitTerms(termIndex)%entry%update(varCont,hostModel=this)
        end if
    else
        termIndex = this%getGeneralTermIndex(name)
        if (allocated(this%modelData)) then 
            call this%generalTerms(termIndex)%entry%update(varCont,this%modelData,hostModel=this)
        else
            call this%generalTerms(termIndex)%entry%update(varCont,hostModel=this)
        end if
    end if

end subroutine updateTermByName
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule model_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
