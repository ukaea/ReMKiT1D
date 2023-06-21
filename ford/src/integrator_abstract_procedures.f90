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
submodule (integrator_abstract_class) integrator_abstract_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains getters and setters for the abstract Integrator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setTimestep(this,timestep)
    !! Setter for dt

    class(Integrator)              ,intent(inout)  :: this
    real(rk)                       ,intent(in)     :: timestep

    this%dt = timestep

end subroutine setTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTimestep(this) result(timestep)
    !! Getter for dt
        
    class(Integrator)    ,intent(in) :: this
    real(rk)                         :: timestep

    if (assertions) call assertPure(this%isDefined(),"Attempted to get timestep from undefined Integrator objects")

    timestep = this%dt

end function getTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setTermGroups(this,groups)
    !! Setter for termGroups

    class(Integrator)               ,intent(inout)  :: this
    type(IntArray)   ,dimension(:)  ,intent(in)     :: groups

    this%termGroups = groups

end subroutine setTermGroups
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTermGroups(this) result(groups)
    !! Getter for termGroups

    class(Integrator)                        ,intent(in) :: this
    type(IntArray) ,allocatable ,dimension(:)            :: groups

    if (assertions) call assertPure(this%isDefined(),"Attempted to get term groups from undefined Integrator objects")

    groups = this%termGroups

end function getTermGroups
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setModelIndices(this,indices)
    !! Setter for modelIndices

    class(Integrator)            ,intent(inout)  :: this
    integer(ik)   ,dimension(:)  ,intent(in)     :: indices

    this%modelIndices = indices 

end subroutine setModelIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getModelIndices(this) result(indices)
    !! Getter for modelIndices

    class(Integrator)                     ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)            :: indices

    if (assertions) call assertPure(this%isDefined(),"Attempted to get model indices from undefined Integrator objects")

    indices = this%modelIndices

end function getModelIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setTimeEvolving(this,evo)
    !! Setter for evolvesTimeVar

    class(Integrator)          ,intent(inout)  :: this
    logical                    ,intent(in)     :: evo

    this%evolvesTimeVar = evo

end subroutine setTimeEvolving
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isTimeEvolving(this) result(evo)
    !! Check if this Integrator is allowed to evolve a time variable if present

    class(Integrator) ,intent(in) :: this
    logical                       :: evo

    if (assertions) call assertPure(this%isDefined(),"Attempted to check whether undefined Integrator evolves time")

    evo = this%evolvesTimeVar

end function isTimeEvolving
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setTimestepController(this,controller)
    !! Setter for dtController

    class(Integrator)          ,intent(inout)  :: this
    class(TimestepController)  ,intent(in)     :: controller

    if(allocated(this%dtController)) deallocate(this%dtController)
    allocate(this%dtController,source=controller)

end subroutine setTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function hasTimestepController(this) result(hasController)
    !! Check if this Integrator has an allocated timestep controller

    class(Integrator) ,intent(in) :: this
    logical                       :: hasController

    hasController = allocated(this%dtController)

end function hasTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
module function getTimestepFromController(this,inputVars) result(timestep)
    !! Get the individual timestep size if Integrator has a timestep controller

    class(Integrator)                     ,intent(inout) :: this 
    class(VariableContainer)              ,intent(in)    :: inputVars
    real(rk)                                             :: timestep

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to get timestep from timestep controller in undefined Integrator")
        call assertPure(inputVars%isDefined(),"Attempted to get timestep from timestep controller in Integrator by passing &
        &undefined inputVars")
        call assertPure(allocated(this%dtController),"getTimestepFromController called from Integrator which doesn't have an&
        & allocated timestep controller")
    end if

    timestep = this%dtController%evaluateTimestep(inputVars,this%dt)

end function getTimestepFromController
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setCommunicationNeeded(this,commNeeded)
    !! Setter for communicationNeeded

    class(Integrator)          ,intent(inout)  :: this
    logical                    ,intent(in)     :: commNeeded

    this%communicationNeeded = commNeeded

end subroutine setCommunicationNeeded
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isCommunicationNeeded(this) result(commNeeded)
    !! Check whether this Integrator requires MPI communication

    class(Integrator) ,intent(in) :: this
    logical                       :: commNeeded

    if (assertions) call assertPure(this%isDefined(),"Attempted to check whether communication is needed by undefined Integrator")

    commNeeded = this%communicationNeeded

end function isCommunicationNeeded
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setCommunicationData(this,commData)
    !! Setter for commData

    class(Integrator)            ,intent(inout)  :: this
    type(CommunicationData)      ,intent(in)     :: commData

    this%commData = commData

end subroutine setCommunicationData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCommunicationData(this) result(commData)
    !! Getter for commData

    class(Integrator)                     ,intent(in) :: this
    type(CommunicationData) ,allocatable              :: commData

    if (assertions) call assertPure(this%isDefined(),"Attempted to get communication data from undefined Integrator")

    commData = this%commData

end function getCommunicationData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setUpdateRules(this,updateRules)
    !! Setter for updateOnInternalIteration

    class(Integrator)                ,intent(inout)  :: this
    type(LogicalArray) ,dimension(:) ,intent(in)     :: updateRules

    this%updateOnInternalIteration = updateRules

end subroutine setUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getUpdateRules(this) result(updateRules)
    !! Getter for updateOnInternalIteration

    class(Integrator)                            ,intent(in) :: this
    type(LogicalArray) ,dimension(:) ,allocatable            :: updateRules

    if (assertions) call assertPure(this%isDefined(),"Attempted to get update rules from undefined Integrator")

    updateRules = this%updateOnInternalIteration 
    
end function getUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNonTrivialUpdate(this,nonTrivialUpdate)
    !! Setter for nonTrivialUpdate

    class(Integrator)          ,intent(inout)  :: this
    logical                    ,intent(in)     :: nonTrivialUpdate

    this%nonTrivialUpdate = nonTrivialUpdate

end subroutine setNonTrivialUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function hasNonTrivialUpdate(this) result(nonTrivialUpdate)
    !! Getter for nonTrivialUpdate

    class(Integrator) ,intent(in) :: this
    logical                       :: nonTrivialUpdate

    if (assertions) call assertPure(this%isDefined(),"Attempted to check whether undefined Integrator has non-trivial update rules")

    nonTrivialUpdate = this%nonTrivialUpdate 
    
end function hasNonTrivialUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setModelDataUpdateRules(this,updateRules)
    !! Setter for updateModelDataOnInternalIteration

    class(Integrator)     ,intent(inout)  :: this
    logical ,dimension(:) ,intent(in)     :: updateRules

    this%updateModelDataOnInternalIteration = updateRules

end subroutine setModelDataUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getModelDataUpdateRules(this) result(updateRules)
    !! Getter for updateModelDataOnInternalIteration

    class(Integrator)                 ,intent(in) :: this
    logical  ,dimension(:) ,allocatable           :: updateRules

    if (assertions) call assertPure(this%isDefined(),"Attempted to get model data update rules from undefined Integrator")

    updateRules = this%updateModelDataOnInternalIteration

end function getModelDataUpdateRules
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNonTrivialModelDataUpdate(this,nonTrivialUpdate)
    !! Setter for nonTrivialModelDataUpdate

    class(Integrator)          ,intent(inout)  :: this
    logical                    ,intent(in)     :: nonTrivialUpdate

    this%nonTrivialModelDataUpdate = nonTrivialUpdate

end subroutine setNonTrivialModelDataUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function hasNonTrivialModelDataUpdate(this) result(nonTrivialUpdate)
    !! Getter for nonTrivialModelDataUpdate

    class(Integrator) ,intent(in) :: this
    logical                       :: nonTrivialUpdate

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to check whether undefined Integrator has non-trivial model data update rules")

    nonTrivialUpdate = this%nonTrivialModelDataUpdate 

end function hasNonTrivialModelDataUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule integrator_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
