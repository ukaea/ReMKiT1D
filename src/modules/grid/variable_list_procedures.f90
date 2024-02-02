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
submodule (variable_list_class) variable_list_procedures
    !! author: Stefan Mijin 
    !!
    !! Contains module procedures associated with the variable list class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initVarList(this) 
    !! Variable list initialization routine

    class(VariableList)           ,intent(inout)  :: this

    call this%makeDefined()
    allocate(this%names(0))
    allocate(this%distf(0))
    allocate(this%isSingleHarmonic(0))
    allocate(this%isScalar(0))
    allocate(this%isOnDualGrid(0))
    allocate(this%isStationary(0))
    allocate(this%priority(0))

end subroutine initVarList  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumVars(this) result(numVars)
    !! Returns number of variables in list

    class(VariableList)  ,intent(in) :: this
    integer(ik)                      :: numVars

    if (assertions) call assertPure(this%isDefined(),"Requested number of variables from undefined variable list")

    numVars = size(this%names)

end function getNumVars
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarNames(this) result(names)
    !! Getter of names

    class(VariableList)                          ,intent(in) :: this
    type(StringArray) ,allocatable ,dimension(:)             :: names

    if (assertions) call assertPure(this%isDefined(),"Requested variable name list from undefined variable list")

    names = this%names

end function getVarNames
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarName(this,ind) result(name)
    !! Return variable name at index ind

    class(VariableList)                  ,intent(in) :: this
    integer(ik)                          ,intent(in) :: ind
    character(:) ,allocatable                        :: name
    
    if (assertions) call assertPure(this%isDefined(),"Requested variable name from undefined variable list")

    name = this%names(ind)%string

end function getVarName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarPriority(this,ind) result(priority)
    !! Return priority of variable with given index 

    class(VariableList)  ,intent(in) :: this
    integer(ik)          ,intent(in) :: ind
    integer(ik)                      :: priority

    if (assertions) call assertPure(this%isDefined(),"Requested variable priority from undefined variable list")

    priority = this%priority(ind)

end function getVarPriority
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarNameRegistered(this,name) result(reg)
    !! Check whether variable with given name is registered

    class(VariableList) ,intent(in)  :: this
    character(*)        ,intent(in)  :: name
    logical                          :: reg

    integer(ik) :: i 

    if (assertions) call assertPure(this%isDefined(),"Attempted to get variable name registration status from undefined variable&
                                    & list object")

    reg = .false.
    do i = 1,size(this%names)
        if (this%names(i)%string == name) then 
            reg = .true. 
            exit 
        end if
    end do

end function isVarNameRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addVar(this,name,isDist,isSingleHarmonic,isScalar,isOnDualGrid,isStationary,priority) 
    !! Add variable with given name to list; isDist determines if variable is a distribution function and is .false. by default;
    !! isSingleHarmonic sets whether a given distribution is only a single harmonic - i.e. a function of just x and v
    !! isScalar tags variable as a scalar (will be stored as a dimension 1 array)
    !! isOnDualGrid marks variable as living on the dual/staggered grid (or having staggered harmonics if it's a distribution)
    !! isStationary marks variable as having d/dt=0
    !! priority is an integer governing operations such as variable derivation

    class(VariableList) ,intent(inout)  :: this
    character(*)        ,intent(in)     :: name
    logical ,optional   ,intent(in)     :: isDist
    logical ,optional   ,intent(in)     :: isSingleHarmonic
    logical ,optional   ,intent(in)     :: isScalar
    logical ,optional   ,intent(in)     :: isOnDualGrid
    logical ,optional   ,intent(in)     :: isStationary
    integer(ik) ,optional ,intent(in)   :: priority

    logical :: distf ,singleH ,scal ,dual,stat
    integer(ik) :: varPriority

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(this%isDefined(),"Attempted to add variable to undefined variable list")
        call assertPure(.not.(this%isVarNameRegistered(name)),&
            "Attempted to add variable with same name as one already in variable list")
    end if

    this%names = [this%names,stringArray(name)]
    distf = .false. 
    if (present(isDist)) distf = isDist
    singleH = .false.
    if (present(isSingleHarmonic)) singleH = isSingleHarmonic

    scal = .false. 
    if (present(isScalar)) scal = isScalar
    !Handle reserved name "time": 
    if (name == "time") scal = .true.

    dual = .false. 
    if (present(isOnDualGrid)) dual = isOnDualGrid

    stat = .false. 
    if (present(isStationary)) stat = isStationary
    
    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(.not. (distf .and. singleH),&
        "Cannot add variable to variable list that is both a distribution and a single harmonic")
        call assertPure(.not. ((distf .or. singleH) .and. scal),&
        "Cannot add variable to variable list that is a scalar and a distribution/single harmonic")
        call assertPure(.not. (scal .and. dual),"Cannotd add variable to variable list that is scalar and on dual grid")
    end if

    varPriority = 0 
    if (present(priority)) varPriority = priority
    this%distf = [this%distf,distf]
    this%isSingleHarmonic = [this%isSingleHarmonic,singleH]
    this%isScalar = [this%isScalar,scal]
    this%priority = [this%priority,varPriority]
    this%isOnDualGrid = [this%isOnDualGrid,dual]
    this%isStationary = [this%isStationary,stat]

end subroutine addVar
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarDist(this,ind) result(distf)
    !! Check whether variable with given index is a full distribution function

    class(VariableList) ,intent(in)  :: this
    integer(ik)         ,intent(in)  :: ind
    logical                          :: distf

    if (assertions) call assertPure(this%isDefined(),"isVarDist called for undefined variable list")

    distf = this%distf(ind)

end function isVarDist
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarSingleHarmonic(this,ind) result(singleH)
    !! Check whether variable with given index is a  single harmonic

    class(VariableList) ,intent(in)  :: this
    integer(ik)         ,intent(in)  :: ind
    logical                          :: singleH

    if (assertions) call assertPure(this%isDefined(),"isVarSingleHarmonic called for undefined variable list")

    singleH = this%isSingleHarmonic(ind)

end function isVarSingleHarmonic
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarScalar(this,ind) result(scal)
    !! Check whether variable with given index is a scalar

    class(VariableList) ,intent(in)  :: this
    integer(ik)         ,intent(in)  :: ind
    logical                          :: scal

    if (assertions) call assertPure(this%isDefined(),"isVarScalar called for undefined variable list")

    scal = this%isScalar(ind)

end function isVarScalar
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarOnDualGrid(this,ind) result(dual)
    !! Check whether variable with given index is a on dual grid

    class(VariableList) ,intent(in)  :: this
    integer(ik)         ,intent(in)  :: ind
    logical                          :: dual
    
    if (assertions) call assertPure(this%isDefined(),"isVarOnDualGrid called for undefined variable list")

    dual = this%isOnDualGrid(ind)

end function isVarOnDualGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarStationary(this,ind) result(stat)
    !! Check whether variable with given index is stationary

    class(VariableList) ,intent(in)  :: this
    integer(ik)         ,intent(in)  :: ind
    logical                          :: stat

    if (assertions) call assertPure(this%isDefined(),"isVarStationary called for undefined variable list")

    stat = this%isStationary(ind)

end function isVarStationary
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarIndex(this,name) result(ind)
    !! Get index of variable with given name

    class(VariableList)  ,intent(in) :: this
    character(*)        ,intent(in)  :: name
    integer(ik)                      :: ind

    logical :: found
    integer(ik) :: i

    if (assertions) call assertPure(this%isDefined(),"getVarIndex called for undefined variable list")
    
    found = .false.

    do i = 1,size(this%names)
        if (this%names(i)%string == name) then 
            found = .true. 
            ind = i 
            exit 
        end if
    end do

    if (assertions) call assertPure(found,"Attempted to getVarIndex with name not in variable list")

end function getVarIndex
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function combineWith(this,other) result(res)
    !! Combine two lists into one

    class(VariableList)  ,intent(in) :: this
    type(VariableList)   ,intent(in) :: other
    type(VariableList)               :: res

    if (assertions) then
        call assertPure(this%isDefined(),"combineWith called for undefined variable list")
        call assertPure(other%isDefined(),"combineWith called for undefined variable list")
    end if

    res%distf = [this%distf,other%distf]
    res%names = [this%names,other%names]
    res%isSingleHarmonic = [this%isSingleHarmonic,other%isSingleHarmonic]
    res%isScalar = [this%isScalar,other%isScalar]
    res%isOnDualGrid = [this%isOnDualGrid,other%isOnDualGrid]
    res%isStationary = [this%isStationary,other%isStationary]
    res%priority = [this%priority,other%priority]

    call res%makeDefined()

end function combineWith
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule variable_list_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
