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
submodule (variable_container_class) variable_container_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the variable container class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initCalculationRule(this,deriv,names) 
    !! Calculation rule object initialization routine

    class(CalculationRule)                    ,intent(inout)  :: this
    class(Derivation) ,optional               ,intent(in)     :: deriv !! Derivation component
    type(StringArray) ,optional ,dimension(:) ,intent(in)     :: names !! Required variable names

    if (assertions) then 
        call assertPure(present(deriv) .eqv. present(names),"If derivation object is passed to calculation rule constructor names&
        & of required variables must also be passed")
        if (present(deriv)) call assertPure(deriv%isDefined(),"Undefined derivation object passed to calculation rule constructor")
    end if
    if (present(names)) this%requiredVarNames = names 
    if (present(deriv)) allocate(this%derivationMethod,source=deriv)

    call this%makeDefined()

end subroutine initCalculationRule
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initVarContainer(this,&
                                        implicitVars,&
                                        derivedVars,&
                                        derivationRules,&
                                        indexingObj,&
                                        partitionObj,&
                                        xHaloWidth,&
                                        procRank) 
    !! Variable container initialization routine

    class(VariableContainer)            ,intent(inout)  :: this
    type(VariableList)                  ,intent(in)     :: implicitVars !! List of implicit variables
    type(VariableList)                  ,intent(in)     :: derivedVars !! List of derived variables
    type(CalculationRule) ,dimension(:) ,intent(in)     :: derivationRules  !! Derivation rules corresponding to derived variables
    type(Indexing)                      ,intent(in)     :: indexingObj !! Reference indexing object used to map between variables and implicit vectors
    type(Partition)                     ,intent(in)     :: partitionObj !! Reference partition object 
    integer(ik)                         ,intent(in)     :: xHaloWidth !! Width of halo in x-direction 
    integer(ik)                         ,intent(in)     :: procRank !! Current processor rank

    integer(ik) :: i ,j ,k ,l
    integer(ik) :: numH ,numV ,minX ,maxX ,minH ,maxH ,locNumX ,locNumH

    logical :: found 

    logical ,allocatable ,dimension(:) :: depthDetermined

    if (assertions) then 

        call assertPure(implicitVars%isDefined(),"Undefined implicit variable list passed to variable container constructor")
        call assertPure(derivedVars%isDefined(),"Undefined derived variable list passed to variable container constructor")
        call assertPure(indexingObj%isDefined(),"Undefined indexing object passed to variable container constructor")
        call assertPure(partitionObj%isDefined(),"Undefined partition object passed to variable container constructor")

        call assertPure(size(derivationRules) == derivedVars%getNumVars(),"Derivation rules passed to variable container &
        &constructor must conform in size to number of variables in derived variable list")

        do i = 1,size(derivationRules)
            call assertPure(derivationRules(i)%isDefined(),"Undefined derivation rule passed to variable container constructor")
            if (allocated(derivationRules(i)%requiredVarNames)) then
            do j = 1, size(derivationRules(i)%requiredVarNames)
                found = implicitVars%isVarNameRegistered(derivationRules(i)%requiredVarNames(j)%string) .or. &
                        derivedVars%isVarNameRegistered(derivationRules(i)%requiredVarNames(j)%string)
                call assertPure(found,"A required &
                &variable name -"//derivationRules(i)%requiredVarNames(j)%string//&
                "- in a derivation rule passed to variable container constructor was not found in the passed implicit &
                &or derived variable lists")
            end do
            end if
        end do

    end if

    this%implicitVarList = implicitVars
    this%derivedVarList = derivedVars

    allocate(this%derivationRules,source=derivationRules)
    allocate(this%requiredDerivationIndices(size(derivationRules)))

    do i = 1 ,size(derivationRules)
        if (allocated(derivationRules(i)%requiredVarNames)) then 
            allocate(this%requiredDerivationIndices(i)%entry(size(derivationRules(i)%requiredVarNames)))
            do j = 1 ,size(derivationRules(i)%requiredVarNames)
                if (implicitVars%isVarNameRegistered(derivationRules(i)%requiredVarNames(j)%string)) then

                    this%requiredDerivationIndices(i)%entry(j) = &
                    implicitVars%getVarIndex(derivationRules(i)%requiredVarNames(j)%string)

                else 

                    this%requiredDerivationIndices(i)%entry(j) = &
                    derivedVars%getVarIndex(derivationRules(i)%requiredVarNames(j)%string) + implicitVars%getNumVars()

                end if
            end do
        else 
            allocate(this%requiredDerivationIndices(i)%entry(0))
        end if
    end do

    allocate(this%implicitVarIndices(implicitVars%getNumVars()))
    allocate(this%implicitToLocIndex(implicitVars%getNumVars()))

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    minH = partitionObj%getMinHAtInd(procRank+1)
    maxH = partitionObj%getMaxHAtInd(procRank+1)
    locNumX = maxX - minX + 1
    locNumH = maxH - minH + 1
    numH = indexingObj%getNumH()
    numV = indexingObj%getNumV()

    allocate(this%variables(implicitVars%getNumVars()+derivedVars%getNumVars()))
    allocate(this%varLens(size(this%variables)))
    do i = 1,size(this%implicitVarIndices)
        if (implicitVars%isVarDist(i)) then 
            allocate(this%implicitVarIndices(i)%entry(locNumX*locNumH*numV))
            allocate(this%implicitToLocIndex(i)%entry(locNumX*locNumH*numV))
            allocate(this%variables(i)%entry(1-xHaloWidth*numH*numV:(locNumX+xHaloWidth)*numH*numV))
            this%variables(i)%entry = 0 
            this%variables(i)%entry(1-xHaloWidth*numH*numV:0) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells
            this%variables(i)%entry(locNumX*numH*numV+1:(locNumX+xHaloWidth)*numH*numV) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells 

            do j = minX,maxX
                do k = minH,maxH
                    do l = 1,numV 
                        this%implicitVarIndices(i)%entry(indexingObj%findDistIndex(j,k,l,.false.)) &
                        = indexingObj%findIndex(implicitVars%getVarName(i),j,k,l,local=.true.)
                        this%implicitToLocIndex(i)%entry(indexingObj%findDistIndex(j,k,l,.false.)) &
                        = indexingObj%findDistIndex(j,k,l,.true.)
                    end do
                end do
            end do

            this%varLens(i) = locNumX*numH*numV

        else if (implicitVars%isVarScalar(i)) then
            error stop "Implicit scalar variables in variable container are not supported"
        else if (minH == 1) then
            allocate(this%implicitVarIndices(i)%entry(locNumX))
            allocate(this%implicitToLocIndex(i)%entry(locNumX))
            allocate(this%variables(i)%entry(1-xHaloWidth:locNumX+xHaloWidth))
            this%variables(i)%entry = 0 
            this%variables(i)%entry(1-xHaloWidth:0) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells
            this%variables(i)%entry(locNumX+1:locNumX+xHaloWidth) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells 
            do j = 1,locNumX
                this%implicitVarIndices(i)%entry(j) = indexingObj%findIndex(implicitVars%getVarName(i),j+minX-1,local=.true.)
                this%implicitToLocIndex(i)%entry(j) = j
            end do

            this%varLens(i) = locNumX

        else 

            allocate(this%implicitVarIndices(i)%entry(0))
            allocate(this%implicitToLocIndex(i)%entry(0))
            allocate(this%variables(i)%entry(1-xHaloWidth:locNumX+xHaloWidth))
            this%variables(i)%entry = 0 
            this%variables(i)%entry(1-xHaloWidth:0) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells
            this%variables(i)%entry(locNumX+1:locNumX+xHaloWidth) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells 

            this%varLens(i) = locNumX

        end if
    end do

     do i = 1,size(derivationRules)
        if (derivedVars%isVarDist(i)) then 
            allocate(this%variables(i+size(this%implicitVarIndices))%entry(1-xHaloWidth*numH*numV:(locNumX+xHaloWidth)*numH*numV))
            this%variables(i+size(this%implicitVarIndices))%entry = 0 
            this%variables(i+size(this%implicitVarIndices))%entry(1-xHaloWidth*numH*numV:0) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells
            
            this%variables(i+size(this%implicitVarIndices))%entry(locNumX*numH*numV+1:(locNumX+xHaloWidth)*numH*numV) &
            = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells 

            this%varLens(i+size(this%implicitVarIndices)) = locNumX*numH*numV
            
        else if (derivedVars%isVarScalar(i)) then
            allocate(this%variables(i+size(this%implicitVarIndices))%entry(1))
            this%variables(i+size(this%implicitVarIndices))%entry = 0 

            this%varLens(i+size(this%implicitVarIndices)) = 1
        else
            allocate(this%variables(i+size(this%implicitVarIndices))%entry(1-xHaloWidth:locNumX+xHaloWidth))
            this%variables(i+size(this%implicitVarIndices))%entry = 0 
            this%variables(i+size(this%implicitVarIndices))%entry(1-xHaloWidth:0) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells
            this%variables(i+size(this%implicitVarIndices))%entry(locNumX+1:locNumX+xHaloWidth) = real(1,kind=rk) !Initialize halos to avoid NaNs in ghost cells 

            this%varLens(i+size(this%implicitVarIndices)) = locNumX
        end if
    end do

    this%maxDerivPriority = 0

    do i = 1,size(derivationRules)
        if (this%derivedVarList%getVarPriority(i) > this%maxDerivPriority) &
            this%maxDerivPriority = this%derivedVarList%getVarPriority(i)
    end do 

    allocate(this%derivationDepth(size(this%variables)))
    this%derivationDepth = -1

    allocate(depthDetermined(size(this%variables)))
    depthDetermined = .false.

    !Set depth to determined for all implicit variables

    depthDetermined(1:implicitVars%getNumVars()) = .true.

    !Loop until all depths are determined
    do 
        do i = 1,size(derivationRules)

            if (.not. depthDetermined(i+implicitVars%getNumVars())) then 
                if (size(this%requiredDerivationIndices(i)%entry) == 0) then
                    !If no variables are required set depth to 0
                    this%derivationDepth(i+implicitVars%getNumVars()) = 0 
                    depthDetermined(i+implicitVars%getNumVars()) = .true.
                else

                    if (all(depthDetermined(this%requiredDerivationIndices(i)%entry))) then

                        this%derivationDepth(i+implicitVars%getNumVars()) = &
                        maxval(this%derivationDepth(this%requiredDerivationIndices(i)%entry)) + 1

                        depthDetermined(i+implicitVars%getNumVars()) = .true.
                    end if

                end if
            end if
        end do

        if (all(depthDetermined)) exit
    end do
    call this%makeDefined()

end subroutine initVarContainer
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarIndex(this,name) result(ind)
    !! Get index of variable with given name

    class(VariableContainer)  ,intent(in) :: this
    character(*)              ,intent(in) :: name
    integer(ik)                           :: ind

    if (assertions) then 
        call assertPure(this%isDefined(),"Variable index requested from undefined variable container")
        call assertPure(this%implicitVarList%isVarNameRegistered(name) .or. this%derivedVarList%isVarNameRegistered(name),&
        "Variable name "//name//" requested from variable container not registered in lists associated with the container")
    end if
    ind = -1
    if (this%implicitVarList%isVarNameRegistered(name)) then 
        ind = this%implicitVarList%getVarIndex(name)
    else if (this%derivedVarList%isVarNameRegistered(name)) then 
        ind = this%implicitVarList%getNumVars()+this%derivedVarList%getVarIndex(name)
    end if

end function getVarIndex
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calculateDerivedVars(this,derivPriority,derivDepth)
    !! Calculate derived variables from implicit variables using derivation rules. If derivPriority is supplied only those variables with
    !! priority <= derivPriority will be derived. If derivDepth is present only variables at that derivation depth are calculated,
    !! otherwise all variables are calculated in derivation depth order

    class(VariableContainer) ,intent(inout)  :: this
    integer(ik) ,optional    ,intent(in)     :: derivPriority
    integer(ik) ,optional    ,intent(in)     :: derivDepth

    integer(ik) :: i ,j ,usedPriority ,currentVarPriority ,ind
    integer(ik) ,allocatable ,dimension(:) :: indicesAtDepth
    logical :: calculateVar

    if (assertions) call assertPure(this%isDefined(),"Attempted to calculate derived variables using undefined container")

    usedPriority = this%maxDerivPriority 
    if (present(derivPriority)) usedPriority = derivPriority

    if (present(derivDepth)) then 

        indicesAtDepth = findIndices(this%derivationDepth==derivDepth)

        do i = 1,size(indicesAtDepth)
            ind = indicesAtDepth(i)  - this%implicitVarList%getNumVars()

            currentVarPriority = this%derivedVarList%getVarPriority(ind)
            calculateVar = allocated(this%derivationRules(ind)%derivationMethod) .and. currentVarPriority <= usedPriority

            if (calculateVar) &
            this%variables(ind + this%implicitVarList%getNumVars())%entry &
            = this%derivationRules(ind)%derivationMethod%calculate(this%variables,this%requiredDerivationIndices(ind)%entry)
        end do

    else

        do j = 0,this%getMaxDepth()
            indicesAtDepth = findIndices(this%derivationDepth==j)

            do i = 1,size(indicesAtDepth)
                ind = indicesAtDepth(i) - this%implicitVarList%getNumVars()

                currentVarPriority = this%derivedVarList%getVarPriority(ind)
                calculateVar = allocated(this%derivationRules(ind)%derivationMethod) .and. currentVarPriority <= usedPriority

                if (calculateVar) &
                this%variables(ind + this%implicitVarList%getNumVars())%entry &
                = this%derivationRules(ind)%derivationMethod%calculate(this%variables,this%requiredDerivationIndices(ind)%entry)
            end do
        end do

    end if

end subroutine calculateDerivedVars
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine extractImplicitVars(this,vec)
    !! Extract implicit variables from a locally linearly indexed real vector

    class(VariableContainer)    ,intent(inout)  :: this
    real(rk)      ,dimension(:) ,intent(in)     :: vec

    integer(ik) :: i 

    if (assertions) call assertPure(this%isDefined(),"Attempted to extract implicit variables using undefined container")

    do i = 1,this%implicitVarList%getNumVars()
        this%variables(i)%entry(this%implicitToLocIndex(i)%entry) = vec(this%implicitVarIndices(i)%entry)
    end do

end subroutine extractImplicitVars
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine copyImplicitVarsToVec(this,vec,ignoreStationary)
    !! Copy variables into a locally linearly indexed real vector. If ignoreStationary is true copies any stationary variables as 0

    class(VariableContainer)    ,intent(inout)  :: this
    real(rk) ,dimension(:)      ,intent(inout)  :: vec
    logical ,optional           ,intent(in)     :: ignoreStationary

    integer(ik) :: i 

    logical :: ignoreStat

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to copy implicit vars to vector variables using undefined container")

    ignoreStat = .false.
    if (present(ignoreStationary)) ignoreStat = ignoreStationary
    do i = 1,this%implicitVarList%getNumVars()
        if (ignoreStat .and. this%implicitVarList%isVarStationary(i)) then 
            vec(this%implicitVarIndices(i)%entry) = 0
        else
            vec(this%implicitVarIndices(i)%entry) = this%variables(i)%entry(this%implicitToLocIndex(i)%entry)
        end if
    end do

end subroutine copyImplicitVarsToVec
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarDist(this,ind) result(distf)
    !! Check whether variable with given index is a distribution function

    class(VariableContainer) ,intent(in)  :: this
    integer(ik)              ,intent(in)  :: ind
    logical                               :: distf

    if (assertions) call assertPure(this%isDefined(),"isVarDist called from undefined container")

    if (ind <= this%implicitVarList%getNumVars()) then 
        distf = this%implicitVarList%isVarDist(ind)
    else
        distf = this%derivedVarList%isVarDist(ind-this%implicitVarList%getNumVars())
    end if

end function isVarDist
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarScalar(this,ind) result(scal)
    !! Check whether variable with given index is a scalar 

    class(VariableContainer) ,intent(in)  :: this
    integer(ik)              ,intent(in)  :: ind
    logical                               :: scal

    if (assertions) call assertPure(this%isDefined(),"isVarScalar called from undefined container")

    if (ind <= this%implicitVarList%getNumVars()) then 
        scal = this%implicitVarList%isVarScalar(ind)
    else
        scal = this%derivedVarList%isVarScalar(ind-this%implicitVarList%getNumVars())
    end if

end function isVarScalar
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarNameRegistered(this,name) result(reg)
    !! Check whether variable with given name is registered in either the implicit or derived list

    class(VariableContainer) ,intent(in)  :: this
    character(*)             ,intent(in)  :: name
    logical                               :: reg

    if (assertions) call assertPure(this%isDefined(),"isVarNameRegistered called from undefined variable container")

    reg = this%implicitVarList%isVarNameRegistered(name)

    if (.not. reg) reg = this%derivedVarList%isVarNameRegistered(name)

end function isVarNameRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarImplicit(this,name) result(imp)
    !! Check whether variable with given name is registered in either the implicit list

    class(VariableContainer) ,intent(in)  :: this
    character(*)             ,intent(in)  :: name
    logical                               :: imp

    if (assertions) call assertPure(this%isDefined(),"isVarImplict called from undefined variable container")

    imp = this%implicitVarList%isVarNameRegistered(name)

end function isVarImplicit
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isVarOnDualGrid(this,name) result(dual)
    !! Check whether variable with given name is on dual grid

    class(VariableContainer) ,intent(in)  :: this
    character(*)             ,intent(in)  :: name
    logical                               :: dual

    integer(ik) :: ind

    if (assertions) call assertPure(this%isDefined(),"isVarOnDualGrid called from undefined variable container")

    ind = this%getVarIndex(name)

    if (ind <= this%implicitVarList%getNumVars()) then 
        dual = this%implicitVarList%isVarOnDualGrid(ind)
    else
        dual = this%derivedVarList%isVarOnDualGrid(ind-this%implicitVarList%getNumVars())
    end if 

end function isVarOnDualGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarName(this,ind) result(name)
    !! Return variable name at index ind

    class(VariableContainer)             ,intent(in) :: this
    integer(ik)                          ,intent(in) :: ind
    character(:) ,allocatable                        :: name

    if (assertions) call assertPure(this%isDefined(),"getVarName called from undefined variable container")

    if (ind <= this%implicitVarList%getNumVars()) then 
        name = this%implicitVarList%getVarName(ind)
    else
        name = this%derivedVarList%getVarName(ind-this%implicitVarList%getNumVars())
    end if 
end function getVarName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getAllVarNames(this) result(names)
    !! Return all variable names in this container

    class(VariableContainer)             ,intent(in) :: this
    type(StringArray) ,allocatable ,dimension(:)     :: names

    integer(ik) :: i

    if (assertions) call assertPure(this%isDefined(),"getAllVarNames called from undefined variable container")

    allocate(names(size(this%variables)))

    do i = 1,size(names)
        names(i)%string = this%getVarName(i)
    end do

end function getAllVarNames
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getImplicitVarNames(this) result(names)
    !! Return all implicit variable names in this container

    class(VariableContainer)             ,intent(in) :: this
    type(StringArray) ,allocatable ,dimension(:)     :: names

    integer(ik) :: i

    if (assertions) call assertPure(this%isDefined(),"getImplicitVarNames called from undefined variable container")

    allocate(names(this%implicitVarList%getNumVars()))

    do i = 1,size(names)
        names(i)%string = this%getVarName(i)
    end do

end function getImplicitVarNames
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarDepth(this,name) result(depth)
    !! Get depth of variable with given name

    class(VariableContainer)  ,intent(in) :: this
    character(*)              ,intent(in) :: name
    integer(ik)                           :: depth

    if (assertions) call assertPure(this%isDefined(),"getVarDepth called from undefined variable container")

    depth = this%derivationDepth(this%getVarIndex(name))

end function getVarDepth
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxDepth(this) result(depth)
    !! Get max of derivationDepth

    class(VariableContainer)  ,intent(in) :: this
    integer(ik)                           :: depth

    if (assertions) call assertPure(this%isDefined(),"getMaxDepth called from undefined variable container")

    depth = maxval(this%derivationDepth)

end function getMaxDepth
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isStationary(this,name) result(stationary)
    !! Check whether variable with given name is stationary

    class(VariableContainer) ,intent(in)  :: this
    character(*)             ,intent(in)  :: name
    logical                               :: stationary

    integer(ik) :: ind

    if (assertions) call assertPure(this%isDefined(),"isStationary called from undefined variable container")

    ind = this%getVarIndex(name)

    if (ind <= this%implicitVarList%getNumVars()) then 
        stationary = this%implicitVarList%isVarStationary(ind)
    else
        stationary = this%derivedVarList%isVarStationary(ind-this%implicitVarList%getNumVars())
    end if 

end function isStationary
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine copyNamedVarsToVec(this,vec,names)
    !! Copy variables into locally indexed vector by name 
    
    class(VariableContainer)              ,intent(inout)  :: this
    real(rk) ,allocatable, dimension(:)   ,intent(inout)  :: vec
    type(StringArray) ,dimension(:)       ,intent(in) :: names 
    
    integer(ik) :: i, totalLen, offset
    integer(ik) ,allocatable ,dimension(:) :: varIndices, varLens

    if (assertions) call assertPure(this%isDefined(),"copyNamedVarsToVec called from undefined variable container")
    totalLen = 0 
    allocate(varIndices(size(names)))
    allocate(varLens(size(names)))
    do i=1,size(names)
        
        if (assertions) call assertPure(this%isVarNameRegistered(names(i)%string),&
                                        "copyNamedVarsToVec called with name not registered in container")

        varIndices(i) = this%getVarIndex(names(i)%string)
        varLens(i) = this%varLens(varIndices(i))

        totalLen = totalLen + varLens(i)
    end do

    if (.not. allocated(vec)) then 
        allocate(vec(totalLen))
    else if (size(vec) /= totalLen) then
        deallocate(vec)
        allocate(vec(totalLen))
    end if 

    offset = 0
    do i=1,size(names) 
        vec(offset+1:offset+varLens(i)) = this%variables(varIndices(i))%entry(1:varLens(i))
        offset = offset + varLens(i)
    end do 
end subroutine copyNamedVarsToVec
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine copyNamedVarsFromVec(this,vec,names)
    !! Copy variables from locally indexed vector by name 
    
    class(VariableContainer)                 ,intent(inout)  :: this
    real(rk)  ,dimension(:)                  ,intent(in)  :: vec
    type(StringArray) ,dimension(:)          ,intent(in) :: names 

    integer(ik) :: i, totalLen, offset
    integer(ik) ,allocatable ,dimension(:) :: varIndices, varLens

    if (assertions) call assertPure(this%isDefined(),"copyNamedVarsFromVec called from undefined variable container")
    totalLen = 0 
    allocate(varIndices(size(names)))
    allocate(varLens(size(names)))
    do i=1,size(names)
        
        if (assertions) call assertPure(this%isVarNameRegistered(names(i)%string),&
            "copyNamedVarsFromVec called with name not registered in container")

        varIndices(i) = this%getVarIndex(names(i)%string)
        varLens(i) = this%varLens(varIndices(i))

        totalLen = totalLen + varLens(i)
    end do

    if (assertions) call assertPure(size(vec) == totalLen, "copyNamedVarsFromVec called with vector of non-conforming size")
    
    offset = 0
    do i=1,size(names) 
        this%variables(varIndices(i))%entry(1:varLens(i)) = vec(offset+1:offset+varLens(i)) 
        offset = offset + varLens(i)
    end do 
end subroutine copyNamedVarsFromVec
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVarLens(this,names) result(lens)
    !! Get lengths of variable data vectors (not including halos) based on a list of names
    
    class(VariableContainer)                 ,intent(in)  :: this
    type(StringArray) ,dimension(:)          ,intent(in)  :: names 
    integer(ik) ,allocatable ,dimension(:)                :: lens 

    integer(ik) :: i
    allocate(lens(size(names)))
    do i=1,size(names)
        
        if (assertions) call assertPure(this%isVarNameRegistered(names(i)%string),&
                                        "copyNamedVarsToVec called with name not registered in container")
        
        lens(i) = this%varLens(this%getVarIndex(names(i)%string))

    end do
end function getVarLens
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule variable_container_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
