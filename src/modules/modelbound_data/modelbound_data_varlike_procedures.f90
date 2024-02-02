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
submodule (modelbound_data_varlike_class) modelbound_data_varlike_procedures
    !! author: Stefan Mijin
    !! 
    !! Contains module procedures associated with the variable-like modelbound data class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initModelboundDataVarlike(this,&
                                                dataList,&
                                                derivationRules,&
                                                partitionObj,&
                                                indexingObj,&
                                                xHaloWidth,&
                                                varCont, &
                                                procRank, &
                                                dataDerivIndices) 
                                                !! Varlike modelbound data initialization routine

    class(ModelboundDataVarlike)        ,intent(inout)  :: this
    type(VariableList)                  ,intent(in)     :: dataList !! Variable list object storing names of data 
    type(CalculationRule) ,dimension(:) ,intent(in)     :: derivationRules !! Calculation rules for each 1D data
    type(Partition)                     ,intent(in)     :: partitionObj !! Partition object used to initialize arrays
    type(Indexing)                      ,intent(in)     :: indexingObj !! Indexing object used to get numV and numH
    integer(ik)                         ,intent(in)     :: xHaloWidth !! Halo width in the x direction
    type(VariableContainer)             ,intent(in)     :: varCont !! Reference variable container for required derivation vars
    integer(ik)                         ,intent(in)     :: procRank !! Rank of the current process
    integer(ik) ,optional ,dimension(:) ,intent(in)     :: dataDerivIndices !! Data indices for which derivations require other modelbound data 

    integer(ik) :: i ,j 
    integer(ik) :: minX ,maxX ,numH ,numV ,locNumX

    if (assertions .or. assertionLvl >= 0) then 

        call assertPure(dataList%isDefined(),"Undefined data variable list passed to variable-like modelbound data constructor")
        call assertPure(partitionObj%isDefined(),"Undefined partition object passed to variable-like modelbound data constructor")
        call assertPure(indexingObj%isDefined(),"Undefined indexing object passed to variable-like modelbound data constructor")

        call assertPure(size(derivationRules) == dataList%getNumVars(),"Derivation rules passed to variable-like modelbound data &
        &constructor must conform in size to number of variables in data variable list")

        call assertPure(varCont%isDefined(),&
        "Undefined reference variable container passed to variable-like modelbound data constructor")

        do i = 1,size(derivationRules)
            call assertPure(derivationRules(i)%isDefined(),&
            "Undefined derivation rule passed to variable-like modelbound data constructor")
        end do

    end if

    this%dataVarList = dataList

    allocate(this%derivationRules,source=derivationRules)
    allocate(this%requiredDerivationIndices(size(derivationRules)))
    allocate(this%data(size(derivationRules)))
    allocate(this%derivedFromMBData(size(derivationRules)))

    this%derivedFromMBData = .false. 
    
    if (present(dataDerivIndices)) this%derivedFromMBData(dataDerivIndices) = .true.

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    locNumX = maxX - minX + 1
    numH = indexingObj%getNumH()
    numV = indexingObj%getNumV()

    do i = 1 ,size(derivationRules)
        if (dataList%isVarDist(i)) then
            allocate(this%data(i)%entry(1-xHaloWidth*numH*numV:(locNumX +xHaloWidth)*numH*numV))
        else if (dataList%isVarSingleHarmonic(i)) then
            allocate(this%data(i)%entry(1-xHaloWidth*numV:(locNumX +xHaloWidth)*numV))
        else if (dataList%isVarScalar(i)) then
            allocate(this%data(i)%entry(1))
        else
            allocate(this%data(i)%entry(1-xHaloWidth:locNumX +xHaloWidth))
        end if 
        this%data(i)%entry = 0
        if (allocated(derivationRules(i)%requiredVarNames)) then 
            allocate(this%requiredDerivationIndices(i)%entry(size(derivationRules(i)%requiredVarNames)))

            if (this%derivedFromMBData(i)) then 

                do j = 1 ,size(derivationRules(i)%requiredVarNames)
                    this%requiredDerivationIndices(i)%entry(j) = dataList%getVarIndex(derivationRules(i)%requiredVarNames(j)%string)
                end do

            else
                do j = 1 ,size(derivationRules(i)%requiredVarNames)
                    this%requiredDerivationIndices(i)%entry(j) = varCont%getVarIndex(derivationRules(i)%requiredVarNames(j)%string)
                end do
            end if
        else 
            allocate(this%requiredDerivationIndices(i)%entry(0))
        end if
    end do

    this%maxDataPriority = 0

    do i = 1,size(derivationRules)
        if (this%dataVarList%getVarPriority(i) > this%maxDataPriority) &
            this%maxDataPriority = this%dataVarList%getVarPriority(i)
    end do 

    call this%makeDefined()

end subroutine initModelboundDataVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateDataVarlike(this,hostModel,inputVars,updatePriority) 
    !!  Update modelbound data based on input variable container

    class(ModelboundDataVarlike)          ,intent(inout) :: this 
    class(ModelSurrogate)                 ,intent(in)    :: hostModel !! Host model - unused
    class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate modelbound data
    integer(ik) ,optional                 ,intent(in)    :: updatePriority !! Priority for this update call (determines which variables are updated)

    integer(ik) :: i ,usedPriority

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to updated undefined variable-like modelbound data")
        call assert(inputVars%isDefined(),"Attempted to updated variable-like modelbound data using undefined inputVars")
    end if

    usedPriority = this%maxDataPriority
    if (present(updatePriority)) usedPriority = updatePriority

    do i = 1, size(this%data)
        if (this%dataVarList%getVarPriority(i) <= usedPriority) then 
            if (this%derivedFromMBData(i)) then 
                this%data(i)%entry &
                = this%derivationRules(i)%derivationMethod%calculate(this%data,this%requiredDerivationIndices(i)%entry)
            else
                this%data(i)%entry &
                = this%derivationRules(i)%derivationMethod%calculate(inputVars%variables,this%requiredDerivationIndices(i)%entry)
            end if
        end if
    end do

end subroutine updateDataVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine copyDataVarlike(this,name,container) 
    !! Copy named modelbound data to passed container 

    class(ModelboundDataVarlike)          ,intent(in)    :: this 
    character(*)                          ,intent(in)    :: name !! Name of data
    real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into - must be rank 1

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to copy data from undefined variable-like modelbound data object")
        call assert(this%dataVarList%isVarNameRegistered(name),&
        "Name "//name//" passed to copyDataVarlike routine of variable-like modelbound data object&
        & is not registered")
    end if

    select rank (container)
    rank (1)
        !gfortran bug workaround with minimal allocation
        if (allocated(container)) then 
            if (ubound(container,1) /= ubound(this%data(this%dataVarList%getVarIndex(name))%entry,1) .or. &
               lbound(container,1) /= lbound(this%data(this%dataVarList%getVarIndex(name))%entry,1)) then 
                deallocate(container)
                allocate(container, source=this%data(this%dataVarList%getVarIndex(name))%entry)
               else
                container = this%data(this%dataVarList%getVarIndex(name))%entry
            end if
        else 
            allocate(container, source=this%data(this%dataVarList%getVarIndex(name))%entry)
        end if
    rank default 
        error stop "container passed to copyDataVarlike is not rank 1"
    end select

end subroutine copyDataVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
module function getDataDimVarlike(this,name) result(dim)
    !! Get data dimensionality (0 if scalar, 1 if fluid, 2 if single harmonic, 3 if distribution)

    class(ModelboundDataVarlike)          ,intent(in)    :: this 
    character(*)                          ,intent(in)    :: name !! Name of data
    integer(ik)                                          :: dim

    integer(ik) :: ind 

    ind = this%dataVarList%getVarIndex(name)

    if (this%dataVarList%isVarScalar(ind)) then 
        dim = 0
    else if (this%dataVarList%isVarSingleHarmonic(ind)) then
        dim = 2
    else if (this%dataVarList%isVarDist(ind)) then 
        dim = 3
    else
        dim = 1
    end if

end function getDataDimVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule modelbound_data_varlike_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
