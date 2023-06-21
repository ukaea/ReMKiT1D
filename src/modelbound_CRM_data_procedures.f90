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
submodule (modelbound_CRM_data_class) modelbound_CRM_data_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the modelbound CRM data 

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initCRMData(this,numTransitions) 
    !! Collisional-radiative modelbound data initialization

    class(ModelboundCRMData)           ,intent(inout)  :: this
    integer(ik)                        ,intent(in)     :: numTransitions !! Expected number of transitions

    call this%makeDefined()

    allocate(this%transitions(numTransitions))
    this%numAddedTransitions = 0 
    this%allTransitionsAdded = .false.

end subroutine initCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addTransition(this,tr) 
    !! Add transition to CRM modelbound data

    class(ModelboundCRMData)           ,intent(inout)  :: this
    class(Transition)                  ,intent(in)     :: tr

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add transition to undefined modelbound CRM data object")
        call assertPure(tr%isDefined(),"Attempted to add undefined transition to modelbound CRM data object")
        call assertPure(.not. this%allTransitionsAdded,"Attempted to add transition to modelbound CRM object when no free &
        &transition slots available")
    end if

    this%numAddedTransitions = this%numAddedTransitions + 1
    allocate(this%transitions(this%numAddedTransitions)%entry,source=tr)
    if (this%numAddedTransitions == size(this%transitions)) this%allTransitionsAdded = .true.

end subroutine addTransition
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setInelData(this,inelData) 
    !! Setter for inelData

    class(ModelboundCRMData)        ,intent(inout)  :: this
    type(InelasticGridData)         ,intent(in)     :: inelData

    if (assertions) call assertPure(inelData%isDefined(),&
    "Undefined inelastic data passed to setInelData on modelbound CRM data object")

    allocate(this%inelData,source=inelData)

end subroutine setInelData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFixedW(this,ind) result(wMat)
    !! Return fixed mapping matrix with given index if inelastic data component is allocated

    class(ModelboundCRMData)       ,intent(in) :: this 
    integer(ik)                    ,intent(in) :: ind
    type(SparseRowData)                        :: wMat

    if (assertions) then 
        call assertPure(this%isDefined(),"getFixedW called from undefined modelbound CRM data object")
        call assertPure(allocated(this%inelData),&
        "getFixedW called from modelbound CRM data object with no allocated inelastic data")
    end if

    wMat = this%inelData%getFixedW(ind)

end function getFixedW
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFixedEmissionVector(this,ind) result(emit)
    !! Return fixed emission vector for mapping matrix with given index if inelastic data component is allocated

    class(ModelboundCRMData)       ,intent(in) :: this 
    integer(ik)                    ,intent(in) :: ind
    real(rk) ,allocatable ,dimension(:)        :: emit

    if (assertions) then 
        call assertPure(this%isDefined(),"getFixedEmissionVector called from undefined modelbound CRM data object")
        call assertPure(allocated(this%inelData),&
        "getFixedEmissionVector called from modelbound CRM data object with no allocated inelastic data")
    end if

    emit = this%inelData%getFixedEmissionVector(ind)

end function getFixedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine interpolateW(this,E,wRes) 
    !!  Interpolate mapping matrices for given energy and store in passed triangular matrix if inelastic data component is allocated. 
    !! Assumes upper triangular structure for wRes if E is positive and lower if it's negative. 

    class(ModelboundCRMData)     ,intent(in)    :: this
    real(rk)                     ,intent(in)    :: E !! Transition energy to interpolate for
    type(SparseRowData)          ,intent(inout) :: wRes !! Lower/upper triangular matrix to store the interpolated weights

    if (assertions) then 
        call assertPure(this%isDefined(),"interpolateW called from undefined modelbound CRM data object")
        call assertPure(allocated(this%inelData),&
        "interpolateW called from modelbound CRM data object with no allocated inelastic data")
    end if

    call this%inelData%interpolateW(E,wRes)

end subroutine interpolateW
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpolatedEmissionVector(this,E) result(emit)
!! Return interpolated emission vector for given input energy E if inelastic data component is allocated

    class(ModelboundCRMData)       ,intent(in) :: this 
    real(rk)                       ,intent(in) :: E
    real(rk) ,allocatable ,dimension(:)        :: emit

    if (assertions) then 
        call assertPure(this%isDefined(),"getInterpolatedEmissionVector called from undefined modelbound CRM data object")
        call assertPure(allocated(this%inelData),&
        "getInterpolatedEmissionVector called from modelbound CRM data object with no allocated inelastic data")
    end if

    emit = this%inelData%getInterpolatedEmissionVector(E)

end function getInterpolatedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTransitionIngoingStates(this,ind) result(inStates)
    !! Get ingoing states of transition with given index 

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                ,intent(in)  :: ind
    integer(ik) ,allocatable ,dimension(:)  :: inStates

    if (assertions) then 
        call assertPure(this%isDefined(),"getTransitionIngoingStates called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"getTransitionIngoingStates called on modelbound CRM data object before all &
        &transitions have been added")
    end if

    inStates = this%transitions(ind)%entry%getIngoingStates()

end function getTransitionIngoingStates  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTransitionOutgoingStates(this,ind) result(outStates)
    !! Get outgoing states of transition with given index

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                ,intent(in)  :: ind
    integer(ik) ,allocatable ,dimension(:)  :: outStates

    if (assertions) then 
        call assertPure(this%isDefined(),"getTransitionOutgoingStates called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"getTransitionOutgoingStates called on modelbound CRM data object before all &
        &transitions have been added")
    end if

    outStates = this%transitions(ind)%entry%getOutgoingStates()

end function getTransitionOutgoingStates 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTransitionRate(this,ind) result(rate)
    !! Get transition rate from transition with given index

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                ,intent(in)  :: ind
    real(rk)    ,allocatable ,dimension(:)  :: rate

    if (assertions) then 
        call assertPure(this%isDefined(),"getTransitionRate called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"getTransitionRate called on modelbound CRM data object before all &
        &transitions have been added")
    end if

    rate = this%transitions(ind)%entry%getRate()

end function getTransitionRate  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTransitionRateMomentum(this,ind) result(rate)
    !! Get momentum transfer rate from transition with given index

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                ,intent(in)  :: ind
    real(rk)    ,allocatable ,dimension(:)  :: rate

    if (assertions) then 
        call assertPure(this%isDefined(),"getTransitionRateMomentum called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"getTransitionRateMomentum called on modelbound CRM data object before all &
        &transitions have been added")
    end if

    rate = this%transitions(ind)%entry%getRateMomentum()

end function getTransitionRateMomentum  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTransitionRateEnergy(this,ind) result(rate)
    !! Get energy transfer rate from transition with given index

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                ,intent(in)  :: ind
    real(rk)    ,allocatable ,dimension(:)  :: rate

    if (assertions) then 
        call assertPure(this%isDefined(),"getTransitionRateEnergy called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"getTransitionRateEnergy called on modelbound CRM data object before all &
        &transitions have been added")
    end if

    rate = this%transitions(ind)%entry%getRateEnergy()

end function getTransitionRateEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionCrossSection(this,col,ind) result(crossSection)
    !! Get cross-section values from column col of the cross-section data of transition with given index

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)              ,intent(in)    :: col
    integer(ik)                ,intent(in)  :: ind
    real(rk)    ,allocatable ,dimension(:)  :: crossSection

    if (assertions) then 
        call assertPure(this%isDefined(),"getTransitionCrossSection called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"getTransitionCrossSection called on modelbound CRM data object before all &
        &transitions have been added")
    end if

    crossSection = this%transitions(ind)%entry%getCrossSectionCol(col)

end function getTransitionCrossSection  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getTransitionEnergy(this,ind) result(energy)
    !! Get transition energy of transition with given index

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                ,intent(in)  :: ind
    real(rk)    ,allocatable ,dimension(:)  :: energy

    if (assertions) then 
        call assertPure(this%isDefined(),"getTransitionEnergy called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"getTransitionEnergy called on modelbound CRM data object before all &
        &transitions have been added")
    end if

    energy = this%transitions(ind)%entry%getTransitionEnergy()

end function getTransitionEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumTransitions(this) result(numTrans)
    !! Get number of transitions registered in this object

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                             :: numTrans

    if (assertions) call assertPure(this%isDefined(),"getNumTransitions called on undefined modelbound CRM data object")

    numTrans = size(this%transitions)
    
end function getNumTransitions  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function ratesIncludeElDensity(this,ind) result(includesElDens)
    !! Check whether given transition includes electron density in the rate by default

    class(ModelboundCRMData)   ,intent(in)  :: this
    integer(ik)                ,intent(in)  :: ind
    logical                                 :: includesElDens

    if (assertions) call assertPure(this%isDefined(),"ratesIncludeElDensity called on undefined modelbound CRM data object")

    includesElDens = this%transitions(ind)%entry%includesElDensity()

end function ratesIncludeElDensity 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getPopulationChangeMatrix(this,ids,transitionIndices) result(popChangeMat)
    !! For a given set of species indices and transition indices returns a matrix whose entries are the change in particle number
    !! of a given species (by index) in transition processes determined by passed indices. Will not provide warnings if any
    !! particular ID is not found in any transition.

    class(ModelboundCRMData)    ,intent(in)  :: this
    integer(ik) ,dimension(:)   ,intent(in)  :: ids
    integer(ik) ,dimension(:)   ,intent(in)  :: transitionIndices
    integer(ik) ,dimension(:,:) ,allocatable :: popChangeMat

    integer(ik) :: i ,j

    integer(ik) ,allocatable ,dimension(:) :: inStates ,outStates

    if (assertions) then 
        call assertPure(this%isDefined(),"getPopulationChangeMatrix called on undefined modelbound CRM data object")
        call assertPure(all(transitionIndices>0),"Out of bounds (lower) transition index detected in getPopulationChangeMatrix")
        call assertPure(all(transitionIndices<=size(this%transitions)),&
                           "Out of bounds (upper) transition index detected in getPopulationChangeMatrix")
    end if

    allocate(popChangeMat(size(ids),size(transitionIndices)))

    do i = 1,size(transitionIndices)
        inStates = this%getTransitionIngoingStates(transitionIndices(i))
        outStates = this%getTransitionOutgoingStates(transitionIndices(i))
        do j = 1,size(ids)
            popChangeMat(j,i) = count(outStates==ids(j)) - count(inStates==ids(j))
        end do
    end do

end function getPopulationChangeMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRequiredDensityData(this,transitionIndex,removeLastState) result(densDataMat)
    !! For a given transition index returns ingoingState data as a (:,2) matrix, where the first column is the list of participating
    !! states and the second column the number of times that state appears. If an ID is zero (electrons) and the transition
    !! includes electron density in the rate, the second column value is reduced by one. If removeLastState is true, the corresponding
    !! number in the second column of the result is reduced by one (useful in implicit terms). 
    !! If any value in the second column drops to 0, the corresponding row will be removed. 
    !! The result can then be used to determine density variables and their powers required to convert the transitions rate values into a source.

    class(ModelboundCRMData)    ,intent(in)  :: this
    integer(ik)                 ,intent(in)  :: transitionIndex
    logical(ik) ,optional       ,intent(in)  :: removeLastState
    integer(ik) ,dimension(:,:) ,allocatable :: densDataMat

    logical :: rmLastState 
    integer(ik) ,allocatable ,dimension(:) :: inStates ,uniqueStates ,stateCount

    integer(ik) :: i ,k

    if (assertions) then 
        call assertPure(this%isDefined(),"getRequiredDensityData called on undefined modelbound CRM data object")
        call assertPure(transitionIndex>0,"Out of bounds (lower) transition index detected in getRequiredDensityData")
        call assertPure(transitionIndex<=size(this%transitions),&
                           "Out of bounds (upper) transition index detected in getRequiredDensityData")
    end if

    rmLastState = .false. 

    if (present(removeLastState)) rmLastState = removeLastState

    inStates = this%getTransitionIngoingStates(transitionIndex)

    uniqueStates = removeDupeInts(inStates)
    allocate(stateCount(size(uniqueStates)))

    do i = 1,size(stateCount)
        stateCount(i) = count(inStates==uniqueStates(i))
        if (uniqueStates(i) == 0 .and. this%ratesIncludeElDensity(transitionIndex)) stateCount(i) = stateCount(i) - 1
    end do

    if (rmLastState) stateCount(size(stateCount)) = stateCount(size(stateCount)) - 1

    allocate(densDataMat(count(stateCount>0),2))

    k = 1 

    do i = 1,size(stateCount)
        if (stateCount(i)>0) then 
            densDataMat(k,1) = uniqueStates(i)
            densDataMat(k,2) = stateCount(i)
            k = k + 1
        end if
    end do

end function getRequiredDensityData  
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateCRMData(this,hostModel,inputVars,updatePriority) 
    !!  Update modelbound data based on input variable container

    class(ModelboundCRMData)              ,intent(inout) :: this 
    class(ModelSurrogate)                 ,intent(in)    :: hostModel !! Host model passed to transitions
    class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate modelbound data
    integer(ik) ,optional                 ,intent(in)    :: updatePriority !! Priority for this update call 

    integer(ik) :: i

    if (assertions) then 
        call assertPure(this%isDefined(),"updateCRMData called on undefined modelbound CRM data object")
        call assertPure(this%allTransitionsAdded,"updateCRMData calledon modelbound CRM data object before all &
        &transitions have been added")
        call assertPure(inputVars%isDefined(),"Undefined inputVars passed to updateCRMData routine on modelbound CRM object")
    end if

    if (present(updatePriority)) then 
        do i = 1,size(this%transitions)
            call this%transitions(i)%entry%update(inputVars,hostModel,this,updatePriority=updatePriority)
        end do
    else
        do i = 1,size(this%transitions)
            call this%transitions(i)%entry%update(inputVars,hostModel,this)
        end do
    end if
end subroutine updateCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine crmCopy(this,name,container) 
    !! Retrieves data based on name format - assumes name of format "rate"//dataSpec//"index"//transIndex, where dataSpec is 0,1,2
    !! and corresponds to rate type (0-particles,1-momentum,2-energy), and transIndex is the transition index in the data or 
    !! "cssl//dataSpec//"index"//transIndex for the cross-section single harmonics where dataSpec is now the cross-section harmoinic index/column
    
    class(ModelboundCRMData)              ,intent(in)    :: this 
    character(*)                          ,intent(in)    :: name 
    real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container 

    integer(ik) :: i ,dataSpec ,transIndex ,indexStart ,csDim ,rateSize

    real(rk) ,allocatable ,dimension(:) :: csVals

    character(4) :: dataPrefix

    dataPrefix = name(1:4)
    
    select rank(container)
    rank (1)

        if (dataPrefix == "rate") then 

            read(name(5:5),*) dataSpec
            read(name(11:len(name)),*) transIndex

            select case(dataSpec)

            case (0)
                
                if (.not. allocated(container)) then 
                    allocate(container,source=this%getTransitionRate(transIndex))
                else
                    if (ubound(container,1) /= ubound(this%getTransitionRate(transIndex),1) .or. &
                        lbound(container,1) /= lbound(this%getTransitionRate(transIndex),1)) then 
                        deallocate(container)
                        allocate(container, source=this%getTransitionRate(transIndex))
                    else

                        container = this%getTransitionRate(transIndex)
                    end if
                end if
                
            case (1)
                if (.not. allocated(container)) then 
                    allocate(container,source=this%getTransitionRateMomentum(transIndex))
                else
                    if (ubound(container,1) /= ubound(this%getTransitionRateMomentum(transIndex),1) .or. &
                        lbound(container,1) /= lbound(this%getTransitionRateMomentum(transIndex),1)) then 
                        deallocate(container)
                        allocate(container, source=this%getTransitionRateMomentum(transIndex))
                    else

                        container = this%getTransitionRateMomentum(transIndex)
                    end if
                end if

            case (2)
                if (.not. allocated(container)) then 
                    allocate(container,source=this%getTransitionRateEnergy(transIndex))
                else
                    if (ubound(container,1) /= ubound(this%getTransitionRateEnergy(transIndex),1) .or. &
                        lbound(container,1) /= lbound(this%getTransitionRateEnergy(transIndex),1)) then 
                        deallocate(container)
                        allocate(container, source=this%getTransitionRateEnergy(transIndex))
                    else

                        container = this%getTransitionRateEnergy(transIndex)
                    end if
                end if

            case default 
                error stop "Unsupported dataSpec integer in crmCopy routine for rate case"
            end select

        else if (dataPrefix == "cssl") then

            indexStart = index(name, "index")

            call assert(indexStart>0,"unsupported data name in crmCopy routine - does not contain 'index' substring")

            read(name(5:indexStart-1),*) dataSpec
            read(name(indexStart+5:len(name)),*) transIndex

            rateSize = this%transitions(transIndex)%entry%getRateSize()
            csDim = this%transitions(transIndex)%entry%getCSDim()

            select case (csDim)

            case (1)

                csVals = this%getTransitionCrossSection(dataSpec,transIndex)

                if (.not. allocated(container)) then 
                    allocate(container(size(csVals)*rateSize))
                else
                    if (size(container) /= size(csVals)*rateSize) then 
                        deallocate(container)
                        allocate(container(size(csVals)*rateSize))
                    end if
                end if

                do i = 1,rateSize
                    container((i-1)*size(csVals)+1:i*size(csVals)) = csVals
                end do

            case (2)

                if (.not. allocated(container)) then 
                    allocate(container,source=this%getTransitionCrossSection(dataSpec,transIndex))
                else
                    if (ubound(container,1) /= ubound(this%getTransitionCrossSection(dataSpec,transIndex),1) .or. &
                        lbound(container,1) /= lbound(this%getTransitionCrossSection(dataSpec,transIndex),1)) then 
                        deallocate(container)
                        allocate(container, source=this%getTransitionCrossSection(dataSpec,transIndex))
                    else

                        container = this%getTransitionCrossSection(dataSpec,transIndex)
                    end if
                end if

            case (-1)

                error stop "crmCopy attempted to retrieve cross-section from transition without one"

            case default

                error stop "Unsupported dataSpec integer in crmCopy routine for cross-section case"

            end select


        else

            error stop "Unsupported data prefix detected in crmCopy routine"

        end if

    rank default
        error stop "container passed to crmCopy is not rank 1"

    end select

end subroutine crmCopy
!-----------------------------------------------------------------------------------------------------------------------------------
module function getDataDimCRM(this,name) result(dim)
    !! Get data dimensionality - currently always returns 1, assuming that the name is associated with a rate

    class(ModelboundCRMData)              ,intent(in)    :: this 
    character(*)                          ,intent(in)    :: name !! Name of data
    integer(ik)                                          :: dim

    character(4) :: dataPrefix

    dataPrefix = name(1:4)


    select case (dataPrefix)

    case ("rate")
        dim = 1
    case ("cssl")
        dim = 2
    case default 

        error stop "Unsupported data prefix detected in getDataDimCRM routine"

    end select

end function getDataDimCRM
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule modelbound_CRM_data_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
