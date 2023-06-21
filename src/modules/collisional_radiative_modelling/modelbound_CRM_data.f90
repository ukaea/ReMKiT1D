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
module modelbound_CRM_data_class
    !! author: Stefan Mijin 
    !!
    !! Houses derived modelbound data class responsible for collisional-radiative data

    use data_kinds                           ,only: rk, ik
    use runtime_constants                    ,only: debugging, assertions
    use god_objects                          ,only: Object
    use assertion_utility                    ,only: assert, assertIdentical, assertPure
    use sparse_row_data_class                ,only: SparseRowData      
    use modelbound_data_abstract_class       ,only: ModelboundData
    use inelastic_grid_data_class            ,only: InelasticGridData 
    use transition_abstract_class            ,only: Transition
    use model_surrogate_class                ,only: ModelSurrogate 
    use variable_container_class             ,only: VariableContainer
    use physics_functions
    use support_functions

    implicit none
    private

    type ,public :: TransitionContainer
        !! Container allowing for heterogeneous Transition object arrays 
        class(Transition) ,allocatable :: entry
    end type TransitionContainer

    type ,public ,extends(ModelboundData) :: ModelboundCRMData
        !! Object responsible for storing and providing centralized access to transition objects and, optionally, inelastic grid data

        type(TransitionContainer) ,allocatable ,dimension(:) ,private :: transitions 
        type(InelasticGridData)   ,allocatable               ,private :: inelData 

        integer(ik)                                          ,private :: numAddedTransitions !! Tracker for number of added transitions
        logical                                              ,private :: allTransitionsAdded !! True if all transitions are added and the data can be used 

        contains

        procedure ,public :: addTransition 
        procedure ,public :: setInelData

        procedure ,public :: getFixedW
        procedure ,public :: getFixedEmissionVector

        procedure ,public :: interpolateW
        procedure ,public :: getInterpolatedEmissionVector

        procedure ,public :: getTransitionIngoingStates
        procedure ,public :: getTransitionOutgoingStates

        procedure ,public :: getTransitionRate 

        procedure ,public :: getTransitionRateMomentum

        procedure ,public :: getTransitionRateEnergy 

        procedure ,public :: getTransitionCrossSection

        procedure ,public :: getTransitionEnergy

        procedure ,public :: getNumTransitions

        procedure ,public :: ratesIncludeElDensity

        procedure ,public :: getPopulationChangeMatrix
        procedure ,public :: getRequiredDensityData

        procedure ,public :: update => updateCRMData 

        procedure ,public :: copyData => crmCopy

        procedure ,public :: getDataDim => getDataDimCRM

        procedure ,public :: init => initCRMData

    end type ModelboundCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initCRMData(this,numTransitions) 
        !! Collisional-radiative modelbound data initialization

        class(ModelboundCRMData)           ,intent(inout)  :: this
        integer(ik)                        ,intent(in)     :: numTransitions !! Expected number of transitions

    end subroutine initCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addTransition(this,tr) 
        !! Add transition to CRM modelbound data

        class(ModelboundCRMData)           ,intent(inout)  :: this
        class(Transition)                  ,intent(in)     :: tr

    end subroutine addTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setInelData(this,inelData) 
        !! Setter for inelData

        class(ModelboundCRMData)        ,intent(inout)  :: this
        type(InelasticGridData)         ,intent(in)     :: inelData

    end subroutine setInelData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getFixedW(this,ind) result(wMat)
        !! Return fixed mapping matrix with given index if inelastic data component is allocated

        class(ModelboundCRMData)       ,intent(in) :: this 
        integer(ik)                    ,intent(in) :: ind
        type(SparseRowData)                        :: wMat
        
    end function getFixedW
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getFixedEmissionVector(this,ind) result(emit)
        !! Return fixed emission vector for mapping matrix with given index if inelastic data component is allocated

        class(ModelboundCRMData)       ,intent(in) :: this 
        integer(ik)                    ,intent(in) :: ind
        real(rk) ,allocatable ,dimension(:)        :: emit
        
    end function getFixedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine interpolateW(this,E,wRes) 
        !!  Interpolate mapping matrices for given energy and store in passed triangular matrix if inelastic data component is allocated. 
        !! Assumes upper triangular structure for wRes if E is positive and lower if it's negative. 

        class(ModelboundCRMData)     ,intent(in)    :: this
        real(rk)                     ,intent(in)    :: E !! Transition energy to interpolate for
        type(SparseRowData)          ,intent(inout) :: wRes !! Lower/upper triangular matrix to store the interpolated weights

    end subroutine interpolateW
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getInterpolatedEmissionVector(this,E) result(emit)
        !! Return interpolated emission vector for given input energy E if inelastic data component is allocated

        class(ModelboundCRMData)       ,intent(in) :: this 
        real(rk)                       ,intent(in) :: E
        real(rk) ,allocatable ,dimension(:)        :: emit
        
    end function getInterpolatedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionIngoingStates(this,ind) result(inStates)
        !! Get ingoing states of transition with given index 

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                ,intent(in)  :: ind
        integer(ik) ,allocatable ,dimension(:)  :: inStates

    end function getTransitionIngoingStates  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionOutgoingStates(this,ind) result(outStates)
        !! Get outgoing states of transition with given index

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                ,intent(in)  :: ind
        integer(ik) ,allocatable ,dimension(:)  :: outStates

    end function getTransitionOutgoingStates 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionRate(this,ind) result(rate)
        !! Get transition rate from transition with given index

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                ,intent(in)  :: ind
        real(rk)    ,allocatable ,dimension(:)  :: rate

    end function getTransitionRate  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionRateMomentum(this,ind) result(rate)
        !! Get momentum transfer rate from transition with given index

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                ,intent(in)  :: ind
        real(rk)    ,allocatable ,dimension(:)  :: rate

    end function getTransitionRateMomentum  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionRateEnergy(this,ind) result(rate)
        !! Get energy transfer rate from transition with given index

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                ,intent(in)  :: ind
        real(rk)    ,allocatable ,dimension(:)  :: rate

    end function getTransitionRateEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionCrossSection(this,col,ind) result(crossSection)
        !! Get cross-section values from column col of the cross-section data of transition with given index

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)              ,intent(in)    :: col
        integer(ik)                ,intent(in)  :: ind
        real(rk)    ,allocatable ,dimension(:)  :: crossSection

    end function getTransitionCrossSection  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getTransitionEnergy(this,ind) result(energy)
        !! Get transition energy of transition with given index

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                ,intent(in)  :: ind
        real(rk)    ,allocatable ,dimension(:)  :: energy

    end function getTransitionEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNumTransitions(this) result(numTrans)
        !! Get number of transitions registered in this object

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                             :: numTrans

    end function getNumTransitions  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function ratesIncludeElDensity(this,ind) result(includesElDens)
        !! Check whether given transition includes electron density in the rate by default

        class(ModelboundCRMData)   ,intent(in)  :: this
        integer(ik)                ,intent(in)  :: ind
        logical                                 :: includesElDens

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

    end function getRequiredDensityData  
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateCRMData(this,hostModel,inputVars,updatePriority) 
        !!  Update modelbound data based on input variable container

        class(ModelboundCRMData)              ,intent(inout) :: this 
        class(ModelSurrogate)                 ,intent(in)    :: hostModel !! Host model passed to transitions
        class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate modelbound data
        integer(ik) ,optional                 ,intent(in)    :: updatePriority !! Priority for this update call 

    end subroutine updateCRMData
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine crmCopy(this,name,container) 
        !! Retrieves data based on name format - assumes name of format "rate"//dataSpec//"index"//transIndex, where dataSpec is 0,1,2
        !! and corresponds to rate type (0-particles,1-momentum,2-energy), and transIndex is the transition index in the data

        class(ModelboundCRMData)              ,intent(in)    :: this 
        character(*)                          ,intent(in)    :: name 
        real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container 

    end subroutine crmCopy
!-----------------------------------------------------------------------------------------------------------------------------------
    module function getDataDimCRM(this,name) result(dim)
        !! Get data dimensionality - currently always returns 1, assuming that the name is associated with a rate

        class(ModelboundCRMData)              ,intent(in)    :: this 
        character(*)                          ,intent(in)    :: name !! Name of data
        integer(ik)                                          :: dim

    end function getDataDimCRM
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module modelbound_CRM_data_class
!-----------------------------------------------------------------------------------------------------------------------------------
 