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
module modelbound_data_varlike_class
    !! author: Stefan Mijin 
    !!
    !! Houses variable-like modelbound data object 

    use data_kinds                            ,only: rk ,ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use variable_container_class              ,only: VariableContainer ,CalculationRule
    use modelbound_data_abstract_class        ,only: ModelboundData
    use support_types                         ,only: RealArray ,IntArray
    use variable_list_class                   ,only: VariableList
    use partition_class                       ,only: Partition
    use indexing_class                        ,only: Indexing

    implicit none
    private

    type ,public ,extends(ModelboundData) :: ModelboundDataVarlike
        !! Modelbound data object containing variable-like data updated using a variable container

        type(RealArray)       ,allocatable ,dimension(:) ,private :: data        !! Values of modelbound data for each variable-like entry
        type(VariableList)                               ,private :: dataVarList !! List of names of data
        type(CalculationRule) ,allocatable ,dimension(:) ,private :: derivationRules !! Calculation rules for contained data
        type(IntArray)        ,allocatable ,dimension(:) ,private :: requiredDerivationIndices !! Lists of required variable indices for each calculation

        logical               ,allocatable ,dimension(:) ,private :: derivedFromMBData !! True for those modelbound variables where derivation 
                                                                                       !! rules use other modelbound data
        integer(ik) ,private :: maxDataPriority !! Highest priority value among data (lowest priority)
        contains

        procedure ,public :: init => initModelboundDataVarlike

        procedure ,public :: update => updateDataVarlike
        procedure ,public :: copyData => copyDataVarlike

        procedure ,public :: getDataDim => getDataDimVarlike

    end type ModelboundDataVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
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

        end subroutine initModelboundDataVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine updateDataVarlike(this,hostModel,inputVars,updatePriority) 
            !!  Update modelbound data based on input variable container

            class(ModelboundDataVarlike)          ,intent(inout) :: this 
            class(ModelSurrogate)                 ,intent(in)    :: hostModel !! Host model - unused
            class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate modelbound data
            integer(ik) ,optional                 ,intent(in)    :: updatePriority !! Priority for this update call (determines which variables are updated)

        end subroutine updateDataVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine copyDataVarlike(this,name,container) 
            !! Copy named modelbound data to passed container 

            class(ModelboundDataVarlike)          ,intent(in)    :: this 
            character(*)                          ,intent(in)    :: name !! Name of data
            real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into - must be rank 1

        end subroutine copyDataVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
        module function getDataDimVarlike(this,name) result(dim)
            !! Get data dimensionality (1 if fluid, 2 if single harmonic, 3 if distribution)

            class(ModelboundDataVarlike)          ,intent(in)    :: this 
            character(*)                          ,intent(in)    :: name !! Name of data
            integer(ik)                                          :: dim

        end function getDataDimVarlike
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module modelbound_data_varlike_class
!-----------------------------------------------------------------------------------------------------------------------------------
 