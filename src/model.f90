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
module model_class
    !! author: Stefan Mijin
    !! 
    !! Houses general model class responsible for storing, directly evaluating and manipulating terms

    use data_kinds                     ,only: rk, ik
    use runtime_constants              ,only: debugging, assertions
    use god_objects                    ,only: Object
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use term_abstract_class            ,only: Term 
    use matrix_term_abstract_class     ,only: MatrixTerm ,MatrixTermIndexingData
    use variable_container_class       ,only: VariableContainer
    use petsc_controller_class         ,only: PETScController
    use support_types                  ,only: IntArray ,StringArray
    use model_surrogate_class          ,only: ModelSurrogate
    use modelbound_data_abstract_class ,only: ModelboundData
    use sparse_row_data_class          ,only: SparseRowData

    implicit none
    private

    type ,public :: MatTermContainer
        !! Container allowing for heterogeneous matrix term arrays
        class(MatrixTerm) ,allocatable :: entry
    end type MatTermContainer

    type ,public :: TermContainer
        !! Container allowing for heterogeneous general term arrays
        class(Term) ,allocatable :: entry
    end type TermContainer


    type ,public ,extends(ModelSurrogate) :: Model
        !! Object storing, evaluating, and manipulating various term objects, with terms categorized as either implicit (MatrixTerm objects)
        !! or general (Term objects). 

        integer(ik)                                       ,private :: numAddedMatrixTerms !! Tracker for number of allocated matrix terms
        integer(ik)                                       ,private :: numAddedGeneralTerms !! Tracker for number of allocated general terms
        type(IntArray)         ,allocatable ,dimension(:) ,private :: implicitTermGroup !! Groups of matrix terms used to identify terms for evaluation/manipulation
        type(IntArray)         ,allocatable ,dimension(:) ,private :: generalTermGroup !! Groups of general terms used to identify terms for evaluation/manipulation

        logical ,allocatable ,dimension(:) ,private :: implicitGroupMixed !! True for each implicit group that is either empty or 
                                                                          !! contains terms that evolve/evaluate more than one
                                                                          !! variable - used to forbid explicit evaluation of such terms
        logical ,allocatable ,dimension(:) ,private :: generalGroupMixed  !! True for each general group that that is either empty 
                                                                          !! or contains terms that evolve/evaluate more than one 
                                                                          !! variable - used to forbid explicit evaluation of such terms

        type(MatTermContainer) ,allocatable ,dimension(:) ,private :: implicitTerms !! Container for MatrixTerm objects living in this model
        type(TermContainer)    ,allocatable ,dimension(:) ,private :: generalTerms !! Container for general Term objects living in this model
 
        logical                                           ,private :: assembled !! True if the model is assembled and ready for use
        logical ,dimension(4)                             ,private :: setupCounter !! All true if the model is ready to accept terms

        class(ModelboundData) ,allocatable                ,private :: modelData !! Optional model data that can be used by contained terms

        type(StringArray) ,allocatable ,dimension(:)      ,private :: implicitTermNames !! Names of implicit terms
        type(StringArray) ,allocatable ,dimension(:)      ,private :: generalTermNames !! Names of general terms

        logical ,allocatable ,dimension(:)                ,private :: skipPattern !! Flag to skip pattern of particular implicit term when adding (for performance on startup)

        contains

        procedure ,public :: isAssembled 
        procedure ,public :: addImplicitTerm
        procedure ,public :: addGeneralTerm

        procedure ,public :: updateTermGroup
        procedure ,public :: evaluateTermGroup 
        procedure ,public :: evaluateTermByName
        procedure ,public :: calculateMatGroupValues
        procedure ,public :: addMatGroupValsToPETSc
        procedure ,public :: assemble

        procedure ,public :: isGroupMixed 
        procedure ,public :: getGroupVarName

        procedure ,public :: setModelData
        procedure ,public :: copyModelData
        procedure ,public :: updateModelData

        procedure ,public :: copyModelboundDataEntry

        procedure ,public :: setNumImplicitTerms
        procedure ,public :: setNumGeneralTerms 
        procedure ,public :: setNumImplicitGroups
        procedure ,public :: setNumGeneralGroups

        procedure ,public :: init => initModel

        procedure ,public :: getImplicitTermRowData
        procedure ,public :: getImplicitTermIndexingData

        procedure ,public :: isTermNameRegistered
        procedure ,public :: isTermNameImplicit

        procedure ,public :: getImplicitTermIndex
        procedure ,public :: getGeneralTermIndex

    end type Model
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initModel(this,numImplicitTerms,numGeneralTerms,numImplicitGroups,numGeneralGroups) 
        !! Model initialization routine

        class(Model)             ,intent(inout)  :: this
        integer(ik) ,optional    ,intent(in)     :: numImplicitTerms !! Number of MatrixTerm objects this model expects to be added
        integer(ik) ,optional    ,intent(in)     :: numGeneralTerms !! Number of general Term objects this model expects to be added
        integer(ik) ,optional    ,intent(in)     :: numImplicitGroups !! Number of implicit/matrix term groups registered with this model 
        integer(ik) ,optional    ,intent(in)     :: numGeneralGroups !! Number of general term groups registered with this model

    end subroutine initModel
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNumImplicitTerms(this,numImplicitTerms) 
        !! Set number of implicit terms and perform allocation

        class(Model)             ,intent(inout)  :: this
        integer(ik)              ,intent(in)     :: numImplicitTerms !! Number of MatrixTerm objects this model expects to be added

    end subroutine setNumImplicitTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNumGeneralTerms(this,numGeneralTerms) 
        !! Set number of general terms and perform allocation

        class(Model)             ,intent(inout)  :: this
        integer(ik)              ,intent(in)     :: numGeneralTerms !! Number of general Term objects this model expects to be added

    end subroutine setNumGeneralTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNumImplicitGroups(this,numImplicitGroups) 
        !! Set number of implicit groups and perform allocation

        class(Model)             ,intent(inout)  :: this
        integer(ik)              ,intent(in)     :: numImplicitGroups !! Number of implicit/matrix term groups registered with this model 

    end subroutine setNumImplicitGroups
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNumGeneralGroups(this,numGeneralGroups) 
        !! Set number of general groups and perform allocation

        class(Model)             ,intent(inout)  :: this
        integer(ik)              ,intent(in)     :: numGeneralGroups !! Number of general term groups registered with this model

    end subroutine setNumGeneralGroups
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isAssembled(this) result(assembled)
        !! Check if model is assembled and ready to use

        class(Model)  ,intent(in) :: this
        logical                   :: assembled

    end function isAssembled
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addImplicitTerm(this,impTerm,implicitGroups,generalGroups,termName,skipPattern)
        !! Add a MatrixTerm object to the model (deallocating the source!), and specify which implicit and general groups it belongs to

        class(Model)                          ,intent(inout)  :: this
        class(MatrixTerm) ,allocatable        ,intent(inout)  :: impTerm !! MatrixTerm object to be reallocated to this model
        integer(ik)             ,dimension(:) ,intent(in)     :: implicitGroups !! Implicit groups the added term should belong to
        integer(ik)             ,dimension(:) ,intent(in)     :: generalGroups !! General groups the added term should belong to
        character(*)                          ,intent(in)     :: termName !! Name of added term for indexing purposes
        logical  ,optional                    ,intent(in)     :: skipPattern !! True if the matrix term pattern should not be added to PETSc preallocation

    end subroutine addImplicitTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addGeneralTerm(this,genTerm,generalGroups,termName)
        !! Add a Term object to the model (deallocating the source!), and specify which general groups it belongs to

        class(Model)                   ,intent(inout)  :: this
        class(Term)       ,allocatable ,intent(inout)  :: genTerm !! General Term object to be reallocated to this model
        integer(ik)      ,dimension(:) ,intent(in)     :: generalGroups !! General groups the added term should belong to
        character(*)                   ,intent(in)     :: termName !! Name of added term for indexing purposes

    end subroutine addGeneralTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateTermGroup(this,groupIndex,varCont)
        !! Update a term group - if groupIndex > size(implicitGroup) it is taken to be in the general group

        class(Model)            ,intent(inout)  :: this
        integer(ik)             ,intent(in)     :: groupIndex !! Group index to be updated 
        type(VariableContainer) ,intent(in)     :: varCont !! Variable container used to update this group

    end subroutine updateTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function evaluateTermGroup(this,groupIndex,varCont) result(res)
        !! Evaluate a term group, returning the sum of all explicit results from the term group 
        !! - if groupIndex > size(implicitGroup) it is taken to be in the general group

        class(Model)                         ,intent(in) :: this
        integer(ik)                          ,intent(in) :: groupIndex !! Group index to evaluate
        type(VariableContainer)              ,intent(in) :: varCont !! Variable container used to evaluate this group
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine calculateMatGroupValues(this,groupIndex,varCont)
        !! Calculate matrix value in implicit term group given by groupIndex

        class(Model)            ,intent(inout)  :: this
        integer(ik)             ,intent(in)     :: groupIndex !! Group index of terms whose matrices should be calculated
        type(VariableContainer) ,intent(in)     :: varCont !! Variable container used for the calculation

    end subroutine calculateMatGroupValues
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addMatGroupValsToPETSc(this,groupIndex,petscCont,mult,petscGroup)
        !! Send off matrix values of given term group to the PETSc controller, multiplied by mult
    
        class(Model)            ,intent(in)    :: this
        integer(ik)             ,intent(in)    :: groupIndex !! Group index of terms whose matrices should be sent to PETSc
        type(PETScController)   ,intent(inout) :: petscCont !! PETScController object housing PETSc matrices 
        real(rk)                ,intent(in)    :: mult !! Multiplier used when adding matrices to PETSc - usually -dt 
        integer(ik) ,optional   ,intent(in)    :: petscGroup

    end subroutine addMatGroupValsToPETSc
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine assemble(this,petscCont)
        !! Assemble all of the matrix terms, preallocate PETScController objects and make sure model is ready for use

        class(Model)                     ,intent(inout) :: this
        type(PETScController)  ,optional ,intent(inout) :: petscCont !! Optional PETScController - should be present if the model has any implicitly evaluated terms

    end subroutine assemble
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isGroupMixed(this,groupIndex) result(mixed)
        !! Return true if group with given index is mixed

        class(Model)  ,intent(in) :: this
        integer(ik)   ,intent(in) :: groupIndex
        logical                   :: mixed

    end function isGroupMixed
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getGroupVarName(this,groupIndex) result(name)
        !! Return evolved variable name for given group index - works only for groups that are not mixed
 
        class(Model)              ,intent(in) :: this
        integer(ik)               ,intent(in) :: groupIndex
        character(:) ,allocatable             :: name

    end function getGroupVarName
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateModelData(this,varCont,updatePriority)
        !! Update this model's modelbound data if allocated

        class(Model)            ,intent(inout)  :: this
        type(VariableContainer) ,intent(in)     :: varCont  !! Variable container used in the update
        integer(ik) ,optional   ,intent(in)     :: updatePriority !! Priority for this update call 

    end subroutine updateModelData
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setModelData(this,modelData)
        !! Setter for modelData

        class(Model)            ,intent(inout)  :: this
        class(ModelboundData)   ,intent(in)     :: modelData

    end subroutine setModelData
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine copyModelData(this,modelData)
        !! Copy model data from this model into modelData

        class(Model)                        ,intent(in)     :: this
        class(ModelboundData) ,allocatable  ,intent(inout)  :: modelData

    end subroutine copyModelData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getImplicitTermRowData(this,termIndex) result(rowData)
        !! Return row data of implicit term with given term index
 
        class(Model)              ,intent(in) :: this
        integer(ik)               ,intent(in) :: termIndex
        type(SparseRowData)                   :: rowData

    end function getImplicitTermRowData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getImplicitTermIndexingData(this,termIndex) result(indexingData)
        !! Return indexing data of implicit term with given term index
 
        class(Model)              ,intent(in) :: this
        integer(ik)               ,intent(in) :: termIndex
        type(MatrixTermIndexingData)          :: indexingData

    end function getImplicitTermIndexingData
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine copyModelboundDataEntry(this,varName,container) 
        !! Copy modelbound variable data with given name. Will throw error if no modelbound data is found

        class(Model)                          ,intent(in)    :: this 
        character(*)                          ,intent(in)    :: varName !! Name of data
        real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into 

    end subroutine copyModelboundDataEntry
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isTermNameRegistered(this,name) result(reg)
        !! Check whether term with given name is registered

        class(Model)        ,intent(in)  :: this
        character(*)        ,intent(in)  :: name
        logical                          :: reg

    end function isTermNameRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isTermNameImplicit(this,name) result(reg)
        !! Check whether term with given name is implicit

        class(Model)        ,intent(in)  :: this
        character(*)        ,intent(in)  :: name
        logical                          :: reg

    end function isTermNameImplicit
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getImplicitTermIndex(this,name) result(ind)
        !! Get index of implicit term with given name

        class(Model)          ,intent(in) :: this
        character(*)         ,intent(in) :: name
        integer(ik)                      :: ind

    end function getImplicitTermIndex
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getGeneralTermIndex(this,name) result(ind)
        !! Get index of general term with given name

        class(Model)          ,intent(in) :: this
        character(*)         ,intent(in) :: name
        integer(ik)                      :: ind

    end function getGeneralTermIndex
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function evaluateTermByName(this,name,varCont) result(res)
        !! Evaluate a term by name

        class(Model)                         ,intent(in) :: this
        character(*)                         ,intent(in) :: name
        type(VariableContainer)              ,intent(in) :: varCont !! Variable container used to evaluate this term
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateTermByName
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module model_class
!-----------------------------------------------------------------------------------------------------------------------------------
 