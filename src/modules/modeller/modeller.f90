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
module modeller_class
    !! author: Stefan Mijin 
    !! 
    !! Houses modeller class responsible for controlling data manipulation and integration

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use variable_container_class    ,only: VariableContainer
    use petsc_controller_class      ,only: PETScController
    use mpi_controller_class        ,only: MPIController ,CommunicationData
    use model_class                 ,only: Model
    use manipulator_abstract_class  ,only: Manipulator 
    use composite_manipulator_class ,only: CompositeManipulator
    use modeller_surrogate_class    ,only: ModellerSurrogate
    use support_types               ,only: StringArray
    use composite_integrator_class  ,only: CompositeIntegrator
    use sparse_row_data_class       ,only: SparseRowData
    use indexing_class              ,only: Indexing
    use support_types
    use status_printing

    implicit none
    private

    type ,public :: ModelContainer
        !! Container for Model objects allowing heterogeneous arrays
        class(Model) ,allocatable :: entry
    end type ModelContainer

    type ,public ,extends(ModellerSurrogate) :: Modeller
        !! Central object responsible for storing variables and models, as well as doing integration and data manipulation

        type(VariableContainer)                          ,private :: vars !! Main variable container
        type(ModelContainer)  ,allocatable ,dimension(:) ,private :: models !! Array holding all models this modeller has access to
        type(PETScController) ,allocatable               ,private :: petscCont !! PETScController object used for solving linear systems
        type(MPIController)                              ,private :: MPICont !! MPIController used for all communication
        type(CommunicationData) ,allocatable             ,private :: commData !! Default communication data containing MPI communication instruction

        type(SparseRowData)   ,allocatable               ,private :: identityMat !! Sparse identity matrix used for PETSc calls

        class(Manipulator)    ,allocatable               ,private :: integ !! Main integrator object 

        class(CompositeManipulator) ,allocatable         ,private :: manip !! Secondary manipulator object for optional data manipulation

        integer(ik)                                      ,private :: numModelsAdded !! Tracked of number of models allocated

        logical(ik)                                      ,private :: assembled !! True if modeller has been assembled and is ready for use

        contains

        procedure ,public :: isAssembled 

        procedure ,public :: addModel

        procedure ,public :: copyVarValuesTo
        procedure ,public :: copyVarValuesFrom
        procedure ,public :: callManipulator
        procedure ,public :: integrate

        procedure ,public :: setIntegrator
        procedure ,public :: setManipulator

        procedure ,public :: calculateIdentityMat

        procedure ,public :: updateModelTermGroup
        procedure ,public :: updateModelTermByName
        procedure ,public :: updateModelData
        procedure ,public :: calculateMatGroupValsInModel
        procedure ,public :: calculateMatValsByTermName
        procedure ,public :: addModelMatGroupToPETSc
        procedure ,public :: linearSolvePETSc
        procedure ,public :: evaluateModelTermGroup
        procedure ,public :: evaluateModelTermByName

        procedure ,public :: copyDataFromModel

        procedure ,public :: getEvolvedVarInTermGroup
        procedure ,public :: isModelTermGroupMixed
        procedure ,public :: getCurrentTime

        procedure ,public :: performComm
        procedure ,public :: safeCommAndDeriv
        procedure ,public :: isTrueEverywhere

        procedure ,public :: assemble

        procedure ,public :: init => initModeller 

    end type Modeller
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initModeller(this,numModels,initVars,mpiCont,petscCont,commData) 
        !! Modeller initialization routine

        class(Modeller)                   ,intent(inout) :: this
        integer(ik)                       ,intent(in)    :: numModels !! Number of models this modeller expects to be added
        type(VariableContainer)           ,intent(in)    :: initVars !! Initial variable container
        type(MPIController)               ,intent(in)    :: mpiCont !! Main MPIController 
        type(PETScController)   ,optional ,intent(in)    :: petscCont !! Optional PETSc controller - should be supplied if any integration/manipulation routine uses PETSc
        type(CommunicationData) ,optional ,intent(in)    :: commData !! Default MPI communication data

    end subroutine initModeller
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addModel(this,newModel)
        !! Add a Model to this modeller

        class(Modeller)                  ,intent(inout)  :: this
        class(Model) ,allocatable        ,intent(inout)  :: newModel

    end subroutine addModel
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isAssembled(this) result(assembled)
        !! Return true if modeller is assembled and ready to evolve variables

        class(Modeller)  ,intent(in) :: this
        logical                      :: assembled

    end function isAssembled
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine copyVarValuesFrom(this,varCont)
        !! Copy values of variables from outside variable container into this modeller's container

        class(Modeller)          ,intent(inout)  :: this
        type(VariableContainer)  ,intent(in)     :: varCont

    end subroutine copyVarValuesFrom
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine copyVarValuesTo(this,varCont)
        !! Copy values of variables to outside variable container from this modeller's container

        class(Modeller)          ,intent(inout)  :: this
        type(VariableContainer)  ,intent(inout)  :: varCont

    end subroutine copyVarValuesTo
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setIntegrator(this,integ)
        !! Setter for integrator object

        class(Modeller)                          ,intent(inout)  :: this
        class(Manipulator)                       ,intent(in)     :: integ

    end subroutine setIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setManipulator(this,manip)
        !! Setter for manipulator object

        class(Modeller)                          ,intent(inout)  :: this
        class(CompositeManipulator)              ,intent(in)     :: manip

    end subroutine setManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateModelTermGroup(this,modelIndex,groupIndex,varCont)
        !! Call the update routine of model with given index for the given term group - optionally use variable container other than the
        !! one stored in the modeller

        class(Modeller)                   ,intent(inout)  :: this
        integer(ik)                       ,intent(in)     :: modelIndex
        integer(ik)                       ,intent(in)     :: groupIndex
        type(VariableContainer) ,optional ,intent(in)     :: varCont

    end subroutine updateModelTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateModelTermByName(this,modelIndex,termName,varCont)
        !! Call the update routine of model with given index for the given term name- optionally use variable container other than the
        !! one stored in the modeller

        class(Modeller)                   ,intent(inout)  :: this
        integer(ik)                       ,intent(in)     :: modelIndex
        character(*)                      ,intent(in)     :: termName
        type(VariableContainer) ,optional ,intent(in)     :: varCont

    end subroutine updateModelTermByName
!-----------------------------------------------------------------------------------------------------------------------------------
    module function evaluateModelTermGroup(this,modelIndex,groupIndex,varCont) result(res)
        !! Call the evaluateTermGroup routine on model with given index and for given term group - optionally use variable container other
        !! than the one stored in the modeller

        class(Modeller)                      ,intent(in) :: this
        integer(ik)                          ,intent(in) :: modelIndex
        integer(ik)                          ,intent(in) :: groupIndex
        type(VariableContainer)  ,optional   ,intent(in) :: varCont
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateModelTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
    module function evaluateModelTermByName(this,modelIndex,name,varCont) result(res)
        !! Call the evaluateTermByName routine on model with given index and for given term name - optionally use variable container other
        !! than the one stored in the modeller

        class(Modeller)                      ,intent(in) :: this
        integer(ik)                          ,intent(in) :: modelIndex
        character(*)                         ,intent(in) :: name
        type(VariableContainer)  ,optional   ,intent(in) :: varCont
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateModelTermByName
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine calculateMatGroupValsInModel(this,modelIndex,groupIndex,varCont)
        !! Calculate matrix value in implicit term group given by groupIndex in model given by modelIndex, and optionally using variable
        !! container other than the one stored in the modeller

        class(Modeller)                   ,intent(inout)  :: this
        integer(ik)                       ,intent(in)     :: modelIndex
        integer(ik)                       ,intent(in)     :: groupIndex
        type(VariableContainer) ,optional ,intent(in)     :: varCont

    end subroutine calculateMatGroupValsInModel
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine calculateMatValsByTermName(this,modelIndex,termName,varCont)
        !! Calculate matrix value in implicit term in model given by modelIndex, and optionally using variable
        !! container other than the one stored in the modeller

        class(Modeller)                   ,intent(inout)  :: this
        integer(ik)                       ,intent(in)     :: modelIndex
        character(*)                      ,intent(in)     :: termName
        type(VariableContainer) ,optional ,intent(in)     :: varCont

    end subroutine calculateMatValsByTermName
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addModelMatGroupToPETSc(this,modelIndex,groupIndex,mult,petscGroup)
        !! Send off matrix values of given term group and model to the PETSc controller, multiplied by mult
    
        class(Modeller)         ,intent(inout)    :: this
        integer(ik)             ,intent(in)       :: modelIndex
        integer(ik)             ,intent(in)       :: groupIndex
        real(rk)                ,intent(in)       :: mult
        integer(ik) ,optional   ,intent(in)       :: petscGroup

    end subroutine addModelMatGroupToPETSc
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine assemble(this,withIdentityMat)
        !! Assemble all of the matrix terms, preallocate PETScController objects and make sure modeller is ready for use. 
        !! If withIdentityMat is true preallocates diagonal elements for the identity matrix.

        class(Modeller)          ,intent(inout) :: this
        logical ,optional        ,intent(in)    :: withIdentityMat

    end subroutine assemble
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine performComm(this,commData,varCont,onlyDepth)
        !! Perform communications using the modeller's MPI controller - optionally uses external CommunicationData instructions or 
        !! exchanges data in external variable container
    
        class(Modeller)                   ,intent(inout) :: this
        type(CommunicationData) ,optional ,intent(in)    :: commData !! Optional non-default communication data 
        type(VariableContainer) ,optional ,intent(inout) :: varCont !! Optional variable container to perform communications on instead of the modeller's 
        integer(ik)             ,optional ,intent(in)    :: onlyDepth !! Only communicate variables at given derivation depth

    end subroutine performComm
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine safeCommAndDeriv(this,commData,varCont,derivPriority)
        !! Perform communications interlaced with derivation calls of increasing derivation depth
    
        class(Modeller)                   ,intent(inout) :: this
        type(CommunicationData) ,optional ,intent(in)    :: commData !! Optional non-default communication data 
        type(VariableContainer) ,optional ,intent(inout) :: varCont !! Optional variable container to perform communications on instead of the modeller's 
        integer(ik)             ,optional ,intent(in)    :: derivPriority !! Derivation priority for interlaced calls

    end subroutine safeCommAndDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine callManipulator(this,priority,inVars,outVars)
        !! Call optional CompositeManipulator manipulate routine with inVars and outVars - the default for the optional VariableContainers is the modeller's VariableContainer

        class(Modeller)                   ,intent(inout) :: this
        integer(ik)                       ,intent(in)    :: priority
        type(VariableContainer) ,optional ,intent(in)    :: inVars
        type(VariableContainer) ,optional ,intent(inout) :: outVars

    end subroutine callManipulator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrate(this,inVars,outVars)
        !! Call ianipulator affect routine with inVars and outVars - the default for the optional VariableContainers is the modeller's VariableContainer

        class(Modeller)                   ,intent(inout) :: this
        type(VariableContainer) ,optional ,intent(in)    :: inVars
        type(VariableContainer) ,optional ,intent(inout) :: outVars

    end subroutine integrate
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isModelTermGroupMixed(this,groupIndex,modelIndex) result(mixed)
    !! Return true if group with given index is mixed in model with given index

        class(Modeller)  ,intent(in) :: this
        integer(ik)      ,intent(in) :: groupIndex
        integer(ik)      ,intent(in) :: modelIndex
        logical                      :: mixed

    end function isModelTermGroupMixed
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getEvolvedVarInTermGroup(this,groupIndex,modelIndex) result(name)
        !! Return variable name for given group index and model index
 
        class(Modeller)           ,intent(in) :: this
        integer(ik)               ,intent(in) :: groupIndex
        integer(ik)               ,intent(in) :: modelIndex
        character(:) ,allocatable             :: name

    end function getEvolvedVarInTermGroup
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCurrentTime(this) result(time)
        !! Return current value of the integrator time variable if the integrator object is a composite integrator or the value of the 
        !! "time" variable if it isn't
 
        class(Modeller)           ,intent(in) :: this
        real(rk)                              :: time

    end function getCurrentTime
!-----------------------------------------------------------------------------------------------------------------------------------
    module function isTrueEverywhere (this,input) result(isTrue)
        !! Return true if input is true on every processor using the MPI controller of the modeller

        class(Modeller) ,intent(inout) :: this
        logical         ,intent(inout) :: input
        logical                        :: isTrue

    end function isTrueEverywhere
!-----------------------------------------------------------------------------------------------------------------------------------
    module function getGlobalMin (this,input) result(min)
        !! Return minimum value of real input computed on all processes using the MPI controller of the modeller

        class(Modeller)        ,intent(inout) :: this
        real(rk)               ,intent(inout) :: input
        real(rk)                              :: min

    end function getGlobalMin
!-----------------------------------------------------------------------------------------------------------------------------------
    module function getGlobalMax (this,input) result(max)
        !! Return maximum value of real input computed on all processes using the MPI controller of the modeller

        class(Modeller)        ,intent(inout) :: this
        real(rk)               ,intent(inout) :: input
        real(rk)                              :: max

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

    end subroutine linearSolvePETSc
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calculateIdentityMat(this,indexingObj)
        !!  Initializes the identity matrix of this modeller based on Indexing object

        class(Modeller)          ,intent(inout)  :: this
        type(Indexing)           ,intent(in)     :: indexingObj

    end subroutine calculateIdentityMat
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateModelData(this,modelIndex,varCont,updatePriority)
        !! Update modelbound data of model with given index if that data is allocated

        class(Modeller)         ,intent(inout)  :: this
        integer(ik)             ,intent(in)     :: modelIndex !! Model index for model whose data is to be update
        type(VariableContainer) ,intent(in)     :: varCont !! Variable container to be used in update
        integer(ik) ,optional   ,intent(in)     :: updatePriority !! Priority for this update call 

    end subroutine updateModelData
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine copyDataFromModel(this,modelIndex,varName,container) 
        !! Copy modelbound variable data with given name from model with given index.

        class(Modeller)                       ,intent(in)    :: this 
        integer(ik)                           ,intent(in)    :: modelIndex
        character(*)                          ,intent(in)    :: varName !! Name of data
        real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into 

    end subroutine copyDataFromModel
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module modeller_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
