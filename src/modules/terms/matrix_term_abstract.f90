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
module matrix_term_abstract_class
    !! author: Stefan Mijin 
    !!
    !! Houses abstract matrix term with default interface

    use data_kinds                     ,only: rk, ik
    use runtime_constants              ,only: debugging, assertions
    use god_objects                    ,only: Object
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use support_types                  ,only: RealArray ,IntArray ,StringArray
    use variable_list_class            ,only: VariableList
    use variable_container_class       ,only: VariableContainer
    use sparse_row_data_class          ,only: SparseRowData
    use petsc_preallocation_data_class ,only: PETScPreallocationData
    use petsc_controller_class         ,only: PETScController
    use basic_interfaces               ,only: coordMapping
    use indexing_class                 ,only: Indexing
    use support_functions              ,only: allCombinations
    use term_abstract_class            ,only: Term
    use modelbound_data_abstract_class ,only: ModelboundData
    use model_surrogate_class          ,only: ModelSurrogate

    implicit none
    private

    type ,public :: DataCoords 
        !! Data type used to store global coordinates of any given matrix row

        integer(ik)    ,allocatable ,dimension(:)   :: rowCoords !! Size 1 or 3 array depending on whether the evolved variable is a distribution 
        integer(ik)    ,allocatable ,dimension(:,:) :: colCoords !! Size (1,:) or (3,:) depending on whether the implicit variable is a distribution
                                                                 !! - the second dimension has the size of the column vector corresponding to given row of matrix

    end type DataCoords

    type ,public :: MatrixTermIndexingData
        !! Indexing data used by the matrix term

        type(DataCoords)    ,allocatable ,dimension(:) :: rowDataCoordsGlobal !! Global coordinates of each row in matrix
        type(DataCoords)    ,allocatable ,dimension(:) :: rowDataCoordsLocal !! Local coordinates of each row in matrix (unflattened)
        integer(ik)         ,allocatable ,dimension(:) :: localRowIndices !! Local row indices corresponding to stored evolved variable (flattened)
        type(IntArray)      ,allocatable ,dimension(:) :: localColIndices !! Local col indices corresponding to stored implicit variable (flattened)
        integer(ik)         ,allocatable ,dimension(:) :: rowReqVarIndices !! Indices of variables entering the matrix as a function of the row coordinates
        integer(ik)         ,allocatable ,dimension(:) :: colReqVarIndices !! Indices of variables entering the matrix as s function of the column coordinates
        character(:)        ,allocatable               :: rowVarName !! Name of the evolved/row variable
        character(:)        ,allocatable               :: colVarName !! Name of the implicit/column variable
        
    end type

    type ,public ,extends(term) ,abstract :: MatrixTerm
        !! Abstract matrix term containing general matrix construction routines and interface with PETSc. 
        !! The values of a matrix are by default constructed by evaluating row and column functions for each element, and multiplying with the corresponding 
        !! multiplicative and normalization constants 

        type(SparseRowData)                            ,private :: rowData  !! Matrix object in sparse row data form
        
        type(RealArray)     ,allocatable ,dimension(:) ,private :: multConst !! Multiplicative constants for each row/column corresponding to the matrix structure
        real(rk)                                       ,private :: normalizationConst = real(1.0d0,kind=rk) !! Normalization constant for this term

        type(MatrixTermIndexingData)                   ,private :: indexingData !! Indexing data used to determing coordinates and required variables

        logical                                        ,private :: fixedMatrix = .false. !! If true the matrix is constructed only once
        logical                                        ,private :: fixedMatrixCalculated = .false. !! True if the fixed matrix has been calculated
        logical                                        ,private :: stationaryEvolvedVar = .false. !! If true, this term is always sent to PETSc with multiplicative constant -1

        logical                                        ,private :: hasNonTrivialRowFun = .true. !! Set to false to not call row function. Defaults to true.
        logical                                        ,private :: hasNonTrivialColFun = .true. !! Set to false to not call col function. Defaults to true.

        type(RealArray)     ,allocatable ,dimension(:) ,private :: colBuffer !Buffer for column vector values
        real(rk)            ,allocatable ,dimension(:) ,private :: rowBuffer !Buffer for row vector values

        contains

        procedure ,public  :: calculateValues 
        procedure ,public  :: setNormalizationConst
        procedure ,public  :: getNormalizationConst
        procedure ,public  :: setReqVars
        procedure ,public  :: setEvolvedAndImplicitVar
        procedure ,public  :: initRowData
        procedure ,public  :: evaluate => evaluateMatTerm

        procedure ,public  :: getMultConst
        procedure ,public  :: setMultConst

        procedure ,public  :: getIndexingData

        procedure ,public  :: setFixedMatrix
        procedure ,public  :: setNontrivialRowFun
        procedure ,public  :: setNontrivialColFun

        procedure ,private :: calculateRows => unityRow
        procedure ,private :: calculateCols => unityCol

        procedure ,public  :: addRowDataToPreallocationData
        procedure ,public  :: addRowDataPatternToController
        procedure ,public  :: addRowValuesToPETScMatrix

        procedure ,public  :: getVarName => getRowVarName

        procedure ,public  :: getRowData
        procedure ,public  :: update => matTermUpdate
        procedure ,public  :: matrixTermUpdate => matrixNoUpdate

    end type MatrixTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine unityRow(this,varCont,rowVals,indexingData) 
        !! Default row function. In general should use the passed variable container, row buffer, and matrix indexing data
        !! to provide values which is only a function of the row 

        class(MatrixTerm)             ,intent(inout)   :: this
        type(VariableContainer)       ,intent(in)      :: varCont
        real(rk) ,dimension(:)        ,intent(inout)   :: rowVals
        type(MatrixTermIndexingData)  ,intent(in)      :: indexingData

    end subroutine unityRow 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine unityCol(this,varCont,colVals,indexingData) 
        !! Default column function. In general should provide values for each column 

        class(MatrixTerm)             ,intent(inout) :: this
        type(VariableContainer)       ,intent(in)    :: varCont
        type(RealArray) ,dimension(:) ,intent(inout) :: colVals
        type(MatrixTermIndexingData)  ,intent(in)    :: indexingData
        
    end subroutine unityCol  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine calculateValues(this,varCont) 
        !! Default matrix entry calculation routine - loops over each row and applies rowFun and colFun to variables in provided container

        class(MatrixTerm)       ,intent(inout)  :: this
        type(VariableContainer) ,intent(in)     :: varCont

    end subroutine calculateValues
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNormalizationConst(this,norm) 
        !! Setter for normalizationConst

        class(MatrixTerm)       ,intent(inout)  :: this
        real(rk)                ,intent(in)     :: norm

    end subroutine setNormalizationConst
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNormalizationConst (this) result(norm)
        !! Getter for normalizationConst

        class(MatrixTerm)  ,intent(in) :: this
        real(rk)                       :: norm

    end function getNormalizationConst
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setReqVars(this,rowReqVars,colReqVars,varCont) 
        !! Set variable names required by the row and column functions and find their indices in variable container

        class(MatrixTerm)               ,intent(inout)  :: this
        type(StringArray) ,dimension(:) ,intent(in)     :: rowReqVars
        type(StringArray) ,dimension(:) ,intent(in)     :: colReqVars
        type(VariableContainer)         ,intent(in)     :: varCont

    end subroutine setReqVars
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setEvolvedAndImplicitVar(this,rowVarName,colVarName,varCont) 
        !! Set evolved (row) and implicit (column) variable names, and check if evolved variable is stationary using varCont

        class(MatrixTerm)         ,intent(inout)  :: this
        character(*)              ,intent(in)     :: rowVarName
        character(*)              ,intent(in)     :: colVarName
        type(VariableContainer)   ,intent(in)     :: varCont

    end subroutine setEvolvedAndImplicitVar
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initRowData(this,rowCoords,rowToColMapping,indexingObj) 
        !! Initialize row data objects based on a set of evolved global coordinates, the row/col var name, and a function that returns 
        !! column coordinates given a row coordinate input. Requires a reference Indexing object.

        class(MatrixTerm)            ,intent(inout)  :: this
        integer(ik)  ,dimension(:,:) ,intent(in)     :: rowCoords
        procedure(coordMapping)                      :: rowToColMapping
        type(Indexing)               ,intent(in)     :: indexingObj

    end subroutine initRowData
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addRowDataToPreallocationData(this,petscPreallocData) 
        !! Add this term's row data to a PETSc preallocation object

        class(MatrixTerm)            ,intent(inout)  :: this
        type(PETScPreallocationData) ,intent(inout)  :: petscPreallocData

    end subroutine addRowDataToPreallocationData
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addRowDataPatternToController(this,petscCont) 
    !! Add this term's row data to a petsc preallocation object in PETSc controller

        class(MatrixTerm)     ,intent(in)     :: this
        type(PETScController) ,intent(inout)  :: petscCont

    end subroutine addRowDataPatternToController
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addRowValuesToPETScMatrix(this,petscCont,mult,petscGroup) 
        !! Add this term's row values to a petsc matrix object in PETSc controller, multiplied by mult. If evolved variable is
        !! stationary mult is ignored and set to -1
    
        class(MatrixTerm)     ,intent(in)     :: this
        type(PETScController) ,intent(inout)  :: petscCont
        real(rk)              ,intent(in)     :: mult
        integer(ik) ,optional ,intent(in)     :: petscGroup

    end subroutine addRowValuesToPETScMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
    module function evaluateMatTerm (this,varCont) result(res)
        !! Get explicit value for the term by crudely dotting each sparse row with the locally stored implicit variable vector in varCont

        class(MatrixTerm)                    ,intent(in) :: this
        type(VariableContainer)              ,intent(in) :: varCont
        real(rk) ,allocatable ,dimension(:)              :: res

    end function evaluateMatTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setMultConst(this,multConst) 
        !! Setter for multConst

        class(MatrixTerm)              ,intent(inout)  :: this
        type(RealArray) ,dimension(:)  ,intent(in)     :: multConst

    end subroutine setMultConst
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getMultConst (this) result(multConst)
        !! Getter for multConst

        class(MatrixTerm)                          ,intent(in) :: this
        type(RealArray) ,allocatable ,dimension(:)             :: multConst

    end function getMultConst
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getRowVarName(this) result(name)
        !! Get name of the evolved variable of this term

        class(MatrixTerm)                    ,intent(in) :: this
        character(:) ,allocatable                        :: name

    end function getRowVarName
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setFixedMatrix(this,isFixed) 
        !! Set whether the matrix is fixed -> only calculated once

        class(MatrixTerm)              ,intent(inout)  :: this
        logical                        ,intent(in)     :: isFixed

    end subroutine setFixedMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNonTrivialRowFun(this,nontriv) 
        !! Set whether the matrix has a non-trivial row function

        class(MatrixTerm)              ,intent(inout)  :: this
        logical                        ,intent(in)     :: nontriv

    end subroutine setNonTrivialRowFun
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNonTrivialColFun(this,nontriv) 
        !! Set whether the matrix has a non-trivial col function

        class(MatrixTerm)              ,intent(inout)  :: this
        logical                        ,intent(in)     :: nontriv

    end subroutine setNonTrivialColFun
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getIndexingData(this) result(indData)
        !! Getter for indexingData

        class(MatrixTerm)                    ,intent(in) :: this
        type(MatrixTermIndexingData)                     :: indData

    end function getIndexingData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getRowData(this) result(rowData)
        !! Getter for rowData

        class(MatrixTerm)           ,intent(in) :: this
        type(SparseRowData)                     :: rowData

    end function 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine matTermUpdate(this,varCont,modelData,hostModel) 
        !! Default matrix term update, call matrixTermUpdate 

        class(MatrixTerm)               ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont
        class(ModelboundData) ,optional ,intent(in)     :: modelData
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel

    end subroutine matTermUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine matrixNoUpdate(this,varCont,indexingData,modelData,hostModel) 
        !! Default MatrixTerm updateMatTerm function - does nothing

        class(MatrixTerm)               ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont
        type(MatrixTermIndexingData)    ,intent(in)     :: indexingData
        class(ModelboundData) ,optional ,intent(in)     :: modelData
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel

    end subroutine matrixNoUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module matrix_term_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 