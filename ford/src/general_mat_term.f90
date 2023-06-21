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
module general_mat_term_class
    !! author: Stefan Mijin 
    !!
    !! Houses general matrix term

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray ,IntArray ,StringArray
    use variable_container_class               ,only: VariableContainer
    use matrix_term_abstract_class             ,only: MatrixTerm ,MatrixTermIndexingData ,DataCoords
    use partition_class                        ,only: Partition
    use grid_class                             ,only: Grid
    use stencil_class                          ,only: Stencil
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use indexing_class                         ,only: Indexing
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate
    use signal_abstract_class                  ,only: TimeSignal ,TimeSignalData
    use model_class                            ,only: Model
    use support_functions                      ,only: removeDupeInts ,allCombinations
    use sparse_row_data_class                  
    
    implicit none
    private

    type ,public :: VarData 
        !! Container for common data used to initialize general matrix terms 
        type(StringArray) ,allocatable ,dimension(:) ,public :: rowVars !! Names of required row variables
        type(StringArray) ,allocatable ,dimension(:) ,public :: colVars !! Names of required column variables
        type(StringArray) ,allocatable ,dimension(:) ,public :: modelboundRowVars !! Named of required modelbound row variables
        type(StringArray) ,allocatable ,dimension(:) ,public :: modelboundColVars !! Named of required modelbound column variables

        real(rk)          ,allocatable ,dimension(:) ,public :: rowVarPowers !! Powers corresponding to row variables
        real(rk)          ,allocatable ,dimension(:) ,public :: colVarPowers !! Powers corresponding to column variables
        real(rk)          ,allocatable ,dimension(:) ,public :: modelboundRowVarPowers !! Powers corresponding to modelbound row variables
        real(rk)          ,allocatable ,dimension(:) ,public :: modelboundColVarPowers !! Powers corresponding to modelbound column variables

    end type VarData

    type ,public :: StencilTemplate 

        class(JaggedArrayGenerator)  ,allocatable   ,public :: stencilGen 
        integer(ik) ,allocatable ,dimension(:,:)    ,public :: rowCoords
        type(Stencil)                               ,public :: defaultStencil
        integer(ik)    ,allocatable ,dimension(:,:) ,public :: overridingStencilCoords
        type(Stencil) ,allocatable ,dimension(:)    ,public :: overridingStencils
        logical                                     ,public :: fixedStencil = .false.

    end type

    type ,public :: CoordProfiles
        !! Container for profiles in the x,h, and v directions 

        real(rk) ,allocatable ,dimension(:) ,public :: xProfile !! Profile in x direction (should be same size as grid)
        real(rk) ,allocatable ,dimension(:) ,public :: hProfile !! Profile in harmonic direction (should conform to number of harmonics)
        real(rk) ,allocatable ,dimension(:) ,public :: vProfile !! Profile in v direction (should be same saize as v grid)

    end type CoordProfiles

    type ,public ,extends(MatrixTerm) :: GeneralMatrixTerm
    !! General matrix term with the capability to take in various stencils, time and spaticoordinate profiles, and evaluate other term groups in the host model.
    !! Required row and column variables passed through the VarData object are raised to corresponding powers.

        type(VarData)                       ,private :: vData !! Common variable data

        type(CoordProfiles)                 ,private :: coordProfile !! Optional multiplicative constants giving x-dependence if size 1 or x,h,v dependence if size 3

        type(TimeSignalData) ,allocatable   ,private :: tData !! Optional time signal data
        real(rk)                            ,private :: timeSignalMult !! Multiplier containing tSignal result. Defaults to 1

        integer(ik) ,allocatable            ,private :: targetTermGroup !! Term group to evaluate and use as a row function
        real(rk) ,allocatable ,dimension(:) ,private :: evalBuffer !! Buffer for evaluating target term group

        class(JaggedArrayGenerator) ,allocatable ,private :: stencilGen !! Optional stencil value generator used to update this term's multConst

        type(RealArray) ,allocatable ,dimension(:) ,private :: stencilVals !! Current values of stencil variables, used as this term's multConst

        type(SparseRowData)                 ,private :: copyTermVals !! Matrix values of the optional copy term

        character(:) ,allocatable           ,private :: copyTermName !! Name of the term whose values should be copied and multiplied with this term's stencil

        real(rk) ,allocatable ,dimension(:) ,private :: modelboundRowVarBuffer !! Buffer holding the product of modelbound row variables raised to their corresponding powers 
        
        real(rk) ,allocatable ,dimension(:) ,private :: modelboundColVarBuffer !! Buffer holding the product of modelbound col variables raised to their corresponding powers 

        logical                             ,private :: xPeriodic !! True if term treats x grid as periodic

        logical                             ,private :: isActive 

        logical                             ,private :: kineticRow !! True if row variable is a distribution
        logical                             ,private :: kineticCol !! True if col variable is a distribution

        logical ,allocatable ,dimension(:)  ,private :: reqRowVarIsDist !! True for all required row variables that are distributions
        logical ,allocatable ,dimension(:)  ,private :: reqColVarIsDist !! True for all required column variables that are distributions

        integer(ik) ,private :: locNumX !! Local number of spatial cells - used when handling modelbound data 
        integer(ik) ,private :: numH !! Local number of spatial cells - used when handling modelbound data 
        integer(ik) ,private :: numV !! Local number of spatial cells - used when handling modelbound data 

        integer(ik) ,allocatable ,dimension(:) ,private :: rowCoordsX !x coordinate associated with each matrix row
        integer(ik) ,allocatable ,dimension(:) ,private :: rowCoordsH !h coordinate associated with each matrix row
        integer(ik) ,allocatable ,dimension(:) ,private :: rowCoordsV !v coordinate associated with each matrix row

        contains

        procedure ,public :: matrixTermUpdate => updateGeneralTerm

        procedure ,public :: calculateRows => generalRow
        procedure ,public :: calculateCols => generalCol
        procedure ,public :: init => initGeneralTerm

    end type GeneralMatrixTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initGeneralTerm(this,gridObj,partitionObj,indexingObj,procRank,&
        evolvedVar,implicitVar,varCont,sTemplate,vData,normConst,&
        coordProfile,tData,evalTermGroup,mbData,copyTermName) 
        !! General matrix term initialization routine 
    
        class(GeneralMatrixTerm)                    ,intent(inout)  :: this
        type(Grid)                                  ,intent(in)     :: gridObj !! Grid object used to initialize stencil
        type(Partition)                             ,intent(in)     :: partitionObj !! Parition object used to determine local number of DoF
        type(Indexing)                              ,intent(in)     :: indexingObj !! Indexing object used to initialize matrix row data
        integer(ik)                                 ,intent(in)     :: procRank !! Current processor rank
        character(*)                                ,intent(in)     :: evolvedVar !! Name of evolved variable
        character(*)                                ,intent(in)     :: implicitVar !! Name of implicit variable
        type(VariableContainer)                     ,intent(in)     :: varCont !! Reference variable container
        type(StencilTemplate)                       ,intent(inout)  :: sTemplate !! StencilTemplate used 
        type(VarData)       ,optional               ,intent(in)     :: vData !! Required variable names and powers
        real(rk)            ,optional               ,intent(in)     :: normConst  !! Normalization constant 
        type(CoordProfiles) ,optional               ,intent(in)     :: coordProfile !! Multiplicative coordinate profiles (globally indexed)
        type(TimeSignalData) ,optional              ,intent(in)     :: tData !! Explicit time dependence 
        integer(ik)       ,optional                 ,intent(in)     :: evalTermGroup !! Optional term group of host model to evaluate and use as additional row variable
        class(ModelboundData) ,optional             ,intent(in)     :: mbData !! Optional modelbound data object used when initializing fixed stancil generators which might require it
        character(*) ,optional                      ,intent(in)     :: copyTermName !! Name of the term whose values should be copied and multiplied with this term's stencil

    end subroutine initGeneralTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine generalRow(this,varCont,rowVals,indexingData) 
        !! Multiplies all required variables and raises them to their corresponding power,using modelbound row variable buffer if allocated, and adding 
        !! spatial and time profiles

        class(GeneralMatrixTerm)      ,intent(inout)   :: this
        type(VariableContainer)       ,intent(in)      :: varCont
        real(rk) ,dimension(:)        ,intent(inout)   :: rowVals
        type(MatrixTermIndexingData)  ,intent(in)      :: indexingData

    end subroutine generalRow 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine generalCol(this,varCont,colVals,indexingData) 
        !! Multiplies all required variables and raises them to their corresponding power,using modelbound row variable buffer if allocated

        class(GeneralMatrixTerm)      ,intent(inout) :: this
        type(VariableContainer)       ,intent(in)    :: varCont
        type(RealArray) ,dimension(:) ,intent(inout) :: colVals
        type(MatrixTermIndexingData)  ,intent(in)    :: indexingData

    end subroutine generalCol 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine updateGeneralTerm(this,varCont,indexingData,modelData,hostModel) 
        !! Update general term if using any modelbound data, updating stencil values, evaluating other terms, or if there is any explicit time dependence 

        class(GeneralMatrixTerm)        ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container - used to update stencil values and signal
        type(MatrixTermIndexingData)    ,intent(in)     :: indexingData
        class(ModelboundData) ,optional ,intent(in)     :: modelData !! Model data used to retrieve modelbound variable
        class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - used to evaluate term group

    end subroutine updateGeneralTerm
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module general_mat_term_class
!-----------------------------------------------------------------------------------------------------------------------------------
 