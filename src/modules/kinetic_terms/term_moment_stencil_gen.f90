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
module term_moment_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for moment terms created using fully kinetic terms

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray
    use variable_container_class               ,only: VariableContainer
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use partition_class                        ,only: Partition
    use v_space_class                          ,only: VSpace
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate
    use matrix_term_abstract_class             ,only: MatrixTermIndexingData, DataCoords
    use sparse_row_data_class                  ,only: SparseRowData
    use model_class                            ,only: Model
    use physical_constants

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: TermMomentStencilGenerator
        !! JaggedArrayGenerator for calculating moment terms based on taking velocity moments of other kinetic terms. The term whose
        !! moment is taken is expected to be evolving a single harmonic of a distribution variable and have a diagonal spatial stencil.
        !! If the implcit variable is a distribution, the spatial stencil must be strictly diagonal (no interpolation allowed).

        integer(ik) ,private :: locNumX !! Local x grid size (including any staggered grid modifications)
        integer(ik) ,private :: numV !! Copy of v grid size
        real(rk) ,allocatable ,dimension(:) ,private :: vVec !! Vector used to calculate the moment (4*pi*v**(2+m)*dv)

        character(:) ,allocatable ,private :: termName !! Name of the term whose moment is to be taken 
        
        type(SparseRowData) ,private :: matBuffer !! Buffer for the matrix values of used term
        type(MatrixTermIndexingData) ,allocatable ,private :: indexingDataBuffer !! Buffer for matrix term indexing data used to take moment
        
        contains

        procedure ,public :: init => initTermMomentGen

        procedure ,public :: calculateInPlace => calcTermMomentVals

    end type TermMomentStencilGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initTermMomentGen(this,partitionObj,vspaceObj,procRank,momentOrder,termName,removeLastCell) 
        !! Term moment stencil value generator initialization routine
    
        class(TermMomentStencilGenerator)         ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        type(VSpace)                              ,intent(in)     :: vspaceObj !! VSpace object used to get interpolation object
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        integer(ik)                               ,intent(in)     :: momentOrder !! Order of moment to be taken
        character(*)                              ,intent(in)     :: termName !! Name of term in host model whose moment should be taken
        logical ,optional                         ,intent(in)     :: removeLastCell !! Set to true if the row variable is staggered and the grid is not periodic. Defaults to false.

    end subroutine initTermMomentGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcTermMomentVals(this,varCont,res,mbData,hostModel)
        !! Calculate term moment stencil values in place. Requires hostModel.

        class(TermMomentStencilGenerator)           ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcTermMomentVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module term_moment_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 