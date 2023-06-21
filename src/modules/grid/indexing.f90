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
module indexing_class
    !! author: Stefan Mijin 
    !! 
    !! Houses Indexing object responsible for handling the indexing of the global variable vector

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use grid_class                  ,only: Grid
    use variable_list_class         ,only: VariableList
    use Partition_class             ,only: Partition 
    use support_types               ,only: IntArray

    implicit none
    private

    type ,public ,extends(Object) :: Indexing
        !! Indexing object containing local and global indexing rules 

        integer(ik) ,allocatable ,dimension(:)   ,private :: procDoFs !! Number of degrees of freedom for each processor
        integer(ik) ,allocatable ,dimension(:,:) ,private :: varOffset !! For indices (i,k) gives offset of variable k in local spatial chunk of processor i
        integer(ik) ,allocatable ,dimension(:)   ,private :: xDoFs !! Degrees of freedom in the x-direction for each processor

        integer(ik)                              ,private :: numV !! Number of velocity cells in grid
        integer(ik)                              ,private :: numX !! Total number of spatial cells in grid
        integer(ik)                              ,private :: numH !! Total numbe of harmonics in grid
 
        type(Partition)                          ,private :: part !! Associated partition object
        type(VariableList)                       ,private :: varList !! Associated implicit variable list 

        integer(ik)                              ,private :: xHaloWidth       !! Halo width in x-direction - only used for periodic indexing
        logical                                  ,private :: periodicIndexing !! Used to make sure that local indexing is correct 
                                                                              !! in periodic case

        contains

        procedure ,public :: getProcDoF 
        procedure ,public :: findIndex
        procedure ,public :: getNumX
        procedure ,public :: getNumV
        procedure ,public :: getNumH
        procedure ,public :: getAllIndicesOfVar

        procedure ,public :: findDistIndex
        procedure ,public :: findLocalXIndex
        procedure ,public :: mapToGlobalIndices

        procedure ,public :: init => initIndexing

    end type Indexing
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initIndexing(this,part,gridData,varList,xHaloWidth) 
            !! Indexing object initialization routine

            class(Indexing)           ,intent(inout)  :: this
            type(Partition)           ,intent(in)     :: part !! Partition component of initialized Indexing object
            type(Grid)                ,intent(in)     :: gridData !! Grid object used to initialize Indexing
            type(VariableList)        ,intent(in)     :: varList !! Implicit variable list
            integer(ik) ,optional     ,intent(in)     :: xHaloWidth !! Width of halo in x-direction - if passed will assume the x-grid is periodic

        end subroutine initIndexing
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getProcDoF (this) result(procDoF)
            !! Getter for procDoFs

            class(Indexing)                       ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)            :: procDoF
 
        end function getProcDoF
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumH (this) result(numH)
            !! Getter for numH

            class(Indexing) ,intent(in) :: this
            integer(ik)                 :: numH
 
        end function getNumH
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumX (this) result(numX)
            !! Getter for numX

            class(Indexing) ,intent(in) :: this
            integer(ik)                 :: numX
 
        end function getNumX
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getNumV (this) result(numV)
            !! Getter for numV

            class(Indexing) ,intent(in) :: this
            integer(ik)                 :: numV
 
        end function getNumV
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function findIndex (this,name,xInd,hInd,vInd,local) result(ind)
            !! Return index of degree of freedom associated with variable with given name and with coordinates given by xInd,hInd,vInd
            !! If local is true return index in the local implicit vector

            class(Indexing)       ,intent(in) :: this
            character(*)          ,intent(in) :: name 
            integer(ik)           ,intent(in) :: xInd 
            integer(ik) ,optional ,intent(in) :: hInd
            integer(ik) ,optional ,intent(in) :: vInd
            logical     ,optional ,intent(in) :: local  
            integer(ik)                       :: ind

        end function findIndex
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function findDistIndex (this,xInd,hInd,vInd,allH,locXInd) result(ind)
            !! Return linear index of distribution corresponding to xInd,hInd,vInd - if allH is true return indexing assuming all harmonics
            !! are indexed locally. If present locXInd is used to identify the processor

            class(Indexing)       ,intent(in) :: this
            integer(ik)           ,intent(in) :: xInd
            integer(ik)           ,intent(in) :: hInd
            integer(ik)           ,intent(in) :: vInd 
            logical               ,intent(in) :: allH
            integer(ik) ,optional ,intent(in) :: locXInd
            integer(ik)                       :: ind
 
        end function findDistIndex
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function findLocalXIndex (this,xInd,locXInd) result(ind)
        !! Return local non-distribution variable index corresponding to xInd,if present, locXInd is used to identify the processor

        class(Indexing)       ,intent(in) :: this
        integer(ik)           ,intent(in) :: xInd
        integer(ik) ,optional ,intent(in) :: locXInd
        integer(ik)                       :: ind
 
        end function findLocalXIndex
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function mapToGlobalIndices (this,name,inputIndices) result(ind)
            !! Returns list of global indices associated with all possible combinations of inputIndices (IntArray of size 1 or 3) and 
            !! given variable name

            class(Indexing)                       ,intent(in) :: this
            character(*)                          ,intent(in) :: name
            type(IntArray) ,dimension(:)          ,intent(in) :: inputIndices
            integer(ik) ,allocatable ,dimension(:)            :: ind
 
        end function mapToGlobalIndices
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getAllIndicesOfVar (this,varInd,procRank) result(ind)
            !! Returns all local indices corresponding to variable with a given index on given processor
        
            class(Indexing)                       ,intent(in) :: this
            integer(ik)                           ,intent(in) :: varInd
            integer(ik)                           ,intent(in) :: procRank 
            integer(ik) ,allocatable ,dimension(:)            :: ind
 
        end function getAllIndicesOfVar
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module indexing_class
!-----------------------------------------------------------------------------------------------------------------------------------
 