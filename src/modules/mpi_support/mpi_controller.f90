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
module mpi_controller_class
    !! author: Stefan Mijin
    !! 
    !! Houses object used to interface with MPI outside of PETSc

    use mpi_f08
    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_types               ,only: IntArray ,RealArray ,StringArray
    use variable_container_class    ,only: VariableContainer 
    use partition_class             ,only: Partition


    implicit none
    private

    type ,public :: CommunicationData
        !! Contains lists of variables to broadcast in processor rows and variables participating in halo exchange
        type(StringArray) ,allocatable ,dimension(:) ,public :: varsToBroadcast  !! Variables to broadcast/exchange in processor row
        type(StringArray) ,allocatable ,dimension(:) ,public :: haloExchangeVars !! Variables that should participate in 
                                                                                 !! halo exchange
        type(StringArray) ,allocatable ,dimension(:) ,public :: scalarsToBroadcast  !! Scalar variables to broadcast everywhere
        integer(ik)       ,allocatable ,dimension(:) ,public :: scalarRoots !! Root processes for each broadcast scalar

    end type CommunicationData

    type ,public ,extends(object) :: MPIController
        !! Provides a centralized interface with MPI and various support routines compatible with ReMKiT1D needs

        type(MPI_Comm) ,private :: worldComm !! Main MPI_Comm_World
        type(MPI_Comm) ,private :: rowComm !! Communicator for processors rows (those with the same x-domain but different h-domains)
        type(MPI_Comm) ,private :: colComm !! Communicator for processor columns (those with the same h-domain but different x-domains)
        integer(ik)    ,private :: worldRank !! Current process' world rank
        integer(ik)    ,private :: worldSize !! Size of worldComm
        integer(ik)    ,private :: rowRank !! Current process' row rank
        integer(ik)    ,private :: rowSize !! Size of rowComm
        integer(ik)    ,private :: colRank !! Current process' column rank
        integer(ik)    ,private :: colSize !! Size of colComm

        integer(ik)                            ,private :: rowNumX !! Number of x-grid points in local processor row
        integer(ik) ,allocatable ,dimension(:) ,private :: colNumX !! Number of x-grid points in each column processor
        integer(ik) ,allocatable ,dimension(:) ,private :: xDispls !! Displacement in x-direction used for MPI gather
        integer(ik)                            ,private :: xHaloWidth !! Halo width in x-direction 
        integer(ik)                            ,private :: numV !! Number of velocity cells in grid
        integer(ik)                            ,private :: numH !! Number of harmonics in grid

        integer(ik)     ,allocatable ,dimension(:) :: rowNumH !! locNumH for each processor in row
        integer(ik)     ,allocatable ,dimension(:) :: rowHOffset !! harmonic offset for each processor in row
        type(RealArray) ,allocatable ,dimension(:) :: distBuffer !! Distribution variable buffer used in halo exchange

        type(IntArray)  ,allocatable ,dimension(:) :: neighbourPairs !! Int array for each x-grid neighbour pair with entry of length 2
                                                                 !! first element is left neighbour rank, second is right - allows 
                                                                 !! for periodic domain 

        logical :: rowsSetUp !! True if row data is set up and ready for communication
        contains

        procedure ,public :: getWorldRank
        procedure ,public :: getWorldSize
        procedure ,public :: getRowRank
        procedure ,public :: getRowSize
        procedure ,public :: getColRank
        procedure ,public :: getColSize
        procedure ,public :: getXHaloWidth
        procedure ,public :: getComm

        procedure ,public :: calculateRowDistData
        procedure ,public :: initializeNeighbourPairs

        procedure ,public :: broadcastVarInRow
        procedure ,public :: exchangeDistVarInRow
        procedure ,public :: exchangeVarXHalos

        procedure ,public :: broadcastReal
        procedure ,public :: broadcastInt
        procedure ,public :: broadcastLogical
        procedure ,public :: broadcastCharacter

        procedure ,public :: barrier

        procedure ,public :: gatherVar
        procedure ,public :: scatterVar

        procedure ,public :: setUpRows

        procedure ,public :: isTrueEverywhere
        procedure ,public :: allreduceMax
        procedure ,public :: allreduceMin

        procedure ,public :: init => initMPIController


    end type MPIController
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine initMPIController(this,numProcsX,numProcsH) 
            !! MPI controller initialization routine - initializes communicators used outside of PETSc. Optionally sets up row/column communicators.

            class(MPIController)      ,intent(inout)  :: this
            integer(ik)   ,optional   ,intent(in)     :: numProcsX !! Number of processes in the x-direction
            integer(ik)   ,optional   ,intent(in)     :: numProcsH !! Number of processes in the h-direction

        end subroutine initMPIController
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine setUpRows(this,numProcsX,numProcsH) 
            !! Set up row/column communication in controller

            class(MPIController)      ,intent(inout)  :: this
            integer(ik)               ,intent(in)     :: numProcsX !! Number of processes in the x-direction
            integer(ik)               ,intent(in)     :: numProcsH !! Number of processes in the h-direction

        end subroutine setUpRows
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getWorldRank (this) result(rank)
            !! Getter for worldRank

            class(MPIController) ,intent(in) :: this
            integer(ik)                      :: rank
 
        end function getWorldRank
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getWorldSize (this) result(size)
            !! Getter for worldSize

            class(MPIController) ,intent(in) :: this
            integer(ik)                      :: size
 
        end function getWorldSize
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getRowRank (this) result(rank)
            !! Getter for rowRank

            class(MPIController) ,intent(in) :: this
            integer(ik)                      :: rank
 
        end function getRowRank
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getRowSize (this) result(size)
            !! Getter for rowSize

            class(MPIController) ,intent(in) :: this
            integer(ik)                      :: size
 
        end function getRowSize
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getColRank (this) result(rank)
            !! Getter for columnRank

            class(MPIController) ,intent(in) :: this
            integer(ik)                      :: rank
 
        end function getColRank
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getColSize (this) result(size)
            !! Getter for columnSize

            class(MPIController) ,intent(in) :: this
            integer(ik)                      :: size
 
        end function getColSize
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getXHaloWidth (this) result(xHaloWidth)
            !! Getter for xHaloWidth

            class(MPIController) ,intent(in) :: this
            integer(ik)                      :: xHaloWidth
 
        end function getXHaloWidth
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getComm (this) result(mpiComm)
            !! Getter for xHaloWidth

            class(MPIController) ,intent(in) :: this
            type(MPI_Comm)                   :: mpiComm
 
        end function getComm
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine calculateRowDistData(this,partitionObj,xHaloWidth,numV) 
            !! Initialize distribution DoF data used in exchanging distribution data in a row

            class(MPIController) ,intent(inout)  :: this
            type(Partition)      ,intent(in)     :: partitionObj !! Partition object used to retrieve DoF information
            integer(ik)          ,intent(in)     :: xHaloWidth !! Halo width in x-direction
            integer(ik)          ,intent(in)     :: numV !! Number of cells in v-grid

        end subroutine calculateRowDistData
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine initializeNeighbourPairs(this,periodic) 
            !! Initialize neighbour pairs - if periodic in X adds additional pair to handle this

            class(MPIController) ,intent(inout)  :: this
            logical              ,intent(in)     :: periodic

        end subroutine initializeNeighbourPairs
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine broadcastVarInRow(this,varContainer,name) 
            !! Broadcast variable with given name from row rank 0 process to other row processes - to be used only with non-distribution 
            !! variables 

            class(MPIController)      ,intent(inout) :: this
            type(VariableContainer)   ,intent(inout) :: varContainer !! Variable container in which to broadcasty
            character(*)              ,intent(in)    :: name !! Name of variable to broadcast

        end subroutine broadcastVarInRow
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine exchangeDistVarInRow(this,varContainer,name) 
            !! Exchanges a distribution variable in a processor row 

            class(MPIController)      ,intent(inout) :: this
            type(VariableContainer)   ,intent(inout) :: varContainer !! Variable container in which to broadcast
            character(*)              ,intent(in)    :: name !! Name of variable to broadcast

        end subroutine exchangeDistVarInRow
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine exchangeVarXHalos(this,varContainer,name,varIsDist) 
            !! Exchanges a halos in x direction for given variable 

            class(MPIController)      ,intent(inout) :: this
            type(VariableContainer)   ,intent(inout) :: varContainer !! Variable container in which to perform halo exchange
            character(*)              ,intent(in)    :: name !! Name of variable to exchange
            logical                   ,intent(in)    :: varIsDist !! Set to true if variable is a distribution

        end subroutine exchangeVarXHalos
!-----------------------------------------------------------------------------------------------------------------------------------
        module function isTrueEverywhere (this,input) result(isTrue)
            !! Return true if input is true on every processor

            class(MPIController) ,intent(inout) :: this
            logical              ,intent(inout) :: input
            logical                             :: isTrue
 
        end function isTrueEverywhere
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine gatherVar(this,localVec,globVec,varIsDist) 
            !! Gather variable values into global vector for output on rank 0. If globVec in not allocated it will be allocated to the correct
            !! size, otherwise it is assumed it is of the correct size.

            class(MPIController)                   ,intent(inout) :: this
            real(rk)              ,dimension(:)    ,intent(in)    :: localVec !! Local vector to gather
            real(rk) ,allocatable ,dimension(:)    ,intent(inout) :: globVec !! Global vector on rank 0 to gather into
            logical                                ,intent(in)    :: varIsDist !! Set to true if the gathered variable is a distribution

        end subroutine gatherVar
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine scatterVar(this,globVec,localVec,varIsDist) 
            !! Scatter global variable values into local vector from rank 0. 

            class(MPIController)                   ,intent(inout) :: this
            real(rk)              ,dimension(:)    ,intent(inout) :: globVec !! Global vector on rank 0 to scatter from
            real(rk) ,allocatable ,dimension(:)    ,intent(inout) :: localVec !! Local vector to scatter into
            logical                                ,intent(in)    :: varIsDist !! Set to true if the scattered variable is a distribution

        end subroutine scatterVar
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine broadcastReal(this,vec,rootProc) 
            !! Broadcast real vector from rank rootProc

            class(MPIController)      ,intent(inout) :: this
            real(rk) ,dimension(:)    ,intent(inout) :: vec
            integer(ik) ,optional     ,intent(in)    :: rootProc

        end subroutine broadcastReal
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine broadcastInt(this,vec) 
            !! Broadcast integer vector from rank 0

            class(MPIController)      ,intent(inout) :: this
            integer(ik) ,dimension(:) ,intent(inout) :: vec

        end subroutine broadcastInt
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine broadcastLogical(this,vec) 
            !! Broadcast logical vector from rank 0

            class(MPIController)      ,intent(inout) :: this
            logical     ,dimension(:) ,intent(inout) :: vec

        end subroutine broadcastLogical
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine broadcastCharacter(this,vec) 
            !! Broadcast character vector from rank 0

            class(MPIController)      ,intent(inout) :: this
            character(:) ,allocatable ,intent(inout) :: vec

        end subroutine broadcastCharacter
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine barrier(this) 
            !! Calls MPI barrier on commWorld

            class(MPIController)      ,intent(inout) :: this

        end subroutine barrier
!-----------------------------------------------------------------------------------------------------------------------------------
        module function allreduceMin (this,input) result(min)
            !! Return min value on all processors
        
            class(MPIController)   ,intent(inout) :: this
            real(rk)               ,intent(inout) :: input
            real(rk)                              :: min
        
        end function allreduceMin
!-----------------------------------------------------------------------------------------------------------------------------------
        module function allreduceMax (this,input) result(max)
            !! Return max value on all processors
        
            class(MPIController)   ,intent(inout) :: this
            real(rk)               ,intent(inout) :: input
            real(rk)                              :: max
        
        end function allreduceMax
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module mpi_controller_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
