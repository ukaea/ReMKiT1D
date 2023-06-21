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
submodule (mpi_controller_class) mpi_controller_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the mpi controller class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMPIController(this,numProcsX,numProcsH) 
    !! MPI controller initialization routine - initializes communicators used outside of PETSc. Optionally sets up row/column communicators.

    class(MPIController)      ,intent(inout)  :: this
    integer(ik)   ,optional   ,intent(in)     :: numProcsX !! Number of processes in the x-direction
    integer(ik)   ,optional   ,intent(in)     :: numProcsH !! Number of processes in the h-direction

    logical :: mpiInitialized

    call MPI_Initialized(mpiInitialized)
    if (.not. mpiInitialized) call MPI_init()
    call MPI_Comm_dup(MPI_Comm_World,this%worldComm)
    call MPI_Comm_size(this%worldComm,this%worldSize)
    call MPI_Comm_rank(this%worldComm,this%worldRank)

    this%rowsSetUp = .false.

    if (present(numProcsX) .and. present(numProcsH)) call this%setUpRows(numProcsX,numProcsH)

    call this%makeDefined()

end subroutine initMPIController 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getWorldRank (this) result(rank)
    !! Getter for worldRank

    class(MPIController) ,intent(in) :: this
    integer(ik)                      :: rank

    if (assertions) call assertPure(this%isDefined(),"getWorldRank called on undefined mpi controller")

    rank = this%worldRank

end function getWorldRank
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getWorldSize (this) result(size)
    !! Getter for worldSize

    class(MPIController) ,intent(in) :: this
    integer(ik)                      :: size

    if (assertions) call assertPure(this%isDefined(),"getWorldSize called on undefined mpi controller")

    size = this%worldSize

end function getWorldSize
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRowRank (this) result(rank)
    !! Getter for rowRank

    class(MPIController) ,intent(in) :: this
    integer(ik)                      :: rank

    if (assertions) then 
        call assertPure(this%isDefined(),"getRowRank called on undefined mpi controller")
        call assertPure(this%rowsSetUp,"getRowRank called from MPI controller before rows set up")
    end if

    rank = this%rowRank

end function getRowRank
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRowSize (this) result(size)
    !! Getter for rowSize

    class(MPIController) ,intent(in) :: this
    integer(ik)                      :: size

    if (assertions) then 
        call assertPure(this%isDefined(),"getRowSize called on undefined mpi controller")
        call assertPure(this%rowsSetUp,"getRowSize called from MPI controller before rows set up")
    end if

    size = this%rowSize

end function getRowSize
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getColRank (this) result(rank)
    !! Getter for columnRank

    class(MPIController) ,intent(in) :: this
    integer(ik)                      :: rank


    if (assertions) then 
        call assertPure(this%isDefined(),"getColRank called on undefined mpi controller")
        call assertPure(this%rowsSetUp,"getColRank called from MPI controller before rows set up")
    end if

    rank = this%colRank

end function getColRank
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getColSize (this) result(size)
    !! Getter for columnSize    

    class(MPIController) ,intent(in) :: this
    integer(ik)                      :: size

    if (assertions) then 
        call assertPure(this%isDefined(),"getColSize called on undefined mpi controller")
        call assertPure(this%rowsSetUp,"getColSize called from MPI controller before rows set up")
    end if

    size = this%colSize

end function getColSize
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getXHaloWidth (this) result(xHaloWidth)
    !! Getter for xHaloWidth

    class(MPIController) ,intent(in) :: this
    integer(ik)                      :: xHaloWidth

    if (assertions) then 
        call assertPure(this%isDefined(),"getXHaloWidth called on undefined mpi controller")
    end if

    xHaloWidth = this%xHaloWidth

end function getXHaloWidth
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine broadcastVarInRow(this,varContainer,name) 
    !! Broadcast variable with given name from row rank 0 process to other row processes - to be used only with non-distribution 
    !! variables 

    class(MPIController)      ,intent(inout) :: this
    type(VariableContainer)   ,intent(inout) :: varContainer !! Variable container in which to broadcasty
    character(*)              ,intent(in)    :: name !! Name of variable to broadcast
    
    real(rk) ,allocatable ,dimension(:) :: buffer

    integer(ik) :: varIndex ,count

    if (assertions) then 
        call assert(this%isDefined(),"broadcastVarInRow called from undefined mpi controller")
        call assert(this%rowsSetUp,"broadcastVarInRow called from MPI controller before rows set up")
        call assert(varContainer%isDefined(),"Undefined variable container passeed to boradcastVarInRow routine")
    end if

    varIndex = varContainer%getVarIndex(name)
    count = size(varContainer%variables(varIndex)%entry)
    buffer = varContainer%variables(varIndex)%entry
    call MPI_Bcast(buffer,count,MPI_Real8,0,this%rowComm)
    varContainer%variables(varIndex)%entry = buffer

end subroutine broadcastVarInRow
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calculateRowDistData(this,partitionObj,xHaloWidth,numV) 
    !! Initialize distribution DoF data used in exchanging distribution data in a row

    class(MPIController) ,intent(inout)  :: this
    type(Partition)      ,intent(in)     :: partitionObj !! Partition object used to retrieve DoF information
    integer(ik)          ,intent(in)     :: xHaloWidth !! Halo width in x-direction
    integer(ik)          ,intent(in)     :: numV !! Number of cells in v-grid

    integer(ik) ,allocatable ,dimension(:) :: rowNumX ,locNumH

    integer(ik) :: i

    if (assertions) then 
        call assert(this%isDefined(),"calculateRowDistDat called from undefined mpi controller")
        call assert(this%rowsSetUp,"calculateRowDist called from MPI controller before rows set up")
        call assert(partitionObj%isDefined(),"Undefined partition passeed to calculateRowDistDat routine")
    end if

    this%xHaloWidth = xHaloWidth
    rowNumX = partitionObj%getLocNumX()
    locNumH = partitionObj%getLocNumH()

    !Allocate column data for MPI_Gatherv
    allocate(this%colNumX(this%colSize))
    allocate(this%xDispls(this%colSize))
    this%xDispls = 0
    do i = 1, this%colSize
        this%colNumX(i) = rowNumX((i-1)*this%rowSize + 1)
    end do
    
    do i = 2,this%colSize
        this%xDispls(i) = sum(this%colNumX(1:i-1))
    end do

    this%rowNumX = rowNumX(this%worldRank+1)
    this%numV = numV 

    allocate(this%rowNumH(0:this%rowSize-1))
    allocate(this%rowHOffset(0:this%rowSize-1))

    this%rowHOffset = 0
    do i = 1, this%rowSize
        this%rowNumH(i-1) = locNumH(i) 
    end do

    this%numH = sum(this%rowNumH)

    do i = 1, this%rowSize-1
        this%rowHOffset(i) = sum(this%rowNumH(0:i-1))
    end do

    allocate(this%distBuffer(0:this%rowSize-1))

    do i = 1, this%rowSize
        allocate(this%distBuffer(i-1)%entry(0:(this%rowNumX+2)*xHaloWidth*locNumH(i)*numV-1))
        this%distBuffer(i-1)%entry = 0
    end do

end subroutine calculateRowDistData
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine exchangeDistVarInRow(this,varContainer,name) 
    !! Exchanges a distribution variable in a processor row 

    class(MPIController)      ,intent(inout) :: this
    type(VariableContainer)   ,intent(inout) :: varContainer !! Variable container in which to broadcast
    character(*)              ,intent(in)    :: name !! Name of variable to broadcast

    integer(ik) :: varIndex ,bufferFirstIndex,bufferLastIndex,varEntryFirstIndex,varEntryLastIndex
    integer(ik) :: i ,j

    if (assertions) then 
        call assert(this%isDefined(),"exchangeDistVarInRow called from undefined mpi controller")
        call assert(this%rowsSetUp,"exchangeDistVarInRow called from MPI controller before rows set up")
        call assert(varContainer%isDefined(),"Undefined variable container passeed to exchangeDistVarInRow routine")
        call assert(allocated(this%rowNumH),"exchangeDistVarInRow called before calculateRowDistData")
    end if
    varIndex = varContainer%getVarIndex(name)
    !Fill local buffer
    do i = 1,this%rowNumX+2*this%xHaloWidth
        bufferFirstIndex = (i-1)*this%rowNumH(this%rowRank)*this%numV
        bufferLastIndex = i*this%rowNumH(this%rowRank)*this%numV-1
        varEntryFirstIndex = (i-this%xHaloWidth-1)*this%numH*this%numV+this%rowHOffset(this%rowRank)*this%numV + 1
        varEntryLastIndex = (i-this%xHaloWidth-1)*this%numH*this%numV+this%rowHOffset(this%rowRank)*this%numV &
        + this%rowNumH(this%rowRank)*this%numV

        this%distBuffer(this%rowRank)%entry(bufferFirstIndex:bufferLastIndex) &
        = varContainer%variables(varIndex)%entry(varEntryFirstIndex:varEntryLastIndex)

    end do
    do i = 0,this%rowSize - 1
        call MPI_Bcast(this%distBuffer(i)%entry,size(this%distBuffer(i)%entry),MPI_Real8,i,this%rowComm)
    end do

    !Copy data back from buffer
    do j = 0,this%rowSize - 1
        do i = 1,this%rowNumX+2*this%xHaloWidth

            bufferFirstIndex = (i-1)*this%rowNumH(j)*this%numV 
            bufferLastIndex = i*this%rowNumH(j)*this%numV-1
            varEntryFirstIndex = (i-this%xHaloWidth-1)*this%numH*this%numV+this%rowHOffset(j)*this%numV + 1
            varEntryLastIndex = (i-this%xHaloWidth-1)*this%numH*this%numV+this%rowHOffset(j)*this%numV &
            + this%rowNumH(j)*this%numV
            
            varContainer%variables(varIndex)%entry(varEntryFirstIndex:varEntryLastIndex) &
            = this%distBuffer(j)%entry(bufferFirstIndex:bufferLastIndex) 

        end do
    end do

end subroutine exchangeDistVarInRow
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initializeNeighbourPairs(this,periodic) 
    !! Initialize neighbour pairs - if periodic in X adds additional pair to handle this

    class(MPIController) ,intent(inout)  :: this
    logical              ,intent(in)     :: periodic

    integer(ik) :: numPairs

    integer(ik) :: i

    if (assertions) then 
        call assert(this%isDefined(),"initializeNeighbourPairs called from undefined mpi controller")
        call assert(this%rowsSetUp,"initializeNeighbourPairs called from MPI controller before rows set up")
    end if

    numPairs = this%worldSize - this%rowSize
    if (periodic) numPairs = this%worldSize
    allocate(this%neighbourPairs(numPairs))
    do i = 1,this%worldSize - this%rowSize
        allocate(this%neighbourPairs(i)%entry(2))
        this%neighbourPairs(i)%entry(1) = i - 1
        this%neighbourPairs(i)%entry(2) = i - 1 + this%rowSize
    end do
   
    if (periodic) then 

        do i = this%worldSize - this%rowSize + 1,this%worldSize
            allocate(this%neighbourPairs(i)%entry(2))
            this%neighbourPairs(i)%entry(1) = i - 1
            this%neighbourPairs(i)%entry(2) = i - 1 + this%rowSize - this%worldSize
        end do

    end if 
end subroutine initializeNeighbourPairs
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine exchangeVarXHalos(this,varContainer,name,varIsDist) 
    !! Exchanges a halos in x direction for given variable 

    class(MPIController)      ,intent(inout) :: this
    type(VariableContainer)   ,intent(inout) :: varContainer !! Variable container in which to perform halo exchange
    character(*)              ,intent(in)    :: name !! Name of variable to exchange
    logical                   ,intent(in)    :: varIsDist !! Set to true if variable is a distribution

    integer(ik) :: haloDataSize ,varSize ,varIndex ,firstBufferIndex,lastBufferIndex

    integer(ik) :: i

    real(rk) ,allocatable ,dimension(:) :: sendbuf, recvbuf

    if (assertions) then 
        call assert(this%isDefined(),"exchangeVarXHalos called from undefined mpi controller")
        call assert(this%rowsSetUp,"exchangeVarXHalos called from MPI controller before rows set up")
        call assert(varContainer%isDefined(),"Undefined variable container passeed to exchangeVarXHalos routine")
        call assert(allocated(this%neighbourPairs),"exchangeVarXHalos called before initializeNeighbourPairs")
    end if

    if (this%xHaloWidth > 0) then 
        varIndex = varContainer%getVarIndex(name)

        haloDataSize = this%xHaloWidth
        if (varIsDist) haloDataSize = haloDataSize * this%numH * this%numV
        firstBufferIndex = size(varContainer%variables(varIndex)%entry)-3*haloDataSize+1
        lastBufferIndex = size(varContainer%variables(varIndex)%entry) - 2*haloDataSize 
        allocate(sendbuf(0:haloDataSize-1))
        allocate(recvbuf(0:haloDataSize-1))
        do i = 1,size(this%neighbourPairs)
            if (this%worldRank == this%neighbourPairs(i)%entry(1)) &
            sendbuf = varContainer%variables(varIndex)%entry(firstBufferIndex:lastBufferIndex)
            if (this%worldRank == this%neighbourPairs(i)%entry(2)) &
            sendbuf =varContainer%variables(varIndex)%entry(1:haloDataSize) 

            if (this%worldRank == this%neighbourPairs(i)%entry(1)) &
            call MPI_Sendrecv(sendbuf,haloDataSize,MPI_Real8,&
                                this%neighbourPairs(i)%entry(2),i,&
                                recvbuf,haloDataSize,MPI_Real8,&
                                this%neighbourPairs(i)%entry(2),i,this%worldComm,MPI_STATUS_IGNORE) 

            if (this%worldRank == this%neighbourPairs(i)%entry(2)) &
            call MPI_Sendrecv(sendbuf,haloDataSize,MPI_Real8,&
                                this%neighbourPairs(i)%entry(1),i,&
                                recvbuf,haloDataSize,MPI_Real8,&
                                this%neighbourPairs(i)%entry(1),i,this%worldComm,MPI_STATUS_IGNORE) 

            if (this%worldRank == this%neighbourPairs(i)%entry(2)) &
            varContainer%variables(varIndex)%entry(1-haloDataSize:0) = recvbuf
            if (this%worldRank == this%neighbourPairs(i)%entry(1)) &
            varContainer%variables(varIndex)%entry(firstBufferIndex+haloDataSize:lastBufferIndex+haloDataSize) = recvbuf

        end do

    end if

end subroutine exchangeVarXHalos
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setUpRows(this,numProcsX,numProcsH) 
    !! Set up row/column communication in controller

    class(MPIController)      ,intent(inout)  :: this
    integer(ik)               ,intent(in)     :: numProcsX !! Number of processes in the x-direction
    integer(ik)               ,intent(in)     :: numProcsH !! Number of processes in the h-direction

    integer(ik) :: colorRows ,colorColumns
    logical :: mpiInitialized

    colorRows = this%worldRank/numProcsH
    colorColumns = mod(this%worldRank,numProcsH)

    if (assertions) then 
        call MPI_Initialized(mpiInitialized)
        call assert(mpiInitialized .and. (.not. this%rowsSetUp),&
        "Attempted to set up MPI rows when MPI is not initialized or rows already set up ")
        call assert(this%worldSize == numProcsX*numProcsH,"The number of processors in the X and H directions supplied&
        & to the mpi controller constructor must multiply to the total number of processors")
    end if

    call MPI_Comm_split(this%worldComm,colorRows,this%worldRank,this%rowComm)
    call MPI_Comm_size(this%rowComm,this%rowSize)
    call MPI_Comm_rank(this%rowComm,this%rowRank)

    call MPI_Comm_split(this%worldComm,colorColumns,this%worldRank,this%colComm)
    call MPI_Comm_size(this%colComm,this%colSize)
    call MPI_Comm_rank(this%colComm,this%colRank)

    this%rowsSetUp = .true.

    call MPI_Barrier(this%worldComm)

end subroutine setUpRows
!-----------------------------------------------------------------------------------------------------------------------------------
module function isTrueEverywhere (this,input) result(isTrue)
    !! Return true if input is true on every processor

    class(MPIController) ,intent(inout) :: this
    logical              ,intent(inout) :: input
    logical                             :: isTrue

    if (assertions) call assert(this%isDefined(),"isTrueEverywhere called on undefined MPI Controller")

    call MPI_Allreduce(input,isTrue,1,MPI_LOGICAL,MPI_LAND,this%worldComm)

end function isTrueEverywhere
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine gatherVar(this,localVec,globVec,varIsDist) 
    !! Gather variable values into global vector for output on rank 0. If globVec in not allocated it will be allocated to the correct
    !! size, otherwise it is assumed it is of the correct size.

    class(MPIController)                   ,intent(inout) :: this
    real(rk)              ,dimension(:)    ,intent(in)    :: localVec !! Local vector to gather
    real(rk) ,allocatable ,dimension(:)    ,intent(inout) :: globVec !! Global vector on rank 0 to gather into
    logical                                ,intent(in)    :: varIsDist !! Set to true if the gathered variable is a distribution

    integer(ik) ,allocatable ,dimension(:) :: countVec,displs !MPI vars
    real(rk) ,allocatable ,dimension(:) :: sendbuf

    integer(ik) :: sendcount ,haloOffset

    if (assertions) then 
        call assert(this%isDefined(),"gatherVar called from undefined MPI controller")
        call assert(this%rowsSetUp,"gatherVar called before row/column data initialized")
    end if

    !Communicate only in the first column
    if (this%rowRank == 0) then 
        countVec = this%colNumX
        displs = this%xDispls
        if (varIsDist) countVec = countVec * this%numV*this%numH
        if (varIsDist) displs = displs * this%numV*this%numH
        if ((.not. allocated(globVec)) .and. (this%colRank == 0)) allocate(globVec(sum(countVec)))
        haloOffset = this%xHaloWidth 
        if (varIsDist) haloOffset = haloOffset * this%numV*this%numH
        sendbuf = localVec(1+haloOffset:countVec(this%colRank+1)+haloOffset)
        sendcount = countVec(this%colRank+1)
        call MPI_Gatherv(sendbuf,sendcount,MPI_Real8,globVec,countVec,displs,MPI_Real8,0,this%colComm)
    end if

end subroutine gatherVar
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine scatterVar(this,globVec,localVec,varIsDist) 
    !! Scatter global variable values into local vectors from rank 0. 

    class(MPIController)                   ,intent(inout) :: this
    real(rk)              ,dimension(:)    ,intent(inout) :: globVec !! Global vector on rank 0 to scatter from
    real(rk) ,allocatable ,dimension(:)    ,intent(inout) :: localVec !! Local vector to scatter into
    logical                                ,intent(in)    :: varIsDist !! Set to true if the scattered variable is a distribution

    integer(ik) ,allocatable ,dimension(:) :: countVec,displs !MPI vars

    integer(ik) :: rcvcount

    if (assertions) then 
        call assert(this%isDefined(),"scatterVar called from undefined MPI controller")
        call assert(this%rowsSetUp,"scatterVar called before row/column data initialized")
    end if
    !Communicate only in the first column
    if (this%rowRank == 0) then 
        countVec = this%colNumX
        displs = this%xDispls
        if (varIsDist) countVec = countVec * this%numV*this%numH
        if (varIsDist) displs = displs * this%numV*this%numH
        rcvCount = countVec(this%colRank+1)
        if ((.not. allocated(localVec))) allocate(localVec(rcvCount))
        
        call MPI_Scatterv(globVec,countVec,displs,MPI_Real8,localVec,rcvCount,MPI_Real8,0,this%colComm)
    end if

end subroutine scatterVar
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine broadcastReal(this,vec,rootProc) 
    !! Broadcast real vector from rank rootProc

    class(MPIController)      ,intent(inout) :: this
    real(rk) ,dimension(:)    ,intent(inout) :: vec
    integer(ik) ,optional     ,intent(in)    :: rootProc

    integer(ik) :: usedRoot

    usedRoot = 0
    if (present(rootProc)) usedRoot = rootProc

    call MPI_Bcast(vec,size(vec),MPI_Real8,usedRoot,this%worldComm)

end subroutine broadcastReal
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine broadcastInt(this,vec) 
    !! Broadcast integer vector from rank 0

    class(MPIController)      ,intent(inout) :: this
    integer(ik) ,dimension(:) ,intent(inout) :: vec

    call MPI_Bcast(vec,size(vec),MPI_Integer,0,this%worldComm)

end subroutine broadcastInt
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine broadcastLogical(this,vec) 
    !! Broadcast logical vector from rank 0

    class(MPIController)      ,intent(inout) :: this
    logical     ,dimension(:) ,intent(inout) :: vec

    call MPI_Bcast(vec,size(vec),MPI_Logical,0,this%worldComm)

end subroutine broadcastLogical
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine broadcastCharacter(this,vec) 
    !! Broadcast character vector from rank 0

    class(MPIController)      ,intent(inout) :: this
    character(:) ,allocatable ,intent(inout) :: vec
    character(len=512)                       :: buffer

    buffer = vec
    call MPI_Bcast(buffer,len(buffer),MPI_Character,0,this%worldComm)
    vec = trim(buffer)

end subroutine broadcastCharacter
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine barrier(this) 
    !! Calls MPI barrier on commWorld


    class(MPIController)      ,intent(inout) :: this

    call MPI_Barrier(this%worldComm)

end subroutine barrier
!-----------------------------------------------------------------------------------------------------------------------------------
module function allreduceMin (this,input) result(min)
    !! Return min value on all processors

    class(MPIController)   ,intent(inout) :: this
    real(rk)               ,intent(inout) :: input
    real(rk)                              :: min

    call MPI_Allreduce(input,min,1,MPI_Real8,MPI_MIN,this%worldComm)

end function allreduceMin
!-----------------------------------------------------------------------------------------------------------------------------------
module function allreduceMax (this,input) result(max)
    !! Return max value on all processors

    class(MPIController)   ,intent(inout) :: this
    real(rk)               ,intent(inout) :: input
    real(rk)                              :: max

    call MPI_Allreduce(input,max,1,MPI_Real8,MPI_MAX,this%worldComm)

end function allreduceMax
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule mpi_controller_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
