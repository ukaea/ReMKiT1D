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
submodule (partition_class) partition_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the partition class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initSimplePartition(this,numProcsX,numProcsH,numX,numH) 
    !! Partition initialization routine - assuming even distributions in x and h directions

    class(Partition)          ,intent(inout)  :: this
    integer(ik)               ,intent(in) :: numProcsX !! Number of processes in x direction
    integer(ik)               ,intent(in) :: numProcsH !! Number of processes in h direction
    integer(ik)               ,intent(in) :: numX !! Total number of x grid points
    integer(ik)               ,intent(in) :: numH !! Total number of h grid points 

    integer(ik)                               :: i, j ,k

    if (assertions .or. assertionLvl >= 0) then 
        call assertPure(mod(numX,numProcsX) == 0,"Number of x grid points given to partition constructor is not divisible by &
        &number of processes in x direction")
        call assertPure(mod(numH,numProcsH) == 0,"Number of harmonics given to partition constructor is not divisible by &
        &number of processes in h direction")
        
    end if

    if (numH > 1 .and. numH/numProcsH == 1) print*,"WARNING: 1 harmonic per processor - crash likely in many kinetic runs"

    allocate(this%minX(numProcsX*numProcsH))
    allocate(this%maxX(numProcsX*numProcsH))
    allocate(this%minH(numProcsX*numProcsH))
    allocate(this%maxH(numProcsX*numProcsH))
    allocate(this%locX(numProcsX*numProcsH))
    allocate(this%locH(numProcsX*numProcsH))

    this%locX = numX/numProcsX
    this%locH = numH/numProcsH

    k = 1
    do i = 1,numProcsX 
        do j = 1,numProcsH 
            this%minX(k) = (i-1) * this%locX(i) + 1
            this%maxX(k) = i * this%locX(i)
            this%minH(k) = (j-1) * this%locH(i) + 1
            this%maxH(k) = j * this%locH(i)
            k = k + 1
        end do 
    end do

    call this%makeDefined()

end subroutine initSimplePartition 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initPartition(this,minX, maxX, minH, maxH) 
    !! Partition initialization routine - general

    class(Partition)          ,intent(inout)  :: this
    integer(ik) ,dimension(:) ,intent(in)     :: minX !! First x-index in each partition
    integer(ik) ,dimension(:) ,intent(in)     :: maxX !! Last x-index in each partition
    integer(ik) ,dimension(:) ,intent(in)     :: minH !! First h-index in each partition
    integer(ik) ,dimension(:) ,intent(in)     :: maxH !! Last x-index in each partition

    if (assertions .or. assertionLvl >= 0) call assertPure(all([size(maxX),size(minH),size(maxH)] == size(minX)),&
    "All arrays passed to generic &
    &partition constructor must be same length")

    call this%makeDefined()

    this%minX = minX 
    this%maxX = maxX 
    this%minH = minH
    this%maxH = maxH

    this%locX = maxX - minX + 1
    this%locH = maxH - minH + 1

end subroutine initPartition
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMinX (this) result(minX)
    !! Getter for minX

    class(Partition)                       ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: minX

    if (assertions) call assertPure(this%isDefined(),"Attempted to get minX vector from undefined partition")

    minX = this%minX

end function getMinX
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMinXAtInd (this,ind) result(minX)
    !! Return minX at index ind

    class(Partition)   ,intent(in) :: this
    integer(ik)        ,intent(in) :: ind
    integer(ik)                    :: minX

    if (assertions) call assertPure(this%isDefined(),"Attempted to get minX value from undefined partition")

    minX = this%minX(ind)

end function getMinXAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxX (this) result(maxX)
    !! Getter for maxX

    class(Partition)                      ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)            :: maxX

    if (assertions) call assertPure(this%isDefined(),"Attempted to get maxX vector from undefined partition")

    maxX = this%maxX

end function getMaxX
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxXAtInd (this,ind) result(maxX)
    !! Return maxX at index ind

    class(Partition)   ,intent(in) :: this
    integer(ik)        ,intent(in) :: ind
    integer(ik)                    :: maxX

    if (assertions) call assertPure(this%isDefined(),"Attempted to get maxX value from undefined partition")

    maxX = this%maxX(ind)
    
end function getMaxXAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMinH (this) result(minH)
    !! Getter for minH

    class(Partition)                      ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)            :: minH

    if (assertions) call assertPure(this%isDefined(),"Attempted to get minH vector from undefined partition")

    minH = this%minH

end function getMinH
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMinHAtInd (this,ind) result(minH)
    !! Return minH at index ind

    class(Partition)   ,intent(in) :: this
    integer(ik)        ,intent(in) :: ind
    integer(ik)                    :: minH

    if (assertions) call assertPure(this%isDefined(),"Attempted to get minH value from undefined partition")

    minH = this%minH(ind)
    
end function getMinHAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxH (this) result(maxH)
    !! Getter for maxH

    class(Partition)                       ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: maxH

    if (assertions) call assertPure(this%isDefined(),"Attempted to get maxX vector from undefined partition")

    maxH = this%maxH

end function getMaxH
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxHAtInd (this,ind) result(maxH)
    !! Return maxH at index ind

    class(Partition)   ,intent(in) :: this
    integer(ik)        ,intent(in) :: ind
    integer(ik)                    :: maxH

    if (assertions) call assertPure(this%isDefined(),"Attempted to get maxH value from undefined partition")

    maxH = this%maxH(ind)
    
end function getMaxHAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLocNumX (this) result(locX)
    !! Getter for locX

    class(Partition)                       ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: locX

    if (assertions) call assertPure(this%isDefined(),"Attempted to get locX from undefined partition")

    locX = this%locX

end function getLocNumX
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLocNumH (this) result(locH)
    !! Getter for locH

    class(Partition)                       ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)             :: locH

    if (assertions) call assertPure(this%isDefined(),"Attempted to get locH from undefined partition")

    locH = this%locH

end function getLocNumH
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function findProc (this,xInd,hInd) result(procInd)
    !! Return partition index which contains xInd and hInd

    class(Partition) ,intent(in) :: this
    integer(ik)      ,intent(in) :: xInd
    integer(ik)      ,intent(in) :: hInd
    integer(ik)                  :: procInd

    integer(ik) ,allocatable ,dimension(:) :: indVec
    logical     ,allocatable ,dimension(:) :: mask

    if (assertions) call assertPure(this%isDefined(),"Called findProc for undefined partition")

    mask = (this%minX <= xInd) .and. (this%maxX >= xInd) .and. (this%minH <= hInd) .and. (this%maxH >= hInd)
    indVec = findIndices(mask)

    if (assertions) call assertPure(size(indVec) == 1, "findProc failed to find processor index for passed x and h indices")

    procInd = indVec(1)

end function findProc
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function filterCoords (this,ind,coords,normalize) result(res)
    !! Filters coordinate list of shape 1,: or 3,: based on partition data, leaving only those present at partition index ind. 
    !! If normalize is true will shift all values in each dimension to so that the smallest value is 1

    class(Partition)             ,intent(in) :: this
    integer(ik)                  ,intent(in) :: ind
    integer(ik) ,dimension(:,:)  ,intent(in) :: coords
    logical  ,optional           ,intent(in) :: normalize
    integer(ik) ,dimension(:,:) ,allocatable :: res

    logical :: norm 

    integer(ik) ,allocatable ,dimension(:) :: presentXCoords ,presentHCoords ,presentVCoords

    integer(ik) :: minX ,maxX ,minH ,maxH ,offsetX ,offsetH ,offsetV
    if (assertions) then 
        call assertPure(this%isDefined(),"filterCoords called for undefined partition")
        call assertPure(size(coords,1)==1 .or. size(coords,1)==3,"coords passed to filter coords must have shape (1,:) or (3,:)")
    end if

    norm = .false. 
    if (present(normalize)) norm = normalize

    minX = this%minX(ind)
    maxX = this%maxX(ind)
    minH = this%minH(ind)
    maxH = this%maxH(ind)

    offsetX = 1 
    offsetH = 1 
    offsetV = 1 

    presentXCoords = removeDupeInts(coords(1,:))
    presentXCoords = pack(presentXCoords,presentXCoords >= minX .and. presentXCoords <= maxX)
    if (norm) offsetX = minval(presentXCoords)
    presentXCoords = presentXCoords - offsetX + 1
    if (size(coords,1) > 1) then 
        presentHCoords = removeDupeInts(coords(2,:))
        presentHCoords = pack(presentHCoords,presentHCoords >= minH .and. presentHCoords <= maxH)
        if (norm) offsetH = minval(presentHCoords)
        presentHCoords = presentHCoords - offsetH + 1

        presentVCoords = removeDupeInts(coords(3,:))
        if (norm) offsetV = minval(presentVCoords)
        presentVCoords = presentVCoords - offsetV + 1
        res = allCombinations([IntArray(presentXCoords),IntArray(presentHCoords),IntArray(presentVCoords)])
    else
        res = allCombinations([IntArray(presentXCoords)])
    end if

end function filterCoords
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule partition_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
