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
submodule (indexing_class) indexing_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the Indexing class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initIndexing(this,part,gridData,varList,xHaloWidth) 
    !! Indexing object initialization routine

    class(Indexing)           ,intent(inout)  :: this
    type(Partition)           ,intent(in)     :: part !! Partition component of initialized Indexing object
    type(Grid)                ,intent(in)     :: gridData !! Grid object used to initialize Indexing
    type(VariableList)        ,intent(in)     :: varList !! Implicit variable list
    integer(ik) ,optional     ,intent(in)     :: xHaloWidth !! Width of halo in x-direction - if passed will assume the x-grid is periodic

    integer(ik) :: numVars, numProcs, numV ,&
                   i ,j ,k 

    integer(ik) ,allocatable ,dimension(:) :: locNumX, locNumH,locMinH

    if (assertions) then 
        call assertPure(part%isDefined(),"Undefined partition passed to indexing object constructor")
        call assertPure(gridData%isDefined(),"Undefined grid passed to indexing object constructor")
        call assertPure(varList%isDefined(),"Undefined variable list passed to indexing object constructor")
    end if

    numVars = varList%getNumVars()
    this%numV = gridData%getNumV()
    this%numX = gridData%getNumX()
    this%numH = gridData%getNumH()
    numProcs = size(part%getMinX())
    this%part = part
    this%varList = varList
    this%xHaloWidth = 0
    this%periodicIndexing = .false. 
    if (present(xHaloWidth)) then 
        this%xHaloWidth = xHaloWidth 
        this%periodicIndexing = .true.
    end if

    allocate(this%procDoFs(numProcs))
    allocate(this%varOffset(numProcs,numVars))
    allocate(this%xDoFs(numProcs))

    this%procDoFs = 0
    this%varOffset = 0 
    this%xDoFs = 0

    locNumX = part%getMaxX() - part%getMinX() + 1
    locNumH = part%getMaxH() - part%getMinH() + 1
    locMinH = part%getMinH()


    do i = 1, numProcs 
        do j = 1, locNumX(i)

            do k = 1, numVars 

                if ((locMinH(i) == 1) .and. (.not. varList%isVarDist(k))) then 

                    this%procDoFs(i) = this%procDoFs(i) + 1

                else if (varList%isVarDist(k)) then 

                    this%procDoFs(i) = this%procDoFs(i) + this%numV * locNumH(i)

                end if

                if ((k < numVars) .and. (j == 1)) this%varOffset(i,k+1) = this%procDoFs(i) 

            end do 

            if (j == 1) this%xDoFs(i) =  this%procDoFs(i)

        end do
    end do

    call this%makeDefined()

end subroutine initIndexing
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getProcDoF (this) result(procDoF)
    !! Getter for procDoFs

    class(Indexing)                       ,intent(in) :: this
    integer(ik) ,allocatable ,dimension(:)            :: procDoF

    if (assertions) call assertPure(this%isDefined(),"getProcDoF called from undefined Indexing object")

    procDoF = this%procDoFs

end function getProcDoF
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumH (this) result(numH)
    !! Getter for numH

    class(Indexing) ,intent(in) :: this
    integer(ik)                 :: numH

    if (assertions) call assertPure(this%isDefined(),"getNumH called for undefined Indexing object")

    numH = this%numH

end function getNumH
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumX (this) result(numX)
    !! Getter for numX

    class(Indexing) ,intent(in) :: this
    integer(ik)             :: numX

    if (assertions) call assertPure(this%isDefined(),"getNumX called for undefined Indexing object")

    numX = this%numX

end function getNumX
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumV (this) result(numV)
    !! Getter for numV

    class(Indexing) ,intent(in) :: this
    integer(ik)             :: numV

    if (assertions) call assertPure(this%isDefined(),"getNumV called for undefined Indexing object")

    numV = this%numV

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

    integer(ik) :: varInd ,procInd ,tempHInd
    logical :: localInd

    localInd = .false. 
    if (present(local)) localInd = local 

    if (assertions) call assertPure(this%isDefined(),"findIndex called from undefined Indexing object")

    varInd = this%varList%getVarIndex(name)

    if (assertions) then 
        
        call assertPure(present(hInd) .eqv. present(vInd),"If called findIndex with supplied hInd/vInd, the other must&
        & be supplied as well")
        call assertPure (this%varList%isVarDist(varInd) .eqv. present(hInd),"If variable passed to findIndex is a&
        & distribution hInd and vInd must be supplied, and if they are supplied the variable must be a distribution")

    end if


    tempHInd = 1 
    if (present(hInd)) tempHInd = hInd

    procInd = this%part%findProc(xInd,tempHInd)

    ind = 0 

    if (.not. localInd) ind = sum(this%procDoFs(1:procInd-1)) 

    ind = ind + (xInd - this%part%getMinXAtInd(procInd))*this%xDoFs(procInd) + this%varOffset(procInd,varInd) + 1

    if (this%varList%isVarDist(varInd)) ind = ind + (hInd - this%part%getMinHAtInd(procInd)) * this%numV + vInd - 1


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

    integer(ik) :: procInd, procIndCheck,usedXInd

    integer(ik) :: checkMinX ,checkMaxX ,thisMinX ,thisMaxX

    if (assertions) call assertPure(this%isDefined(),"findDistIndex called from undefined Indexing object")
    
    if (present(locXInd)) then 
        procInd = this%part%findProc(locXInd,hInd)
    else
        procInd = this%part%findProc(xInd,hInd)
    end if

    usedXInd = xInd - this%part%getMinXAtInd(procInd) + 1
    if (this%periodicIndexing) then 
        procIndCheck = this%part%findProc(xInd,hInd)
        checkMinX = this%part%getMinXAtInd(procIndCheck)
        checkMaxX = this%part%getMaxXAtInd(procIndCheck)
        thisMinX = this%part%getMinXAtInd(procInd)
        thisMaxX = this%part%getMaxXAtInd(procInd)

        if ((checkMinX == 1) .and. (thisMaxX == this%numX)) then 
            if (checkMaxX /= thisMinX - 1) usedXInd = usedXInd + this%numX !puts index into right halo if more than 2
                                                                   !processors in x direction

            !Handle pathological 2 processor case in x direction - might be overkill
            if ((checkMaxX == thisMinX - 1) .and. (usedXInd < 1-this%xHaloWidth)) usedXInd = usedXInd + this%numX
        end if 

        if ((checkMaxX == this%numX) .and. (thisMinX == 1))  then 
            if (checkMinX /= thisMaxX + 1) usedXInd = usedXInd - this%numX !puts index into left halo if more than 2
                                                                   !processors in x direction

            !Handle pathological 2 processor case in x direction - might be overkill
            if ((checkMinX == thisMaxX +1) .and. (usedXInd > thisMaxX+this%xHaloWidth)) usedXInd = usedXInd - this%numX
        end if

    end if
    if (allH) then

        ind = (usedXInd-1)*this%numH*this%numV + (hInd - 1) * this%numV + vInd

    else

        ind = (usedXInd-1)&
            *(this%part%getMaxHAtInd(procInd)-this%part%getMinHAtInd(procInd)+1) * this%numV &
            + (hInd - this%part%getMinHAtInd(procInd)) * this%numV + vInd
    end if


end function findDistIndex
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function mapToGlobalIndices (this,name,inputIndices) result(ind)
    !! Returns list of global indices associated with all possible combinations of inputIndices (IntArray of size 1 or 3) and 
    !! given variable name

    class(Indexing)                       ,intent(in) :: this
    character(*)                          ,intent(in) :: name
    type(intArray) ,dimension(:)          ,intent(in) :: inputIndices
    integer(ik) ,allocatable ,dimension(:)            :: ind

    logical :: compatibleSize 

    integer(ik) :: i ,j ,k ,l

    if (assertions) then
        call assertPure(this%isDefined(),"mapToGlobalIndices called from undefined Indexing object")

        compatibleSize = (this%varList%isVarDist(this%varList%getVarIndex(name)) .and. (size(inputIndices) == 3)) .or. &
                    ((.not. this%varList%isVarDist(this%varList%getVarIndex(name))) .and. (size(inputIndices) == 1))
        call assertPure(compatibleSize,"inputIndices array passed to mapToGlobalIndices must be of size one if the named variable &
        &is not a distribution and 3 otherwise")

        call assertPure(all(inputIndices(1)%entry > 0),"inputIndices(1) - x indices out of bounds (lower)&
        & in mapToGlobalIndices call")
        call assertPure(all(inputIndices(1)%entry <= this%numX),"inputIndices(1) - x indices out of bounds (upper)&
        & in mapToGlobalIndices call")

        if (size(inputIndices) > 1) then 

            call assertPure(all(inputIndices(2)%entry > 0),"inputIndices(2) - h indices out of bounds (lower)&
            & in mapToGlobalIndices call")
            call assertPure(all(inputIndices(2)%entry <= this%numH),"inputIndices(2) - h indices out of bounds (upper)&
            & in mapToGlobalIndices call")

            call assertPure(all(inputIndices(3)%entry > 0),"inputIndices(3) - v indices out of bounds (lower)&
            & in mapToGlobalIndices call")
            call assertPure(all(inputIndices(3)%entry <= this%numV),"inputIndices(3) - v indices out of bounds (upper)&
            & in mapToGlobalIndices call")

        end if
    end if

    if (size(inputIndices) > 1) then 

        allocate(ind(size(inputIndices(1)%entry)*size(inputIndices(2)%entry)*size(inputIndices(3)%entry)))
        l = 1
        do i = 1,size(inputIndices(1)%entry)
            do j = 1,size(inputIndices(2)%entry)
                do k = 1,size(inputIndices(3)%entry)

                    ind(l) = this%findIndex(name,inputIndices(1)%entry(i),inputIndices(2)%entry(j),inputIndices(3)%entry(k))
                    l = l + 1
                end do
            end do
        end do

    else
        
        allocate(ind(size(inputIndices(1)%entry)))
        l = 1
        do i = 1,size(inputIndices(1)%entry)

            ind(l) = this%findIndex(name,inputIndices(1)%entry(i))
            l = l + 1

        end do

    end if

end function mapToGlobalIndices
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function findLocalXIndex (this,xInd,locXInd) result(ind)
    !! Return local non-distribution variable index corresponding to xInd,if present, locXInd is used to identify the processor

    class(Indexing)       ,intent(in) :: this
    integer(ik)           ,intent(in) :: xInd
    integer(ik) ,optional ,intent(in) :: locXInd
    integer(ik)                       :: ind

    integer(ik) :: procInd ,procIndCheck

    integer(ik) :: checkMinX ,checkMaxX ,thisMinX ,thisMaxX

    if (assertions) call assertPure(this%isDefined(),"findLocalXIndex called from undefined Indexing object")

    if (present(locXInd)) then 
        procInd = this%part%findProc(locXInd,1)
    else
        procInd = this%part%findProc(xInd,1)
    end if

    ind = xInd - this%part%getMinXAtInd(procInd) + 1

    !Handle periodic case cleanly 

    if (this%periodicIndexing) then 
        procIndCheck = this%part%findProc(xInd,1)
        checkMinX = this%part%getMinXAtInd(procIndCheck)
        checkMaxX = this%part%getMaxXAtInd(procIndCheck)
        thisMinX = this%part%getMinXAtInd(procInd)
        thisMaxX = this%part%getMaxXAtInd(procInd)

        if ((checkMinX == 1) .and. (thisMaxX == this%numX)) then 
            if (checkMaxX /= thisMinX - 1) ind = ind + this%numX !puts index into right halo if more than 2
                                                                   !processors in x direction

            !Handle pathological 2 processor case in x direction - might be overkill
            if ((checkMaxX == thisMinX - 1) .and. (ind < 1-this%xHaloWidth)) ind = ind + this%numX
        end if 

        if ((checkMaxX == this%numX) .and. (thisMinX == 1))  then 
            if (checkMinX /= thisMaxX + 1) ind = ind - this%numX !puts index into left halo if more than 2
                                                                   !processors in x direction

            !Handle pathological 2 processor case in x direction - might be overkill
            if ((checkMinX == thisMaxX +1) .and. (ind > thisMaxX+this%xHaloWidth)) ind = ind - this%numX
        end if
    end if

end function findLocalXIndex
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getAllIndicesOfVar (this,varInd,procRank) result(ind)
    !! Returns all local indices corresponding to variable with a given index on given processor, accounting for dual grid variables

    class(Indexing)                       ,intent(in) :: this
    integer(ik)                           ,intent(in) :: varInd
    integer(ik)                           ,intent(in) :: procRank 
    integer(ik) ,allocatable ,dimension(:)            :: ind

    logical :: varActive ,isDist
    integer(ik) :: i ,j,varSize ,usedLocNumX
    integer(ik) ,allocatable ,dimension(:) :: locNumX ,locNumH ,locMinH

    if (assertions) call assertPure(this%isDefined(),"getAllIndicesOfVar called from undefined Indexing object")

    locNumX = this%part%getMaxX() - this%part%getMinX() + 1
    locNumH = this%part%getMaxH() - this%part%getMinH() + 1
    locMinH = this%part%getMinH()
    isDist = this%varList%isVarDist(varInd)
    varActive = isDist .or. locMinH(procRank+1) == 1
    if (varActive) then 
        varSize = 1
        if (isDist) varSize =  this%numV * locNumH(procRank+1)
        usedLocNumX = locNumX(procRank+1)
        if (.not. this%periodicIndexing) then 
            if (this%varList%isVarOnDualGrid(varInd) .and. this%part%getMaxXAtInd(procRank+1) == this%numX) &
            usedLocNumX = usedLocNumX - 1
        end if
        allocate(ind(varSize*usedLocNumX))
        do i = 1,usedLocNumX
            ind((i-1)*varSize+1:i*varSize) = [((i-1)*this%xDoFs(procRank+1)+this%varOffset(procRank+1,varInd)+j , j = 1,varSize)]
        end do
    else 
        allocate(ind(0))
    end if
end function getAllIndicesOfVar
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule indexing_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
