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
submodule (central_diff_grad_derivation_class) central_diff_grad_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the central difference derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCentralDiffDeriv(this,refGeometry,refPartition,procRank,varPowers,multConst)
    !! Initialize central difference derivation object

    class(CentralDiffDerivation)        ,intent(inout) :: this
    type(Geometry)                      ,intent(in)    :: refGeometry !! Geometry object used to calculate central difference
    type(Partition)                     ,intent(in)    :: refPartition !! Partition object used to calculate central difference
    integer(ik)                         ,intent(in)    :: procRank !! Current processor rank
    real(rk) ,optional ,dimension(:)    ,intent(in)    :: varPowers !! Optional fluid variable powers
    real(rk) ,optional                  ,intent(in)    :: multConst !! Optional multiplicative constant - default 1

    logical :: pGrid ,containsLeftBoundary ,containsRightBoundary

    integer(ik) :: minX, maxX ,i 

    real(rk) ,allocatable ,dimension(:) :: dx ,fullLinInterp ,linInterp 
    real(rk) :: locDx

    if (assertions) then 

        call assert(refGeometry%isDefined(),"Undefined geometry object passed to initCentralDiffDeriv")
        call assert(refPartition%isDefined(),"Undefined partition object passed to initCentralDiffDeriv")

    end if

    dx = refGeometry%getCellWidths()
    minX = refPartition%getMinXAtInd(procRank+1)
    maxX = refPartition%getMaxXAtInd(procRank+1)
    containsLeftBoundary = minX == 1
    containsRightBoundary = maxX == size(dx)
    pGrid = refGeometry%isPeriodic() 

    allocate(fullLinInterp,source = refGeometry%getLinInterp())
    allocate(linInterp(0:maxX-minX+1))
    linInterp = fullLinInterp(minX - 1:maxX)

    allocate(this%leftMult(0:maxX-minX+2))
    this%leftMult = 0

    allocate(this%rightMult(0:maxX-minX+2))
    this%rightMult = 0

    allocate(this%centralMult(0:maxX-minX+2))
    this%centralMult = 0

    do i = 1, maxX - minX + 1
        locDx = dx(minX-1+i)
        this%centralMult(i) = (- linInterp(i-1) + (real(1,kind=rk)-linInterp(i)))/locDx
        this%leftMult(i-1) = - (real(1,kind=rk)-linInterp(i-1))/locDx
        this%rightMult(i+1) = linInterp(i)/locDx
        if (.not. pGrid) then 
            if (containsLeftBoundary .and. i == 1) then 
                this%leftMult(i-1) = 0
                this%centralMult(i) = -2/(dx(minX-1+i)+dx(minX+i))
                this%rightMult(i+1) = 2/(dx(minX-1+i)+dx(minX+i))
            end if

            if (containsRightBoundary .and. i == maxX - minX+1) then 
                this%leftMult(i-1) = -2/(dx(minX-1+i)+dx(minX+i-2))
                this%centralMult(i) = 2/(dx(minX-1+i)+dx(minX+i-2))
                this%rightMult(i+1) = 0
            end if
        end if

    end do
    

    if (present(varPowers)) then 
        this%varPowers = varPowers 
    else
        allocate(this%varPowers(0))
    end if

    this%multConst = real(1,kind=rk)
    if (present(multConst)) this%multConst = multConst 

    call this%makeDefined()

end subroutine initCentralDiffDeriv 
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateCentralDiffDeriv(this,inputArray,indices) result(output)

    class(CentralDiffDerivation)       ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,upperBound

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateCentralDiffDeriv called on undefined derivation object")
        call assertPure(size(indices) == size(this%varPowers)+1,"Number of indices passed to calculateCentralDiffDeriv does &
        &not conform with derivation's varPowers component")
        call assertPure(all(indices>0),"indices passed to calculateCentralDiffDeriv out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateCentralDiffDeriv out of bounds - upper")

    end if
    allocate(output,mold=inputArray(indices(1))%entry)
    output = 0
    upperBound = size(this%centralMult)-2
    output(1:upperBound) = this%multConst * (inputArray(indices(1))%entry(1:upperBound) * this%centralMult(1:upperBound) &
                                             +inputArray(indices(1))%entry(0:upperBound-1) * this%leftMult(0:upperBound-1) &
                                             + inputArray(indices(1))%entry(2:upperBound+1) * this%rightMult(2:upperBound+1))

    do i = 2, size(indices)
        output = output * inputArray(indices(i))%entry ** this%varPowers(i-1)
    end do

end function calculateCentralDiffDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule central_diff_grad_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
