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
submodule (d2dv2_derivation_class) d2dv2_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the d2dv2 derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initD2DV2Derivation(this,vSpaceObj,outerV,innerV,vidfdvAtZero,targetH)
    !! Initialize second order velocity derivative derivation

    class(D2DV2Derivation)                  ,intent(inout) :: this
    type(VSpace)                            ,intent(in)    :: vSpaceObj
    type(RealArray) ,optional ,dimension(:) ,intent(in)    :: outerV !! Outer velocity vector corresponding to cell centres for each included harmonic. Defaults to ones.
    type(RealArray) ,optional ,dimension(:) ,intent(in)    :: innerV !! Inner velocity vector corresponding to right cell boundaries for each included harmonic. Defaults to ones.
    type(RealArray) ,optional ,dimension(:) ,intent(in)    :: vidfdvAtZero !! Extrapolation of v_i*df/dv at zero in the form A1*f(v1)+A2*f(v2) where A's correspond to included harmonics (default = 0)
    integer(ik) ,optional                   ,intent(in)    :: targetH !! Harmonic to take derivative of. If not present will return full distribution. 

    integer(ik) :: i ,numV

    real(rk) ,allocatable ,dimension(:) :: vGrid

    if (assertions) call assert (vSpaceObj%isDefined(),"Undefined VSpace object passed to initD2DV2Derivation")

    this%numH = vSpaceObj%getNumH()
    this%dvCopy = vSpaceObj%getVCellWidths()

    numV = size(this%dvCopy)
    allocate(this%dvPlus(numV))
    this%dvPlus = 0
    vGrid = vSpaceObj%getVGrid()
    this%dvPlus(1:numV-1) = vGrid(2:numV) -  vGrid(1:numV-1)

    if (present(targetH)) then 
        this%includedHs = [targetH]
    else
        this%includedHs = [(i,i=1,this%numH)]
    end if

    allocate(this%vidfdvAtZero(size(this%includedHs)))
    allocate(this%outerV(size(this%includedHs)))
    allocate(this%innerV(size(this%includedHs)))

    do i = 1 ,size(this%includedHs)
        this%vidfdvAtZero(i)%entry = real([0,0],kind=rk)
        this%outerV(i)%entry = [(real(1,kind=rk),i=1,numV)]
        this%innerV(i)%entry = [(real(1,kind=rk),i=1,numV)]
    end do

    if (present(outerV)) then 
        if (assertions) then 
            call assert(size(outerV) == size(this%includedHs),&
            "outerV passed to initD2DV2Derivation does not conform with expected number of harmonics")
            do i = 1,size(outerV)
                call assert(size(outerV(i)%entry) == numV,"outerV entry passed to initD2DV2Derivation not of size numV")
            end do
        end if
        this%outerV = outerV
    end if

    if (present(innerV)) then 
        if (assertions) then 
            call assert(size(innerV) == size(this%includedHs),&
            "innerV passed to initD2DV2Derivation does not conform with expected number of harmonics")
            do i = 1,size(innerV)
                call assert(size(innerV(i)%entry) == numV,"innerV entry passed to initD2DV2Derivation not of size numV")
            end do
        end if
        this%innerV = innerV
    end if

    if (present(vidfdvAtZero)) then 
        if (assertions) then 
            call assert(size(vidfdvAtZero) == size(this%includedHs),&
            "vidfdvAtZero passed to initD2DV2Derivation does not conform with expected number of harmonics")
            do i = 1,size(vidfdvAtZero)
                call assert(size(vidfdvAtZero(i)%entry) == 2,"vidfdvAtZero entry passed to initD2DV2Derivation not of size 2")
            end do
        end if
        this%vidfdvAtZero = vidfdvAtZero
    end if
    
    call this%makeDefined()

end subroutine initD2DV2Derivation  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateD2DV2(this,inputArray,indices) result(output)

    class(D2DV2Derivation)             ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,j ,inferredNumX ,numV ,indOffset ,indOffsetOut ,lBoundInput ,hInd

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateD2DV2 called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateD2DV2 must be 1")
        call assertPure(all(indices>0),"indices passed to calculateD2DV2 out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateD2DV2 out of bounds - upper")

    end if

    numV = size(this%dvCopy)
    inferredNumX = size(inputArray(indices(1))%entry)/(this%numH*numV)
    lBoundInput = lbound(inputArray(indices(1))%entry,1)

    allocate(output(inferredNumX*size(this%includedHs)*numV))
    output = 0

    do i = 1,inferredNumX
        do j = 1,size(this%includedHs)
        indOffset = (i-1)*this%numH*numV+(this%includedHs(j)-1)*numV + lBoundInput - 1
        hInd = this%includedHs(j)-minval(this%includedHs) + 1
        indOffsetOut = (i-1)*size(this%includedHs)*numV+(hInd-1)*numV
        
        output(indOffsetOut + 1) = this%innerV(hInd)%entry(1)*(inputArray(indices(1))%entry(indOffset+2)&
                                    - inputArray(indices(1))%entry(indOffset+1))/this%dvPlus(1)&
                                    - inputArray(indices(1))%entry(indOffset+1) * this%vidfdvAtZero(hInd)%entry(1) &
                                    - inputArray(indices(1))%entry(indOffset+2) * this%vidfdvAtZero(hInd)%entry(2)

        output(indOffsetOut + 2 : indOffsetOut + numV - 1) = this%innerV(hInd)%entry(2:numV-1) *&
                                            (inputArray(indices(1))%entry(indOffset+3:indOffset+numV)&
                                            - inputArray(indices(1))%entry(indOffset+2:indOffset+numV-1) )/this%dvPlus(2:numV-1) &
                                - this%innerV(hInd)%entry(1:numV-2) * (inputArray(indices(1))%entry(indOffset+2:indOffset+numV-1) &
                                            - inputArray(indices(1))%entry(indOffset+1:indOffset+numV-2) )/this%dvPlus(1:numV-2)

        output(indOffsetOut + numV) = - this%innerV(hInd)%entry(numV-1) * (inputArray(indices(1))%entry(indOffset+numV) &
                                        - inputArray(indices(1))%entry(indOffset+numV-1) )/this%dvPlus(numV-1)

        output(indOffsetOut + 1:indOffsetOut + numV) = output(indOffsetOut + 1:indOffsetOut+numV) * &
                                                      this%outerV(hInd)%entry / this%dvCopy
        
        end do
    end do

end function calculateD2DV2
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule d2dv2_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
