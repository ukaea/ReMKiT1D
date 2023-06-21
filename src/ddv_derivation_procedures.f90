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
submodule (ddv_derivation_class) ddv_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the ddv derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDDVDerivation(this,vSpaceObj,outerV,innerV,vifAtZero,targetH)
    !! Initialize first order velocity derivative derivation

    class(DDVDerivation)                    ,intent(inout) :: this
    type(VSpace)                            ,intent(in)    :: vSpaceObj
    type(RealArray) ,optional ,dimension(:) ,intent(in)    :: outerV !! Outer velocity vector corresponding to cell centres for each included harmonic. Defaults to ones.
    type(RealArray) ,optional ,dimension(:) ,intent(in)    :: innerV !! Inner velocity vector corresponding to right cell boundaries for each included harmonic. Defaults to ones.
    type(RealArray) ,optional ,dimension(:) ,intent(in)    :: vifAtZero !! Extrapolation of v_i*f at zero in the form A1*f(v1)+A2*f(v2) where A's are given for each included harmonic(default = 0)
    integer(ik) ,optional                   ,intent(in)    :: targetH !! Harmonic to take derivative of. If not present will return full distribution/include all harmonics. 

    integer(ik) :: i ,numV

    if (assertions) call assert (vSpaceObj%isDefined(),"Undefined VSpace object passed to initDDVDerivation")

    this%numH = vSpaceObj%getNumH()
    this%dvCopy = vSpaceObj%getVCellWidths()
    this%vInterpCopy = vSpaceObj%getVLinInterp()

    numV = size(this%dvCopy)

    if (present(targetH)) then 
        this%includedHs = [targetH]
    else
        this%includedHs = [(i,i=1,this%numH)]
    end if
    allocate(this%vifAtZero(size(this%includedHs)))
    allocate(this%outerV(size(this%includedHs)))
    allocate(this%innerV(size(this%includedHs)))

    do i = 1 ,size(this%includedHs)
        this%vifAtZero(i)%entry = real([0,0],kind=rk)
        this%outerV(i)%entry = [(real(1,kind=rk),i=1,numV)]
        this%innerV(i)%entry = [(real(1,kind=rk),i=1,numV)]
    end do

    if (present(outerV)) then 
        if (assertions) then 
            call assert(size(outerV) == size(this%includedHs),&
            "outerV passed to initDDVDerivation does not conform with expected number of harmonics")
            do i = 1,size(outerV)
                call assert(size(outerV(i)%entry) == numV,"outerV entry passed to initDDVDerivation not of size numV")
            end do
        end if
        this%outerV = outerV
    end if

    if (present(innerV)) then 
        if (assertions) then 
            call assert(size(innerV) == size(this%includedHs),&
            "innerV passed to initDDVDerivation does not conform with expected number of harmonics")
            do i = 1,size(innerV)
                call assert(size(innerV(i)%entry) == numV,"innerV entry passed to initDDVDerivation not of size numV")
            end do
        end if
        this%innerV = innerV
    end if

    if (present(vifAtZero)) then 
        if (assertions) then 
            call assert(size(vifAtZero) == size(this%includedHs),&
            "vifAtZero passed to initDDVDerivation does not conform with expected number of harmonics")
            do i = 1,size(vifAtZero)
                call assert(size(vifAtZero(i)%entry) == 2,"vifAtZero entry passed to initDDVDerivation not of size 2")
            end do
        end if
        this%vifAtZero = vifAtZero
    end if
    call this%makeDefined()

end subroutine initDDVDerivation  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateDDV(this,inputArray,indices) result(output)

    class(DDVDerivation)               ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,j ,inferredNumX ,numV ,indOffset ,indOffsetOut ,lBoundInput ,hInd

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateDDV called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateDDV must be 1")
        call assertPure(all(indices>0),"indices passed to calculateDDV out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateDDV out of bounds - upper")

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

            output(indOffsetOut + 1) = inputArray(indices(1))%entry(indOffset+1) &
                                        * ((real(1,kind=rk) - this%vInterpCopy(1))*this%innerV(hInd)%entry(1) &
                                        - this%vifAtZero(hInd)%entry(1))+&
                                        inputArray(indices(1))%entry(indOffset+2) &
                                        * (this%vInterpCopy(1)*this%innerV(hInd)%entry(1) - this%vifAtZero(hInd)%entry(2))

            output(indOffsetOut + 2 : indOffsetOut + numV - 1) = - inputArray(indices(1))%entry(indOffset+1:indOffset+numV-2) &
                                                * (real(1,kind=rk) - this%vInterpCopy(1:numV-2))*this%innerV(hInd)%entry(1:numV-2)&
                                                + inputArray(indices(1))%entry(indOffset+2:indOffset+numV-1) * &
                                                ((real(1,kind=rk) - this%vInterpCopy(2:numV-1))*this%innerV(hInd)%entry(2:numV-1) -&
                                                this%vInterpCopy(1:numV-2)*this%innerV(hInd)%entry(1:numV-2)) & 
                                                + inputArray(indices(1))%entry(indOffset+3:indOffset+numV) &
                                                * this%vInterpCopy(2:numV-1)*this%innerV(hInd)%entry(2:numV-1)

            output(indOffsetOut + numV) = - inputArray(indices(1))%entry(indOffset+numV-1) &
                            * (real(1,kind=rk) - this%vInterpCopy(numV-1))*this%innerV(hInd)%entry(numV-1) &
                            - inputArray(indices(1))%entry(indOffset + numV) * &
                            this%vInterpCopy(numV-1)*this%innerV(hInd)%entry(numV-1)

            output(indOffsetOut + 1:indOffsetOut + numV) = output(indOffsetOut + 1:indOffsetOut+numV) * &
                                                            this%outerV(hInd)%entry/this%dvCopy
        
        end do
    end do

end function calculateDDV
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule ddv_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
