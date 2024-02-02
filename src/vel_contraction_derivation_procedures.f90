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
submodule (vel_contraction_derivation_class) vel_contraction_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the velocity contraction derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initContracDeriv(this,h,g,refVSpace,expH)
    !! Initialize velocity contraction derivation object

    class(VelContracDerivation)         ,intent(inout) :: this
    integer(ik)                         ,intent(in)    :: h !! Harmonic to contract with
    real(rk)           ,dimension(:)    ,intent(in)    :: g !! Velocity space vector 
    type(VSpace)                        ,intent(in)    :: refVSpace !! Reference velocity space
    integer(ik)   ,optional             ,intent(in)    :: expH !! Expected number of harmonics

    if (assertions .or. assertionLvl >= 0) then 

        call assert(refVSpace%isDefined(),"Undefined velocity space object passed to initContracDeriv")
        call assert(h > 0,"Harmonic index passed to initContracDeriv out of bounds-lower")
        call assert(h <= refVSpace%getNumH(),"Harmonic index passed to initContracDeriv out of bounds - upper")

        call assert(size(g) == refVSpace%getNumV(),"Velocity space vector g passed to initContracDeriv&
        & does not conform to velocity grid")

    end if

    this%h = h 
    this%g = g
    this%numH = refVSpace%getNumH()
    if (present(expH)) this%numH = expH

    call this%makeDefined()

end subroutine initContracDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateContracDeriv(this,inputArray,indices) result(output)

    class(VelContracDerivation)        ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,inferredLocNumX,lboundInput,offset

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateContracDeriv called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateContracDeriv must be 1")
        call assertPure(all(indices>0),"indices passed to calculateContracDeriv out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateContracDeriv out of bounds - upper")

    end if

    inferredLocNumX = size(inputArray(indices(1))%entry) / (this%numH*size(this%g))

    allocate(output(inferredLocNumX))
    lboundInput = lbound(inputArray(indices(1))%entry,1)
    do i = 1,inferredLocNumX
        offset = lboundInput + (i-1)*this%numH*size(this%g) + (this%h-1)*size(this%g) - 1
        output(i) = dot_product(inputArray(indices(1))%entry(offset+1:offset+size(this%g)),this%g)
    end do
end function calculateContracDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule vel_contraction_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
