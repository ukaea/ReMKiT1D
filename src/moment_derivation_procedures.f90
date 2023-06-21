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
submodule (moment_derivation_class) moment_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the moment derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMomentDeriv(this,mOrder,h,refVSpace,varPowers,g,multConst)
    !! Initialize moment derivation object

    class(MomentDerivation)             ,intent(inout) :: this
    integer(ik)                         ,intent(in)    :: mOrder !! Moment order 
    integer(ik)                         ,intent(in)    :: h !! Harmonic to take moment of
    type(VSpace)                        ,intent(in)    :: refVSpace !! Velocity space object to make copy of
    real(rk) ,optional ,dimension(:)    ,intent(in)    :: varPowers !! Optional fluid variable powers
    real(rk) ,optional ,dimension(:)    ,intent(in)    :: g !! Optional velocity space vector
    real(rk) ,optional                  ,intent(in)    :: multConst !! Optional multiplicative constant - default 1

    if (assertions) then 

        call assert(refVSpace%isDefined(),"Undefined velocity space object passed to initMomentDeriv")
        call assert(h > 0,"Harmonic index passed to initMomentDeriv out of bounds-lower")
        call assert(h <= refVSpace%getNumH(),"Harmonic index passed to initMomentDeriv out of bounds - upper")

        if (present(g)) call assert(size(g) == refVSpace%getNumV(),"Optional velocity space vector g passed to initMomentDeriv&
        & does not conform to velocity grid")

    end if

    this%mOrder = mOrder 
    this%h = h 
    this%vSpaceCopy = refVSpace

    allocate(this%g(refVSpace%getNumV()))
    this%g = real(1,kind=rk)

    if (present(g)) this%g = g

    if (present(varPowers)) then 
        this%varPowers = varPowers 
    else
        allocate(this%varPowers(0))
    end if

    this%multConst = real(1,kind=rk)
    if (present(multConst)) this%multConst = multConst 

    call this%makeDefined()

end subroutine initMomentDeriv 
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateMomentDeriv(this,inputArray,indices) result(output)

    class(MomentDerivation)            ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,inferredLocNumX

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateMomentDeriv called on undefined derivation object")
        call assertPure(size(indices) == size(this%varPowers)+1,"Number of indices passed to calculateMomentDeriv does not conform&
        & with derivation's varPowers component")
        call assertPure(all(indices>0),"indices passed to calculateMomentDeriv out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateMomentDeriv out of bounds - upper")

        inferredLocNumX = size(inputArray(indices(1))%entry) / (this%vSpaceCopy%getNumH()*this%vSpaceCopy%getNumV())

        do i = 1,size(this%varPowers)
            call assertPure(size(inputArray(indices(i+1))%entry) == inferredLocNumX,"An expected fluid variable passed to&
            & calculateMomentDeriv does not conform to local size inferred from first passed variable &
            &(expected to be a distribution")
        end do

    end if

    output = this%multConst * this%vSpaceCopy%calculateMoment(inputArray(indices(1))%entry,this%h,this%mOrder,this%g)

    do i = 2, size(indices)
        output = output * inputArray(indices(i))%entry ** this%varPowers(i-1)
    end do

end function calculateMomentDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule moment_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
