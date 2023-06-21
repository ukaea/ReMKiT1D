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
submodule (simple_derivation_class) simple_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the simple derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initSimpleDeriv(this,varPowers,multConst)
    !! Initialize simple derivation object

    class(SimpleDerivation)   ,intent(inout) :: this
    real(rk) ,dimension(:)    ,intent(in)    :: varPowers
    real(rk) ,optional        ,intent(in)    :: multConst

    this%multConst = real(1,kind=rk)
    if (present(multConst)) this%multConst = multConst 

    this%varPowers = varPowers

    call this%makeDefined()

end subroutine initSimpleDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateSimple(this,inputArray,indices) result(output)

    class(SimpleDerivation)            ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i 

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateSimple called on undefined derivation object")
        call assertPure(size(indices) == size(this%varPowers),"Number of indices passed to calculateSimple does not conform with&
        & derivation's varPowers component")
        call assertPure(all(indices>0),"indices passed to calculateSimple out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateSimple out of bounds - upper")

    end if

    output = this%multConst * inputArray(indices(1))%entry ** this%varPowers(1)
    
    do i = 2, size(indices)
        output = output * inputArray(indices(i))%entry ** this%varPowers(i)
    end do

end function calculateSimple
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule simple_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
