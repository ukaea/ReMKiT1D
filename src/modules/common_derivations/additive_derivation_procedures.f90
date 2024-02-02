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
submodule (additive_derivation_class) additive_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the additive derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initAdditiveDeriv(this,numDerivs,expectedNumIndices, resultPower, linearCoefficients)
    !! Initialize additive derivation object, allocating the expected number of derivation objects

    class(AdditiveDerivation)          ,intent(inout) :: this
    integer(ik)                        ,intent(in)    :: numDerivs
    integer(ik)                        ,intent(in)    :: expectedNumIndices
    real(rk)   ,optional               ,intent(in)    :: resultPower
    real(rk)   ,optional ,dimension(:) ,intent(in)    :: linearCoefficients

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(numDerivs > 0,"numDerivs passed to initAdditiveDeriv must be greated than 0")
    end if

    allocate(this%derivs(numDerivs))
    allocate(this%derivIndices(numDerivs))
    this%numAddedDerivations = 0
    this%expectedNumIndices = expectedNumIndices

    this%resultPower = real(1,kind=rk)
    if (present(resultPower)) this%resultPower = resultPower

    if (present(linearCoefficients)) then
        if (assertions .or. assertionLvl >= 0) call assert(size(linearCoefficients) == numDerivs,&
        "Passed linearCoefficients in initAdditiveDeriv must conform to number of derivations")

        this%linearCoefficients = linearCoefficients
    else
        this%linearCoefficients = [(real(1,kind=rk),i=1,numDerivs)]
    end if

    call this%makeDefined()

end subroutine initAdditiveDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addDerivation(this,deriv,activeIndices)
    !! Adds derivation to additive derivation object

    class(AdditiveDerivation)   ,intent(inout) :: this
    class(Derivation)           ,intent(in)    :: deriv
    integer(ik) ,dimension(:)   ,intent(in)    :: activeIndices !! Entries of indices object passed to calculateAdditive to be passed to added derivation 

    if (assertions) then 

        call assert(this%isDefined(),"Attempted to add derivation to undefined additive derivation object")
        call assert(deriv%isDefined(),"Attempted to add undefined derivation to additive derivation object")

        call assert(all(activeIndices > 0), &
        "activeIndices passed to addDerivation routine of additiveDerivation class out of bounds - lower")

        call assert(all(activeIndices <= this%expectedNumIndices), &
        "activeIndices passed to addDerivation routine of additiveDerivation class out of bounds - upper &
        &(above number of expected indices)")

        call assert(this%numAddedDerivations < size(this%derivs),&
        "Attempted to add derivation to additive derivation object when no free slot available")

    end if

    this%numAddedDerivations = this%numAddedDerivations + 1 

    allocate(this%derivs(this%numAddedDerivations)%entry,source=deriv)
    this%derivIndices(this%numAddedDerivations)%entry = activeIndices


end subroutine addDerivation  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateAdditive(this,inputArray,indices) result(output)

    class(AdditiveDerivation)          ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i 

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateAdditive called on undefined derivation object")
        call assertPure(size(this%derivs) == this%numAddedDerivations,"calculateAdditive called before all derivations added")

        call assertPure(size(indices) == this%expectedNumIndices,&
        "Number of indices passed to calculateAdditive does not conform with&
        & derivation's expectedNumIndices component")
        
        call assertPure(all(indices>0),"indices passed to calculateAdditive out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateAdditive out of bounds - upper")

    end if

    output = this%derivs(1)%entry%calculate(inputArray,indices(this%derivIndices(1)%entry)) * this%linearCoefficients(1)

    do i = 2,size(this%derivs)
        output = output + this%derivs(i)%entry%calculate(inputArray,indices(this%derivIndices(i)%entry))*this%linearCoefficients(i)
    end do

    output = output**this%resultPower

end function calculateAdditive
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule additive_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
