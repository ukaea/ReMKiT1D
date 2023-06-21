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
submodule (vel_tensor_prod_derivation_class) vel_tensor_prod_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the velocity tensor product derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVTProdDeriv(this,refVSpace,velVec,power)
    !! Initialize velocity tensor product derivation object

    class(VelTProdDerivation)           ,intent(inout) :: this
    type(VSpace)                        ,intent(in)    :: refVSpace !! Reference velocity space
    real(rk) ,optional ,dimension(:)    ,intent(in)    :: velVec !! Optional velocity space vector. Defaults to velocity grid values
    real(rk) ,optional                  ,intent(in)    :: power !! Optional power to raise the shifted velocity vector to. Defaults to 1

    if (assertions) then
        call assert(refVSpace%isDefined(),"Undefined refVSpace passed to initVTProdDeriv")
        if (present(velVec)) call assert(size(velVec)==refVSpace%getNumV(),&
        "velVec passed to initVTProdDeriv does not conform to number of velocity grid points")
    end if

    this%v = refVSpace%getVGrid()
    if (present(velVec)) this%v = velVec

    this%power = real(1,kind=rk)
    if (present(power)) this%power = power

    call this%makeDefined()

end subroutine initVTProdDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateVTProdDeriv(this,inputArray,indices) result(output)

    class(VelTProdDerivation)          ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output
    
    if (assertions) then 

        call assertPure(this%isDefined(),"calculateVTProdDeriv called on undefined derivation object")
        call assertPure(any(size(indices) == [1,2]),"Number of indices passed to calculateVTProdDeriv must be 1 or 2")
        call assertPure(all(indices>0),"indices passed to calculateVTProdDeriv out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateVTProdDeriv out of bounds - upper")

    end if

    if (size(indices) == 2) then
        output = shiftedFlatTensorProduct(this%v,inputArray(indices(1))%entry,inputArray(indices(2))%entry,this%power)
    else
        output = shiftedFlatTensorProduct(this%v,inputArray(indices(1))%entry,power=this%power)
    end if


end function calculateVTProdDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule vel_tensor_prod_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
