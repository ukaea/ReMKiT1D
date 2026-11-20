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
submodule (polynomial_fun_derivation_class) polynomial_fun_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the polynomial function derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initPolyFunDeriv(this,polyPowers,polyCoeffs,constCoeff)
    !! Initialize polynomial function derivation object

    class(PolyFunDeriv)       ,intent(inout) :: this
    real(rk) ,dimension(:)    ,intent(in)    :: polyPowers
    real(rk) ,dimension(:)    ,intent(in)    :: polyCoeffs
    real(rk) ,optional        ,intent(in)    :: constCoeff

    if (assertions .or. assertionLvl >= 0) then 
        call assert(size(polyPowers) > 0 ,"polyPowers passed to initPolyFunDeriv must be of size 1 or greater")
        call assert(size(polyPowers) == size(polyCoeffs),"polyPowers and polyCoeffs passed to initPolyDeriv must be of&
                                                                    & same size")
    end if
    this%constCoeff = real(0,kind=rk)
    if (present(constCoeff)) this%constCoeff = constCoeff 

    this%polyPowers = polyPowers

    this%polyCoeffs = polyCoeffs

    call this%makeDefined()

end subroutine initPolyFunDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculatePolyFun(this,inputArray,indices) result(output)

    class(PolyFunDeriv)                ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,srcIdx
    real(rk)    :: p

    if (assertions) then 

        call assertPure(this%isDefined(),"calculatePolyFun called on undefined derivation object")
        call assertPure(size(indices) == 1 .or. size(indices) == size(this%polyPowers),&
        "Number of indices passed to calculatePolyFun must be 1 or equal to the number of non-constant polynomial powers")
        call assertPure(all(indices>0),"indices passed to calculatePolyFun out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculatePolyFun out of bounds - upper")

    end if

    allocate(output(size(inputArray(indices(1))%entry)))
    output(:) = this%constCoeff

    do i = 1, size(this%polyPowers)
        if (size(indices) > 1) then
            srcIdx = indices(i)
        else
            srcIdx = indices(1)
        end if
        p = this%polyPowers(i)
        if (abs(p - 1.0d0) <= unitPowerTolerance) then
            output(:) = output(:) + this%polyCoeffs(i) * inputArray(srcIdx)%entry(:)
        else
            output(:) = output(:) + this%polyCoeffs(i) * (inputArray(srcIdx)%entry(:) ** p)
        end if
    end do

end function calculatePolyFun
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule polynomial_fun_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
