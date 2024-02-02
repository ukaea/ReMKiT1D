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
submodule (gen_int_polynomial_fun_derivation_class) gen_int_polynomial_fun_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the generalized integer powered polynomial function derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initGenIntPolyFunDeriv(this,polyPowers,polyCoeffs,maxPowers,funcName,multConst)
    !! Initialize generalized integer powered polynomial function derivation object

    class(GenIntPolyFunDeriv)       ,intent(inout) :: this
    type(IntArray) ,dimension(:)    ,intent(in)    :: polyPowers
    real(rk)       ,dimension(:)    ,intent(in)    :: polyCoeffs
    integer(ik)    ,dimension(:)    ,intent(in)    :: maxPowers
    character(*) ,optional          ,intent(in)    :: funcName 
    real(rk) ,optional              ,intent(in)    :: multConst

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then

        call assert(size(polyPowers) == size(polyCoeffs),&
        "polyPowers and polyCoeffs passed to initGenIntPolyFunDeriv must be of same size")

        do i = 1,size(polyPowers)
            call assert(size(polyPowers(i)%entry)==size(maxPowers),&
            "All polyPowers passed to initGenIntPolyFunDeriv must be of same size as maxPowers")
            call assert(all(polyPowers(i)%entry >= 0),"Negative polyPowers are not supported by GenIntPolyFunDeriv")
            call assert(all(polyPowers(i)%entry <= maxPowers),&
            "All polyPowers must be less than or equal to corresponding maxPowers")
        end do

    end if

    this%polyPowers = polyPowers
    this%polyCoeffs = polyCoeffs
    this%maxPowers = maxPowers

    this%multConst = real(1,kind=rk)
    if (present(multConst)) this%multConst = multConst

    if (present(funcName)) this%funcName = funcName

    call this%makeDefined()

end subroutine initGenIntPolyFunDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateGenIntPolyFun(this,inputArray,indices) result(output)

    class(GenIntPolyFunDeriv)          ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i ,j

    type(RealArrayD2) ,allocatable ,dimension(:) :: varBuffer

    real(rk) ,allocatable ,dimension(:)          :: monomialBuffer


    if (assertions) then 

        call assertPure(this%isDefined(),"calculateGenIntPolyFun called on undefined derivation object")
        call assertPure(size(indices) == size(this%maxPowers),&
        "calculateGenIntPolyFun called with indices array not conforming to maxPowers")
        call assertPure(all(indices>0),"indices passed to calculateGenIntPolyFun out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateGenIntPolyFun out of bounds - upper")

    end if

    allocate(varBuffer(size(indices)))

    do i = 1,size(varBuffer)
        allocate(varBuffer(i)%entry(this%maxPowers(i)))

        do j = 1,this%maxPowers(i)
            varBuffer(i)%entry(j)%entry = inputArray(indices(i))%entry**j
        end do
    end do

    allocate(output,mold=inputArray(indices(1))%entry)

    output = 0

    monomialBuffer = output

    do i = 1,size(this%polyCoeffs)
        monomialBuffer = this%polyCoeffs(i)
        do j = 1,size(varBuffer)
            if (this%polyPowers(i)%entry(j) > 0) monomialBuffer = monomialBuffer &
                                                                 * varBuffer(j)%entry(this%polyPowers(i)%entry(j))%entry
        end do
        output = output + monomialBuffer
    end do

    if (allocated(this%funcName)) then 

        select case(this%funcName)
        case("exp")
            output = exp(output)
        case("log")
            where (abs(output) < epsilon(output))
                output = real(1,kind=rk)
            end where
            output = log(output)
        case("sin")
            output = sin(output)
        case("cos")
            output = cos(output)
        case("abs")
            output = abs(output)
        case("tan")
            output = tan(output)
        case("atan")
            output = atan(output)
        case("asin")
            output = asin(output)
        case("acos")
            output = acos(output)
        case("sign")
            output = sign(real(1,kind=rk),output)
        case("erf")
            output = erf(output)
        case("erfc")
            output = erfc(output)
        case default 
            error stop "unsupported function name passed to generalized integer polynomial derivation"
        end select 

    end if

    output = output * this%multConst

end function calculateGenIntPolyFun
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule gen_int_polynomial_fun_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
