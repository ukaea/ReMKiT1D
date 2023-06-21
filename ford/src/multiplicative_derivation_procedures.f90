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
submodule (multiplicative_derivation_class) multiplicative_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the multiplicative derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains 
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMultDeriv(this,innerDeriv,innerIndices,outerDeriv,outerIndices,innerPower,outerPower,innerFuncName)
    !! Initialize multiplicative derivation object

    class(MultiplicativeDerivation)     ,intent(inout) :: this
    class(Derivation)                   ,intent(in)    :: innerDeriv
    integer(ik) ,dimension(:)           ,intent(in)    :: innerIndices
    class(Derivation) ,optional         ,intent(in)    :: outerDeriv 
    integer(ik) ,optional ,dimension(:) ,intent(in)    :: outerIndices
    real(rk) ,optional                  ,intent(in)    :: innerPower
    real(rk) ,optional                  ,intent(in)    :: outerPower
    character(*) ,optional              ,intent(in)    :: innerFuncName 

    if (assertions) then 
        call assert(innerDeriv%isDefined(),"innerDeriv passed to initMultDeriv not defined")
        if (present(outerDeriv)) then 
            call assert(outerDeriv%isDefined(),"outerDeriv passed to initMultDeriv not defined")
            call assert(present(outerIndices),"If outerDeriv passed to initMultDeriv outerIndices must also be passed")
        end if

    end if

    allocate(this%innerDeriv,source=innerDeriv)
    allocate(this%innerIndices,source=innerIndices)

    if (present(outerDeriv)) then 
        allocate(this%outerDeriv,source=outerDeriv)
        allocate(this%outerIndices,source=outerIndices)
    end if

    this%innerPower = real(1,kind=rk)
    this%outerPower = real(1,kind=rk)

    if (present(innerPower)) this%innerPower = innerPower
    if (present(outerPower)) this%outerPower = outerPower

    if (present(innerFuncName)) allocate(this%innerFuncName,source=innerFuncName)

    call this%makeDefined()

end subroutine initMultDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateMultiplicative(this,inputArray,indices) result(output)

    class(MultiplicativeDerivation)    ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateMultiplicative called on undefined derivation object")
        
        call assertPure(all(indices>0),"indices passed to calculateMultiplicative out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateMultiplicative out of bounds - upper")

    end if

    output = this%innerDeriv%calculate(inputArray,indices(this%innerIndices))

    if (allocated(this%innerFuncName)) then 

        select case(this%innerFuncName)
        case("exp")
            output = exp(output)**this%innerPower
        case("log")
            output = log(output)**this%innerPower
        case("sin")
            output = sin(output)**this%innerPower
        case("cos")
            output = cos(output)**this%innerPower
        case("abs")
            output = abs(output)**this%innerPower
        case("tan")
            output = tan(output)**this%innerPower
        case("atan")
            output = atan(output)**this%innerPower
        case("asin")
            output = asin(output)**this%innerPower
        case("acos")
            output = acos(output)**this%innerPower
        case("sign")
            output = sign(real(1,kind=rk),output)
        case("erf")
            output = erf(output)**this%innerPower
        case("erfc")
            output = erfc(output)**this%innerPower
        case default 
            error stop "unsupported function name passed to multiplicative derivation"
        end select 

    else
        output = output**this%innerPower
    end if

    if (allocated(this%outerDeriv)) then 
        output = output*this%outerDeriv%calculate(inputArray,indices(this%outerIndices))**this%outerPower
    end if
    
end function calculateMultiplicative
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule multiplicative_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
