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
submodule (unary_transforms) unary_transforms_procedures
!! author: Stefan Mijin  
!!
!! Contains the implementations of various unary transform routines

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine associateFunctionPointer(name,funcPointer) 
    !! Returns a function pointer to a unary transformation based on a string name

    character(*)                                  ,intent(in)    :: name
    procedure(realArrayFunctionGenParam) ,pointer ,intent(inout) :: funcPointer

    select case(name)
    case("exp")
        funcPointer => unaryExp
    case("log")
        funcPointer => unaryLog
    case("sin")
        funcPointer => unarySin
    case("cos")
        funcPointer => unaryCos
    case("abs")
        funcPointer => unaryAbs
    case("tan")
        funcPointer => unaryTan
    case("atan")
        funcPointer => unaryAtan
    case("asin")
        funcPointer => unaryAsin
    case("acos")
        funcPointer => unaryAcos
    case("sign")
        funcPointer => unarySign
    case("erf")
        funcPointer => unaryErf
    case("erfc")
        funcPointer => unaryErfc
    case("ipow")
        funcPointer => unaryPowInt
    case("rpow")
        funcPointer => unaryPowReal
    case("shift")
        funcPointer => unaryShift
    case("cont")
        funcPointer => unaryContract
    case("expand")
        funcPointer => unaryExpand
    case("none")
        funcPointer => null()
    case default 
        error stop "unknown unary transform name passed to associateFunctionPointer"
    end select

end subroutine associateFunctionPointer
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryLog(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for natural log. String name "log".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = log(input)

end function unaryLog
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryExp(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for exp. String name "exp".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = exp(input)

end function unaryExp
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unarySin(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for sin. String name "sin".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = sin(input)

end function unarySin
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryCos(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for cos. String name "cos".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = cos(input)

end function unaryCos
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryAbs(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for absolute value. String name "abs".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = abs(input)

end function unaryAbs
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unarySign(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for sign function. String name "sign".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = sign(real(1,kind=rk),input)

end function unarySign
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryTan(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for tan. String name "tan".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output
    
    output = tan(input)

end function unaryTan
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryAsin(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for arcsin. String name "asin".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = asin(input)

end function unaryAsin
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryAcos(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for arccos. String name "acos".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = acos(input)

end function unaryAcos
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryAtan(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for arctan. String name "atan".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = atan(input)

end function unaryAtan
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryErf(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for the error function. String name "erf".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = erf(input)

end function unaryErf
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryErfc(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for the complementary error function. String name "erfc".

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = erfc(input)
    
end function unaryErfc
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryPowReal(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for the raising to a real power. String name "rpow". Uses realParams(1) as the power.

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = input ** realParams(1)

end function unaryPowReal
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryPowInt(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for the raising to an integer power. String name "ipow". Uses intParams(1) as the power.

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    output = input ** intParams(1)

end function unaryPowInt
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryShift(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper cyclically permuting input. String name "shift". Uses intParams(1) as the shift amount.

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    allocate(output(size(input)))

    if (intParams(1)<0) then
        output(:size(output)+intParams(1)) = input(1-intParams(1):)
        output(size(output)+intParams(1)+1:) = input(:-intParams(1))
    else 
        output(:intParams(1)) = input(size(input)-intParams(1)+1:)
        output(intParams(1)+1:) = input(:size(input)-intParams(1))
    end if

end function unaryShift
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryContract(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper contracting an array with a smaller array. String name "cont". Uses realParams as the contracting array.
    !! After contraction it uses intParams(1) to determine the expected output size. The input array should be evenly divided by 
    !! size(realParams) and intParams(1). intParams(2) is then used to determine which sub-array of size intParams(1) is returned

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    integer(ik) :: i ,offset
    real(rk) ,allocatable ,dimension(:)              :: contractionResult

    allocate(output(intParams(1)))

    do i = 1,intParams(1)
        offset = (i-1)*size(realParams)*size(input)/intParams(1) + (intParams(2)-1)*size(realParams)
        output(i) = dot_product(input(offset+1:offset+size(realParams)),realParams)
    end do
end function unaryContract
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unaryExpand(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper expanding an array with another array. String name "expand". Uses realParams as the expanding array.
    !! After expansion it uses intParams(1) to determine how many times the result should be copied. The result is
    !! output((i-1)*size(input)*size(realParams)+(j-1)*size(realParams)+k) = input(j)*realParams(k) where j=1,size(input) and i=1,intParams(1)

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    integer(ik) i,j,offset
    allocate(output(intParams(1)*size(input)*size(realParams)))

    do i = 1,intParams(1)
        do j = 1,size(input)
            offset = (i-1)*size(realParams)*size(input) + (j-1)*size(realParams)
            output(offset+1:offset+size(realParams)) = input(i)*realParams
        end do
    end do

end function unaryExpand
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule unary_transforms_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
