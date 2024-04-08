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
    case("slopeRatio")
        funcPointer => unarySlopeRatio
    case("superbeeLimiter")
        funcPointer => unarySuperbee
    case("minmodLimiter")
        funcPointer => unaryMinmod
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
    !! The input size is assumed to be evenly divided by size(realParams) and by intParams(1). The input does not have to have size 
    !! equal to size(realParams)*intParams(1). intParams(1) determines the expected size of the output, while intParams(2)
    !! determines which of the strided slices is returnd. See below for example.
    !! 
    !! Contraction is performed on strided slices of the input of length size(realParams). The strides are determined by intParams(1). 
    !! For example, contracting [1,2,3,4,5,6,7,8] with [1,2] yields the intermediate
    !! result [5,11,17,23]. Then if intParams(1) is , we are looking to get a result of length 1. If intParams(2) is then 2, we get [11].
    !! If intParams(1) is 2, we want a result of length 2, sampled evenly from the contracted array (this is useful for ReMKiT1D's
    !! flattened represenation of distributions). If intParams(2) is 1, we get [5,17], and if intParams(2) is 2, we get [11,23]

    real(rk)               ,dimension(:) ,intent(in) :: input 
    real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
    integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
    logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
    real(rk) ,allocatable ,dimension(:)              :: output

    integer(ik) :: i ,offset

    allocate(output(intParams(1)))

    do i = 1,intParams(1)
        offset = (i-1)*size(input)/intParams(1) + (intParams(2)-1)*size(realParams)
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
            output(offset+1:offset+size(realParams)) = input(j)*realParams
        end do
    end do

end function unaryExpand
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function unarySlopeRatio(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper calculating r_i = (u_i - u_{i-n})/(u_{i+n}-u_i) where n is intParams[1]. If the absolute value of the denominator is less than
    !! realParams[1], the result will be realParams[2] (with the appropriate sign) if the numerator is not less than realParams[1],
    !! and 1 otherwise

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

        real(rk) ,allocatable ,dimension(:) :: denominator, numerator 

        numerator = input - unaryShift(input,intParams=intParams)
        denominator = unaryShift(input,intParams=-intParams) - input 

        output = numerator/denominator 

        where (abs(denominator) < realParams(1) .and. abs(numerator - denominator) > realParams(1))
            
            output = realParams(2) * sign(real(1,kind=rk),numerator)*sign(real(1,kind=rk),denominator)
        else where (abs(denominator) < realParams(1) .and. abs(numerator - denominator) < realParams(1))
            output = real(1,kind=rk)

        end where

    end function unarySlopeRatio
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unarySuperbee(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for the superbee limiter. 

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

        output = max(real(0,kind=rk),min(2*input,real(1,kind=rk)),min(input,real(2,kind=rk)))

    end function unarySuperbee
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryMinmod(input,realParams,intParams,logicalParams) result(output)
    !! Unary wrapper for the minmod limiter. 

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

        output = max(real(0,kind=rk),min(input,real(1,kind=rk)))

    end function unaryMinmod
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule unary_transforms_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
