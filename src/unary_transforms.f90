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
module unary_transforms
    !! Contains various unary transforms for the calculation tree

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging ,assertions
    use basic_interfaces            ,only: realArrayFunctionGenParam
    use assertion_utility       
    use support_types               

    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine associateFunctionPointer(name,funcPointer) 
        !! Returns a function pointer to a unary transformation based on a string name

        character(*)                                  ,intent(in)    :: name
        procedure(realArrayFunctionGenParam) ,pointer ,intent(inout) :: funcPointer

    end subroutine associateFunctionPointer
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryLog(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for natural log. String name "log".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryLog
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryExp(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for exp. String name "exp".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryExp
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unarySin(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for sin. String name "sin".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unarySin
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryCos(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for cos. String name "cos".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryCos
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryAbs(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for absolute value. String name "abs".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryAbs
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unarySign(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for sign function. String name "sign".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unarySign
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryTan(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for tan. String name "tan".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryTan
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryAsin(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for arcsin. String name "asin".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryAsin
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryAcos(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for arccos. String name "acos".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryAcos
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryAtan(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for arctan. String name "atan".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryAtan
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryErf(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for the error function. String name "erf".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryErf
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryErfc(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for the complementary error function. String name "erfc".

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryErfc
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryPowReal(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for the raising to a real power. String name "rpow". Uses realParams(1) as the power.

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryPowReal
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryPowInt(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper for the raising to an integer power. String name "ipow". Uses intParams(1) as the power.

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

    end function unaryPowInt
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function unaryShift(input,realParams,intParams,logicalParams) result(output)
        !! Unary wrapper cyclically permuting input. String name "shift". Uses intParams(1) as the shift amount.

        real(rk)               ,dimension(:) ,intent(in) :: input 
        real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
        integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
        logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
        real(rk) ,allocatable ,dimension(:)              :: output

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

    end function unaryExpand
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module unary_transforms
!-----------------------------------------------------------------------------------------------------------------------------------