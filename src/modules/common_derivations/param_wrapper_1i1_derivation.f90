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
module param_wrapper_1i1_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation wrapping a real function with one real array input and one integer parameter

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation
    use basic_interfaces            ,only: realArrayFunctionIntParam

    implicit none
    private

    type ,public ,extends(Derivation) :: FunWrapperDerivation1I1
        !! Wraps single integer parameter and single real array function, passing first index variable as input, 
        !! optionally multiplying it with a constant before passing

        integer(ik)                                  ,private :: param !! Integer parameter
        real(rk)                                     ,private :: multConst !! Multiplicative constant - default 1
        real(rk)                                     ,private :: multConstInner !! Multiplicative constant before applying function - default 1
        procedure(realArrayFunctionIntParam) ,pointer ,nopass :: funPtr => null() !! Wrapped function pointer

        contains

        procedure ,public :: init => initWrapper1I1

        procedure ,public :: calculate => calculateWrapper1I1

        final             :: finalizeWrapper1I1

    end type FunWrapperDerivation1I1
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initWrapper1I1(this,fun,param,multConst,multConstInner)
        !! Initialize 1I1 function wrapper derivation

        class(FunWrapperDerivation1I1)   ,intent(inout) :: this
        procedure(realArrayFunctionIntParam)            :: fun
        integer(ik)                      ,intent(in)    :: param
        real(rk) ,optional               ,intent(in)    :: multConst
        real(rk) ,optional               ,intent(in)    :: multConstInner

    end subroutine initWrapper1I1  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateWrapper1I1(this,inputArray,indices) result(output)

        class(FunWrapperDerivation1I1)     ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateWrapper1I1
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module subroutine finalizeWrapper1I1(this) 
        !! Deallocate pointer component

        type(FunWrapperDerivation1I1) ,intent(inout) :: this

    end subroutine finalizeWrapper1I1 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module param_wrapper_1i1_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 