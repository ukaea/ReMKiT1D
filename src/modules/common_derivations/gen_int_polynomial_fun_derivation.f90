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
module gen_int_polynomial_fun_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation taking in multiple variables and performing a generalized integer powered polynomial calculation

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray ,IntArray ,RealArrayD2
    use derivation_abstract_class   ,only: Derivation

    implicit none
    private

    type ,public ,extends(Derivation) :: GenIntPolyFunDeriv
        !! Derivation taking in multiple variables and calculating multConst*fun(sum coef_i prod vars^powers) where prod vars^powers is
        !! shorthand for a product of multiple variables raised to corresponding powers. Basically behaves like a multiplicative derivation where the
        !! inner derivation is a massive sum of simple derivation outputs with integer powers. Implemented to allow AMJUEL rate fit derivations efficiently

        type(IntArray) ,allocatable ,dimension(:) ,private :: polyPowers !! Variable powers for each monomial
        real(rk)       ,allocatable ,dimension(:) ,private :: polyCoeffs !! Polynomial coefficients for each monomial
        real(rk)                                  ,private :: multConst !! Constant coefficient - default 1
        integer(ik)    ,allocatable ,dimension(:) ,private :: maxPowers !! Maximum powers to be precalculated for each passed variable

        character(:) ,allocatable ,private :: funcName !! Intrinsic function name optionally applied to polynomial 
        contains

        procedure ,public :: init => initGenIntPolyFunDeriv

        procedure ,public :: calculate => calculateGenIntPolyFun

    end type GenIntPolyFunDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initGenIntPolyFunDeriv(this,polyPowers,polyCoeffs,maxPowers,funcName,multConst)
        !! Initialize generalized integer powered polynomial function derivation object
    
        class(GenIntPolyFunDeriv)       ,intent(inout) :: this
        type(IntArray) ,dimension(:)    ,intent(in)    :: polyPowers
        real(rk)       ,dimension(:)    ,intent(in)    :: polyCoeffs
        integer(ik)    ,dimension(:)    ,intent(in)    :: maxPowers
        character(*) ,optional          ,intent(in)    :: funcName 
        real(rk) ,optional              ,intent(in)    :: multConst

    end subroutine initGenIntPolyFunDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateGenIntPolyFun(this,inputArray,indices) result(output)

        class(GenIntPolyFunDeriv)          ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateGenIntPolyFun
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module gen_int_polynomial_fun_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 