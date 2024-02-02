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
module polynomial_fun_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation taking in a single variable and uses it in a polynomial function

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation

    implicit none
    private

    type ,public ,extends(Derivation) :: PolyFunDeriv
        !! Derivation taking in a single variable v and returning c0 + sum(c_i*v**p_i), where p are powers and c polynomial coefficents,
        !! Can optionally take as many conforming variables as there are coefficients and return c0 + sum(c_i*v_i**p_i)

        real(rk) ,allocatable ,dimension(:) ,private :: polyPowers !! Polynomial powers
        real(rk) ,allocatable ,dimension(:) ,private :: polyCoeffs !! Polynomial coefficients
        real(rk)                            ,private :: constCoeff !! Constant coefficient - default 0

        contains

        procedure ,public :: init => initPolyFunDeriv

        procedure ,public :: calculate => calculatePolyFun

    end type PolyFunDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initPolyFunDeriv(this,polyPowers,polyCoeffs,constCoeff)
        !! Initialize polynomial function derivation object
    
        class(PolyFunDeriv)       ,intent(inout) :: this
        real(rk) ,dimension(:)    ,intent(in)    :: polyPowers
        real(rk) ,dimension(:)    ,intent(in)    :: polyCoeffs
        real(rk) ,optional        ,intent(in)    :: constCoeff

    end subroutine initPolyFunDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculatePolyFun(this,inputArray,indices) result(output)

        class(PolyFunDeriv)                ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculatePolyFun
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module polynomial_fun_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 