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
module multiplicative_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses composite derivation class where two derivation results get multiplied together, optionally applying an elementary function (exp, log, sin, cos) to one of them and raising them to corresponding powers

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray ,IntArray
    use derivation_abstract_class   ,only: Derivation ,DerivationContainer

    implicit none
    private

    type ,public ,extends(Derivation) :: MultiplicativeDerivation
        !! Composite derivation class containing two calculation rules applied additively. Both derivations must return 
        !! same length array. Optionally applies one of several intrinsic functions to "inner" derivation result before multiplication.

        class(Derivation) ,allocatable ,private :: outerDeriv !! Outer multiplicative derivation - optional
        class(Derivation) ,allocatable ,private :: innerDeriv !! Inner multiplicative derivation - this one has function applied to it

        integer(ik) ,allocatable ,dimension(:) ,private :: outerIndices !! Subset of index vector passed to outer derivation
        integer(ik) ,allocatable ,dimension(:) ,private :: innerIndices !! Subset of index vector passed to inner derivation

        character(:) ,allocatable ,private :: innerFuncName !! Intrinsic function name optionally applied to result of inner derivation before raising to power and multiplying with outer result 
        real(rk) ,private :: innerPower !! Power to raise result of inner derivation (after applying optional function). Defaults to 1. 
        real(rk) ,private :: outerPower !! Power to raise result of outer derivation . Defaults to 1.

        contains

        procedure ,public :: init => initMultDeriv

        procedure ,public :: calculate => calculateMultiplicative

    end type MultiplicativeDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
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

    end subroutine initMultDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateMultiplicative(this,inputArray,indices) result(output)

        class(MultiplicativeDerivation)    ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateMultiplicative
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module multiplicative_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 