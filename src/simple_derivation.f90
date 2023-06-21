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
module simple_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses simple derivation class handling derivations of the form const*product(variables**powers)

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation

    implicit none
    private

    type ,public ,extends(Derivation) :: SimpleDerivation
        !! Simple derivation object assuming all variables are confroming. Just multiplies the input array entries and raises them to
        !! set powers. Optionally multiplies the result with a constant. 

        real(rk) ,allocatable ,dimension(:) ,private :: varPowers !! Powers corresponding to each variable - indices in calculate must conform to this
        real(rk)                            ,private :: multConst !! Multiplicative constant - default 1

        contains

        procedure ,public :: init => initSimpleDeriv

        procedure ,public :: calculate => calculateSimple

    end type SimpleDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initSimpleDeriv(this,varPowers,multConst)
        !! Initialize simple derivation object

        class(SimpleDerivation)   ,intent(inout) :: this
        real(rk) ,dimension(:)    ,intent(in)    :: varPowers
        real(rk) ,optional        ,intent(in)    :: multConst

    end subroutine initSimpleDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateSimple(this,inputArray,indices) result(output)

        class(SimpleDerivation)            ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateSimple
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module simple_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 