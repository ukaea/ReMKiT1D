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
module additive_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses composite derivation class where multiple derivation results get added

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray ,IntArray
    use derivation_abstract_class   ,only: Derivation ,DerivationContainer

    implicit none
    private

    type ,public ,extends(Derivation) :: AdditiveDerivation
        !! Composite derivation class containing a number of calculation rules applied additively. All derivations must return 
        !! same length array

        type(DerivationContainer) ,allocatable ,dimension(:) ,private :: derivs !! Derivation whose results to be added
        type(IntArray)            ,allocatable ,dimension(:) ,private :: derivIndices !! Subset of indices vector passed to each derivation
        integer(ik)                                          ,private :: numAddedDerivations !! Counter for keeping track of added derivations
        integer(ik)                                          ,private :: expectedNumIndices !! Expected size of indices array - for errorChecking
        real(rk)                                             ,private :: resultPower !! Optional power to which the result is raised = default 1
        real(rk) ,allocatable ,dimension(:)                  ,private :: linearCoefficients !! Additive linear combination coefficients associated with each derivation. Defaults to ones.

        contains

        procedure ,public :: init => initAdditiveDeriv
        procedure ,public :: addDerivation

        procedure ,public :: calculate => calculateAdditive

    end type AdditiveDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initAdditiveDeriv(this,numDerivs,expectedNumIndices, resultPower, linearCoefficients)
        !! Initialize additive derivation object, allocating the expected number of derivation objects

        class(AdditiveDerivation)          ,intent(inout) :: this
        integer(ik)                        ,intent(in)    :: numDerivs
        integer(ik)                        ,intent(in)    :: expectedNumIndices
        real(rk)   ,optional               ,intent(in)    :: resultPower
        real(rk)   ,optional ,dimension(:) ,intent(in)    :: linearCoefficients

    end subroutine initAdditiveDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addDerivation(this,deriv,activeIndices)
        !! Adds derivation to additive derivation object

        class(AdditiveDerivation)   ,intent(inout) :: this
        class(Derivation)           ,intent(in)    :: deriv
        integer(ik) ,dimension(:)   ,intent(in)    :: activeIndices !! Entries of indices object passed to calculateAdditive to be passed to added derivation 

    end subroutine addDerivation  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateAdditive(this,inputArray,indices) result(output)

        class(AdditiveDerivation)            ,intent(inout) :: this
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateAdditive
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module additive_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 