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
module vel_contraction_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that takes a general contraction in velocity space of a distribution function harmonic

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation
    use v_space_class               ,only: VSpace

    implicit none
    private

    type ,public ,extends(Derivation) :: VelContracDerivation
        !! Derivation that calculates the contraction of velocity vector g with h'th harmonic of single passed distribution variable

        integer(ik)                         ,private :: h !! Harmonic to contract with
        real(rk) ,allocatable ,dimension(:) ,private :: g !! Velocity space vector     
        integer(ik)                         ,private :: numH !! Local copy of number of harmonics
        
        contains

        procedure ,public :: init => initContracDeriv

        procedure ,public :: calculate => calculateContracDeriv

    end type VelContracDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initContracDeriv(this,h,g,refVSpace,expH)
        !! Initialize velocity contraction derivation object
    
        class(VelContracDerivation)         ,intent(inout) :: this
        integer(ik)                         ,intent(in)    :: h !! Harmonic to contract with
        real(rk)           ,dimension(:)    ,intent(in)    :: g !! Velocity space vector 
        type(VSpace)                        ,intent(in)    :: refVSpace !! Reference velocity space
        integer(ik)   ,optional             ,intent(in)    :: expH !! Expected number of harmonics (defaults to total number of harmonics)

    end subroutine initContracDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateContracDeriv(this,inputArray,indices) result(output)

        class(VelContracDerivation)        ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateContracDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module vel_contraction_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 