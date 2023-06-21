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
module moment_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses moment derivation class handling derivations of the form const*product(fluidVariables**powers)* m-th moment of h-th
    !! harmonic of f * g, where g is an optional constant velocity vector 

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation
    use v_space_class               ,only: VSpace

    implicit none
    private

    type ,public ,extends(Derivation) :: MomentDerivation
        !! Moment derivation object. First step in derivation is taking the mOrder-th moment of the h-th harmonic of the first passed
        !! variable (assumes that it is a distribution), optionally multiplied by constant velocity space vector g. The results is then
        !! optionally multiplied by a constant and a product of fluid variables raised to powers (as in SimpleDerivation).  

        real(rk) ,allocatable ,dimension(:) ,private :: varPowers !! Powers corresponding to each fluid variable - indices in calculate must conform size of this + 1
        real(rk)                            ,private :: multConst !! Multiplicative constant - default 1
        integer(ik)                         ,private :: mOrder !! Moment order 
        integer(ik)                         ,private :: h !! Harmonic to take moment of
        real(rk) ,allocatable ,dimension(:) ,private :: g !! Optional velocity space vector
        type(VSpace)                        ,private :: vSpaceCopy !! Copy of velocity space object used to take the moment 
                                                                   !! (NOTE: this might cause performance issues if many momentDerivations are present at once)        
        
        contains

        procedure ,public :: init => initMomentDeriv

        procedure ,public :: calculate => calculateMomentDeriv

    end type MomentDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initMomentDeriv(this,mOrder,h,refVSpace,varPowers,g,multConst)
        !! Initialize moment derivation object

        class(MomentDerivation)             ,intent(inout) :: this
        integer(ik)                         ,intent(in)    :: mOrder !! Moment order 
        integer(ik)                         ,intent(in)    :: h !! Harmonic to take moment of
        type(VSpace)                        ,intent(in)    :: refVSpace !! Velocity space object to make copy of
        real(rk) ,optional ,dimension(:)    ,intent(in)    :: varPowers !! Optional fluid variable powers
        real(rk) ,optional ,dimension(:)    ,intent(in)    :: g !! Optional velocity space vector
        real(rk) ,optional                  ,intent(in)    :: multConst !! Optional multiplicative constant - default 1

    end subroutine initMomentDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateMomentDeriv(this,inputArray,indices) result(output)

        class(MomentDerivation)            ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateMomentDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module moment_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 