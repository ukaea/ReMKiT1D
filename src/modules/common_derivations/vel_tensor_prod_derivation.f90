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
module vel_tensor_prod_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that takes calculates the flattened velocity vector and x vector tensor product, resulting in a single harmonic variable

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions, assertionLvl
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use derivation_abstract_class   ,only: Derivation
    use v_space_class               ,only: VSpace
    use support_functions

    implicit none
    private

    type ,public ,extends(Derivation) :: VelTProdDerivation
        !! Derivation that calculates the flattened tensor product of vector v with a fluid variable, optionally shifting v by a second fluid variable
        !! and raising the shifted values to a power

        real(rk) ,allocatable ,dimension(:) ,private :: v !! Velocity space vector     
        real(rk)                            ,private :: power !! Power to raise the shifted velocity values to
        
        contains

        procedure ,public :: init => initVTProdDeriv

        procedure ,public :: calculate => calculateVTProdDeriv

    end type VelTProdDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVTProdDeriv(this,refVSpace,velVec,power)
        !! Initialize velocity tensor product derivation object

        class(VelTProdDerivation)           ,intent(inout) :: this
        type(VSpace)                        ,intent(in)    :: refVSpace !! Reference velocity space
        real(rk) ,optional ,dimension(:)    ,intent(in)    :: velVec !! Optional velocity space vector. Defaults to velcity grid values
        real(rk) ,optional                  ,intent(in)    :: power !! Optional power to raise the shifted velocity vector to. Defaults to 1

    end subroutine initVTProdDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateVTProdDeriv(this,inputArray,indices) result(output)

        class(VelTProdDerivation)          ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateVTProdDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module vel_tensor_prod_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 