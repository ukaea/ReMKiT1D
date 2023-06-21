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
module ccl_diff_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses Chang-Cooper-Langdon diffusion coefficient derivation

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use physics_functions
    use derivation_abstract_class   ,only: Derivation
    use v_space_class               ,only: VSpace
    use physical_constants

    implicit none
    private

    type ,public ,extends(Derivation) :: CCLDiffDerivation
        !! Calculates the Chang-Cooper-Langdon diffusion coefficient for f0 Coulomb collisions assuming that the passed variables are
        !! the single harmonic f0 and the ccl interpolation weight

        real(rk) ,allocatable ,dimension(:) ,private :: v2dv !! v^2 dv vector used for outer integration
        real(rk) ,allocatable ,dimension(:) ,private :: vdvPlus !! v_{m+1/2}(v_{m+1}-v_m) vector used for inner integration
        real(rk) ,allocatable ,dimension(:) ,private :: vStar !! (v_m+v_{m+1})/2 vector 

        integer(ik) ,private :: numV
        
        contains

        procedure ,public :: init => initCCLDiff

        procedure ,public :: calculate => calculateCCLDiff

    end type CCLDiffDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCCLDiff(this,vSpaceObj)
        !! Initialize Chang-Cooper-Langdon diffusion coefficient derivation

        class(CCLDiffDerivation)     ,intent(inout) :: this
        type(VSpace)                 ,intent(in)    :: vSpaceObj

    end subroutine initCCLDiff  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateCCLDiff(this,inputArray,indices) result(output)

        class(CCLDiffDerivation)           ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateCCLDiff
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module ccl_diff_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 