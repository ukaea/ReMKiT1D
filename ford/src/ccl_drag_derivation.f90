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
module ccl_drag_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses Chang-Cooper-Langdon drag coefficient derivation

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

    type ,public ,extends(Derivation) :: CCLDragDerivation
        !! Calculates the Chang-Cooper-Langdon drag coefficient for f0 Coulomb collisions assuming that the passed variable is
        !! the single harmonic f0

        real(rk) ,allocatable ,dimension(:) ,private :: v2dv !! v^2 dv vector used for integration

        integer(ik) ,private :: numV
        
        contains

        procedure ,public :: init => initCCLDrag

        procedure ,public :: calculate => calculateCCLDrag

    end type CCLDragDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCCLDrag(this,vSpaceObj)
        !! Initialize Chang-Cooper-Langdon drag coefficient derivation

        class(CCLDragDerivation)     ,intent(inout) :: this
        type(VSpace)                 ,intent(in)    :: vSpaceObj

    end subroutine initCCLDrag  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateCCLDrag(this,inputArray,indices) result(output)

        class(CCLDragDerivation)           ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateCCLDrag
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module ccl_drag_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 