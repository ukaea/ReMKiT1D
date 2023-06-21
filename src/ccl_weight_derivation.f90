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
module ccl_weight_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses Chang-Cooper-Langdon interpolation weight derivation

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use physics_functions
    use derivation_abstract_class   ,only: Derivation
    use v_space_class               ,only: VSpace
    use ieee_arithmetic

    implicit none
    private

    type ,public ,extends(Derivation) :: CCLWeightDerivation
        !! Calculates the Chang-Cooper-Langdon interpolation weight for f0 Coulomb collisions assuming that the two passed variables
        !! are single harmonic C and D integrals (see 4.5 in Mijin thesis)

        real(rk) ,allocatable ,dimension(:) ,private :: dvPlus !! v_{n+1} - v_n vector 

        integer(ik) ,private :: numV
        
        contains

        procedure ,public :: init => initCCLWeights

        procedure ,public :: calculate => calculateCCLWeights

    end type CCLWeightDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCCLWeights(this,vSpaceObj)
        !! Initialize Chang-Cooper-Langdon interpolation weight derivation

        class(CCLWeightDerivation)     ,intent(inout) :: this
        type(VSpace)                   ,intent(in)    :: vSpaceObj

    end subroutine initCCLWeights  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateCCLWeights(this,inputArray,indices) result(output)

        class(CCLWeightDerivation)         ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateCCLWeights
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module ccl_weight_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 