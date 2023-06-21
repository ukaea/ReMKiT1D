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
module maxwellian_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that calculates a Maxwellian based on temperature and density and the velocity grid

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use physics_functions
    use derivation_abstract_class   ,only: Derivation
    use v_space_class               ,only: VSpace

    implicit none
    private

    type ,public ,extends(Derivation) :: MaxwellianDerivation
        !! Calculates Maxwellian distribution. Expects temperature and density indices (in that order), and assumes velocity is normalized
        !! to electron thermal speed.

        integer(ik)                          ,private :: numH !! Total number of harmonics 
        real(rk) ,allocatable ,dimension(:)  ,private :: vGridCopy !! Local copy of velocity grid for easier calculations

        contains

        procedure ,public :: init => initMaxwellianDeriv

        procedure ,public :: calculate => calculateMaxwellian

    end type MaxwellianDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initMaxwellianDeriv(this,vSpaceObj)
        !! Initialize Maxwellian derivation

        class(MaxwellianDerivation)      ,intent(inout) :: this
        type(VSpace)                     ,intent(in)    :: vSpaceObj

    end subroutine initMaxwellianDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateMaxwellian(this,inputArray,indices) result(output)

        class(MaxwellianDerivation)        ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateMaxwellian
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module maxwellian_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 