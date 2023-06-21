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
module cold_ion_ij_int_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses Shkarofsky I and J integral for a drifting delta function derivation object

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use physics_functions
    use derivation_abstract_class   ,only: Derivation
    use grid_class                  ,only: Grid

    implicit none
    private

    type ,public ,extends(Derivation) :: ColdIonIJIntDerivation
        !! Calculates the Shkarofsky I or J integrals for cold drifting ions (normalized to ion density). The passed variable is
        !! assumed to be the ion flow velocity 

        real(rk) ,allocatable ,dimension(:) ,private :: vGridCopy !! Local copy of velocity grid

        integer(ik) ,private :: ind !! Integral index
        integer(ik) ,allocatable ,dimension(:) ,private :: lGridCopy !! Local copy of l harmonic grid
        
        logical ,private :: isJInt !! True if the calculated integral is the J integral instead of the I integral.
        
        contains

        procedure ,public :: init => initColdIJInt

        procedure ,public :: calculate => calculateColdIJInt

    end type ColdIonIJIntDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initColdIJInt(this,gridObj,ind,isJInt)
        !! Initialize col dion Shkarofsky I/J integral derivation derivation

        class(ColdIonIJIntDerivation)     ,intent(inout) :: this
        type(Grid)                        ,intent(in)    :: gridObj
        integer(ik)                       ,intent(in)    :: ind !! Index of integral
        logical ,optional                 ,intent(in)    :: isJInt !! If true the lower triangular J integral is calculated instead of the I integral. Defaults to false.

    end subroutine initColdIJInt  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateColdIJInt(this,inputArray,indices) result(output)

        class(ColdIonIJIntDerivation)      ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateColdIJInt
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module cold_ion_ij_int_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 