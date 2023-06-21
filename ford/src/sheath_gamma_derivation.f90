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
module sheath_gamma_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation returning electron sheath heat transmission coefficient with given electron and ion temperatures

    use data_kinds                   ,only: rk ,ik
    use runtime_constants            ,only: debugging, assertions
    use assertion_utility            ,only: assert, assertIdentical, assertPure
    use god_objects                  ,only: Object
    use support_types                ,only: RealArray 
    use derivation_abstract_class    ,only: Derivation
    use physics_functions

    implicit none
    private

    type ,public ,extends(Derivation) :: ElectronSheathGammaDerivation
        !! Returns the sheath heat transmission coefficient for electrons depending on electron and ion temperatures in the boundary
        !! cell. Assumes the indices are in Te,Ti order

        real(rk) ,private :: massRatio !! e-i mass ratio
        integer(ik) ,private :: boundaryIndex !! Index of the local cell used for this boundary calculation

        contains

        procedure ,public :: init => initElectronSheathGamma

        procedure ,public :: calculate => calculateElectronGamma

    end type ElectronSheathGammaDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initElectronSheathGamma(this,massRatio,boundaryIndex)
        !! Initialize electron sheath gamma derivation object

        class(ElectronSheathGammaDerivation)   ,intent(inout) :: this
        real(rk)                               ,intent(in)    :: massRatio
        integer(ik)                            ,intent(in)    :: boundaryIndex

    end subroutine initElectronSheathGamma  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateElectronGamma(this,inputArray,indices) result(output)

        class(ElectronSheathGammaDerivation),intent(inout)    :: this 
        type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:)  ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                 :: output

    end function calculateElectronGamma
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module sheath_gamma_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 