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
module derivation_abstract_class
    !! author: Stefan Mijin 
    !!
    !! Houses abstract Derivation object definitions used to build rules for transforming evolved into derived variables

    use data_kinds                  ,only: rk ,ik
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 

    implicit none
    private

    type ,public :: DerivationContainer 
        class(Derivation) ,allocatable :: entry 
    end type DerivationContainer

    type ,public ,extends(Object), abstract :: Derivation
        !! Abstract derivation object defining the interface for calculating derived variables and data

        contains

        procedure(calculation) ,deferred :: calculate

    end type Derivation
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        function calculation(this,inputArray,indices) result(output)

            import :: ik ,rk ,RealArray ,Derivation

            class(Derivation)                  ,intent(inout)    :: this 
            type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
            integer(ik)           ,dimension(:) ,intent(in)    :: indices           
            real(rk) ,allocatable ,dimension(:)                :: output

        end function calculation
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module derivation_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 