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
module ij_int_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses Shkarofsky I and J integral derivation object

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use physics_functions
    use derivation_abstract_class   ,only: Derivation
    use v_space_class               ,only: VSpace
    use sparse_row_data_class       ,only: SparseRowData

    implicit none
    private

    type ,public ,extends(Derivation) :: IJIntDerivation
        !! Calculates the Shkarofsky I or J integrals for each passed harmonic. 

        type(SparseRowData) ,private :: ijMatBuffer !! Sparse buffer used to calculate the integrals

        integer(ik) ,private :: numV
        
        contains

        procedure ,public :: init => initIJInt

        procedure ,public :: calculate => calculateIJInt

    end type IJIntDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initIJInt(this,vSpaceObj,ind,isJInt)
        !! Initialize Shkarofsky I/J integral derivation derivation

        class(IJIntDerivation)     ,intent(inout) :: this
        type(VSpace)               ,intent(in)    :: vSpaceObj
        integer(ik)                ,intent(in)    :: ind !! Index of indegral
        logical ,optional          ,intent(in)    :: isJInt !! If true the lower triangular J integral is calculated instead of the I integral. Defaults to false.

    end subroutine initIJInt  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateIJInt(this,inputArray,indices) result(output)

        class(IJIntDerivation)             ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateIJInt
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module ij_int_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 