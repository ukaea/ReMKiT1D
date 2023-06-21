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
module ddv_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that calculates velocity space first order derivatives

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

    type ,public ,extends(Derivation) :: DDVDerivation
        !! Calculates v_o*d(v_i f)/dv, where v_o and v_i are velocity vectors, and f is a distribution. Can return either a distribution
        !! or a single harmonic. Assumes that f at v->inf is 0.

        integer(ik)                          ,private :: numH !! Total number of harmonics 
        real(rk) ,allocatable ,dimension(:)  ,private :: dvCopy !! Local copy of velocity grid widths for easier calculations
        real(rk) ,allocatable ,dimension(:)  ,private :: vInterpCopy !! Local copy of velocity interpolation factors for easier calculations

        type(RealArray) ,allocatable ,dimension(:)  ,private :: outerV !! Outer velocity vector corresponding to cell centres for each included harmonic. Defaults to ones.
        type(RealArray) ,allocatable ,dimension(:)  ,private :: innerV !! Inner velocity vector corresponding to right cell boundaries for each included harmonic. Defaults to ones.

        type(RealArray) ,allocatable ,dimension(:)  ,private :: vifAtZero !! Extrapolation of v_i*f at zero in the form A1*f(v1)+A2*f(v2) where A's are given for each included harmonic(default = 0)

        integer(ik) ,allocatable ,dimension(:) ,private :: includedHs !! Harmonics included in the output (either all or just the targeted harmonics). Defaults to all.
        contains

        procedure ,public :: init => initDDVDerivation

        procedure ,public :: calculate => calculateDDV

    end type DDVDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initDDVDerivation(this,vSpaceObj,outerV,innerV,vifAtZero,targetH)
        !! Initialize first order velocity derivative derivation

        class(DDVDerivation)                    ,intent(inout) :: this
        type(VSpace)                            ,intent(in)    :: vSpaceObj
        type(RealArray) ,optional ,dimension(:) ,intent(in)    :: outerV !! Outer velocity vector corresponding to cell centres for each included harmonic. Defaults to ones.
        type(RealArray) ,optional ,dimension(:) ,intent(in)    :: innerV !! Inner velocity vector corresponding to right cell boundaries for each included harmonic. Defaults to ones.
        type(RealArray) ,optional ,dimension(:) ,intent(in)    :: vifAtZero !! Extrapolation of v_i*f at zero in the form A1*f(v1)+A2*f(v2) where A's are given for each included harmonic(default = 0)
        integer(ik) ,optional                   ,intent(in)    :: targetH !! Harmonic to take derivative of. If not present will return full distribution/include all harmonics. 

    end subroutine initDDVDerivation  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateDDV(this,inputArray,indices) result(output)

        class(DDVDerivation)               ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateDDV
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module ddv_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 