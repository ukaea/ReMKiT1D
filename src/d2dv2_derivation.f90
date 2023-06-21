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
module d2dv2_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that calculates velocity space second order derivatives

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

    type ,public ,extends(Derivation) :: D2DV2Derivation
        !! Calculates v_o*d(v_i df/dv)/dv, where v_o and v_i are velocity vectors, and f is a distribution. Can return either a distribution
        !! or a single harmonic. Assumes that df/dv at v->inf is 0.

        integer(ik)                          ,private :: numH !! Total number of harmonics 
        real(rk) ,allocatable ,dimension(:)  ,private :: dvCopy !! Local copy of velocity grid widths for easier calculations
        real(rk) ,allocatable ,dimension(:)  ,private :: dvPlus !! Local copy of v_n+1 - v_n

        type(RealArray) ,allocatable ,dimension(:)  ,private :: outerV !! Outer velocity vector corresponding to cell centres for each included harmonic. Defaults to ones.
        type(RealArray) ,allocatable ,dimension(:)  ,private :: innerV !! Inner velocity vector corresponding to right cell boundaries for each included harmonic. Defaults to ones.

        type(RealArray) ,allocatable ,dimension(:)  ,private :: vidfdvAtZero !! Extrapolation of v_i*df/dv at zero in the form A1*f(v1)+A2*f(v2) where A's correspond to included harmonics (default = 0)

        integer(ik) ,allocatable ,dimension(:) ,private :: includedHs !! Harmonics included in the output (either all or just the targeted harmonics). Defaults to all.

        contains

        procedure ,public :: init => initD2DV2Derivation

        procedure ,public :: calculate => calculateD2DV2

    end type D2DV2Derivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initD2DV2Derivation(this,vSpaceObj,outerV,innerV,vidfdvAtZero,targetH)
        !! Initialize second order velocity derivative derivation

        class(D2DV2Derivation)                  ,intent(inout) :: this
        type(VSpace)                            ,intent(in)    :: vSpaceObj
        type(RealArray) ,optional ,dimension(:) ,intent(in)    :: outerV !! Outer velocity vector corresponding to cell centres for each included harmonic. Defaults to ones.
        type(RealArray) ,optional ,dimension(:) ,intent(in)    :: innerV !! Inner velocity vector corresponding to right cell boundaries for each included harmonic. Defaults to ones.
        type(RealArray) ,optional ,dimension(:) ,intent(in)    :: vidfdvAtZero !! Extrapolation of v_i*df/dv at zero in the form A1*f(v1)+A2*f(v2) where A's correspond to included harmonics (default = 0)
        integer(ik) ,optional                   ,intent(in)    :: targetH !! Harmonic to take derivative of. If not present will return full distribution. 

    end subroutine initD2DV2Derivation  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateD2DV2(this,inputArray,indices) result(output)

        class(D2DV2Derivation)             ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateD2DV2
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module d2dv2_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 