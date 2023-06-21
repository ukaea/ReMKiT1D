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
module harmonic_extractor_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that extracts single harmonics from distribution variables

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

    type ,public ,extends(Derivation) :: HExtractorDerivation
        !! Extracts target harmonic from distribution variable

        integer(ik)                          ,private :: numH !! Total number of harmonics 
        integer(ik)                          ,private :: numV !! Velocity grid size

        integer(ik)                          ,private :: targetH !! Harmonic to extract
        contains

        procedure ,public :: init => initHEDerivation

        procedure ,public :: calculate => calculateHE

    end type HExtractorDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initHEDerivation(this,vSpaceObj,targetH)
        !! Initialize harmonic extractor derivation

        class(HExtractorDerivation)    ,intent(inout) :: this
        type(VSpace)                   ,intent(in)    :: vSpaceObj
        integer(ik)                    ,intent(in)    :: targetH !! Harmonic to extract

    end subroutine initHEDerivation  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateHE(this,inputArray,indices) result(output)

        class(HExtractorDerivation)        ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateHE
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module harmonic_extractor_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 