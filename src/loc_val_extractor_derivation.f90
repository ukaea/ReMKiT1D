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
module loc_val_extractor_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that extracts a single positional value for a fluid variable

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray 
    use physics_functions
    use derivation_abstract_class   ,only: Derivation
    use partition_class             ,only: Partition

    implicit none
    private

    type ,public ,extends(Derivation) :: LocValExtractorDerivation
        !! Extracts target location value from a fluid variable as a scalar variable. Expects a single input.

        logical                              ,private :: isActive

        integer(ik)                          ,private :: targetX !! x location to extract from - local indexing
        integer(ik)                          ,private :: locNumX
        contains

        procedure ,public :: init => initLocValDeriv

        procedure ,public :: calculate => calculateLocVal

    end type LocValExtractorDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initLocValDeriv(this,partObj,numProc,targetX)
        !! Initialize harmonic extractor derivation
    
        class(LocValExtractorDerivation)    ,intent(inout) :: this
        type(Partition)                     ,intent(in)    :: partObj
        integer(ik)                         ,intent(in)    :: numProc
        integer(ik)                         ,intent(in)    :: targetX !! x location to extract from - global indexing

    end subroutine initLocValDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateLocVal(this,inputArray,indices) result(output)

        class(LocValExtractorDerivation)        ,intent(inout)    :: this 
        type(RealArray)            ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)                ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable      ,dimension(:)                :: output

    end function calculateLocVal
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module loc_val_extractor_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 