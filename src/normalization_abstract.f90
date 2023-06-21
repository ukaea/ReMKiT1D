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
module normalization_abstract_class
    !! author: Stefan Mijin
    !! 
    !! Houses abstract object responsible for calculating and determining normalization constants

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use support_types                         ,only: NamedReal ,StringArray

    implicit none
    private

    type ,public ,extends(Object), abstract :: Normalization
        !! Abstract normalization object providing access routines

        type(NamedReal) ,allocatable ,dimension(:) ,private :: normalizationVals

        contains

        procedure ,public :: getNormalizationValue
        procedure ,public :: getCustomNormalization
        procedure ,public :: setNormalizationVals

    end type Normalization
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine setNormalizationVals(this,normVals)
        !! Set normalization values array of normalization object

        class(Normalization)          ,intent(inout)  :: this
        type(NamedReal) ,dimension(:) ,intent(in)     :: normVals

    end subroutine setNormalizationVals  
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNormalizationValue(this,name) result(val)
        !! Get normalization value associated with passed name

        class(Normalization)      ,intent(in) :: this
        character(*)              ,intent(in) :: name
        real(rk)                              :: val

    end function getNormalizationValue   
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getCustomNormalization(this,names,powers,multConst) result(val)
        !! Get normalization value calculated as the product of normalization values with passed names raised to corresponding 
        !! powers and optionally multiplied by a constant

        class(Normalization)              ,intent(in) :: this
        type(StringArray) ,dimension(:)   ,intent(in) :: names
        real(rk)          ,dimension(:)   ,intent(in) :: powers
        real(rk) ,optional                ,intent(in) :: multConst
        real(rk)                                      :: val

    end function getCustomNormalization   
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module normalization_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 