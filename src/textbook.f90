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
module textbook_class
    !! author: Stefan Mijin 
    !!
    !! Houses class containing named derivation rules 

    use data_kinds                           ,only: rk, ik
    use runtime_constants                    ,only: debugging, assertions
    use god_objects                          ,only: Object
    use assertion_utility                    ,only: assert, assertIdentical, assertPure
    use derivation_abstract_class            ,only: Derivation, DerivationContainer
    use mat_derivation_abstract_class        ,only: MatDerivation, MatDerivationContainer
    use support_types
    use key_names

    implicit none
    private

    type ,public ,extends(Object) :: Textbook
        !! Object storing named derivations 

        type(DerivationContainer) ,allocatable ,dimension(:) ,private :: derivations
        type(StringArray)         ,allocatable ,dimension(:) ,private :: derivationNames

        type(MatDerivationContainer) ,allocatable ,dimension(:) ,private :: matDerivations
        type(StringArray)            ,allocatable ,dimension(:) ,private :: matDerivationNames

        contains

        procedure ,public :: addDerivation 
        procedure ,public :: isDerivationRegistered
        procedure ,public :: copyDerivation

        procedure ,public :: addMatDerivation 
        procedure ,public :: isMatDerivationRegistered
        procedure ,public :: copyMatDerivation

        procedure ,public :: init => initTextbook

    end type Textbook
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initTextbook(this) 
        !! Textbook object initialization

        class(Textbook)           ,intent(inout)  :: this

    end subroutine initTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addDerivation(this,deriv,name) 
       !! Add derivation object to textbook

        class(Textbook)          ,intent(inout)  :: this
        class(Derivation)        ,intent(in)     :: deriv
        character(*)             ,intent(in)     :: name

    end subroutine addDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isDerivationRegistered(this,name) result(reg)
        !! Check whether derivation with given name is registered in the textbook

        class(Textbook)          ,intent(in)  :: this
        character(*)             ,intent(in)  :: name
        logical                               :: reg

    end function isDerivationRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine copyDerivation(this,name,deriv) 
        !! Copy derivation with given name into passed deriv object, overwriting any existing derivation

        class(Textbook)                ,intent(in)    :: this
        character(*)                   ,intent(in)    :: name
        class(Derivation) ,allocatable ,intent(inout) :: deriv

    end subroutine copyDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine addMatDerivation(this,deriv,name) 
       !! Add matrix derivation object to textbook

        class(Textbook)          ,intent(inout)  :: this
        class(MatDerivation)     ,intent(in)     :: deriv
        character(*)             ,intent(in)     :: name

    end subroutine addMatDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isMatDerivationRegistered(this,name) result(reg)
        !! Check whether matrix derivation with given name is registered in the textbook

        class(Textbook)          ,intent(in)  :: this
        character(*)             ,intent(in)  :: name
        logical                               :: reg

    end function isMatDerivationRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine copyMatDerivation(this,name,deriv) 
        !! Copy matrix derivation with given name into passed deriv object, overwriting any existing derivation

        class(Textbook)                   ,intent(in)    :: this
        character(*)                      ,intent(in)    :: name
        class(MatDerivation) ,allocatable ,intent(inout) :: deriv

    end subroutine copyMatDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module textbook_class
!-----------------------------------------------------------------------------------------------------------------------------------
 