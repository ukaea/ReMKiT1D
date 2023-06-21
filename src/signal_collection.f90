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
module signal_collection_class
    !! author: Stefan Mijin 
    !!
    !! Houses class containing named signal wrappers 

    use data_kinds                           ,only: rk, ik
    use runtime_constants                    ,only: debugging, assertions
    use god_objects                          ,only: Object
    use assertion_utility                    ,only: assert, assertIdentical, assertPure
    use signal_abstract_class                ,only: TimeSignal
    use support_types

    implicit none
    private

    type ,private :: SignalContainer 
        !! Container object for signal types
        class(TimeSignal) ,allocatable :: entry
    end type

    type ,public ,extends(Object) :: SignalCollection
        !! Object storing named signals

        type(SignalContainer) ,allocatable ,dimension(:) ,private :: signals
        type(StringArray)     ,allocatable ,dimension(:) ,private :: signalNames

        contains

        procedure ,public :: addSignal
        procedure ,public :: isSignalRegistered
        procedure ,public :: copySignal

        procedure ,public :: init => initSignalCollection

    end type SignalCollection
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initSignalCollection(this) 
        !! SignalCollection object initialization

        class(SignalCollection)           ,intent(inout)  :: this

    end subroutine initSignalCollection
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addSignal(this,sig,name) 
       !! Add named signal to collection

        class(SignalCollection)  ,intent(inout)  :: this
        class(TimeSignal)            ,intent(in)     :: sig
        character(*)             ,intent(in)     :: name

    end subroutine addSignal
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function isSignalRegistered(this,name) result(reg)
        !! Check whether signal with given name is registered in the collection

        class(SignalCollection)  ,intent(in)  :: this
        character(*)             ,intent(in)  :: name
        logical                               :: reg

    end function isSignalRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine copySignal(this,name,sig) 
        !! Copy signal with given name into passed sig object, overwriting any existing allocation

        class(SignalCollection)    ,intent(in)    :: this
        character(*)               ,intent(in)    :: name
        class(TimeSignal) ,allocatable ,intent(inout) :: sig

    end subroutine copySignal
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module signal_collection_class
!-----------------------------------------------------------------------------------------------------------------------------------
 