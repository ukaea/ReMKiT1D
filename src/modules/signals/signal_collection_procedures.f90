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
submodule (signal_collection_class) signal_collection_procedure
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the signal collection class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initSignalCollection(this) 
    !! SignalCollection object initialization

    class(SignalCollection)           ,intent(inout)  :: this

    call this%makeDefined()

    allocate(this%signals(0))
    allocate(this%signalNames(0))

end subroutine initSignalCollection
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addSignal(this,sig,name) 
    !! Add named signal to collection

    class(SignalCollection)  ,intent(inout)  :: this
    class(TimeSignal)        ,intent(in)     :: sig
    character(*)             ,intent(in)     :: name

    type(SignalContainer) ,allocatable ,dimension(:) :: tempSignals 

    integer(ik) :: i

    if (assertions) then 

        call assertPure(this%isDefined(),"addSignal called from undefined SignalCollection")
        call assertPure(.not. this%isSignalRegistered(name),"TimeSignal name passed to addSignal is already registered")

    end if

    allocate(tempSignals(size(this%signals)+1))

    do i = 1, size(this%signals)
        allocate(tempSignals(i)%entry,source=this%signals(i)%entry)
    end do

    allocate(tempSignals(size(this%signals)+1)%entry,source=sig)
    call move_alloc(tempSignals,this%signals)

    this%signalNames = [this%signalNames,StringArray(name)]

end subroutine addSignal
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isSignalRegistered(this,name) result(reg)
    !! Check whether signal with given name is registered in the collection

    class(SignalCollection)  ,intent(in)  :: this
    character(*)             ,intent(in)  :: name
    logical                               :: reg

    integer(ik) :: i

    if (assertions) &
    call assertPure(this%isDefined(),"Attempted to get signal name registration status from undefined SignalCollection")

    reg = .false.
    do i = 1,size(this%signalNames)
        if (this%signalNames(i)%string == name) then 
            reg = .true. 
            exit 
        end if
    end do

end function isSignalRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine copySignal(this,name,sig) 
    !! Copy signal with given name into passed sig object, overwriting any existing allocation

    class(SignalCollection)        ,intent(in)    :: this
    character(*)                   ,intent(in)    :: name
    class(TimeSignal) ,allocatable ,intent(inout) :: sig

    logical :: found
    integer(ik) :: i,ind

    if (assertions) &
    call assertPure(this%isDefined(),"Attempted to copy signal from undefined SignalCollection")

    found = .false.

    do i = 1,size(this%signalNames)
        if (this%signalNames(i)%string == name) then 
            found = .true. 
            ind = i 
            exit 
        end if
    end do

    if (assertions) call assertPure(found,"Attempted to copy signal name not registered in SignalCollection")

    if (allocated(sig)) deallocate(sig)
    allocate(sig,source=this%signals(ind)%entry)

end subroutine copySignal
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule signal_collection_procedure
!-----------------------------------------------------------------------------------------------------------------------------------
