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
submodule (textbook_class) textbook_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the textbook class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initTextbook(this) 
    !! Textbook object initialization

    class(Textbook)           ,intent(inout)  :: this

    call this%makeDefined()

    allocate(this%derivations(0))
    allocate(this%derivationNames(0))

    allocate(this%matDerivations(0))
    allocate(this%matDerivationNames(0))
    
end subroutine initTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addDerivation(this,deriv,name) 
    !! Add derivation object to textbook

    class(Textbook)          ,intent(inout)  :: this
    class(Derivation)        ,intent(in)     :: deriv
    character(*)             ,intent(in)     :: name

    type(DerivationContainer) ,allocatable ,dimension(:) :: tempDeriv

    integer(ik) :: i

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add derivation to undefined textbook")
        call assertPure(deriv%isDefined(),"Attempted to add undefined derivation to textbook")

        call assertPure(name /= keyNone,"The name "//keyNone//" for derivations is reserved and cannot be used in textbooks")
        call assertPure(.not. this%isDerivationRegistered(name),"Attempted to add derivation to textbook with same name as an &
        &already added derivation")

    end if

    allocate(tempDeriv(size(this%derivations)+1))

    do i = 1, size(this%derivations)
        allocate(tempDeriv(i)%entry,source=this%derivations(i)%entry)
    end do

    allocate(tempDeriv(size(this%derivations)+1)%entry,source=deriv)

    call move_alloc(tempDeriv,this%derivations)
    this%derivationNames = [this%derivationNames,StringArray(name)] 

end subroutine addDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isDerivationRegistered(this,name) result(reg)
    !! Check whether derivation with given name is registered in the textbook

    class(Textbook)          ,intent(in)  :: this
    character(*)             ,intent(in)  :: name
    logical                               :: reg

    integer(ik) :: i

    if (assertions) call assertPure(this%isDefined(),"Attempted to get derivation name registration status from undefined textbook")

    reg = .false.
    do i = 1,size(this%derivationNames)
        if (this%derivationNames(i)%string == name) then 
            reg = .true. 
            exit 
        end if
    end do

end function isDerivationRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine copyDerivation(this,name,deriv) 
    !! Copy derivation with given name into passed deriv object, overwriting any existing derivation

    class(Textbook)                ,intent(in)    :: this
    character(*)                   ,intent(in)    :: name
    class(Derivation) ,allocatable ,intent(inout) :: deriv

    logical :: found
    integer(ik) :: i,ind

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to copy derivation from undefined textbook")
    end if

    found = .false.

    do i = 1,size(this%derivationNames)
        if (this%derivationNames(i)%string == name) then 
            found = .true. 
            ind = i 
            exit 
        end if
    end do

    if (assertions) call assertPure(found,"Attempted to copy derivation name "//name// " not registered in textbook")

    if (allocated(deriv)) deallocate(deriv)
    allocate(deriv,source=this%derivations(ind)%entry)

end subroutine copyDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine addMatDerivation(this,deriv,name) 
    !! Add matrix derivation object to textbook

    class(Textbook)          ,intent(inout)  :: this
    class(MatDerivation)     ,intent(in)     :: deriv
    character(*)             ,intent(in)     :: name

    type(MatDerivationContainer) ,allocatable ,dimension(:) :: tempDeriv

    integer(ik) :: i

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to add matrix derivation to undefined textbook")
        call assertPure(deriv%isDefined(),"Attempted to add undefined matrix derivation to textbook")

        call assertPure(name /= keyNone,"The name "//keyNone//" for derivations is reserved and cannot be used in textbooks")

        call assertPure(.not. this%isDerivationRegistered(name),"Attempted to add matrix derivation to textbook with same name as &
        &an already added matrix derivation")

    end if

    allocate(tempDeriv(size(this%matDerivations)+1))

    do i = 1, size(this%matDerivations)
        allocate(tempDeriv(i)%entry,source=this%matDerivations(i)%entry)
    end do

    allocate(tempDeriv(size(this%matDerivations)+1)%entry,source=deriv)

    call move_alloc(tempDeriv,this%matDerivations)
    this%matDerivationNames = [this%matDerivationNames,StringArray(name)] 


end subroutine addMatDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isMatDerivationRegistered(this,name) result(reg)
    !! Check whether matrix derivation with given name is registered in the textbook

    class(Textbook)          ,intent(in)  :: this
    character(*)             ,intent(in)  :: name
    logical                               :: reg

    integer(ik) :: i

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to get matrix derivation name registration status from undefined textbook")

    reg = .false.
    do i = 1,size(this%matDerivationNames)
        if (this%matDerivationNames(i)%string == name) then 
            reg = .true. 
            exit 
        end if
    end do

end function isMatDerivationRegistered
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine copyMatDerivation(this,name,deriv) 
    !! Copy matrix derivation with given name into passed deriv object, overwriting any existing derivation

    class(Textbook)                   ,intent(in)    :: this
    character(*)                      ,intent(in)    :: name
    class(MatDerivation) ,allocatable ,intent(inout) :: deriv

    logical :: found
    integer(ik) :: i,ind

    if (assertions) then 
        call assertPure(this%isDefined(),"Attempted to copy matrix derivation from undefined textbook")
    end if

    found = .false.

    do i = 1,size(this%matDerivationNames)
        if (this%matDerivationNames(i)%string == name) then 
            found = .true. 
            ind = i 
            exit 
        end if
    end do

    if (assertions) call assertPure(found,"Attempted to copy matrix derivation name "//name// "not registered in textbook")

    if (allocated(deriv)) deallocate(deriv)
    allocate(deriv,source=this%matDerivations(ind)%entry)

end subroutine copyMatDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule textbook_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
