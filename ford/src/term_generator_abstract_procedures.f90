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
submodule (term_generator_abstract_class) term_generator_abstract_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the abstract TermGenerator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setImplicitTerms(this,impTerms) 
    !! Move allocation of impTerms to this%implicitTerms

    class(TermGenerator)                              ,intent(inout) :: this
    type(MatTermContainer) ,allocatable ,dimension(:) ,intent(inout) :: impTerms

    call move_alloc(impTerms,this%implicitTerms)

end subroutine setImplicitTerms
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setGeneralTerms(this,genTerms) 
    !! Move allocation of genTerms to this%generalTerms

    class(TermGenerator)                           ,intent(inout) :: this
    type(TermContainer) ,allocatable ,dimension(:) ,intent(inout) :: genTerms
    
    call move_alloc(genTerms,this%generalTerms)

end subroutine setGeneralTerms
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setGeneratorPrefix(this,prefix) 
    !! Set prefix for added term names

    class(TermGenerator)    ,intent(inout) :: this
    character(*)            ,intent(in)    :: prefix

    this%generatorPrefix = prefix

end subroutine setGeneratorPrefix
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumImplicitTerms(this) result(numTerms)
    !! Get size of this%implicitTerms

    class(TermGenerator)   ,intent(in)  :: this
    integer(ik)                         :: numTerms

    if (assertions) then 
        call assertPure(this%isDefined(),"getNumImplicitTerms called from undefined term generator")
        call assertPure(allocated(this%implicitTerms),"getNumImplicitTerms called from term generator with unallocated implicit &
                                                      &term container")
    end if

    numTerms = size(this%implicitTerms)

end function getNumImplicitTerms
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumGeneralTerms(this) result(numTerms)
    !! Get size of this%generalTerms

    class(TermGenerator)   ,intent(in)  :: this
    integer(ik)                         :: numTerms

    if (assertions) then 
        call assertPure(this%isDefined(),"getNumGeneralTerms called from undefined term generator")
        call assertPure(allocated(this%generalTerms),"getNumGeneralTerms called from term generator with unallocated general &
                                                      &term container")
    end if

    numTerms = size(this%generalTerms)

end function getNumGeneralTerms
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine moveTerms(this,modelObj,impTermImpGroups,impTermGenGroups,genTermGroups) 
    !! Move terms to modelObj 

    class(TermGenerator)                   ,intent(inout) :: this
    class(Model)                           ,intent(inout) :: modelObj
    type(IntArray) ,optional ,dimension(:) ,intent(in)    :: impTermImpGroups
    type(IntArray) ,optional ,dimension(:) ,intent(in)    :: impTermGenGroups
    type(IntArray) ,optional ,dimension(:) ,intent(in)    :: genTermGroups

    integer(ik) :: i 

    type(IntArray) ,allocatable ,dimension(:) :: iTermIGroups, iTermGGroups, gTermGroups

    character(len=30) :: intToStrBuffer

    if (present(impTermImpGroups)) then 
        iTermIGroups = impTermImpGroups
    else
        allocate(iTermIGroups(size(this%implicitTerms)))
        do i = 1, size(iTermIGroups)
            iTermIGroups(i)%entry = [1]
        end do
    end if

    if (present(impTermGenGroups)) then 
        iTermGGroups = impTermGenGroups
    else
        allocate(iTermGGroups(size(this%implicitTerms)))
        do i = 1, size(iTermGGroups)
            iTermGGroups(i)%entry = [1]
        end do
    end if

    if (present(genTermGroups)) then 
        gTermGroups = genTermGroups
    else
        allocate(gTermGroups(size(this%generalTerms)))
        do i = 1, size(gTermGroups)
            gTermGroups(i)%entry = [1]
        end do
    end if

    if (assertions) then 

        call assert(this%isDefined(),"moveTerms called from undefined term generator")
        call assert(modelObj%isDefined(),"moveTerms called with undefined model")
        call assert(allocated(this%implicitTerms),"moveTerms called before implicit term container of term generator is allocated")
        call assert(allocated(this%generalTerms),"moveTerms called before general term container of term generator is allocated")

        call assert(size(this%implicitTerms) == size(iTermIGroups),&
        "moveTerms called with non-conforming implicit term implict groups")

        call assert(size(this%implicitTerms) == size(iTermGGroups),&
        "moveTerms called with non-conforming implicit term general groups")

        call assert(size(this%generalTerms) == size(gTermGroups),&
        "moveTerms called with non-conforming general term groups")

    end if

    do i = 1, size(this%implicitTerms)
        intToStrBuffer = ""
        write(intToStrBuffer,'(I0)') i
        call modelObj%addImplicitTerm(this%implicitTerms(i)%entry,iTermIGroups(i)%entry,iTermGGroups(i)%entry,&
                                     this%generatorPrefix//"_implicit_"//trim(intToStrBuffer))

    end do

    do i = 1,size(this%generalTerms)
        intToStrBuffer = ""
        write(intToStrBuffer,'(I0)') i
        call modelObj%addGeneralTerm(this%generalTerms(i)%entry,gTermGroups(i)%entry,&
                                      this%generatorPrefix//"_general_"//trim(intToStrBuffer))

    end do
end subroutine moveTerms
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule term_generator_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
