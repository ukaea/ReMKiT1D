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
module term_generator_abstract_class
    !! author: Stefan Mijin
    !!
    !! Houses abstract TermGenerator object, used to automatically generate and add terms to models 

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use god_objects                           ,only: Object
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use support_types                         
    use model_class                           ,only: Model, MatTermContainer ,TermContainer
    use modelbound_data_abstract_class        ,only: ModelboundData

    implicit none
    private

    type ,public :: TermGeneratorContainer
        !! Container object for term generators
        class(TermGenerator) ,allocatable :: entry
    end type TermGeneratorContainer

    type ,public ,extends(object), abstract :: TermGenerator
        !! Abstract TermGenerator object, responsible for automatic generation of terms

        type(MatTermContainer) ,allocatable ,dimension(:) ,private :: implicitTerms 
        type(TermContainer)    ,allocatable ,dimension(:) ,private :: generalTerms

        character(:) ,allocatable ,private :: generatorPrefix

        contains

        procedure ,public :: setImplicitTerms
        procedure ,public :: setGeneralTerms
        procedure ,public :: setGeneratorPrefix

        procedure ,public :: getNumImplicitTerms
        procedure ,public :: getNumGeneralTerms
        procedure ,public :: getGeneratorPrefix

        procedure ,public :: moveTerms

        procedure(generateTerms) ,deferred :: generate

    end type TermGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        subroutine generateTerms(this,mbData) 
            !! Abstract routine for generating terms 

            import :: TermGenerator ,ModelboundData

            class(TermGenerator)  ,intent(inout) :: this 
            class(ModelboundData) ,optional ,intent(in) :: mbData
 
        end subroutine generateTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setImplicitTerms(this,impTerms) 
        !! Move allocation of impTerms to this%implicitTerms

        class(TermGenerator)                              ,intent(inout) :: this
        type(MatTermContainer) ,allocatable ,dimension(:) ,intent(inout) :: impTerms
    
    end subroutine setImplicitTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setGeneralTerms(this,genTerms) 
        !! Move allocation of genTerms to this%generalTerms

        class(TermGenerator)                           ,intent(inout) :: this
        type(TermContainer) ,allocatable ,dimension(:) ,intent(inout) :: genTerms
    
    end subroutine setGeneralTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setGeneratorPrefix(this,prefix) 
        !! Set prefix for added term names

        class(TermGenerator)    ,intent(inout) :: this
        character(*)            ,intent(in)    :: prefix
    
    end subroutine setGeneratorPrefix
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNumImplicitTerms(this) result(numTerms)
        !! Get size of this%implicitTerms

        class(TermGenerator)   ,intent(in)  :: this
        integer(ik)                         :: numTerms
    
    end function getNumImplicitTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getNumGeneralTerms(this) result(numTerms)
        !! Get size of this%generalTerms

        class(TermGenerator)   ,intent(in)  :: this
        integer(ik)                         :: numTerms
    
    end function getNumGeneralTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getGeneratorPrefix(this) result(prefix)
        !! Get size of this%generalTerms

        class(TermGenerator)   ,intent(in)  :: this
        character(:) ,allocatable           :: prefix
    
    end function getGeneratorPrefix
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine moveTerms(this,modelObj,impTermImpGroups,impTermGenGroups,genTermGroups) 
        !! Move terms to modelObj 

        class(TermGenerator)                   ,intent(inout) :: this
        class(Model)                           ,intent(inout) :: modelObj
        type(IntArray) ,optional ,dimension(:) ,intent(in)    :: impTermImpGroups
        type(IntArray) ,optional ,dimension(:) ,intent(in)    :: impTermGenGroups
        type(IntArray) ,optional ,dimension(:) ,intent(in)    :: genTermGroups
    
    end subroutine moveTerms
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module term_generator_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 
