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
module multiplicative_generator_core_class
    !! author: Stefan Mijin
    !!
    !! Core of stencil generators based on tensor products

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray ,IntArray
    use partition_class                        ,only: Partition
    use grid_class                             ,only: Grid
    use support_functions                      ,only: removeDupeInts ,flatTensorProduct ,allCombinations

    implicit none
    private
   
    type ,public ,extends(Object) :: MultiplicativeGeneratorCore
        !! Provides capability of building stencil value objects based on tensor products of individual stencil value array

        integer(ik) ,allocatable ,dimension(:,:) :: localRowCoords !! Local values of row coordinates, used to access passed values

        logical ,private :: fluidCol !! If true will only use x values 
        logical ,private :: vValsDependOnX !! True if the v stencil is expected to have x dependence (assumes to be in a form similar to single harmonic vars)

        real(rk) ,allocatable ,dimension(:) :: fixedVVals !! Fixed velocity stencil values used when the column variable is a distribution and the row variable isn't

        integer(ik) ,private :: expXSize !! Expected size of x stencil
        integer(ik) ,private :: expVSize !! Expected size of v stencil (when it does not depend on x)

        contains

        procedure ,public :: init => initCore

        procedure ,public :: calculateInPlace => calcCoreVals

    end type MultiplicativeGeneratorCore
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCore(this,gridObj,partitionObj,procRank,rowCoords,fluidCol,fixedVVals,vValsDependOnX) 
        !! Multiplicative generator core initialization routine
    
        class(MultiplicativeGeneratorCore)        ,intent(inout)  :: this
        type(Grid)                                ,intent(in)     :: gridObj !! Grid object used to determine local number of rows
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        integer(ik) ,dimension(:,:)               ,intent(in)     :: rowCoords !! Global row coordinate values
        logical ,optional                         ,intent(in)     :: fluidCol !! True if column variable is fluid. Defaults to false.
        real(rk) ,optional ,dimension(:)          ,intent(in)     :: fixedVVals !! Fixed velocity stencil values when the row variable is fluid and the column variable is kinetic
        logical ,optional                         ,intent(in)     :: vValsDependOnX !! True if v stencil has spatial dependence. Defaults to false.

    end subroutine initCore
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcCoreVals(this,res,xVals,hVals,vVals)
        !! Calculate multiplicative core values. All passed stencil values should be indexed starting at the first evolved point. If there are
        !! any gaps in the row values those must also be included for indexing to work.

        class(MultiplicativeGeneratorCore)          ,intent(inout) :: this
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        type(RealArray)           ,dimension(:)     ,intent(in)    :: xVals !! x stencil values 
        type(RealArray) ,optional ,dimension(:)     ,intent(in)    :: hVals !! h stencil values
        type(RealArray) ,optional ,dimension(:)     ,intent(in)    :: vVals !! v stencil values

    end subroutine calcCoreVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module multiplicative_generator_core_class
!-----------------------------------------------------------------------------------------------------------------------------------
 