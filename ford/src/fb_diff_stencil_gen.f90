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
module fb_diff_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for fixed forward or backwards differencing

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray
    use variable_container_class               ,only: VariableContainer
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use partition_class                        ,only: Partition
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: FBDiffStencilValGenerator
        !! JaggedArrayGenerator for calculating forward or backwards difference stencils. These are two-point stencils of the form
        !! [-Jout(i)*Jin(i),Jout(i)*Jin(i+1)], or [-Jout(i)*Jin(i-1),Jout(i)*Jin(i)], depending on whether the difference is forward 
        !! or backwards. Jout and Jin should be defined on the x grid, and should act like the (inverse) Jacobian/Hodge star operator. 

        logical ,private :: backwardsDiff !! If true the difference is backwards. Defaults to .false.

        real(rk) ,allocatable ,dimension(:) ,private :: innerJ !! Inner/column jacobian/hodge star (left jacobian if forward diff and vice versa)
        real(rk) ,allocatable ,dimension(:) ,private :: outerJ !! Outer/row inverse jacobian/hodge star

        logical ,private :: periodicGrid !! True if stencil is calculated on periodic grid
        logical ,private :: staggeredGridMode !! If true will zero out the rightmost contribution to backwards differencing or remove last row for forward differencing - defaults to false

        integer(ik) ,private :: minX !! Local minX (obtained from partition data)
        integer(ik) ,private :: maxX !! Local maxX

        contains

        procedure ,public :: init => initFBValGen

        procedure ,public :: calculateInPlace => calcFBVals

    end type FBDiffStencilValGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initFBValGen(this,partitionObj,procRank,innerJ,outerJ,xPeriodic,backwardsDiff,staggeredGridMode) 
        !! Forward/backward differentiation stencil value generator initialization routine
    
        class(FBDiffStencilValGenerator)          ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        real(rk) ,dimension(:)                    ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star (should conform to x-grid)
        real(rk) ,dimension(:)                    ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star (should conform to x-grid)
        logical ,optional                         ,intent(in)     :: xPeriodic !! Used to determine if outer processors should have their stencils trimmed.  Defaults to .false.  
        logical ,optional                         ,intent(in)     :: backwardsDiff !! If true the difference is backwards. Defaults to .false.
        logical ,optional                         ,intent(in)     :: staggeredGridMode !! If true will zero out the rightmost contribution to backwards differencing or remove last row for forward differencing

    end subroutine initFBValGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcFBVals(this,varCont,res,mbData,hostModel)
        !! Calculate forward/backwards diff stencil values in place (does not depend on varCont)

        class(FBDiffStencilValGenerator)            ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcFBVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module fb_diff_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 