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
module uwc_diff_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for upwinding and central differencing

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
   
    type ,public ,extends(JaggedArrayGenerator) :: UWCDiffStencilValGenerator
        !! JaggedArrayGenerator for calculating an upwinding or central difference stencil based on interpolation to cell faces. Optionally
        !! interpolates a variable with given index in variable container. Interpolation variable required in upwinding mode.

        real(rk) ,allocatable ,dimension(:) ,private :: outerJ !! Outer/row inverse jacobian/hodge star
        real(rk) ,allocatable ,dimension(:) ,private :: linInterp !! linear interpolation coefficients to right cell faces
        real(rk) ,allocatable ,dimension(:) ,private :: linInterpEffective !! linear interpolation coefficients to right cell faces, taking into account potential upwinding
        real(rk) ,allocatable ,dimension(:) ,private :: innerJ !! Inner/column jacobian/hodge star (right cell faces)

        logical ,private :: periodicGrid !! True if stencil is calculated on periodic grid

        integer(ik) ,private :: upwindingMode !! 0 - no upwinding, 1 - upwinding with stencil containing the flux jacobian, 2- upwinding with stencil not including flux jacobian
        integer(ik) ,private :: minX !! Local minX (obtained from partition data)
        integer(ik) ,private :: maxX !! Local maxX

        integer(ik) ,allocatable ,private :: interpVarIndex !! Variable to optionally interpolate as part of stencil
        real(rk) ,allocatable ,dimension(:) ,private :: interpVarBuffer
        logical ,private :: staggeredGridMode !! Remove last row if true

        contains

        procedure ,public :: init => initUWCDiffStencil

        procedure ,public :: calculateInPlace => calcUWCDiffVals

    end type UWCDiffStencilValGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initUWCDiffStencil(this,partitionObj,procRank,innerJ,outerJ,linInterp,xPeriodic,interpVarIndex,upwindingMode,&
                                        staggeredGridMode) 
        !! Central/upwind differentiation stencil value generator initialization routine
    
        class(UWCDiffStencilValGenerator)         ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        real(rk) ,dimension(:)                    ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star (should conform to x-grid)
        real(rk) ,dimension(:)                    ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star (should conform to x-grid)
        real(rk) ,dimension(:)                    ,intent(in)     :: linInterp !! Linear interpolation coefficients to right cell faces (should be size(xGrid)+1)
        logical ,optional                         ,intent(in)     :: xPeriodic !! Used to determine if outer processors should have their stencils trimmed. Defaults to .false.  
        integer(ik) ,optional                     ,intent(in)     :: interpVarIndex !! Optional interpolated variable index 
        integer(ik) ,optional                     ,intent(in)     :: upwindingMode !! 0 - no upwinding, 1 - upwinding with stencil containing the flux jacobian, 2- upwinding with stencil not including flux jacobian. Defaults to 0.
        logical ,optional                         ,intent(in)     :: staggeredGridMode !! Removes last row if true
        
    end subroutine initUWCDiffStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcUWCDiffVals(this,varCont,res,mbData,hostModel)
        !! Calculate central/upwind diff stencil values in place 

        class(UWCDiffStencilValGenerator)           ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcUWCDiffVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module uwc_diff_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 