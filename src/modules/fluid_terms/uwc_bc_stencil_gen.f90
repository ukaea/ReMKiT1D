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
module uwc_bc_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for flux-like boundary conditions
    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions, assertionLvl
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
   
    type ,public ,extends(JaggedArrayGenerator) :: UWCBCStencilValGenerator
        !! JaggedArrayGenerator for calculating a flux-like boundary condition stencil. Can take 2 variable indices, one for a flux jacobian
        !! variable and the other for a lower jacobian bound variable. If the flux jacobian is present it is extrapolated to the boundary and the stencil 
        !! multiplied by it. If the lower jacobian bound variable is present it is extrapolated to the boundary and substituted for the flux jacobian if 
        !! the jacobian's projection to the boundary surface normal is smaller than the lower bound variable. If setting extrapolate to true the stencils
        !! will extrapolate the implicit/column variables as well. The stencil is assumed to either be [0,-1] for right boundary or [0,1] for left boundary.
        !! Will return size(0) stencil if the current processor does not contain the corresponding boundary.

        real(rk) ,private :: outerJ !! Outer/row inverse jacobian/hodge star
        real(rk) ,private :: linInterp !! linear interpolation coefficient for face before the boundary - used for flux and lower bound variable
        real(rk) ,private :: linExterp !! Linear extrapolation coefficient for column variables
        real(rk) ,private :: innerJ !! Inner/column jacobian/hodge star at boundary face

        logical ,private :: leftBoundary !! True if the stencil values are for the left boundary condition
        logical ,private :: extrapolate !! True if the stencil should extrapolate implicit/column variables in addition to the jacobian

        integer(ik) ,private :: minX !! Local minX (obtained from partition data)
        integer(ik) ,private :: maxX !! Local maxX
        integer(ik) ,private :: inferredGridSize

        integer(ik) ,allocatable ,private :: interpVarIndex !! Variable to optionally interpolate as part of stencil
        integer(ik) ,allocatable ,private :: lowerBoundVarIndex !! Variable index of lower jacobian bound variable

        real(rk) ,allocatable ,private :: fixedLowerBound !! A constant lower bound value (ignored if lowerBoundVarIndex is allocated)
        contains

        procedure ,public :: init => initUWCBCStencil

        procedure ,public :: calculateInPlace => calcUWCBCVals

    end type UWCBCStencilValGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initUWCBCStencil(this,partitionObj,procRank,innerJ,outerJ,linInterp,extrapolate,&
                                       interpVarIndex,leftBoundary,lowerBoundVarIndex,fixedLowerBound,linExterp) 
        !! Flux-like boundary condition stencil value generator initialization routine
    
        class(UWCBCStencilValGenerator)           ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        real(rk)                                  ,intent(in)     :: innerJ !! Inner/column jacobian/hodge star at boundary face
        real(rk)                                  ,intent(in)     :: outerJ !! Outer/row inverse jacobian/hodge star at boundary cell
        real(rk)                                  ,intent(in)     :: linInterp !! Linear interpolation coefficient at face before boundary
        logical ,optional                         ,intent(in)     :: extrapolate !! Extrapolate implicit/column variables. Defaults to .false.
        logical ,optional                         ,intent(in)     :: leftBoundary !! True if the stencil values are for the left boundary condition
        integer(ik) ,optional                     ,intent(in)     :: interpVarIndex !! Variable to optionally interpolate as part of stencil (flux jacobian)
        integer(ik) ,optional                     ,intent(in)     :: lowerBoundVarIndex !! Variable index of lower jacobian bound variable
        real(rk)    ,optional                     ,intent(in)     :: fixedLowerBound !! A constant lower bound value (ignored if lowerBoundVarIndex is allocated)
        real(rk)    ,optional                     ,intent(in)     :: linExterp !! Linear extrapolation coefficient for column variables. Computed from linInterp by default
        
    end subroutine initUWCBCStencil
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcUWCBCVals(this,varCont,res,mbData,hostModel)
        !! Calculate flux-like boundary stencil values in place 

        class(UWCBCStencilValGenerator)             ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcUWCBCVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module uwc_bc_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 