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
module interp_stencil_gen_class
    !! author: Stefan Mijin
    !!
    !! Stencil generator for implicit interpolation/extrapolation between regular and staggered grids

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions, assertionLvl
    use god_objects                            ,only: Object
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use support_types                          ,only: RealArray
    use variable_container_class               ,only: VariableContainer
    use jagged_array_generator_class           ,only: JaggedArrayGenerator
    use partition_class                        ,only: Partition
    use geometry_class                         ,only: Geometry
    use modelbound_data_abstract_class         ,only: ModelboundData
    use model_surrogate_class                  ,only: ModelSurrogate

    implicit none
    private
   
    type ,public ,extends(JaggedArrayGenerator) :: InterpStencilGenerator
        !! JaggedArrayGenerator for calculating implicit interpolation stencil between regular and staggered grids. Expects that the
        !! xStencil is [0,1] for regular to staggered grid interpolation, and [-1,0] for staggered to regular. If non-periodic and
        !! going from staggered to regular the values are implicitly extrapolated to the first and last cell centres from the edges 
        !! leading up to the boundary. In this case the generator expects that the left boundary stencil is overriden with [0,1] and the right
        !! stencil with [-2,-1]

        real(rk) ,allocatable ,dimension(:) ,private :: linInterp !! Linear interpolation coefficients for inner points
        real(rk)                            ,private :: linExterpL !! Linear extrapolation coefficient for the left boundary point during inverse interpolation on non-periodic grid
        real(rk)                            ,private :: linExterpR !! Linear extrapolation coefficient for the right boundary point during inverse interpolation on non-periodic grid

        logical ,private :: periodicGrid !! True if stencil is calculated on periodic grid
        logical ,private :: staggeredGridMode !!  If true will interpolate from staggered to regular grid. Defaults to false

        logical ,private :: containsLeftBoundary !! True if the local x-grid contains the left boundary
        logical ,private :: containsRightBoundary !! True if the local x-grid contains the right boundary

        integer(ik) ,private :: locNumX !! Size of local grid chunk (without any halos)

        contains

        procedure ,public :: init => initInterpValGen

        procedure ,public :: calculateInPlace => calcInterpVals

    end type InterpStencilGenerator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initInterpValGen(this,partitionObj,geometryObj,procRank,staggeredGridMode) 
        !! Central differentiation stencil value generator initialization routine
    
        class(InterpStencilGenerator)             ,intent(inout)  :: this
        type(Partition)                           ,intent(in)     :: partitionObj !! Partition object used to determine local number of rows
        type(Geometry)                            ,intent(in)     :: geometryObj !! Geometry object used to get interpolation/extrapolation coefficients
        integer(ik)                               ,intent(in)     :: procRank !! Current processor rank
        logical ,optional                         ,intent(in)     :: staggeredGridMode !!  If true will interpolate from staggered to regular grid.

    end subroutine initInterpValGen
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine calcInterpVals(this,varCont,res,mbData,hostModel)
        !! Calculate interpolation stencil values in place (does not depend on varCont)

        class(InterpStencilGenerator)               ,intent(inout) :: this
        type(VariableContainer)                     ,intent(in)    :: varCont
        type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
        class(ModelboundData) ,optional             ,intent(in)    :: mbData
        class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    end subroutine calcInterpVals
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module interp_stencil_gen_class
!-----------------------------------------------------------------------------------------------------------------------------------
 