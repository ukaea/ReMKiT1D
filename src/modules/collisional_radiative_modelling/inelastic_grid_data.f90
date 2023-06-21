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
module inelastic_grid_data_class
    !! author: Stefan Mijin 
    !!
    !! Houses  InelasticGridData class 
    !! @todo: Determine whether separate W and emit interpolations are valid and incorporate class into test suite

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use v_space_class               ,only: VSpace
    use sparse_row_data_class       ,only: SparseRowData      
    use inelastic_mapping 
    use lin_interp1D_class          ,only: Interpolation1D
    use support_types

    implicit none
    private

    type ,public ,extends(Object) :: InelasticGridData
        !! Object responsible for storing inelastic mapping matrices, determining active emitter cels, and interpolating variable
        !! energy mapping matrices

        real(rk)            ,allocatable ,dimension(:) ,private :: eInterpGrid !! Energy grid for interpolating variable energy mappings
        type(SparseRowData) ,allocatable ,dimension(:) ,private :: interpW !! Mappings evaluated on the interpolation grid
        real(rk)            ,allocatable ,dimension(:) ,private :: fixedE !! Energies used to evaluate fixed mappings
        type(SparseRowData) ,allocatable ,dimension(:) ,private :: fixedW !! Mappings with fixed transition energies
        type(RealArray)     ,allocatable ,dimension(:) ,private :: fixedEmissionVecs !! Emission vectors for fixed mappings

        integer(ik) :: numV !! Local copy of velocity grid size used for building emission vectors

        contains

        procedure ,public :: getFixedW
        procedure ,public :: getFixedEmissionVector

        procedure ,public :: interpolateW
        procedure ,public :: getInterpolatedEmissionVector

        procedure ,public :: init => initInelData

    end type InelasticGridData
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine initInelData(this,space,fixedEnergies,interpolationGrid) 
        !! Inelastic grid data initialization routine

        class(InelasticGridData)           ,intent(inout)  :: this
        type(VSpace)                       ,intent(in)     :: space !! Velocity space on which weights should be calculated
        real(rk) ,optional ,dimension(:)   ,intent(in)     :: fixedEnergies !! Energy values for fixed energy mappings
        real(rk) ,optional ,dimension(:)   ,intent(in)     :: interpolationGrid !! Energy values for mappings on the interpolation grid

    end subroutine initInelData
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getFixedW(this,ind) result(wMat)
        !! Return fixed mapping matrix with given index

        class(InelasticGridData)       ,intent(in) :: this 
        integer(ik)                    ,intent(in) :: ind
        type(SparseRowData)                        :: wMat
        
    end function getFixedW
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getFixedEmissionVector(this,ind) result(emit)
        !! Return fixed emission vector for mapping matrix with given index

        class(InelasticGridData)       ,intent(in) :: this 
        integer(ik)                    ,intent(in) :: ind
        real(rk) ,allocatable ,dimension(:)        :: emit
        
    end function getFixedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine interpolateW(this,E,wRes) 
        !!  Interpolate mapping matrices for given energy and store in passed triangular matrix. Assumes upper triangular structure for wRes
        !! if E is positive and lower if it's negative

        class(InelasticGridData)     ,intent(in)    :: this
        real(rk)                     ,intent(in)    :: E !! Transition energy to interpolate for
        type(SparseRowData)          ,intent(inout) :: wRes !! Lower/upper triangular matrix to store the interpolated weights

    end subroutine interpolateW
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getInterpolatedEmissionVector(this,E) result(emit)
        !! Return interpolated emission vector for given input energy E

        class(InelasticGridData)       ,intent(in) :: this 
        real(rk)                       ,intent(in) :: E
        real(rk) ,allocatable ,dimension(:)        :: emit
        
    end function getInterpolatedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module inelastic_grid_data_class
!-----------------------------------------------------------------------------------------------------------------------------------
 