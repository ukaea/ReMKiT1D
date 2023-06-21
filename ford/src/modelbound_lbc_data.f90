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
module modelbound_lbc_data_class
    !! author: Stefan Mijin 
    !!
    !! Houses modelbound data for kinetic logical boundary condition 

    use data_kinds                            ,only: rk ,ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use variable_container_class              ,only: VariableContainer ,CalculationRule
    use modelbound_data_abstract_class        ,only: ModelboundData
    use support_types                         ,only: RealArray ,IntArray
    use variable_list_class                   ,only: VariableList
    use support_functions                     ,only: allPl
    use mat_derivation_abstract_class         ,only: MatDerivation
    use v_space_class                         ,only: VSpace
    use physical_constants

    implicit none
    private

    type ,public ,extends(ModelboundData) :: ModelboundLBCData
        !! Modelbound data object containing data used to construct kinetic logical boundary condition. Assumes no m>0 harmonics.
        !! Assumes standard normalization.
        !! Data accessible through standard data copy routine: "gamma" - sheath heat transmission coefficient
        !!                                                     "potential" - potential drop at sheath normalized to boundary temperature
        !!                                                     "coVel" - cut-off velocity
        !!                                                     "shTemp" - sheath temperature obtained using the cut-off distribution

        real(rk) ,allocatable ,dimension(:,:)   ,private :: fixedPLL !! Fixed Pll tensor component
        real(rk) ,allocatable ,dimension(:,:,:) ,private :: bufferPLL !! Total Pll tensor (for each harmonic combination and velocity cell)

        integer(ik) ,private :: coCell !! Index of cell containing cut-off

        real(rk) ,dimension(2) ,private :: interpdv !! Cell widths of new cells around cut-off
        real(rk) ,dimension(2) ,private :: interpCoord !! Cell centre coordinates of new cells around cut-off
        real(rk) ,dimension(2) ,private :: vInterp !! Interpolation coefficients around cut-off

        real(rk) ,private :: coVel !! Cut-off velocity value
        real(rk) ,private :: gamma !! Sheath heat transmission coeffient 
        real(rk) ,private :: potential !! Sheath potential normalized to boundary temperature
        real(rk) ,private :: shTemp !! Sheath entrance temperature

        real(rk) ,allocatable ,dimension(:,:) ,private :: fext !! Extrapolated distribution function at cell edge
        real(rk) ,allocatable ,dimension(:)   ,private :: vGridCopy !! Local copy of velocity grid
        real(rk) ,allocatable ,dimension(:)   ,private :: dvCopy !! Local copy of velocity cell widths
        real(rk) ,allocatable ,dimension(:)   ,private :: vBoundaries !! Velocity values at left cell boundaries

        logical ,private :: isLeftBoundary !! True if this data is associated with the left boundary

        real(rk) ,private :: bisTol !! False position bisection tolerance 

        class(MatDerivation) ,allocatable      ,private :: fextDeriv !! Derivation of the extrapolated distribution function
        integer(ik) ,allocatable ,dimension(:) ,private :: fextReqIndices !! Indices of variables required for fext derivation

        integer(ik) ,private :: ionCurrentVarIndex !! Ion current variable index
        integer(ik) ,allocatable ,private :: totalCurrentVarIndex !! Total current variable index. Defaults to 0 current.

        integer(ik) ,private :: maxL !! Max l harmonic in decomposition

        logical ,private :: isActive
        contains

        procedure ,public :: init => initMBLBC

        procedure ,public :: update => updateDataLBC
        procedure ,public :: copyData => copyDataLBC

        procedure ,public :: getDataDim => getDataDimLBC

        procedure ,public :: getPll
        procedure ,public :: getCoCell
        procedure ,public :: getInterpCoords
        procedure ,public :: getInterpWidths
        procedure ,public :: getInterpCoeffs

        procedure ,private :: calculatePll
        procedure ,private :: calculateInterps
        procedure ,private :: boundaryHarmonic
        procedure ,private :: interpMom
        procedure ,private :: calculateCutOff

    end type ModelboundLBCData
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine initMBLBC(this,vSpaceObj,distExtDerivation,distExtReqVarIndices,&
                                    ionCurrentVarIndex,isActive,totalCurrentVarIndex,bisTol,isLeftBoundary) 
            !! Varlike modelbound data initialization routine

            class(ModelboundLBCData)        ,intent(inout)  :: this
            type(VSpace)                    ,intent(in)     :: vSpaceObj !! Velocity space object used to get grid data
            class(MatDerivation)            ,intent(in)     :: distExtDerivation !! Extrapolated distribution derivation object
            integer(ik) ,dimension(:)       ,intent(in)     :: distExtReqVarIndices !! Required variable indices for the distribution derivation 
            integer(ik)                     ,intent(in)     :: ionCurrentVarIndex !! Scalar variable index representing ion current at boundary
            logical                         ,intent(in)     :: isActive !! True if the modelbound data should be updated (use to avoid updates on processors with no boundary)
            integer(ik) ,optional           ,intent(in)     :: totalCurrentVarIndex !! Scalar variable index representing total current through boundary. Defaults to 0 current.
            real(rk)    ,optional           ,intent(in)     :: bisTol !! Bisection tolerance for false position method used to calculate cut-off velocity. Defaults to 1e-12
            logical     ,optional           ,intent(in)     :: isLeftBoundary !! True if boundary this data refers to is the left boundary. Defaults to false. 

        end subroutine initMBLBC
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine updateDataLBC(this,hostModel,inputVars,updatePriority) 
            !!  Update modelbound data based on input variable container

            class(ModelboundLBCData)              ,intent(inout) :: this 
            class(ModelSurrogate)                 ,intent(in)    :: hostModel !! Host model - unused
            class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate modelbound data
            integer(ik) ,optional                 ,intent(in)    :: updatePriority !! Priority for this update call (determines which variables are updated) - unused here

        end subroutine updateDataLBC
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine copyDataLBC(this,name,container) 
            !! Copy named modelbound data to passed container. 

            class(ModelboundLBCData)              ,intent(in)    :: this 
            character(*)                          ,intent(in)    :: name !! Name of data
            real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into 

        end subroutine copyDataLBC
!-----------------------------------------------------------------------------------------------------------------------------------
        module function getDataDimLBC(this,name) result(dim)
            !! Get data dimensionality - will return 0 for scalars

            class(ModelboundLBCData)              ,intent(in)    :: this 
            character(*)                          ,intent(in)    :: name !! Name of data
            integer(ik)                                          :: dim

        end function getDataDimLBC
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine calculatePll(this,rowsToUpdate)
            !! Calculate full Pll tensor based on modelbound data 
        
            class(ModelboundLBCData)              ,intent(inout)    :: this 
            integer(ik) ,optional ,dimension(:)   ,intent(in)       :: rowsToUpdate

        end subroutine calculatePll
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine calculateInterps(this)
            !! Calculate interpolation quantities based on modelbound data 

            class(ModelboundLBCData)              ,intent(inout)    :: this 

        end subroutine calculateInterps
!-----------------------------------------------------------------------------------------------------------------------------------
        module function interpMom(this,order,f) result(res)
            !! Calculate moment of single local harmonic with interpolated grid points near cut-off

            class(ModelboundLBCData)              ,intent(inout)    :: this 
            integer(ik)                           ,intent(in)       :: order 
            real(rk) ,dimension(:)                ,intent(in)       :: f

            real(rk) :: res
 
        end function interpMom
!-----------------------------------------------------------------------------------------------------------------------------------
        module function boundaryHarmonic(this,lNum) result(res)
            !! Calculate lNum harmonic at boundary (assuming m=0) using Pll tensor and interpolation

            class(ModelboundLBCData)              ,intent(inout)    :: this 
            integer(ik)                           ,intent(in)       :: lNum 

            real(rk) ,allocatable ,dimension(:) :: res
 
        end function boundaryHarmonic
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine calculateCutOff(this,ionCurrent,totalCurrent)
            !! Calculate cut-off cell and velocity based on modelbound data and currentd
        
            class(ModelboundLBCData)              ,intent(inout)    :: this 
            real(rk)                              ,intent(in)       :: ionCurrent !! ion current into sheath 
            real(rk)                              ,intent(in)       :: totalCurrent !! total current into sheath to be matched

        end subroutine calculateCutOff

!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getCoCell(this) result(res)
            !! Getter for the cut-off cell coordinate

            class(ModelboundLBCData)              ,intent(in)    :: this
            integer(ik)                                          :: res

        end function getCoCell
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getPll(this,lNum) result(res)
            !! Getter for the cut-off (lNum,:,:) decomposition tensor

            class(ModelboundLBCData)             ,intent(in)    :: this
            integer(ik)                          ,intent(in)    :: lNum 
            real(rk) ,allocatable ,dimension(:,:)               :: res

        end function getPll
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getInterpCoords(this) result(res)
            !! Getter for the interpolated cut-off cell coordinates

            class(ModelboundLBCData)  ,intent(in)    :: this
            real(rk) ,dimension(2)                   :: res

        end function getInterpCoords
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getInterpWidths(this) result(res)
            !! Getter for the interpolated cut-off cell widths

            class(ModelboundLBCData)  ,intent(in)    :: this
            real(rk) ,dimension(2)                   :: res

        end function getInterpWidths
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getInterpCoeffs(this) result(res)
            !! Getter for the interpolation coefficients

            class(ModelboundLBCData)  ,intent(in)    :: this
            real(rk) ,dimension(2)                   :: res

        end function getInterpCoeffs
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module modelbound_lbc_data_class
!-----------------------------------------------------------------------------------------------------------------------------------
 