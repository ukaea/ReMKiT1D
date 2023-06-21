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
module f_scaling_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation that extrapolates a distribution function to a boundary or last cell before boundary based on density scaling

    use data_kinds                    ,only: rk ,ik
    use runtime_constants             ,only: debugging, assertions
    use assertion_utility             ,only: assert, assertIdentical, assertPure
    use god_objects                   ,only: Object
    use support_types                 ,only: RealArray 
    use mat_derivation_abstract_class ,only: MatDerivation
    use partition_class               ,only: Partition

    implicit none
    private

    type ,public ,extends(MatDerivation) :: FScalingDerivation
        !! Extrapolates a distribution function to the boundary or last cell centre as a (numH,numV) matrix. Assumes that there are no m>0 harmonics for staggered grids.
        !! Expects 1-4 variables. The first is the extrapolated distribution, the second the density vector. The third is the staggered density
        !! if the variables are staggered. The value of the density at the boundary (scalar) is either the third or the fourth variable, and should be present
        !! only if extrapolating to the boundary.

        integer(ik) ,dimension(2) :: exterpCoords !! Coordinates used for extrapolation. Should be the indices of the cells closest to the boundary on the regular and dual grids, respectively.

        logical ,private :: hasBoundary !! True if the processor this derivation is calculated on has the corresponding external boundary. If it doesn't, this derivation returns 0.

        logical ,private :: leftBoundary 

        integer(ik) ,private :: numV !! Number of velocity grid points
        integer(ik) ,private :: numH !! Number of harmonics

        logical ,private :: staggeredVars !! True if the distribution has staggered harmonics. Defaults to false.
        logical ,private :: extrapolateToBoundary !! True if the extrapolation should be performed to the cell boundary and not the last cell centre before the boundary.

        contains

        procedure ,public :: init => initFScaling

        procedure ,public :: calculate => calculateFScaling

        procedure ,public :: getScalingFactors

    end type FScalingDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initFScaling(this,partitionObj,procRank,numV,leftBoundary,staggeredVars,extrapolateToBoundary)
        !! Initialize distribution scaling extrapolation derivation object

        class(FScalingDerivation)     ,intent(inout)  :: this
        type(Partition)               ,intent(in)     :: partitionObj !! Partition object used to determine local processor grid chunk
        integer(ik)                   ,intent(in)     :: procRank !! Current processor rank
        integer(ik)                   ,intent(in)     :: numV !! Number of expected velocity grid points
        logical  ,optional            ,intent(in)     :: leftBoundary !! True if extrapolating to left boundary. Defaults to false.
        logical  ,optional            ,intent(in)     :: staggeredVars !! True if the distribution has staggered harmonics. Defaults to false.
        logical  ,optional            ,intent(in)     :: extrapolateToBoundary !! True if the extrapolation should be performed to the cell boundary and not the last cell centre before the boundary.

    end subroutine initFScaling  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateFScaling(this,inputArray,indices) result(output)

        class(FScalingDerivation)            ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:)  ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:,:)               :: output

    end function calculateFScaling
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getScalingFactors(this,inputArray,indices) result(output)

        class(FScalingDerivation)            ,intent(in)    :: this 
        type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:)  ,intent(in)    :: indices           
        real(rk)              ,dimension(2)                 :: output

    end function getScalingFactors
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module f_scaling_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 