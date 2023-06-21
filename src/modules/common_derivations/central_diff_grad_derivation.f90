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
module central_diff_grad_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation object that returns values proportional to centrally differenced gradient of given variable

    use data_kinds                  ,only: rk ,ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use support_types               ,only: RealArray
    use derivation_abstract_class   ,only: Derivation
    use geometry_class              ,only: Geometry
    use partition_class             ,only: Partition

    implicit none
    private

    type ,public ,extends(Derivation) :: CentralDiffDerivation
        !! Returns equal to the central difference of the first passed (fluid) variable, optionally multiplied by a constant and 
        !! product of fluid variables raised to corresponding powers

        real(rk) ,allocatable ,dimension(:) ,private :: varPowers !! Powers corresponding to each fluid variable - indices in calculate must conform to 
                                                                  !! size of this + 1
        real(rk)                            ,private :: multConst !! Multiplicative constant - default 1
        
        real(rk) ,allocatable ,dimension(:) ,private :: leftMult
        real(rk) ,allocatable ,dimension(:) ,private :: centralMult
        real(rk) ,allocatable ,dimension(:) ,private :: rightMult


        contains

        procedure ,public :: init => initCentralDiffDeriv

        procedure ,public :: calculate => calculateCentralDiffDeriv

    end type CentralDiffDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCentralDiffDeriv(this,refGeometry,refPartition,procRank,varPowers,multConst)
        !! Initialize central difference derivation object
    
        class(CentralDiffDerivation)        ,intent(inout) :: this
        type(Geometry)                      ,intent(in)    :: refGeometry !! Geometry object used to calculate central difference
        type(Partition)                     ,intent(in)    :: refPartition !! Partition object used to calculate central difference
        integer(ik)                         ,intent(in)    :: procRank !! Current processor rank
        real(rk) ,optional ,dimension(:)    ,intent(in)    :: varPowers !! Optional fluid variable powers
        real(rk) ,optional                  ,intent(in)    :: multConst !! Optional multiplicative constant - default 1

    end subroutine initCentralDiffDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateCentralDiffDeriv(this,inputArray,indices) result(output)

        class(CentralDiffDerivation)          ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateCentralDiffDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module central_diff_grad_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 