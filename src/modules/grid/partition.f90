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
module partition_class
    !! Houses Partition object responsible for decomposing the x-h domain

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions ,assertionLvl
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use support_functions           ,only: findIndices ,allCombinations ,removeDupeInts
    use support_types               ,only: IntArray

    implicit none
    private

    type ,public ,extends(Object) :: Partition
        !! Object responsible for storing the x-h domain decomposition

        integer(ik) ,allocatable ,dimension(:) ,private :: minX !! First x-index in each partition
        integer(ik) ,allocatable ,dimension(:) ,private :: maxX !! Last x-index in each partition 
        integer(ik) ,allocatable ,dimension(:) ,private :: minH !! First h-index in each partition
        integer(ik) ,allocatable ,dimension(:) ,private :: maxH !! Last h-index in each partition

        integer(ik) ,allocatable ,dimension(:) ,private :: locX !! Number of x-grid points for each partition
        integer(ik) ,allocatable ,dimension(:) ,private :: locH !! Number of h-grid points for each partition

        contains

        procedure ,public :: getMinX
        procedure ,public :: getMinXAtInd
        procedure ,public :: getMinH
        procedure ,public :: getMinHAtInd
        procedure ,public :: getMaxX
        procedure ,public :: getMaxXAtInd
        procedure ,public :: getMaxH
        procedure ,public :: getMaxHAtInd

        procedure ,public :: getLocNumX
        procedure ,public :: getLocNumH

        procedure ,public :: findProc

        procedure ,public :: filterCoords

        procedure ,public :: initSimplePartition
        procedure ,public :: init => initPartition

    end type Partition
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine initSimplePartition(this,numProcsX,numProcsH,numX,numH) 
            !! Partition initialization routine - assuming even distributions in x and h directions

            class(Partition)          ,intent(inout)  :: this
            integer(ik)               ,intent(in) :: numProcsX !! Number of processes in x direction
            integer(ik)               ,intent(in) :: numProcsH !! Number of processes in h direction
            integer(ik)               ,intent(in) :: numX !! Total number of x grid points
            integer(ik)               ,intent(in) :: numH !! Total number of h grid points 

        end subroutine initSimplePartition
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initPartition(this,minX, maxX, minH, maxH) 
            !! Partition initialization routine - general

            class(Partition)          ,intent(inout)  :: this
            integer(ik) ,dimension(:) ,intent(in)     :: minX !! First x-index in each partition
            integer(ik) ,dimension(:) ,intent(in)     :: maxX !! Last x-index in each partition
            integer(ik) ,dimension(:) ,intent(in)     :: minH !! First h-index in each partition
            integer(ik) ,dimension(:) ,intent(in)     :: maxH !! Last x-index in each partition

        end subroutine initPartition
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMinX (this) result(minX)
            !! Getter for minX

            class(Partition)                       ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)             :: minX
 
        end function getMinX
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMinXAtInd (this,ind) result(minX)
            !! Return minX at index ind

            class(Partition)   ,intent(in) :: this
            integer(ik)        ,intent(in) :: ind
            integer(ik)                    :: minX
 
        end function getMinXAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMaxX (this) result(maxX)
            !! Getter for maxX

            class(Partition)                      ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)            :: maxX
 
        end function getMaxX
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMaxXAtInd (this,ind) result(maxX)
            !! Return maxX at index ind

            class(Partition)   ,intent(in) :: this
            integer(ik)        ,intent(in) :: ind
            integer(ik)                    :: maxX
 
        end function getMaxXAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMinH (this) result(minH)
            !! Getter for minH

            class(Partition)                      ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)            :: minH
 
        end function getMinH
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMinHAtInd (this,ind) result(minH)
            !! Return minH at index ind

            class(Partition)   ,intent(in) :: this
            integer(ik)        ,intent(in) :: ind
            integer(ik)                    :: minH
 
        end function getMinHAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMaxH (this) result(maxH)
            !! Getter for maxH

            class(Partition)                       ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)             :: maxH
 
        end function getMaxH
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getMaxHAtInd (this,ind) result(maxH)
            !! Return maxH at index ind

            class(Partition)   ,intent(in) :: this
            integer(ik)        ,intent(in) :: ind
            integer(ik)                    :: maxH
 
        end function getMaxHAtInd
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getLocNumX (this) result(locX)
            !! Getter for locX

            class(Partition)                       ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)             :: locX
 
        end function getLocNumX
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getLocNumH (this) result(locH)
            !! Getter for locH

            class(Partition)                       ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)             :: locH
 
        end function getLocNumH
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function findProc (this,xInd,hInd) result(procInd)
            !! Return partition index which contains xInd and hInd

            class(Partition) ,intent(in) :: this
            integer(ik)      ,intent(in) :: xInd
            integer(ik)      ,intent(in) :: hInd
            integer(ik)                  :: procInd
 
        end function findProc
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function filterCoords (this,ind,coords,normalize) result(res)
            !! Filters coordinate list of shape 1,: or 3,: based on partition data. If normalize is true will shift all values in each
            !! dimension to so that the smallest value is 1

            class(Partition)             ,intent(in) :: this
            integer(ik)                  ,intent(in) :: ind
            integer(ik) ,dimension(:,:)  ,intent(in) :: coords
            logical  ,optional           ,intent(in) :: normalize
            integer(ik) ,dimension(:,:) ,allocatable :: res
 
        end function filterCoords
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module partition_class
!-----------------------------------------------------------------------------------------------------------------------------------
 