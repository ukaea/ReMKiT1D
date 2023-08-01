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
module flat_nd_data
    !! author: Stefan Mijin 
    !! 
    !! Houses flattened data represenation with multidimensional indexing

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure

    implicit none
    private

    type ,public ,extends(Object) :: FlatNDData
        !! Flat representation of multidimensional data allowing for multidimensional indexing using vector notation - i.e. dataObj%getValue([1,2,3]) = someData(1,2,3). Allows for dimensionality agnostic implementations of various algorithms.

        integer(ik)     ,allocatable ,dimension(:) ,private :: dims !! Dimension sizes of each dimension
        real(rk)        ,allocatable ,dimension(:) ,public  :: data !! Flattened data array

        contains

        procedure ,public :: getValue 
        procedure ,public :: get1DSlice 
        procedure ,public :: get1DSliceIndices

        procedure ,public :: init => initFlatNDData
        procedure ,public :: directInit

    end type FlatNDData
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initFlatNDData(this,array) 
            !! Initializes the data based on a deferred rank input array

            class(FlatNDData)                ,intent(inout)  :: this
            real(rk)              ,dimension(..) ,intent(in) :: array
        end subroutine initFlatNDData   
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine directInit(this,array,dims) 
            !! Initializes the data based on an already flat array. Requires dimensions/shape to be explicitly passed

            class(FlatNDData)                ,intent(inout)  :: this
            real(rk)               ,dimension(:) ,intent(in) :: array
            integer(ik)            ,dimension(:) ,intent(in) :: dims
        end subroutine directInit   
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getValue (this,indexSet) result(val)
            !! Return multidimensional value for given indexSet (should match array dimension)
 
            class(FlatNDData)                ,intent(in) :: this
            integer(ik) ,dimension(:)        ,intent(in) :: indexSet
            real(rk)                                     :: val
 
        end function getValue
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function get1DSlice (this,indexSet,sliceIndex) result(val)
            !! Return 1D slice of data at dimension given by sliceIndex and with the other indices set to indexSet

            class(FlatNDData)                ,intent(in) :: this
            integer(ik) ,dimension(:)        ,intent(in) :: indexSet
            integer(ik)                      ,intent(in) :: sliceIndex
            real(rk) ,allocatable ,dimension(:)          :: val

        end function get1DSlice
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function get1DSliceIndices (this,indexSet,sliceIndex) result(val)
            !! Return 1D slice indices of data at dimension given by sliceIndex and with the other indices set to indexSet

            class(FlatNDData)                ,intent(in) :: this
            integer(ik) ,dimension(:)        ,intent(in) :: indexSet
            integer(ik)                      ,intent(in) :: sliceIndex
            integer(ik) ,allocatable ,dimension(:)       :: val

        end function get1DSliceIndices
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module flat_nd_data
!-----------------------------------------------------------------------------------------------------------------------------------
 