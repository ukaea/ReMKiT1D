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
submodule (flat_nd_data_class) flat_nd_data_procedures
!! author: Stefan Mijin 
!! 
!!  Contains module procedures associated with FlatNDData

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine initFlatNDData(this,array) 
            !! Initializes the data based on a deferred rank input array

            class(FlatNDData)                ,intent(inout)  :: this
            real(rk)             ,dimension(..)  ,intent(in) :: array

            this%dims = shape(array)
            select rank (array)
            rank(1) ; this%data = pack(array,.true.)
            rank(2) ; this%data = pack(array,.true.)
            rank(3) ; this%data = pack(array,.true.)
            rank(4) ; this%data = pack(array,.true.)
            rank(5) ; this%data = pack(array,.true.)
            rank(6) ; this%data = pack(array,.true.)
            rank(7) ; this%data = pack(array,.true.)
            rank default 
                error stop "rank of input array to FlatNDData not supported"
            end select

            call this%makeDefined()
        end subroutine initFlatNDData   
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module subroutine directInit(this,array,dims) 
            !! Initializes the data based on an already flat array. Requires dimensions/shape to be explicitly passed

            class(FlatNDData)                ,intent(inout)  :: this
            real(rk)               ,dimension(:) ,intent(in) :: array
            integer(ik)            ,dimension(:) ,intent(in) :: dims

            if (assertions) call assertPure(size(array) == product(dims),&
                            "directInit called on FlatNDData with non-conforming dims and array")

            this%dims = dims 
            this%data = array 

            call this%makeDefined()

        end subroutine directInit 
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getDims (this) result(dims)
            !! Return shape of data stored
 
            class(FlatNDData)                ,intent(in) :: this
            integer(ik) ,allocatable ,dimension(:)       :: dims

            if (assertions) then 
                call assertPure(this%isDefined(),"getValue called on undefined FlatNDData object")
            end if

            dims = this%dims
 
        end function getDims
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function getValue (this,indexSet) result(val)
            !! Return multidimensional value for given indexSet (should match array dimension)
 
            class(FlatNDData)                ,intent(in) :: this
            integer(ik) ,dimension(:)        ,intent(in) :: indexSet
            real(rk)                                     :: val

            integer(ik) :: i ,index

            if (assertions) then 
                call assertPure(this%isDefined(),"getValue called on undefined FlatNDData object")
                call assertPure(size(indexSet) == size(this%dims),"indexSet passed to getValue for FlatNDData has wrong size")

                call assertPure(all(this%dims >= indexSet),&
                                "Index out of range for indexSet passed to getValue of FlatNDData - upper")
                call assertPure(all(0 < indexSet),&
                                "Index out of range for indexSet passed to getValue of FlatNDData - lower")


            end if
            index = indexSet(1)
            do i = 2, size(this%dims)
                index = index + (indexSet(i)-1)*product(this%dims(:i-1))
            end do

            val = this%data(index)
 
     end function getValue
!-----------------------------------------------------------------------------------------------------------------------------------
     pure module function get1DSlice (this,indexSet,sliceIndex) result(val)
        !! Return 1D slice of data at dimension given by sliceIndex and with the other indices set to indexSet

        class(FlatNDData)                ,intent(in) :: this
        integer(ik) ,dimension(:)        ,intent(in) :: indexSet
        integer(ik)                      ,intent(in) :: sliceIndex
        real(rk) ,allocatable ,dimension(:)          :: val

        val = this%data(this%get1DSliceIndices(indexSet,sliceIndex))

    end function get1DSlice
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function get1DSliceIndices (this,indexSet,sliceIndex) result(val)
        !! Return 1D slice indices of data at dimension given by sliceIndex and with the other indices set to indexSet

        class(FlatNDData)                ,intent(in) :: this
        integer(ik) ,dimension(:)        ,intent(in) :: indexSet
        integer(ik)                      ,intent(in) :: sliceIndex
        integer(ik) ,allocatable ,dimension(:)       :: val

        integer(ik) :: offset ,i

        if (assertions) then 
            call assertPure(this%isDefined(),"get1DSliceIndices called on undefined FlatNDData object")
            call assertPure(size(indexSet) == size(this%dims)-1,&
                                "indexSet passed to get1DSliceIndices for FlatNDData has wrong size")

            call assertPure(all(pack(this%dims,[(i .ne. sliceIndex, i=1,size(this%dims))]) >= indexSet),&
                            "Index out of range for indexSet passed to get1DSliceIndices of FlatNDData - upper")
            call assertPure(all(0 < indexSet),&
                            "Index out of range for indexSet passed to get1DSliceIndices of FlatNDData - lower")


        end if

        offset = indexSet(1)
        if (sliceIndex == 1) then 
            offset = 0
            do i = 2, size(this%dims)
                offset = offset + (indexSet(i-1)-1)*product(this%dims(:i-1))
            end do

            val = offset + [(i,i=1,this%dims(sliceIndex))] 

        else 
            do i = 2, sliceIndex-1
                offset = offset + (indexSet(i)-1)*product(this%dims(:i-1))
            end do

            do i = sliceIndex+1, size(this%dims)
                offset = offset + (indexSet(i-1)-1)*product(this%dims(:i-1))
            end do

            val = offset + ([(i,i=1,this%dims(sliceIndex))] - 1) * product(this%dims(:sliceIndex-1))
        end if

    end function get1DSliceIndices
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule flat_nd_data_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
