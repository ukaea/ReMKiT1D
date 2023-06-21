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
module support_functions
    !! Contains various support and math functions

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging ,assertions
    use assertion_utility       
    use support_types               

    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface findIndices
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function findIndicesMatrix(mask) result(indices)
        !! Find locations of .true. values in logical mask of rank 2

        logical     ,dimension(:,:) ,intent(in)   :: mask
        integer(ik) ,dimension(:,:) ,allocatable  :: indices

    end function findIndicesMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function findIndicesVector(mask) result(indices)
        !! Find locations of .true. values in logical mask of rank 1

        logical     ,dimension(:) ,intent(in)   :: mask
        integer(ik) ,dimension(:) ,allocatable  :: indices

    end function findIndicesVector
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    interface removeName
!-----------------------------------------------------------------------------------------------------------------------------------
    pure elemental module function removeNameIntArray(input) result(output)
        !! Remove name from named integer array and return IntArray

        type(NamedIntegerArray) ,intent(in) :: input 
        type(IntArray)                      :: output

    end function removeNameIntArray
!-----------------------------------------------------------------------------------------------------------------------------------
    pure elemental module function removeNameRealArray(input) result(output)
        !! Remove name from named real array and return RealArray

        type(NamedRealArray) ,intent(in) :: input 
        type(RealArray)                  :: output

    end function removeNameRealArray
!-----------------------------------------------------------------------------------------------------------------------------------
    pure elemental module function removeNameLogicalArray(input) result(output)
        !! Remove name from named logical array and return LogicalArray

        type(NamedLogicalArray) ,intent(in) :: input 
        type(LogicalArray)                  :: output

    end function removeNameLogicalArray
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function allCombinations(array) result(res)
        !! Convert IntArray(:) vector into a 2D array of shape (size(IntArray),:) containing all possible combination of IntArray(:)%entries.
        !! Works only for sizes 1 and 3

        type(IntArray) ,dimension(:)   ,intent(in)   :: array
        integer(ik)    ,dimension(:,:) ,allocatable  :: res

    end function allCombinations
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function findNearestPointsInArray(array,point) result(pair)
        !! Finds indices of two of the values nearest to a point in a monotonic array

        real(rk)       ,dimension(:)   ,intent(in)   :: array
        real(rk)                       ,intent(in)   :: point
        integer(ik)    ,dimension(2)                 :: pair

    end function findNearestPointsInArray
!-----------------------------------------------------------------------------------------------------------------------------------
    elemental module function expInt1(val) result(res)
        !! Allan and Hastings approximation of E1 exponential integral, with coefficients from routine in:
        !! Shanjie Zhang, Jianming Jin,
        !! Computation of Special Functions,
        !! Wiley, 1996,
        !! ISBN: 0-471-11963-6,
        !! LC: QA351.C45.

        real(rk) ,intent(in)   :: val
        real(rk)               :: res 

    end function expInt1
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function removeDupeInts(array) result(uniques)
        !! Removes duplicates from integer array - unoptimized

        integer(ik)    ,dimension(:)   ,intent(in)   :: array
        integer(ik) ,allocatable   ,dimension(:)     :: uniques

    end function removeDupeInts
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function withinBounds(array,lowerBound,upperBound) result(res)
        !! Returns map array >= lowerBound .and. array <= upperBound

        integer(ik)    ,dimension(:)   ,intent(in)   :: array
        integer(ik)                    ,intent(in)   :: lowerBound
        integer(ik)                    ,intent(in)   :: upperBound
        logical     ,allocatable   ,dimension(:)     :: res

    end function withinBounds
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function flatTensorProduct(array1,array2,array3) result(res)
        !! Takes a 3-fold tensor product of rank-1 real arrays and flattens it

        real(rk) ,dimension(:)   ,intent(in)  :: array1
        real(rk) ,dimension(:)   ,intent(in)  :: array2
        real(rk) ,dimension(:)   ,intent(in)  :: array3
        real(rk) ,dimension(:)   ,allocatable :: res

    end function flatTensorProduct
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function jaggedArray(array,mask) result(res)
        !! Takes rank-2 array and optional mask and produces rank-1 array of realArrays - a jagged array - by masking the first dimension

        real(rk)           ,dimension(:,:)   ,intent(in)  :: array
        logical  ,optional ,dimension(:,:)   ,intent(in)  :: mask
        type(RealArray) ,dimension(:)   ,allocatable :: res

    end function jaggedArray
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function allPl(points,maxL) result(res)
        !! Return rank-2 array of Legendre polynomials evaluated at given points. Result shape is (size(points),0:maxL). Uses recursion
        !! formula.

        real(rk)           ,dimension(:)   ,intent(in)  :: points
        integer(ik)                        ,intent(in)  :: maxL
        real(rk) ,allocatable ,dimension(:,:)           :: res

    end function allPl
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function shiftedFlatTensorProduct(array1,array2,shiftArray,power) result(res)
        !! Takes 2 arrays and produces a flattened tensor product, such that array1's dimension is the inner dimension. If present,
        !! shift array should conform to array2 and will be added to array1 as the product is taken
        !! The shifted array1 is then optionally raised to a power.

        real(rk)           ,dimension(:)   ,intent(in)  :: array1
        real(rk)           ,dimension(:)   ,intent(in)  :: array2
        real(rk) ,optional ,dimension(:)   ,intent(in)  :: shiftArray
        real(rk) ,optional                 ,intent(in)  :: power
        real(rk)           ,dimension(:)   ,allocatable :: res

    end function shiftedFlatTensorProduct
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function triangularIntArray(dim,lower) result(res)
        !! Returns upper or lower triangular pattern as IntArrays with dimension dim

        integer(ik)                    ,intent(in)  :: dim
        logical  ,optional             ,intent(in)  :: lower
        type(IntArray) ,dimension(:)   ,allocatable :: res

    end function triangularIntArray
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module support_functions
!-----------------------------------------------------------------------------------------------------------------------------------