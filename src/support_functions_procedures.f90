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
submodule (support_functions) support_functions_procedures
    !! author: Stefan Mijin  
    !!
    !! Contains the implementations of various support routines

    implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
    contains
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function findIndicesMatrix(mask) result(indices)
        !! Find locations of .true. values in logical mask of rank 2

        logical     ,dimension(:,:) ,intent(in)   :: mask           
        integer(ik) ,dimension(:,:) ,allocatable  :: indices       
        integer(ik)                               :: i, & 
                                                     numIndices     
        logical     ,dimension(:,:) ,allocatable  :: maskTemp

        numIndices = count(mask)                                    !Count number of found indices

        allocate(indices(numIndices,2))
        maskTemp = mask

        do i = 1, numIndices                                        
            indices(i,:) = findloc(maskTemp,.true.)            
            maskTemp(indices(i,1),indices(i,2)) = .false.
        end do

    end function findIndicesMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function findIndicesVector(mask) result(indices)
        !! Find locations of .true. values in logical mask of rank 1

        logical     ,dimension(:) ,intent(in)   :: mask
        integer(ik) ,dimension(:) ,allocatable  :: indices
        integer(ik)                             :: i

        indices = pack([(i,i=1,size(mask))],mask)                   

    end function findIndicesVector
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function allCombinations(array) result(res)
        !! Convert IntArray(:) vector into a 2D array of shape (size(IntArray),:) containing all possible combination of IntArray(:)%entries.
        !! Works only for sizes 1 and 3

        type(IntArray) ,dimension(:)   ,intent(in)   :: array
        integer(ik)    ,dimension(:,:) ,allocatable  :: res
        
        integer(ik) :: i ,j ,k ,l
        integer(ik) ,allocatable ,dimension(:) :: countElements 

        allocate(countElements(size(array)))
        do i = 1,size(array)
            countElements(i) = size(array(i)%entry)
        end do

        allocate(res(size(array),product(countElements)))
        res = 0
        if (size(countElements) == 1) then 
            res(1,:) = array(1)%entry
        end if
        l = 1
        if(size(countElements) == 3) then 
            do i = 1,size(array(1)%entry)
                do j = 1,size(array(2)%entry)
                    do k = 1,size(array(3)%entry)
                        
                        res(1,l) = array(1)%entry(i)
                        res(2,l) = array(2)%entry(j)
                        res(3,l) = array(3)%entry(k)

                        l = l + 1

                    end do 
                end do
            end do

        end if

    end function allCombinations
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function findNearestPointsInArray(array,point) result(pair)
    !! Finds indices of two of the values nearest to a point in a monotonic array

        real(rk)       ,dimension(:)   ,intent(in)   :: array
        real(rk)                       ,intent(in)   :: point
        integer(ik)    ,dimension(2)                 :: pair

        integer(ik) :: closestPoint

        pair = 0 

        if (point < array(1)) then 
            pair(2) = 1
        else if (point > array(size(array))) then 
            pair(1) = size(array)
        else 
            closestPoint = minloc(abs(array - point),1)
            if (array(closestPoint) - point > 0) then 
                pair(1) = closestPoint - 1
                pair(2) = closestPoint 
            else 
                pair(1) = closestPoint 
                pair(2) = closestPoint + 1
            end if
        end if

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

        real(rk) ,parameter :: a(1:6) = real([-0.57721566d0,&
                                               0.99999193d0,&
                                               -0.24991055d0,&
                                               5.519968d-02, &
                                               -9.76004d-03, &
                                               1.07857d-03],kind=rk), &
                               b(1:4) = real([0.2677737343d0, &
                                              8.6347608925d0, &
                                              18.059016973d0,&
                                              8.5733287401d0],kind=rk), &
                               c(1:4) = real([3.9584969228d0, &
                                              21.0996530827d0, &
                                              25.6329561486d0,&
                                              9.5733223454d0],kind=rk)
        real(rk) :: x5(6),x3(4)

        if (val <= 0) error stop "negative value passed to expInt1"

        if (val < real(1,kind=rk)) then 
            x5 = [real(1,kind=rk),val,val**2,val**3,val**4,val**5]
            res = - log(val) + dot_product(a,x5)
        else
            x3 = [real(1,kind=rk),val,val**2,val**3]
            res = exp(val) * dot_product(b,x3)/(val*dot_product(c,x3))
        end if

    end function expInt1
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function removeDupeInts(array) result(uniques)
        !! Returns sorted array of unique integer values

        integer(ik) ,dimension(:) ,intent(in) :: array(:)

        integer(ik) ,dimension(:) , allocatable  :: uniques
        integer(ik) ,dimension(:) , allocatable  :: buffer
        logical ,dimension(:) ,allocatable :: mask

        buffer = array
        uniques = array

        call mergeSort(buffer,1,size(array)+1,uniques) 

        allocate(mask(size(array)))
        mask=.false.
        mask(:size(array)-1) = uniques(:size(array)-1) == uniques(2:)

        uniques=pack(uniques,.not. mask)

    end function removeDupeInts
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function withinBounds(array,lowerBound,upperBound) result(res)
        !! Returns map array >= lowerBound .and. array <= upperBound

        integer(ik)    ,dimension(:)   ,intent(in)   :: array
        integer(ik)                    ,intent(in)   :: lowerBound
        integer(ik)                    ,intent(in)   :: upperBound
        logical     ,allocatable   ,dimension(:)     :: res

        res = (array >= lowerBound) .and. (array <= upperBound)

    end function withinBounds
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function flatTensorProduct(array1,array2,array3) result(res)
        !! Takes a 3-fold tensor product of rank-1 real arrays and flattens it

        real(rk) ,dimension(:)   ,intent(in)  :: array1
        real(rk) ,dimension(:)   ,intent(in)  :: array2
        real(rk) ,dimension(:)   ,intent(in)  :: array3
        real(rk) ,dimension(:)   ,allocatable :: res

        integer(ik) :: i ,j ,k ,l

        allocate(res(size(array1)*size(array2)*size(array3)))
        res = 0
        l = 1
        do i = 1,size(array1)
            do j = 1,size(array2)
                do k = 1,size(array3)
                    
                    res(l) = array1(i)*array2(j)*array3(k)

                    l = l + 1

                end do 
            end do
        end do

    end function flatTensorProduct
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function jaggedArray(array,mask) result(res)
        !! Takes rank-2 array and optional mask and produces rank-1 array of realArrays - a jagged array - by masking the first dimension

        real(rk)           ,dimension(:,:)   ,intent(in)  :: array
        logical  ,optional ,dimension(:,:)   ,intent(in)  :: mask
        type(RealArray) ,dimension(:)   ,allocatable :: res

        integer(ik) :: i

        allocate(res(size(array,2)))

        if (present(mask)) then 
            do i = 1,size(res)
                res(i)%entry = pack(array(:,i),mask(:,i))
            end do
        else
            do i = 1,size(res)
                res(i)%entry = array(:,i)
            end do
        end if

    end function jaggedArray
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function allPl(points,maxL) result(res)
        !! Return rank-2 array of Legendre polynomials evaluated at given points. Result shape is (size(points),0:maxL). Uses recursion
        !! formula.

        real(rk)           ,dimension(:)   ,intent(in)  :: points
        integer(ik)                        ,intent(in)  :: maxL
        real(rk) ,allocatable ,dimension(:,:)           :: res

        integer(ik) :: i

        if (maxL >= 0) then 

            allocate(res(size(points),0:maxL))

            res(:,0) = real(1,kind=rk)

            if (maxL > 0) res(:,1) = points 

            do i = 2,maxL 
                res(:,i) = ((2*i-1)*points*res(:,i-1) - (i-1) * res(:,i-2)) / real(i,kind=rk)
            end do
        else
            allocate(res(size(points),0))
        end if

    end function allPl
!-----------------------------------------------------------------------------------------------------------------------------------
    pure elemental module function removeNameIntArray(input) result(output)
        !! Remove name from named integer array and return IntArray

        type(NamedIntegerArray) ,intent(in) :: input 
        type(IntArray)                      :: output

        output = IntArray(input%values)

    end function removeNameIntArray
!-----------------------------------------------------------------------------------------------------------------------------------
    pure elemental module function removeNameRealArray(input) result(output)
        !! Remove name from named real array and return RealArray

        type(NamedRealArray) ,intent(in) :: input 
        type(RealArray)                  :: output

        output = RealArray(input%values)

    end function removeNameRealArray
!-----------------------------------------------------------------------------------------------------------------------------------
    pure elemental module function removeNameLogicalArray(input) result(output)
        !! Remove name from named logical array and return LogicalArray

        type(NamedLogicalArray) ,intent(in) :: input 
        type(LogicalArray)                  :: output

        output = LogicalArray(input%values)

    end function removeNameLogicalArray
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

        integer(ik) :: i 

        real(rk) ,allocatable ,dimension(:) :: usedShift 
        real(rk)                            :: usedPower

        allocate(res(size(array1)*size(array2)))

        if (assertions) then
            if (present(shiftArray)) call assertPure(size(shiftArray)==size(array2),&
            "shiftArray in shiftedFlatTensorProduct must conform to array2")
        end if

        if (present(shiftArray)) then
            usedShift = shiftArray
        else
            allocate(usedShift(size(array2)))
            usedShift = 0
        end if

        usedPower = real(1,kind=rk)
        if (present(power)) usedPower = power
        res = 0
        do i = 1,size(array2)
                    res((i-1)*size(array1)+1:i*size(array1)) = (array1 + usedShift(i))**power*array2(i)
        end do

    end function shiftedFlatTensorProduct
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function triangularIntArray(dim,lower) result(res)
        !! Returns upper or lower triangular array pattern as IntArrays with dimension dim

        integer(ik)                    ,intent(in)  :: dim
        logical  ,optional             ,intent(in)  :: lower
        type(IntArray) ,dimension(:)   ,allocatable :: res

        logical :: low 
        integer(ik) :: i ,j

        allocate(res(dim))

        low = .false.
        if (present(lower)) low = lower

        if (low) then
            do i = 1,dim
                res(i)%entry = [(j,j=1,i)]
            end do
        else
            do i = 1,dim
                res(i)%entry = [(j,j=i,dim)]
            end do
        end if

    end function triangularIntArray
!-----------------------------------------------------------------------------------------------------------------------------------
    pure recursive subroutine mergeSort(buffer, start, end, array)

        integer(ik) ,dimension(:),intent(inout) :: array, buffer
        integer(ik), intent(in) :: start,end
        integer(ik) :: mid

        if (end-start<2) return   
                
        mid = (end+start)/2
        call mergeSort(array, start, mid, buffer)
        call mergeSort(array, mid, end, buffer)
        call merge(buffer, start, mid, end, array)

    end subroutine mergeSort
!-----------------------------------------------------------------------------------------------------------------------------------
  pure subroutine merge(inputArray, start, mid, end, outputArray)

    integer(ik) ,dimension(:) ,intent(inout) :: inputArray,outputArray
    integer(ik) ,intent(in) :: start,mid,end
    integer(ik)    :: i,leftIndex,rightIndex
    logical :: rightIndexOOB
    leftIndex = start
    rightIndex = mid

    do i = start,end-1
       rightIndexOOB = rightIndex >= end
       if (.not. rightIndexOOB) rightIndexOOB = inputArray(leftIndex) <= inputArray(rightIndex) ! Avoids error in some compiler debug modes
       if (leftIndex < mid .and. rightIndexOOB) then
          outputArray(i)=inputArray(leftIndex)
          leftIndex=leftIndex+1
       else
          outputArray(i)=inputArray(rightIndex)
          rightIndex = rightIndex + 1
       end if
    end do

end subroutine merge
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule support_functions_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
