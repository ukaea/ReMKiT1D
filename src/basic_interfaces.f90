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
module basic_interfaces
    !! author: Stefan Mijin
    !!
    !! Contains a set of basic (abstract) interfaces that do not depend on packages other than basic support used by other modules 

    use data_kinds                  ,only: rk, ik 
    use support_types               ,only: RealArray ,IntArray

    implicit none
    public
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function realFunction(input) result(output)

            import :: rk

            real(rk) ,intent(in) :: input 
            real(rk)             :: output

        end function realFunction
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function realArrayFunction(input) result(output)

            import :: rk

            real(rk)              ,dimension(:) ,intent(in) :: input 
            real(rk) ,allocatable ,dimension(:)             :: output

        end function realArrayFunction
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function realArrayFunctionIntParam(input,param) result(output)

            import :: rk ,ik

            real(rk)              ,dimension(:) ,intent(in) :: input 
            integer(ik)                         ,intent(in) :: param
            real(rk) ,allocatable ,dimension(:)             :: output

        end function realArrayFunctionIntParam
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function realArrayFunctionGenParam(input,realParams,intParams,logicalParams) result(output)

            import :: rk ,ik

            real(rk)               ,dimension(:) ,intent(in) :: input 
            real(rk)     ,optional ,dimension(:) ,intent(in) :: realParams
            integer(ik)  ,optional ,dimension(:) ,intent(in) :: intParams
            logical      ,optional ,dimension(:) ,intent(in) :: logicalParams
            real(rk) ,allocatable ,dimension(:)             :: output

        end function realArrayFunctionGenParam
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function indexedTransform(inputArray,indices) result(output)

            import :: ik ,rk ,RealArray

            type(RealArray)       ,dimension(:) ,intent(in) :: inputArray 
            integer(ik)           ,dimension(:) ,intent(in) :: indices
            real(rk) ,allocatable ,dimension(:)             :: output

        end function indexedTransform
!-----------------------------------------------------------------------------------------------------------------------------------
        pure function coordMapping(inputArray) result(output)

            import :: ik ,IntArray

            integer(ik)    ,dimension(:)        ,intent(in) :: inputArray 
            type(IntArray) ,allocatable ,dimension(:)       :: output

        end function coordMapping
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module basic_interfaces
!-----------------------------------------------------------------------------------------------------------------------------------
 