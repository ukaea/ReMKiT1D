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
module support_types
    !! Contains common useful datatypes

    use data_kinds                  ,only: rk, ik

    implicit none
    public

    type ,public :: IntArray
        !! Array of integer values

        integer(ik) ,allocatable ,dimension(:) ,public :: entry 

    end type IntArray

    type ,public :: RealArray
        !! Array of real values 

        real(rk) ,allocatable ,dimension(:) ,public :: entry 

    end type RealArray

    type ,public :: RealArrayD2
        !! Real array of depth 2 

        type(RealArray) ,allocatable ,dimension(:) ,public :: entry 

    end type RealArrayD2

    type ,public :: StringArray
        !! Array of character arrays (a string)

        character(:) ,allocatable ,public :: string
  
    end type StringArray

    type ,public :: LogicalArray
        !! Array of logical values

        logical ,allocatable ,dimension(:) ,public :: entry
  
    end type LogicalArray

    type ,public :: NamedInteger
        !! Named integer value 

        character(:) ,allocatable ,public :: name
        integer(ik)               ,public :: value

    end type NamedInteger

    type ,public :: NamedReal
        !! Named real value

        character(:) ,allocatable ,public :: name
        real(rk)                  ,public :: value

    end type NamedReal

    type ,public :: NamedLogical
        !! Named logical value

        character(:) ,allocatable ,public :: name
        logical                   ,public :: value

    end type NamedLogical

    type ,public :: NamedString
        !! Named string 

        character(:) ,allocatable ,public :: name
        character(:) ,allocatable ,public :: value

    end type NamedString

    type ,public :: NamedScalarContainer 
        !! Container with arrays on named integers, reals, logicals, and strings

        type(NamedInteger) ,allocatable ,dimension(:) ,public :: intData
        type(NamedReal)    ,allocatable ,dimension(:) ,public :: realData
        type(NamedLogical) ,allocatable ,dimension(:) ,public :: logicalData
        type(NamedString)  ,allocatable ,dimension(:) ,public :: stringData

    end type NamedScalarContainer

    type ,public :: NamedIntegerArray
        !! Named integer array 

        character(:) ,allocatable               ,public :: name
        integer(ik)  ,allocatable, dimension(:) ,public :: values

    end type NamedIntegerArray

    type ,public :: NamedRealArray
        !! Named real value array

        character(:) ,allocatable ,public               :: name
        real(rk)     ,allocatable, dimension(:) ,public :: values

    end type NamedRealArray

    type ,public :: NamedLogicalArray
        !! Named logical value array 

        character(:) ,allocatable               ,public :: name
        logical      ,allocatable, dimension(:) ,public :: values

    end type NamedLogicalArray

    type ,public :: NamedStringArray
        !! Named array of strings

        character(:)                            ,allocatable ,public :: name
        type(StringArray) ,allocatable, dimension(:)         ,public :: values

    end type NamedStringArray

    type ,public :: NamedArrayContainer 
        !! Container with named arrays of integers, reals, logicals, and strings

        type(NamedIntegerArray) ,allocatable ,dimension(:) ,public :: intData
        type(NamedRealArray)    ,allocatable ,dimension(:) ,public :: realData
        type(NamedLogicalArray) ,allocatable ,dimension(:) ,public :: logicalData
        type(NamedStringArray)  ,allocatable ,dimension(:) ,public :: stringData
        
    end type NamedArrayContainer
!-----------------------------------------------------------------------------------------------------------------------------------
 end module support_types
!-----------------------------------------------------------------------------------------------------------------------------------
 