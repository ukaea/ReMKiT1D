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
module signal_abstract_class
    !! author: Stefan Mijin
    !! 
    !! Houses abstract signal class - a wrapper for explicit time dependence

    use data_kinds                            ,only: rk

    implicit none
    private

    type ,public :: TimeSignalData
        !! Container for time signal data used tu inject explicit time dependence
        real(rk) ,allocatable ,dimension(:) ,public :: tParams !! Parameters for tSignal
        real(rk)                            ,public :: tPeriod !! Period for tSignal
        class(TimeSignal) ,allocatable      ,public :: tSignal !! Optional signal object used to inject explicit time dependence

    end type TimeSignalData

    type ,public ,abstract :: TimeSignal
        !! Abstract signal class wrapping explicit time dependence. NOTE: Not meant to be initializable, hence not an object

        contains

        procedure(signalCalculation) ,deferred :: calculate

    end type TimeSignal
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        function signalCalculation(this,time,period,params) result(res)

            import :: TimeSignal ,rk

            class(TimeSignal)                    ,intent(inout) :: this 
            real(rk)                         ,intent(in)    :: time
            real(rk)                         ,intent(in)    :: period
            real(rk) ,optional ,dimension(:) ,intent(in)    :: params
            real(rk)                                        :: res


        end function signalCalculation
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module signal_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 