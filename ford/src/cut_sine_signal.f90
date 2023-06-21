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
module cut_sine_signal_class
    !! author: Stefan Mijin
    !! 
    !! Houses cut sine signal class - periodic signal with passed period, interpeting params(1) and params(2) as the point in the period 
    !! when the signal turns on/off, and is the shape of a half-sine

    use data_kinds                           ,only: rk
    use signal_abstract_class                ,only: TimeSignal
    use runtime_constants                    ,only: debugging, assertions
    use assertion_utility                    ,only: assert, assertIdentical, assertPure
    use physical_constants                   ,only: pi

    implicit none
    private

    type ,public ,extends(TimeSignal) :: CutSineSignal
        !! Periodic signal where in params(1) and params(2) are between 0 and 1. Returns half-sine if params(1)<mod(time,period)<params(2).
        !! Otherwise returns 0.

        contains

        procedure ,public :: calculate => calculateCutSine

    end type CutSineSignal
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
        module function calculateCutSine(this,time,period,params) result(res)

            class(CutSineSignal)             ,intent(inout) :: this 
            real(rk)                         ,intent(in)    :: time
            real(rk)                         ,intent(in)    :: period
            real(rk) ,optional ,dimension(:) ,intent(in)    :: params
            real(rk)                                        :: res


        end function calculateCutSine
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module cut_sine_signal_class
!-----------------------------------------------------------------------------------------------------------------------------------
 