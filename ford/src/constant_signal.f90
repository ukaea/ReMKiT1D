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
module constant_signal_class
    !! author: Stefan Mijin
    !! 
    !! Houses constant signal class - always returns 1

    use data_kinds                            ,only: rk
    use signal_abstract_class                 ,only: TimeSignal
    

    implicit none
    private

    type ,public ,extends(TimeSignal) :: ConstSignal
        !! Signal wrapper always returning 1

        contains

        procedure ,public :: calculate => calculateConst

    end type ConstSignal
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
        module function calculateConst(this,time,period,params) result(res)

            class(ConstSignal)               ,intent(inout) :: this 
            real(rk)                         ,intent(in)    :: time
            real(rk)                         ,intent(in)    :: period
            real(rk) ,optional ,dimension(:) ,intent(in)    :: params
            real(rk)                                        :: res


        end function calculateConst
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module constant_signal_class
!-----------------------------------------------------------------------------------------------------------------------------------
 