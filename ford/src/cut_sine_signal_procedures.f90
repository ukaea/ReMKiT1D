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
submodule (cut_sine_signal_class) cut_sine_signal_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the cut sine signal class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateCutSine(this,time,period,params) result(res)

    class(CutSineSignal)             ,intent(inout) :: this 
    real(rk)                         ,intent(in)    :: time
    real(rk)                         ,intent(in)    :: period
    real(rk) ,optional ,dimension(:) ,intent(in)    :: params
    real(rk)                                        :: res

    if (assertions) then 
        call assert(present(params),"calculateCutSine reguires parameter params to be present")
        call assert(size(params)>1,"calculateCutSine expects at least 2 params")
        call assert(params(1)>=0,"params(1) passed to calculateCutSine must be non-negative")
        call assert(params(2)>params(1),"params(2) passed to calculateCutSine must be greater than params(1)")
        call assert(params(2)<=real(1,kind=rk),"params(2) passed to calculateCutSine must be less than or equal to 1")
        call assert(period>0,"period argument passed to calculateCutSine must be positive")
    end if

    res = 0

    if (mod(time,period)>=params(1)*period .and. params(2)*period>mod(time,period)) &
    res = sin(pi*(mod(time,period)-params(1)*period)/(period*(params(2)-params(1))))
    
end function calculateCutSine
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule cut_sine_signal_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
