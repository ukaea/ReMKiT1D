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
submodule (simple_transition_class) simple_transition_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the simple transition class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initSimpleTransition(this,locNumX,inState,outState,energy,fixedRate) 
    !! Initialization routine for SimpleTransition object

    class(SimpleTransition)            ,intent(inout)  :: this
    integer(ik)                        ,intent(in)     :: locNumX !! Local number of spatial cells
    integer(ik)                        ,intent(in)     :: inState !! Pre-transition state
    integer(ik)                        ,intent(in)     :: outState !! Post-transition state
    real(rk)                           ,intent(in)     :: energy !! Transition energy
    real(rk)                           ,intent(in)     :: fixedRate !! Fixed transition rate

    real(rk) ,allocatable ,dimension(:) :: rateVec

    allocate(rateVec(locNumX))
    rateVec = fixedRate 
    this%transitionEnergy = energy

    call this%setStates([inState],[outState])
    call this%setRate(rateVec)
    call this%setRateEnergy(rateVec*energy)
    call this%setRateMomentum(rateVec*0)

    call this%makeDefined()

end subroutine initSimpleTransition
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getEnergy(this) result(energyCost)
    !! Returns array representing energy cost of this transition

    class(SimpleTransition)            ,intent(in) :: this 
    real(rk) ,allocatable ,dimension(:)            :: energyCost

    if (assertions) call assertPure(this%isDefined(),"getEnergy called on undefined simple transition object")

    allocate(energyCost(1))
    energyCost = this%transitionEnergy

end function getEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule simple_transition_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
