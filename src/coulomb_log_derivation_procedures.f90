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
submodule (coulomb_log_derivation_class) coulomb_log_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the Coulomb Log derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCoulombLogDeriv(this,ionZ,locNumX,densNorm,tempNorm,electronLog,&
        ionLog,ionZ2,ionMassRatio,removeLogLeiDiscontinuity)
    !! Initialize Coulomb Log derivation object

    class(CoulombLogDerivation)   ,intent(inout) :: this
    real(rk)                      ,intent(in)    :: ionZ 
    integer(ik)                   ,intent(in)    :: locNumX
    real(rk)                      ,intent(in)    :: densNorm
    real(rk)                      ,intent(in)    :: tempNorm
    logical ,optional             ,intent(in)    :: electronLog
    logical ,optional             ,intent(in)    :: ionLog
    real(rk) ,optional            ,intent(in)    :: ionZ2
    real(rk) ,optional            ,intent(in)    :: ionMassRatio
    logical ,optional             ,intent(in)    :: removeLogLeiDiscontinuity

    this%locNumX = locNumX
    this%ionZ = ionZ 
    this%densNorm = densNorm
    this%tempNorm = tempNorm 
    this%electronLog = .false. 
    this%ionLog = .false. 
    this%removeLogLeiDiscontinuity = .false. 
    this%ionZ2 = ionZ 
    this%ionMassRatio = real(1,kind=rk)

    if (present(electronLog)) this%electronLog = electronLog
    if (present(ionLog)) this%ionLog = ionLog
    if (present(removeLogLeiDiscontinuity)) this%removeLogLeiDiscontinuity = removeLogLeiDiscontinuity

    if (this%ionLog) call assert(.not. this%electronLog,"Ambiguous ionLog/electronLog passed to initCoulombLogDeriv")

    if (present(ionZ2)) this%ionZ2 = ionZ2
    if (present(ionMassRatio)) this%ionMassRatio = ionMassRatio
    
    call this%makeDefined()

end subroutine initCoulombLogDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateCoulombLog(this,inputArray,indices) result(output)

    class(CoulombLogDerivation)        ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    logical :: passedZ

    integer(ik) :: i

    passedZ = size(indices) == 3

    if (assertions) then 

        call assertPure(indices>0,"indices passed to calculateCoulombLog out of bounds - lower")
        call assertPure(indices<=size(inputArray),"indices passed to calculateCoulombLog out of bounds - upper")

        if (this%ionLog) then 
            call assertPure(size(indices) == 4,"indices passed to calculateCoulombLog must be size 4 when calculating ion-ion log")
        else 

            call assertPure(any(size(indices) == [2,3]),"indices passed to calculateCoulombLog must be either size 2 or 3 when &
                                                         &calculating electon logs")
        end if
    end if

    allocate(output,mold=inputArray(indices(1))%entry)
    output = 0

    if (this%electronLog) then 

        do i = 1,this%locNumX
            output(i) = logLee(inputArray(indices(1))%entry(i)*this%tempNorm,&
                        inputArray(indices(2))%entry(i)*this%densNorm)
        end do

    else if (this%ionLog) then 

        do i = 1,this%locNumX
            output(i) = logLii(this%ionZ,this%ionZ2,this%ionMassRatio,inputArray(indices(1))%entry(i)*this%densNorm,&
                        inputArray(indices(2))%entry(i)*this%densNorm,inputArray(indices(3))%entry(i)*this%tempNorm,&
                        inputArray(indices(4))%entry(i)*this%tempNorm)
        end do

    else
        
        if (passedZ) then 
            do i = 1,this%locNumX
                output(i) = logLei(inputArray(indices(1))%entry(i)*this%tempNorm,&
                            inputArray(indices(2))%entry(i)*this%densNorm,&
                            inputArray(indices(3))%entry(i),removeDisc=this%removeLogLeiDiscontinuity)
            end do
            
        else
            do i = 1,this%locNumX
                output(i) = logLei(inputArray(indices(1))%entry(i)*this%tempNorm,&
                            inputArray(indices(2))%entry(i)*this%densNorm,&
                            this%ionZ,removeDisc=this%removeLogLeiDiscontinuity)
            end do
    
        end if

    end if

end function calculateCoulombLog
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule coulomb_log_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
