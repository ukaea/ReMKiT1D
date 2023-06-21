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
submodule (interpolation_derivation_class) interpolation_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the interpolation derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initInterpDeriv(this,geometryObj,gridObj,minX,maxX,inverseInterp,distInterp)
    !! Initialize interpolation derivation object

    class(InterpolationDerivation)   ,intent(inout) :: this
    type(Geometry)                   ,intent(in)    :: geometryObj
    type(Grid)                       ,intent(in)    :: gridObj
    integer(ik)                      ,intent(in)    :: minX 
    integer(ik)                      ,intent(in)    :: maxX
    logical ,optional                ,intent(in)    :: inverseInterp
    logical ,optional                ,intent(in)    :: distInterp

    real(rk) ,allocatable ,dimension(:) :: dx ,lInterp,lInterpDual

    integer(ik) :: numX 

    if (assertions) then 
        call assert(geometryObj%isDefined(),"Undefined geometry object passed to initInterpDeriv")
        call assert(gridObj%isDefined(),"Undefined grid object passed to initInterpDeriv")
    end if
    
    this%inverseInterp = .false. 
    if (present(inverseInterp)) this%inverseInterp = inverseInterp

    this%periodicGrid = geometryObj%isPeriodic()
    numX = gridObj%getNumX()
    dx = geometryObj%getCellWidths()
    this%locNumX = maxX - minX + 1
    this%containsLeftBoundary = minX == 1 
    this%containsRightBoundary = maxX == numX

    allocate(lInterp(0:numX))
    lInterp = geometryObj%getLinInterp()
    allocate(lInterpDual(0:numX))
    lInterpDual = geometryObj%getLinInterp(dualGrid=.true.)

    this%linInterp = lInterp(minX:maxX)
    this%linInterpDual = lInterpDual

    this%linExterp = real(1,kind=rk) - lInterp(numX-1)

    this%linExterpLDual = dx(1)/(2*dx(2))
    this%linExterpRDual = dx(numX)/(2*dx(numX-1))

    this%distInterp = .false. 
    if (present(distInterp)) this%distInterp = .true. 

    this%numH = gridObj%getNumH()
    this%numV = gridObj%getNumV()

    allocate(this%oddHarmonic(this%numH))
    this%oddHarmonic = mod(gridObj%getLGrid(),2) /= 0

    call this%makeDefined()

end subroutine initInterpDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateInterp(this,inputArray,indices) result(output)

    class(InterpolationDerivation)     ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    integer(ik) :: i,j

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateInterp called on undefined derivation object")
        call assertPure(size(indices) == 1,"Only one index should be passed to calculateInterp")

    end if
    allocate(output,mold=inputArray(indices(1))%entry)
    output= real(1,kind=rk) !Avoiding zeros in halo 

    if (this%distInterp) then 

        do i = 1, this%numH

            if (this%oddHarmonic(i)) then

                do j = 1,this%locNumX

                    output((j-1)*this%numH*this%numV+(i-1)*this%numV+1:(j-1)*this%numH*this%numV+i*this%numV) = &
                    inputArray(indices(1))%entry((j-2)*this%numH*this%numV+(i-1)*this%numV+1:(j-2)*this%numH*this%numV+i*this%numV)&
                    * (real(1,kind=rk) - this%linInterpDual(j))&
                    + inputArray(indices(1))%entry((j-1)*this%numH*this%numV+(i-1)*this%numV+1&
                    :(j-1)*this%numH*this%numV+i*this%numV) * this%linInterpDual(j)

                end do

                if (.not. this%periodicGrid) then

                    if (this%containsLeftBoundary) &
                    output((i-1)*this%numV+1:i*this%numV) = inputArray(indices(1))%entry((i-1)*this%numV+1:i*this%numV) &
                    * (real(1,kind=rk)+this%linExterpLDual) &
                    - this%linExterpLDual*inputArray(indices(1))%entry(this%numH*this%numV+(i-1)*this%numV+1&
                    :this%numH*this%numV+i*this%numV)

                    if (this%containsRightBoundary) &
                    output((this%locNumX-1)*this%numH*this%numV+(i-1)*this%numV+1&
                    :(this%locNumX-1)*this%numH*this%numV+i*this%numV) = &
                    inputArray(indices(1))%entry((this%locNumX-2)*this%numH*this%numV+(i-1)*this%numV+1&
                    :(this%locNumX-2)*this%numH*this%numV+i*this%numV) * (real(1,kind=rk)+this%linExterpRDual) &
                    - this%linExterpRDual*inputArray(indices(1))%entry((this%locNumX-3)*this%numH*this%numV+(i-1)*this%numV+1&
                    :(this%locNumX-3)*this%numH*this%numV+i*this%numV)

                end if

            else

                do j = 1,this%locNumX
                    output((j-1)*this%numH*this%numV+(i-1)*this%numV+1:(j-1)*this%numH*this%numV+i*this%numV) = &
                    inputArray(indices(1))%entry((j-1)*this%numH*this%numV+(i-1)*this%numV+1:(j-1)*this%numH*this%numV+i*this%numV)&
                    * (real(1,kind=rk) - this%linInterp(j))&
                    + inputArray(indices(1))%entry(j*this%numH*this%numV+(i-1)*this%numV+1:j*this%numH*this%numV+i*this%numV)&
                    * this%linInterp(j)
                end do

                if (.not. this%periodicGrid) then
                    if (this%containsRightBoundary) &
                    output((this%locNumX-1)*this%numH*this%numV+(i-1)*this%numV+1&
                    :(this%locNumX-1)*this%numH*this%numV+i*this%numV) = &
                    inputArray(indices(1))%entry((this%locNumX-1)*this%numH*this%numV+(i-1)*this%numV+1&
                    :(this%locNumX-1)*this%numH*this%numV+i*this%numV) * (real(1,kind=rk)+this%linExterp) &
                    - this%linExterp*inputArray(indices(1))%entry((this%locNumX-2)*this%numH*this%numV+(i-1)*this%numV+1&
                    :(this%locNumX-2)*this%numH*this%numV+i*this%numV)
                end if

            end if

           

        end do

    else

        if (this%inverseInterp) then 
            output(1:this%locNumX) = inputArray(indices(1))%entry(0:this%locNumX-1) * (real(1,kind=rk)&
                                    -this%linInterpDual(1:this%locNumX))&
                                    +  inputArray(indices(1))%entry(1:this%locNumX) * this%linInterpDual(1:this%locNumX)
        else
            output(1:this%locNumX) = inputArray(indices(1))%entry(1:this%locNumX) *(real(1,kind=rk)-this%linInterp(1:this%locNumX))&
                                    +  inputArray(indices(1))%entry(2:this%locNumX+1) * this%linInterp(1:this%locNumX)
        end if
        
        if (.not. this%periodicGrid) then 
            if (this%inverseInterp .and.  this%containsLeftBoundary) &
            output(1) = inputArray(indices(1))%entry(1) * (real(1,kind=rk)+this%linExterpLDual) &
                    - this%linExterpLDual*inputArray(indices(1))%entry(2)

            if (.not. this%inverseInterp .and. this%containsRightBoundary) &
            output(this%locNumX) = inputArray(indices(1))%entry(this%locNumX) * (real(1,kind=rk)+this%linExterp) &
                                - this%linExterp*inputArray(indices(1))%entry(this%locNumX-1)

            if (this%inverseInterp .and. this%containsRightBoundary) &
            output(this%locNumX) = inputArray(indices(1))%entry(this%locNumX-1) * (real(1,kind=rk)+this%linExterpRDual) &
                                - this%linExterpRDual*inputArray(indices(1))%entry(this%locNumX-2)
        end if
    end if

end function calculateInterp
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule interpolation_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
