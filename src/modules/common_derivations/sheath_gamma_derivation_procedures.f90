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
submodule (sheath_gamma_derivation_class) sheath_gamma_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the sheath gamma derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initElectronSheathGamma(this,massRatio,boundaryIndex)
    !! Initialize electron sheath gamma derivation object

    class(ElectronSheathGammaDerivation)   ,intent(inout) :: this
    real(rk)                               ,intent(in)    :: massRatio
    integer(ik)                            ,intent(in)    :: boundaryIndex

    this%massRatio = massRatio
    this%boundaryIndex = boundaryIndex
    
    call this%makeDefined()

end subroutine initElectronSheathGamma  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateElectronGamma(this,inputArray,indices) result(output)

    class(ElectronSheathGammaDerivation),intent(inout)    :: this 
    type(RealArray)       ,dimension(:)  ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:)  ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                 :: output

    if (assertions) then 

        call assertPure(indices>0,"indices passed to calculateElectronGamma out of bounds - lower")
        call assertPure(indices<=size(inputArray),"indices passed to calculateElectronGamma out of bounds - upper")

        call assertPure(size(indices) == 2,"indices passed to calculateElectronGamma must be size 2")
       
    end if

    allocate(output(1))
    output = 0

    output(1) = real(2,kind=rk) - shPotDrop(this%massRatio,inputArray(indices(2))%entry(this%boundaryIndex)/&
                                         inputArray(indices(1))%entry(this%boundaryIndex))

end function calculateElectronGamma
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule sheath_gamma_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
