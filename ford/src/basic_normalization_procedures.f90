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
submodule (basic_normalization_class) basic_normalization_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the basic normalization class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initNormalizationFromJSON(this,jsonCont,mpiCont)
    !! Initialize normalization based on config.json, loading temperature, density, and reference ion charge. 
    !! The resulting object will have the following normalization quantities:
    !! 
    !! 1. density (in m^{-3})
    !! 2. temperature (in eV)
    !! 3. reference ion Z
    !! 4. velocity (used for the velocity grid)
    !! 5. speed (here equal to the velocity)
    !! 6. time (normalized to e-i collision time)
    !! 7. length (velocity * time)
    !! 8. EField (here elMass * velocity/(elCharge * time))
    !! 9. heatFlux (here elMass * density * velocity^3 /2)
    !! 10. crossSection (here 1/(time*density*velocity))
    !!
    !! All names taken from key_names module

    class(BasicNormalization) ,intent(inout) :: this
    type(JSONController)      ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)       ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedReal) ,dimension(3) :: inputNorms
    type(NamedReal) ,dimension(10) :: finalNorms 
    type(NamedReal) ,dimension(10) :: finalNormsNameCropped

    real(rk) :: Te ,ne, Z ,E0

    integer(ik) :: i

    if (assertions) call assert(mpiCont%isDefined(),"Undefined MPIController passed to initNormalizationFromJSON")

    inputNorms(1) = NamedReal(keyNormalization//"."//keyDensNorm,real(1.0d19,kind=rk))
    inputNorms(2) = NamedReal(keyNormalization//"."//keyTempEVNorm,real(10,kind=rk))
    inputNorms(3) = NamedReal(keyNormalization//"."//keyIonZNorm,real(1,kind=rk))

    call jsonCont%load(inputNorms)
    
    finalNorms(1:3) = inputNorms 

    ne = inputNorms(1)%value 
    Te = inputNorms(2)%value 
    Z = inputNorms(3)%value

    finalNorms(4) = NamedReal(keyNormalization//"."//keyVelGridNorm,elVthermal(Te))
    finalNorms(5) = NamedReal(keyNormalization//"."//keySpeedNorm,finalNorms(4)%value)
    finalNorms(6) = NamedReal(keyNormalization//"."//keyTimeNorm,collTimeei(Te,ne,Z,logLei(Te,ne,Z)))
    finalNorms(7) = NamedReal(keyNormalization//"."//keyLengthNorm,finalNorms(4)%value*finalNorms(6)%value)

    E0 = elMass * finalNorms(4)%value /(elCharge * finalNorms(6)%value)

    finalNorms(8) = NamedReal(keyNormalization//"."//keyEFieldNorm,E0)

    finalNorms(9) = NamedReal(keyNormalization//"."//keyHeatFluxNorm,elMass * ne * finalNorms(4)%value**3/2)
    finalNorms(10) = NamedReal(keyNormalization//"."//keyCrossSectionNorm,&
                               real(1,kind=rk)/(ne*finalNorms(7)%value))
    finalNormsNameCropped = finalNorms

    do i = 1,10
        finalNormsNameCropped(i)%name = &
        finalNormsNameCropped(i)%name(len(keyNormalization//".")+1:len(finalNormsNameCropped(i)%name))
    end do
    

    call jsonCont%output(finalNorms)

    call this%setNormalizationVals(finalNormsNameCropped)
    call this%makeDefined()

end subroutine initNormalizationFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule basic_normalization_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
