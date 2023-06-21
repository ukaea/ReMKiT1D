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
submodule (species_class) species_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the species class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initFromJSON(this,name,jsonCont,mpiCont)
    !! Initialize species from JSON file. Here for extensibility 

    class(Species)          ,intent(inout) :: this
    character(*)            ,intent(in)    :: name 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedStringArray) ,dimension(1) :: vars 
    type(NamedReal)        ,dimension(2) :: realParams 
    type(NamedInteger)     ,dimension(1) :: idParam

    if (assertions) call assert(mpiCont%isDefined(),&
    "Undefined mpi controller passed to initFromJSON constructor for species object")

    vars(1)%name = keySpecies//"."//name//"."//keyAssociatedVars
    allocate(vars(1)%values(0))
    idParam(1) = NamedInteger(keySpecies//"."//name//"."//keyID,-1)
    realParams(1) = NamedReal(keySpecies//"."//name//"."//keyCharge,real(1,kind=rk))
    realParams(2) = NamedReal(keySpecies//"."//name//"."//keyMassA,real(1,kind=rk))

    call jsonCont%load(idParam)
    call jsonCont%load(vars)
    call jsonCont%load(realParams)

    call this%init(idParam(1)%value,name,realParams(1)%value,realParams(2)%value,vars(1)%values)

    realParams(1)%value = this%getCharge()
    realParams(2)%value = this%getMass()/amu

    call jsonCont%output(idParam)
    call jsonCont%output(realParams)
    call jsonCont%output(vars)
    
end subroutine initFromJSON  
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initSpecies(this,id,name,Z,A,vars)
    !! Species initialization. Z is the charge and A the atomic mass (in amus). A/Z are ignored if id = 0, which is reserved for electrons.  
    !! vars is a StringArray containing names of associated variables for easy access 

    class(Species)                  ,intent(inout) :: this
    integer(ik)                     ,intent(in)    :: id 
    character(*)                    ,intent(in)    :: name 
    real(rk)                        ,intent(in)    :: Z 
    real(rk)                        ,intent(in)    :: A 
    type(StringArray) ,dimension(:) ,intent(in)    :: vars 

    this%id = id
    this%name = name  

    !Electron case
    if (id == 0) then 

        this%charge = -real(1,kind=rk)
        this%mass = elMass

    else 

        this%charge = Z 
        this%mass = A * amu

    end if

    this%associatedVars = vars 

    call this%makeDefined()
    
end subroutine initSpecies 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getID (this) result(id)
    !! Getter for id

    class(Species)       ,intent(in) :: this
    integer(ik)                      :: id

    if (assertions) call assertPure(this%isDefined(),"getID called on undefined species object")

    id = this%id

end function getID
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMass (this) result(mass)
    !! Getter for mass

    class(Species)    ,intent(in) :: this
    real(rk)                      :: mass

    if (assertions) call assertPure(this%isDefined(),"getMass called on undefined species object")

    mass = this%mass

end function getMass
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCharge (this) result(charge)
    !! Getter for charge

    class(Species)    ,intent(in) :: this
    real(rk)                      :: charge

    if (assertions) call assertPure(this%isDefined(),"getCharge called on undefined species object")

    charge = this%charge

end function getCharge
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getName (this) result(name)
    !! Getter for name

    class(Species)    ,intent(in) :: this
    character(:) ,allocatable     :: name

    if (assertions) call assertPure(this%isDefined(),"getName called on undefined species object")

    name = this%name

end function getName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getAssociatedVars (this) result(associatedVars)
    !! Getter for associatedVars

    class(Species)                      ,intent(in) :: this
    type(StringArray) ,allocatable ,dimension(:)    :: associatedVars

    if (assertions) call assertPure(this%isDefined(),"getAssociatedVars called on undefined species object")

    associatedVars = this%associatedVars 
    
end function getAssociatedVars
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule species_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
