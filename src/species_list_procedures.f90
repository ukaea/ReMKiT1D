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
submodule (species_list_class) species_list_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the species list class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initSpeciesList(this,jsonCont,mpiCont)
    !! SpeciesList initialization from config.json.

    class(SpeciesList)      ,intent(inout) :: this
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedStringArray) ,dimension(1) :: speciesNames 

    integer(ik) :: i ,j 

    if (assertions .or. assertionLvl >= 0) call assert(mpiCont%isDefined(),&
    "Undefined mpi controller passed to initSpeciesList")

    speciesNames(1)%name = keySpecies//"."//keyNames
    allocate(speciesNames(1)%values(0))

    call jsonCont%load(speciesNames)
    call jsonCont%output(speciesNames)

    allocate(this%speciesEntries(size(speciesNames(1)%values)))
    allocate(this%speciesIDs(size(speciesNames(1)%values)))
    allocate(this%speciesNames(size(speciesNames(1)%values)))

    do i = 1,size(this%speciesEntries)

        this%speciesNames(i) = speciesNames(1)%values(i)
        call this%speciesEntries(i)%initFromJSON(this%speciesNames(i)%string,jsonCont,mpiCont)
        this%speciesIDs(i) = this%speciesEntries(i)%getID()

        if (assertions .or. assertionLvl >= 0) then 

            do j = 1,i-1
                call assert(this%speciesNames(i)%string /= this%speciesNames(j)%string,"Duplicate species name detected")
                call assert(this%speciesIDs(i) /= this%speciesIDs(j),"Duplicate species id detected")
            end do

        end if

    end do

    call this%makeDefined()
    
end subroutine initSpeciesList 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getSpeciesFromID (this,id) result(spec)
    !! Getter for species based on id
    
    class(SpeciesList)   ,intent(in) :: this
    integer(ik)          ,intent(in) :: id 
    type(Species) ,allocatable       :: spec

    integer(ik) :: i ,ind

    logical :: found

    if (assertions) call assertPure(this%isDefined(),"getSpeciesFromID called from undefined species list")

    found = .false. 

    do i = 1, size(this%speciesIDs)

        if (this%speciesIDs(i) == id) then 
            ind = i 
            found = .true.
            exit 
        end if

    end do  

    if (assertions) call assertPure(found,"getSpeciesFromID failed to find species with given id")

    spec = this%speciesEntries(ind)

end function getSpeciesFromID
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getSpeciesIDs (this) result(ids)
    !! Getter for speciesIDs

    class(SpeciesList)                      ,intent(in) :: this
    integer(ik)  ,allocatable, dimension(:)             :: ids

    if (assertions) call assertPure(this%isDefined(),"getSpeciesIDs called from undefined species list")

    ids = this%speciesIDs

end function getSpeciesIDs
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getSpeciesFromName (this,name) result(spec)
    !! Getter for species based on name
    
    class(SpeciesList)   ,intent(in) :: this
    character(*)         ,intent(in) :: name 
    type(Species) ,allocatable       :: spec

    integer(ik) :: i ,ind

    logical :: found

    if (assertions) call assertPure(this%isDefined(),"getSpeciesFromName called from undefined species list")

    found = .false. 

    do i = 1, size(this%speciesNames)

        if (this%speciesNames(i)%string == name) then 
            ind = i 
            found = .true.
            exit 
        end if

    end do  

    if (assertions) call assertPure(found,"getSpeciesFromName failed to find species with given name")

    spec = this%speciesEntries(ind)

end function getSpeciesFromName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getSpeciesVarFromName (this,name,ind) result(var)
    !! Return associated var with index ind of species with given name

    class(SpeciesList)   ,intent(in) :: this
    character(*)         ,intent(in) :: name
    integer(ik)          ,intent(in) :: ind 
    character(:) ,allocatable        :: var

    type(Species) ,allocatable       :: spec
    type(StringArray) ,allocatable ,dimension(:) :: associatedVars
    
    if (assertions) call assertPure(this%isDefined(),"getSpeciesVarFromName called from undefined species list")

    spec = this%getSpeciesFromName(name)

    associatedVars = spec%getAssociatedVars()

    if (assertions) call assertPure(size(associatedVars) >= ind .and. ind > 0,"ind and name combination passed to &
                                    &getSpeciesVarFromName does not result in associated variable")

    var = associatedVars(ind)%string

end function getSpeciesVarFromName
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getSpeciesVarFromID (this,id,ind) result(var)
    !! Return associated var with index ind of species with given id

    class(SpeciesList)   ,intent(in) :: this
    integer(ik)          ,intent(in) :: id 
    integer(ik)          ,intent(in) :: ind 
    character(:) ,allocatable        :: var

    type(Species) ,allocatable       :: spec
    type(StringArray) ,allocatable ,dimension(:) :: associatedVars
    
    if (assertions) call assertPure(this%isDefined(),"getSpeciesVarFromID called from undefined species list")

    spec = this%getSpeciesFromID(id)

    associatedVars = spec%getAssociatedVars()

    if (assertions) call assertPure(size(associatedVars) >= ind .and. ind > 0,"ind and id combination passed to &
                                    &getSpeciesVarFromID does not result in associated variable")

    var = associatedVars(ind)%string

end function getSpeciesVarFromID
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule species_list_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
