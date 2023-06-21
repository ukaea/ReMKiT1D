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
module basic_normalization_class
    !! author: Stefan Mijin
    !! 
    !! Houses simple normalization class based on electron quantities and using JSON file data

    use data_kinds                            ,only: rk, ik
    use runtime_constants                     ,only: debugging, assertions
    use assertion_utility                     ,only: assert, assertIdentical, assertPure
    use god_objects                           ,only: Object
    use support_types                         ,only: NamedReal
    use mpi_controller_class                  ,only: MPIController 
    use json_controller_class                 ,only: JSONController
    use normalization_abstract_class          ,only: Normalization
    use physical_constants
    use physics_functions
    use key_names 
    

    implicit none
    private

    type ,public ,extends(Normalization) :: BasicNormalization
        !! Basic normalization based on electron quantities and containing most common values

        contains

        procedure ,public :: init => initNormalizationFromJSON

    end type BasicNormalization
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
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

    end subroutine initNormalizationFromJSON  
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module basic_normalization_class
!-----------------------------------------------------------------------------------------------------------------------------------
 