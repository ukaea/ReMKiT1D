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
module manipulator_abstract_class
    !! author: Stefan Mijin
    !!
    !! Houses abstract Manipulator object 

    use god_objects                           ,only: Object
    use modeller_surrogate_class              ,only: ModellerSurrogate
    use variable_container_class              ,only: VariableContainer

    implicit none
    private

    type , public :: ManipulatorContainer 
        class(Manipulator) ,allocatable :: entry 
    end type

    type ,public ,extends(Object), abstract :: Manipulator
        !! Abstract Manipulator object used to manipulate data through Modeller callback

        contains

        procedure(manipulation) ,deferred :: affect

    end type Manipulator
!-----------------------------------------------------------------------------------------------------------------------------------
!Abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        subroutine manipulation(this,manipulatedModeller,outputVars,inputVars) 
            !! Transform inputVars data to outputVars data based on this Manipulator and passed Modeller callback
            import :: Manipulator ,ModellerSurrogate ,VariableContainer

            class(Manipulator)                    ,intent(inout) :: this 
            class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller
            !! Modeller to be used in callback data manipulation
            class(VariableContainer)              ,intent(inout) :: outputVars
            !! VariableContainer object to store the manipulation output 
            class(VariableContainer)              ,intent(in)    :: inputVars
            !! VariableContainer object housing input data for the manipulation routine

        end subroutine manipulation
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module manipulator_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 