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
module timestep_controller_abstract_class
    !! author: Stefan Mijin 
    !! 
    !! Houses abstract interface for optional timestep controller object

    use data_kinds                  ,only: rk
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use variable_container_class    ,only: VariableContainer

    implicit none
    private

    type ,public ,extends(Object) ,abstract:: TimestepController
        !! Abstract timestep controller object

        contains

        procedure(timestepCalc) ,deferred :: evaluateTimestep

    end type TimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface 
        function timestepCalc(this,inputVars,currentTimestep) result(timestep)
            !! Calculate timestep in accordance with timestep controller rules and using current timestep and input variables

            import :: TimestepController ,VariableContainer ,rk

            class(TimestepController)             ,intent(inout) :: this 
            class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate timestep
            real(rk)                              ,intent(in)    :: currentTimestep !! Current timestep to be used if the controller rescales timesteps
            real(rk)                                             :: timestep

        end function timestepCalc
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module timestep_controller_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 