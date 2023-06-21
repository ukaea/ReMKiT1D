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
module operator_abstract_class
    !! author: Stefan Mijin 
    !!
    !! Houses abstract Operator class determining the update and actOn interfaces

    use data_kinds                     ,only: rk
    use god_objects                    ,only: Object
    use variable_container_class       ,only: VariableContainer
    use runtime_constants              ,only: debugging, assertions
    use assertion_utility              ,only: assert, assertIdentical, assertPure
    use basic_interfaces               ,only: realArrayFunction

    implicit none
    private

    type ,public ,extends(Object) ,abstract :: Operator
        !! Abstract operator class for use in explicit terms to transform real arrays

        contains

        procedure                            ,public   :: update => noUpdate
        procedure(realArrayFunction) ,nopass ,deferred :: actOn

    end type Operator
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine noUpdate(this,varCont) 
            !! Default operator update function - does nothing 

            class(Operator)         ,intent(inout)  :: this
            type(VariableContainer) ,intent(in)     :: varCont

        end subroutine noUpdate
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module operator_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 