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
module modelbound_data_abstract_class
    !! author: Stefan Mijin
    !! 
    !! Houses abstract model-bound data object used for updating terms

    use data_kinds                            ,only: rk ,ik
    use god_objects                           ,only: Object
    use model_surrogate_class                 ,only: ModelSurrogate
    use variable_container_class              ,only: VariableContainer

    implicit none
    private

    type ,public ,extends(Object), abstract :: ModelboundData
        !! Abstract general modelbound data object used to build classes that contain data for updating terms

        contains

        procedure(dataCalculation)          ,deferred :: update
        procedure(retrieveData)             ,deferred :: copyData
        procedure(getDataDimensionality)    ,deferred :: getDataDim

    end type ModelboundData
!-----------------------------------------------------------------------------------------------------------------------------------
    abstract interface
!-----------------------------------------------------------------------------------------------------------------------------------
        subroutine dataCalculation(this,hostModel,inputVars,updatePriority) 

            import :: ModelboundData ,ModelSurrogate ,VariableContainer ,ik

            class(ModelboundData)                 ,intent(inout) :: this 
            class(ModelSurrogate)                 ,intent(in)    :: hostModel 
            class(VariableContainer)              ,intent(in)    :: inputVars 
            integer(ik) ,optional                 ,intent(in)    :: updatePriority

        end subroutine dataCalculation
!-----------------------------------------------------------------------------------------------------------------------------------
        subroutine retrieveData(this,name,container) 

            import :: rk, ModelboundData

            class(ModelboundData)                 ,intent(in)    :: this 
            character(*)                          ,intent(in)    :: name
            real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container

        end subroutine retrieveData
!-----------------------------------------------------------------------------------------------------------------------------------
        function getDataDimensionality(this,name) result(dim) 

            import :: ik, ModelboundData

            class(ModelboundData)                 ,intent(in)    :: this 
            character(*)                          ,intent(in)    :: name
            integer(ik)                                          :: dim

        end function getDataDimensionality
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module modelbound_data_abstract_class
!-----------------------------------------------------------------------------------------------------------------------------------
 