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
module extrapolation_initialization_support
    !! author: Stefan Mijin
    !!
    !! Contains support for initializing extrapolation objects based on JSON data

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use geometry_class                         ,only: Geometry
    use grid_class                             ,only: Grid 
    use partition_class                        ,only: Partition
    use mpi_controller_class                   ,only: MPIController
    use json_controller_class                  ,only: JSONController
    use extrapolation_abstract_class           ,only: Extrapolation
    use lin_extrapolation_class                ,only: LinExtrapolation
    use log_extrapolation_class                ,only: LogExtrapolation
    use support_types
    use key_names

    implicit none 
    private 

    public :: initExtrapolationFromJSON

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initExtrapolationFromJSON(extObj,jsonPrefix,partObj,gridObj,geometryObj,jsonCont,mpiCont)
    !! Initialize extrapolation from JSON data - requires an olready opened JSON file

        class(Extrapolation) ,allocatable ,intent(inout) :: extObj
        character(*)                      ,intent(in)    :: jsonPrefix 
        type(Partition)                   ,intent(in)    :: partObj
        type(Grid)                        ,intent(in)    :: gridObj
        type(Geometry)                    ,intent(in)    :: geometryObj
        type(JSONController)              ,intent(inout) :: jsonCont
        type(MPIController)               ,intent(inout) :: mpiCont
        
    end subroutine initExtrapolationFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initLinExtrapolationFromJSON(extObj,jsonPrefix,partObj,gridObj,geometryObj,jsonCont,mpiCont)
        !! Initialize linear extrapolation from JSON data

        class(Extrapolation) ,allocatable ,intent(inout) :: extObj
        character(*)                      ,intent(in)    :: jsonPrefix 
        type(Partition)                   ,intent(in)    :: partObj
        type(Grid)                        ,intent(in)    :: gridObj
        type(Geometry)                    ,intent(in)    :: geometryObj
        type(JSONController)              ,intent(inout) :: jsonCont
        type(MPIController)               ,intent(inout) :: mpiCont
        
    end subroutine initLinExtrapolationFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initLogExtrapolationFromJSON(extObj,jsonPrefix,partObj,gridObj,geometryObj,jsonCont,mpiCont,linInterp)
        !! Initialize logarithmic extrapolation from JSON data

        class(Extrapolation) ,allocatable ,intent(inout) :: extObj
        character(*)                      ,intent(in)    :: jsonPrefix 
        type(Partition)                   ,intent(in)    :: partObj
        type(Grid)                        ,intent(in)    :: gridObj
        type(Geometry)                    ,intent(in)    :: geometryObj
        type(JSONController)              ,intent(inout) :: jsonCont
        type(MPIController)               ,intent(inout) :: mpiCont
        logical                           ,intent(in)    :: linInterp
        
    end subroutine initLogExtrapolationFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module extrapolation_initialization_support
!-----------------------------------------------------------------------------------------------------------------------------------