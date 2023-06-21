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
submodule (extrapolation_initialization_support) extrapolation_initialization_support_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedure implementations for extrapolation initializastion from JSON data

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
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

    type(NamedString) ,dimension(1) :: extType 

    extType(1) = NamedString(jsonPrefix//"."//keyExtrapolation//"."//keyType,"")

    call jsonCont%load(extType)
    call jsonCont%output(extType)

    select case (extType(1)%value)
    case (keyLinExtrapolation)
        call initLinExtrapolationFromJSON(extObj,jsonPrefix//"."//keyExtrapolation,partObj,gridObj,geometryObj,jsonCont,mpiCont)
    case (keyLogExtrapolation)
        call initLogExtrapolationFromJSON(extObj,jsonPrefix//"."//keyExtrapolation,&
                                          partObj,gridObj,geometryObj,jsonCont,mpiCont,.false.)
    case (keyLinlogExtrapolation)
        call initLogExtrapolationFromJSON(extObj,jsonPrefix//"."//keyExtrapolation,&
                                          partObj,gridObj,geometryObj,jsonCont,mpiCont,.true.)
    case default 
        error stop "Unsupported extrapolation strategy type detected"
    end select
    
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

    type(NamedLogical) ,dimension(1) :: leftBoundary ,staggeredVars

    type(NamedInteger) ,dimension(1) :: usedHaloWidth 
    type(LinExtrapolation) :: linExterpObj

    leftBoundary(1) = NamedLogical(jsonPrefix//"."//keyLeftBoundary,.false.)
    staggeredVars(1) = NamedLogical(jsonPrefix//"."//keyStaggeredVars,.false.)
    usedHaloWidth(1) = NamedInteger(jsonPrefix//"."//keyExpectedHaloWidth,mpiCont%getXHaloWidth())

    call jsonCont%load(leftBoundary)
    call jsonCont%output(leftBoundary)
    call jsonCont%load(staggeredVars)
    call jsonCont%output(staggeredVars)
    call jsonCont%load(usedHaloWidth)
    call jsonCont%output(usedHaloWidth)

    call linExterpObj%init(partObj,gridObj,mpiCont%getWorldRank(),usedHaloWidth(1)%value,&
                           geometryObj,leftBoundary(1)%value,staggeredVars(1)%value)

    allocate(extObj,source=linExterpObj)

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

    type(NamedLogical) ,dimension(1) :: leftBoundary ,staggeredVars

    type(NamedInteger) ,dimension(1) :: usedHaloWidth 
    type(LogExtrapolation) :: logExtrapObj

    leftBoundary(1) = NamedLogical(jsonPrefix//"."//keyLeftBoundary,.false.)
    staggeredVars(1) = NamedLogical(jsonPrefix//"."//keyStaggeredVars,.false.)
    usedHaloWidth(1) = NamedInteger(jsonPrefix//"."//keyExpectedHaloWidth,mpiCont%getXHaloWidth())

    call jsonCont%load(leftBoundary)
    call jsonCont%output(leftBoundary)
    call jsonCont%load(staggeredVars)
    call jsonCont%output(staggeredVars)
    call jsonCont%load(usedHaloWidth)
    call jsonCont%output(usedHaloWidth)

    call logExtrapObj%init(partObj,gridObj,mpiCont%getWorldRank(),usedHaloWidth(1)%value,&
                           geometryObj,leftBoundary(1)%value,staggeredVars(1)%value,linInterp)

    allocate(extObj,source=logExtrapObj)
    
end subroutine initLogExtrapolationFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule extrapolation_initialization_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
