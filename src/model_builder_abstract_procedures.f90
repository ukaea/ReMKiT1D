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
submodule (model_builder_abstract_class) model_builder_abstract_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains procedures associated with the abstract ModelBuilder class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadParams(this,jsonCont,mpiCont) 
    !! Load parameters from "./config.json"

    class(ModelBuilder)   ,intent(inout)  :: this
    type(JSONController)  ,intent(inout)  :: jsonCont
    !! JSONController object responsible for reading the config file
    type(MPIController)   ,intent(inout)  :: mpiCont
    !! MPIController object to be used with JSON IO

    if (assertions) then 
        call assert(mpiCont%isDefined(),"Undefined mpiCont passed to loadParams on model builder")
    end if

    call jsonCont%load(this%scalarParams)
    call jsonCont%load(this%arrayParams)

end subroutine loadParams
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputUsedParams(this,jsonCont,mpiCont) 
    !! Output used parameters to "./used_config.json"

    class(ModelBuilder)   ,intent(inout)  :: this
    type(JSONController)  ,intent(inout)  :: jsonCont
    !! JSONController object responsible for writing the config file
    type(MPIController)   ,intent(inout)  :: mpiCont
    !! MPIController object to be used with JSON IO

    if (assertions) then 
        call assert(mpiCont%isDefined(),"Undefined mpiCont passed to outputUsedParams on model builder")
    end if
    call jsonCont%output(this%scalarParams)
    call jsonCont%output(this%arrayParams)

end subroutine outputUsedParams
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getScalarParams(this) result(params)
    !! Getter for ModelBuilder scalarParams

    class(ModelBuilder)   ,intent(in)  :: this
    type(NamedScalarContainer)         :: params

    params = this%scalarParams

end function getScalarParams
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getArrayParams(this) result(params)
    !! Getter for ModelBuilder arrayParams

    class(ModelBuilder)  ,intent(in)  :: this
    type(NamedArrayContainer)         :: params

    params = this%arrayParams

end function getArrayParams
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setScalarParams(this,params) 
    !! Setter for ModelBuilder scalarParams

    class(ModelBuilder)        ,intent(inout)  :: this
    type(NamedScalarContainer) ,intent(in)     :: params

    this%scalarParams = params

end subroutine setScalarParams
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setArrayParams(this,params) 
    !! Setter for ModelBuilder scalarParams

    class(ModelBuilder)        ,intent(inout)  :: this
    type(NamedArrayContainer)  ,intent(in)     :: params

    this%arrayParams = params

end subroutine setArrayParams
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule model_builder_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
