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
module hdf5_controller_class
    !! Contains serial hdf5 routines for creating files containing variable and grid data as well as parallel restart dump routines

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: Object
    use mpi_controller_class        ,only: MPIController
    use variable_container_class    ,only: VariableContainer
    use grid_class                  ,only: Grid
    use support_types               ,only: StringArray ,RealArray
    use hdf5
    use iso_c_binding
    implicit none 
    public

    type ,extends(Object) ,public :: HDF5Controller
        !! Controller object responsible for HDF5 IO

        type(StringArray) ,allocatable ,dimension(:) ,private :: IOVarNames !! Variable names this controller is responsible for
        integer(ik)       ,allocatable ,dimension(:) ,private :: IOVarIndices !! Indices of variables this controller is responsible for
        character(:)      ,allocatable               ,private :: filepath !! IO filepath
 
        type(RealArray)   ,allocatable ,dimension(:) ,private :: IOBuffer !! Buffer used in HDF5 IO
        
        contains 

        procedure ,public :: outputVarsSerial
        procedure ,public :: loadVarsSerial
        procedure ,public :: outputGridDataSerial
        procedure ,public :: dumpRestartFiles
        procedure ,public :: loadRestartFiles
        procedure ,public :: init => initHDF5Controller

    end type HDF5Controller

!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initHDF5Controller(this,varCont,varNames,filepath) 
        !! HDF5 controller initialization routine - determines variables this controller is responsible for and the output path. The default
        !! path is "./", corresponding to the executable location in the build tree.

        class(HDF5Controller)           ,intent(inout)  :: this
        type(VariableContainer)         ,intent(in)     :: varCont !! Reference variable container
        type(StringArray) ,dimension(:) ,intent(in)     :: varNames !! Names of variables this controller is responsible for
        character(*)      ,optional     ,intent(in)     :: filepath !! Filepath this controller will use 

    end subroutine initHDF5Controller
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine outputVarsSerial(this,mpiCont,varCont,filename,IDNum) 
        !!Gather and output variables this controller is responsible for. Default filename is "ReMKiT1DVarOutput". If IDnum is present it is
        !! appended to the filename. The output is serial from processor 0.

        class(HDF5Controller)           ,intent(inout)  :: this
        type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used for output
        character(*) ,optional          ,intent(in)     :: filename  !! Name of output file
        integer(ik)  ,optional          ,intent(in)     :: IDNum !! IDNum to be appended to filename

    end subroutine outputVarsSerial
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadVarsSerial(this,mpiCont,varCont,filename) 
        !!Load and scatter variables this controller is responsible for. Default filename is "ReMKiT1DVarInput". If a variable is not
        !!found in the file it is set to 0. 

        class(HDF5Controller)           ,intent(inout)  :: this
        type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
        type(VariableContainer)         ,intent(inout)  :: varCont !! Variable container to be loaded into
        character(*) ,optional          ,intent(in)     :: filename  !! Name of input file

    end subroutine loadVarsSerial
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine outputGridDataSerial(this,mpiCont,gridObj,filename) 
        !! Gather and output grid data. Default filename is "ReMKiT1DGridOutput". The output is serial from processor 0.

        class(HDF5Controller)           ,intent(inout)  :: this
        type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
        type(Grid)                      ,intent(in)     :: gridObj !! Output grid object
        character(*) ,optional          ,intent(in)     :: filename  !! Name of output file

    end subroutine outputGridDataSerial
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine dumpRestartFiles(this,mpiCont,varCont,filename) 
        !! Each processor dumps all of its variables into separate files. Default filename "restart", appended with processor rank.

        class(HDF5Controller)           ,intent(inout)  :: this
        type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
        type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used for output
        character(*) ,optional          ,intent(in)     :: filename !! Name of output file

    end subroutine dumpRestartFiles
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadRestartFiles(this,mpiCont,varCont,filename) 
        !! Each processor loads all of its variables from files of the form filename_rank.h5. It is assumed that the files exist and that
        !! they are compatible with the variable container

        class(HDF5Controller)           ,intent(inout)  :: this
        type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
        type(VariableContainer)         ,intent(inout)  :: varCont !! Variable container used for output
        character(*) ,optional          ,intent(in)     :: filename !! Name of output file

    end subroutine loadRestartFiles
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module hdf5_controller_class
!-----------------------------------------------------------------------------------------------------------------------------------