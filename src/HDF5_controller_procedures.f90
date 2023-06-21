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
submodule (HDF5_controller_class) HDF5_controller_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the HDF5 controller class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initHDF5Controller(this,varCont,varNames,filepath) 
    !! HDF5 controller initialization routine - determines variables this controller is responsible for and the output path. The default
    !! path is "./", corresponding to the executable location in the build tree.

    class(HDF5Controller)           ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Reference variable container
    type(StringArray) ,dimension(:) ,intent(in)     :: varNames !! Names of variables this controller is responsible for
    character(*)      ,optional     ,intent(in)     :: filepath !! Filepath this controller will use 

    integer(ik) :: i

    if (assertions) call assert(varCont%isDefined(),"Undefined variable container passed to HDF5 controller constructor")
    
    this%IOVarNames = varNames
    allocate(this%IOVarIndices(size(varNames)))
    allocate(this%IOBuffer(size(varNames)))

    do i = 1, size(varNames)
        this%IOVarIndices(i) = varCont%getVarIndex(varNames(i)%string)
    end do

    if (present(filepath)) then
        this%filepath = filepath
    else
        this%filepath = "./"
    end if

    call this%makeDefined()

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

    character(:) ,allocatable                       :: usedFilename 
    character(len = 80)                             :: tmpstring
    integer(HID_T)                                  :: file
    integer(HID_T)                                  :: dset ,space
    character(len=17) ,parameter                    :: defaultFilename = "ReMKiT1DVarOutput"

    integer(ik)                                     :: i ,hdferr

    integer(HSIZE_T)                  ,dimension(1) :: dsetsize

    if (assertions) then 
        call assert(this%isDefined(),"outputVarsSerial called from undefined HDF5 controller")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to outputVarsSerial")
        call assert(varCont%isDefined(),"Undefined variable container passed to outputVarsSerial")
    end if

    !First gather all data to rank 0
    do i = 1, size(this%IOVarIndices)
        if (varCont%isVarScalar(this%IOVarIndices(i))) then  !Handle scalar var case
            if (mpiCont%getWorldRank() == 0) this%IOBuffer(i)%entry = varCont%variables(this%IOVarIndices(i))%entry
        else
        call mpiCont%gatherVar(varCont%variables(this%IOVarIndices(i))%entry,&
                               this%IOBuffer(i)%entry,&
                               varCont%isVarDist(this%IOVarIndices(i)))
        end if
    end do  

    !Start output from rank 0
    if (mpiCont%getWorldRank() == 0) then 

        call h5open_f(hdferr)

        !Handle file naming
        usedFilename = this%filepath//defaultFilename
        if (present(filename)) usedFilename = this%filepath//filename

        if (present(IDNum)) then 
            write(tmpstring,'(I0)') IDNum 
            usedFilename = usedFilename//"_"//trim(tmpstring)//".h5"
        else 
            usedFilename = usedFilename//".h5"
        end if

        call h5fcreate_f(usedFilename, H5F_ACC_TRUNC_F, file, hdferr)

        do i = 1,size(this%IOVarIndices)
            dsetsize(1) = size(this%IOBuffer(i)%entry)
            call h5screate_simple_f(1, dsetsize, space, hdferr)
            call h5dcreate_f(file, this%IOVarNames(i)%string, int(H5T_IEEE_F64LE, HID_T), space, dset, hdferr)
            call h5dwrite_f(dset, h5kind_to_type(rk,H5_REAL_KIND), this%IOBuffer(i)%entry, dsetsize,hdferr)
            call h5dclose_f(dset , hdferr)
            call h5sclose_f(space, hdferr)
        end do

        call h5fclose_f(file , hdferr)
        call h5close_f(hdferr)
    end if

    !Maybe need a Barrier here?
end subroutine outputVarsSerial
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadVarsSerial(this,mpiCont,varCont,filename) 
    !!Load and scatter variables this controller is responsible for. Default filename is "ReMKiT1DVarInput". If a variable is not
    !!found in the file it is set to 0. 

    class(HDF5Controller)           ,intent(inout)  :: this
    type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
    type(VariableContainer)         ,intent(inout)  :: varCont !! Variable container to be loaded into
    character(*) ,optional          ,intent(in)     :: filename  !! Name of input file

    character(:) ,allocatable                       :: usedFilename 
    character(len = 80)                             :: tmpstring
    integer(HID_T)                                  :: file ,rootID
    integer(HID_T)                                  :: dset ,space
    character(len=16) ,parameter                    :: defaultFilename = "ReMKiT1DVarInput"

    integer(ik)                                     :: i ,hdferr

    integer(HSIZE_T)                  ,dimension(1) :: dsetsize,maxsize
    real(rk) ,allocatable ,dimension(:) :: locBuffer
    logical ,allocatable ,dimension(:) :: varFound

    if (assertions) then 
        call assert(this%isDefined(),"loadVarsSerial called from undefined HDF5 controller")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to loadVarsSerial")
        call assert(varCont%isDefined(),"Undefined variable container passed to loadVarsSerial")
    end if

    !Deallocate allocated buffer variables
    do i = 1,size(this%IOBuffer)
        if (allocated(this%IOBuffer(i)%entry)) deallocate(this%IOBuffer(i)%entry)
    end do
    allocate(varFound(size(this%IOVarIndices)))
    varFound = .false.
    !Start output from rank 0
    if (mpiCont%getWorldRank() == 0) then 

        call h5open_f(hdferr)

        !Handle file naming
        usedFilename = this%filepath//defaultFilename
        if (present(filename)) usedFilename = this%filepath//filename
        usedFilename = usedFilename//".h5"

        call h5fopen_f(usedFilename, H5F_ACC_RDONLY_F, file, hdferr)
        
        do i = 1,size(this%IOVarIndices)
            call h5lexists_f(file,this%IOVarNames(i)%string,varFound(i),hdferr)
            !Add more sofisticated handling of missing variables
            if (varFound(i)) then 
                call h5dopen_f(file,this%IOVarNames(i)%string,dset,hdferr)
                call h5dget_space_f(dset,space,hdferr)
                call h5sget_simple_extent_dims_f(space,dsetsize,maxsize,hdferr)
                allocate(this%IOBuffer(i)%entry(dsetsize(1)))

                call h5dread_f(dset, h5kind_to_type(rk,H5_REAL_KIND), this%IOBuffer(i)%entry, dsetsize,hdferr,H5S_ALL_F,space)

                call h5sclose_f(space, hdferr)
                call h5dclose_f(dset , hdferr)

            end if
        end do

        call h5fclose_f(file , hdferr)
        call h5close_f(hdferr)
    end if

    call mpiCont%barrier()

    do i = 1,size(this%IOVarIndices)
        if (allocated(locBuffer)) deallocate(locBuffer)
        call mpiCont%broadcastLogical(varFound)
        if (varFound(i)) then
            if (varCont%isVarScalar(this%IOVarIndices(i))) then  !Handle scalar var case
                allocate(locBuffer(1))
                if (.not. allocated(this%IOBuffer(i)%entry)) allocate(this%IOBuffer(i)%entry(1))
                if (mpiCont%getWorldRank() == 0) locBuffer = this%IOBuffer(i)%entry
                call mpiCont%broadcastReal(locBuffer)
            else
                call mpiCont%scatterVar(this%IOBuffer(i)%entry,locBuffer,varCont%isVarDist(this%IOVarIndices(i)))
            end if
            if (mpiCont%getRowRank() == 0) then
                varCont%variables(this%IOVarIndices(i))%entry(1:size(locBuffer)) = locBuffer
            end if
            call mpiCont%broadcastVarInRow(varCont,this%IOVarNames(i)%string)
        else 
            varCont%variables(this%IOVarIndices(i))%entry = 0
        end if
    end do

end subroutine loadVarsSerial
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputGridDataSerial(this,mpiCont,gridObj,filename) 
    !! Gather and output grid data. Default filename is "ReMKiT1DGridOutput". The output is serial from processor 0.

    class(HDF5Controller)           ,intent(inout)  :: this
    type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
    type(Grid)                      ,intent(in)     :: gridObj !! Output grid object
    character(*) ,optional          ,intent(in)     :: filename  !! Name of output file

    character(:) ,allocatable                       :: usedFilename 
    integer(HID_T)                                  :: file
    integer(HID_T)                                  :: dset ,space
    character(len=18) ,parameter                    :: defaultFilename = "ReMKiT1DGridOutput"

    integer(ik)                                     :: i ,hdferr

    real(rk)             ,allocatable ,dimension(:) :: xGrid, vGrid
    integer(ik)          ,allocatable ,dimension(:) :: lGrid ,mGrid ,imaginary
    logical              ,allocatable ,dimension(:) :: isImaginary

    integer(HSIZE_T)                  ,dimension(1) :: dsetsize

    if (assertions) then 
        call assert(this%isDefined(),"outputGridDataSerial called from undefined HDF5 controller")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to outputGridDataSerial")
        call assert(gridObj%isDefined(),"Undefined variable container passed to outputGridDataSerial")
    end if

    !Start output from rank 0
    if (mpiCont%getWorldRank() == 0) then 

        call h5open_f(hdferr)

        !Handle file naming
        usedFilename = this%filepath//defaultFilename
        if (present(filename)) usedFilename = this%filepath//filename

        usedFilename = usedFilename//".h5"

        xGrid = gridObj%getXGrid()
        vGrid = gridObj%getVGrid()
        lGrid = gridObj%getLGrid()
        mGrid = gridObj%getMGrid()
        isImaginary = gridObj%getHarmonicIm()

        allocate(imaginary(size(isImaginary)))

        imaginary = 0

        where (isImaginary)
            imaginary = 1
        end where

        call h5fcreate_f(usedFilename, H5F_ACC_TRUNC_F, file, hdferr)

        !Write xgrid data
        dsetsize(1) = size(xGrid)
        call h5screate_simple_f(1, dsetsize, space, hdferr)
        call h5dcreate_f(file, "x", int(H5T_IEEE_F64LE, HID_T), space, dset, hdferr)
        call h5dwrite_f(dset, h5kind_to_type(rk,H5_REAL_KIND), xGrid, dsetsize,hdferr)
        call h5dclose_f(dset , hdferr)
        call h5sclose_f(space, hdferr)

        !Write vgrid data
        dsetsize(1) = size(vGrid)
        call h5screate_simple_f(1, dsetsize, space, hdferr)
        call h5dcreate_f(file, "v", int(H5T_IEEE_F64LE, HID_T), space, dset, hdferr)
        call h5dwrite_f(dset, h5kind_to_type(rk,H5_REAL_KIND), vGrid, dsetsize,hdferr)
        call h5dclose_f(dset , hdferr)
        call h5sclose_f(space, hdferr)

        !Write harmonic data
        dsetsize(1) = size(lGrid)
        call h5screate_simple_f(1, dsetsize, space, hdferr)

        call h5dcreate_f(file, "l", int(H5T_STD_I32LE, HID_T), space, dset, hdferr)
        call h5dwrite_f(dset, h5kind_to_type(ik,H5_INTEGER_KIND), lGrid, dsetsize,hdferr)
        call h5dclose_f(dset , hdferr)

        call h5dcreate_f(file, "m", int(H5T_STD_I32LE, HID_T), space, dset, hdferr)
        call h5dwrite_f(dset, h5kind_to_type(ik,H5_INTEGER_KIND), mGrid, dsetsize,hdferr)
        call h5dclose_f(dset , hdferr)

        call h5dcreate_f(file, "im", int(H5T_STD_I32LE, HID_T), space, dset, hdferr)
        call h5dwrite_f(dset, h5kind_to_type(ik,H5_INTEGER_KIND), imaginary, dsetsize,hdferr)
        call h5dclose_f(dset , hdferr)

        call h5sclose_f(space, hdferr)

        call h5fclose_f(file , hdferr)
        call h5close_f(hdferr)
    end if

end subroutine outputGridDataSerial
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine dumpRestartFiles(this,mpiCont,varCont,filename) 
    !! Each processor dumps all of its variables into separate files. Default filename "restart", appended with processor rank.

    class(HDF5Controller)           ,intent(inout)  :: this
    type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used for output
    character(*) ,optional          ,intent(in)     :: filename !! Name of output file

    character(:) ,allocatable                       :: usedFilename 
    character(len = 80)                             :: tmpstring
    integer(HID_T)                                  :: file
    integer(HID_T)                                  :: dset ,space
    character(len=7) ,parameter                    :: defaultFilename = "restart"

    integer(ik)                                     :: i ,hdferr 

    integer(HSIZE_T)                  ,dimension(1) :: dsetsize

    if (assertions) then 
        call assert(this%isDefined(),"dumpRestartFiles called from undefined HDF5 controller")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to dumpRestartFiles")
        call assert(varCont%isDefined(),"Undefined variable container passed to dumpRestartFiles")
    end if

    call h5open_f(hdferr)

    !Handle file naming
    usedFilename = this%filepath//defaultFilename
    if (present(filename)) usedFilename = this%filepath//filename

    write(tmpstring,'(I0)') mpiCont%getWorldRank()
    usedFilename = usedFilename//"_"//trim(tmpstring)//".h5"

    call h5fcreate_f(usedFilename, H5F_ACC_TRUNC_F, file, hdferr)

    do i = 1,size(varCont%variables)
        dsetsize(1) = size(varCont%variables(i)%entry)
        call h5screate_simple_f(1, dsetsize, space, hdferr)
        call h5dcreate_f(file, varCont%getVarName(i), int(H5T_IEEE_F64LE, HID_T), space, dset, hdferr)
        call h5dwrite_f(dset, h5kind_to_type(rk,H5_REAL_KIND), varCont%variables(i)%entry, dsetsize,hdferr)
        call h5dclose_f(dset , hdferr)
        call h5sclose_f(space, hdferr)
    end do

    call h5fclose_f(file , hdferr)
    call h5close_f(hdferr)

end subroutine dumpRestartFiles
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadRestartFiles(this,mpiCont,varCont,filename) 
    !! Each processor loads all of its variables from files of the form filename_rank.h5. It is assumed that the files exist and that
    !! they are compatible with the variable container

    class(HDF5Controller)           ,intent(inout)  :: this
    type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController used by the HDF5 controller for communication
    type(VariableContainer)         ,intent(inout)  :: varCont !! Variable container used for output
    character(*) ,optional          ,intent(in)     :: filename !! Name of output file

    character(:) ,allocatable                       :: usedFilename 
    character(len = 80)                             :: tmpstring
    integer(HID_T)                                  :: file
    integer(HID_T)                                  :: dset
    character(len=7) ,parameter                     :: defaultFilename = "restart"

    integer(ik)                                     :: i ,hdferr 

    integer(HSIZE_T)                  ,dimension(1) :: dsetsize

    if (assertions) then 
        call assert(this%isDefined(),"dumpRestartFiles called from undefined HDF5 controller")
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to dumpRestartFiles")
        call assert(varCont%isDefined(),"Undefined variable container passed to dumpRestartFiles")
    end if

    call h5open_f(hdferr)

    !Handle file naming
    usedFilename = this%filepath//defaultFilename
    if (present(filename)) usedFilename = this%filepath//filename

    write(tmpstring,'(I0)') mpiCont%getWorldRank()
    usedFilename = usedFilename//"_"//trim(tmpstring)//".h5"

    call h5fopen_f(usedFilename, H5F_ACC_RDONLY_F, file, hdferr)

    do i = 1,size(varCont%variables)
        dsetsize(1) = size(varCont%variables(i)%entry)
        call h5dopen_f(file, varCont%getVarName(i), dset, hdferr)
        call h5dread_f(dset, int(H5T_IEEE_F64LE, HID_T), varCont%variables(i)%entry, dsetsize,hdferr)
        call h5dclose_f(dset , hdferr)
    end do

    call h5fclose_f(file , hdferr)
    call h5close_f(hdferr)

end subroutine loadRestartFiles
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule HDF5_controller_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
