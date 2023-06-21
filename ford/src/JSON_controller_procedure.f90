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
submodule (JSON_controller_class) JSON_controller_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the JSON controller class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadFile(this,mpiCont,filepath) 
    !! Loads json file on rank 0. The default filepath here is "./config.json".

    class(JSONController)           ,intent(inout)  :: this
    type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController object to be used by the JSONController for communication
    character(*)      ,optional     ,intent(in)     :: filepath !! Non-default filepath

    character(:) ,allocatable                       :: usedFilepath 
    character(len=13) ,parameter                    :: defaultFilepath = "./config.json"
    type(json_core)                                 :: core
    type(json_value)                       ,pointer :: p

    integer(ik) :: currentRank

    character(len=65)   :: intToStrBuffer


    if (assertions) call assert(mpiCont%isDefined(),"Undefined mpi controller passed to loadFile routine of json controller")

    currentRank = 0 

    usedFilepath = defaultFilepath
    if (allocated(this%alternativeJSONFilepath)) usedFilepath = this%alternativeJSONFilepath
    if (present(filepath)) usedFilepath = filepath

    if (mpiCont%getWorldRank() == 0) then 
        print*,"Loading config file "//usedFilepath//" on MPI rank 0"
            
        call this%file%initialize(comment_char='/')
        call this%file%load(filename=usedFilepath)

        if (this%file%failed()) then 
            call core%initialize()
            call core%create_object(p,'')
            call core%print(p,usedFilepath)
            call core%destroy(p)
            call this%file%load(filename=usedFilepath)
        end if

        this%fileOpen = .true.
    end if
    call mpiCont%barrier() 

    if (mpiCont%getWorldRank() /= 0) then 
        intToStrBuffer = ""
        write(intToStrBuffer,'(I0)') mpiCont%getWorldRank()
        print*,"Loading config file "//usedFilepath//" on MPI rank "//trim(intToStrBuffer)
        
        call this%file%initialize(comment_char='/')
        call this%file%load(filename=usedFilepath)

        this%fileOpen = .true.
    end if
    call mpiCont%barrier() 
end subroutine loadFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine closeFile(this,mpiCont,saveFile,filepath) 
    !! Closes currently open json file. If saveFile is true, saves the file being worked on before closing. 
    !! The default filepath here is "./config.json".

    class(JSONController)           ,intent(inout)  :: this
    type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController object to be used by the JSONController for communication
    logical  ,optional              ,intent(in)     :: saveFile !! True if the file should be saved before closing
    character(*)      ,optional     ,intent(in)     :: filepath !! Non-default filepath for saving

    character(:) ,allocatable                       :: usedFilepath 
    character(len=13) ,parameter                    :: defaultFilepath = "./config.json"

    if (assertions) call assert(mpiCont%isDefined(),"Undefined mpi controller passed to closeFile routine of json controller")

    if (this%fileOpen) then

        if (mpiCont%getWorldRank() == 0) then

            usedFilepath = defaultFilepath
            if (allocated(this%alternativeJSONFilepath)) usedFilepath = this%alternativeJSONFilepath
            if (present(filepath)) usedFilepath = filepath

            if (present(saveFile)) then
                if (saveFile) call this%file%print(usedFilepath)
            end if

        end if

        call this%file%destroy()

        this%fileOpen = .false.

    end if
    call mpiCont%barrier()

end subroutine closeFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadScalarParams(this,vars) 
    !! Load named scalar parameters by calling individual load routines. 

    class(JSONController)           ,intent(inout)  :: this
    type(NamedScalarContainer)      ,intent(inout)  :: vars !! Values to load

    call this%load(vars%realData)
    call this%load(vars%intData)
    call this%load(vars%logicalData)
    call this%load(vars%stringData)

end subroutine loadScalarParams
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setAlternativeJSONPath(this,filepath) 
    !! Set alternative default JSON filepath

    class(JSONController)           ,intent(inout)  :: this
    character(*)                    ,intent(in)     :: filepath 

    this%alternativeJSONFilepath = filepath

end subroutine setAlternativeJSONPath
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getAlternativeJSONPath(this) result(filepath)
    !! Get alternative default JSON filepath

    class(JSONController)           ,intent(in)     :: this
    character(:) ,allocatable                       :: filepath 

    allocate(filepath,source=this%alternativeJSONFilepath)

end function getAlternativeJSONPath
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadArrayParams(this,vars) 
    !!Load named array parameters by calling individual load routines.  

    class(JSONController)           ,intent(inout)  :: this
    type(NamedArrayContainer)       ,intent(inout)  :: vars !! Values to load

    call this%load(vars%realData)
    call this%load(vars%intData)
    call this%load(vars%logicalData)
    call this%load(vars%stringData)

end subroutine loadArrayParams
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedReals(this,vars) 
    !! Load named reals from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't found, the
    !! passed value is not modified. 

    class(JSONController)           ,intent(inout)  :: this
    type(NamedReal) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    real(rk) ,allocatable ,dimension(:)   :: buffer 
    integer(ik)                           :: i
    logical                               :: found
    
    allocate(buffer(size(vars)))

    if (assertions) call assert(this%fileOpen,"loadNamedReals called with no open json file")

    do i = 1, size(vars)

        buffer(i) = vars(i)%value   !Save previous value in case the variable is not found
        call this%file%get(vars(i)%name,vars(i)%value,found)
        if (.not. found) vars(i)%value = buffer(i)

    end do

end subroutine loadNamedReals
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedInts(this,vars) 
    !! Load named int from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't found, the
    !! passed value is not modified. 

    class(JSONController)            ,intent(inout)  :: this
    type(NamedInteger) ,dimension(:) ,intent(inout)  :: vars !! Values to load

    integer(ik) ,allocatable ,dimension(:)   :: buffer 
    integer(ik)                              :: i
    logical                                  :: found
    
    allocate(buffer(size(vars)))

    if (assertions) call assert(this%fileOpen,"loadNamedIntegers called with no open json file")

    do i = 1, size(vars)

        buffer(i) = vars(i)%value   !Save previous value in case the variable is not found
        call this%file%get(vars(i)%name,vars(i)%value,found)
        if (.not. found) vars(i)%value = buffer(i)

    end do

end subroutine loadNamedInts
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedLogicals(this,vars) 
    !! Load named logicals from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't found, 
    !! the passed value is not modified. 

    class(JSONController)              ,intent(inout)  :: this
    type(NamedLogical) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    logical ,allocatable ,dimension(:)   :: buffer 
    integer(ik)                          :: i
    logical                              :: found
    
    allocate(buffer(size(vars)))

    if (assertions) call assert(this%fileOpen,"loadNamedLogicals called with no open json file")

    do i = 1, size(vars)

        buffer(i) = vars(i)%value   !Save previous value in case the variable is not found
        call this%file%get(vars(i)%name,vars(i)%value,found)
        if (.not. found) vars(i)%value = buffer(i)

    end do

end subroutine loadNamedLogicals
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedStrings(this,vars) 
    !! Load named strings from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't found, the
    !! passed value is not modified. 

    class(JSONController)           ,intent(inout)  :: this
    type(NamedString) ,dimension(:) ,intent(inout)  :: vars !! Values to load

    character(:) ,allocatable            :: buffer 
    integer(ik)                          :: i
    logical                              :: found
    
    call assert(this%fileOpen,"loadNamedStrings called with no open json file")

    do i = 1, size(vars)

        buffer = vars(i)%value   !Save previous value in case the variable is not found
        call this%file%get(vars(i)%name,vars(i)%value,found)
        if (.not. found) vars(i)%value = buffer

    end do

end subroutine loadNamedStrings
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedRealArrays(this,vars) 
    !! Load named real arrays from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't found, 
    !! the passed value is not modified. 

    class(JSONController)                ,intent(inout)  :: this
    type(NamedRealArray) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    type(realArray) ,allocatable ,dimension(:)  :: buffer 
    integer(ik)                                 :: i
    logical                                     :: found
    
    allocate(buffer(size(vars)))

    if (assertions) call assert(this%fileOpen,"loadNamedRealArrays called with no open json file")

    do i = 1, size(vars)

        buffer(i)%entry = vars(i)%values   !Save previous value in case the variable is not found
        call this%file%get(vars(i)%name,vars(i)%values,found)
        if (.not. found) vars(i)%values = buffer(i)%entry

    end do

end subroutine loadNamedRealArrays
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedIntArrays(this,vars) 
    !! Load named int arrays from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't found, 
    !! the passed value is not modified. 

    class(JSONController)                    ,intent(inout)  :: this
    type(NamedIntegerArray)  ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    type(intArray) ,allocatable ,dimension(:)   :: buffer 
    integer(ik)                                 :: i
    logical                                     :: found
    
    allocate(buffer(size(vars)))

    if (assertions) call assert(this%fileOpen,"loadNamedIntArrays called with no open json file")

    do i = 1, size(vars)

        buffer(i)%entry = vars(i)%values   !Save previous value in case the variable is not found
        call this%file%get(vars(i)%name,vars(i)%values,found)
        if (.not. found) vars(i)%values = buffer(i)%entry

    end do


end subroutine loadNamedIntArrays
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedLogicalArrays(this,vars) 
    !! Load named logical arrays from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't 
    !! found, the passed value is not modified. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedLogicalArray) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    type(logicalArray) ,allocatable ,dimension(:)  :: buffer 
    integer(ik)                                    :: i
    logical                                        :: found
    
    allocate(buffer(size(vars)))

    if (assertions) call assert(this%fileOpen,"loadNamedLogicalArrays called with no open json file")

    do i = 1, size(vars)

        buffer(i)%entry = vars(i)%values   !Save previous value in case the variable is not found
        call this%file%get(vars(i)%name,vars(i)%values,found)
        if (.not. found) vars(i)%values = buffer(i)%entry

    end do

end subroutine loadNamedLogicalArrays
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine loadNamedStringArrays(this,vars)
    !!  Load named string arrays from currently open json file on rank 0, then broadcast them to all processors. If a variable isn't 
    !! found, the passed value is not modified. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedStringArray) ,dimension(:)    ,intent(inout)  :: vars !! Values to load

    integer(ik)                                    :: i,j
    logical                                        :: found
    character(len=65) ,allocatable ,dimension(:)   :: readBuffer
    
    if (assertions) call assert(this%fileOpen,"loadNamedStringArrays called with no open json file")
    
    do i = 1, size(vars) 
        call this%file%get(vars(i)%name,readBuffer,found)
        if (found) then 
            if (allocated(vars(i)%values)) deallocate(vars(i)%values)
            allocate(vars(i)%values(size(readBuffer)))
            do j = 1, size(readBuffer)
                vars(i)%values(j)%string = trim(readBuffer(j))
            end do
        end if
    end do

end subroutine loadNamedStringArrays
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputScalarParamsToFile(this,vars) 
    !! Outputs named scalar parameters to json file from rank 0. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedScalarContainer)              ,intent(inout)  :: vars !! Values to output

    call this%output(vars%realData)
    call this%output(vars%intData)
    call this%output(vars%logicalData)
    call this%output(vars%stringData)

end subroutine outputScalarParamsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputArrayParamsToFile(this,vars) 
    !! Outputs named array parameters to json file from rank 0. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedArrayContainer)               ,intent(inout)  :: vars !! Values to output

    call this%output(vars%realData)
    call this%output(vars%intData)
    call this%output(vars%logicalData)
    call this%output(vars%stringData)

end subroutine outputArrayParamsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedRealsToFile(this,vars) 
    !! Outputs named reals to json file from rank 0. First attempts to update the variable, and if it is not found it is added.

    class(JSONController)           ,intent(inout)  :: this
    type(NamedReal) ,dimension(:)   ,intent(inout)  :: vars !! Values to output

    integer(ik)                                     :: i

    if (assertions) call assert(this%fileOpen,"outputNamedRealsToFile called with no open json file")

    do i = 1, size(vars)

        call this%file%add(vars(i)%name,vars(i)%value)

    end do
        

end subroutine outputNamedRealsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedIntsToFile(this,vars) 
    !! Outputs named ints to json file from rank 0. First attempts to update the variable, and if it is not found it is added.

    class(JSONController)              ,intent(inout)  :: this
    type(NamedInteger) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

    integer(ik)                                     :: i

    if (assertions) call assert(this%fileOpen,"outputNamedIntsToFile called with no open json file")

    do i = 1, size(vars)

        
        call this%file%add(vars(i)%name,vars(i)%value)

    end do

end subroutine outputNamedIntsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedLogicalsToFile(this,vars) 
    !! Outputs named logicals to json file from rank 0. First attempts to update the variable, and if it is not found it is added.

    class(JSONController)              ,intent(inout)  :: this
    type(NamedLogical) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

    integer(ik)                                     :: i

    if (assertions) call assert(this%fileOpen,"outputNamedLogicalsToFile called with no open json file")

    do i = 1, size(vars)

        call this%file%add(vars(i)%name,vars(i)%value)

    end do

end subroutine outputNamedLogicalsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedStringsToFile(this,vars) 
    !! Outputs named strings to json file from rank 0. First attempts to update the variable, and if it is not found it is added.

    class(JSONController)             ,intent(inout)  :: this
    type(NamedString) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

    integer(ik)                                     :: i

    if (assertions) call assert(this%fileOpen,"outputNamedStringsToFile called with no open json file")

    do i = 1, size(vars)

        call this%file%add(vars(i)%name,vars(i)%value)

    end do

end subroutine outputNamedStringsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedRealArraysToFile(this,vars) 
    !! Outputs named real arrays to json file from rank 0.

    class(JSONController)                ,intent(inout)  :: this
    type(NamedRealArray) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

    integer(ik)                                     :: i

    if (assertions) call assert(this%fileOpen,"outputNamedRealArraysToFile called with no open json file")

    do i = 1, size(vars)

        call this%file%add(vars(i)%name,vars(i)%values)

    end do

end subroutine outputNamedRealArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedIntArraysToFile(this,vars) 
    !! Outputs named int arrays to json file from rank 0. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedIntegerArray) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

    integer(ik)                                     :: i

    if (assertions) call assert(this%fileOpen,"outputNamedIntArraysToFile called with no open json file")

    do i = 1, size(vars)

        call this%file%add(vars(i)%name,vars(i)%values)

    end do

end subroutine outputNamedIntArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedLogicalArraysToFile(this,vars) 
    !! Outputs named logical arrays to json file from rank 0. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedLogicalArray) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

    integer(ik)                                     :: i

    if (assertions) call assert(this%fileOpen,"outputNamedLogicalArraysToFile called with no open json file")

    do i = 1, size(vars)

        call this%file%add(vars(i)%name,vars(i)%values)

    end do

end subroutine outputNamedLogicalArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedStringArraysToFile(this,vars) 
    !! Outputs named string arrays to json file from rank 0. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedStringArray) ,dimension(:)   ,intent(inout)  :: vars 

    integer(ik)                                             :: i,j

    character(len = 80)                             :: tmpstring

    if (assertions) call assert(this%fileOpen,"outputNamedStringArraysToFile called with no open json file")

    do i = 1, size(vars)
        do j = 1,size(vars(i)%values)
            write(tmpstring,'(I0)') j
            call this%file%add(vars(i)%name//"["//trim(tmpstring)//"]",vars(i)%values(j)%string)
        end do
    end do

end subroutine outputNamedStringArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule JSON_controller_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
