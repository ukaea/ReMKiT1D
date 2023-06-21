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
module json_controller_class
    !! Contains json-fortran routines for reading/writing config files

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use god_objects                 ,only: object
    use support_types               ! Using all named scalars and arrays
    use mpi_controller_class        ,only: MPIController
    use json_module
    implicit none 
    public

    type ,public :: JSONController
        !! Object responsible for reading and writing to JSON files

        type(json_file)           ,private :: file !! json file handle
        logical                   ,private :: fileOpen = .false. !! True if this controller has opened a file
        character(:) ,allocatable ,private :: alternativeJSONFilepath !! If allocated will be used instead of default config.json

        contains 

        procedure ,public :: loadFile 
        procedure ,public :: closeFile 

        procedure ,public :: setAlternativeJSONPath
        procedure ,public :: getAlternativeJSONPath

        procedure ,public :: loadArrayParams
        procedure ,public :: loadScalarParams 

        procedure ,public :: outputArrayParamsToFile
        procedure ,public :: outputScalarParamsToFile

        generic ,public :: load => loadNamedIntArrays, loadNamedInts, loadNamedLogicals, loadNamedLogicalArrays, &
                                    loadNamedRealArrays, loadNamedReals, loadNamedStringArrays, loadNamedStrings,loadArrayParams,&
                                    loadScalarParams

        generic ,public :: output => outputNamedIntArraysToFile, outputNamedIntsToFile, outputNamedLogicalsToFile, &
        outputNamedLogicalArraysToFile, outputNamedRealArraysToFile, outputNamedRealsToFile, outputNamedStringArraysToFile,&
         outputNamedStringsToFile,outputArrayParamsToFile,outputScalarParamsToFile

        procedure ,public :: loadNamedReals
        procedure ,public :: loadNamedInts
        procedure ,public :: loadNamedLogicals
        procedure ,public :: loadNamedStrings

        procedure ,public :: loadNamedRealArrays
        procedure ,public :: loadNamedIntArrays
        procedure ,public :: loadNamedLogicalArrays
        procedure ,public :: loadNamedStringArrays

        procedure ,public :: outputNamedRealsToFile
        procedure ,public :: outputNamedIntsToFile
        procedure ,public :: outputNamedLogicalsToFile
        procedure ,public :: outputNamedStringsToFile

        procedure ,public :: outputNamedRealArraysToFile
        procedure ,public :: outputNamedIntArraysToFile
        procedure ,public :: outputNamedLogicalArraysToFile
        procedure ,public :: outputNamedStringArraysToFile

    end type JSONController
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadFile(this,mpiCont,filepath) 
        !! Loads json file on rank 0. The default filepath here is "./config.json".

        class(JSONController)           ,intent(inout)  :: this
        type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController object to be used by the JSONController for communication
        character(*)      ,optional     ,intent(in)     :: filepath !! Non-default filepath

    end subroutine loadFile
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine closeFile(this,mpiCont,saveFile,filepath) 
        !! Closes currently open json file. If saveFile is true, saves the file being worked on before closing. 
        !! The default filepath here is "./config.json".

        class(JSONController)           ,intent(inout)  :: this
        type(MPIController)             ,intent(inout)  :: mpiCont !! MPIController object to be used by the JSONController for communication
        logical  ,optional              ,intent(in)     :: saveFile !! True if the file should be saved before closing
        character(*)      ,optional     ,intent(in)     :: filepath !! Non-default filepath for saving

    end subroutine closeFile
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setAlternativeJSONPath(this,filepath) 
        !! Set alternative default JSON filepath

        class(JSONController)           ,intent(inout)  :: this
        character(*)                    ,intent(in)     :: filepath 

    end subroutine setAlternativeJSONPath
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function getAlternativeJSONPath(this) result(filepath)
        !! Get alternative default JSON filepath

        class(JSONController)           ,intent(in)     :: this
        character(:) ,allocatable                       :: filepath 

    end function getAlternativeJSONPath
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadScalarParams(this,vars) 
        !! Load named scalar parameters by calling individual load routines.  

        class(JSONController)           ,intent(inout)  :: this
        type(NamedScalarContainer)      ,intent(inout)  :: vars !! Values to load

    end subroutine loadScalarParams
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadArrayParams(this,vars) 
        !!Load named array parameters by calling individual load routines.  

        class(JSONController)           ,intent(inout)  :: this
        type(NamedArrayContainer)       ,intent(inout)  :: vars !! Values to load

    end subroutine loadArrayParams
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedReals(this,vars) 
        !! Load named reals from currently open json file. If a variable isn't found, the
        !! passed value is not modified. 

        class(JSONController)           ,intent(inout)  :: this
        type(NamedReal) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedReals
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedInts(this,vars) 
        !! Load named int from currently open json file. If a variable isn't found, the
        !! passed value is not modified. 

        class(JSONController)            ,intent(inout)  :: this
        type(NamedInteger) ,dimension(:) ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedInts
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedLogicals(this,vars) 
        !! Load named logicals from currently open json file. If a variable isn't found, 
        !! the passed value is not modified. 

        class(JSONController)              ,intent(inout)  :: this
        type(NamedLogical) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedLogicals
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedStrings(this,vars) 
        !! Load named strings from currently open json file. If a variable isn't found, the
        !! passed value is not modified. 

        class(JSONController)           ,intent(inout)  :: this
        type(NamedString) ,dimension(:) ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedStrings
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedRealArrays(this,vars) 
        !! Load named real arrays from currently open json file. If a variable isn't found, 
        !! the passed value is not modified. 

        class(JSONController)                ,intent(inout)  :: this
        type(NamedRealArray) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedRealArrays
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedIntArrays(this,vars) 
        !! Load named int arrays from currently open json file. If a variable isn't found, 
        !! the passed value is not modified. 

        class(JSONController)                    ,intent(inout)  :: this
        type(NamedIntegerArray)  ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedIntArrays
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedLogicalArrays(this,vars) 
        !! Load named logical arrays from currently open json file. If a variable isn't 
        !! found, the passed value is not modified. 

        class(JSONController)                   ,intent(inout)  :: this
        type(NamedLogicalArray) ,dimension(:)   ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedLogicalArrays
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine loadNamedStringArrays(this,vars)
        !!  Load named string arrays from currently open json file. If a variable isn't 
        !! found, the passed value is not modified. 

        class(JSONController)                   ,intent(inout)  :: this
        type(NamedStringArray) ,dimension(:)    ,intent(inout)  :: vars !! Values to load

    end subroutine loadNamedStringArrays
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedRealsToFile(this,vars) 
    !! Outputs named reals to json file.

    class(JSONController)           ,intent(inout)  :: this
    type(NamedReal) ,dimension(:)   ,intent(inout)  :: vars !! Values to output

end subroutine outputNamedRealsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputScalarParamsToFile(this,vars) 
    !! Outputs named scalar parameters to json file. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedScalarContainer)              ,intent(inout)  :: vars !! Values to output

end subroutine outputScalarParamsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputArrayParamsToFile(this,vars) 
    !! Outputs named array parameters to json file. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedArrayContainer)               ,intent(inout)  :: vars !! Values to output

end subroutine outputArrayParamsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedIntsToFile(this,vars) 
    !! Outputs named ints to json file.

    class(JSONController)              ,intent(inout)  :: this
    type(NamedInteger) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

end subroutine outputNamedIntsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedLogicalsToFile(this,vars) 
    !! Outputs named logicals to json file.

    class(JSONController)              ,intent(inout)  :: this
    type(NamedLogical) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

end subroutine outputNamedLogicalsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedStringsToFile(this,vars) 
    !! Outputs named strings to json file.

    class(JSONController)             ,intent(inout)  :: this
    type(NamedString) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

end subroutine outputNamedStringsToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedRealArraysToFile(this,vars) 
    !! Outputs named real arrays to json file.

    class(JSONController)                ,intent(inout)  :: this
    type(NamedRealArray) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

end subroutine outputNamedRealArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedIntArraysToFile(this,vars) 
    !! Outputs named int arrays to json file. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedIntegerArray) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

end subroutine outputNamedIntArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedLogicalArraysToFile(this,vars) 
    !! Outputs named logical arrays to json file. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedLogicalArray) ,dimension(:)   ,intent(inout)  :: vars !! Values ot output

end subroutine outputNamedLogicalArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine outputNamedStringArraysToFile(this,vars) 
    !! Outputs named string arrays to json file. 

    class(JSONController)                   ,intent(inout)  :: this
    type(NamedStringArray) ,dimension(:)   ,intent(inout)  :: vars 

end subroutine outputNamedStringArraysToFile
!-----------------------------------------------------------------------------------------------------------------------------------
end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module json_controller_class
!-----------------------------------------------------------------------------------------------------------------------------------