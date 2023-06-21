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
module status_printing
    !! author: Stefan Mijin 
    !!
    !! Contains a simple status printing utility for basic console output

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging ,assertions
    use assertion_utility           ,only: assert
    use iso_fortran_env             ,only: error_unit, output_unit

    implicit none
    private

    integer(ik) ,private :: rank !! MPI rank for printing purposes

    logical     ,private :: printerReady !! True if the printer can be used

    public :: preparePrinter ,printNamedValue ,printError ,printMessage

!-----------------------------------------------------------------------------------------------------------------------------------
    interface printNamedValue

        module subroutine printNamedInt(name,input,allProcs)

            character(*)      ,intent(in) :: name
            integer(ik)       ,intent(in) :: input

            logical ,optional ,intent(in) :: allProcs

        end subroutine printNamedInt

        module subroutine printNamedReal(name,input,allProcs)
    
            character(*) ,intent(in) :: name
            real(rk)     ,intent(in) :: input

            logical ,optional ,intent(in) :: allProcs

        end subroutine printNamedReal

        module subroutine printNamedLogical(name,input,allProcs)
    
            character(*) ,intent(in) :: name
            logical      ,intent(in) :: input

            logical ,optional ,intent(in) :: allProcs

        end subroutine printNamedLogical

    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
    contains
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine preparePrinter(rankIn)
        !! Prepare printer utility for use - currently only sets rank

        integer(ik) ,intent(in) :: rankIn 

        rank = rankIn
        printerReady = .true.

    end subroutine preparePrinter
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine printNamedInt(name,input,allProcs)
        !! Generic printNamedValue interface component - prints name: input for integer input. Prefaces with rank number if allProcs is true

        character(*) ,intent(in) :: name
        integer(ik)  ,intent(in) :: input

        logical ,optional ,intent(in) :: allProcs

        logical :: aP

        if (assertions) call assert(printerReady,"Attempted to use printer module before printer was prepared")

        aP  = .false. 
        if (present(allProcs)) aP = allProcs
        if (aP) then 
            write(output_unit,*) "PID",rank," - ",name,": ",input
            write(output_unit,*)
        else if (rank == 0) then
            write(output_unit,*) name,": ",input
            write(output_unit,*)
        end if

    end subroutine printNamedInt
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine printNamedReal(name,input,allProcs)
        !! Generic printNamedValue interface component - prints name: input for real input. Prefaces with rank number if allProcs is true

        character(*) ,intent(in) :: name
        real(rk)     ,intent(in) :: input

        logical ,optional ,intent(in) :: allProcs

        logical :: aP

        if (assertions) call assert(printerReady,"Attempted to use printer module before printer was prepared")

        aP  = .false. 
        if (present(allProcs)) aP = allProcs
        if (aP) then 
            write(output_unit,*) "PID",rank," - ",name,": ",input
            write(output_unit,*)
        else if (rank == 0) then
            write(output_unit,*) name,": ",input
            write(output_unit,*)
        end if

    end subroutine printNamedReal
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine printNamedLogical(name,input,allProcs)
        !! Generic printNamedValue interface component - prints name: input for logical input. Prefaces with rank number if allProcs is true

        character(*) ,intent(in) :: name
        logical      ,intent(in) :: input

        logical ,optional ,intent(in) :: allProcs

        logical :: aP

        if (assertions) call assert(printerReady,"Attempted to use printer module before printer was prepared")

        aP  = .false. 
        if (present(allProcs)) aP = allProcs
        if (aP) then 
            write(output_unit,*) "PID",rank," - ",name,": ",input
            write(output_unit,*)
        else if (rank == 0) then
            write(output_unit,*) name,": ",input
            write(output_unit,*)
        end if

    end subroutine printNamedLogical
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine printMessage(msg,allProcs)
        !! Print message to output_unit - if allProcs true prefaces with rank number

        character(*) ,intent(in) :: msg

        logical ,optional ,intent(in) :: allProcs

        logical :: aP

        if (assertions) call assert(printerReady,"Attempted to use printer module before printer was prepared")

        aP  = .false. 
        if (present(allProcs)) aP = allProcs
        if (aP) then 
            write(output_unit,*) "PID",rank," - ",msg
            write(output_unit,*)
        else if (rank == 0) then
            write(output_unit,*) msg
            write(output_unit,*)
        end if

    end subroutine printMessage
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine printError(msg)
        !! Print error to error_unit on all processes

        character(*) ,intent(in) :: msg

        if (assertions) call assert(printerReady,"Attempted to use printer module before printer was prepared")

        write(output_unit,*) "PID",rank," - ERROR:",msg
        write(output_unit,*)

    end subroutine printError
!-----------------------------------------------------------------------------------------------------------------------------------
 end module status_printing
!-----------------------------------------------------------------------------------------------------------------------------------
 