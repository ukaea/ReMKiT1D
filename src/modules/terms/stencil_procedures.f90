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
submodule (stencil_class) stencil_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the stencil class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initStencil(this,xStencil,vStencil,hStencil,mapToDist,xPeriodic,vStencilFixed) 
    !! Stencil object initialization routine

    class(stencil)                      ,intent(inout)  :: this
    integer(ik) ,optional ,dimension(:) ,intent(in)     :: xStencil !! Stencil in x-direction
    integer(ik) ,optional ,dimension(:) ,intent(in)     :: vStencil !! Stencil in v-direction 
    integer(ik) ,optional ,dimension(:) ,intent(in)     :: hStencil !! l harmonic stencil
    logical     ,optional               ,intent(in)     :: mapToDist !! Set to true if stencil maps to harmonic/velocity space
    logical     ,optional               ,intent(in)     :: xPeriodic !! Set to true if x-grid is periodic
    type(IntArray) ,optional ,dimension(:) ,intent(in)  :: vStencilFixed !! Optional fixed vStencil 

    logical :: toDist
    if (assertions) then 

        if (present(vStencilFixed)) then 
                call assertPure(present(mapToDist),"If vStencilFixed is present mapToDist must be set to true")
                call assertPure(mapToDist,"If vStencilFixed is present mapToDist must be set to true")

        end if
    end if
    if (present(xStencil)) then 
        call this%xStencil%init(xStencil) 
    else
        call this%xStencil%init() 
    end if

    toDist = .false.
    if (present(mapToDist)) toDist = mapToDist 

    if (toDist) then 
        allocate(this%vStencil)
        allocate(this%hStencil)

        if (present(vStencilFixed)) then
            call this%vStencil%init(fixedStencil=vStencilFixed)
        else if (present(vStencil)) then 
            call this%vStencil%init(vStencil)
        else
            call this%vStencil%init()
        end if

        if (present(hStencil)) then 
            call this%hStencil%init(hStencil)
        else
            call this%hStencil%init()
        end if

    end if

    this%xPeriodic = .false. 
    if (present(xPeriodic)) this%xPeriodic = xPeriodic
    
    call this%makeDefined()

end subroutine initStencil
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function mapCoords(this,gridObj,inputArray) result(output)
    !! Stencil coordinate mapping routine

    class(stencil)                      ,intent(in) :: this
    type(Grid)                          ,intent(in) :: gridObj !! Grid used to construct the mapping
    integer(ik)    ,dimension(:)        ,intent(in) :: inputArray  !! Input array of coordinates (size 1 - [x], or size 3 [x,h,v])
    type(IntArray) ,allocatable ,dimension(:)       :: output 

    if (assertions) call assertPure(gridObj%isDefined(),"Undefined grid object passed to stencil mapping routine")

    if (allocated(this%vStencil)) then    
        allocate(output(3))
    else 
        allocate(output(1))
    end if

    output(1)%entry = this%xStencil%mapCoords(inputArray(1),gridObj%getNumX(),this%xPeriodic)

    if (size(output)>1) then   

        if (size(inputArray) == 1) then 
            output(2)%entry = this%hStencil%mapCoords(0,gridObj%getNumH())
            output(3)%entry = this%vStencil%mapCoords(0,gridObj%getNumV())
        else
            output(2)%entry = this%hStencil%mapCoords(inputArray(2),gridObj%getNumH())
            output(3)%entry = this%vStencil%mapCoords(inputArray(3),gridObj%getNumV())
        end if

    end if

end function mapCoords
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule stencil_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
