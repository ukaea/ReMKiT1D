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
submodule (stencil1d_class) stencil1d_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the stencil1d class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initStencil(this,rawStencil,fixedStencil) 
    !! Stencil1D object initialization routine

    class(Stencil1D)                    ,intent(inout)  :: this
    integer(ik) ,optional ,dimension(:) ,intent(in)     :: rawStencil   !! Optional raw stencil, defaults to [0] - a diagonal stencil.
    type(IntArray) ,optional ,dimension(:) ,intent(in)  :: fixedStencil !! Optional fixed stencil 

    if (present(rawStencil)) then 
        this%rawStencil = rawStencil 
    else
        this%rawStencil = [0]
    end if

    if (present(fixedStencil)) then 
        this%fixedStencil = fixedStencil
    end if

call this%makeDefined()

end subroutine initStencil
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function mapCoords(this,inputCoord,dimSize,periodic) result(output)
    !! Stencil1D coordinate mapping routine

    class(Stencil1D)       ,intent(in) :: this
    integer(ik)            ,intent(in) :: inputCoord  !! Input coordinate value
    integer(ik)            ,intent(in) :: dimSize !! Size of dimension in which mapping is done
    logical      ,optional ,intent(in) :: periodic !! True if dimension is periodic
    integer(ik) ,allocatable ,dimension(:)       :: output 

    logical :: isPeriodic 

    if (this%isStencilFixed()) then 
        output = this%fixedStencil(inputCoord)%entry
        return 
    end if 
    
    isPeriodic = .false. 
    if (present(periodic)) isPeriodic = periodic 

    output = pack(inputCoord + this%rawStencil,this%getMask(inputCoord,dimSize,isPeriodic))

    !Correct for periodic dimension
    if (isPeriodic) then
        where (output <= 0)
            output = output + dimSize
        else where (output > dimSize)
            output = output - dimSize
        end where
    end if

end function mapCoords
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMask(this,coord,dimSize,periodic) result(res)
    !! Get logical mask for included stencil points for given coordinate and dimension size. If periodic, the stencil dimension is treated as being periodic with periodicity equal to dimSize.

    class(Stencil1D)                    ,intent(in) :: this
    integer(ik)                         ,intent(in) :: coord  
    integer(ik)                         ,intent(in) :: dimSize
    logical      ,optional              ,intent(in) :: periodic
    logical ,allocatable ,dimension(:)              :: res 

    logical :: isPeriodic 

    isPeriodic = .false. 
    if (present(periodic)) isPeriodic = periodic 

    allocate(res(size(this%rawStencil)))

    res = withinBounds(this%rawStencil+coord,1,dimSize)
    res = res .or. (isPeriodic .and. (coord + this%rawStencil <=0 ))
    res = res .or. (isPeriodic .and. (coord + this%rawStencil > dimSize))

end function getMask
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getStencilDims(this) result(dim)
    !! Returns size of stencil

    class(Stencil1D)       ,intent(in) :: this
    integer(ik)                        :: dim 

    dim = size(this%rawStencil)

end function getStencilDims
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isStencilFixed(this) result(stencilIsFixed)
    !! Check if stencil is fixed

    class(Stencil1D) ,intent(in) :: this
    logical                       :: stencilIsFixed 

    if (assertions) call assertPure(this%isDefined(),"isStencilFixed called on undefined Stencil1D")

    stencilIsFixed = allocated(this%fixedStencil)

end function isStencilFixed
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFixedStencil(this) result(fixedStencil)
    !! Return values of fixed stencil

    class(Stencil1D)                          ,intent(in) :: this
    type(IntArray) ,allocatable ,dimension(:)             :: fixedStencil  

    if (assertions) then 
        call assertPure(this%isDefined(),"getFixedStencil called on undefined Stencil1D")
        call assertPure(this%isStencilFixed(),"getFixedStencil called on Stencil1D with no fixed stencil")
    end if

    fixedStencil = this%fixedStencil

end function getFixedStencil
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule stencil1d_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
