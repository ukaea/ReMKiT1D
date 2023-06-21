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
submodule (grid_class) grid_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the grid class


implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initGrid(this,x,v,maxL,maxM) 
    !! Grid initialization routine

    class(Grid)           ,intent(inout)  :: this
    real(rk)    ,dimension(:) ,intent(in) :: x !! Positions of x-grid cell centres
    real(rk)    ,dimension(:) ,intent(in) :: v !! Positions of v-grid cell centres
    integer(ik)               ,intent(in) :: maxL !! Highest resolved l-harmonic
    integer(ik)               ,intent(in) :: maxM !! Highest resolved m-harmonic

    integer(ik)                                      :: i ,j

    if (assertions) then 
        call assertPure(maxL >= 0, "Negative maxL passed to grid constructor")
        call assertPure(maxM >= 0,"Negative maxM passed to grid constructor")
        call assertPure(maxM <= maxL,"Max m number passed to grid constructor must be less than or equal to the max l number")
    end if

    this%xGrid = x 
    this%vGrid = v 
    this%maxL = maxL 
    this%maxM = maxM

    allocate(this%lGrid(0))
    allocate(this%mGrid(0))
    allocate(this%imaginaryHarmonic(0))

    do i = 0,maxL 
        do j = 0,min(i,maxM)
            this%lGrid = [this%lGrid,i]
            this%mGrid = [this%mGrid,j]
            this%imaginaryHarmonic = [this%imaginaryHarmonic,.false.]

            if (j>0) then
                this%lGrid = [this%lGrid,i]
                this%mGrid = [this%mGrid,j]
                this%imaginaryHarmonic = [this%imaginaryHarmonic,.true.]
            end if
        end do
    end do

    call this%makeDefined()

end subroutine initGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getXGrid (this) result(x)
!! Getter for xGrid

    class(Grid)                           ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)             :: x

    if (assertions) call assertPure(this%isDefined(),"getXGrid called for undefined grid")

    x = this%xGrid

end function getXGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVGrid (this) result(v)
    !! Getter for vGrid
    
    class(Grid)                           ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)             :: v
    
    if (assertions) call assertPure(this%isDefined(),"getVGrid called for undefined grid")

    v = this%vGrid

end function getVGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxL (this) result(maxL)
    !! Getter for maxL

    class(Grid) ,intent(in) :: this
    integer(ik)             :: maxL

    if (assertions) call assertPure(this%isDefined(),"getMaxL called for undefined grid")

    maxL = this%maxL

end function getMaxL
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMaxM (this) result(maxM)
    !! Getter for maxM

    class(Grid) ,intent(in) :: this
    integer(ik)             :: maxM

    if (assertions) call assertPure(this%isDefined(),"getMaxM called for undefined grid")

    maxM = this%maxM

end function getMaxM
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumH (this) result(numH)
    !! Return total number of resolved harmonic on grid

    class(Grid) ,intent(in) :: this
    integer(ik)             :: numH

    if (assertions) call assertPure(this%isDefined(),"getNumH called for undefined grid")

    numH = size(this%lGrid)

end function getNumH
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumX (this) result(numX)
    !! Return number of x points on grid

    class(Grid) ,intent(in) :: this
    integer(ik)             :: numX

    if (assertions) call assertPure(this%isDefined(),"getNumX called for undefined grid")

    numX = size(this%xGrid)

end function getNumX
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumV (this) result(numV)
    !! Return number of v points on grid

    class(Grid) ,intent(in) :: this
    integer(ik)             :: numV

    if (assertions) call assertPure(this%isDefined(),"getNumV called for undefined grid")

    numV = size(this%vGrid)

end function getNumV
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLGrid (this) result(l)
    !! Getter for lGrid

    class(Grid)                           ,intent(in) :: this
    integer(ik)   ,allocatable ,dimension(:)          :: l

    if (assertions) call assertPure(this%isDefined(),"getLGrid called for undefined grid")

    l = this%lGrid

end function getLGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getMGrid (this) result(m)
    !! Getter for mGrid

    class(Grid)                           ,intent(in) :: this
    integer(ik)   ,allocatable ,dimension(:)          :: m

    if (assertions) call assertPure(this%isDefined(),"getMGrid called for undefined grid")

    m = this%mGrid

end function getMGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getHarmonicIm (this) result(im)
    !! Getter for imaginaryHarmonic

    class(Grid)                       ,intent(in) :: this
    logical   ,allocatable ,dimension(:)          :: im

    if (assertions) call assertPure(this%isDefined(),"getHarmonicIm called for undefined grid")

    im = this%imaginaryHarmonic

end function getHarmonicIm
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getH (this,l,m,im) result(h)
    !! Return index of harmonic l,m, (if im=true returns the imaginary component)

    class(Grid)  ,intent(in) :: this
    integer(ik)  ,intent(in) :: l,m
    logical      ,intent(in) :: im
    integer(ik)              :: h

    integer(ik) ,allocatable ,dimension(:) :: hVec 

    if (assertions) then 
        call assertPure(this%isDefined(),"getH called for undefined grid")
        call assertPure(this%maxL >= l, "l passed to getH out of bounds - upper (maxL)")
        call assertPure(0 <= l, "l passed to getH out of bounds - lower (0)")
        call assertPure(min(this%maxM,l) >= m, "m passed to getH out of bounds - upper (min(this%maxM,l))")
        call assertPure(0 <= m, "m passed to getH out of bounds - lower (0)")
    end if

    hVec = findIndices((this%lGrid == l) .and. (this%mGrid == m) .and. (this%imaginaryHarmonic .eqv. im))

    if (assertions) call assertPure(size(hVec) > 0,"l,m,im triple not found in grid using getH")

    h = hVec(1)

end function getH
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getL (this,ind) result(l)
    !! Return l number for given harmonic index

    class(Grid)  ,intent(in) :: this
    integer(ik)  ,intent(in) :: ind
    integer(ik)              :: l

    if (assertions) call assertPure(this%isDefined(),"getL called for undefined grid")

    l = this%lGrid(ind)

end function getL
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getM (this,ind) result(m)
    !! Return m number for given harmonic index

    class(Grid)  ,intent(in) :: this
    integer(ik)  ,intent(in) :: ind
    integer(ik)              :: m

    if (assertions) call assertPure(this%isDefined(),"getM called for undefined grid")

    m = this%mGrid(ind)

end function getM
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isImaginary (this,ind) result(im)
    !! Return true if harmonic with given index is imaginary

    class(Grid)  ,intent(in) :: this
    integer(ik)  ,intent(in) :: ind
    logical                  :: im

    if (assertions) call assertPure(this%isDefined(),"gisImaginary called for undefined grid")

    im = this%imaginaryHarmonic(ind)

end function isImaginary
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule grid_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
