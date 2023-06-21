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
submodule (v_space_class) v_space_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the VSpace class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initVSpace(this,gridObj) 
    !! VSpace initialization routine
    
    class(VSpace)           ,intent(inout)  :: this
    type(Grid)              ,intent(in)     :: gridObj !! Grid object used to initialize VSpace

    integer(ik)                             :: i 

    if (assertions) call assertPure(gridObj%isDefined(),"Attempted to initialize vSpace object using undefined grid")

    this%vGrid = gridObj%getVGrid()
    this%numH = gridObj%getNumH()
    this%numV = size(this%vGrid)
    allocate(this%vWidths(this%numV))
    this%vWidths(1) = 2 * this%vGrid(1)
    do i = 2, this%numV 
        this%vWidths(i) = 2 * (this%vGrid(i) - sum(this%vWidths(1:i-1)))
    end do

    allocate(this%linInterp(this%numV))

    this%linInterp(1:this%numV-1) = this%vWidths(1:this%numV-1)/(this%vWidths(1:this%numV-1)+this%vWidths(2:this%numV))
    this%linInterp(this%numV) = real(1.0d0,kind=rk)
    call this%makeDefined()

end subroutine initVSpace
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVGrid (this) result(v)
    !! Getter for vGrid

    class(VSpace)                         ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)             :: v

    if (assertions) call assertPure(this%isDefined(),"Velocity grid requested from undefined vSpace object")

    v = this%vGrid

end function getVGrid
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumH (this) result(numH)
    !! Return total number of resolved harmonic on grid

    class(VSpace) ,intent(in) :: this
    integer(ik)               :: numH

    if (assertions) call assertPure(this%isDefined(),"numH requested from undefined vSpace object")

    numH = this%numH 

end function getNumH
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNumV (this) result(numV)
    !! Return number of v points on grid

    class(VSpace) ,intent(in) :: this
    integer(ik)               :: numV

    if (assertions) call assertPure(this%isDefined(),"numV requested from undefined vSpace object")

    numV = this%numV

end function getNumV
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVCellWidths (this) result(dv)
    !! Getter for vWidths

    class(VSpace)                       ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)           :: dv

    if (assertions) call assertPure(this%isDefined(),"Velocity grid cell widths requested from undefined vSpace object")

    dv = this%vWidths

end function getVCellWidths
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getVLinInterp (this) result(linInterp)
    !! Getter for linInterp

    class(VSpace)                       ,intent(in) :: this
    real(rk)   ,allocatable ,dimension(:)           :: linInterp

    if (assertions) call assertPure(this%isDefined(),&
    "Velocity linear interpolation coefficients requested from undefined vSpace object")

    linInterp = this%linInterp

end function getVLinInterp
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNearestPoints (this,v) result(pair)
    !! Return two nearest points for given velocity value v. If the first point is 0, the velocity is below the lowest velocity in the 
    !! grid, and if the second point is 0 the velocity is above the greatest v in the grid.

    class(VSpace)             ,intent(in) :: this
    real(rk)                  ,intent(in) :: v
    integer(ik)   ,dimension(2)           :: pair

    if (assertions) call assertPure(this%isDefined(),"getNearestPoints called from undefined vSpace object")

    pair = findNearestPointsInArray(this%vGrid,v)

end function getNearestPoints
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getContainingVCell (this,v) result(ind)
    !! Return index of cell which containes the given velocity value v. If the returned index is 0, the point is outside of the grid. 

    class(VSpace)             ,intent(in) :: this
    real(rk)                  ,intent(in) :: v
    integer(ik)                           :: ind

    if (assertions) call assertPure(this%isDefined(),"getContainingVCell called from undefined vSpace object")

    ind = 0 

    if ((v >= 0) .and. (v < this%vGrid(this%numV) + this%vWidths(this%numV)/2)) &
    ind = findloc((v >= this%vGrid - this%vWidths/2) .and. (v < this%vGrid + this%vWidths/2),.true.,1)

end function getContainingVCell
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function calculateMoment (this,f,h,mOrder,g,gDependsOnX) result(res)
    !! Calculate moment of the h-th harmonic of passed local distribution function f, optionally multiplied by g(v). If gDependsOnX 
    !! is .true. assumes that g is given as a strided array with stride numV, and it is used to allocate the result - 
    !! this allows g to not have a halo while f does. The moment is of order mOrder. 
    
    class(VSpace)                       ,intent(in) :: this
    real(rk)    ,dimension(:)           ,intent(in) :: f
    integer(ik)                         ,intent(in) :: h, mOrder
    real(rk)    ,optional ,dimension(:) ,intent(in) :: g
    logical     ,optional               ,intent(in) :: gDependsOnX

    real(rk)   ,allocatable ,dimension(:)           :: res

    integer(ik)                                     :: i ,fOffset ,gOffset ,fghaloDiff
    logical                                         :: gX

    if (assertions) call assertPure(this%isDefined(),&
    "Attempted to calculate moment of distribution function using undefined vSpace")

    if (present(g)) then 

        gX = .false.
        if (present(gDependsOnX)) gX = gDependsOnX

        if (gX) then  
            allocate(res(size(g)/(this%numV))) !If g depends on x use g's x dependence for the allocation of the result
            do i = 1, size(res) 
                fghaloDiff = (size(f)/(this%numV*this%numH) - size(g)/(this%numV))/2
                gOffset = (i-1)*this%numV
                fOffset = (i-1+fghaloDiff)*this%numH*this%numV + (h-1)*this%numV 
                res(i) = simpleMoment(mOrder,f(fOffset+1:fOffset+this%numV),this%vGrid,this%vWidths,g(gOffset+1:gOffset+this%numV))
            end do
        else
            allocate(res(size(f)/(this%numV*this%numH)))
            do i = 1, size(res) 
                fOffset = (i-1)*this%numH*this%numV + (h-1)*this%numV 
                res(i) = simpleMoment(mOrder,f(fOffset+1:fOffset+this%numV),this%vGrid,this%vWidths,g)
            end do

        end if

    else
        allocate(res(size(f)/(this%numV*this%numH)))
        do i = 1, size(res) 
            fOffset = (i-1)*this%numH*this%numV + (h-1)*this%numV 
            res(i) = simpleMoment(mOrder,f(fOffset+1:fOffset+this%numV),this%vGrid,this%vWidths)
        end do

    end if
end function calculateMoment
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getShkarofskyIMat (this,index) result(res)
    !! Return sparse row data format for the lower triangular Shkarofsky I_index integral.

    class(VSpace)                       ,intent(in) :: this
    integer(ik)                         ,intent(in) :: index
    type(SparseRowData)                             :: res

    integer(ik) :: i ,j
    real(rk) ,allocatable ,dimension(:) :: vdvArr
    if (assertions) call assertPure(this%isDefined(),"getShkarofskyIMat called from undefined velocity space")

    vdvArr = 4*pi*this%vGrid**(2+index)*this%vWidths

    call res%init()
    do i = 1,this%getNumV()
        call res%addRow(i,[(j,j=1,i)],vdvArr(1:i)/this%vGrid(i)**index)
        res%values(i)%entry(i) = res%values(i)%entry(i)/2
    end do
end function getShkarofskyIMat
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getShkarofskyJMat (this,index) result(res)
    !! Return sparse row data format for the upper triangular Shkarofsky J_index integral.

    class(VSpace)                       ,intent(in) :: this
    integer(ik)                         ,intent(in) :: index
    type(SparseRowData)                             :: res

    integer(ik) :: i ,j
    real(rk) ,allocatable ,dimension(:) :: vdvArr
    if (assertions) call assertPure(this%isDefined(),"getShkarofskyJMat called from undefined velocity space")

    vdvArr = 4*pi*this%vGrid**(2+index)*this%vWidths

    call res%init()
    do i = 1,this%getNumV()
        call res%addRow(i,[(j,j=i,this%getNumV())],vdvArr(i:this%getNumV())/this%vGrid(i)**index)
        res%values(i)%entry(1) = res%values(i)%entry(1)/2
    end do
end function getShkarofskyJMat
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule v_space_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
