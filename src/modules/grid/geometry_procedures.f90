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
submodule (geometry_class) geometry_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the geometry class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initGeometry(this,cellWidths,jLeft,jRight,periodicGrid) 
    !! Geometry initialization routine

    class(Geometry)           ,intent(inout)  :: this
    real(rk)    ,dimension(:) ,intent(in)     :: cellWidths !! Spatial cell widths
    real(rk)    ,dimension(:) ,intent(in)     :: jLeft !! Left face jacobian values
    real(rk)    ,dimension(:) ,intent(in)     :: jRight !! Right face jacobian values
    logical     ,optional     ,intent(in)     :: periodicGrid !! Set to true if the grid is periodic

    integer(ik) :: i

    logical :: periodic

    if (assertions) then
        call assertPure((size(cellWidths) == size(jLeft)) .and. (size(jLeft) == size(jRight)),&
        "Cell widths and face jacobians passed to geometry initializations must all be of same size")
    end if

    this%cellWidths = cellWidths 
    this%jacobianLeft = jLeft 
    this%jacobianRight = jRight 
    this%jacobianCentre = (jLeft + jRight)/2 

    periodic = .false. 
    if (present(periodicGrid)) periodic = periodicGrid

   
    allocate(this%linInterp(0:size(cellWidths)))
    
    this%linInterp = 0

    do i = 1,size(cellWidths)-1
        this%linInterp(i) = cellWidths(i)/(cellWidths(i)+cellWidths(i+1)) 
    end do 

    this%linInterp(size(cellWidths)) = real(1.0d0,kind=rk) !Handle non-periodic interpolation without explicit checks

    if (periodic) then 
        this%linInterp(size(cellWidths)) = cellWidths(size(cellWidths))/(cellWidths(size(cellWidths))+cellWidths(1)) 
        this%linInterp(0) = this%linInterp(size(cellWidths))
    end if

    this%periodicGrid = periodic 
    
    call this%makeDefined()
    
end subroutine initGeometry
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCellWidths (this,dualGrid,extendedBoundaryCells) result(dx)
    !! Getter for cellWidths, if dualGrid is true returns dual grid values based off of cellWidths. If extendedBoundaryCells is
    !! true will extend dual grid cells at boundaries to the boundaries themselves if not periodic (defaults to true)

    class(Geometry)                       ,intent(in) :: this
    logical ,optional                     ,intent(in) :: dualGrid
    logical ,optional                     ,intent(in) :: extendedBoundaryCells
    real(rk)   ,allocatable ,dimension(:)             :: dx

    logical :: dGrid,exBoundaryCells

    if (assertions) call assertPure(this%isDefined(),"getCellWidths called from undefined geometry object")

    dGrid = .false.
    if (present(dualGrid)) dGrid = dualGrid

    exBoundaryCells = .true.
    if (present(extendedBoundaryCells)) exBoundaryCells = extendedBoundaryCells

    dx = this%cellWidths

    if (dGrid) then 
        dx(1:size(dx)-1) = (this%cellWidths(1:size(dx)-1) + this%cellWidths(2:size(dx)))/2
        if (this%periodicGrid) then 
            dx(size(dx)) = (this%cellWidths(1) + this%cellWidths(size(dx)))/2
        else if (exBoundaryCells) then
            dx(1) = this%cellWidths(1) + this%cellWidths(2)/2
            dx(size(dx)-1) = this%cellWidths(size(dx)) + this%cellWidths(size(dx)-1)/2
            dx(size(dx)) = 10*epsilon(dx) ! Set to this to avoid divide by 0 errors
        end if
    end if
end function getCellWidths
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getJacobianLeft (this,dualGrid) result(jLeft)
    !! Getter for jacobianLeft, if dualGrid is true returns dual grid values 

    class(Geometry)                       ,intent(in) :: this
    logical ,optional                     ,intent(in) :: dualGrid
    real(rk)   ,allocatable ,dimension(:)             :: jLeft
    logical :: dGrid

    if (assertions) call assertPure(this%isDefined(),"getJacobianLeft called from undefined geometry object")

    dGrid = .false.
    if (present(dualGrid)) dGrid = dualGrid

    jLeft = this%jacobianLeft

    if (dGrid) then 
        jLeft = this%jacobianCentre
        if (.not. this%periodicGrid) then 
            jLeft(1) = this%jacobianLeft(1)
            jLeft(size(jLeft)) = this%jacobianRight(size(jLeft))
        end if
    end if

end function getJacobianLeft
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getJacobianRight (this,dualGrid,extendedBoundaryCells) result(jRight)
    !! Getter for jacobianRight, if dualGrid is true returns dual grid values 

    class(Geometry)                       ,intent(in) :: this
    logical ,optional                     ,intent(in) :: dualGrid
    logical ,optional                     ,intent(in) :: extendedBoundaryCells
    real(rk)   ,allocatable ,dimension(:)             :: jRight
    logical :: dGrid, exBoundaryCells

    if (assertions) call assertPure(this%isDefined(),"getJacobianRight called from undefined geometry object")

    dGrid = .false.
    if (present(dualGrid)) dGrid = dualGrid

    exBoundaryCells = .true.
    if (present(extendedBoundaryCells)) exBoundaryCells = extendedBoundaryCells

    jRight = this%jacobianRight

    if (dGrid) then 

        jRight(1:size(jRight)-1) = this%jacobianCentre(2:size(jRight))

        if (this%periodicGrid) then 
            jRight(size(jRight)) = this%jacobianCentre(1)
        else if (exBoundaryCells) then
            jRight(size(jRight)-1) = this%jacobianRight(size(jRight))
            jRight(size(jRight)) = this%jacobianRight(size(jRight))
        end if
    end if

end function getJacobianRight
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getJacobianCentre (this,dualGrid,extendedBoundaryCells) result(jCentre)
    !! Getter for jacobianCentre, if dualGrid is true returns dual grid values 

    class(Geometry)                       ,intent(in) :: this
    logical ,optional                     ,intent(in) :: dualGrid
    logical ,optional                     ,intent(in) :: extendedBoundaryCells
    real(rk)   ,allocatable ,dimension(:)             :: jCentre
    logical :: dGrid ,exBoundaryCells

    if (assertions) call assertPure(this%isDefined(),"getJacobianCentre called from undefined geometry object")

    dGrid = .false.
    if (present(dualGrid)) dGrid = dualGrid

    exBoundaryCells = .true.
    if (present(extendedBoundaryCells)) exBoundaryCells = extendedBoundaryCells

    jCentre = this%jacobianCentre

    if (dGrid) then 

        !The straightforward thing here is to take jCentre on the dual grid to be jRight on the original grid. However, since 
        !jCentre is used to get the (dual) cell volume as jCentre*dx, this would not be accurate for the dual grid, since the 
        !dual cell cross-section is in general two connected trapezoids. Hence the following is done to make sure that the dual cell volume is correct.
        !NOTE: In the non-periodic case, the first and second to last dual cells extend to the actual system boundaries so that the dual grid 
        !covers the original grid 

        jCentre(1:size(jCentre)-1) =  this%jacobianRight(1:size(jCentre)-1)/2 &
                                    + (this%jacobianCentre(1:size(jCentre)-1)*this%cellWidths(1:size(jCentre)-1) &
                                    + this%jacobianCentre(2:size(jCentre))*this%cellWidths(2:size(jCentre)))&
                                    / (2*(this%cellWidths(1:size(jCentre)-1)+this%cellWidths(2:size(jCentre))))
        

        if (this%periodicGrid) then 
            jCentre(size(jCentre)) = this%jacobianRight(size(jCentre))/2 &
                                    + (this%jacobianCentre(size(jCentre))*this%cellWidths(size(jCentre)) &
                                    + this%jacobianCentre(1)*this%cellWidths(1))&
                                    / (2*(this%cellWidths(size(jCentre))+this%cellWidths(1)))
        else if (exBoundaryCells) then

            jCentre(1) = this%cellWidths(1)*this%jacobianCentre(1) &
                        + (this%jacobianRight(1)+this%jacobianCentre(2))*this%cellWidths(2)/4
            jCentre(1) = jCentre(1)/(this%cellWidths(1)+this%cellWidths(2)/2)

            jCentre(size(jCentre)-1) = this%cellWidths(size(jCentre))*this%jacobianCentre(size(jCentre)) &
                                        + (this%jacobianRight(size(jCentre)-1)+this%jacobianCentre(size(jCentre)-1))&
                                        *this%cellWidths(size(jCentre)-1)/4

            jCentre(size(jCentre)-1) = jCentre(size(jCentre)-1)/(this%cellWidths(size(jCentre))&
                                      +this%cellWidths(size(jCentre)-1)/2)

            jCentre(size(jCentre)) = 10*epsilon(jCentre) ! Set to this to avoid divide by 0 errors

        end if
    end if

end function getJacobianCentre
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getLinInterp (this,dualGrid) result(linInterp)
    !! Getter for linInterp, if dualGrid is true returns dual grid values (0.5 everywhere)

    class(Geometry)                       ,intent(in) :: this
    logical ,optional                     ,intent(in) :: dualGrid
    real(rk)   ,allocatable ,dimension(:)             :: linInterp
    logical :: dGrid

    if (assertions) call assertPure(this%isDefined(),"getLinInterp called from undefined geometry object")

    dGrid = .false.
    if (present(dualGrid)) dGrid = dualGrid

    linInterp = this%linInterp

    if (dGrid) linInterp = real(0.5d0,kind=rk)

end function getLinInterp
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function isPeriodic (this) result(periodic)
    !! Getter for periodicGrid

    class(Geometry) ,intent(in) :: this
    logical                     :: periodic

    if (assertions) call assertPure(this%isDefined(),"isPeriodic called from undefined geometry object")

    periodic = this%periodicGrid

end function isPeriodic
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule geometry_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
