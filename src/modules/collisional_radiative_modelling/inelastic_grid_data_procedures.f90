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
submodule (inelastic_grid_data_class) inelastic_grid_data_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the inelastic grid data class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initInelData(this,space,fixedEnergies,interpolationGrid) 
    !! Inelastic grid data initialization routine

    class(InelasticGridData)           ,intent(inout)  :: this
    type(VSpace)                       ,intent(in)     :: space !! Velocity space on which weights should be calculated
    real(rk) ,optional ,dimension(:)   ,intent(in)     :: fixedEnergies !! Energy values for fixed energy mappings
    real(rk) ,optional ,dimension(:)   ,intent(in)     :: interpolationGrid !! Energy values for mappings on the interpolation grid

    integer(ik) :: i

    if (assertions) call assertPure(space%isDefined(),"Undefined vSpace object passed to inelastic grid data constructor")

    this%numV = size(space%getVGrid())
    if (present(fixedEnergies)) then 
        this%fixedE = fixedEnergies 
        allocate(this%fixedW(size(fixedEnergies)))
        allocate(this%fixedEmissionVecs(size(fixedEnergies)))
        do i = 1,size(fixedEnergies)
            this%fixedW(i) = inelWeights(space,fixedEnergies(i))
            allocate(this%fixedEmissionVecs(i)%entry(this%numV))
            call fillEmissionVector(this%fixedW(i),this%fixedEmissionVecs(i)%entry)
        end do
    else 
        allocate(this%fixedE(0))
        allocate(this%fixedW(0))
        allocate(this%fixedEmissionVecs(0))
    end if

    if (present(interpolationGrid)) then 
        if (assertions .and. size(interpolationGrid) > 0) then
            call assertPure(size(interpolationGrid) > 1,&
                            "interpolationGrid in InelasticGridData must have length > 1")
            call assertPure(all(interpolationGrid(2:size(interpolationGrid)) > interpolationGrid(1:size(interpolationGrid)-1)),&
                            "interpolationGrid in InelasticGridData must be strictly monotonically increasing")
            call assertPure(all(interpolationGrid /= 0),&
                            "interpolationGrid in InelasticGridData cannot contain zero energy")
        end if
        this%eInterpGrid = interpolationGrid 
        allocate(this%interpW(size(interpolationGrid)))
        do i = 1,size(interpolationGrid)
            this%interpW(i) = inelWeights(space,interpolationGrid(i))
        end do
    else 
        allocate(this%eInterpGrid(0))
        allocate(this%interpW(0))
    end if

    call this%makeDefined()

end subroutine initInelData
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFixedW(this,ind) result(wMat)
    !! Return fixed mapping matrix with given index

    class(InelasticGridData)       ,intent(in) :: this 
    integer(ik)                    ,intent(in) :: ind
    type(SparseRowData)                        :: wMat

    if (assertions) then 
        call assertPure(this%isDefined(),"getFixedW called on undefined inelastic grid data object")
        call assertPure(ind <= size(this%fixedW),"Index passed to getFixedW out of range - upper")
        call assertPure(ind > 0,"Index passed to getFixedW out of range - lower")

    end if

    wMat = this%fixedW(ind)

end function getFixedW
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getFixedEmissionVector(this,ind) result(emit)
    !! Return fixed emission vector for mapping matrix with given index

    class(InelasticGridData)       ,intent(in) :: this 
    integer(ik)                    ,intent(in) :: ind
    real(rk) ,allocatable ,dimension(:)        :: emit

    if (assertions) then 
        call assertPure(this%isDefined(),"getFixedEmissionVector called on undefined inelastic grid data object")
        call assertPure(ind <= size(this%fixedW),"Index passed to getFixedEmissionVector out of range - upper")
        call assertPure(ind > 0,"Index passed to getFixedEmissionVector out of range - lower")

    end if

    emit = this%fixedEmissionVecs(ind)%entry

end function getFixedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine interpolateW(this,E,wRes) 
    !!  Interpolate mapping matrices for given energy and store in passed triangular matrix. Assumes upper triangular structure for wRes
    !! if E is positive and lower if it's negative

    class(InelasticGridData)     ,intent(in)    :: this
    real(rk)                     ,intent(in)    :: E !! Transition energy to interpolate for
    type(SparseRowData)          ,intent(inout) :: wRes !! Lower/upper triangular matrix to store the interpolated weights

    type(Interpolation1D) :: interpObj

    integer(ik) ,allocatable ,dimension(:) :: ind 
    real(rk) ,allocatable ,dimension(:) :: weight, tempEmit

    integer(ik) :: i ,j ,rowInd

    if (assertions) then
        call assertPure(this%isDefined(),&
                        "interpolateW called for undefined inelastic grid data object")
        call assertPure(size(this%eInterpGrid) > 1,&
                        "interpolateW requires energy interpolation grid of length > 1")
    end if

    call interpObj%init(this%eInterpGrid,[E])
    ind = interpObj%getFirstDataIndices()
    weight = interpObj%getInterpWeights()

    if (assertions) then
        call assertPure(weight(1) >= 0,"Transition energy E outside energy interpolation grid in interpolateW")
        call assertPure(this%eInterpGrid(ind(1)) * this%eInterpGrid(ind(1)+1) > 0,&
                        "Transition energy E falls in sign-changing interval of energy interpolation grid in interpolateW")
    end if

    !Handle upper triangular wRes with positive E 
    if (E > 0) then 
        do i = 1, size(wRes%rowIndex)
            wRes%values(i)%entry = 0
        end do

        do i = 1, size(this%interpW(ind(1))%rowIndex)
            rowInd = this%interpW(ind(1))%rowIndex(i)
            do j = 1,size(this%interpW(ind(1))%columnVector(i)%entry)
                wRes%values(rowInd)%entry(this%interpW(ind(1))%columnVector(i)%entry(j)-(rowInd-1)) = &
                this%interpW(ind(1))%values(i)%entry(j) * (real(1,kind=rk) - weight(1))
            end do
        end do
 
        do i = 1, size(this%interpW(ind(1)+1)%rowIndex)
            rowInd = this%interpW(ind(1)+1)%rowIndex(i)
            do j = 1,size(this%interpW(ind(1)+1)%columnVector(i)%entry)
                wRes%values(rowInd)%entry(this%interpW(ind(1)+1)%columnVector(i)%entry(j)-(rowInd-1)) = &
                wRes%values(rowInd)%entry(this%interpW(ind(1)+1)%columnVector(i)%entry(j)-(rowInd-1)) + &
                this%interpW(ind(1)+1)%values(i)%entry(j) * weight(1)
            end do
        end do
    else !Handle lower triangular wRes with negative E
        do i = 1, size(wRes%rowIndex)
            wRes%values(i)%entry = 0
        end do

        do i = 1, size(this%interpW(ind(1))%rowIndex)
            rowInd = this%interpW(ind(1))%rowIndex(i)
            do j = 1,size(this%interpW(ind(1))%columnVector(i)%entry)
                wRes%values(rowInd)%entry(this%interpW(ind(1))%columnVector(i)%entry(j)) = &
                this%interpW(ind(1))%values(i)%entry(j) * (real(1,kind=rk) - weight(1))
            end do
        end do
        do i = 1, size(this%interpW(ind(1)+1)%rowIndex)
            rowInd = this%interpW(ind(1)+1)%rowIndex(i)
            do j = 1,size(this%interpW(ind(1)+1)%columnVector(i)%entry)
                wRes%values(rowInd)%entry(this%interpW(ind(1)+1)%columnVector(i)%entry(j)) = &
                wRes%values(rowInd)%entry(this%interpW(ind(1)+1)%columnVector(i)%entry(j)) + &
                this%interpW(ind(1)+1)%values(i)%entry(j) * weight(1)
            end do
        end do
    end if 

end subroutine interpolateW
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpolatedEmissionVector(this,E) result(emit)
    !! Return interpolated emission vector for given input energy E

    class(InelasticGridData)       ,intent(in) :: this 
    real(rk)                       ,intent(in) :: E
    real(rk) ,allocatable ,dimension(:)        :: emit

    type(Interpolation1D) :: interpObj

    integer(ik) ,allocatable ,dimension(:) :: ind 
    real(rk) ,allocatable ,dimension(:) :: weight, tempEmit
 
    if (assertions) then
        call assertPure(this%isDefined(),&
                        "getInterpolatedEmissionVector called for undefined inelastic grid data object")
        call assertPure(size(this%eInterpGrid) > 1,&
                        "getInterpolatedEmissionVector requires energy interpolation grid of length > 1")
    end if
    allocate(emit(this%numV))
    allocate(tempEmit(this%numV))
    call interpObj%init(this%eInterpGrid,[E])
    ind = interpObj%getFirstDataIndices()
    weight = interpObj%getInterpWeights()
    if (assertions) then
        call assertPure(weight(1) >= 0,&
                        "Transition energy E outside energy interpolation grid in getInterpolatedEmissionVector")
        call assertPure(this%eInterpGrid(ind(1)) * this%eInterpGrid(ind(1)+1) > 0,&
                        "Transition energy E falls in sign-changing interval of energy interpolation grid in&
                        & getInterpolatedEmissionVector")
    end if
    call fillEmissionVector(this%interpW(ind(1)),emit)
    call fillEmissionVector(this%interpW(ind(1)+1),tempEmit)

    emit = emit * (real(1,kind=rk) - weight(1)) + weight(1) * tempEmit

end function getInterpolatedEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule inelastic_grid_data_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
