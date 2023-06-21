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
submodule (modelbound_lbc_data_class) modelbound_lbc_data_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the modelbound logical boundary condition data class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMBLBC(this,vSpaceObj,distExtDerivation,distExtReqVarIndices,&
                                ionCurrentVarIndex,isActive,totalCurrentVarIndex,bisTol,isLeftBoundary) 
    !! Varlike modelbound data initialization routine

    class(ModelboundLBCData)        ,intent(inout)  :: this
    type(VSpace)                    ,intent(in)     :: vSpaceObj !! Velocity space object used to get grid data
    class(MatDerivation)            ,intent(in)     :: distExtDerivation !! Extrapolated distribution derivation object
    integer(ik) ,dimension(:)       ,intent(in)     :: distExtReqVarIndices !! Required variable indices for the distribution derivation 
    integer(ik)                     ,intent(in)     :: ionCurrentVarIndex !! Scalar variable index representing ion current at boundary
    logical                         ,intent(in)     :: isActive !! True if the modelbound data should be updated (use to avoid updates on processors with no boundary)
    integer(ik) ,optional           ,intent(in)     :: totalCurrentVarIndex !! Scalar variable index representing total current through boundary. Defaults to 0 current.
    real(rk)    ,optional           ,intent(in)     :: bisTol !! Bisection tolerance for false position method used to calculate cut-off velocity. Defaults to 1e-12
    logical     ,optional           ,intent(in)     :: isLeftBoundary !! True if boundary this data refers to is the left boundary. Defaults to false. 

    if (assertions) then 

            call assertPure(vSpaceObj%isDefined(),"Undefined VSpace object passed to initMBLBC")
            call assertPure(distExtDerivation%isDefined(),"Undefined distribution derivation passed to init MBLBC")

    end if

    this%isActive = isActive

    if (isActive) then 
        allocate(this%fextDeriv,source=distExtDerivation)
        this%fextReqIndices = distExtReqVarIndices
        this%ionCurrentVarIndex = ionCurrentVarIndex
        if (present(totalCurrentVarIndex)) this%totalCurrentVarIndex = totalCurrentVarIndex
        
        this%bisTol = real(1d-12,kind=rk)
        if (present(bisTol)) this%bisTol = bisTol

        this%isLeftBoundary = .false.
        if (present(isLeftBoundary)) this%isLeftBoundary = isLeftBoundary

        this%vGridCopy = vSpaceObj%getVGrid()
        this%dvCopy = vSpaceObj%getVCellWidths()

        allocate(this%vBoundaries(1:size(this%vGridCopy)+1))
        this%vBoundaries(1:size(this%vGridCopy)) = this%vGridCopy - this%dvCopy/2
        this%vBoundaries(size(this%vGridCopy)+1) = this%vGridCopy(size(this%vGridCopy)) + this%dvCopy(size(this%vGridCopy))/2

        this%maxL = vSpaceObj%getNumH() - 1 !Assuming no m>0 harmonics 

        allocate(this%fext(this%maxL+1,size(this%dvCopy)))
        this%fext = 0 

        allocate(this%fixedPLL(0:this%maxL,0:this%maxL))
        this%fixedPLL = fixedPll(this%maxL)
    end if

    call this%makeDefined()
end subroutine initMBLBC
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateDataLBC(this,hostModel,inputVars,updatePriority) 
    !!  Update modelbound data based on input variable container

    class(ModelboundLBCData)              ,intent(inout) :: this 
    class(ModelSurrogate)                 ,intent(in)    :: hostModel !! Host model - unused
    class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate modelbound data
    integer(ik) ,optional                 ,intent(in)    :: updatePriority !! Priority for this update call (determines which variables are updated) - unused here

    real(rk) :: ionCurrent, totCurrent ,boundaryWe ,boundaryDens ,boundaryFlux ,boundaryHeatflux

    real(rk) ,allocatable ,dimension(:) :: f0, f1 

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to updated undefined logical boundary modelbound data")
        call assert(inputVars%isDefined(),"Attempted to updated logical boundary modelbound data using undefined inputVars")
    end if

    if (this%isActive) then 
        this%fext = this%fextDeriv%calculate(inputVars%variables,this%fextReqIndices)

        ionCurrent = inputVars%variables(this%ionCurrentVarIndex)%entry(1)
        totCurrent = 0 

        if (allocated(this%totalCurrentVarIndex)) &
        totCurrent = inputVars%variables(this%totalCurrentVarIndex)%entry(1)

        call this%calculateCutOff(ionCurrent,totCurrent)
        call this%calculatePll()

        f0 = this%boundaryHarmonic(0)
        f1 = this%boundaryHarmonic(1)
        boundaryDens = this%interpMom(0,f0)
        boundaryWe = this%interpMom(2,f0)
        boundaryFlux = this%interpMom(1,f1)/real(3,kind=rk)
        boundaryHeatflux = this%interpMom(3,f1)/real(3,kind=rk)
        this%shTemp = 2 * boundaryWe/(3*boundaryDens) - elMass * boundaryFlux**2/boundaryDens**2

        this%gamma = boundaryHeatflux/(boundaryFlux*this%shTemp)
        this%potential = this%coVel**2/this%shTemp
    end if

end subroutine updateDataLBC
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine copyDataLBC(this,name,container) 
    !! Copy named modelbound data to passed container.

    class(ModelboundLBCData)              ,intent(in)    :: this 
    character(*)                          ,intent(in)    :: name !! Name of data
    real(rk) ,allocatable ,dimension(..)  ,intent(inout) :: container !! Container to copy into 

    if (assertions) call assert(this%isDefined(),"Attempted to copy data from undefined logical boundary modelbound data object")

    select rank (container)
    rank (1)
        if (allocated(container)) then 
            if (ubound(container,1) /= 1 .or. &
               lbound(container,1) /= 1) then 
                deallocate(container)
                allocate(container(1))
            end if
        else
            allocate(container(1))
        end if

        select case (name)
        case ("gamma")
            container = this%gamma
        case ("potential")
            container = this%potential
        case ("coVel")
            container = this%coVel
        case ("shTemp")
            container = this%shTemp
        case default
            error stop "unregistered name requested from ModelboundLBCData object"
        end select
    rank default 
        error stop "container passed to copyDataLBC is not rank 1"
    end select

end subroutine copyDataLBC
!-----------------------------------------------------------------------------------------------------------------------------------
module function getDataDimLBC(this,name) result(dim)
    !! Get data dimensionality - will return 0 for scalars

    class(ModelboundLBCData)              ,intent(in)    :: this 
    character(*)                          ,intent(in)    :: name !! Name of data
    integer(ik)                                          :: dim

    if (assertions) &
    call assert(this%isDefined(),"Attempted to get data dimensionality from undefined logical boundary modelbound data object")

    select case (name)
    case ("gamma")
        dim = 0
    case ("potential")
        dim = 0
    case ("coVel")
        dim = 0
    case ("shTemp")
        dim = 0
    case default
        error stop "unregistered name requested from ModelboundLBCData object"
    end select

end function getDataDimLBC
!-----------------------------------------------------------------------------------------------------------------------------------
function fixedPll(maxL) result(res)
    !! Calculate fixed component of Pll tensor 

    integer(ik)           ,intent(in)     :: maxL
    real(rk) ,allocatable ,dimension(:,:) :: legPl 
    real(rk) ,allocatable ,dimension(:,:) :: res 

    integer(ik) :: i ,j
    
    allocate(legPl(1,0:maxL))
    legPl = allPl([real(0,kind=rk)],maxL)

    allocate(res(0:maxL,0:maxL))
     
    res = 0 
    do i = 0,maxL 

        if (mod(i,2)==1) then

            do j = 0,maxL 

                if (mod(j,2)==1 .and. i == j) res(j,i) = real(1,kind=rk) 

                if (mod(j,2)==0) then 

                    if (j==0) then 
                        res(j,i) = - i * legPl(1,i-1)
                    else
                        res(j,i) = j * legPl(1,j-1)*legPl(1,i) - i * legPl(1,j)*legPl(1,i-1)
                    end if

                    res(j,i) = res(j,i) * real(2*j+1,kind=rk)/(real((j-i)*(i+j+1),kind=rk))
                end if 
            end do 

        end if

    end do

end function fixedPll
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calculatePll(this,rowsToUpdate)
    !! Calculate full Pll tensor based on modelbound data 

    class(ModelboundLBCData)              ,intent(inout)    :: this 
    integer(ik) ,optional ,dimension(:)   ,intent(in)       :: rowsToUpdate

    integer(ik) :: maxL ,numV ,i ,j ,k

    real(rk) ,allocatable ,dimension(:) :: xmin

    real(rk) ,allocatable ,dimension(:,:) :: legPl

    integer(ik) ,allocatable ,dimension(:) :: lRows
    if (assertions) call assert(this%isDefined(),"calculatePll called on undefined modelbound LBC data")

    maxL = size(this%fixedPLL,1)-1
    numV = size(this%vGridCopy)

    if (.not. allocated(this%bufferPLL)) then 
        allocate(this%bufferPLL(numV,0:maxL,0:maxL))
    end if
    xmin = - this%coVel/this%vGridCopy
    xmin(this%coCell) = - this%coVel/this%interpCoord(2)

    xmin = max(xmin,real(-1,kind=rk))
    if (present(rowsToUpdate)) then 
        lRows = rowsToUpdate
    else
        lRows = [(i,i=0,maxL)]
    end if
    allocate(legPl(numV,0:maxL))
    legPl = allPl(xmin,maxL)

    do i = 0,maxL
        do j = 1,size(lRows) 
            this%bufferPLL(:,lRows(j),i) = 0

            if (i==lRows(j)) then 
                this%bufferPLL(1:this%coCell-1,lRows(j),i) = real(1,kind=rk)
                this%bufferPLL(this%coCell:numV,lRows(j),i) = real(0.5d0,kind=rk) - xmin(this%coCell:numV)&
                                                                                   *legPl(this%coCell:numV,lRows(j))**2/2

                do k = 1,i-1
                    this%bufferPLL(this%coCell:numV,lRows(j),i) = this%bufferPLL(this%coCell:numV,lRows(j),i) &
                                                - legPl(this%coCell:numV,k)*(xmin(this%coCell:numV)*legPl(this%coCell:numV,k) &
                                                                            - legPl(this%coCell:numV,k+1))
                                                                            
                end do

                this%bufferPLL(:,lRows(j),i) = this%bufferPLL(:,lRows(j),i) * (-1) ** i
            else

                this%bufferPLL(this%coCell:numV,lRows(j),i) = - xmin(this%coCell:numV) * legPl(this%coCell:numV,i) &
                                                      * legPl(this%coCell:numV,lRows(j)) 

                if (i > 0) this%bufferPLL(this%coCell:numV,lRows(j),i) = this%bufferPLL(this%coCell:numV,lRows(j),i) &
                                                                - i*legPl(this%coCell:numV,lRows(j))*legPl(this%coCell:numV,i-1)&
                                                                /real(lRows(j)-i,kind=rk)

                if (lRows(j) > 0) this%bufferPLL(this%coCell:numV,lRows(j),i) = this%bufferPLL(this%coCell:numV,lRows(j),i) &
                                                                 + lRows(j)*legPl(this%coCell:numV,lRows(j)-1)&
                                                                 *legPl(this%coCell:numV,i) &
                                                                 /real(lRows(j)-i,kind=rk)

                this%bufferPLL(this%coCell:numV,lRows(j),i) = this%bufferPLL(this%coCell:numV,lRows(j),i) &
                                                       * (-1)**i * real(2*lRows(j)+1,kind=rk)/real(2*(i+lRows(j)+1),kind=rk)
            end if

            this%bufferPLL(:,lRows(j),i) = this%bufferPLL(:,lRows(j),i) + this%fixedPLL(lRows(j),i)
            if (this%isLeftBoundary) this%bufferPLL(:,lRows(j),i) = this%bufferPLL(:,lRows(j),i) * (-1)**(i+lRows(j))

        end do 

    end do

end subroutine calculatePll
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calculateInterps(this)
    !! Calculate interpolation quantities based on modelbound data 

    class(ModelboundLBCData)              ,intent(inout)    :: this 

    if (assertions) call assert(this%isDefined(),"calculateInterps called on undefined modelbound LBC data")

    this%interpCoord(1) = (this%vBoundaries(this%coCell-1)+this%coVel)/2 !Left value
    this%interpCoord(2) = (this%vBoundaries(this%coCell+1)+this%coVel)/2 !Right value

    this%vInterp(1) = real(1,kind=rk) - (this%interpCoord(1)-this%vGridCopy(this%coCell-1))/this%dvCopy(this%coCell-1)
    this%vInterp(2) = real(1,kind=rk) - (this%interpCoord(2)-this%vGridCopy(this%coCell))/this%dvCopy(this%coCell)

    this%interpdv(1) = this%coVel - this%vBoundaries(this%coCell-1)
    this%interpdv(2) = this%vBoundaries(this%coCell+1) - this%coVel 

end subroutine calculateInterps
!-----------------------------------------------------------------------------------------------------------------------------------
module function interpMom(this,order,f) result(res)
    !! Calculate moment of single local harmonic with interpolated grid points near cut-off

    class(ModelboundLBCData)              ,intent(inout)    :: this 
    integer(ik)                           ,intent(in)       :: order 
    real(rk) ,dimension(:)                ,intent(in)       :: f

    real(rk) :: res
    
    if (assertions) call assert(this%isDefined(),"interpMom called on undefined modelbound LBC data")

    res = 0 

    res = res + dot_product(f(1:this%coCell-2),this%vGridCopy(1:this%coCell-2)**(order+2)*this%dvCopy(1:this%coCell-2))

    res = res + f(this%coCell-1) * this%interpCoord(1)**2*this%interpdv(1)*this%vGridCopy(this%coCell-1)**order

    res = res + f(this%coCell) * this%interpCoord(2)**2*this%interpdv(2)*this%vGridCopy(this%coCell)**order

    res = res + dot_product(f(this%coCell+1:size(f)),this%vGridCopy(this%coCell+1:size(f))**(order+2)&
                *this%dvCopy(this%coCell+1:size(f)))

    res = res * 4*pi

end function interpMom
!-----------------------------------------------------------------------------------------------------------------------------------
module function boundaryHarmonic(this,lNum) result(res)
    !! Calculate lNum harmonic at boundary (assuming m=0) using Pll tensor and interpolation

    class(ModelboundLBCData)              ,intent(inout)    :: this 
    integer(ik)                           ,intent(in)       :: lNum 

    real(rk) ,allocatable ,dimension(:) :: res

    integer(ik) :: i

    if (assertions) call assert(this%isDefined(),"boundaryHarmonic called on undefined modelbound LBC data")

    allocate(res(size(this%vGridCopy)))

    res = 0 

    do i = 1,this%coCell-2
        res(i) = dot_product(this%bufferPLL(i,lNum,:),this%fext(:,i))
    end do

    !Left interpolated cell
    res(this%coCell-1) = dot_product(this%bufferPLL(this%coCell-1,lNum,:),this%vInterp(1) * this%fext(:,this%coCell-1) &
                                                                    + (real(1,kind=rk)-this%vInterp(1))*this%fext(:,this%coCell))
    
    !Right interpolated cell
    if (this%coCell == size(res)) then 
        res(this%coCell) = dot_product(this%bufferPLL(this%coCell,lNum,:),this%vInterp(2) * this%fext(:,this%coCell))

    else
        res(this%coCell) = dot_product(this%bufferPLL(this%coCell,lNum,:),this%vInterp(2) * this%fext(:,this%coCell) &
        + (real(1,kind=rk)-this%vInterp(2))*this%fext(:,this%coCell+1))
    end if

    do i = this%coCell+1,size(res)
        res(i) = dot_product(this%bufferPLL(i,lNum,:),this%fext(:,i))
    end do
                                
end function boundaryHarmonic
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calculateCutOff(this,ionCurrent,totalCurrent)
    !! Calculate cut-off cell and velocity based on modelbound data and currents

    class(ModelboundLBCData)              ,intent(inout)    :: this 
    real(rk)                              ,intent(in)       :: ionCurrent !! ion current into sheath 
    real(rk)                              ,intent(in)       :: totalCurrent !! total current into sheath to be matched

    integer(ik) :: i 

    real(rk) :: elFlux ,v1,v2,y1, y2
    if (assertions) call assert(this%isDefined(),"calculateCutOff called on undefined modelbound LBC data")

    do i = size(this%vGridCopy),2,-1
        this%coVel = this%vBoundaries(i)
        this%coCell = i
        call this%calculateInterps()
        call this%calculatePll([1])

        elFlux = this%interpMom(1,this%boundaryHarmonic(1))/real(3,kind=rk)
        if (abs(elFlux) > abs(ionCurrent - totalCurrent)) exit
    end do

    if (this%coCell == 2) error stop "calculateCutOff cannot match logical boundary condition fluxes"

    !Initialize false position method variables 

    v1 = this%coVel
    y1 = elFlux - ionCurrent + totalCurrent
    v2 = this%vBoundaries(this%coCell+1)
    this%coVel = v2
    call this%calculateInterps()
    call this%calculatePll([1])

    elFlux = this%interpMom(1,this%boundaryHarmonic(1))/real(3,kind=rk)

    y2 = elFlux - ionCurrent + totalCurrent

    do i = 1,30
        this%coVel = (v1*y2-v2*y1)/(y2-y1)
        call this%calculateInterps()
        call this%calculatePll([1])
        elFlux = this%interpMom(1,this%boundaryHarmonic(1))/real(3,kind=rk)
        if (abs(elFlux-ionCurrent+ totalCurrent)/abs(ionCurrent- totalCurrent) < this%bisTol) exit 

        if ((elFlux-ionCurrent+ totalCurrent) * y1 < 0 ) then 

            v2 = this%coVel
            y2 = elFlux-ionCurrent+ totalCurrent

        else
            v1 = this%coVel
            y1 = elFlux-ionCurrent+ totalCurrent
        end if

    end do

end subroutine 
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCoCell(this) result(res)
    !! Getter for the cut-off cell coordinate

    class(ModelboundLBCData)              ,intent(in)    :: this
    integer(ik)                                          :: res

    if (assertions) call assertPure(this%isDefined(),"getCoCell called on undefined modelbound LBC data")

    res = this%coCell

end function getCoCell
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getPll(this,lNum) result(res)
    !! Getter for the cut-off (lNum,:,:) decomposition tensor

    class(ModelboundLBCData)             ,intent(in)    :: this
    integer(ik)                          ,intent(in)    :: lNum 
    real(rk) ,allocatable ,dimension(:,:)               :: res

    if (assertions) call assertPure(this%isDefined(),"getPll called on undefined modelbound LBC data")

    res = this%bufferPLL(:,lNum,:)

end function getPll
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpCoords(this) result(res)
    !! Getter for the interpolated cut-off cell coordinates

    class(ModelboundLBCData)  ,intent(in)    :: this
    real(rk) ,dimension(2)                   :: res

    if (assertions) call assertPure(this%isDefined(),"getInterpCoords called on undefined modelbound LBC data")
    
    res = this%interpCoord

end function getInterpCoords
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpWidths(this) result(res)
    !! Getter for the interpolated cut-off cell widths

    class(ModelboundLBCData)  ,intent(in)    :: this
    real(rk) ,dimension(2)                   :: res

    if (assertions) call assertPure(this%isDefined(),"getInterpWidths called on undefined modelbound LBC data")

    res = this%interpdv 
    
end function getInterpWidths
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getInterpCoeffs(this) result(res)
    !! Getter for the interpolation coefficients

    class(ModelboundLBCData)  ,intent(in)    :: this
    real(rk) ,dimension(2)                   :: res 

    if (assertions) call assertPure(this%isDefined(),"getInterpCoeffs called on undefined modelbound LBC data")

    res = this%vInterp

end function getInterpCoeffs
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule modelbound_lbc_data_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
