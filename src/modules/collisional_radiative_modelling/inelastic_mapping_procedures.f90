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
submodule (inelastic_mapping) inelastic_mapping_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with inelastic collision mappings

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function inelWeights(space,E) result(wMat)
    !! Calculate velocity space inelastic mapping weight matrix for transition with energy E on given velocity space. Assumes that 
    !! velocity and energy are suitably normalized (i.e. v^2 is an energy)

    type(VSpace)        ,intent(in) :: space !! Velocity space object to calculate the mapping on
    real(rk)            ,intent(in) :: E    !! Mapping transition energy
    type(SparseRowData)             :: wMat

    real(rk) ,allocatable ,dimension(:) :: vGrid ,vWidths ,eGrid ! vGrid data, eGrid is just vGrid**2

    logical ,allocatable ,dimension(:,:) :: wMap ! Logical map of emitters to absorbers, true where there is an allowed transition

    integer(ik) ,allocatable ,dimension(:,:) :: idealAbsorbers

    integer(ik) :: i, j, numV ,idealEmitter ,firstCell,secondCell ,numPairs 

    integer(ik) :: deltaFirst ,deltaSecond !Numbers of absorber pairs in which first/second ideal absorber are present

    type(IntArray) ,allocatable ,dimension(:) :: absorberList

    real(rk) :: aListEnergyWidth, & !Sum of energy widths of all cells in particular absorber list
                gammaP, & !Pair emission fraction
                eInterp 

    real(rk) ,allocatable ,dimension(:,:) :: denseWMat

    if (assertions) call assertPure(space%isDefined(),"Undefined vSpace passed to inelWeights function")

    vGrid = space%getVGrid()
    vWidths = space%getVCellWidths()
    eGrid = vGrid**2
    numV = size(vGrid)
    allocate(wMap(numV,numV))
    allocate(denseWMat(numV,numV))
    denseWMat = 0
    wMap = .false.
    
    !Find ideal absorbers for each eligible point and set wMap to true where both absorbers exist
    allocate(idealAbsorbers(2,numV))
    idealAbsorbers = 0
    do i = 1,numV 
        if (eGrid(i) >= E) idealAbsorbers(:,i) = space%getNearestPoints(sqrt(eGrid(i) - E))
        if (all(idealAbsorbers(:,i) > 0)) wMap(idealAbsorbers(:,i),i) = .true.
    end do

    !Find ideal emitters for each point and update wMap accordingly 
    
    do i  = 1,numV
        idealEmitter = 0
        if (eGrid(i) >= - E) then 
            idealEmitter = space%getContainingVCell(sqrt(eGrid(i) + E))
            ! If ideal emitter is found, make sure that it has enough energy to emit to cell i
            if (idealEmitter > 0) then 
                if (eGrid(idealEmitter)-eGrid(i) > E) then !if superelastic emitter doesn't have enough energy
                    idealEmitter = idealEmitter + (int(sign(real(1,kind=rk),E),kind=ik) - 1)/2 !if E is negative shift emitter to left 
                else !if inelastic emitter doesn't have enough energy
                    idealEmitter = idealEmitter + (int(sign(real(1,kind=rk),E),kind=ik) + 1)/2 !if E is positive shift emitter to
                end if
            end if
            !Correct for idealEmitter > numV case 
            if (idealEmitter > numV) idealEmitter = 0
        end if
        if (idealEmitter > 0) then 
            !Register i as an absorber of its ideal emitter
            wMap(i,idealEmitter) = .true.
        else
            !If no ideal emitter is found, remove cell i from all absorber pairs
            wMap(i,:) = .false. 
        end if
    end do

    !Final pass through wMap removes all emitters that have lost one of their ideal absorbers after last step

    do i = 1,numV 
        if (all(idealAbsorbers(:,i) > 0)) then 
            if (.not. all(wMap(idealAbsorbers(:,i),i))) wMap(:,i) = .false.
        else
            wMap(:,i) = .false.
        end if
    end do

    allocate(absorberList(numV))
    do i = 1,numV 
        !Fill out total absorber list
        absorberList(i)%entry = findIndices(wMap(:,i))
        numPairs = size(absorberList(i)%entry)
        if (numPairs > 0) numPairs = numPairs - 1
        if (numPairs > 0) then 
            aListEnergyWidth = 2*sum(vGrid(absorberList(i)%entry)*vWidths(absorberList(i)%entry))
            deltaFirst = count(wMap(idealAbsorbers(2,i):numV,i)) !Number of times first ideal absorber appears in a pair
            deltaSecond = count(wMap(1:idealAbsorbers(1,i),i))   !Number of times second ideal absorber appearts in a pair
            do j = 1, numPairs 
                if (absorberList(i)%entry(j) <= idealAbsorbers(1,i)) then !Pairs which contain points <= to the first ideal absorber
                    firstCell = absorberList(i)%entry(j)
                    secondCell = idealAbsorbers(2,i)
                    gammaP = 2 * vGrid(firstCell) * vWidths(firstCell)  !Contribution from first (non-ideal) absorber
                                                         
                    if (absorberList(i)%entry(j) == idealAbsorbers(1,i)) gammaP = gammaP / deltaFirst !Handle ideal absorber pair case
                    
                    gammaP = gammaP + & !Contribution from second (ideal) absorber
                    2 * vGrid(secondCell) * vWidths(secondCell) / deltaSecond
                    
                else !Pairs which include points greater than the second ideal absorber
                    firstCell = idealAbsorbers(1,i)
                    secondCell = absorberList(i)%entry(j+1)

                    gammaP = 2 * vGrid(firstCell) & !Contribution from first (ideal) absorber
                                                         * vWidths(firstCell) / deltaFirst

                    gammaP = gammaP + & !Contribution from second (non-ideal) absorber
                    2 * vGrid(secondCell) * vWidths(secondCell)
                end if

                gammaP = gammaP / aListEnergyWidth

                eInterp = (eGrid(i) - eGrid(firstCell) - E)/(eGrid(secondCell) - eGrid(firstCell))
                denseWMat(firstCell,i) = denseWMat(firstCell,i) & !Add contribution of first cell in this pair to total matrix
                + gammaP * eGrid(i)*vWidths(i) * (real(1.0d0,kind=rk) - eInterp) / (eGrid(firstCell)*vWidths(firstCell))

                denseWMat(secondCell,i) = denseWMat(secondCell,i) & !Add contribution of second cell in this pair to total matrix
                + gammaP * eGrid(i)*vWidths(i) * eInterp / (eGrid(secondCell)*vWidths(secondCell))
            end do
        end if
    end do

    call wMat%init()

    !Add rows to sparse matrix
    do i = 1, numV 
        if (any(wMap(i,:))) then 
            call wMat%addRow(i,findIndices(wMap(i,:)),denseWMat(i,findIndices(wMap(i,:))))
        end if
    end do

end function inelWeights
!-----------------------------------------------------------------------------------------------------------------------------------
!> Fills passed velocity space vector emit with zeros or ones, depending on whether there are any entries in the corresponding 
!> column of passed mapping weight matrix
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine fillEmissionVector(wMat,emit) 

    type(sparseRowData)    ,intent(in) :: wMat
    real(rk) ,dimension(:) ,intent(inout) :: emit

    integer(ik) :: i ,j

    if (assertions) call assertPure(wMat%isDefined(),"Undefined sparse matrix passed to fillEmissionVector")

    emit = 0
    do i = 1,size(emit)
        do j = 1,size(wMat%rowIndex)
            if (any(wMat%columnVector(j)%entry == i)) then 
                emit(i) = real(1,kind=rk) 
                exit
            end if
        end do
    end do

end subroutine fillEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule inelastic_mapping_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
