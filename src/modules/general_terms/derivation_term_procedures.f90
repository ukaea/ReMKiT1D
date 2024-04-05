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
submodule (derivation_explicit_term_class) derivation_term_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the derivation term class 

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDerivationTerm(this,gridObj,partitionObj,procRank,evolvedVar,varCont,derivObj,derivIndices,mbVarName)
    !! Explicit derivation term initialization routine 

    class(DerivationTerm)                       ,intent(inout)  :: this
    type(Grid)                                  ,intent(in)     :: gridObj
    type(Partition)                             ,intent(in)     :: partitionObj !! Parition object used to determine local number of DoF
    integer(ik)                                 ,intent(in)     :: procRank !! Current processor rank
    character(*)                                ,intent(in)     :: evolvedVar !! Name of evolved variable
    type(VariableContainer)                     ,intent(in)     :: varCont !! Reference variable container
    class(Derivation)                           ,intent(in)     :: derivObj !! Derivation object used by the term 
    integer(ik)              ,dimension(:)      ,intent(in)     :: derivIndices !! Required variable indices for the derivation object
    character(*)          ,optional             ,intent(in)     :: mbVarName !! Optional modelbound variable with which to multiply the derivation result 

    integer(ik) :: minX, maxX 

    if (assertions .or. assertionLvl >= 0) then 
        call assert(partitionObj%isDefined(),"Undefined partition object passed to derivation term constructor")
        call assert(varCont%isDefined(),"Undefined variable container object passed to derivation term constructor")
    end if 

    minX = partitionObj%getMinXAtInd(procRank+1)
    maxX = partitionObj%getMaxXAtInd(procRank+1)
    this%locNumX = maxX - minX + 1
    this%numV = gridObj%getNumV()
    this%numH = gridObj%getNumH()

    this%kineticRow = varCont%isVarDist(varCont%getVarIndex(evolvedVar))

    if (assertions .or. assertionLvl >=0) &
        call assert(.not. this%kineticRow, &
        "Derivation terms do not currently support evolving distribution variables")

    this%isActive = partitionObj%getMinHAtInd(procRank+1) == 1

    call this%setEvolvedVar(evolvedVar)

    allocate(this%derivObj,source=derivObj)
    this%derivReqIndices = derivIndices

    allocate(this%resultBuffer,source=varCont%variables(varCont%getVarIndex(evolvedVar))%entry)
    this%resultBuffer = real(0,kind=rk)
    if (present(mbVarName)) then 
        this%mbVarName = mbVarName
    end if

    call this%makeDefined()

end subroutine initDerivationTerm
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine updateDerivationTerm(this,varCont,modelData,hostModel) 
    !! Update function, does not call update on operator, assuming it is never set. Only updates the modelbound data buffer if
    !! required

    class(DerivationTerm)           ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update 
    class(ModelboundData) ,optional ,intent(in)     :: modelData !! Reference model data - unused
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Host model - unused

    real(rk)        ,allocatable ,dimension(:)      :: modelboundDataVals
    integer(ik)                                     :: mbHalo ,lboundBuffer ,mbDataDim ,bufferHalo ,haloDiff ,mbLBound

    if (assertions) then 
        call assert(this%isDefined(),"Attempted to update undefined derivation term object")
        call assert(varCont%isDefined(),"Attempted to update derivation term object using undefined variable container")

        if (allocated(this%mbVarName))&
        call assert(present(modelData),&
        "Attempted to updated derivation term object with modelbound dependencies without passing modelData")

    end if
    if (this%isActive) then 

        this%resultBuffer = this%derivObj%calculate(varCont%variables,this%derivReqIndices)

        if (allocated(this%mbVarName)) then
            !Determine buffer halo size
            lboundBuffer = lbound(this%resultBuffer,1)
            bufferHalo = 1 - lboundBuffer

            call modelData%copyData(this%mbVarName,modelboundDataVals)
            mbDataDim = modelData%getDataDim(this%mbVarName)
            mbLBound = lbound(modelboundDataVals,1)
            select case (mbDataDim)
            case (0)
                this%resultBuffer = this%resultBuffer *  modelboundDataVals(1)
            case (1) 
                mbHalo = (size(modelboundDataVals) - this%locNumX)/2
                haloDiff = bufferHalo - mbHalo
                this%resultBuffer(lboundBuffer+haloDiff:lboundBuffer&
                +haloDiff+size(modelboundDataVals)-1) = &
                this%resultBuffer(lboundBuffer+haloDiff:lboundBuffer&
                +haloDiff+size(modelboundDataVals)-1) *  modelboundDataVals

            case default
                error stop "Unsupported dimensionality of rank 1 modelbound data detected in derivation term"
            end select
        end if

    end if


end subroutine updateDerivationTerm
!-----------------------------------------------------------------------------------------------------------------------------------
module function derivationOuterFun(this,varCont) result(res)
    !! Outer function calling the contained derivation object, optionally multiplying the result with the modelbound variable
    !! buffer

    class(DerivationTerm)           ,intent(in)   :: this
    type(VariableContainer)         ,intent(in)   :: varCont

    real(rk) ,allocatable           ,dimension(:) :: res  

    if (assertions) call assert(this%isDefined(),"derivationOuterFun called on undefined term")

    res = this%resultBuffer

end function derivationOuterFun
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule derivation_term_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
