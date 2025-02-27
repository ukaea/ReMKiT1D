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
submodule (fluid_stencil_templates) fluid_stencil_templates_procedures
!! author: Stefan Mijin 
!! 
!! Contains fluid stencil template procedures

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initFluidStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar 
    character(*)               ,intent(in)    :: implicitVar 


    type(NamedString) ,dimension(1) :: stencilType 

    stencilType(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyStencilType,"")
    call envObj%jsonCont%load(stencilType)
    call envObj%jsonCont%output(stencilType)

    select case(stencilType(1)%value)
    case (keyDiagonal)
        call initDiagonalStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyCentralDiffInterp)
        call initCentralDifferenceInterpTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyStaggeredDiff)
        call initStaggeredDifferenceTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyUpwindedDiff)
        call initUpwindingDifferenceTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyBoundaryStencil)
        call initBCTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyDiffusionStencil)
        call initDiffusionStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyTermMomentStencil)
        call initTermMomentStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyCustomFluid1DStencil)
        call initCustomFluid1DStencil(stencilTemplateObj,envObj,jsonPrefix)
    case default
        error stop "Unsupported stencil type detected by initFluidStencilTemplate"
    end select

end subroutine initFluidStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDiagonalStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize diagonal tencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedIntegerArray) ,dimension(1) :: evolvedXCells


    evolvedXCells(1)%name = jsonPrefix//"."//keyStencilData//"."//keyEvolvedXCells
    allocate(evolvedXCells(1)%values(0))

    call envObj%jsonCont%load(evolvedXCells)
    call envObj%jsonCont%output(evolvedXCells)

    if (size(evolvedXCells(1)%values) > 0) then 
        call initDiagonalStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,evolvedXCells(1)%values)
    else 
        call initDiagonalStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar)
    end if

end subroutine initDiagonalStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDiagonalStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,evolvedXCells)
    !! Initialize diagonal stencil template based on direct inputs

    type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)         ,intent(inout) :: envObj
    character(*)                     ,intent(in)    :: evolvedVar
    character(*)                     ,intent(in)    :: implicitVar
    integer(ik) ,optional ,dimension(:) ,intent(in) :: evolvedXCells

    type(InterpStencilGenerator) :: interpStencilGen

    integer(ik) :: i
    logical     :: pGrid, staggeredRowVar ,staggeredColVar 

    integer(ik) ,allocatable ,dimension(:) :: validXCoords ,usedCoords

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    pGrid = envObj%geometryObj%isPeriodic()

    validXCoords = [(i,i=1,envObj%gridObj%getNumX())]
    if (staggeredRowVar .and. (.not. pGrid)) &
    validXCoords = [(i,i=1,envObj%gridObj%getNumX()-1)]

    if (present(evolvedXCells)) then 
        usedCoords = evolvedXCells
    else
        usedCoords = validXCoords
    end if 

    if (assertions .or. assertionLvl >= 0) then 
        do i = 1,size(usedCoords)
            call assert(any(validXCoords == usedCoords(i)),&
            "Unallowed value detected in evolvedXCells in initDiagonalStencilTemplateDirect")
        end do
    end if

    stencilTemplateObj%fixedStencil = .true.
    stencilTemplateObj%rowCoords = reshape(usedCoords,[1,size(usedCoords)])

    if (staggeredRowVar .eqv. staggeredColVar) then 

        call stencilTemplateObj%defaultStencil%init(xStencil=[0],xPeriodic=pGrid)

        ! No stencil gen since default multConst is all ones
    else

        if (staggeredColVar) then 

            call stencilTemplateObj%defaultStencil%init(xStencil=[-1,0],xPeriodic=pGrid)
            if (.not. pGrid) then !Handle extrapolation points
                
                allocate(stencilTemplateObj%overridingStencils(2))
                call stencilTemplateObj%overridingStencils(1)%init(xStencil=[0,1],xPeriodic=pGrid)
                call stencilTemplateObj%overridingStencils(2)%init(xStencil=[-2,-1],xPeriodic=pGrid)
                stencilTemplateObj%overridingStencilCoords = reshape([1,envObj%gridObj%getNumX()],[1,2])

            end if
        else

            call stencilTemplateObj%defaultStencil%init(xStencil=[0,1],xPeriodic=pGrid)

        end if

        call interpStencilGen%init(envObj%partitionObj,envObj%geometryObj,envObj%mpiCont%getWorldRank(),&
                                   staggeredGridMode=staggeredColVar)

        allocate(stencilTemplateObj%stencilGen,source=interpStencilGen)
        
    end if

end subroutine initDiagonalStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCentralDifferenceInterpTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize central difference with face interpolation stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedLogical) ,dimension(1) :: ignoreJacobian 
    type(NamedString)  ,dimension(1) :: interpolatedVarName

    ignoreJacobian(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyIgnoreJacobian,.false.)
    interpolatedVarName(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyInterpolatedVar,keyNone)

    call envObj%jsonCont%load(ignoreJacobian)
    call envObj%jsonCont%output(ignoreJacobian)
    call envObj%jsonCont%load(interpolatedVarName)
    call envObj%jsonCont%output(interpolatedVarName)

    if (interpolatedVarName(1)%value /= keyNone) then 
        call initCentralDifferenceInterpTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                        ignoreJacobian(1)%value,interpolatedVarName(1)%value)
    else
        call initCentralDifferenceInterpTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian(1)%value)
    end if

    
end subroutine initCentralDifferenceInterpTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCentralDifferenceInterpTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                            ignoreJacobian,interpolatedVarName)
    !! Initialize central difference with face interpolation stencil template based on direct inputs

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    logical                    ,intent(in)    :: ignoreJacobian
    character(*) ,optional     ,intent(in)    :: interpolatedVarName

    type(UWCDiffStencilValGenerator) :: stencilGen

    integer(ik) :: i
    logical     :: pGrid, staggeredRowVar ,staggeredColVar 

    real(rk) ,allocatable ,dimension(:) :: dx ,centreJ ,rightJ ,linInterp ,outerJ,innerJ

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions .or. assertionLvl >= 0) then 
        call assert(staggeredRowVar .eqv. staggeredColVar, "centralDifferenceInterpolated stencil requires both the&
        & evolved and implicit variables to be defined on the same grid")
        if (present(interpolatedVarName)) then 
        call assert(envObj%externalVars%isVarNameRegistered(interpolatedVarName),interpolatedVarName//&
        " not registered in passed environment wrapper")
        end if
    end if
    pGrid = envObj%geometryObj%isPeriodic()

    call stencilTemplateObj%defaultStencil%init(xStencil = [-1,0,1],xPeriodic=pGrid)

    stencilTemplateObj%fixedStencil = .not. present (interpolatedVarName) 
    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX())],[1,envObj%gridObj%getNumX()])


    if (staggeredRowVar .and. (.not. pGrid)) &
    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX()-1)],[1,envObj%gridObj%getNumX() - 1]) !Remove last row in staggered non-periodic case

    allocate(linInterp,source=envObj%geometryObj%getLinInterp(dualGrid=staggeredRowVar))
    dx = envObj%geometryObj%getCellWidths(dualGrid=staggeredRowVar)
    centreJ = envObj%geometryObj%getJacobianCentre(dualGrid=staggeredRowVar)
    rightJ = envObj%geometryObj%getJacobianRight(dualGrid=staggeredRowVar)

    if (ignoreJacobian) then 
        outerJ = 1/dx
        innerJ = [(real(1,kind=rk),i=1,size(outerJ))]
    else
        outerJ = 1/(dx*centreJ)
        innerJ = rightJ
    end if

    if (present(interpolatedVarName)) then
        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,linInterp,xPeriodic=pGrid,&
                        interpVarIndex=envObj%externalVars%getVarIndex(interpolatedVarName),upwindingMode=0,&
                        staggeredGridMode=staggeredRowVar) 
    else
        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,linInterp,xPeriodic=pGrid,&
                         staggeredGridMode=staggeredRowVar) 
    end if

    allocate(stencilTemplateObj%stencilGen,source=stencilGen)

end subroutine initCentralDifferenceInterpTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStaggeredDifferenceTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize staggered difference with stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedLogical) ,dimension(1) :: ignoreJacobian 
    
    ignoreJacobian(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyIgnoreJacobian,.false.)

    call envObj%jsonCont%load(ignoreJacobian)
    call envObj%jsonCont%output(ignoreJacobian)

    call initStaggeredDifferenceTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian(1)%value)

end subroutine initStaggeredDifferenceTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStaggeredDifferenceTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian)
    !! Initialize staggered difference with stencil template based on direct input 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    logical                    ,intent(in)    :: ignoreJacobian

    type(FBDiffStencilValGenerator) :: stencilGen

    integer(ik) :: i

    logical     :: pGrid, staggeredRowVar ,staggeredColVar 

    real(rk) ,allocatable ,dimension(:) :: dx ,centreJ ,rightJ ,outerJ,innerJ

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions .or. assertionLvl >= 0) call assert(staggeredRowVar .neqv. staggeredColVar,&
     "staggeredDifference stencil requires that the evolved and&
    & implicit variables be defined on different grids")

    pGrid = envObj%geometryObj%isPeriodic()

    stencilTemplateObj%fixedStencil = .true.
    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX())],[1,envObj%gridObj%getNumX()])

    if (staggeredRowVar .and. (.not. pGrid)) &
    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX()-1)],[1,envObj%gridObj%getNumX()-1]) !Remove last row in staggered non-periodic case

    dx = envObj%geometryObj%getCellWidths(dualGrid=staggeredRowVar,extendedBoundaryCells=.false.)
    centreJ = envObj%geometryObj%getJacobianCentre(dualGrid=staggeredRowVar,extendedBoundaryCells=.false.)
    rightJ = envObj%geometryObj%getJacobianRight(dualGrid=staggeredRowVar,extendedBoundaryCells=.false.)

    if (ignoreJacobian) then 
        outerJ = 1/dx
        innerJ = [(real(1,kind=rk),i=1,size(outerJ))]
    else
        outerJ = 1/(dx*centreJ)
        innerJ = rightJ
    end if

    if (staggeredRowVar) then

        call stencilTemplateObj%defaultStencil%init(xStencil=[0,1],xPeriodic=pGrid)

        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,&
                             xPeriodic=pGrid,staggeredGridMode=.true.) 

    else
        call stencilTemplateObj%defaultStencil%init(xStencil=[-1,0],xPeriodic=pGrid)

        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,&
                             xPeriodic=pGrid,backwardsDiff=.true.,staggeredGridMode=.true.) 
        
    end if

    allocate(stencilTemplateObj%stencilGen,source=stencilGen)


end subroutine initStaggeredDifferenceTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initUpwindingDifferenceTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize upwinding difference with stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedLogical) ,dimension(1) :: ignoreJacobian 
    type(NamedString)  ,dimension(1) :: fluxJacVar

    ignoreJacobian(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyIgnoreJacobian,.false.)
    fluxJacVar(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyFluxJacVar,"")

    call envObj%jsonCont%load(ignoreJacobian)
    call envObj%jsonCont%load(fluxJacVar)
    call envObj%jsonCont%output(ignoreJacobian)
    call envObj%jsonCont%output(fluxJacVar)

    call initUpwindingDifferenceTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                               ignoreJacobian(1)%value,fluxJacVar(1)%value)

end subroutine initUpwindingDifferenceTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initUpwindingDifferenceTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian,fluxJacVar)
    !! Initialize upwinding difference with stencil template based on direct inputs

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    logical                    ,intent(in)    :: ignoreJacobian
    character(*)               ,intent(in)    :: fluxJacVar

    type(UWCDiffStencilValGenerator) :: stencilGen

    integer(ik) :: i
    logical     :: pGrid, staggeredRowVar ,staggeredColVar 

    real(rk) ,allocatable ,dimension(:) :: dx ,centreJ ,rightJ ,linInterp ,outerJ,innerJ

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions .or. assertionLvl >= 0) call assert(staggeredRowVar .eqv. staggeredColVar,&
     "upwindingDifference stencil requires both the&
    & evolved and implicit variables to be defined on the same grid")

    if (assertions .or. assertionLvl >= 0) then 
            call assert(envObj%externalVars%isVarNameRegistered(fluxJacVar),fluxJacVar//&
                                    " not registered in passed environment wrapper")
    end if

    pGrid = envObj%geometryObj%isPeriodic()

    call stencilTemplateObj%defaultStencil%init(xStencil = [-1,0,1],xPeriodic=pGrid)

    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX())],[1,envObj%gridObj%getNumX()])

    call stencilTemplateObj%defaultStencil%init(xStencil = [-1,0,1],xPeriodic=pGrid)

    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX())],[1,envObj%gridObj%getNumX()])

    if (staggeredRowVar .and. (.not. pGrid)) &
    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX()-1)],[1,envObj%gridObj%getNumX() - 1]) !Remove last row in staggered non-periodic case

    allocate(linInterp,source=envObj%geometryObj%getLinInterp(dualGrid=staggeredRowVar))
    dx = envObj%geometryObj%getCellWidths(dualGrid=staggeredRowVar)
    centreJ = envObj%geometryObj%getJacobianCentre(dualGrid=staggeredRowVar)
    rightJ = envObj%geometryObj%getJacobianRight(dualGrid=staggeredRowVar)

    if (ignoreJacobian) then 
        outerJ = 1/dx
        innerJ = [(real(1,kind=rk),i=1,size(outerJ))]
    else
        outerJ = 1/(dx*centreJ)
        innerJ = rightJ
    end if

        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,linInterp,xPeriodic=pGrid,&
                            interpVarIndex=envObj%externalVars%getVarIndex(fluxJacVar),upwindingMode=1,&
                            staggeredGridMode=staggeredRowVar) 

    allocate(stencilTemplateObj%stencilGen,source=stencilGen)

end subroutine initUpwindingDifferenceTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initBCTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize flux-like boundary condition stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedLogical) ,dimension(1) :: ignoreJacobian
    type(NamedLogical) ,dimension(1) :: leftBoundary
    type(NamedLogical) ,dimension(1) :: dontExtrapolate
    type(NamedLogical) ,dimension(1) :: noLowerBound

    type(NamedString)  ,dimension(1) :: fluxJacVar
    type(NamedString)  ,dimension(1) :: lowerBoundVar

    type(NamedReal)    ,dimension(1) :: fixedLowerBound

    ignoreJacobian(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyIgnoreJacobian,.false.)
    leftBoundary(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyLeftBoundary,.false.)
    dontExtrapolate(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyDontExtrapolate,.false.)
    noLowerBound(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyNoLowerBound,.false.)

    fluxJacVar(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyFluxJacVar,keyNone)
    lowerBoundVar(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyLowerBoundVar,keyNone)

    fixedLowerBound(1) = NamedReal(jsonPrefix//"."//keyStencilData//"."//keyFixedLowerBound,0)

    call envObj%jsonCont%load(ignoreJacobian)
    call envObj%jsonCont%load(leftBoundary)
    call envObj%jsonCont%load(dontExtrapolate)
    call envObj%jsonCont%load(noLowerBound)
    call envObj%jsonCont%load(fluxJacVar)
    call envObj%jsonCont%load(lowerBoundVar)
    call envObj%jsonCont%load(fixedLowerBound)
    call envObj%jsonCont%output(ignoreJacobian)
    call envObj%jsonCont%output(leftBoundary)
    call envObj%jsonCont%output(dontExtrapolate)
    call envObj%jsonCont%output(noLowerBound)
    call envObj%jsonCont%output(fluxJacVar)
    call envObj%jsonCont%output(lowerBoundVar)
    call envObj%jsonCont%output(fixedLowerBound)

    if (fluxJacVar(1)%value /= keyNone) then 

        if (lowerBoundVar(1)%value /= keyNone) then
            call initBCTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian(1)%value,&
                                    leftBoundary(1)%value,dontExtrapolate(1)%value &
                                    ,noLowerBound(1)%value,fixedLowerBound(1)%value,fluxJacVar(1)%value,lowerBoundVar(1)%value)
        else 
            call initBCTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian(1)%value,&
                                    leftBoundary(1)%value,dontExtrapolate(1)%value &
                                    ,noLowerBound(1)%value,fixedLowerBound(1)%value,fluxJacVar(1)%value)
        end if
    else
        call initBCTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian(1)%value,&
                                    leftBoundary(1)%value,dontExtrapolate(1)%value &
                                    ,noLowerBound(1)%value,fixedLowerBound(1)%value)
    end if


end subroutine initBCTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initBCTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,ignoreJacobian,leftBoundary,dontExtrapolate&
                                      ,noLowerBound,fixedLowerBound,fluxJacVar,lowerBoundVar)
    !! Initialize flux-like boundary condition stencil template based on direct inputs 
                                      
    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    logical                    ,intent(in)    :: ignoreJacobian
    logical                    ,intent(in)    :: leftBoundary
    logical                    ,intent(in)    :: dontExtrapolate
    logical                    ,intent(in)    :: noLowerBound
    real(rk)                   ,intent(in)    :: fixedLowerBound
    character(*) ,optional     ,intent(in)    :: fluxJacVar
    character(*) ,optional     ,intent(in)    :: lowerBoundVar

    type(UWCBCStencilValGenerator) :: stencilGen

    integer(ik) :: i

    real(rk) ,allocatable ,dimension(:) :: dx ,dxReg ,centreJ ,faceJ ,linInterp 

    real(rk) :: outerJ ,innerJ ,lInterp, lExterp

    logical  :: staggeredRowVar ,staggeredColVar ,staggeredFluxJac,staggeredLowerBound

    integer(ik) ,allocatable ,dimension(:) :: xStencil 

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions .or. assertionLvl >= 0) call assert(staggeredRowVar .eqv. staggeredColVar, "boundaryStencil requires both the&
    & evolved and implicit variables to be defined on the same grid")

    if (assertions .or. assertionLvl >= 0) then 
        if (present(fluxJacVar)) then 
            call assert(envObj%externalVars%isVarNameRegistered(fluxJacVar),fluxJacVar//&
                                    " not registered in passed environment wrapper")

            call assert(.not. envObj%externalVars%isVarOnDualGrid(fluxJacVar),"fluxJacVar in boundaryStencil must be &
            &defined on regular grid")
        end if

        if (present(lowerBoundVar)) then 
            call assert(envObj%externalVars%isVarNameRegistered(lowerBoundVar),lowerBoundVar//&
                                    " not registered in passed environment wrapper")
            call assert(.not. envObj%externalVars%isVarOnDualGrid(fluxJacVar),"lowerBoundVar in boundaryStencil must be &
                                    &defined on regular grid")
        end if

    end if


    if (leftBoundary) then 
        xStencil = [0,1]
        stencilTemplateObj%rowCoords = reshape([1],[1,1])  
    else
        xStencil = [0,-1]
        stencilTemplateObj%rowCoords = reshape([envObj%gridObj%getNumX()],[1,1])
        if (staggeredRowVar) stencilTemplateObj%rowCoords = reshape([envObj%gridObj%getNumX()-1],[1,1])
    end if

    if (dontExtrapolate) xStencil = [0] 

    call stencilTemplateObj%defaultStencil%init(xStencil = xStencil)

    dx = envObj%geometryObj%getCellWidths(dualGrid=staggeredRowVar)
    dxReg = envObj%geometryObj%getCellWidths()
    centreJ = envObj%geometryObj%getJacobianCentre(dualGrid=staggeredRowVar)
    faceJ = envObj%geometryObj%getJacobianRight(dualGrid=staggeredRowVar)
    allocate(linInterp,source=envObj%geometryObj%getLinInterp())
    if (leftBoundary) faceJ = envObj%geometryObj%getJacobianLeft(dualGrid=staggeredRowVar)

    lInterp = linInterp(envObj%gridObj%getNumX()-1)
    if (leftBoundary) lInterp = linInterp(1)
    lExterp = real(1,kind=rk) - lInterp
    if (leftBoundary) lExterp = lInterp

    if (staggeredRowVar) then 
        lExterp = dxReg(envObj%gridObj%getNumX()-1)/dxReg(envObj%gridObj%getNumX()-2)
        if (leftBoundary) lExterp = dxReg(1)/dxReg(2)
    end if

    if (ignoreJacobian) then 
        if (leftBoundary) then
            outerJ = 1/dx(1)
        else
            outerJ = 1/dx(envObj%gridObj%getNumX())
            if (staggeredRowVar) outerJ = 1/dx(envObj%gridObj%getNumX()-1)
        end if
        innerJ =real(1,kind=rk)
    else
        if (leftBoundary) then
            outerJ = 1/(dx(1)*centreJ(1))
            innerJ = faceJ(1)
        else
            outerJ = 1/(dx(envObj%gridObj%getNumX())*centreJ(envObj%gridObj%getNumX()))
            innerJ = faceJ(envObj%gridObj%getNumX())
            if (staggeredRowVar) then 
                outerJ = 1/(dx(envObj%gridObj%getNumX()-1)*centreJ(envObj%gridObj%getNumX()-1))
                innerJ = faceJ(envObj%gridObj%getNumX()-1)
            end if
        end if
    end if

    if (present(fluxJacVar)) then 

        if (present(lowerBoundVar)) then
            call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,lInterp,&
            extrapolate = .not. dontExtrapolate, interpVarIndex = envObj%externalVars%getVarIndex(fluxJacVar),&
            leftBoundary=leftBoundary,&
            lowerBoundVarIndex=envObj%externalVars%getVarIndex(lowerBoundVar),linExterp=lExterp) 
        else if (.not. noLowerBound) then 
            call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,lInterp,&
            extrapolate = .not. dontExtrapolate, interpVarIndex = envObj%externalVars%getVarIndex(fluxJacVar),&
            leftBoundary=leftBoundary,&
            fixedLowerBound=fixedLowerBound,linExterp=lExterp)
        else 
            call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,lInterp,&
            extrapolate = .not. dontExtrapolate, interpVarIndex = envObj%externalVars%getVarIndex(fluxJacVar),&
            leftBoundary=leftBoundary,linExterp=lExterp)
        end if
    else
        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,lInterp,&
        extrapolate = .not. dontExtrapolate,leftBoundary=leftBoundary,linExterp=lExterp)
    end if

    allocate(stencilTemplateObj%stencilGen,source=stencilGen)

end subroutine initBCTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDiffusionStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize diffusion stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedString)      ,dimension(1) :: ruleName
    type(NamedStringArray)  ,dimension(1) :: reqVars
    type(NamedLogical)     ,dimension(1) :: notInterp

    class(Derivation) ,allocatable :: method
    type(CalculationRule) :: calcRule

    integer(ik) :: i

    ruleName(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyRuleName,keyNone)
    reqVars(1)%name = jsonPrefix//"."//keyStencilData//"."//keyReqVarNames
    allocate(reqVars(1)%values(0))

    notInterp = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyDoNotInterpolateD,.false.)

    call envObj%jsonCont%load(ruleName)
    call envObj%jsonCont%load(reqVars)
    call envObj%jsonCont%load(notInterp)
    call envObj%jsonCont%output(ruleName)
    call envObj%jsonCont%output(reqVars)
    call envObj%jsonCont%output(notInterp)

    if (assertions .or. assertionLvl >= 0) then 
        if (ruleName(1)%value /= keyNone) then 
            call assert(envObj%textbookObj%isDerivationRegistered(ruleName(1)%value),ruleName(1)%name//&
                                    " not registered in passed environment wrapper")

            
            do i = 1,size(reqVars(1)%values)
                call assert(envObj%externalVars%isVarNameRegistered(reqVars(1)%values(i)%string),reqVars(1)%values(i)%string//&
                " not registered in passed environment wrapper")
            end do
        end if

    end if

    if (ruleName(1)%value /= keyNone) then 
        call envObj%textbookObj%copyDerivation(ruleName(1)%value,method)
        call calcRule%init(method,reqVars(1)%values)

        call initDiffusionStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,calcRule,&
                                        doNotInterpolateD=notInterp(1)%value)
    else
        call initDiffusionStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,doNotInterpolateD=notInterp(1)%value)
    end if

end subroutine initDiffusionStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDiffusionStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,calcRule,doNotInterpolateD)
    !! Initialize diffusion stencil template based on direct inputs

    type(StencilTemplate)           ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)        ,intent(inout) :: envObj
    character(*)                    ,intent(in)    :: evolvedVar
    character(*)                    ,intent(in)    :: implicitVar
    type(CalculationRule) ,optional ,intent(in)    :: calcRule
    logical ,optional               ,intent(in)    :: doNotInterpolateD

    type(DiffusionStencilValGenerator) :: stencilGen 

    logical     :: pGrid, staggeredRowVar ,staggeredColVar ,notInterp

    integer(ik) :: i ,haloSize
    real(rk) ,allocatable ,dimension(:) :: dx ,dxp ,centreJ ,rightJ ,linInterp ,outerJ,innerJ

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions .or. assertionLvl >= 0) then 
        call assert(.not. staggeredRowVar, "diffusionStencil requires both the&
        & evolved and implicit variables to be defined on the regular grid")
        call assert(.not. staggeredColVar, "diffusionStencil requires both the&
        & evolved and implicit variables to be defined on the regular grid")
    end if

    pGrid = envObj%geometryObj%isPeriodic()
    haloSize = envObj%mpiCont%getXHaloWidth()

    dx = envObj%geometryObj%getCellWidths()
    centreJ = envObj%geometryObj%getJacobianCentre()
    rightJ = envObj%geometryObj%getJacobianRight()
    allocate(linInterp,source=envObj%geometryObj%getLinInterp())
    outerJ = 1/(dx*centreJ)
    dxp = (dx + [dx(2:size(dx)),real(0,kind=rk)])/2 
    innerJ = rightJ/dxp

    notInterp = .false.
    if (present(doNotInterpolateD)) notInterp = doNotInterpolateD

    call stencilTemplateObj%defaultStencil%init(xStencil = [-1,0,1],xPeriodic = pGrid)

    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX())],[1,envObj%gridObj%getNumX()])

    if (present(calcRule)) then 
        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,linInterp&
        ,xPeriodic=pGrid,diffCoeffDerivRule=calcRule,xHaloWidth=haloSize,doNotInterpolateD=notInterp)
    else
        call stencilGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,linInterp&
        ,xPeriodic=pGrid,xHaloWidth=haloSize,doNotInterpolateD=notInterp)
    end if

    allocate(stencilTemplateObj%stencilGen,source=stencilGen)

end subroutine initDiffusionStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCustomFluid1DStencil(stencilTemplateObj,envObj,jsonPrefix)
    !! Initialize custom fluid 1D stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template

    type(NamedIntegerArray)             ,dimension(1) :: xStencil 
    type(NamedStringArray)              ,dimension(1) :: varContColVars ,mbDataColVars
    type(NamedRealArray) ,dimension(:) ,allocatable   :: colVecs 
    character(len=30) :: intToStrBuffer

    integer(ik) :: i

    xStencil(1) = NamedIntegerArray(jsonPrefix//"."//keyStencilData//"."//keyXStencil,[0])
    varContColVars(1)%name = jsonPrefix//"."//keyStencilData//"."//keyColumnVarContVars
    allocate(varContColVars(1)%values(0))
    mbDataColVars(1)%name = jsonPrefix//"."//keyStencilData//"."//keyColumnMBDataVars
    allocate(mbDataColVars(1)%values(0))

    call envObj%jsonCont%load(xStencil)
    call envObj%jsonCont%load(varContColVars)
    call envObj%jsonCont%load(mbDataColVars)
    call envObj%jsonCont%output(xStencil)
    call envObj%jsonCont%output(varContColVars)
    call envObj%jsonCont%output(mbDataColVars)

    allocate(colVecs(size(xStencil(1)%values)))

    do i = 1,size(colVecs)
        intToStrBuffer=""
        write(intToStrBuffer,'(I0)') i
        colVecs(i)%name = jsonPrefix//"."//keyStencilData//"."//keyColumnVec//trim(intToStrBuffer)
        allocate(colVecs(i)%values(0))
    end do
    
    call envObj%jsonCont%load(colVecs)
    call envObj%jsonCont%output(colVecs)

    ! This needs to be done for the stencil generator to correctly understand there are no required variables
    do i = 1,size(varContColVars(1)%values)
        if (varContColVars(1)%values(i)%string == keyNone) deallocate(varContColVars(1)%values(i)%string)
    end do

    do i = 1,size(mbDataColVars(1)%values)
        if (mbDataColVars(1)%values(i)%string == keyNone) deallocate(mbDataColVars(1)%values(i)%string)
    end do

    call initCustomFluid1DStencilDirect(stencilTemplateObj,envObj,xStencil(1)%values,&
                                       removeName(colVecs),varContColVars(1)%values,mbDataColVars(1)%values)
end subroutine initCustomFluid1DStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCustomFluid1DStencilDirect(stencilTemplateObj,envObj&
        ,xStencil,columnVecs,varContColVarNames,mbColVarNames)
    !! Initialize custom fluid 1D stencil template based on direct inputs

    type(StencilTemplate)           ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)        ,intent(inout) :: envObj
    integer(ik) ,dimension(:)       ,intent(in)    :: xStencil 
    type(RealArray) ,dimension(:)   ,intent(in)    :: columnVecs
    type(StringArray) ,dimension(:) ,intent(in)    :: varContColVarNames 
    type(StringArray) ,dimension(:) ,intent(in)    :: mbColVarNames 

    type(FluidStencilGen1D) :: stencilGen 
    type(StencilGenerator1D) :: stencilGenSimple
    type(Stencil1D) :: stencil1DObj

    logical     :: pGrid ,noVariableDependence

    integer(ik) :: i 

    pGrid = envObj%geometryObj%isPeriodic()

    if (assertions .or. assertionLvl >= 0) then 
        call assert(size(xStencil)==size(varContColVarNames),&
        "varContColVarNames in CustomFluid1DStencil initialization must have the same size as the stencil")
        call assert(size(xStencil)==size(mbColVarNames),&
        "mbColVarNames in CustomFluid1DStencil initialization must have the same size as the stencil")
        call assert(size(xStencil)==size(columnVecs),&
        "columnVecs in CustomFluid1DStencil initialization must have the same size as the stencil")
        do i = 1,size(columnVecs) 
            call assert(size(columnVecs(i)%entry)==envObj%gridObj%getNumX(),&
            "Column vectors in CustomFluid1DStencil must conform to x grid")
        end do
    end if

    call stencilTemplateObj%defaultStencil%init(xStencil = xStencil,xPeriodic = pGrid)

    stencilTemplateObj%rowCoords = reshape([(i,i=1,envObj%gridObj%getNumX())],[1,envObj%gridObj%getNumX()])

    call stencil1DObj%init(xStencil)

    noVariableDependence = .true. 

    do i = 1,size(xStencil)
        if (allocated(varContColVarNames(i)%string)) noVariableDependence = .false.
        if (allocated(mbColVarNames(i)%string)) noVariableDependence = .false.
    end do

    if (noVariableDependence) then
        call stencilGenSimple%init(stencil1DObj,columnVecs,periodicDim=pGrid,&
                                    coordInterval=[envObj%partitionObj%getMinXAtInd(envObj%mpiCont%getWorldRank()+1),&
                                                envObj%partitionObj%getMaxXAtInd(envObj%mpiCont%getWorldRank()+1)])
        stencilTemplateObj%fixedStencil = .true.
        allocate(stencilTemplateObj%stencilGen,source=stencilGenSimple)
    else
        call stencilGen%init(stencil1DObj,columnVecs,varContColVarNames,mbColVarNames,periodicDim=pGrid,&
                            coordInterval=[envObj%partitionObj%getMinXAtInd(envObj%mpiCont%getWorldRank()+1),&
                                        envObj%partitionObj%getMaxXAtInd(envObj%mpiCont%getWorldRank()+1)])

        allocate(stencilTemplateObj%stencilGen,source=stencilGen)
    end if

end subroutine initCustomFluid1DStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule fluid_stencil_templates_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
