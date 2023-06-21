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
submodule (kinetic_stencil_templates) kinetic_stencil_templates_procedures
!! author: Stefan Mijin 
!! 
!! Contains kinetic stencil template procedures

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initKineticStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
    !! Initialize kineti stencil template based on existing model, environment object, and JSON file

    type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)         ,intent(inout) :: envObj
    character(*)                     ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)                     ,intent(in)    :: evolvedVar
    character(*)                     ,intent(in)    :: implicitVar
    class(ModelboundData) ,optional ,intent(in)     :: mbData


    type(NamedString) ,dimension(1) :: stencilType 

    stencilType(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyStencilType,"")
    call envObj%jsonCont%load(stencilType)
    call envObj%jsonCont%output(stencilType)

    select case(stencilType(1)%value)
    case (keyDiagonal)
        call initKinDiagonalStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyMomentStencil)
        call initMomentStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyKinSpatialDiffStencil)
        call initSpatialDiffStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyDDVStencil)
        call initDDVStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyVDiffStencil)
        call initVelDiffusionStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyShkarofskyIJ)
        call initIJStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyFixedBoltzmann)
        call assert(present(mbData),"Modelbound data required for "//keyFixedBoltzmann)

        select type(mbData)
        type is (ModelboundCRMData)
            call initFixedBoltzmannStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
        class default
            error stop "Unexpected modelbound data type detected for fixed Boltzmann stencil template"
        end select
    case (keyVariableBoltzmann)
        call assert(present(mbData),"Modelbound data required for "//keyVariableBoltzmann)

        select type(mbData)
        type is (ModelboundCRMData)
            call initVariableBoltzmannStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
        class default
            error stop "Unexpected modelbound data type detected for fixed Boltzmann stencil template"
        end select

    case (keyScalingLBC)
        call initScalingLBCStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case (keyTermMomentStencil)
        call initTermMomentStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    case default
        error stop "Unsupported stencil type detected by initKineticStencilTemplate"
    end select

end subroutine initKineticStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initKinDiagonalStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize diagonal stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedIntegerArray) ,dimension(1) :: evolvedHarmonics,evolvedVelCells,evolvedXCells

    integer(ik) :: i

 
    evolvedXCells(1)%name = jsonPrefix//"."//keyStencilData//"."//keyEvolvedXCells
    allocate(evolvedXCells(1)%values(0))
    evolvedHarmonics(1) = NamedIntegerArray(jsonPrefix//"."//keyStencilData//"."//keyEvolvedHarmonics,&
                                            [(i,i=1,envObj%gridObj%getNumH())])
    evolvedVelCells(1) = NamedIntegerArray(jsonPrefix//"."//keyStencilData//"."//keyEvolvedVCells,&
                                            [(i,i=1,envObj%gridObj%getNumV())])

    call envObj%jsonCont%load(evolvedXCells)
    call envObj%jsonCont%output(evolvedXCells)
    call envObj%jsonCont%load(evolvedHarmonics)
    call envObj%jsonCont%output(evolvedHarmonics)
    call envObj%jsonCont%load(evolvedVelCells)
    call envObj%jsonCont%output(evolvedVelCells)

    call initKinDiagonalStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                              evolvedXCells(1)%values,evolvedHarmonics(1)%values,evolvedVelCells(1)%values)

end subroutine initKinDiagonalStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initKinDiagonalStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                       evolvedXCells,evolvedHarmonics,evolvedVCells)
    !! Initialize diagonal stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    integer(ik) ,dimension(:)  ,intent(in)    :: evolvedXCells
    integer(ik) ,dimension(:)  ,intent(in)    :: evolvedHarmonics
    integer(ik) ,dimension(:)  ,intent(in)    :: evolvedVCells

    type(InterpStencilGenerator) :: interpStencilGen

    integer(ik) :: i
    logical     :: pGrid, staggeredRowVar ,staggeredColVar , interpolationRequired

    logical ,allocatable ,dimension(:) :: oddL

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,validXCoords ,usedCoords ,usedHCoords, usedVCoords

    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen

    if (assertions) call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initKinDiagonalStencilTemplateDirect must be invoked with evolvedVar being a distribution")

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions) then
        if (.not. staggeredRowVar) call assert(.not. staggeredColVar,"If the evolved distribution in &
        &initKinDiagonalStencilTemplateDirect is not staggered the column variable cannot be staggered either")
    end if

    lGrid = envObj%gridObj%getLGrid()
    
    oddL = mod(lGrid(evolvedHarmonics),2) == 1

    if (staggeredRowVar) then 
        call assert(all(oddL) .or. all(.not. oddL),"If the evolved distribution is staggered then the evolved harmonics in &
        &initKinDiagonalStencilTemplate must either all have odd l or all have even l")
    end if

    pGrid = envObj%geometryObj%isPeriodic()

    validXCoords = [(i,i=1,envObj%gridObj%getNumX())]
    if (staggeredRowVar .and. (.not. pGrid) .and. all(oddL)) &
    validXCoords = [(i,i=1,envObj%gridObj%getNumX()-1)]

    if (size(evolvedXCells)>0) then 
        usedCoords = evolvedXCells
    else
        usedCoords = validXCoords
    end if 

    usedHCoords = [(i,i=1,envObj%gridObj%getNumH())]
    if (size(evolvedHarmonics)>0) usedHCoords = evolvedHarmonics

    usedVCoords = [(i,i=1,envObj%gridObj%getNumV())]
    if (size(evolvedVCells)>0) usedVCoords = evolvedVCells   

    if (assertions) then 
        do i = 1,size(evolvedXCells)
            call assert(any(validXCoords == usedCoords(i)),&
            "Unallowed value detected in evolvedXCells in initDiagonalStencilTemplateDirect")
        end do
    end if

    stencilTemplateObj%fixedStencil = .true.
    stencilTemplateObj%rowCoords = allCombinations([IntArray(usedCoords),&
                                                   IntArray(usedHCoords),&
                                                   IntArray(usedVCoords)])

    interpolationRequired = .not. envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)) &
                            .and. staggeredRowVar .and. (all(.not. oddL) .eqv. staggeredColVar)
    if (.not. interpolationRequired) then 

        call stencilTemplateObj%defaultStencil%init(xPeriodic=pGrid,&
             mapToDist=envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)))

        ! No stencil gen since default multConst is all ones
    else

        if (staggeredColVar) then 

            call stencilTemplateObj%defaultStencil%init(xStencil=[-1,0],xPeriodic=pGrid)
            if (.not. pGrid) then !Handle extrapolation points

                stencilTemplateObj%overridingStencilCoords = allCombinations([IntArray([1,envObj%gridObj%getNumX()]),&
                                                             IntArray(usedHCoords),&
                                                             IntArray(usedVCoords)])
                
                allocate(stencilTemplateObj%overridingStencils(size(stencilTemplateObj%overridingStencilCoords,2)))
                do i = 1, size(stencilTemplateObj%overridingStencils)
                    if (stencilTemplateObj%overridingStencilCoords(1,i) == 1) then
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[0,1],xPeriodic=pGrid)
                    else
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[-2,-1],xPeriodic=pGrid)

                    end if
                end do

            end if
        else

            call stencilTemplateObj%defaultStencil%init(xStencil=[0,1],xPeriodic=pGrid)

        end if

        call interpStencilGen%init(envObj%partitionObj,envObj%geometryObj,envObj%mpiCont%getWorldRank(),&
                                   staggeredGridMode=staggeredColVar)

        allocate(genCore)
        call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                         ,stencilTemplateObj%rowCoords,fluidCol=.true.)

        call multStencilGen%init(genCore,fluidCol=.true.)

        call multStencilGen%setXGen(interpStencilGen)

        allocate(stencilTemplateObj%stencilGen,source=multStencilGen)
        
    end if

end subroutine initKinDiagonalStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMomentStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize distribution moment stencil template based on environment object and JSON file

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedInteger) ,dimension(1) :: momentHarmonic, momentOrder


    momentHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyMomentHarmonic,0)
    momentOrder(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyMomentOrder,0)

    call envObj%jsonCont%load(momentHarmonic)
    call envObj%jsonCont%output(momentHarmonic)
    call envObj%jsonCont%load(momentOrder)
    call envObj%jsonCont%output(momentOrder)

    call initMomentStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                         momentHarmonic(1)%value,momentOrder(1)%value)

end subroutine initMomentStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMomentStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,momentHarmonic,momentOrder)
    !! Initialize distribution moment stencil template based on direct input data

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    integer(ik)                ,intent(in)    :: momentHarmonic
    integer(ik)                ,intent(in)    :: momentOrder

    type(InterpStencilGenerator) :: interpStencilGen

    integer(ik) :: i ,usedLocNumX
    logical     :: pGrid, oddL, staggeredRowVar ,staggeredColVar ,interpolationRequired

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,locNumX

    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen

    if (assertions) then 
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)),&
                                "initMomentStencilTemplateDirect must be invoked with implicitVar being a distribution")
        call assert(.not. envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initMomentStencilTemplateDirect must be invoked with evolvedVar being a fluid variable")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions) then
        if (.not. staggeredColVar) call assert(.not. staggeredRowVar,"If the implicit distribution in &
        &initMomentStencilTemplateDirect is not staggered the evolved variable cannot be staggered either")
    end if

    lGrid = envObj%gridObj%getLGrid()

    oddL = mod(lGrid(momentHarmonic),2)==1

    stencilTemplateObj%fixedStencil = .true.
    stencilTemplateObj%rowCoords = allCombinations([IntArray([(i,i=1,envObj%gridObj%getNumX())])])

    pGrid = envObj%geometryObj%isPeriodic()
    if (staggeredRowVar .and. (.not. pGrid)) &
    stencilTemplateObj%rowCoords = allCombinations([IntArray([(i,i=1,envObj%gridObj%getNumX()-1)])]) !Remove last spatial row in staggered non-periodic case

    interpolationRequired = staggeredColVar .and. (oddL .neqv. staggeredRowVar)

    locNumX = envObj%partitionObj%getLocNumX()
    usedLocNumX = locNumX(envObj%mpiCont%getWorldRank()+1)

    if (envObj%partitionObj%getMaxHAtInd(envObj%mpiCont%getWorldRank()+1)==envObj%gridObj%getNumX()) then 
        if (staggeredRowVar .and. (.not. pGrid)) usedLocNumX = usedLocNumX - 1
    end if
    allocate(genCore)
    call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                         ,stencilTemplateObj%rowCoords,&
                         fixedVVals=4*pi*envObj%gridObj%getVGrid()**(2+momentOrder)*envObj%vSpaceObj%getVCellWidths())

    call multStencilGen%init(genCore,&
                        initXVals=jaggedArray(reshape(real([(1,i=1,usedLocNumX)],kind=rk),[1,usedLocNumX])),&
                        initHVals=jaggedArray(reshape(real([(1,i=1,1)],kind=rk),[1,1])),&
                    initVVals=jaggedArray(reshape(real([(1,i=1,envObj%gridObj%getNumV())],kind=rk),[1,envObj%gridObj%getNumV()])))

    if (interpolationRequired) then 

        if (oddL) then 

            call stencilTemplateObj%defaultStencil%init(xStencil=[-1,0],hStencil=[momentHarmonic],&
                                                       vStencil=[(i,i=1,envObj%gridObj%getNumV())],xPeriodic=pGrid,mapToDist=.true.)
            if (.not. pGrid) then !Handle extrapolation points

                stencilTemplateObj%overridingStencilCoords = allCombinations([IntArray([1,envObj%gridObj%getNumX()])])
                
                allocate(stencilTemplateObj%overridingStencils(size(stencilTemplateObj%overridingStencilCoords,2)))
                do i = 1, size(stencilTemplateObj%overridingStencils)
                    if (stencilTemplateObj%overridingStencilCoords(1,i) == 1) then
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[0,1],hStencil=[momentHarmonic],&
                                                    vStencil=[(i,i=1,envObj%gridObj%getNumV())],xPeriodic=pGrid,mapToDist=.true.)
                    else
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[-2,-1],hStencil=[momentHarmonic],&
                                                    vStencil=[(i,i=1,envObj%gridObj%getNumV())],xPeriodic=pGrid,mapToDist=.true.)
                    end if
                end do

            end if
        else

            call stencilTemplateObj%defaultStencil%init(xStencil=[0,1],hStencil=[momentHarmonic],&
                                                    vStencil=[(i,i=1,envObj%gridObj%getNumV())],xPeriodic=pGrid,mapToDist=.true.)

        end if

        call interpStencilGen%init(envObj%partitionObj,envObj%geometryObj,envObj%mpiCont%getWorldRank(),&
                                   staggeredGridMode=oddL)

        call multStencilGen%setXGen(interpStencilGen)
        
    else

        call stencilTemplateObj%defaultStencil%init(xStencil=[0],hStencil=[momentHarmonic],&
                                                    vStencil=[(i,i=1,envObj%gridObj%getNumV())],xPeriodic=pGrid,mapToDist=.true.)

    end if

    allocate(stencilTemplateObj%stencilGen,source=multStencilGen)


end subroutine initMomentStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initSpatialDiffStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize d/dx kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedInteger) ,dimension(1) :: rowHarmonic, colHarmonic

    rowHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyRowHarmonic,0)
    colHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyColHarmonic,0)

    call envObj%jsonCont%load(rowHarmonic)
    call envObj%jsonCont%output(rowHarmonic)
    call envObj%jsonCont%load(colHarmonic)
    call envObj%jsonCont%output(colHarmonic)

    call initSpatialDiffStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                              rowHarmonic(1)%value,colHarmonic(1)%value)

end subroutine initSpatialDiffStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initSpatialDiffStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,colHarmonic)
    !! Initialize d/dx kinetic stencil template based on direct input. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    integer(ik)                ,intent(in)    :: rowHarmonic
    integer(ik)                ,intent(in)    :: colHarmonic

    integer(ik) :: i  ,spatialStencilCase
    logical     :: pGrid, oddRowL ,oddColL, staggeredRowVar ,staggeredColVar

    integer(ik) ,allocatable ,dimension(:) :: lGrid 

    real(rk) ,allocatable ,dimension(:) :: dx ,linInterp ,outerJ ,innerJ

    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen

    type(UWCDiffStencilValGenerator) :: centralDiffGen
    type(FBDiffStencilValGenerator)  :: staggeredDiffGen

    if (assertions) then 
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)),&
                                "initSpatialDiffStencilTemplateDirect must be invoked with implicitVar being a distribution")
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initSpatialDiffStencilTemplateDirect must be invoked with evolvedVar being a distributione")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions) call assert(staggeredColVar .eqv. staggeredRowVar,&
    "initSpatialDiffStencilTemplateDirect expects both row and column variables to have the same staggered status")

    lGrid = envObj%gridObj%getLGrid()

    oddRowL = mod(lGrid(rowHarmonic),2)==1
    oddColL = mod(lGrid(colHarmonic),2)==1

    stencilTemplateObj%fixedStencil = .true.
    stencilTemplateObj%rowCoords = allCombinations([IntArray([(i,i=1,envObj%gridObj%getNumX())]),&
                                                    IntArray([rowHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])])

    pGrid = envObj%geometryObj%isPeriodic()
    if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) &
    stencilTemplateObj%rowCoords = allCombinations([IntArray([(i,i=1,envObj%gridObj%getNumX()-1)]),&
                                                    IntArray([rowHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])]) !Remove last spatial row in staggered non-periodic case

    allocate(linInterp,source=envObj%geometryObj%getLinInterp(dualGrid=staggeredRowVar))
    dx = envObj%geometryObj%getCellWidths(dualGrid=staggeredRowVar .and. oddRowL,extendedBoundaryCells=.false.)
    
    !Ignore Jacobian data by default
    outerJ = 1/dx
    innerJ = [(real(1,kind=rk),i=1,size(outerJ))]

    allocate(genCore)
    call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                         ,stencilTemplateObj%rowCoords)

    call multStencilGen%init(genCore,&
                        initHVals=jaggedArray(reshape(real([(1,i=1,1)],kind=rk),[1,1])),&
                    initVVals=jaggedArray(reshape(real([(1,i=1,envObj%gridObj%getNumV())],kind=rk),[1,envObj%gridObj%getNumV()])))

    spatialStencilCase = 0 !Central differencing 

    if (staggeredRowVar .and. (oddColL .neqv. oddRowL)) then 
        spatialStencilCase = 1 !Forward diff
        if (oddColL) spatialStencilCase = 2 !Backwards differencing
    end if

    select case (spatialStencilCase) 

    case (0)
        call stencilTemplateObj%defaultStencil%init(xStencil = [-1,0,1],hStencil=[colHarmonic-rowHarmonic]&
                                                   ,xPeriodic=pGrid,mapToDist=.true.)

        call centralDiffGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,linInterp,xPeriodic=pGrid,&
                                                   staggeredGridMode=staggeredRowVar) 

        call multStencilGen%setXGen(centralDiffGen)
        
    case (1)
        call stencilTemplateObj%defaultStencil%init(xStencil = [0,1],hStencil=[colHarmonic-rowHarmonic]&
                                                   ,xPeriodic=pGrid,mapToDist=.true.)

        call staggeredDiffGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,&
                             xPeriodic=pGrid,staggeredGridMode=.true.) 

        call multStencilGen%setXGen(staggeredDiffGen)

    case (2)
        call stencilTemplateObj%defaultStencil%init(xStencil = [-1,0],hStencil=[colHarmonic-rowHarmonic]&
                                                   ,xPeriodic=pGrid,mapToDist=.true.)

        call staggeredDiffGen%init(envObj%partitionObj,envObj%mpiCont%getWorldRank(),innerJ,outerJ,&
                                  xPeriodic=pGrid,staggeredGridMode=.true.,backwardsDiff=.true.) 

        call multStencilGen%setXGen(staggeredDiffGen)

    end select

    allocate(stencilTemplateObj%stencilGen,source=multStencilGen)

end subroutine initSpatialDiffStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDDVStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize d/dv kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedRealArray) ,dimension(1) :: fixedC, fixedInterp ,cfAtZero 
    type(NamedString) ,dimension(1) :: modelboundC ,modelboundInterp
    type(NamedInteger) ,dimension(1) :: rowHarmonic, colHarmonic

    integer(ik) :: i

    fixedC(1) = NamedRealArray(jsonPrefix//"."//keyStencilData//"."//keyFixedC,[(real(1,kind=rk),i=1,envObj%gridObj%getNumV())])
    fixedInterp(1) = NamedRealArray(jsonPrefix//"."//keyStencilData//"."//keyFixedInterp,envObj%vSpaceObj%getVLinInterp())

    cfAtZero(1) = NamedRealArray(jsonPrefix//"."//keyStencilData//"."//keyCFAtZero,real([0,0],kind=rk))

    modelboundC(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyModelboundC,keyNone)
    modelboundInterp(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyModelboundInterp,keyNone)

    rowHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyRowHarmonic,0)
    colHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyColHarmonic,0)

    call envObj%jsonCont%load(rowHarmonic)
    call envObj%jsonCont%output(rowHarmonic)
    call envObj%jsonCont%load(colHarmonic)
    call envObj%jsonCont%output(colHarmonic)
    call envObj%jsonCont%load(fixedC)
    call envObj%jsonCont%output(fixedC)
    call envObj%jsonCont%load(fixedInterp)
    call envObj%jsonCont%output(fixedInterp)
    call envObj%jsonCont%load(cfAtZero)
    call envObj%jsonCont%output(cfAtZero)
    call envObj%jsonCont%load(modelboundC)
    call envObj%jsonCont%output(modelboundC)
    call envObj%jsonCont%load(modelboundInterp)
    call envObj%jsonCont%output(modelboundInterp)

    if (modelboundC(1)%value /= keyNone) then 

        if (modelboundInterp(1)%value /= keyNone) then
            call initDDVStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                              rowHarmonic(1)%value,colHarmonic(1)%value,fixedC=fixedC(1)%values,&
                                              fixedInterp=fixedInterp(1)%values,cfAtZero=cfAtZero(1)%values,&
                                              modelboundC=modelboundC(1)%value,modelboundInterp=modelboundInterp(1)%value)
        else
            call initDDVStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                rowHarmonic(1)%value,colHarmonic(1)%value,fixedC=fixedC(1)%values,&
                                                fixedInterp=fixedInterp(1)%values,cfAtZero=cfAtZero(1)%values,&
                                                modelboundC=modelboundC(1)%value)
        end if 

    else

        if (modelboundInterp(1)%value /= keyNone) then
            call initDDVStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                rowHarmonic(1)%value,colHarmonic(1)%value,fixedC=fixedC(1)%values,&
                                                fixedInterp=fixedInterp(1)%values,cfAtZero=cfAtZero(1)%values,&
                                                modelboundInterp=modelboundInterp(1)%value)
        else
            call initDDVStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                                rowHarmonic(1)%value,colHarmonic(1)%value,fixedC=fixedC(1)%values,&
                                                fixedInterp=fixedInterp(1)%values,cfAtZero=cfAtZero(1)%values)
        end if 

    end if


end subroutine initDDVStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initDDVStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,colHarmonic,&
                                              fixedC,fixedInterp,cfAtZero,modelboundC,modelboundInterp)
    !! Initialize d/dv kinetic stencil template based on direct input. 

    type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)         ,intent(inout) :: envObj
    character(*)                     ,intent(in)    :: evolvedVar
    character(*)                     ,intent(in)    :: implicitVar
    integer(ik)                      ,intent(in)    :: rowHarmonic
    integer(ik)                      ,intent(in)    :: colHarmonic
    real(rk) ,optional ,dimension(:) ,intent(in)    :: fixedC
    real(rk) ,optional ,dimension(:) ,intent(in)    :: fixedInterp
    real(rk) ,optional ,dimension(2) ,intent(in)    :: cfAtZero
    character(*) ,optional           ,intent(in)    :: modelboundC
    character(*) ,optional           ,intent(in)    :: modelboundInterp

    integer(ik) :: i  
    logical     :: pGrid, oddRowL ,oddColL, staggeredRowVar ,staggeredColVar ,interpolationRequired

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,activeXGrid

    real(rk) ,allocatable ,dimension(:) :: usedFixedC ,usedFixedInterp 
    real(rk) ,dimension(2)              :: usedCFAtZero

    type(InterpStencilGenerator) :: interpStencilGen
    type(DDVStencilGenerator) :: ddvStencilGen
    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen

    if (assertions) then 
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)),&
                                "initDDVStencilTemplateDirect must be invoked with implicitVar being a distribution")
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initDDVStencilTemplateDirect must be invoked with evolvedVar being a distributione")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions) call assert(staggeredColVar .eqv. staggeredRowVar,&
    "initDDVStencilTemplateDirect expects both row and column variables to have the same staggered status")

    if (present(fixedC)) then 
        usedFixedC = fixedC
    else
        usedFixedC = [(real(1,kind=rk),i=1,envObj%gridObj%getNumV())]
    end if

    if (present(fixedInterp)) then 
        usedFixedInterp = fixedInterp
    else
        usedFixedInterp = envObj%vSpaceObj%getVLinInterp()
    end if

    if (present(cfAtZero)) then 
        usedCFAtZero = cfAtZero
    else
        usedCFAtZero = real([0,0],kind=rk)
    end if

    if (assertions) then 

        call assert(size(usedFixedC)==envObj%gridObj%getNumV(),&
                    "fixedC passed to initDDVStencilTemplateDirect does not conform with velocity grid size")
        call assert(size(usedFixedInterp)==envObj%gridObj%getNumV(),&
                    "fixedInterp passed to initDDVStencilTemplateDircet does not conform with velocity grid size")

        call assert(size(usedCFAtZero)==2,&
                    "cfAtZero passed to initDDVStencilTemplateDirect should have length 2")

    end if

    pGrid = envObj%geometryObj%isPeriodic()

    lGrid = envObj%gridObj%getLGrid()

    oddRowL = mod(lGrid(rowHarmonic),2)==1
    oddColL = mod(lGrid(colHarmonic),2)==1

    stencilTemplateObj%fixedStencil = .true.

    if (present(modelboundC) .or. present(modelboundInterp)) stencilTemplateObj%fixedStencil = .false.

    interpolationRequired = staggeredRowVar .and. (oddColL .neqv. oddRowL)

    activeXGrid = [(i,i=1,envObj%gridObj%getNumX())]

    if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) &
    activeXGrid = [(i,i=1,envObj%gridObj%getNumX()-1)]

    stencilTemplateObj%rowCoords = allCombinations([IntArray(activeXGrid),&
                                                    IntArray([rowHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])])

    if (present(modelboundC)) then 

        if (present(modelboundInterp)) then
            call ddvStencilGen%init(envObj%partitionObj,envObj%vSpaceObj,envObj%mpiCont%getWorldRank(),activeXGrid,&
                                    fixedC=usedFixedC,fixedInterp=usedFixedInterp,cfAtZero=usedCFAtZero,&
                                    mbC=modelboundC,mbInterp=modelboundInterp)
        else
            call ddvStencilGen%init(envObj%partitionObj,envObj%vSpaceObj,envObj%mpiCont%getWorldRank(),activeXGrid,&
                                    fixedC=usedFixedC,fixedInterp=usedFixedInterp,cfAtZero=usedCFAtZero,&
                                    mbC=modelboundC)
        end if 

    else

        if (present(modelboundInterp)) then
            call ddvStencilGen%init(envObj%partitionObj,envObj%vSpaceObj,envObj%mpiCont%getWorldRank(),activeXGrid,&
                                    fixedC=usedFixedC,fixedInterp=usedFixedInterp,cfAtZero=usedCFAtZero,&
                                    mbInterp=modelboundInterp)
        else
            call ddvStencilGen%init(envObj%partitionObj,envObj%vSpaceObj,envObj%mpiCont%getWorldRank(),activeXGrid,&
            fixedC=usedFixedC,fixedInterp=usedFixedInterp,cfAtZero=usedCFAtZero)
        end if 

    end if

    if (interpolationRequired) then 

        allocate(genCore)
        call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                            ,stencilTemplateObj%rowCoords,vValsDependOnX=.true.)

        call multStencilGen%init(genCore,&
                            initHVals=jaggedArray(reshape(real([(1,i=1,1)],kind=rk),[1,1])))

        if (oddColL) then 

            call stencilTemplateObj%defaultStencil%init(xStencil=[-1,0],hStencil=[colHarmonic-rowHarmonic]&
                                                       ,vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)
            if (.not. pGrid) then !Handle extrapolation points

                stencilTemplateObj%overridingStencilCoords = allCombinations([IntArray([1,envObj%gridObj%getNumX()]),&
                                                                            IntArray([rowHarmonic]),&
                                                                            IntArray([(i,i=1,envObj%gridObj%getNumV())])])
                
                allocate(stencilTemplateObj%overridingStencils(size(stencilTemplateObj%overridingStencilCoords,2)))
                do i = 1, size(stencilTemplateObj%overridingStencils)
                    if (stencilTemplateObj%overridingStencilCoords(1,i) == 1) then
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[0,1],&
                                                    hStencil=[colHarmonic-rowHarmonic],&
                                                    vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)
                    else
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[-2,-1],&
                                                    hStencil=[colHarmonic-rowHarmonic],&
                                                    vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)
                    end if
                end do

            end if
        else

            call stencilTemplateObj%defaultStencil%init(xStencil=[0,1],hStencil=[colHarmonic-rowHarmonic],&
                                                        vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)

        end if

        call interpStencilGen%init(envObj%partitionObj,envObj%geometryObj,envObj%mpiCont%getWorldRank(),&
                                   staggeredGridMode=oddColL)

        call multStencilGen%setXGen(interpStencilGen)
        call multStencilGen%setVGen(ddvStencilGen)

        allocate(stencilTemplateObj%stencilGen,source=multStencilGen)


    else
        call stencilTemplateObj%defaultStencil%init(hStencil=[colHarmonic-rowHarmonic],&
                                                    vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)

        allocate(stencilTemplateObj%stencilGen,source=ddvStencilGen)
        
    end if

end subroutine initDDVStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVelDiffusionStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize d^2/d^2v kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    integer(ik) :: i  

    type(NamedRealArray) ,dimension(1) :: fixedA ,adfAtZero 
    type(NamedString) ,dimension(1) :: modelboundA 
    type(NamedInteger) ,dimension(1) :: rowHarmonic, colHarmonic

    fixeda(1) = NamedRealArray(jsonPrefix//"."//keyStencilData//"."//keyFixedA,[(real(1,kind=rk),i=1,envObj%gridObj%getNumV())])

    adfAtZero(1) = NamedRealArray(jsonPrefix//"."//keyStencilData//"."//keyADFAtZero,real([0,0],kind=rk))

    modelboundA(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyModelboundA,keyNone)

    rowHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyRowHarmonic,0)
    colHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyColHarmonic,0)

    call envObj%jsonCont%load(rowHarmonic)
    call envObj%jsonCont%output(rowHarmonic)
    call envObj%jsonCont%load(colHarmonic)
    call envObj%jsonCont%output(colHarmonic)
    call envObj%jsonCont%load(fixedA)
    call envObj%jsonCont%output(fixedA)
    call envObj%jsonCont%load(adfAtZero)
    call envObj%jsonCont%output(adfAtZero)
    call envObj%jsonCont%load(modelboundA)
    call envObj%jsonCont%output(modelboundA)

    if (modelboundA(1)%value /= keyNone) then 

        call initVelDiffusionStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic(1)%value,&
                                                   colHarmonic(1)%value,fixedA=fixedA(1)%values,adfAtZero=adfAtZero(1)%values,&
                                                   modelboundA=modelboundA(1)%value)


    else

        call initVelDiffusionStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic(1)%value,&
                                                   colHarmonic(1)%value,fixedA=fixedA(1)%values,adfAtZero=adfAtZero(1)%values)

    end if

end subroutine initVelDiffusionStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVelDiffusionStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,colHarmonic,&
                                                        fixedA,adfAtZero,modelboundA)
    !! Initialize d^2/d^2v kinetic stencil template based on direct inputs. 

    type(StencilTemplate)            ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)         ,intent(inout) :: envObj
    character(*)                     ,intent(in)    :: evolvedVar
    character(*)                     ,intent(in)    :: implicitVar
    integer(ik)                      ,intent(in)    :: rowHarmonic
    integer(ik)                      ,intent(in)    :: colHarmonic
    real(rk) ,optional ,dimension(:) ,intent(in)    :: fixedA
    real(rk) ,optional ,dimension(2) ,intent(in)    :: adfAtZero
    character(*) ,optional           ,intent(in)    :: modelboundA

    integer(ik) :: i  
    logical     :: pGrid, oddRowL ,oddColL, staggeredRowVar ,staggeredColVar ,interpolationRequired

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,activeXGrid

    real(rk) ,allocatable ,dimension(:) :: usedFixedA
    real(rk) ,dimension(2)              :: usedADFAtZero

    type(InterpStencilGenerator) :: interpStencilGen
    type(VDiffStencilGen) :: vdStencilGEn
    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen

    if (assertions) then 
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)),&
                                "initVelDiffusionStencilTemplateDirect must be invoked with implicitVar being a distribution")
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initVelDiffusionStencilTemplateDirect must be invoked with evolvedVar being a distributione")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions) call assert(staggeredColVar .eqv. staggeredRowVar,&
    "initVelDiffusionStencilTemplateDirect expects both row and column variables to have the same staggered status")

    if (present(fixedA)) then 
        usedFixedA = fixedA
    else
        usedFixedA = [(real(1,kind=rk),i=1,envObj%gridObj%getNumV())]
    end if 

    if (present(adfAtZero)) then 
        usedADFAtZero = adfAtZero
    else
        usedADFAtZero = real([0,0],kind=rk)
    end if
    
    if (assertions) then 

        call assert(size(usedFixedA)==envObj%gridObj%getNumV(),& 
                    "fixedA passed to initVelDiffusionStencilTemplateDirect does not conform with velocity grid size")

        call assert(size(usedADFAtZero)==2,&
                    "adfAtZero passed to initVelDiffusionStencilTemplateDirect should have length 2")

    end if
    pGrid = envObj%geometryObj%isPeriodic()

    lGrid = envObj%gridObj%getLGrid()

    oddRowL = mod(lGrid(rowHarmonic),2)==1
    oddColL = mod(lGrid(colHarmonic),2)==1

    stencilTemplateObj%fixedStencil = .true.

    if (present(modelboundA)) stencilTemplateObj%fixedStencil = .false.

    interpolationRequired = staggeredRowVar .and. (oddColL .neqv. oddRowL)

    activeXGrid = [(i,i=1,envObj%gridObj%getNumX())]

    if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) &
    activeXGrid = [(i,i=1,envObj%gridObj%getNumX()-1)]

    stencilTemplateObj%rowCoords = allCombinations([IntArray(activeXGrid),&
                                                    IntArray([rowHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])])

    if (present(modelboundA)) then 


        call vdStencilGEn%init(envObj%partitionObj,envObj%vSpaceObj,envObj%mpiCont%getWorldRank(),activeXGrid,&
                                    fixedA=usedFixedA,adfAtZero=usedADFAtZero,&
                                    mbA=modelboundA)


    else

        call vdStencilGEn%init(envObj%partitionObj,envObj%vSpaceObj,envObj%mpiCont%getWorldRank(),activeXGrid,&
                                    fixedA=usedFixedA,adfAtZero=usedADFAtZero)

    end if

    if (interpolationRequired) then 

        allocate(genCore)
        call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                            ,stencilTemplateObj%rowCoords,vValsDependOnX=.true.)

        call multStencilGen%init(genCore,&
                            initHVals=jaggedArray(reshape(real([(1,i=1,1)],kind=rk),[1,1])))

        if (oddColL) then 

            call stencilTemplateObj%defaultStencil%init(xStencil=[-1,0],hStencil=[colHarmonic-rowHarmonic]&
                                                        ,vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)
            if (.not. pGrid) then !Handle extrapolation points

                stencilTemplateObj%overridingStencilCoords = allCombinations([IntArray([1,envObj%gridObj%getNumX()]),&
                                                                            IntArray([rowHarmonic]),&
                                                                            IntArray([(i,i=1,envObj%gridObj%getNumV())])])
                
                allocate(stencilTemplateObj%overridingStencils(size(stencilTemplateObj%overridingStencilCoords,2)))
                do i = 1, size(stencilTemplateObj%overridingStencils)
                    if (stencilTemplateObj%overridingStencilCoords(1,i) == 1) then
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[0,1],&
                                                    hStencil=[colHarmonic-rowHarmonic],&
                                                    vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)
                    else
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[-2,-1],&
                                                    hStencil=[colHarmonic-rowHarmonic],&
                                                    vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)
                    end if
                end do

            end if
        else

            call stencilTemplateObj%defaultStencil%init(xStencil=[0,1],hStencil=[colHarmonic-rowHarmonic],&
                                                        vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)

        end if

        call interpStencilGen%init(envObj%partitionObj,envObj%geometryObj,envObj%mpiCont%getWorldRank(),&
                                   staggeredGridMode=oddColL)

        call multStencilGen%setXGen(interpStencilGen)
        call multStencilGen%setVGen(vdStencilGEn)

        allocate(stencilTemplateObj%stencilGen,source=multStencilGen)


    else
        call stencilTemplateObj%defaultStencil%init(hStencil=[colHarmonic-rowHarmonic]&
                                                    ,vStencil=[-1,0,1],xPeriodic=pGrid,mapToDist=.true.)

        allocate(stencilTemplateObj%stencilGen,source=vdStencilGEn)
        
    end if

end subroutine initVelDiffusionStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initIJStencilTemplate(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize Shkarofsky I/J integral kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedInteger) ,dimension(1) :: rowHarmonic, colHarmonic ,integralIndex
    type(NamedLogical) ,dimension(1) :: jIntegral

    rowHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyRowHarmonic,0)
    colHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyColHarmonic,0)
    integralIndex(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyIntegralIndex,0)

    jIntegral(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyJInt,.false.)

    call envObj%jsonCont%load(rowHarmonic)
    call envObj%jsonCont%output(rowHarmonic)
    call envObj%jsonCont%load(colHarmonic)
    call envObj%jsonCont%output(colHarmonic)
    call envObj%jsonCont%load(integralIndex)
    call envObj%jsonCont%output(integralIndex)
    call envObj%jsonCont%load(jIntegral)
    call envObj%jsonCont%output(jIntegral)

    call initIJStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic(1)%value,colHarmonic(1)%value,&
                                     integralIndex(1)%value,jIntegral(1)%value)

end subroutine initIJStencilTemplate
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initIJStencilTemplateDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                              rowHarmonic,colHarmonic,integralIndex,jIntegral)
    !! Initialize Shkarofsky I/J integral kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    integer(ik)                ,intent(in)    :: rowHarmonic
    integer(ik)                ,intent(in)    :: colHarmonic
    integer(ik)                ,intent(in)    :: integralIndex
    logical                    ,intent(in)    :: jIntegral

    integer(ik) :: i ,usedLocNumX
    logical     :: pGrid, oddRowL ,oddColL, staggeredRowVar ,staggeredColVar ,interpolationRequired

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,activeXGrid ,locNumX

    type(InterpStencilGenerator) :: interpStencilGen
    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen

    type(SparseRowData) :: stencilMat

    if (assertions) then 
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)),&
                                "initIJStencilTemplateDirect must be invoked with implicitVar being a distribution")
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initIJStencilTemplateDirect must be invoked with evolvedVar being a distributione")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    if (assertions) call assert(staggeredColVar .eqv. staggeredRowVar,&
    "initIJStencilTemplateDirect expects both row and column variables to have the same staggered status")

    pGrid = envObj%geometryObj%isPeriodic()

    lGrid = envObj%gridObj%getLGrid()

    oddRowL = mod(lGrid(rowHarmonic),2)==1
    oddColL = mod(lGrid(colHarmonic),2)==1

    stencilTemplateObj%fixedStencil = .true.

    interpolationRequired = staggeredRowVar .and. (oddColL .neqv. oddRowL)

    activeXGrid = [(i,i=1,envObj%gridObj%getNumX())]

    if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) &
    activeXGrid = [(i,i=1,envObj%gridObj%getNumX()-1)]

    stencilTemplateObj%rowCoords = allCombinations([IntArray(activeXGrid),&
                                                    IntArray([rowHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])])

    allocate(genCore)
    call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                        ,stencilTemplateObj%rowCoords)

    if (jIntegral) then
        stencilMat = envObj%vSpaceObj%getShkarofskyJMat(integralIndex)
    else
        stencilMat = envObj%vSpaceObj%getShkarofskyIMat(integralIndex)
    end if

    locNumX = envObj%partitionObj%getLocNumX()
    usedLocNumX = locNumX(envObj%mpiCont%getWorldRank()+1)

    if (envObj%partitionObj%getMaxHAtInd(envObj%mpiCont%getWorldRank()+1)==envObj%gridObj%getNumX()) then 
        if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) usedLocNumX = usedLocNumX - 1
    end if

    call multStencilGen%init(genCore,&
                        initXVals=jaggedArray(reshape(real([(1,i=1,usedLocNumX)],kind=rk),[1,usedLocNumX])),&
                        initHVals=jaggedArray(reshape(real([(1,i=1,1)],kind=rk),[1,1])),&
                        initVVals = stencilMat%values)

    if (interpolationRequired) then 

        if (oddColL) then 
            call stencilTemplateObj%defaultStencil%init(xStencil=[-1,0],hStencil=[colHarmonic-rowHarmonic]&
                                                        ,vStencilFixed=stencilMat%columnVector,xPeriodic=pGrid,mapToDist=.true.)
            if (.not. pGrid) then !Handle extrapolation points
                stencilTemplateObj%overridingStencilCoords = allCombinations([IntArray([1,envObj%gridObj%getNumX()]),&
                                                                            IntArray([rowHarmonic]),&
                                                                            IntArray([(i,i=1,envObj%gridObj%getNumV())])])
                
                allocate(stencilTemplateObj%overridingStencils(size(stencilTemplateObj%overridingStencilCoords,2)))
                do i = 1, size(stencilTemplateObj%overridingStencils)
                    if (stencilTemplateObj%overridingStencilCoords(1,i) == 1) then
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[0,1],&
                                                    hStencil=[colHarmonic-rowHarmonic],&
                                                    vStencilFixed=stencilMat%columnVector,xPeriodic=pGrid,mapToDist=.true.)
                    else
                        call stencilTemplateObj%overridingStencils(i)%init(xStencil=[-2,-1],&
                                                    hStencil=[colHarmonic-rowHarmonic],&
                                                    vStencilFixed=stencilMat%columnVector,xPeriodic=pGrid,mapToDist=.true.)
                    end if
                end do

            end if
        else

            call stencilTemplateObj%defaultStencil%init(xStencil=[0,1],hStencil=[colHarmonic-rowHarmonic],&
                                                        vStencilFixed=stencilMat%columnVector,xPeriodic=pGrid,mapToDist=.true.)

        end if

        call interpStencilGen%init(envObj%partitionObj,envObj%geometryObj,envObj%mpiCont%getWorldRank(),&
                                   staggeredGridMode=oddColL)

        call multStencilGen%setXGen(interpStencilGen)


    else
        call stencilTemplateObj%defaultStencil%init(hStencil=[colHarmonic-rowHarmonic]&
                                                    ,vStencilFixed=stencilMat%columnVector,xPeriodic=pGrid,mapToDist=.true.)
        
    end if

    allocate(stencilTemplateObj%stencilGen,source=multStencilGen)


end subroutine initIJStencilTemplateDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initFixedBoltzmannStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
    !! Initialize fixed mapping Boltzmann kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    class(ModelboundCRMData)   ,intent(in)    :: mbData

    type(NamedInteger) ,dimension(1) :: rowHarmonic ,transitionIndex, fixedEnergyIndex
    type(NamedLogical) ,dimension(1) :: absorptionTerm ,dbTerm

    rowHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyRowHarmonic,0)
    transitionIndex(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyTransitionIndex,0)
    fixedEnergyIndex(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyFixedEnergyIndex,0)

    absorptionTerm(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyAbsorptionTerm,.false.)
    dbTerm(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyDBTerm,.false.)

    call envObj%jsonCont%load(rowHarmonic)
    call envObj%jsonCont%output(rowHarmonic)
    call envObj%jsonCont%load(transitionIndex)
    call envObj%jsonCont%output(transitionIndex)
    call envObj%jsonCont%load(fixedEnergyIndex)
    call envObj%jsonCont%output(fixedEnergyIndex)
    call envObj%jsonCont%load(absorptionTerm)
    call envObj%jsonCont%output(absorptionTerm)
    call envObj%jsonCont%load(dbTerm)
    call envObj%jsonCont%output(dbTerm)

    call initFixedBoltzmannStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,mbData,rowHarmonic(1)%value,&
                                         transitionIndex(1)%value,fixedEnergyIndex(1)%value,absorptionTerm(1)%value,dbTerm(1)%value)

end subroutine initFixedBoltzmannStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initFixedBoltzmannStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,mbData,&
                                           evolvedHarmonic,transitionIndex,fixedEnergyIndex,absorptionTerm,dbTerm)
    !! Initialize fixed mapping Boltzmann kinetic stencil template based on environment non-JSON inputs. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    class(ModelboundCRMData)   ,intent(in)    :: mbData
    integer(ik)                ,intent(in)    :: evolvedHarmonic
    integer(ik)                ,intent(in)    :: transitionIndex
    integer(ik)                ,intent(in)    :: fixedEnergyIndex
    logical                    ,intent(in)    :: absorptionTerm
    logical                    ,intent(in)    :: dbTerm

    integer(ik) :: i ,usedLocNumX
    logical     :: pGrid, oddRowL, staggeredRowVar

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,activeXGrid ,locNumX

    type(FixedBoltzmannStencilGen) :: boltzStencilGen
    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen
    type(IntArray) ,allocatable ,dimension(:) :: fixedVStencil

    type(SparseRowData) :: stencilMat

    if (assertions) then 
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)),&
                                "initFixedBoltzmannStencilDirect must be invoked with implicitVar being a distribution")
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initFixedBoltzmannStencilDirect must be invoked with evolvedVar being a distributione")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)

    pGrid = envObj%geometryObj%isPeriodic()
    
    lGrid = envObj%gridObj%getLGrid()

    oddRowL = mod(lGrid(evolvedHarmonic),2)==1

    stencilTemplateObj%fixedStencil = .not. dbTerm

    activeXGrid = [(i,i=1,envObj%gridObj%getNumX())]

    if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) &
    activeXGrid = [(i,i=1,envObj%gridObj%getNumX()-1)]

    stencilTemplateObj%rowCoords = allCombinations([IntArray(activeXGrid),&
                                                    IntArray([evolvedHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])])

    allocate(genCore)
    call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                        ,stencilTemplateObj%rowCoords,vValsDependOnX=dbTerm)

    if (absorptionTerm) then
        stencilMat = mbData%getFixedW(fixedEnergyIndex)
        allocate(fixedVStencil(envObj%gridObj%getNumV()))
        do i = 1,envObj%gridObj%getNumV()
            allocate(fixedVStencil(i)%entry(0))
        end do

        do i = 1,size(stencilMat%rowIndex)
            fixedVStencil(stencilMat%rowIndex(i))%entry = stencilMat%columnVector(i)%entry
        end do
        call stencilTemplateObj%defaultStencil%init(vStencilFixed=fixedVStencil,mapToDist=.true.)  
    else
        call stencilTemplateObj%defaultStencil%init(mapToDist=.true.)  
    end if

    locNumX = envObj%partitionObj%getLocNumX()
    usedLocNumX = locNumX(envObj%mpiCont%getWorldRank()+1)

    if (envObj%partitionObj%getMaxHAtInd(envObj%mpiCont%getWorldRank()+1)==envObj%gridObj%getNumX()) then 
        if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) usedLocNumX = usedLocNumX - 1
    end if

    call multStencilGen%init(genCore,&
                        initXVals=jaggedArray(reshape(real([(1,i=1,usedLocNumX)],kind=rk),[1,usedLocNumX])),&
                        initHVals=jaggedArray(reshape(real([(1,i=1,1)],kind=rk),[1,1])))

    call boltzStencilGen%init(envObj%vSpaceObj,transitionIndex,fixedEnergyIndex,&
                            lNum=lGrid(evolvedHarmonic),absorptionTerm=absorptionTerm,dbTerm=dbTerm)

    call multStencilGen%setVGen(boltzStencilGen)
    allocate(stencilTemplateObj%stencilGen,source=multStencilGen)

end subroutine initFixedBoltzmannStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initScalingLBCStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize scaling logical boundary condition kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    integer(ik) :: i 

    type(NamedInteger) ,dimension(1) :: rowHarmonic ,colHarmonic
    type(NamedLogical) ,dimension(1) :: leftBoundary
    type(NamedIntegerArray) ,dimension(1) :: decompHarmonics
    type(NamedStringArray) ,dimension(1) :: derivReqVars
    type(NamedString) ,dimension(1) :: derivName

    rowHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyRowHarmonic,0)
    colHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyColHarmonic,1)

    leftBoundary(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyLeftBoundary,.false.)
    derivName(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyRuleName,"")
    decompHarmonics(1) = NamedIntegerArray(jsonPrefix//"."//keyStencilData//"."//keyDecompHarmonics,&
                        [(i,i=1,envObj%gridObj%getNumH())])
    derivReqVars(1)%name = jsonPrefix//"."//keyStencilData//"."//keyReqVarNames
    allocate(derivReqVars(1)%values(0))

    call envObj%jsonCont%load(rowHarmonic)
    call envObj%jsonCont%output(rowHarmonic)
    call envObj%jsonCont%load(colHarmonic)
    call envObj%jsonCont%output(colHarmonic)
    call envObj%jsonCont%load(leftBoundary)
    call envObj%jsonCont%output(leftBoundary)
    call envObj%jsonCont%load(decompHarmonics)
    call envObj%jsonCont%output(decompHarmonics)
    call envObj%jsonCont%load(derivReqVars)
    call envObj%jsonCont%output(derivReqVars)
    call envObj%jsonCont%load(derivName)
    call envObj%jsonCont%output(derivName)

    call initScalingLBCStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic(1)%value,colHarmonic(1)%value,&
                                    decompHarmonics(1)%values,derivName(1)%value,derivReqVars(1)%values,leftBoundary(1)%value)
    
end subroutine initScalingLBCStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initScalingLBCStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,rowHarmonic,colHarmonic,&
                                              decompHarmonics,derivName,derivReqVars,leftBoundary)
    !! Initialize scaling logical boundary condition kinetic stencil template based on direct input. 

    type(StencilTemplate)           ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)        ,intent(inout) :: envObj
    character(*)                    ,intent(in)    :: evolvedVar
    character(*)                    ,intent(in)    :: implicitVar
    integer(ik)                     ,intent(in)    :: rowHarmonic
    integer(ik)                     ,intent(in)    :: colHarmonic
    integer(ik) ,dimension(:)       ,intent(in)    :: decompHarmonics
    character(*)                    ,intent(in)    :: derivName
    type(StringArray) ,dimension(:) ,intent(in)    :: derivReqVars
    logical ,optional               ,intent(in)    :: leftBoundary

    integer(ik) :: i ,xCoord ,xStencil
    logical     :: oddRowL, staggeredRowVar ,isActive ,lBoundary

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,derivReqIndices

    real(rk) ,allocatable ,dimension(:) :: dx 

    type(ScalingLBCStencilGen) :: stencilGen

    class(MatDerivation) ,allocatable :: derivObj

    if (assertions) &
    call assert(evolvedVar==implicitVar,"initScalingLBCStencilDirect requires that the evolved and implicit variables be the same")

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)

    lBoundary = .false. 
    if (present(leftBoundary)) lBoundary = leftBoundary

    allocate(derivReqIndices(size(derivReqVars)))
    do i = 1,size(derivReqIndices)
        derivReqIndices(i) = envObj%externalVars%getVarIndex(derivReqVars(i)%string)
    end do

    lGrid = envObj%gridObj%getLGrid()

    oddRowL = mod(lGrid(rowHarmonic),2)==1

    call envObj%textbookObj%copyMatDerivation(derivName,derivObj)

    select type (derivObj)
    type is (FScalingDerivation)

        if (staggeredRowVar) then
            call assert(.not. oddRowL,"initScalingLBCStencilDirect cannot have odd row l numbers on a staggered grid")
            if (.not. lBoundary) &
            call assert(all(mod(lGrid(decompHarmonics),2)==1) .or. all(mod(lGrid(decompHarmonics),2)==0),&
            "If the right boundary is treated using initScalingLBCStencilDirect on a staggered grid, included decomposition &
            &harmonics must either all be even or all odd")
        end if

        xCoord = envObj%gridObj%getNumX()
        if (lBoundary) xCoord = 1
        stencilTemplateObj%rowCoords = allCombinations([IntArray([xCoord]),&
                                                    IntArray([rowHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])])

        isActive = size(envObj%partitionObj%filterCoords(envObj%mpiCont%getWorldRank()+1,stencilTemplateObj%rowCoords),2) > 0

        xStencil = 0
        if (.not. lBoundary .and. all(mod(lGrid(decompHarmonics),2)==1)) xStencil = -1 
        call stencilTemplateObj%defaultStencil%init(xStencil=[xStencil], &
                                                    hStencil=decompHarmonics-rowHarmonic,& 
                                                    vStencil=[0,1],mapToDist=.true.)

        dx = envObj%geometryObj%getCellWidths()

        call stencilGen%init(envObj%vSpaceObj,isActive,derivObj,derivReqIndices,lGrid(colHarmonic),&
                                dx(xCoord),includedDecompHarmonics=decompHarmonics)

        allocate(stencilTemplateObj%stencilGen,source=stencilGen)
    class default
        error stop &
        "initScalingLBCStencilDirect requires passed matrix derivation &
        &name to correspond to a derivation of type FScalingDerivation"
    end select
    
end subroutine initScalingLBCStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initTermMomentStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar)
    !! Initialize term moment kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar

    type(NamedString)  ,dimension(1) :: termName
    type(NamedInteger) ,dimension(1) :: momentOrder ,colHarmonic
    
    termName(1) = NamedString(jsonPrefix//"."//keyStencilData//"."//keyTermName,"")
    momentOrder(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyMomentOrder,0)
    colHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyColHarmonic,1)

    call envObj%jsonCont%load(colHarmonic)
    call envObj%jsonCont%output(colHarmonic)
    call envObj%jsonCont%load(momentOrder)
    call envObj%jsonCont%output(momentOrder)
    call envObj%jsonCont%load(termName)
    call envObj%jsonCont%output(termName)

    call initTermMomentStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                     momentOrder(1)%value,colHarmonic(1)%value,termName(1)%value)

end subroutine initTermMomentStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initTermMomentStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,&
                                              momentOrder,colHarmonic,termName)
    !! Initialize term moment kinetic stencil template based on direct input.

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    integer(ik)                ,intent(in)    :: momentOrder
    integer(ik)                ,intent(in)    :: colHarmonic
    character(*)               ,intent(in)    :: termName

    type(TermMomentStencilGenerator) :: stencilGen

    integer(ik) :: i ,usedLocNumX
    logical     :: pGrid, staggeredRowVar ,staggeredColVar ,oddL

    integer(ik) ,allocatable ,dimension(:) :: locNumX ,xStencil ,lGrid

    if (assertions) then 
        call assert(.not. envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initTermMomentStencilDirect must be invoked with evolvedVar being a fluid variable")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)
    staggeredColVar = envObj%externalVars%isVarOnDualGrid(implicitVar)

    stencilTemplateObj%rowCoords = allCombinations([IntArray([(i,i=1,envObj%gridObj%getNumX())])])

    pGrid = envObj%geometryObj%isPeriodic()
    if (staggeredRowVar .and. (.not. pGrid)) &
    stencilTemplateObj%rowCoords = allCombinations([IntArray([(i,i=1,envObj%gridObj%getNumX()-1)])]) !Remove last spatial row in staggered non-periodic case
    xStencil = [0]

    if (.not. envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar))) then
        if (staggeredRowVar .neqv. staggeredColVar) then
            xStencil = [-1,0]
            if (staggeredRowVar) xStencil = [0,1]
        end if
    else 
        lGrid = envObj%gridObj%getLGrid()
    
        oddL = mod(lGrid(colHarmonic),2) == 1
        if (assertions) call assert(oddL .eqv. staggeredRowVar,"If initTermMomentStencilDirect is called with an implcit &
        &distribution variable the row variable must be staggered if the implicit harmonic has odd l")
    end if
    call stencilTemplateObj%defaultStencil%init(xStencil=xStencil,hStencil=[colHarmonic],&
                                                    vStencil=[(i,i=1,envObj%gridObj%getNumV())],xPeriodic=pGrid,&
                                            mapToDist=envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)))

    if ((.not. pGrid) .and. size(xStencil) > 1) then !Handle extrapolation points
        stencilTemplateObj%overridingStencilCoords = allCombinations([IntArray([1,envObj%gridObj%getNumX()])])
        allocate(stencilTemplateObj%overridingStencils(size(stencilTemplateObj%overridingStencilCoords,2)))
        do i = 1, size(stencilTemplateObj%overridingStencils)
            if (stencilTemplateObj%overridingStencilCoords(1,i) == 1) then
                call stencilTemplateObj%overridingStencils(i)%init(xStencil=[0,1],xPeriodic=pGrid)
            else
                call stencilTemplateObj%overridingStencils(i)%init(xStencil=[-2,-1],xPeriodic=pGrid)

            end if
        end do
    end if
    call stencilGen%init(envObj%partitionObj,envObj%vSpaceObj,envObj%mpiCont%getWorldRank(),&
                         momentOrder,termName,removeLastCell=staggeredRowVar .and. (.not. pGrid))
    allocate(stencilTemplateObj%stencilGen,source=stencilGen)

end subroutine initTermMomentStencilDirect
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVariableBoltzmannStencil(stencilTemplateObj,envObj,jsonPrefix,evolvedVar,implicitVar,mbData)
    !! Initialize variable mapping Boltzmann kinetic stencil template based on environment object and JSON file. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: jsonPrefix !! prefix for JSON keys for this template
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    class(ModelboundCRMData)   ,intent(in)    :: mbData
    
    type(NamedInteger) ,dimension(1) :: rowHarmonic ,transitionIndex
    type(NamedLogical) ,dimension(1) :: absorptionTerm ,superelasticTerm

    rowHarmonic(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyRowHarmonic,0)
    transitionIndex(1) = NamedInteger(jsonPrefix//"."//keyStencilData//"."//keyTransitionIndex,0)

    absorptionTerm(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keyAbsorptionTerm,.false.)
    superelasticTerm(1) = NamedLogical(jsonPrefix//"."//keyStencilData//"."//keySuperelasticTerm,.false.)

    call envObj%jsonCont%load(rowHarmonic)
    call envObj%jsonCont%output(rowHarmonic)
    call envObj%jsonCont%load(transitionIndex)
    call envObj%jsonCont%output(transitionIndex)
    call envObj%jsonCont%load(absorptionTerm)
    call envObj%jsonCont%output(absorptionTerm)
    call envObj%jsonCont%load(superelasticTerm)
    call envObj%jsonCont%output(superelasticTerm)

    call initVariableBoltzmannStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,mbData,rowHarmonic(1)%value,&
                                         transitionIndex(1)%value,absorptionTerm(1)%value,superelasticTerm(1)%value)

end subroutine initVariableBoltzmannStencil
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initVariableBoltzmannStencilDirect(stencilTemplateObj,envObj,evolvedVar,implicitVar,mbData,&
    evolvedHarmonic,transitionIndex,absorptionTerm,superelasticTerm)
    !! Initialize variable mapping Boltzmann kinetic stencil template based on environment non-JSON inputs. 

    type(StencilTemplate)      ,intent(inout) :: stencilTemplateObj
    type(EnvironmentWrapper)   ,intent(inout) :: envObj
    character(*)               ,intent(in)    :: evolvedVar
    character(*)               ,intent(in)    :: implicitVar
    class(ModelboundCRMData)   ,intent(in)    :: mbData
    integer(ik)                ,intent(in)    :: evolvedHarmonic
    integer(ik)                ,intent(in)    :: transitionIndex
    logical                    ,intent(in)    :: absorptionTerm
    logical                    ,intent(in)    :: superelasticTerm

    integer(ik) :: i ,usedLocNumX
    logical     :: pGrid, oddRowL, staggeredRowVar

    integer(ik) ,allocatable ,dimension(:) :: lGrid ,activeXGrid ,locNumX

    type(VariableBoltzmannStencilGen) :: boltzStencilGen
    type(MultiplicativeGeneratorCore) ,allocatable :: genCore 
    type(MultiplicativeStencilGen) :: multStencilGen
    type(IntArray) ,allocatable ,dimension(:) :: fixedVStencil

    if (assertions) then 
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(implicitVar)),&
                                "initVariableBoltzmannStencilDirect must be invoked with implicitVar being a distribution")
        call assert(envObj%externalVars%isVarDist(envObj%externalVars%getVarIndex(evolvedVar)),&
                                "initVariableBoltzmannStencilDirect must be invoked with evolvedVar being a distributione")
    end if

    staggeredRowVar = envObj%externalVars%isVarOnDualGrid(evolvedVar)

    pGrid = envObj%geometryObj%isPeriodic()
    
    lGrid = envObj%gridObj%getLGrid()

    oddRowL = mod(lGrid(evolvedHarmonic),2)==1

    activeXGrid = [(i,i=1,envObj%gridObj%getNumX())]

    if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) &
    activeXGrid = [(i,i=1,envObj%gridObj%getNumX()-1)]

    stencilTemplateObj%rowCoords = allCombinations([IntArray(activeXGrid),&
                                                    IntArray([evolvedHarmonic]),&
                                                    IntArray([(i,i=1,envObj%gridObj%getNumV())])])

    allocate(genCore)
    call genCore%init(envObj%gridObj,envObj%partitionObj,envObj%mpiCont%getWorldRank()&
                        ,stencilTemplateObj%rowCoords,vValsDependOnX=.true.)

    if (absorptionTerm) then
        fixedVStencil = triangularIntArray(envObj%gridObj%getNumV(),lower=superelasticTerm)
        call stencilTemplateObj%defaultStencil%init(vStencilFixed=fixedVStencil,mapToDist=.true.)  
    else
        call stencilTemplateObj%defaultStencil%init(mapToDist=.true.)  
    end if

    locNumX = envObj%partitionObj%getLocNumX()
    usedLocNumX = locNumX(envObj%mpiCont%getWorldRank()+1)

    if (envObj%partitionObj%getMaxHAtInd(envObj%mpiCont%getWorldRank()+1)==envObj%gridObj%getNumX()) then 
        if (staggeredRowVar .and. (.not. pGrid) .and. oddRowL) usedLocNumX = usedLocNumX - 1
    end if

    call multStencilGen%init(genCore,&
                        initXVals=jaggedArray(reshape(real([(1,i=1,usedLocNumX)],kind=rk),[1,usedLocNumX])),&
                        initHVals=jaggedArray(reshape(real([(1,i=1,1)],kind=rk),[1,1])))

    call boltzStencilGen%init(envObj%vSpaceObj,transitionIndex,&
                            lNum=lGrid(evolvedHarmonic),absorptionTerm=absorptionTerm,superelasticTerm=superelasticTerm)

    call multStencilGen%setVGen(boltzStencilGen)
    allocate(stencilTemplateObj%stencilGen,source=multStencilGen)

end subroutine initVariableBoltzmannStencilDirect 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule kinetic_stencil_templates_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
