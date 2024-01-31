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
submodule (custom_derivation_support) custom_derivation_support_procedures
!! author: Stefan Mijin 
!! 
!! Contains custom derivation support procedure implementations

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addCustomDerivationsToTextbook(textbookObj,gridObj,geometryObj,partObj,&
    vspaceObj,normObj,speciesListObj,varList,jsonCont,mpiCont)
    !! Adds custom derivations to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    type(Grid)              ,intent(in)    :: gridObj   
    type(Geometry)          ,intent(in)    :: geometryObj 
    type(Partition)         ,intent(in)    :: partObj
    type(VSpace)            ,intent(in)    :: vSpaceObj
    class(Normalization)    ,intent(in)    :: normObj 
    type(SpeciesList)       ,intent(in)    :: speciesListObj  
    type(VariableList)      ,intent(in)    :: varList
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedStringArray) ,dimension(1) :: customDerivationTags 
    type(NamedString) ,dimension(1) :: derivType

    integer(ik) :: i

    character(:) ,allocatable :: derivTag,derivTypeName 

    customDerivationTags(1)%name = keyCustomDerivations//"."//keyTags
    allocate(customDerivationTags(1)%values(0))

    call jsonCont%load(customDerivationTags)

    call jsonCont%output(customDerivationTags)

    do i = 1,size(customDerivationTags(1)%values)
        derivTag = customDerivationTags(1)%values(i)%string
        derivType(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyType,"")
        
        call jsonCont%load(derivType)
        call jsonCont%output(derivType)

        derivTypeName = derivType(1)%value

        select case(derivTypeName)
        case (keySimpleDeriv)
            call addSimpleDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        case (keyPolyFunDeriv)
            call addPolyDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        case (keyAdditiveDeriv)
            call addAdditiveDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        case (keyMultDeriv)
            call addMultiplicativeDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        case (keyBoundedExtDeriv)
            call addBoundedExtDerivationToTextbook(textbookObj,derivTag,gridObj,partObj,geometryObj,jsonCont,mpiCont)
        case (keyCIIJDeriv)
            call addColdIonIJIntDerivationToTextbook(textbookObj,derivTag,gridObj,jsonCont,mpiCont)
        case (keyIJDeriv)
            call addIJIntDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        case (keyHarmonicExtractor)
            call addHarmonicExtractorDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        case (keyFScalingExtrapolation)
            call addFScalingDerivationToTextbook(textbookObj,derivTag,vSpaceObj,partObj,jsonCont,mpiCont)
        case (keyDDVDeriv)
            call addDDVDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        case (keyD2DV2Deriv)
            call addD2DV2DerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        case (keyMomentDerivation)
            call addMomentDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        case (keyLocValExtractor)
            call addLocalValExtractorDerivationToTextbook(textbookObj,derivTag,partObj,jsonCont,mpiCont)
        case (keyVelContracDeriv)
            call addVelContracDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        case (keyVelTProdDeriv)
            call addVelTProdDeivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        case (keyGenIntPolyFunDeriv)
            call addGenIntPolyDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        case (keyRangeFilterDerivation)
            call addRangeFilterDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        case (keyCalculationTreeDerivation)
            call addCalculationTreeDerivationToTextbook(textbookObj,derivTag,varList,jsonCont,mpiCont)
        case (keynDLinInterpDerivation)
            call addnDLinInterpDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        case default 
            error stop "Unsupported custom derivation type detected in addCustomDerivationsToTextbook"
        end select
    end do


end subroutine addCustomDerivationsToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addSimpleDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
    !! Add simple derivation to textbook based on JSON data

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedReal) ,dimension(1) :: multConst
    type(NamedRealArray) ,dimension(1) :: varPowers

    type(SimpleDerivation) :: simpleDeriv

    multConst(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyMultConst,real(1,kind=rk))
    varPowers(1)%name = keyCustomDerivations//"."//derivTag//"."//keyVarPowers
    allocate(varPowers(1)%values(0))

    call jsonCont%load(multConst)
    call jsonCont%load(varPowers)

    call jsonCont%output(multConst)
    call jsonCont%output(varPowers)

    call simpleDeriv%init(varPowers(1)%values,multConst(1)%value)

    call textbookObj%addDerivation(simpleDeriv,derivTag)

end subroutine addSimpleDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addPolyDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
    !! Add a polynomial derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedReal) ,dimension(1) :: constCoeff
    type(NamedRealArray) ,dimension(1) :: polyPowers, polyCoeffs

    type(PolyFunDeriv) :: polyDeriv 

    constCoeff(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyConstCoeff,real(0,kind=rk))
    polyPowers(1)%name = keyCustomDerivations//"."//derivTag//"."//keyPolyPowers
    allocate(polyPowers(1)%values(0))

    polyCoeffs(1)%name = keyCustomDerivations//"."//derivTag//"."//keyPolyCoeffs
    allocate(polyCoeffs(1)%values(0))

    call jsonCont%load(constCoeff)
    call jsonCont%load(polyCoeffs)
    call jsonCont%load(polyPowers)

    call jsonCont%output(constCoeff)
    call jsonCont%output(polyCoeffs)
    call jsonCont%output(polyPowers)

    call polyDeriv%init(polyPowers(1)%values,polyCoeffs(1)%values,constCoeff(1)%value)

    call textbookObj%addDerivation(polyDeriv,derivTag)

end subroutine addPolyDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addAdditiveDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
    !! Add an additive derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedReal) ,dimension(1) :: resultPower
    type(NamedStringArray) ,dimension(1) :: derivTags 
    type(NamedIntegerArray) ,allocatable ,dimension(:) :: derivIndices 
    type(NamedRealArray) ,dimension(1) :: linCoeffs

    integer(ik) :: i ,maxIndex

    type(AdditiveDerivation) :: additiveDeriv

    class(Derivation) ,allocatable :: derivBuffer

    resultPower(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyResultPower,real(1,kind=rk))

    call jsonCont%load(resultPower)
    call jsonCont%output(resultPower)

    derivTags(1)%name = keyCustomDerivations//"."//derivTag//"."//keyDerivTags
    allocate(derivTags(1)%values(0))

    call jsonCont%load(derivTags)
    call jsonCont%output(derivTags)

    linCoeffs(1)%name = keyCustomDerivations//"."//derivTag//"."//keyLinearCoefficients
    allocate(linCoeffs(1)%values(0))

    call jsonCont%load(linCoeffs)
    call jsonCont%output(linCoeffs)

    allocate(derivIndices(size(derivTags(1)%values)))
    do i = 1,size(derivTags(1)%values)
        derivIndices(i)%name = &
        keyCustomDerivations//"."//derivTag//"."//derivTags(1)%values(i)%string//"."//keyDerivIndices
        allocate(derivIndices(i)%values(0))
    end do

    call jsonCont%load(derivIndices)
    call jsonCont%output(derivIndices)

    maxIndex = 0

    do i = 1,size(derivIndices)
        maxIndex = max(maxIndex,maxval(derivIndices(i)%values))
    end do

    if (size(linCoeffs(1)%values) > 0) then
        call assert(size(linCoeffs(1)%values) == size(derivTags(1)%values),linCoeffs(1)%name//" must have same size as "&
        // derivTags(1)%name)
    else
        linCoeffs(1)%values = [(real(1,kind=rk),i=1,size(derivTags(1)%values))]
    end if

    call additiveDeriv%init(size(derivIndices),maxIndex,resultPower(1)%value,linearCoefficients=linCoeffs(1)%values)

    do i = 1,size(derivTags(1)%values)
        call textbookObj%copyDerivation(derivTags(1)%values(i)%string,derivBuffer)
        call additiveDeriv%addDerivation(derivBuffer,derivIndices(i)%values)
    end do

    call textbookObj%addDerivation(additiveDeriv,derivTag)

end subroutine addAdditiveDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addMultiplicativeDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
    !! Add a multiplicative derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedReal) ,dimension(1) :: innerDerivPower ,outerDerivPower
    type(NamedString) ,dimension(1) :: innerDerivation, outerDerivation ,funcName
    type(NamedIntegerArray) ,dimension(1) :: innerDerivIndices ,outerDerivIndices

    class(Derivation) ,allocatable :: derivBuffer, outerDerivBuffer

    type(MultiplicativeDerivation) :: multiplicativeDeriv

    innerDerivPower(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyInnerDerivPower,real(1,kind=rk))
    outerDerivPower(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyOuterDerivPower,real(1,kind=rk))

    call jsonCont%load(innerDerivPower)
    call jsonCont%output(innerDerivPower)
    call jsonCont%load(outerDerivPower)
    call jsonCont%output(outerDerivPower)

    innerDerivation(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyInnerDerivation,"")
    outerDerivation(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyOuterDerivation,keyNone)
    funcName(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyInnerDerivFuncName,keyNone)

    call jsonCont%load(innerDerivation)
    call jsonCont%output(innerDerivation)
    call jsonCont%load(outerDerivation)
    call jsonCont%output(outerDerivation)
    call jsonCont%load(funcName)
    call jsonCont%output(funcName)
    
    innerDerivIndices(1)%name = keyCustomDerivations//"."//derivTag//"."//keyInnerDerivIndices
    allocate(innerDerivIndices(1)%values(0))

    outerDerivIndices(1)%name = keyCustomDerivations//"."//derivTag//"."//keyOuterDerivIndices
    allocate(outerDerivIndices(1)%values(0))

    call jsonCont%load(innerDerivIndices)
    call jsonCont%output(innerDerivIndices)
    call jsonCont%load(outerDerivIndices)
    call jsonCont%output(outerDerivIndices)

    call textbookObj%copyDerivation(innerDerivation(1)%value,derivBuffer)
    if (outerDerivation(1)%value /= keyNone) then  
        call textbookObj%copyDerivation(outerDerivation(1)%value,outerDerivBuffer)

        if (funcName(1)%value /=keyNone) then 

            call multiplicativeDeriv%init(derivBuffer,innerDerivIndices(1)%values,outerDeriv=outerDerivBuffer,&
                                            outerIndices=outerDerivIndices(1)%values,&
                                            innerPower=innerDerivPower(1)%value,outerPower=outerDerivPower(1)%value,&
                                            innerFuncName=funcName(1)%value)
        else

            call multiplicativeDeriv%init(derivBuffer,innerDerivIndices(1)%values,outerDeriv=outerDerivBuffer,&
                                            outerIndices=outerDerivIndices(1)%values,&
                                            innerPower=innerDerivPower(1)%value,outerPower=outerDerivPower(1)%value)

        end if

    else

        if (funcName(1)%value /=keyNone) then 

            call multiplicativeDeriv%init(derivBuffer,innerDerivIndices(1)%values,&
                                            innerPower=innerDerivPower(1)%value,&
                                            innerFuncName=funcName(1)%value)
        else

            call multiplicativeDeriv%init(derivBuffer,innerDerivIndices(1)%values,&
                                            innerPower=innerDerivPower(1)%value)

        end if

    end if

    call textbookObj%addDerivation(multiplicativeDeriv,derivTag)

end subroutine addMultiplicativeDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addBoundedExtDerivationToTextbook(textbookObj,derivTag,gridObj,partObj,geometryObj,jsonCont,mpiCont)
    !! Add a bounded extrapolation derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(Grid)              ,intent(in)    :: gridObj   
    type(Partition)         ,intent(in)    :: partObj
    type(Geometry)          ,intent(in)    :: geometryObj 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedLogical) ,dimension(1) :: expectUpperBoundVar ,expectLowerBoundVar &
                                       ,ignoreUpperBound, ignoreLowerBound
    type(NamedReal) ,dimension(1) :: fixedLowerBound ,fixedUpperBound

    real(rk) ,allocatable ,dimension(:) :: dx ,linInterp 
    real(rk)                            :: lInterp, lExterp
    type(BoundedExtDerivation)          :: boundedExtDeriv
    class(Extrapolation) ,allocatable   :: extrapObj

    expectUpperBoundVar(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyExpectUpperBoundVar,.false.)
    expectLowerBoundVar(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyExpectLowerBoundVar,.false.)
    ignoreUpperBound(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyIgnoreUpperBound,.false.)
    ignoreLowerBound(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyIgnoreLowerBound,.false.)

    fixedLowerBound(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyFixedLowerBound,real(0,kind=rk))
    fixedUpperBound(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyFixedUpperBound,huge(real(0,kind=rk)))

    call jsonCont%load(expectLowerBoundVar)
    call jsonCont%output(expectLowerBoundVar)
    call jsonCont%load(expectUpperBoundVar)
    call jsonCont%output(expectUpperBoundVar)
    call jsonCont%load(ignoreLowerBound)
    call jsonCont%output(ignoreLowerBound)
    call jsonCont%load(ignoreUpperBound)

    call jsonCont%load(fixedLowerBound)
    call jsonCont%output(fixedLowerBound)

    call jsonCont%load(fixedUpperBound)
    call jsonCont%output(fixedUpperBound)
    
    call initExtrapolationFromJSON(extrapObj,keyCustomDerivations//"."//derivTag,partObj,gridObj,geometryObj,jsonCont,mpiCont)

    if (ignoreLowerBound(1)%value) then 
        if (ignoreUpperBound(1)%value) then 
            call boundedExtDeriv%init(partObj,mpiCont%getWorldRank(),extrapObj)
        else
            call boundedExtDeriv%init(partObj,mpiCont%getWorldRank(),extrapObj,fixedUpperBound=fixedUpperBound(1)%value,&
                                    expectUpperBoundVar=expectUpperBoundVar(1)%value)
        end if

    else

        if (ignoreUpperBound(1)%value) then 
            call boundedExtDeriv%init(partObj,mpiCont%getWorldRank(),extrapObj,fixedLowerBound=fixedLowerBound(1)%value,&
            expectLowerBoundVar=expectLowerBoundVar(1)%value)
        else
            call boundedExtDeriv%init(partObj,mpiCont%getWorldRank(),extrapObj,fixedUpperBound=fixedUpperBound(1)%value,&
            expectUpperBoundVar=expectUpperBoundVar(1)%value,fixedLowerBound=fixedLowerBound(1)%value,&
            expectLowerBoundVar=expectLowerBoundVar(1)%value)
        end if

    end if

    call textbookObj%addDerivation(boundedExtDeriv,derivTag)

end subroutine addBoundedExtDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addColdIonIJIntDerivationToTextbook(textbookObj,derivTag,gridObj,jsonCont,mpiCont)
    !! Add a cold ion I/J integral derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(Grid)              ,intent(in)    :: gridObj   
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedLogical) ,dimension(1) :: isJInt
    type(NamedInteger) ,dimension(1) :: ind

    type(ColdIonIJIntDerivation) :: ciIJDeriv

    isJInt(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyIsJIntegral,.false.)
    ind(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyIndex,0)

    call jsonCont%load(isJInt)
    call jsonCont%output(isJInt)
    call jsonCont%load(ind)
    call jsonCont%output(ind)
    
    call ciIJDeriv%init(gridObj,ind(1)%value,isJInt=isJInt(1)%value)

    call textbookObj%addDerivation(ciIJDeriv,derivTag)

end subroutine addColdIonIJIntDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addIJIntDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
    !! Add a I/J integral derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj   
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedLogical) ,dimension(1) :: isJInt
    type(NamedInteger) ,dimension(1) :: ind

    type(IJIntDerivation) :: ijDeriv

    isJInt(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyIsJIntegral,.false.)
    ind(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyIndex,0)

    call jsonCont%load(isJInt)
    call jsonCont%output(isJInt)
    call jsonCont%load(ind)
    call jsonCont%output(ind)
    
    call ijDeriv%init(vSpaceObj,ind(1)%value,isJInt=isJInt(1)%value)

    call textbookObj%addDerivation(ijDeriv,derivTag)

end subroutine addIJIntDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addHarmonicExtractorDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
    !! Add a harmonic extractor derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj   
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedInteger) ,dimension(1) :: ind

    type(HExtractorDerivation) :: heDeriv

    ind(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyIndex,0)

    call jsonCont%load(ind)
    call jsonCont%output(ind)
    
    call heDeriv%init(vSpaceObj,ind(1)%value)

    call textbookObj%addDerivation(heDeriv,derivTag)

end subroutine addHarmonicExtractorDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addFScalingDerivationToTextbook(textbookObj,derivTag,vSpaceObj,partObj,jsonCont,mpiCont)
    !! Add a distribution scaling extrapolatio matrix derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj
    type(Partition)         ,intent(in)    :: partObj   
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedLogical) ,dimension(1) :: extToBound, staggeredVars ,leftBoundary

    type(FScalingDerivation) :: fScalDeriv

    extToBound(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyExtrapolateToBoundary,.false.)
    staggeredVars(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyStaggeredVars,.false.)
    leftBoundary(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyLeftBoundary,.false.)

    call jsonCont%load(extToBound)
    call jsonCont%output(extToBound)
    call jsonCont%load(staggeredVars)
    call jsonCont%output(staggeredVars)
    call jsonCont%load(leftBoundary)
    call jsonCont%output(leftBoundary)

    call fScalDeriv%init(partObj,mpiCont%getWorldRank(),vSpaceObj%getNumV(),leftBoundary=leftBoundary(1)%value,&
                         staggeredVars=staggeredVars(1)%value,extrapolateToBoundary=extToBound(1)%value)

    call textbookObj%addMatDerivation(fScalDeriv,derivTag)

end subroutine addFScalingDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addDDVDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
    !! Add a d/dv derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedRealArray) ,dimension(:) , allocatable:: outerV, innerV  ,vifAtZero 
    type(NamedInteger) ,dimension(1) :: targetH 

    type(DDVDerivation) :: ddvDeriv 

    character(len=30) :: intToStrBuffer


    integeR(ik) :: i 

    targetH(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyTargetH,0)
    

    call jsonCont%load(targetH)
    call jsonCont%output(targetH)

    if (targetH(1)%value > 0) then
        allocate(outerV(1))
        allocate(innerV(1))
        allocate(vifAtZero(1))
        
        outerV(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyOuterV,[(real(1,kind=rk),i=1,vSpaceObj%getNumV())])
        innerV(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyInnerV,[(real(1,kind=rk),i=1,vSpaceObj%getNumV())])
        vifAtZero(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyVIFAtZero,real([0,0],kind=rk))

        call jsonCont%load(outerV)
        call jsonCont%output(outerV)
        call jsonCont%load(innerV)
        call jsonCont%output(innerV)
        call jsonCont%load(vifAtZero)
        call jsonCont%output(vifAtZero)

        if (assertions .or. assertionLvl >= 0) then 
            call assert(size(outerV(1)%values) == vSpaceObj%getNumV(),outerV(1)%name//" must be of size numV")
            call assert(size(innerV(1)%values) == vSpaceObj%getNumV(),innerV(1)%name//" must be of size numV")
            call assert(size(vifAtZero(1)%values) == 2,vifAtZero(1)%name//" must be size 2")
        end if

        call ddvDeriv%init(vSpaceObj,outerV=removeName(outerV),innerV=removeName(innerV),&
                           vifAtZero=removeName(vifAtZero),targetH=targetH(1)%value)
    else
        allocate(outerV(vSpaceObj%getNumH()))
        allocate(innerV(vSpaceObj%getNumH()))
        allocate(vifAtZero(vSpaceObj%getNumH()))

        do i = 1,vSpaceObj%getNumH()
            intToStrBuffer=""
            write(intToStrBuffer,'(I0)') i
            outerV(i)%name = keyCustomDerivations//"."//derivTag//"."//keyOuterV//"."//keyHIs//trim(intToStrBuffer)
            outerV(i)%values = real([(1,i=1,vSpaceObj%getNumV())],kind=rk)
            innerV(i)%name = keyCustomDerivations//"."//derivTag//"."//keyInnerV//"."//keyHIs//trim(intToStrBuffer)
            innerV(i)%values = real([(1,i=1,vSpaceObj%getNumV())],kind=rk)
            vifAtZero(i) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyVIFAtZero//"."//keyHIs//trim(intToStrBuffer)&
                                         ,real([0,0],kind=rk))
        end do

        call jsonCont%load(outerV)
        call jsonCont%output(outerV)
        call jsonCont%load(innerV)
        call jsonCont%output(innerV)
        call jsonCont%load(vifAtZero)
        call jsonCont%output(vifAtZero)

        do i = 1,vSpaceObj%getNumH()
            if (assertions .or. assertionLvl >= 0) then 
                call assert(size(outerV(i)%values) == vSpaceObj%getNumV(),outerV(i)%name//" must be of size numV")
                call assert(size(innerV(i)%values) == vSpaceObj%getNumV(),innerV(i)%name//" must be of size numV")
                call assert(size(vifAtZero(i)%values) == 2,vifAtZero(i)%name//" must be size 2")
            end if
        end do
        call ddvDeriv%init(vSpaceObj,outerV=removeName(outerV),innerV=removeName(innerV),&
                           vifAtZero=removeName(vifAtZero))
    end if

    call textbookObj%addDerivation(ddvDeriv,derivTag)

end subroutine addDDVDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addD2DV2DerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
    !! Add a d^2/dv^2 derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedRealArray) ,dimension(:) , allocatable:: outerV, innerV  ,vidfdvAtZero 
    type(NamedInteger) ,dimension(1) :: targetH 

    type(D2DV2Derivation) :: d2dv2Deriv 

    character(len=30) :: intToStrBuffer


    integeR(ik) :: i 

    targetH(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyTargetH,0)
    
    call jsonCont%load(targetH)
    call jsonCont%output(targetH)

    if (targetH(1)%value > 0) then
        allocate(outerV(1))
        allocate(innerV(1))
        allocate(vidfdvAtZero(1))
        
        outerV(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyOuterV,[(real(1,kind=rk),i=1,vSpaceObj%getNumV())])
        innerV(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyInnerV,[(real(1,kind=rk),i=1,vSpaceObj%getNumV())])
        vidfdvAtZero(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyVIDFDVAtZero,real([0,0],kind=rk))

        call jsonCont%load(outerV)
        call jsonCont%output(outerV)
        call jsonCont%load(innerV)
        call jsonCont%output(innerV)
        call jsonCont%load(vidfdvAtZero)
        call jsonCont%output(vidfdvAtZero)

        if (assertions .or. assertionLvl >= 0) then 
            call assert(size(outerV(1)%values) == vSpaceObj%getNumV(),outerV(1)%name//" must be of size numV")
            call assert(size(innerV(1)%values) == vSpaceObj%getNumV(),innerV(1)%name//" must be of size numV")
            call assert(size(vidfdvAtZero(1)%values) == 2,vidfdvAtZero(1)%name//" must be size 2")
        end if

        call d2dv2Deriv%init(vSpaceObj,outerV=removeName(outerV),innerV=removeName(innerV),&
                           vidfdvAtZero=removeName(vidfdvAtZero),targetH=targetH(1)%value)
    else
        allocate(outerV(vSpaceObj%getNumH()))
        allocate(innerV(vSpaceObj%getNumH()))
        allocate(vidfdvAtZero(vSpaceObj%getNumH()))

        do i = 1,vSpaceObj%getNumH()
            intToStrBuffer=""
            write(intToStrBuffer,'(I0)') i
            outerV(i)%name = keyCustomDerivations//"."//derivTag//"."//keyOuterV//"."//keyHIs//trim(intToStrBuffer)
            outerV(i)%values = real([(1,i=1,vSpaceObj%getNumV())],kind=rk)
            innerV(i)%name = keyCustomDerivations//"."//derivTag//"."//keyInnerV//"."//keyHIs//trim(intToStrBuffer)
            innerV(i)%values = real([(1,i=1,vSpaceObj%getNumV())],kind=rk)
            vidfdvAtZero(i) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyVIFAtZero//"."&
                                            //keyHIs//trim(intToStrBuffer),real([0,0],kind=rk))
        end do

        call jsonCont%load(outerV)
        call jsonCont%output(outerV)
        call jsonCont%load(innerV)
        call jsonCont%output(innerV)
        call jsonCont%load(vidfdvAtZero)
        call jsonCont%output(vidfdvAtZero)

        do i = 1,vSpaceObj%getNumH()
            if (assertions .or. assertionLvl >= 0) then 
                call assert(size(outerV(i)%values) == vSpaceObj%getNumV(),outerV(i)%name//" must be of size numV")
                call assert(size(innerV(i)%values) == vSpaceObj%getNumV(),innerV(i)%name//" must be of size numV")
                call assert(size(vidfdvAtZero(i)%values) == 2,vidfdvAtZero(i)%name//" must be size 2")
            end if
        end do
        call d2dv2Deriv%init(vSpaceObj,outerV=removeName(outerV),innerV=removeName(innerV),&
                           vidfdvAtZero=removeName(vidfdvAtZero))
    end if

    call textbookObj%addDerivation(d2dv2Deriv,derivTag)

end subroutine addD2DV2DerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addMomentDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
    !! Add a custom moment derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedRealArray) ,dimension(1) :: gVec, varPowers
    type(NamedInteger) ,dimension(1) :: momentHarmonic, momentOrder 
    type(NamedReal) ,dimension(1) :: multConst

    type(MomentDerivation) :: momentDeriv 

    integer(ik) :: i

    momentHarmonic(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyMomentHarmonic,0)
    momentOrder(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyMomentOrder,0)
    
    call jsonCont%load(momentHarmonic)
    call jsonCont%output(momentHarmonic)
    call jsonCont%load(momentOrder)
    call jsonCont%output(momentOrder)

    varPowers(1)%name = keyCustomDerivations//"."//derivTag//"."//keyVarPowers
    allocate(varPowers(1)%values(0))

    gVec(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyGVector,&
                                 [(real(1,kind=rk),i=1,vSpaceObj%getNumV())])


    call jsonCont%load(gVec)
    call jsonCont%output(gVec)
    call jsonCont%load(varPowers)
    call jsonCont%output(varPowers)

    multConst(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyMultConst,real(1,kind=rk))

    call jsonCont%load(multConst)
    call jsonCont%output(multConst)

    call momentDeriv%init(momentOrder(1)%value,momentHarmonic(1)%value,vSpaceObj,&
                          varPowers(1)%values,gVec(1)%values,multConst(1)%value)

    call textbookObj%addDerivation(momentDeriv,derivTag)

end subroutine addMomentDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addLocalValExtractorDerivationToTextbook(textbookObj,derivTag,partObj,jsonCont,mpiCont)
    !! Add a local fluid value extraction derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(Partition)         ,intent(in)    :: partObj
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedInteger) ,dimension(1) :: targetX 

    type(LocValExtractorDerivation) :: locValDeriv

    targetX(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyTargetX,1)

    call jsonCont%load(targetX)
    call jsonCont%output(targetX)

    call locValDeriv%init(partObj,mpiCont%getWorldRank(),targetX(1)%value)

    call textbookObj%addDerivation(locValDeriv,derivTag)
    
end subroutine addLocalValExtractorDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addVelContracDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
    !! Add velocity space contraction derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedRealArray) ,dimension(1) :: gVec
    type(NamedInteger) ,dimension(1) :: targetHarmonic ,expectedNumHarmonics

    integer(ik) :: i

    type(VelContracDerivation) :: velContracDeriv

    targetHarmonic(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyTargetH,0)
    expectedNumHarmonics(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyExpectedNumHarmonics,0)

    gVec(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyGVector,&
    [(real(1,kind=rk),i=1,vSpaceObj%getNumV())])

    call jsonCont%load(gVec)
    call jsonCont%output(gVec)
    call jsonCont%load(targetHarmonic)
    call jsonCont%output(targetHarmonic)
    call jsonCont%load(expectedNumHarmonics)
    call jsonCont%output(expectedNumHarmonics)

    if (expectedNumHarmonics(1)%value > 0) then
        call velContracDeriv%init(targetHarmonic(1)%value,gVec(1)%values,vSpaceObj,expectedNumHarmonics(1)%value)
    else
        call velContracDeriv%init(targetHarmonic(1)%value,gVec(1)%values,vSpaceObj)
    end if

    call textbookObj%addDerivation(velContracDeriv,derivTag)

end subroutine addVelContracDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addVelTProdDeivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
    !! Add velocity vector tensor product derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VSpace)            ,intent(in)    :: vSpaceObj
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedRealArray) ,dimension(1) :: gVec
    type(NamedReal) ,dimension(1) :: shiftPower 

    integer(ik) :: i

    type(VelTProdDerivation) :: tProdDeriv

    shiftPower(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyPower,real(1,kind=rk))

    gVec(1) = NamedRealArray(keyCustomDerivations//"."//derivTag//"."//keyGVector,vSpaceObj%getVGrid())

    call jsonCont%load(gVec)
    call jsonCont%output(gVec)
    call jsonCont%load(shiftPower)
    call jsonCont%output(shiftPower)

    call tProdDeriv%init(vSpaceObj,gVec(1)%values,shiftPower(1)%value)

    call textbookObj%addDerivation(tProdDeriv,derivTag)

end subroutine addVelTProdDeivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addGenIntPolyDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
    !! Add a generalized integer powered polynomial derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedReal) ,dimension(1) :: multConst
    type(NamedRealArray) ,dimension(1) :: polyCoeffs
    type(NamedIntegerArray) ,allocatable ,dimension(:) :: polyPowers 
    type(NamedString) ,dimension(1) :: funcName

    integer(ik) ,allocatable ,dimension(:) :: maxPowers

    integer(ik) :: i ,j
    character(len=30) :: intToStrBuffer

    type(GenIntPolyFunDeriv) :: genIntPolyDeriv 

    multConst(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyMultConst,real(1,kind=rk))

    polyCoeffs(1)%name = keyCustomDerivations//"."//derivTag//"."//keyPolyCoeffs
    allocate(polyCoeffs(1)%values(0))

    funcName(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyFunctionName,"")

    call jsonCont%load(multConst)
    call jsonCont%load(polyCoeffs)
    call jsonCont%load(funcName)

    call jsonCont%output(multConst)
    call jsonCont%output(polyCoeffs)
    call jsonCont%output(funcName)

    allocate(polyPowers(size(polyCoeffs(1)%values)))

    do i = 1, size(polyCoeffs(1)%values)
        intToStrBuffer=""
        write(intToStrBuffer,'(I0)') i
        polyPowers(i)%name = keyCustomDerivations//"."//derivTag//"."//keyPolyPowers//"."//keyIndex//trim(intToStrBuffer)
        allocate(polyPowers(i)%values(0))
    end do

    call jsonCont%load(polyPowers)
    call jsonCont%output(polyPowers)


    do i = 1,size(polyCoeffs(1)%values)
        if (.not. allocated(maxPowers)) then
            allocate(maxPowers(size(polyPowers(i)%values)))
            maxPowers = polyPowers(i)%values
        end if
        if (assertions .or. assertionLvl >= 0) call assert(size(maxPowers) == size(polyPowers(i)%values),&
        "A polyPowers entry in "//derivTag//&
                                    " does not conform to expected number of variables")
        do j = 1,size(maxPowers)
            maxPowers(j) = max(maxPowers(j),polyPowers(i)%values(j))
        end do
    end do

    if (len(funcName(1)%value) > 0) then
        call genIntPolyDeriv%init(removeName(polyPowers),polyCoeffs(1)%values,maxPowers&
                                  ,funcName(1)%value,multConst=multConst(1)%value)
    else
        call genIntPolyDeriv%init(removeName(polyPowers),polyCoeffs(1)%values,maxPowers&
                                  ,multConst=multConst(1)%value)

    end if

    call textbookObj%addDerivation(genIntPolyDeriv,derivTag)

end subroutine addGenIntPolyDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addRangeFilterDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
    !! Add a range filter derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedString) ,dimension(1) :: derivName
    type(NamedIntegerArray) ,dimension(1) :: controlIndices, derivIndices
    type(NamedRealArray) ,allocatable ,dimension(:) :: controlRanges

    class(Derivation) ,allocatable :: derivBuffer

    type(RangeFilterDerivation) :: rangeFilterDeriv

    integer(ik) :: i
    character(len=30) :: intToStrBuffer

    controlIndices(1)%name = keyCustomDerivations//"."//derivTag//"."//keyControlIndices
    allocate(controlIndices(1)%values(0))

    derivIndices(1) = NamedIntegerArray(keyCustomDerivations//"."//derivTag//"."//keyDerivIndices,[0]) 

    derivName(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyRuleName,"")

    call jsonCont%load(derivName)
    call jsonCont%output(derivName)
    
    call jsonCont%load(controlIndices)
    call jsonCont%output(controlIndices)

    call jsonCont%load(derivIndices)
    call jsonCont%output(derivIndices)
    
    allocate(controlRanges(size(controlIndices(1)%values)))

    do i = 1, size(controlRanges)

        intToStrBuffer=""
        write(intToStrBuffer,'(I0)') i
        controlRanges(i)%name = keyCustomDerivations//"."//derivTag//"."//keyControlRanges//"."//keyIndex//trim(intToStrBuffer)
        allocate(controlRanges(i)%values(0))

    end do

    call jsonCont%load(controlRanges)
    call jsonCont%output(controlRanges)

    call textbookObj%copyDerivation(derivName(1)%value,derivBuffer)

    if (all(derivIndices(1)%values == [0])) then 
        call rangeFilterDeriv%init(derivBuffer,controlIndices(1)%values,removeName(controlRanges))
    else
        call rangeFilterDeriv%init(derivBuffer,controlIndices(1)%values,&
                                   removeName(controlRanges),derivIndices=derivIndices(1)%values)
    end if
    
    call textbookObj%addDerivation(rangeFilterDeriv,derivTag)

end subroutine addRangeFilterDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addCalculationTreeDerivationToTextbook(textbookObj,derivTag,varList,jsonCont,mpiCont)
    !! Add a calculation tree derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(VariableList)      ,intent(in)    :: varList
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedInteger) ,dimension(1) :: numNodes ,parentIndex
    type(NamedIntegerArray) ,dimension(1) :: childIndices,unaryIntParams
    type(NamedRealArray) ,dimension(1) :: unaryRealParams 
    type(NamedLogicalArray) ,dimension(1) :: unaryLogicalParams 
    type(NamedLogical) ,dimension(1) :: additiveMode 
    type(NamedString) ,dimension(1) :: unaryTransform, leafVariable 
    type(NamedReal) ,dimension(1) :: constantCoeff

    type(CalculationTreeDerivation) :: calcTreeDeriv
    type(FlatTree) :: fTree

    integer(ik) :: i ,leafVariableIndex
    character(len=30) :: intToStrBuffer

    real(rk) :: defaultConstant

    numNodes(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyNumNodes,0)

    call jsonCont%load(numNodes)
    call jsonCont%output(numNodes)

    allocate(fTree%kernels(numNodes(1)%value))
    allocate(fTree%children(numNodes(1)%value))
    allocate(fTree%parent(numNodes(1)%value))

    do i = 1,numNodes(1)%value

        intToStrBuffer=""
        write(intToStrBuffer,'(I0)') i

        parentIndex(1) = NamedInteger(keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                    //"."//keyParent,0)
        additiveMode(1) = NamedLogical(keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                    //"."//keyAdditiveMode,.false.)
        unaryTransform(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                    //"."//keyUnaryTransform,keyNone)
        leafVariable(1) = NamedString(keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                    //"."//keyLeafVar,keyNone)
        
        childIndices(1)%name = keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                //"."//keyChildren
        if (allocated(childIndices(1)%values)) deallocate(childIndices(1)%values)
        allocate(childIndices(1)%values(0))

        unaryIntParams(1)%name = keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                //"."//keyUnaryIntParams

        if (allocated(unaryIntParams(1)%values)) deallocate(unaryIntParams(1)%values)
        allocate(unaryIntParams(1)%values(0))

        unaryRealParams(1)%name = keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                //"."//keyUnaryRealParams
        if (allocated(unaryRealParams(1)%values)) deallocate(unaryRealParams(1)%values)
        allocate(unaryRealParams(1)%values(0))

        unaryLogicalParams(1)%name = keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                //"."//keyUnaryLogicalParams
        if (allocated(unaryLogicalParams(1)%values)) deallocate(unaryLogicalParams(1)%values)
        allocate(unaryLogicalParams(1)%values(0))


        call jsonCont%load(parentIndex)
        call jsonCont%output(parentIndex)
        call jsonCont%load(additiveMode)
        call jsonCont%output(additiveMode)
        call jsonCont%load(unaryIntParams)
        call jsonCont%output(unaryIntParams)
        call jsonCont%load(unaryRealParams)
        call jsonCont%output(unaryRealParams)
        call jsonCont%load(unaryLogicalParams)
        call jsonCont%output(unaryLogicalParams)
        call jsonCont%load(unaryTransform)
        call jsonCont%output(unaryTransform)
        call jsonCont%load(leafVariable)
        call jsonCont%output(leafVariable)
        call jsonCont%load(childIndices)
        call jsonCont%output(childIndices)


        defaultConstant = real(1,kind=rk)
        if (additiveMode(1)%value) defaultConstant = 0
        constantCoeff(1) = NamedReal(keyCustomDerivations//"."//derivTag//"."//keyNodes//"."//keyIndex//trim(intToStrBuffer)&
                                    //"."//keyConstant,defaultConstant)

        call jsonCont%load(constantCoeff)
        call jsonCont%output(constantCoeff)

        leafVariableIndex = 0 
        if (leafVariable(1)%value /= keyNone) then 
            call assert(varList%isVarNameRegistered(leafVariable(1)%value),&
            "Leaf variable "//leafVariable(1)%value// " in calculation tree not registered in variable list")
            leafVariableIndex =varList%getVarIndex(leafVariable(1)%value)
        end if

        fTree%kernels(i)%additiveMode = additiveMode(1)%value
        fTree%kernels(i)%constant = constantCoeff(1)%value
        fTree%kernels(i)%unaryTransformationTag = unaryTransform(1)%value
        fTree%kernels(i)%leafVarIndex = leafVariableIndex
        fTree%kernels(i)%unaryIntParams = unaryIntParams(1)%values
        fTree%kernels(i)%unaryLogicalParams = unaryLogicalParams(1)%values
        fTree%kernels(i)%unaryRealParams = unaryRealParams(1)%values

        fTree%parent(i) = parentIndex(1)%value
        fTree%children(i)%entry = childIndices(1)%values

    end do

    call calcTreeDeriv%init(fTree)

    call textbookObj%addDerivation(calcTreeDeriv,derivTag)

end subroutine addCalculationTreeDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addnDLinInterpDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
    !! Add an n-D linear interpolation derivation to textbook based on JSON config file

    type(Textbook)          ,intent(inout) :: textbookObj 
    character(*)            ,intent(in)    :: derivTag 
    type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
    type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    type(NamedStringArray)  ,dimension(1) :: gridNames
    type(NamedRealArray)    ,dimension(1) :: gridVals 
    type(NamedRealArray)    ,dimension(1) :: dataVals 
    type(NamedIntegerArray) ,dimension(1) :: dataDims

    type(Interpolation1D) ,allocatable ,dimension(:) :: interps1D
    type(InterpolationND) :: interpND
    type(FlatNDData)      :: flatData

    type(NDInterpolationDerivation) :: ndInterpDerivation

    integer(ik) :: i 

    gridNames(1)%name = keyCustomDerivations//"."//derivTag//"."//keyGrids//"."//keyNames
    allocate(gridNames(1)%values(0))

    call jsonCont%load(gridNames)
    call jsonCont%output(gridNames)
    allocate(interps1D(size(gridNames(1)%values)))
    allocate(gridVals(1)%values(0))
    do i = 1,size(gridNames(1)%values)
        gridVals(1)%name = keyCustomDerivations//"."//derivTag//"."//keyGrids//"."//gridNames(1)%values(i)%string
        call jsonCont%load(gridVals)
        call jsonCont%output(gridVals)

        call interps1D(i)%init(gridVals(1)%values)
    end do

    dataDims(1)%name = keyCustomDerivations//"."//derivTag//"."//keyData//"."//keyDims
    allocate(dataDims(1)%values(0))

    dataVals(1)%name = keyCustomDerivations//"."//derivTag//"."//keyData//"."//keyValues
    allocate(dataVals(1)%values(0))

    call jsonCont%load(dataDims)
    call jsonCont%output(dataDims)
    call jsonCont%load(dataVals)
    call jsonCont%output(dataVals)

    print *, dataDims(1)%values
    call flatData%directInit(dataVals(1)%values,dataDims(1)%values)
    call interpND%init(interps1D)

    call ndInterpDerivation%init(interpND,flatData)

    call textbookObj%addDerivation(ndInterpDerivation,derivTag)

end subroutine addnDLinInterpDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule custom_derivation_support_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
