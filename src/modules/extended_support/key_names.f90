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
module key_names
    !! author: Stefan Mijin 
    !!
    !! Provides centralized housing for JSON and other key names

    ! Grid key names
    character(len=*) ,parameter :: keyXGrid = "xGrid" 
    character(len=*) ,parameter :: keyVGrid = "vGrid" 
    character(len=*) ,parameter :: keyCellWidths = "cellWidths" 
    character(len=*) ,parameter :: keyCellCoords = "cellCentreCoords" 
    character(len=*) ,parameter :: keyPeriodic = "isPeriodic" 
    character(len=*) ,parameter :: keyBuildFromWidths = "constructFromWidths" 
    character(len=*) ,parameter :: keyMaxL = "maxL"
    character(len=*) ,parameter :: keyMaxM = "maxM"
    character(len=*) ,parameter :: keyFaceJacobians = "faceJacobians" 
    character(len=*) ,parameter :: keyLengthInMeters = "isLengthInMeters" 

    ! Variable keys
    character(len=*) ,parameter :: keyVariables = "variables" 
    character(len=*) ,parameter :: keyDerivedVars = "derivedVariables" 
    character(len=*) ,parameter :: keyImplicitVars = "implicitVariables" 
    character(len=*) ,parameter :: keyIsDist = "isDistribution"
    character(len=*) ,parameter :: keyIsScalar = "isScalar"
    character(len=*) ,parameter :: keyIsSingleHarmonic = "isSingleHarmonic"
    character(len=*) ,parameter :: keyIsOnDualGrid = "isOnDualGrid"
    character(len=*) ,parameter :: keyIsStationary = "isStationary"
    character(len=*) ,parameter :: keyInitVals = "initVals"
    character(len=*) ,parameter :: keyPriority = "priority"

    ! Library flags 

    character(len=*) ,parameter :: keyMPI = "MPI" 
    character(len=*) ,parameter :: keyCommData = "commData"
    character(len=*) ,parameter :: keyVarsToBroadcast = "varsToBroadcast"
    character(len=*) ,parameter :: keyHaloVars = "haloExchangeVars"
    character(len=*) ,parameter :: keyScalarsToBroadcast = "scalarVarsToBroadcast"
    character(len=*) ,parameter :: keyScalarsRoots = "scalarBroadcastRoots"
    character(len=*) ,parameter :: keyNumPX = "numProcsX" 
    character(len=*) ,parameter :: keyNumPH = "numProcsH" 
    character(len=*) ,parameter :: keyXHaloWidth = "xHaloWidth" 
    character(len=*) ,parameter :: keyPETSc = "PETSc" 
    character(len=*) ,parameter :: keySolOptions = "solverOptions" 
    character(len=*) ,parameter :: keyMaxSolIters = "maxSolverIters" 
    character(len=*) ,parameter :: keySolRelTol = "solverToleranceRel" 
    character(len=*) ,parameter :: keySolAbsTol = "solverToleranceAbs" 
    character(len=*) ,parameter :: keySolDivTol = "solverToleranceDiv" 
    character(len=*) ,parameter :: keySolKSP = "kspSolverType" 
    character(len=*) ,parameter :: keyHyprePC = "hyprePCType" 
    character(len=*) ,parameter :: keyPETScOpts = "PETScCommandLineOpts" 
    character(len=*) ,parameter :: keyHDF5 = "HDF5"
    character(len=*) ,parameter :: keyOutputVars = "outputVars"
    character(len=*) ,parameter :: keyInputVars = "inputVars" 


    ! Generic keys
    character(len=*) ,parameter :: keyNames = "names" 
    character(len=*) ,parameter :: keyNameSingle = "name" 
    character(len=*) ,parameter :: keyNone = "none"
    character(len=*) ,parameter :: keyFilepath = "filepath"
    character(len=*) ,parameter :: keyActive = "active"
    character(len=*) ,parameter :: keyID = "ID"
    character(len=*) ,parameter :: keyTags = "tags"
    character(len=*) ,parameter :: keyReqVarNames = "requiredVarNames"
    character(len=*) ,parameter :: keyReqVarPowers = "requiredVarPowers"
    character(len=*) ,parameter :: keyMultConst = "multConst"
    character(len=*) ,parameter :: keyType = "type"
    character(len=*) ,parameter :: keyGroupIndices = "groupIndices"
    character(len=*) ,parameter :: keyEvolvedVar = "evolvedVar"
    character(len=*) ,parameter :: keyImplicitVar = "implicitVar"
    character(len=*) ,parameter :: keyMode = "mode"
    character(len=*) ,parameter :: keyFrequency = "frequency"
    character(len=*) ,parameter :: keySave = "save"
    character(len=*) ,parameter :: keyLoad = "load"
    character(len=*) ,parameter :: keyAllTermImpGroups = "allTermsImplicitGroups"
    character(len=*) ,parameter :: keyAllTermGenGroups = "allTermsGeneralGroups"
    character(len=*) ,parameter :: keyImplicitTermGroups = "implicitGroups"
    character(len=*) ,parameter :: keyGeneralTermGroups = "generalGroups"
    character(len=*) ,parameter :: keySpeciesName = "speciesName"
    character(len=*) ,parameter :: keySpeciesNamesPlural = "speciesNames"
    character(len=*) ,parameter :: keyAdditionalColumnVars = "additionalColumnVariables"
    character(len=*) ,parameter :: keyColumnVarPowers = "columnVarPowers"
    character(len=*) ,parameter :: keyMultVar = "multiplierVariable"
    character(len=*) ,parameter :: keyIndex = "index"
    character(len=*) ,parameter :: keyPower = "power"
    character(len=*) ,parameter :: keyFixedMatrix = "fixedMatrix"
    character(len=*) ,parameter :: keyFunctionName = "functionName"
    character(len=*) ,parameter :: keyObjGroups = "objGroups"
    character(len=*) ,parameter :: keyData = "data"
    character(len=*) ,parameter :: keyDims = "dims"
    character(len=*) ,parameter :: keyValues = "values"
    character(len=*) ,parameter :: keyUpdate = "update"
    character(len=*) ,parameter :: keyAccumulate = "accumulate"
    
    ! Rule keys
    character(len=*) ,parameter :: keyDerivRule = "derivationRule" 
    character(len=*) ,parameter :: keyRuleName = "ruleName"

    ! Integrator keys
    character(len=*) ,parameter :: keyIntegrator = "integrator"
    character(len=*) ,parameter :: keyNumImplicitGroups = "numImplicitGroups"
    character(len=*) ,parameter :: keyNumGenGroups = "numGeneralGroups"
    character(len=*) ,parameter :: keyTimestepController = "timestepController"
    character(len=*) ,parameter :: keyRescaleTimestep = "rescaleTimestep"
    character(len=*) ,parameter :: keyUseMaxVal = "useMaxVal"
    character(len=*) ,parameter :: keyBDE = "backwardsEuler"
    character(len=*) ,parameter :: keyRK = "RK"
    character(len=*) ,parameter :: keyOrder = "order"
    character(len=*) ,parameter :: keyNonlinTol = "nonlinTol"
    character(len=*) ,parameter :: keyAbsTol = "absTol"
    character(len=*) ,parameter :: keyConvergenceVars = "convergenceVars"
    character(len=*) ,parameter :: keyMaxNonlinIters = "maxNonlinIters"
    character(len=*) ,parameter :: keyAssociatedPETScGroup = "associatedPETScGroup"
    character(len=*) ,parameter :: keyUse2Norm = "use2Norm"
    character(len=*) ,parameter :: keyInternalStepControl = "internalStepControl"
    character(len=*) ,parameter :: keyStartingNumSteps = "startingNumSteps"
    character(len=*) ,parameter :: keyStepMultiplier = "stepMultiplier"
    character(len=*) ,parameter :: keyMinNumNonlinInters = "minNumNonlinIters"
    character(len=*) ,parameter :: keyMaxBDERestarts = "maxBDERestarts"
    character(len=*) ,parameter :: keyBDEConsolidationInterval = "BDEConsolidationInterval"
    character(len=*) ,parameter :: keyStepDecrament = "stepDecrament"
    character(len=*) ,parameter :: keyStepTags = "stepTags"
    character(len=*) ,parameter :: keyIntegratorTags = "integratorTags"
    character(len=*) ,parameter :: keyIntegratorTagSingle = "integratorTag"
    character(len=*) ,parameter :: keyEvolvedModels = "evolvedModels"
    character(len=*) ,parameter :: keyInternalUpdateGroups = "internallyUpdatedGroups"
    character(len=*) ,parameter :: keyInternalUpdateModelData = "internallyUpdateModelData"
    character(len=*) ,parameter :: keyGlobalStepFraction = "globalStepFraction"
    character(len=*) ,parameter :: keyUseInitialInput = "useInitialInput"
    character(len=*) ,parameter :: keyAllowTimeEvolution = "allowTimeEvolution"
    character(len=*) ,parameter :: keyInitialTimestep = "initialTimestep"
    character(len=*) ,parameter :: keyRelaxationWeight = "relaxationWeight"
    character(len=*) ,parameter :: keyCVODE = "CVODE"
    character(len=*) ,parameter :: keyCVODEGMRESMaxRestarts = "maxRestarts"
    character(len=*) ,parameter :: keyRelTol = "relTol"
    character(len=*) ,parameter :: keyBBDPreParams = "CVODEPreBBDParams"
    character(len=*) ,parameter :: keyCVODEAM = "CVODEUseAdamsMoulton"
    character(len=*) ,parameter :: keyCVODEStabDet = "CVODEUseStabLimDet"
    character(len=*) ,parameter :: keyCVODEMaxOrder = "CVODEMaxOrder"
    character(len=*) ,parameter :: keyCVODEMaxInternalSteps = "CVODEMaxInternalSteps"
    character(len=*) ,parameter :: keyCVODEMaxStep = "CVODEMaxStepSize"
    character(len=*) ,parameter :: keyCVODEMinStep = "CVODEMinStepSize"
    character(len=*) ,parameter :: keyCVODEInitStep = "CVODEInitStepSize"
    

    ! Manipulator keys 

    character(len=*) ,parameter :: keyManipulators = "manipulators"
    character(len=*) ,parameter :: keyGroupEvaluator = "groupEvaluator"
    character(len=*) ,parameter :: keyEvaluatedGroup = "evaluatedTermGroup"
    character(len=*) ,parameter :: keyResultVarName = "resultVarName"
    character(len=*) ,parameter :: keyTermEvaluator = "termEvaluator"
    character(len=*) ,parameter :: keyEvaluatedModelNames = "evaluatedModelNames"
    character(len=*) ,parameter :: keyEvaluatedTermNames = "evaluatedTermNames"
    character(len=*) ,parameter :: keyExtractor = "modelboundDataExtractor"
    character(len=*) ,parameter :: keyModelTag = "modelTag"
    character(len=*) ,parameter :: keyModelboundDataName = "modelboundDataName"

    ! General model keys 

    character(len=*) ,parameter :: keyModels = "models"

    ! Custom model keys 

    character(len=*) ,parameter :: keyCustomModel = "customModel"
    character(len=*) ,parameter :: keyTermTags = "termTags"
    character(len=*) ,parameter :: keyVarData = "varData"
    character(len=*) ,parameter :: keyReqRowVarNames = "requiredRowVarNames"
    character(len=*) ,parameter :: keyReqRowVarPowers = "requiredRowVarPowers"
    character(len=*) ,parameter :: keyReqColVarNames = "requiredColVarNames"
    character(len=*) ,parameter :: keyReqColVarPowers = "requiredColVarPowers"
    character(len=*) ,parameter :: keyReqMBRowVarNames = "requiredMBRowVarNames"
    character(len=*) ,parameter :: keyReqMBRowVarPowers = "requiredMBRowVarPowers"
    character(len=*) ,parameter :: keyReqMBColVarNames = "requiredMBColVarNames"
    character(len=*) ,parameter :: keyReqMBColVarPowers = "requiredMBColVarPowers"
    character(len=*) ,parameter :: keyTermType = "termType"
    character(len=*) ,parameter :: keyCustomNormConst = "customNormConst"
    character(len=*) ,parameter :: keyNormNames = "normNames"
    character(len=*) ,parameter :: keyNormPowers = "normPowers"
    character(len=*) ,parameter :: keyTimeSignalData = "timeSignalData"
    character(len=*) ,parameter :: keyStencilData = "stencilData"
    character(len=*) ,parameter :: keyStencilType = "stencilType"
    character(len=*) ,parameter :: keyModelboundData = "modelboundData"
    character(len=*) ,parameter :: keyModelboundDataType = "modelboundDataType"
    character(len=*) ,parameter :: keyTermGenerators = "termGenerators"
    character(len=*) ,parameter :: keySpatialProfile = "spatialProfile"
    character(len=*) ,parameter :: keyHarmonicProfile = "harmonicProfile"
    character(len=*) ,parameter :: keyVelocityProfile = "velocityProfile"
    character(len=*) ,parameter :: keyTimeSignalPeriod = "timeSignalPeriod"
    character(len=*) ,parameter :: keyTimeSignalParams = "timeSignalParams"
    character(len=*) ,parameter :: keyRealTimePeriod = "realTimePeriod"
    character(len=*) ,parameter :: keySkipPattern = "skipPattern"
    character(len=*) ,parameter :: keyMultCopyTermName = "multCopyTermName"
    character(len=*) ,parameter :: keyReqMBVarName = "requiredMBVarName"
    character(len=*) ,parameter :: keyMatrixTerm = "matrixTerm"
    character(len=*) ,parameter :: keyDerivationTerm = "derivationTerm"

    ! Term generator keys 

    character(len=*) ,parameter :: keyGeneratorTag = "generatorTag"
    
    character(len=*) ,parameter :: keyCRMDensTermGen = "CRMDensityEvolution"
    character(len=*) ,parameter :: keyEvolvedSpeciesIDs = "evolvedSpeciesIDs"
    character(len=*) ,parameter :: keyIncludedTransitionInds = "includedTransitionIndices"

    character(len=*) ,parameter :: keyCRMElEnergyTermGen = "CRMElectronEnergyEvolution"
    character(len=*) ,parameter :: keyElectronEnergyVar = "electronEnergyDensity"

    character(len=*) ,parameter :: keyCRMBoltzTermGen = "CRMFixedBoltzmannCollInt"
    character(len=*) ,parameter :: keyEvolvedHarmonicSingle = "evolvedHarmonic"
    character(len=*) ,parameter :: keyAssociatedVarIndex = "associatedVarIndex"
    character(len=*) ,parameter :: keyFixedEnergyIndices = "fixedEnergyIndices"

    character(len=*) ,parameter :: keyCRMSecElTermGen = "CRMSecondaryElectronTerms"

    character(len=*) ,parameter :: keyCRMVarBoltzTermGen = "CRMVariableBoltzmannCollInt"

    ! Custom modelbound data keys

    !Varlike data
    character(len=*) ,parameter :: keyVarLikeMB = "varlikeData"
    character(len=*) ,parameter :: keyDataNames = "dataNames"
    character(len=*) ,parameter :: keyIsDerivedFromOtherData = "isDerivedFromOtherData"
    character(len=*) ,parameter :: keyDerivationPriority = "derivationPriority"

    !CRM data
    character(len=*) ,parameter :: keyCRMData = "modelboundCRMData"
    character(len=*) ,parameter :: keyElState = "electronStateID"
    character(len=*) ,parameter :: keyInelGrid = "inelasticGridData"
    character(len=*) ,parameter :: keyFixedTransitionEnergies = "fixedTransitionEnergies"
    character(len=*) ,parameter :: keyTransitionTags = "transitionTags"
    character(len=*) ,parameter :: keyTransitions = "transitions"
    character(len=*) ,parameter :: keyIngoingState = "ingoingState"
    character(len=*) ,parameter :: keyIngoingStatePlural = "ingoingStates"
    character(len=*) ,parameter :: keyOutgoingState = "outgoingState"
    character(len=*) ,parameter :: keyOutgoingStatePlural = "outgoingStates"
    character(len=*) ,parameter :: keyRate = "rate"
    character(len=*) ,parameter :: keyFixedEnergyIndex = "fixedEnergyIndex"
    character(len=*) ,parameter :: keyDirectTransitionEnergyIndex = "directTransitionFixedEnergyIndex"
    character(len=*) ,parameter :: keyFixedEnergy = "fixedEnergy"
    character(len=*) ,parameter :: keyTakeMomentumMoment = "takeMomentumMoment"
    character(len=*) ,parameter :: keyDistributionVarName = "distributionVarName"
    character(len=*) ,parameter :: keyTransitionIndex = "transitionIndex"
    character(len=*) ,parameter :: keyDirectTransitionIndex = "directTransitionIndex"
    character(len=*) ,parameter :: keyElectronTemperatureVar = "electronTemperatureVar"
    character(len=*) ,parameter :: keyCrossSectionData = "crossSectionData"
    character(len=*) ,parameter :: keyPresentHarmonics = "presentHarmonics"
    character(len=*) ,parameter :: keyMaxCrossSectionL = "maxCrossSectionL"
    character(len=*) ,parameter :: keyLPrefix = "l="
    character(len=*) ,parameter :: keyDegeneracyRule = "degeneracyRuleName"
    character(len=*) ,parameter :: keyDegeneracyRuleReqVars = "degeneracyRuleReqVars"
    character(len=*) ,parameter :: keyFixedDegenRatio = "fixedDegeneracyRatio"
    character(len=*) ,parameter :: keyCSUpdatePriority = "crossSectionUpdatePriority"
    character(len=*) ,parameter :: keySimpleTransition = "simpleTransition"
    character(len=*) ,parameter :: keyFixedECSTransition = "fixedECSTransition"
    character(len=*) ,parameter :: keyDBTransition = "detailedBalanceTransition"
    character(len=*) ,parameter :: keyDerivedTransition = "derivedTransition"
    character(len=*) ,parameter :: keyMomentumRule = "momentumRateDerivationRule"
    character(len=*) ,parameter :: keyMomentumReqVarNames = "momentumRateDerivationReqVarNames"
    character(len=*) ,parameter :: keyEnergyRule = "energyRateDerivationRule"
    character(len=*) ,parameter :: keyEnergyReqVarNames = "energyRateDerivationReqVarNames"
    character(len=*) ,parameter :: keyJanevRadRecomb = "JanevRadRecomb"
    character(len=*) ,parameter :: keyJanevCollExIon = "JanevCollExIon"
    character(len=*) ,parameter :: keyJanevCollDeexRecomb = "JanevCollDeexRecomb"
    character(len=*) ,parameter :: keyStartHState = "startHState"
    character(len=*) ,parameter :: keyEndHState = "endHState"
    character(len=*) ,parameter :: keyVariableECSTransition = "variableECSTransition"
    character(len=*) ,parameter :: keyCrossSectionDerivations = "crossSectionDerivations"
    character(len=*) ,parameter :: keyCSDerivationHarmonics = "crossSectionDerivationHarmonics"
    character(len=*) ,parameter :: keyEnergyDerivationName = "energyDerivationName"
    character(len=*) ,parameter :: keyEnergyDerivationReqVars = "energyDerivationReqVars"

    ! LBC data

    character(len=*) ,parameter :: keyLBCData = "modelboundLBCData"
    character(len=*) ,parameter :: keyIonCurrentVar = "ionCurrentVarName"
    character(len=*) ,parameter :: keyTotalCurrentVar = "totalCurrentVarName"
    character(len=*) ,parameter :: keyBisTol = "bisectionTolerance"

    ! Stencil template keys 

    character(len=*) ,parameter :: keyInterpolatedVar = "interpolatedVarName"
    character(len=*) ,parameter :: keyDiagonal = "diagonalStencil"
    character(len=*) ,parameter :: keyEvolvedXCells = "evolvedXCells"
    character(len=*) ,parameter :: keyStaggeredDiff = "staggeredDifferenceStencil"
    character(len=*) ,parameter :: keyCentralDiffInterp = "centralDifferenceInterpolated"
    character(len=*) ,parameter :: keyUpwindedDiff = "upwindedDifference"
    character(len=*) ,parameter :: keyIgnoreJacobian = "ignoreJacobian"
    character(len=*) ,parameter :: keyLeftBoundary = "leftBoundary"
    character(len=*) ,parameter :: keyDontExtrapolate = "dontExtrapolate"
    character(len=*) ,parameter :: keyLowerBoundVar = "lowerBoundVar"
    character(len=*) ,parameter :: keyFluxJacVar = "fluxJacVar"
    character(len=*) ,parameter :: keyFixedLowerBound = "fixedLowerBound"
    character(len=*) ,parameter :: keyNoLowerBound = "noLowerBound"
    character(len=*) ,parameter :: keyBoundaryStencil = "boundaryStencil"
    character(len=*) ,parameter :: keyDiffusionStencil = "diffusionStencil"
    character(len=*) ,parameter :: keyDoNotInterpolateD = "doNotInterpolateDiffCoeff"
    character(len=*) ,parameter :: keyEvolvedHarmonics = "evolvedHarmonics"
    character(len=*) ,parameter :: keyEvolvedVCells = "evolvedVCells"
    character(len=*) ,parameter :: keyMomentStencil = "momentStencil"
    character(len=*) ,parameter :: keyMomentHarmonic = "momentHarmonic"
    character(len=*) ,parameter :: keyMomentOrder = "momentOrder"
    character(len=*) ,parameter :: keyKinSpatialDiffStencil = "kineticSpatialDiffStencil"
    character(len=*) ,parameter :: keyRowHarmonic = "rowHarmonic"
    character(len=*) ,parameter :: keyColHarmonic = "colHarmonic"
    character(len=*) ,parameter :: keyDDVStencil = "ddvStencil"
    character(len=*) ,parameter :: keyFixedC = "fixedC"
    character(len=*) ,parameter :: keyFixedInterp = "fixedInterp"
    character(len=*) ,parameter :: keyModelboundC = "modelboundC"
    character(len=*) ,parameter :: keyModelboundInterp = "modelboundInterp"
    character(len=*) ,parameter :: keyCFAtZero = "cfAtZero"
    character(len=*) ,parameter :: keyVDiffStencil = "vDiffusionStencil"
    character(len=*) ,parameter :: keyFixedA = "fixedA"
    character(len=*) ,parameter :: keyModelboundA = "modelboundA"
    character(len=*) ,parameter :: keyADFAtZero = "adfAtZero"
    character(len=*) ,parameter :: keyShkarofskyIJ = "shkarofskyIJStencil"
    character(len=*) ,parameter :: keyJInt = "JIntegral"
    character(len=*) ,parameter :: keyIntegralIndex = "integralIndex"
    character(len=*) ,parameter :: keyFixedBoltzmann = "boltzmannStencil"
    character(len=*) ,parameter :: keyAbsorptionTerm = "absorptionTerm"
    character(len=*) ,parameter :: keyDBTerm = "detailedBalanceTerm"
    character(len=*) ,parameter :: keyScalingLBC = "scalingLogicalBoundaryStencil"
    character(len=*) ,parameter :: keyDecompHarmonics = "includedDecompHarmonics"
    character(len=*) ,parameter :: keyTermMomentStencil = "termMomentStencil"
    character(len=*) ,parameter :: keyTermName = "termName"
    character(len=*) ,parameter :: keyVariableBoltzmann = "variableBoltzmannStencil"
    character(len=*) ,parameter :: keySuperelasticTerm = "superelasticTerm"
    character(len=*) ,parameter :: keyCustomFluid1DStencil = "customFluid1DStencil"
    character(len=*) ,parameter :: keyXStencil = "xStencil"
    character(len=*) ,parameter :: keyColumnVec = "columnVector"
    character(len=*) ,parameter :: keyColumnVarContVars = "columnVarContVars"
    character(len=*) ,parameter :: keyColumnMBDataVars = "columnMBDataVars"

    ! Signal keys 

    character(len=*) ,parameter :: keySignal = "timeSignalType"
    character(len=*) ,parameter :: keyConstSignal = "constant"
    character(len=*) ,parameter :: keyHatSignal = "hat"
    character(len=*) ,parameter :: keyCutSineSignal = "cutSine"

    ! Normalization keys

    character(len=*) ,parameter :: keyNormalization = "normalization"
    character(len=*) ,parameter :: keyTimeNorm = "time"
    character(len=*) ,parameter :: keyTempEVNorm = "eVTemperature"
    character(len=*) ,parameter :: keyDensNorm = "density"
    character(len=*) ,parameter :: keyVelGridNorm = "velGrid"
    character(len=*) ,parameter :: keyLengthNorm = "length"
    character(len=*) ,parameter :: keySpeedNorm = "speed"
    character(len=*) ,parameter :: keyEFieldNorm = "EField"
    character(len=*) ,parameter :: keyHeatFluxNorm = "heatFlux"
    character(len=*) ,parameter :: keyCrossSectionNorm = "crossSection"
    character(len=*) ,parameter :: keyIonZNorm = "referenceIonZ"

    ! Species keys 

    character(len=*) ,parameter :: keySpecies = "species"
    character(len=*) ,parameter :: keyMassA = "atomicMass"
    character(len=*) ,parameter :: keyCharge = "charge"
    character(len=*) ,parameter :: keyAssociatedVars = "associatedVars"

    ! Standard textbook 

    character(len=*) ,parameter :: keyStandardTextbook = "standardTextbook"
    character(len=*) ,parameter :: keyTempDerivIDs = "temperatureDerivSpeciesIDs"
    character(len=*) ,parameter :: keyElectronPolytropicC = "electronPolytropicCoeff"
    character(len=*) ,parameter :: keyIonPolytropicC = "ionPolytropicCoeff"
    character(len=*) ,parameter :: keySheathGammaIonSpeciesID = "electronSheathGammaIonSpeciesID"
    character(len=*) ,parameter :: keyRemoveLogLeiDisc = "removeLogLeiDiscontinuity"

    ! Custom derivation keys 

    character(len=*) ,parameter :: keyCustomDerivations = "customDerivations"
    character(len=*) ,parameter :: keyVarPowers = "varPowers"
    character(len=*) ,parameter :: keySimpleDeriv = "simpleDerivation"
    character(len=*) ,parameter :: keyAdditiveDeriv = "additiveDerivation"
    character(len=*) ,parameter :: keyMultDeriv = "multiplicativeDerivation"
    character(len=*) ,parameter :: keyDerivTags = "derivationTags"
    character(len=*) ,parameter :: keyResultPower = "resultPower"
    character(len=*) ,parameter :: keyLinearCoefficients = "linearCoefficients"
    character(len=*) ,parameter :: keyDerivIndices = "derivationIndices"
    character(len=*) ,parameter :: keyInnerDerivation = "innerDerivation"
    character(len=*) ,parameter :: keyOuterDerivation = "outerDerivation"
    character(len=*) ,parameter :: keyInnerDerivPower = "innerDerivPower"
    character(len=*) ,parameter :: keyOuterDerivPower = "outerDerivPower"
    character(len=*) ,parameter :: keyInnerDerivIndices = "innerDerivIndices"
    character(len=*) ,parameter :: keyOuterDerivIndices = "outerDerivIndices"
    character(len=*) ,parameter :: keyInnerDerivFuncName = "innerDerivFuncName"
    character(len=*) ,parameter :: keyPolyFunDeriv = "polynomialFunctionDerivation"
    character(len=*) ,parameter :: keyPolyPowers = "polynomialPowers"
    character(len=*) ,parameter :: keyPolyCoeffs = "polynomialCoefficients"
    character(len=*) ,parameter :: keyConstCoeff = "constantPolynomialCoefficient"
    character(len=*) ,parameter :: keyBoundedExtDeriv = "boundedExtrapolationDerivation"
    character(len=*) ,parameter :: keyFixedUpperBound = "fixedUpperBound"
    character(len=*) ,parameter :: keyIgnoreUpperBound = "ignoreUpperBound"
    character(len=*) ,parameter :: keyIgnoreLowerBound = "ignoreLowerBound"
    character(len=*) ,parameter :: keyExpectLowerBoundVar = "expectLowerBoundVar"
    character(len=*) ,parameter :: keyExpectUpperBoundVar = "expectUpperBoundVar"
    character(len=*) ,parameter :: keyStaggeredVars = "staggeredVars"
    character(len=*) ,parameter :: keyCIIJDeriv = "coldIonIJIntegralDerivation"
    character(len=*) ,parameter :: keyIsJIntegral = "isJIntegral"
    character(len=*) ,parameter :: keyIJDeriv = "IJIntegralDerivation"
    character(len=*) ,parameter :: keyHarmonicExtractor= "harmonicExtractorDerivation"
    character(len=*) ,parameter :: keyFScalingExtrapolation = "distScalingExtrapDerivation"
    character(len=*) ,parameter :: keyExtrapolateToBoundary = "extrapolateToBoundary"
    character(len=*) ,parameter :: keyDDVDeriv = "ddvDerivation"
    character(len=*) ,parameter :: keyOuterV = "outerV"
    character(len=*) ,parameter :: keyInnerV = "innerV"
    character(len=*) ,parameter :: keyVIFAtZero = "vifAtZero"
    character(len=*) ,parameter :: keyTargetH = "targetH"
    character(len=*) ,parameter :: keyHIs = "h="
    character(len=*) ,parameter :: keyD2DV2Deriv = "d2dv2Derivation"
    character(len=*) ,parameter :: keyVIDFDVAtZero = "vidfdvAtZero"
    character(len=*) ,parameter :: keyMomentDerivation = "momentDerivation"
    character(len=*) ,parameter :: keyGVector = "gVector"
    character(len=*) ,parameter :: keyLocValExtractor= "locValExtractorDerivation"
    character(len=*) ,parameter :: keyTargetX= "targetX"
    character(len=*) ,parameter :: keyVelContracDeriv= "velocityContractionDerivation"
    character(len=*) ,parameter :: keyExpectedNumHarmonics= "expectedNumberOfHarmonics"
    character(len=*) ,parameter :: keyVelTProdDeriv= "velocityTensorProdDerivation"
    character(len=*) ,parameter :: keyGenIntPolyFunDeriv = "generalizedIntPowerPolyDerivation"
    character(len=*) ,parameter :: keyRangeFilterDerivation = "rangeFilterDerivation"
    character(len=*) ,parameter :: keyControlIndices = "controlIndices"
    character(len=*) ,parameter :: keyControlRanges = "controlRanges"
    character(len=*) ,parameter :: keyCalculationTreeDerivation = "calculationTreeDerivation"
    character(len=*) ,parameter :: keyNodes = "nodes"
    character(len=*) ,parameter :: keyNumNodes = "numNodes"
    character(len=*) ,parameter :: keyChildren = "children"
    character(len=*) ,parameter :: keyParent = "parent"
    character(len=*) ,parameter :: keyAdditiveMode = "isAdditiveNode"
    character(len=*) ,parameter :: keyConstant = "constant"
    character(len=*) ,parameter :: keyLeafVar = "leafVariable"
    character(len=*) ,parameter :: keyUnaryTransform = "unaryTransform"
    character(len=*) ,parameter :: keyUnaryRealParams = "unaryRealParams"
    character(len=*) ,parameter :: keyUnaryIntParams = "unaryIntegerParams"
    character(len=*) ,parameter :: keyUnaryLogicalParams = "unaryLogicalParams"
    character(len=*) ,parameter :: keynDLinInterpDerivation = "nDLinInterpDerivation"
    character(len=*) ,parameter :: keyGrids = "grids"


    !Extrapolation keys 

    character(len=*) ,parameter :: keyExtrapolation = "extrapolationStrategy"
    character(len=*) ,parameter :: keyLinExtrapolation = "linExtrapolation"
    character(len=*) ,parameter :: keyLogExtrapolation = "logExtrapolation"
    character(len=*) ,parameter :: keyLinlogExtrapolation = "linlogExtrapolation"
    character(len=*) ,parameter :: keyExpectedHaloWidth = "expectedHaloWidth"

    ! Timeloop keys

    character(len=*) ,parameter :: keyTimeloop = "timeloop"
    character(len=*) ,parameter :: keyFixedSteps = "fixedNumSteps"
    character(len=*) ,parameter :: keyNumTimestep = "numTimesteps"
    character(len=*) ,parameter :: keySaveInterval = "fixedSaveInterval"
    character(len=*) ,parameter :: keyNormTimeTarget = "normalizedTimeTarget"
    character(len=*) ,parameter :: keyRealTimeTarget = "realTimeTarget"
    character(len=*) ,parameter :: keyTimeTarget = "timeValueTarget"
    character(len=*) ,parameter :: keyMinSaveInterval = "minimumSaveInterval"
    character(len=*) ,parameter :: keyOutputMode = "outputMode"
    character(len=*) ,parameter :: keyRestart = "restart"
    character(len=*) ,parameter :: keyResetTime = "resetTime"
    character(len=*) ,parameter :: keyHDF5LoadInit = "loadInitValsFromHDF5"
    character(len=*) ,parameter :: keyInitPath = "initValFilename"
    character(len=*) ,parameter :: keyOutputDrivenMode = "outputDriven"
    character(len=*) ,parameter :: keyOutputPoints = "outputPoints"
    character(len=*) ,parameter :: keyInitialOutputIndex = "initialOutputIndex"

!----------------------------------------------------------------------------------------------------------------------------------
end module key_names
!----------------------------------------------------------------------------------------------------------------------------------
 
