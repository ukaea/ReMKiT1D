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
module initialization_support
    !! author: Stefan Mijin
    !!
    !! Contains various support routines used in data initialization

    use data_kinds                          ,only: rk, ik
    use runtime_constants                   ,only: debugging, assertions, assertionLvl
    use assertion_utility                   ,only: assert, assertIdentical, assertPure
    use variable_container_class            ,only: VariableContainer, CalculationRule
    use variable_list_class                 ,only: VariableList 
    use indexing_class                      ,only: Indexing
    use json_controller_class               ,only: JSONController 
    use mpi_controller_class                ,only: MPIController ,CommunicationData
    use derivation_abstract_class           ,only: Derivation
    use textbook_class                      ,only: Textbook
    use grid_class                          ,only: Grid
    use v_space_class                       ,only: VSpace
    use partition_class                     ,only: Partition
    use geometry_class                      ,only: Geometry
    use petsc_controller_class              ,only: PETScController, SolverOptions
    use hdf5_controller_class               ,only: HDF5Controller
    use species_list_class                  ,only: SpeciesList
    use species_class                       ,only: Species
    use simple_derivation_class             ,only: SimpleDerivation
    use additive_derivation_class           ,only: AdditiveDerivation
    use moment_derivation_class             ,only: MomentDerivation
    use sheath_gamma_derivation_class       ,only: ElectronSheathGammaDerivation
    use normalization_abstract_class        ,only: Normalization
    use implicit_PicardBDE_integrator_class ,only: PicardBDEIntegrator ,InternalControllerOptions
    use composite_integrator_class          ,only: CompositeIntegrator, IntegratorCallStep
    use explicit_rk_integrator_class        ,only: ExplicitRKIntegrator
    use simple_timestep_controller_class    ,only: SimpleTimestepController
    use signal_collection_class             ,only: SignalCollection
    use constant_signal_class               ,only: ConstSignal
    use hat_signal_class                    ,only: HatSignal
    use cut_sine_signal_class               ,only: CutSineSignal
    use interpolation_derivation_class      ,only: InterpolationDerivation
    use central_diff_grad_derivation_class  ,only: CentralDiffDerivation
    use coulomb_log_derivation_class        ,only: CoulombLogDerivation
    use maxwellian_derivation_class         ,only: MaxwellianDerivation
    use ccl_diff_derivation_class           ,only: CCLDiffDerivation
    use ccl_drag_derivation_class           ,only: CCLDragDerivation
    use ccl_weight_derivation_class         ,only: CCLWeightDerivation
    use cvode_integrator_class              ,only: CVODEOptions, CVODEIntegrator
    use status_printing
    use custom_derivation_support 
    use support_types
    use key_names
    use physical_constants

    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVarListFromJSON(varList,jsonCont,mpiCont,isDerivedList)
        !! Initialize variable list based on variable data from a JSON file. If isDerivedList = .true. the list is initialized
        !! based on derived list data, otherwise it uses the implicit list.

        type(VariableList)      ,intent(inout) :: varList
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
        logical ,optional       ,intent(in)    :: isDerivedList

    end subroutine initVarListFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initGridFromJSON(gridObj,jsonCont,mpiCont,lengthNorm)
        !! Initialize grid object based on data from a JSON file. 

        type(Grid)              ,intent(inout) :: gridObj
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
        real(rk) ,optional      ,intent(in)    :: lengthNorm  !! Length normalization if the supplied grid is in meters, defaults to 1

    end subroutine initGridFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initPartFromJSON(partObj,gridObj,jsonCont,mpiCont)
        !! Initialize simple partition object based on data from a JSON file and grid object

        type(Partition)         ,intent(inout) :: partObj
        type(Grid)              ,intent(in)    :: gridObj     !! Grid object used to get numX and numH
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine initPartFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initGeometryFromJSON(geometryObj,gridObj,jsonCont,mpiCont)
        !! Initialize geometry object based on data from a JSON file and grid object

        type(Geometry)          ,intent(inout) :: geometryObj
        type(Grid)              ,intent(in)    :: gridObj     !! Grid object used to infer cell widths for consistency
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine initGeometryFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initPETScContFromJSON(petscCont,indexingObj,jsonCont,mpiCont)
        !! Initialize PETSc controller object based on data from a JSON file and indexing object

        type(PETScController)   ,intent(inout) :: petscCont
        type(Indexing)          ,intent(in)    :: indexingObj 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController and to initialize PETSc

    end subroutine initPETScContFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initHDF5ContFromJSON(hdf5Cont,varCont,jsonCont,mpiCont)
        !! Initialize HDF5 controller object based on data from a JSON file and variable container object

        type(HDF5Controller)    ,intent(inout) :: hdf5Cont
        type(VariableContainer) ,intent(in)    :: varCont 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController and to initialize HDF5

    end subroutine initHDF5ContFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initVarContFromJSON(varCont,indexingObj,partObj,textbookObj,jsonCont,mpiCont)
        !! Initialize variable container based on variable data from a JSON file. 

        type(VariableContainer) ,intent(inout) :: varCont
        type(Indexing)          ,intent(in)    :: indexingObj !! Indexing object used in variable container initialization 
        type(Partition)         ,intent(in)    :: partObj     !! Partition object used in variable container initialization 
        type(Textbook)          ,intent(in)    :: textbookObj !! Textbook object used to retrieve derivation rules for variable container 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine initVarContFromJSON
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initStandardIntegrator(integratorObj,varCont,indexingObj,jsonCont,mpiCont)
        !! Initialize standard integrator based on data from a JSON file. The standard integrator is a composit integrator containing
        !! either RK or BDE integrators. 
    
        type(CompositeIntegrator) ,intent(inout) :: integratorObj 
        type(VariableContainer)   ,intent(in)    :: varCont
        type(Indexing)            ,intent(in)    :: indexingObj !! Indexing object used in BDE integrator initialization
        type(JSONController)      ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)       ,intent(inout) :: mpiCont     !! MPIController used with JSONController

    end subroutine initStandardIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initStandardTextbook(textbookObj,gridObj,geometryObj,partObj,&
                                           vSpaceObj,normObj,speciesListObj,varList,jsonCont,mpiCont)
        ! Initialize standard textbook object based on grid, geometry, partition, species, and normalization data as well as a JSON file. 
        !! Assumes that the velocity grid is normalized to electron thermal speed.
        !! The textbook includes (passed variables are in appearance order unless otherwise specified):
        !!
        !! "flowSpeedFromFlux" = flux/density (normalized to n_0 * u_0)
        !! "sonicSpeed" for each ion species  (negative IDs) = sqrt((poly_e * T_e + poly_i * T_i)/m_i) where i is the index of the ion species (normalized to u_0)
        !! "tempFromEnergy" for each species specified in config file as requiring temperature derivation = 2/3 * energy/density - mass * flux**2/(3*density**2)
        !! "leftElectronGamma" = sheath heat transmission coeffiecient for electrons at the left boundary using T_e and T_i
        !! "rightElectronGamma" = sheath heat transmission coefficient for electrons at the right boundary using T_e and T_i
        !! "densityMoment" = zeroth moment of f_0
        !! "energyMoment" = second order of f_0
        !! "cclDragCoeff" = Chang-Cooper-Langdon drag coefficient - assumes passed variable is the single harmonic f_0
        !! "cclDiffusionCoeff" = Chang-Cooper-Langdon diffusion coefficient - assumes passed variables are the single harmonic f_0 and cclWeight
        !! "cclWeight" = Chang-Cooper-Langdon interpolation weight - assumes passed variables are the single harmonic cclDragCoeff and cclDiffusionCoeff
        !! if l=1 is resolved:
        !! "fluxMoment" = first moment of f_1/3 (normalized to n_0*u_0)
        !! "heatFluxMoment" = third moment of f_1/3 (normalized to m_e * n_0 * v_th**3 / 2)
        !! if l=2 is resolved:
        !! "viscosityTensorxxMoment" = second moment of 2*f_2/15 (normalized to n_0 * e * T_0)
        !! Interpolation for staggered grids:   
        !! "gridToDual" = interpolates one variable from regular to staggered(dual) grid
        !! "dualToGrid" = interpolates one variable from staggered(dual) to regular grid
        !! "distributionInterp" = interpolates one distribution function (even harmonics to dual, odd to regular grid)
        !! Other useful derivations:
        !! "gradDeriv" = calculates gradient of variable on regular grid, interpolating on cell faces and extrapolating at boundaries
        !! "logLei" = calculates electron-ion Coulomb log taking in electron temperature and density (derivations added for each ion species)
        !! "logLee" = calculates electron self collision Coulomb log taking in electron temperature and density 
        !! "logLii" = calculates ion-ion collision frequency for each ion collision combination taking in the two species temperatures and 
        !!            densities (used as"logLiis_S" where s and S are the first and second ion species name) 
        !! "maxwellianDistribution" = calculates distribution function with Maxwellian l=0 harmonic and all other harmonics 0. Takes in temperature and density
        !! Also automatically includes all defined custom derivations
    
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

    end subroutine 
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initStandardSignals(signalCollectionObj)
    !! Initialize signal collection with commonly used signals and their names. The signals are:
    !! 1. constant signal
    !! 2. hat signal
    !! 3. cut sine signal

    type(SignalCollection) ,intent(inout) :: signalCollectionObj 

    end subroutine initStandardSignals
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine setDistHarmonic(dist,h,xArr,vArr)
        !! Set the h-th harmonic of distribution function variable dist to xArr*vArr at corresponding x, and v coordinates. 
        !! vArr should be either rank 1 or 2, and if it is rank 2 it should be size(numV,size(xArr))

        real(rk) ,dimension(:)  ,intent(inout) :: dist 
        integer(ik)             ,intent(in)    :: h
        real(rk) ,dimension(:)  ,intent(in)    :: xArr 
        real(rk) ,dimension(..) ,intent(in)    :: vArr
        
    end subroutine setDistHarmonic
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module initialization_support
!-----------------------------------------------------------------------------------------------------------------------------------
