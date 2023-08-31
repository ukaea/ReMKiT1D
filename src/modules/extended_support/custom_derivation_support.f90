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
module custom_derivation_support
    !! author: Stefan Mijin
    !!
    !! Contains support for adding custom derivations to a textbook

    use data_kinds                             ,only: rk, ik
    use runtime_constants                      ,only: debugging, assertions
    use assertion_utility                      ,only: assert, assertIdentical, assertPure
    use textbook_class                         ,only: Textbook
    use grid_class                             ,only: Grid
    use partition_class                        ,only: Partition
    use v_space_class                          ,only: VSpace
    use normalization_abstract_class           ,only: Normalization
    use variable_list_class                    ,only: VariableList
    use species_list_class                     ,only: SpeciesList
    use json_controller_class                  ,only: JSONController
    use mpi_controller_class                   ,only: MPIController
    use geometry_class                         ,only: Geometry
    use simple_derivation_class                ,only: SimpleDerivation
    use polynomial_fun_derivation_class        ,only: PolyFunDeriv
    use additive_derivation_class              ,only: AdditiveDerivation
    use multiplicative_derivation_class        ,only: MultiplicativeDerivation
    use bounded_ext_derivation_class           ,only: BoundedExtDerivation
    use cold_ion_ij_int_derivation_class       ,only: ColdIonIJIntDerivation
    use ij_int_derivation_class                ,only: IJIntDerivation
    use harmonic_extractor_derivation_class    ,only: HExtractorDerivation
    use f_scaling_derivation_class             ,only: FScalingDerivation
    use ddv_derivation_class                   ,only: DDVDerivation
    use d2dv2_derivation_class                 ,only: D2DV2Derivation
    use derivation_abstract_class              ,only: Derivation
    use extrapolation_abstract_class           ,only: Extrapolation
    use moment_derivation_class                ,only: MomentDerivation
    use loc_val_extractor_derivation_class     ,only: LocValExtractorDerivation
    use vel_contraction_derivation_class       ,only: VelContracDerivation
    use vel_tensor_prod_derivation_class       ,only: VelTProdDerivation
    use gen_int_polynomial_fun_derivation_class,only: GenIntPolyFunDeriv
    use range_filter_derivation_class          ,only: RangeFilterDerivation
    use calculation_tree_class                 ,only: FlatTree
    use calculation_tree_derivation_class      ,only: CalculationTreeDerivation
    use lin_interpnd_class                     ,only: InterpolationND
    use lin_interp1D_class                     ,only: Interpolation1D
    use lin_interpnd_derivation_class          ,only: NDInterpolationDerivation
    use flat_nd_data_class                     ,only: FlatNDData
    use extrapolation_initialization_support
    use support_types
    use support_functions
    use key_names

    implicit none 
    private 

    public :: addCustomDerivationsToTextbook

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
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

    end subroutine addCustomDerivationsToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addSimpleDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        !! Add simple derivation to textbook based on JSON data

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addSimpleDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addPolyDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        !! Add a polynomial derivation to textbook based on JSON config file

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addPolyDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addAdditiveDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        !! Add an additive derivation to textbook based on JSON config file

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addAdditiveDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addMultiplicativeDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        !! Add a multiplicative derivation to textbook based on JSON config file

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

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

    end subroutine addBoundedExtDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addColdIonIJIntDerivationToTextbook(textbookObj,derivTag,gridObj,jsonCont,mpiCont)
        !! Add a cold ion I/J integral derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(Grid)              ,intent(in)    :: gridObj   
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addColdIonIJIntDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addIJIntDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        !! Add a I/J integral derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj   
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addIJIntDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addHarmonicExtractorDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        !! Add a harmonic extractor derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj   
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addHarmonicExtractorDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addFScalingDerivationToTextbook(textbookObj,derivTag,vSpaceObj,partObj,jsonCont,mpiCont)
        !! Add a distribution scaling extrapolation matrix derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj
        type(Partition)         ,intent(in)    :: partObj   
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addFScalingDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addDDVDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        !! Add a d/dv derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addDDVDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addD2DV2DerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        !! Add a d^2/dv^2 derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addD2DV2DerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addMomentDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        !! Add a custom moment derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addMomentDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addLocalValExtractorDerivationToTextbook(textbookObj,derivTag,partObj,jsonCont,mpiCont)
        !! Add a local fluid value extraction derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(Partition)         ,intent(in)    :: partObj
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addLocalValExtractorDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addVelContracDerivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        !! Add velocity space contraction derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
        
    end subroutine addVelContracDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addVelTProdDeivationToTextbook(textbookObj,derivTag,vSpaceObj,jsonCont,mpiCont)
        !! Add velocity vector tensor product derivation to textbook based on JSON config file
    
        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VSpace)            ,intent(in)    :: vSpaceObj
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 
        
    end subroutine addVelTProdDeivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addGenIntPolyDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        !! Add a generalized integer powered polynomial derivation to textbook based on JSON config file

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addGenIntPolyDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addRangeFilterDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        !! Add a range filter derivation to textbook based on JSON config file

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addRangeFilterDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addCalculationTreeDerivationToTextbook(textbookObj,derivTag,varList,jsonCont,mpiCont)
        !! Add a calculation tree derivation to textbook based on JSON config file

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(VariableList)      ,intent(in)    :: varList
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addCalculationTreeDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addnDLinInterpDerivationToTextbook(textbookObj,derivTag,jsonCont,mpiCont)
        !! Add an n-D linear interpolation derivation to textbook based on JSON config file

        type(Textbook)          ,intent(inout) :: textbookObj 
        character(*)            ,intent(in)    :: derivTag 
        type(JSONController)    ,intent(inout) :: jsonCont    !! JSONController used to get parameters from ./config.json 
        type(MPIController)     ,intent(inout) :: mpiCont     !! MPIController used with JSONController 

    end subroutine addnDLinInterpDerivationToTextbook
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module custom_derivation_support
!-----------------------------------------------------------------------------------------------------------------------------------