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
module janev_crm_data_support
    !! author: Stefan Mijin 
    !! 
    !! Contains functions for adding particular processes with Janev fits to modelbound CRM data

    use data_kinds                         ,only: rk, ik
    use runtime_constants                  ,only: debugging, assertions
    use assertion_utility                  ,only: assert, assertIdentical, assertPure
    use modelbound_CRM_data_class          ,only: ModelboundCRMData
    use fixed_ecs_transition_class         ,only: FixedECSTransition
    use db_transition_class                ,only: DBTransition
    use derived_transition_class           ,only: DerivedTransition
    use simple_derivation_class            ,only: SimpleDerivation
    use param_wrapper_1i1_derivation_class ,only: FunWrapperDerivation1I1
    use v_space_class                      ,only: VSpace
    use janev_fits
    use physical_constants

    implicit none 

    public

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addJanevRadRecombTransition(crmData,locNumX,endState,temperatureVarIndex,eVTempNorm,densNorm,timeNorm)
        !! Adds Janev radiative recombination transition to endState of hydrogen. Normalizes the rate to densNorm/timeNorm for 
        !! easier use.
        
        type(ModelboundCRMData) ,intent(inout) :: crmData 
        integer(ik)             ,intent(in)    :: locNumX 
        integer(ik)             ,intent(in)    :: endState 
        integer(ik)             ,intent(in)    :: temperatureVarIndex
        real(rk)                ,intent(in)    :: eVTempNorm
        real(rk)                ,intent(in)    :: densNorm
        real(rk)                ,intent(in)    :: timeNorm 

    end subroutine addJanevRadRecombTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addJanevCollExIonTransition(crmData,locNumX,startState,endState,distVarIndex,&
                                                  eVTempNorm,csNorm,refVSpace,fixedWIndex)
        !! Adds Janev collisional excitation/ionization transition from startState to endState of hydrogen.
        !! Normalizes the cross-section for easier later use. Cross-section only isotropic.
        
        type(ModelboundCRMData) ,intent(inout) :: crmData 
        integer(ik)             ,intent(in)    :: locNumX 
        integer(ik)             ,intent(in)    :: startState 
        integer(ik)             ,intent(in)    :: endState 
        integer(ik)             ,intent(in)    :: distVarIndex
        real(rk)                ,intent(in)    :: eVTempNorm
        real(rk)                ,intent(in)    :: csNorm
        type(VSpace) ,target    ,intent(inout) :: refVSpace !! Target for the transition reference pointer
        integer(ik)             ,intent(in)    :: fixedWIndex !! Index of the inelastic weight array corresponding to this transition

    end subroutine addJanevCollExIonTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine addJanevCollDeexRecombTransition(crmData,locNumX,startState,endState,distVarIndex,&
                                                      tempVarIndex,eVTempNorm,densNorm,refVSpace,&
                                                      directFixedWIndex,fixedWIndex,directTransitionIndex,csUpdatePriority)
        !! Adds Janev collisional deexcitation/recombination transition from startState to endState of hydrogen.
        !! Normalizes the cross-section for easier later use. Cross-section only isotropic.

        type(ModelboundCRMData) ,intent(inout) :: crmData 
        integer(ik)             ,intent(in)    :: locNumX 
        integer(ik)             ,intent(in)    :: startState 
        integer(ik)             ,intent(in)    :: endState 
        integer(ik)             ,intent(in)    :: distVarIndex
        integer(ik)             ,intent(in)    :: tempVarIndex
        real(rk)                ,intent(in)    :: eVTempNorm
        real(rk)                ,intent(in)    :: densNorm
        type(VSpace) ,target    ,intent(inout) :: refVSpace !! Target for the transition reference pointer
        integer(ik)             ,intent(in)    :: directFixedWIndex !! Index of the inelastic weight array corresponding to 
                                                                    !! the direct counterpart of this transition
        integer(ik)             ,intent(in)    :: fixedWIndex !! Index of the inelastic weight array corresponding to this transition
        
        integer(ik)             ,intent(in)    :: directTransitionIndex !! Index of the direct transition in crmData
        integer(ik)             ,intent(in)    :: csUpdatePriority !! Update priority for detailed balance cross-section

end subroutine addJanevCollDeexRecombTransition
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module janev_crm_data_support
!-----------------------------------------------------------------------------------------------------------------------------------