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
module coulomb_log_derivation_class
    !! author: Stefan Mijin 
    !!
    !! Houses derivation returning coulomb log values based on density and temperature (and optionally ionization)

    use data_kinds                   ,only: rk ,ik
    use runtime_constants            ,only: debugging, assertions
    use assertion_utility            ,only: assert, assertIdentical, assertPure
    use god_objects                  ,only: Object
    use support_types                ,only: RealArray 
    use derivation_abstract_class    ,only: Derivation
    use physics_functions

    implicit none
    private

    type ,public ,extends(Derivation) :: CoulombLogDerivation
        !! Returns values of the coulomb log based on electron density and temperature variables (and optionally ionization). 
        !! The passed indices are assumed to be in T,n,Z order if ee or ei logarithm, or in n1,n2,T1,T2 order if ii logarithm. 
        !! Assumes that density is normalized to m^{-3} units, and temperature to eV.

        logical  ,private :: electronLog !! If true returns the electron-electron Coulomb log instead of the electron-ion log. Default is false
        logical  ,private :: ionLog !! If true returns the ion-ion Coulomb log instead of the electron-ion log. Default is false
        real(rk) ,private :: ionZ !! The ion Z used if no ionization variable is passed, or if calculating ion-ion logarithm
        real(rk) ,private :: ionZ2 !! The second ion Z used if calculating ion-ion logarithm (Defaults to ionZ)
        real(rk) ,private :: ionMassRatio !! Mass ratio used if calculating ion-ion logarithm (Defaults to 1)

        integer(ik) ,private :: locNumX !! Local number of cells (used to avoid dividing by zero in ghost/halo cells)

        real(rk) ,private :: densNorm !! Density normalization
        real(rk) ,private :: tempNorm !! Temperature normalization

        contains

        procedure ,public :: init => initCoulombLogDeriv

        procedure ,public :: calculate => calculateCoulombLog

    end type CoulombLogDerivation
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initCoulombLogDeriv(this,ionZ,locNumX,densNorm,tempNorm,electronLog,ionLog,ionZ2,ionMassRatio)
        !! Initialize Coulomb Log derivation object

        class(CoulombLogDerivation)   ,intent(inout) :: this
        real(rk)                      ,intent(in)    :: ionZ 
        integer(ik)                   ,intent(in)    :: locNumX
        real(rk)                      ,intent(in)    :: densNorm
        real(rk)                      ,intent(in)    :: tempNorm
        logical ,optional             ,intent(in)    :: electronLog
        logical ,optional             ,intent(in)    :: ionLog
        real(rk) ,optional            ,intent(in)    :: ionZ2
        real(rk) ,optional            ,intent(in)    :: ionMassRatio

    end subroutine initCoulombLogDeriv  
!-----------------------------------------------------------------------------------------------------------------------------------
    module function calculateCoulombLog(this,inputArray,indices) result(output)

        class(CoulombLogDerivation)        ,intent(inout)    :: this 
        type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
        integer(ik)           ,dimension(:) ,intent(in)    :: indices           
        real(rk) ,allocatable ,dimension(:)                :: output

    end function calculateCoulombLog
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module coulomb_log_derivation_class
!-----------------------------------------------------------------------------------------------------------------------------------
 