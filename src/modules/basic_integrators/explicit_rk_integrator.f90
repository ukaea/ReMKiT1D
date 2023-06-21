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
module explicit_rk_integrator_class
    !! author: Stefan Mijin
    !! 
    !! Contains an explicit Runge-Kutta integrator class with the option to supply user Butcher tableaus

    use data_kinds                         ,only: rk ,ik
    use god_objects                        ,only: Object
    use runtime_constants                  ,only: debugging, assertions
    use assertion_utility                  ,only: assert, assertIdentical, assertPure
    use modeller_surrogate_class           ,only: ModellerSurrogate
    use variable_container_class           ,only: VariableContainer
    use modeller_class                     ,only: Modeller 
    use integrator_abstract_class          ,only: Integrator
    use support_types                      ,only: IntArray ,RealArray ,LogicalArray
    use mpi_controller_class               ,only: CommunicationData
    use timestep_controller_abstract_class ,only: TimestepController

    implicit none
    private

    type ,public :: BTableau
        !! Butcher tableau derived type 

        type(RealArray) ,allocatable ,dimension(:) ,public :: a !! Lower triangular matrix of coefficients
        real(rk)        ,allocatable ,dimension(:) ,public :: b !! Coefficients for summing over intermediate values
        real(rk)        ,allocatable ,dimension(:) ,public :: c !! Time coefficients

    end type

    type ,public ,extends(Integrator) :: ExplicitRKIntegrator
        !! Explicit Runge-Kutta integrator

        type(BTableau)                               ,private :: tableau !! Butcher tableau used by the integrator
        type(IntArray)  ,allocatable ,dimension(:)   ,private :: tableauANonzeros !! Locations of non-zeros in the  a_ij tableau matrix
        type(RealArray) ,allocatable ,dimension(:,:) ,private :: intermediateValues !! Buffer for intermediate RK values for various variables
        type(VariableContainer) ,allocatable         ,private :: buffer !! VariableContainer buffer for passing to Modeller routines

        contains

        procedure ,public :: affect => integrateRK

        procedure ,public :: init => initRKIntegrator

    end type ExplicitRKIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    interface
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrateRK(this,manipulatedModeller,outputVars,inputVars) 
        !! Implementation of abstract manipulate routine for the case of Runge-Kutta integrator

        class(ExplicitRKIntegrator)           ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    end subroutine integrateRK
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine initRKIntegrator(this,modelList,termGroups,order,tableau,evolvesTimeVar,dtController,initialTimestep)
        !! Constructs an RK integrator with initial evaluated model indices and correspoding term groups if those
        !! are present. If neither order or tableau are provided defaults to forward Euler time-stepping. Possible order values are 1,2,3,4
        !! where 2 is the midpoint method, 3 is SSPRK3, and 4 is the standard RK4 method. 

        implicit none 

        class(ExplicitRKIntegrator)               ,intent(inout) :: this
        integer(ik)    ,optional    ,dimension(:) ,intent(in)    :: modelList !! List of models this integrator will be responsible for
        type(IntArray) ,optional    ,dimension(:) ,intent(in)    :: termGroups !! Term groups this integrator is responsible for - should conform with modelList
        integer(ik)    ,optional                  ,intent(in)    :: order !! Runge-Kutta order
        type(BTableau) ,optional                  ,intent(in)    :: tableau !! User-defined Butcher tableau
        logical        ,optional                  ,intent(in)    :: evolvesTimeVar !! Set to true if this integrator is allowed to change the "time" varible (if present in passed variable container)
        class(TimestepController) ,optional       ,intent(in)    :: dtController !! User-supplied timestep controller object
        real(rk)                  ,optional       ,intent(in)    :: initialTimestep !! Default timestep

    end subroutine initRKIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface 
!-----------------------------------------------------------------------------------------------------------------------------------
 end module explicit_rk_integrator_class
!-----------------------------------------------------------------------------------------------------------------------------------
 