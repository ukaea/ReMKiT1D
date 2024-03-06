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
submodule (cvode_integrator_class) cvode_integrator_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains module procedures associated with the CVODE integrator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
    module subroutine integrateCVODE(this,manipulatedModeller,outputVars,inputVars) 
        !! Implementation of abstract manipulate routine for the case of Runge-Kutta integrator

        class(CVODEIntegrator)                ,intent(inout) :: this 
        class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during integration
        class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the integration output 
        class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the integration routine

    end subroutine integrateCVODE
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initCVODEIntegrator(this,mpiCont,modelList,termGroups,evolvesTimeVar,dtController,initialTimestep)

        class(CVODEIntegrator)                    ,intent(inout) :: this
        type(MPIController)                       ,intent(in)    :: mpiCont !! 
        integer(ik)    ,optional    ,dimension(:) ,intent(in)    :: modelList !! List of models this integrator will be responsible for
        type(IntArray) ,optional    ,dimension(:) ,intent(in)    :: termGroups !! Term groups this integrator is responsible for - should conform with modelList
        logical        ,optional                  ,intent(in)    :: evolvesTimeVar !! Set to true if this integrator is allowed to change the "time" varible (if present in passed variable container)
        class(TimestepController) ,optional       ,intent(in)    :: dtController !! User-supplied timestep controller object
        real(rk)                  ,optional       ,intent(in)    :: initialTimestep !! Default timestep


        integer(c_int) :: ierr 

 !       ierr = FSUNContext_Create(mpiCont%getComm(),this%sunctx) 

        if (assertions .or. assertionLvl >= 0) then
            
            if (present(termGroups)) call assert(present(modelList),"Term groups object passed to RK integrator constructor without&
                                                                   & model list")
            if (present(modelList) .and. present(termGroups)) call assert(size(modelList)==size(termGroups),"If both model list and&
            & term groups arguments are passed to initRKIntegrator both must have the same length")
            if (present(dtController)) call assert(dtController%isDefined(),"Undefined timestep controller passed to explicit RK&
            & integrator constructor")
        end if

        if (present(evolvesTimeVar)) call this%setTimeEvolving(evolvesTimeVar) 
        if (present(modelList)) call this%setModelIndices(modelList)
        if (present(termGroups)) call this%setTermGroups(termGroups)
        if (present(initialTimestep)) call this%setTimestep(initialTimestep)
        if (present(dtController)) call this%setTimestepController(dtController)
        call this%setNonTrivialModelDataUpdate(.false.)
        call this%setNonTrivialUpdate(.false.)

        call this%makeDefined()


    end subroutine initCVODEIntegrator
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule cvode_integrator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
