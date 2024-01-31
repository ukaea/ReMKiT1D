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
module simple_timestep_controller_class
    !! author: Stefan Mijin 
    !! 
    !! Houses simple timestep controller which uses the min (max) value of a product of variables, scaled to corresponding powers

    use data_kinds                         ,only: rk ,ik
    use runtime_constants                  ,only: debugging, assertions, assertionLvl
    use god_objects                        ,only: Object
    use assertion_utility                  ,only: assert, assertIdentical, assertPure
    use variable_container_class           ,only: VariableContainer
    use mpi_controller_class               ,only: MPIController
    use timestep_controller_abstract_class ,only: TimestepController
    use support_types

    implicit none
    private

    type ,public ,extends(TimestepController) :: SimpleTimestepController
        !! Timestep controller that calculates the timestep using the min (max) value of a product of variables raised to corresponding
        !! powers. Variables are assumed fluid-like.

        type(MPIController)                    ,private :: mpiCont !! Local copy of MPI controller used for allreduce calls
        integer(ik) ,allocatable ,dimension(:) ,private :: reqVarIndices !! Indices of variables to be raised to a power and multiplied
        real(rk)    ,allocatable ,dimension(:) ,private :: reqVarPowers !! Powers corresponding to required variables 
        real(rk)                               ,private :: multConst !! Optional multiplicative const 
        logical                                ,private :: useMaxVal !! Use max value of product = false by default
        logical                                ,private :: rescaleTimestep !! Multiply the calculated timestep with currentTimestep dummy argument in evaluateTimestep

        contains

        procedure ,public :: init => initSimpleTimestepController

        procedure ,public :: evaluateTimestep => evaluateTimestepSimple

    end type SimpleTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine initSimpleTimestepController(this,mpiCont,varCont,reqVars,reqVarPowers&
                                                    ,multConst,useMaxVal,rescaleTimestep) 
            !! Simple timestep controller initialization routine 
        
            class(SimpleTimestepController)           ,intent(inout)  :: this
            type(MPIController)                       ,intent(in)     :: mpiCont !! Reference MPI controller
            type(VariableContainer)                   ,intent(in)     :: varCont !! Reference variable container
            type(StringArray)           ,dimension(:) ,intent(in)     :: reqVars !! Required variable list 
            real(rk)                    ,dimension(:) ,intent(in)     :: reqVarPowers !! Powers corresponding to required variables
            real(rk)          ,optional               ,intent(in)     :: multConst  !! Normalization constant
            logical           ,optional               ,intent(in)     :: useMaxVal !! True if max value of product is used instead of min
            logical           ,optional               ,intent(in)     :: rescaleTimestep !! True if currentTimestep is rescaled using the calculated timestep value

        end subroutine initSimpleTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
        module function evaluateTimestepSimple(this,inputVars,currentTimestep) result(timestep)
            !! Calculate timestep as min(max) of variable product

            class(SimpleTimestepController)       ,intent(inout) :: this 
            class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate timestep
            real(rk)                              ,intent(in)    :: currentTimestep !! Current timestep to be used if the controller rescales timesteps
            real(rk)                                             :: timestep

        end function evaluateTimestepSimple
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module simple_timestep_controller_class
!-----------------------------------------------------------------------------------------------------------------------------------
 