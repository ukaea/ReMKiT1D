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
submodule (simple_timestep_controller_class) simple_timestep_controller_procedures
!! author: Stefan Mijin
!! 
!! Contains module procedures associated with the simple timestep controller class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
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

    integer(ik) :: i

    if (assertions .or. assertionLvl >= 0) then 
        call assert(mpiCont%isDefined(),"Undefined MPI controller passed to initSimppleTimestepController")
        call assert(varCont%isDefined(),"Undefined reference variable container passed to initSimppleTimestepController")
        call assert(size(reqVars) == size(reqVarPowers), "reqVars and reqVarPowers passed to initSimpleTimestepController &
        &do not conform in size")
    end if

    this%multConst = real(1,kind=rk)
    this%useMaxVal = .false. 
    this%rescaleTimestep = .false. 

    if (present(multConst)) this%multConst = multConst
    if (present(useMaxVal)) this%useMaxVal = useMaxVal
    if (present(rescaleTimestep)) this%rescaleTimestep = rescaleTimestep

    allocate(this%reqVarIndices(size(reqVars)))

    do i = 1, size(reqVars)
        this%reqVarIndices(i) = varCont%getVarIndex(reqVars(i)%string)
    end do

    this%reqVarPowers = reqVarPowers

    this%mpiCont = mpiCont

    call this%makeDefined()
    
end subroutine initSimpleTimestepController
!-----------------------------------------------------------------------------------------------------------------------------------
module function evaluateTimestepSimple(this,inputVars,currentTimestep) result(timestep)
    !! Calculate timestep as min(max) of variable product

    class(SimpleTimestepController)       ,intent(inout) :: this 
    class(VariableContainer)              ,intent(in)    :: inputVars !! Variable container used to calculate timestep
    real(rk)                              ,intent(in)    :: currentTimestep !! Current timestep to be used if the controller rescales timesteps
    real(rk)                                             :: timestep

    real(rk) ,allocatable ,dimension(:) :: localTimestep

    real(rk) :: localReduce

    integer(ik) :: i ,inferredHaloSize ,upperBound

    inferredHaloSize = 1 - lbound(inputVars%variables(this%reqVarIndices(1))%entry,1)
    upperBound = ubound(inputVars%variables(this%reqVarIndices(1))%entry,1) - inferredHaloSize

    localTimestep = this%multConst * inputVars%variables(this%reqVarIndices(1))%entry(1:upperBound) ** this%reqVarPowers(1)

    do i = 2, size(this%reqVarIndices)
        localTimestep = localTimestep * inputVars%variables(this%reqVarIndices(i))%entry(1:upperBound) ** this%reqVarPowers(i)
    end do
    
    if (this%rescaleTimestep) localTimestep = localTimestep * currentTimestep

    if (this%useMaxVal) then 
        localReduce = maxval(localTimestep)
        timestep = this%mpiCont%allreduceMax(localReduce)
    else 
        localReduce = minval(localTimestep)
        timestep = this%mpiCont%allreduceMin(localReduce)
    end if
end function evaluateTimestepSimple
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule simple_timestep_controller_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
