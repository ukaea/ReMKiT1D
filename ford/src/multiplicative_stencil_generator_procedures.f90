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
submodule (multiplicative_stencil_generator_class) multiplicative_stencil_generator_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the multiplicative stencil generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initMultGen(this,coreObj,fluidCol,initXVals,initHVals,initVVals) 
    !! Multiplicative stencil value generator initialization routine

    class(MultiplicativeStencilGen)                ,intent(inout)  :: this
    type(MultiplicativeGeneratorCore) ,allocatable ,intent(inout)  :: coreObj !! Multiplicative core - will be deallocated after call
    logical ,optional                              ,intent(in)     :: fluidCol !! True if the column variable for this stencil is fluid. 
    type(RealArray) ,optional ,dimension(:)        ,intent(in)     :: initXVals !! Optional initial raw x stencil values. Defaults to unallocated. 
    type(RealArray) ,optional ,dimension(:)        ,intent(in)     :: initHVals !! Optional initial raw h stencil values. Defaults to unallocated. 
    type(RealArray) ,optional ,dimension(:)        ,intent(in)     :: initVVals !! Optional initial raw v stencil values. Defaults to unallocated. 

    if (assertions) &
    call assert(coreObj%isDefined(),"Undefined multiplicative core passed to multiplicative stencil generator constructor")

    call move_alloc(coreObj,this%core)

    if (present(initXVals)) allocate(this%xVals,source=initXVals)
    if (present(initHVals)) allocate(this%hVals,source=initHVals)
    if (present(initVVals)) allocate(this%vVals,source=initVVals)

    this%fluidCol = .false.
    if (present(fluidCol)) this%fluidCol = fluidCol 

    call this%makeDefined()

end subroutine initMultGen
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine calcMultVals (this,varCont,res,mbData,hostModel)
    !! Calculate multiplicative stencil values in place

    class(MultiplicativeStencilGen)             ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)  ,intent(inout) :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    if (assertions) then
        call assert(this%isDefined(),"calcMultVals called from undefined multiplicative stencil generator")
        call assert(allocated(this%xVals) .or. allocated(this%xStencilGen), &
        "calcMultVals called when neither xVals nor xStencilGen have been allocated")

        if (.not. this%fluidCol) then 
            call assert(allocated(this%hVals) .or. allocated(this%hStencilGen), &
            "calcMultVals called in kinetic column case when neither hVals nor hStencilGen have been allocated")
            call assert(allocated(this%vVals) .or. allocated(this%vStencilGen), &
            "calcMultVals called in kinetic column case when neither vVals nor vStencilGen have been allocated")
        end if
    end if
    
    if (allocated(this%xStencilGen)) then 
        if (present(mbData)) then 
            if (present(hostModel)) then
                call this%xStencilGen%calculateInPlace(varCont,this%xVals,mbData,hostModel)
            else
                call this%xStencilGen%calculateInPlace(varCont,this%xVals,mbData)
            end if
        else
            if (present(hostModel)) then
                call this%xStencilGen%calculateInPlace(varCont,this%xVals,hostModel=hostModel)
            else
                call this%xStencilGen%calculateInPlace(varCont,this%xVals)
            end if
        end if

    end if

    if (.not. this%fluidCol) then  

        if (allocated(this%hStencilGen)) then 
            if (present(mbData)) then 
                if (present(hostModel)) then
                    call this%hStencilGen%calculateInPlace(varCont,this%hVals,mbData,hostModel)
                else
                    call this%hStencilGen%calculateInPlace(varCont,this%hVals,mbData)
                end if
            else
                if (present(hostModel)) then
                    call this%hStencilGen%calculateInPlace(varCont,this%hVals,hostModel=hostModel)
                else
                    call this%hStencilGen%calculateInPlace(varCont,this%hVals)
                end if
            end if

        end if

        if (allocated(this%vStencilGen)) then 
            if (present(mbData)) then 
                if (present(hostModel)) then
                    call this%vStencilGen%calculateInPlace(varCont,this%vVals,mbData,hostModel)
                else
                    call this%vStencilGen%calculateInPlace(varCont,this%vVals,mbData)
                end if
            else
                if (present(hostModel)) then
                    call this%vStencilGen%calculateInPlace(varCont,this%vVals,hostModel=hostModel)
                else
                    call this%vStencilGen%calculateInPlace(varCont,this%vVals)
                end if
            end if

        end if

    end if

    !Calculate stencil values using core
    
    if (this%fluidCol) then 
        call this%core%calculateInPlace(res,this%xVals)
    else
        call this%core%calculateInPlace(res,this%xVals,this%hVals,this%vVals)
    end if
end subroutine calcMultVals 
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setXGen (this,gen)
    !! Setter for xStencilGen. 

    class(MultiplicativeStencilGen)               ,intent(inout) :: this
    class(JaggedArrayGenerator)                   ,intent(in)    :: gen

    if (assertions) then 
        call assert(this%isDefined(),"setXGen called on undefined multiplicative stencil generator")
        call assert(gen%isDefined(),"Undefined stencil value generator passed to setXGen")
    end if

    if (allocated(this%xStencilGen)) deallocate(this%xStencilGen)
    allocate(this%xStencilGen,source=gen)

end subroutine setXGen 
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setHGen (this,gen)
    !! Setter for hStencilGen. 

    class(MultiplicativeStencilGen)               ,intent(inout) :: this
    class(JaggedArrayGenerator)                   ,intent(in)    :: gen

    if (assertions) then 
        call assert(this%isDefined(),"setHGen called on undefined multiplicative stencil generator")
        call assert(gen%isDefined(),"Undefined stencil value generator passed to setHGen")
    end if

    if (allocated(this%hStencilGen)) deallocate(this%hStencilGen)
    allocate(this%hStencilGen,source=gen)

end subroutine setHGen 
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine setVGen (this,gen)
    !! Setter for vStencilGen. 

    class(MultiplicativeStencilGen)               ,intent(inout) :: this
    class(JaggedArrayGenerator)                   ,intent(in)    :: gen

    if (assertions) then 
        call assert(this%isDefined(),"setVGen called on undefined multiplicative stencil generator")
        call assert(gen%isDefined(),"Undefined stencil value generator passed to setVGen")
    end if

    if (allocated(this%vStencilGen)) deallocate(this%vStencilGen)
    allocate(this%vStencilGen,source=gen)

end subroutine setVGen 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule multiplicative_stencil_generator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
