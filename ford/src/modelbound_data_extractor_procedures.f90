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
submodule (modelbound_data_extractor_class) modelbound_data_extractor_procedures
!! author: Stefan Mijin 
!! 
!! Contains procedures associated with the ModelboundDataExtractor class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine initExtractor(this,resultVarIndex,modelIndex,modelboundDataName) 
    !! ModelboundDataExtractor initialization routine

    class(ModelboundDataExtractor) ,intent(inout)  :: this
    integer(ik)                    ,intent(in)     :: resultVarIndex !! Index of variable to write the result in
    integer(ik)                    ,intent(in)     :: modelIndex     !! Index of model housing required modelbound data
    character(*)                   ,intent(in)     :: modelboundDataName !! Name of data to extract

    this%resultVarIndex = resultVarIndex 
    this%modelIndex = modelIndex 
    this%modelboundDataName = modelboundDataName

    call this%makeDefined()

end subroutine initExtractor
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine extract(this,manipulatedModeller,outputVars,inputVars) 
    !! Implementation of abstract manipulate routine for the Extractor

    class(ModelboundDataExtractor)        ,intent(inout) :: this 
    class(ModellerSurrogate)              ,intent(inout) :: manipulatedModeller !! Modeller to be used in callbacks during manipulation
    class(VariableContainer)              ,intent(inout) :: outputVars !! VariableContainer object to store the manipulation output 
    class(VariableContainer)              ,intent(in)    :: inputVars !! VariableContainer object housing input data for the manipulation routine

    real(rk) ,allocatable ,dimension(:) :: modelboundDataBuffer 

    integer(ik) :: inferredHaloDiff ,lboundVar

     if (assertions) then 
        call assert(this%isDefined(),"extract routine called by undefined ModelboundDataExtractor")
        call assert(manipulatedModeller%isDefined(),"Undefined Modeller object passed to ModelboundDataExtractor extract routine")
        call assert(outputVars%isDefined(),"outputVars passed to extract not defined")
        call assert(inputVars%isDefined(),"inputVars passed to extract not defined")
    end if

    select type (manipulatedModeller)
    type is (Modeller)
  
        ! copy modelbound data to make sure it is of the correct size 
        allocate(modelboundDataBuffer,mold=outputVars%variables(this%resultVarIndex)%entry)
        call manipulatedModeller%copyDataFromModel(this%modelIndex,this%modelboundDataName,modelboundDataBuffer)

        if (assertions) call assert(size(modelboundDataBuffer) <= size(outputVars%variables(this%resultVarIndex)%entry),&
                                   "modelbound data variable requested from ModelboundDataExtractor does not fit into result var") 

        ! Handle case where modelbound data does not have the same halo size as the result variable
        inferredHaloDiff = (size(outputVars%variables(this%resultVarIndex)%entry) - size(modelboundDataBuffer))/2
        lboundVar = lbound(outputVars%variables(this%resultVarIndex)%entry,1)
        outputVars%variables(this%resultVarIndex)%entry(lboundVar+inferredHaloDiff:&
                                                        lboundVar+inferredHaloDiff + size(modelboundDataBuffer) - 1) = &
                                                        modelboundDataBuffer
    class default
        error stop "Unsupported surrogate passed to ModelboundDataExtractor"
    end select

end subroutine extract
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule modelbound_data_extractor_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
