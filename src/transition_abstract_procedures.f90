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
submodule (transition_abstract_class) transition_abstract_procedures
    !! author: Stefan Mijin
    !!
    !! Contains abstract transition procedures

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getIngoingStates(this) result(inStates)
    !! Getter for ingoingStates

    class(transition)          ,intent(in)  :: this
    integer(ik) ,allocatable ,dimension(:)  :: inStates

    if (assertions) call assertPure(this%isDefined(),"getIngoingStates call on unidentified transition object")

    inStates = this%ingoingStates

end function getIngoingStates  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getOutgoingStates(this) result(outStates)
    !! Getter for outgoingStates

    class(transition)          ,intent(in)  :: this
    integer(ik) ,allocatable ,dimension(:)  :: outStates

    if (assertions) call assertPure(this%isDefined(),"getOutgoingStates call on unidentified transition object")

    outStates = this%outgoingStates

end function getOutgoingStates  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setStates(this,inStates,outStates)
    !! Setter for both ingoing and outgoing states

    class(Transition)          ,intent(inout)  :: this
    integer(ik)  ,dimension(:) ,intent(in)     :: inStates !! Ingoing state IDs
    integer(ik)  ,dimension(:) ,intent(in)     :: outStates !! Outgoing state IDs

    this%outgoingStates = outStates
    this%ingoingStates = inStates 

end subroutine setStates  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setRate(this,rate)
    !! Setter for rate values

    class(transition)        ,intent(inout)  :: this
    real(rk)  ,dimension(:)  ,intent(in)     :: rate

    this%rate = rate

end subroutine setRate  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRate(this) result(rate)
    !! Getter for rate values

    class(transition)          ,intent(in)  :: this
    real(rk)    ,allocatable ,dimension(:)  :: rate

    if (assertions) call assertPure(this%isDefined(),"getRate called for undefined transition object")

    rate = this%rate

end function getRate  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setRateMomentum(this,rate)
    !! Setter for rateMomentum values

    class(transition)        ,intent(inout)  :: this
    real(rk)  ,dimension(:)  ,intent(in)     :: rate

    this%rateMomentum = rate

end subroutine setRateMomentum  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRateMomentum(this) result(rate)
    !! Getter for rateMomentum values

    class(transition)          ,intent(in)  :: this
    real(rk)    ,allocatable ,dimension(:)  :: rate

    if (assertions) call assertPure(this%isDefined(),"getRateMomentum called for undefined transition object")

    rate = this%rateMomentum

end function getRateMomentum  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setRateEnergy(this,rate)
    !! Setter for rateEnergy

    class(transition)        ,intent(inout)  :: this
    real(rk)  ,dimension(:)  ,intent(in)     :: rate

    this%rateEnergy = rate

end subroutine setRateEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRateEnergy(this) result(rate)
    !! Getter for rateEnergy

    class(transition)          ,intent(in)  :: this
    real(rk)    ,allocatable ,dimension(:)  :: rate

    if (assertions) call assertPure(this%isDefined(),"getRateEnergy called for undefined transition object")

    rate = this%rateEnergy

end function getRateEnergy  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setCrossSectionCol(this,crossSection,col)
    !! Set cross-section values in column col   

    class(Transition)        ,intent(inout)  :: this
    real(rk)  ,dimension(:)  ,intent(in)     :: crossSection
    integer(ik)              ,intent(in)     :: col

    this%crossSection(:,col) = crossSection

end subroutine setCrossSectionCol  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setCrossSection(this,crossSection)
    !! Setter for crossSection

    class(Transition)         ,intent(inout)  :: this
    real(rk)  ,dimension(:,:) ,intent(in)     :: crossSection

    this%crossSection = crossSection

end subroutine setCrossSection  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCrossSectionCol(this,col) result(crossSection)
    !! Get cross-section values from column col

    class(transition)          ,intent(in)  :: this
    integer(ik)              ,intent(in)    :: col
    real(rk)    ,allocatable ,dimension(:)  :: crossSection

    if (assertions) call assertPure(this%isDefined(),"getCrossSectionCol called for undefined transition object")

    if (col > size(this%crossSection,2)) then 
        allocate(crossSection(size(this%crossSection,1)))
        crossSection = 0
    else 
        crossSection = this%crossSection(:,col)
    end if

end function getCrossSectionCol  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function includesElDensity(this) result(includesDens)
    !! Check whether rates in this transition include an electron density factor

    class(Transition)        ,intent(in)  :: this
    logical                               :: includesDens

    if (assertions) call assertPure(this%isDefined(),"includesElDensity called for undefined transition object")

    includesDens = this%rateContainsElDensity

end function includesElDensity  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setIncludeElectronDensity(this,includeDens)
    !! Setter for rateContainsElDensity

    class(Transition)    ,intent(inout)  :: this
    logical              ,intent(in)     :: includeDens

    this%rateContainsElDensity = includeDens

end subroutine setIncludeElectronDensity  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCSDim(this) result(csDim)
    !! Getter for csDim

    class(Transition)        ,intent(in)  :: this
    integer(ik)                           :: csDim

    if (assertions) call assertPure(this%isDefined(),"getCSDim called for undefined transition object")

    csDim = this%csDim

end function getCSDim  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setCSDim(this,csDim)
    !! Setter for csDim

    class(Transition)    ,intent(inout)  :: this
    integer(ik)          ,intent(in)     :: csDim

    this%csDim = csDim

end subroutine setCSDim  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getRateSize(this) result(rateSize)
    !! Getter for rate array length

    class(Transition)        ,intent(in)  :: this
    integer(ik)                           :: rateSize

    if (assertions) call assertPure(this%isDefined(),"getRateSize called for undefined transition object")

    rateSize = size(this%rate)

end function getRateSize  
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine noUpdate(this,varCont,hostModel,hostData,updatePriority)
    !! Default update routine - does nothing

    class(Transition)               ,intent(inout)  :: this
    type(VariableContainer)         ,intent(in)     :: varCont !! Variable container used in update
    class(ModelSurrogate) ,optional ,intent(in)     :: hostModel !! Optional host model reference for callbacks during update
    class(ModelboundData) ,optional ,intent(in)     :: hostData !! Optional host data reference for callbacks during update
    integer(ik) ,optional           ,intent(in)     :: updatePriority !! Priority for this update call 

end subroutine noUpdate  
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule transition_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
