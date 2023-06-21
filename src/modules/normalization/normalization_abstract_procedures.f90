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
submodule (normalization_abstract_class) normalization_abstract_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains abstract normalization procedures 

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module subroutine setNormalizationVals(this,normVals)
    !! Set normalization values array of normalization object

    class(Normalization)          ,intent(inout)  :: this
    type(NamedReal) ,dimension(:) ,intent(in)     :: normVals

    if (allocated(this%normalizationVals)) deallocate(this%normalizationVals)
    allocate(this%normalizationVals,source=normVals)

end subroutine setNormalizationVals  
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getNormalizationValue(this,name) result(val)
    !! Get normalization value associated with passed name

    class(Normalization)      ,intent(in) :: this
    character(*)              ,intent(in) :: name
    real(rk)                              :: val

    integer(ik) :: i 

    logical :: found

    if (assertions) call assertPure(this%isDefined(),"Normalization value requested from undefined normalization object")

    found = .false. 

    do i = 1, size(this%normalizationVals)
        if (this%normalizationVals(i)%name == name) then 
            found = .true. 
            val = this%normalizationVals(i)%value
            exit
        end if
    end do

    if (.not. found) error stop "getNormalizationValue called for name not in normalization values"

end function getNormalizationValue
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function getCustomNormalization(this,names,powers,multConst) result(val)
    !! Get normalization value calculated as the product of normalization values with passed names raised to corresponding 
    !! powers and optionally multiplied by a constant

    class(Normalization)              ,intent(in) :: this
    type(StringArray) ,dimension(:)   ,intent(in) :: names
    real(rk)          ,dimension(:)   ,intent(in) :: powers
    real(rk) ,optional                ,intent(in) :: multConst
    real(rk)                                      :: val

    integer(ik) :: i 

    if (assertions) call assertPure(this%isDefined(),"Custom normalization value requested from undefined normalization object")

    val = real(1,kind=rk)
    if (present(multConst)) val = multConst 

    do i = 1, size(names)
        val = val * this%getNormalizationValue(names(i)%string)**powers(i)
    end do
    
end function getCustomNormalization 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule normalization_abstract_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
