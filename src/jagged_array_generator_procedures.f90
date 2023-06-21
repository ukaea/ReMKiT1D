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
submodule (jagged_array_generator_class) jagged_array_generator_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the abstract jagged array generator class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculate(this,varCont,mbData,hostModel) result(res)
    !! Use in place version of stencil calculation to return values

    class(JaggedArrayGenerator)                ,intent(inout) :: this
    type(VariableContainer)                     ,intent(in)    :: varCont
    type(RealArray) ,allocatable ,dimension(:)                 :: res
    class(ModelboundData) ,optional             ,intent(in)    :: mbData
    class(ModelSurrogate) ,optional             ,intent(in)    :: hostModel

    if (present(mbData)) then
        if (present(hostModel)) then
            call this%calculateInPlace(varCont,res,mbData,hostModel)
        else
            call this%calculateInPlace(varCont,res,mbData)
        end if
    else
        if (present(hostModel)) then 
            call this%calculateInPlace(varCont,res,hostModel=hostModel)
        else
            call this%calculateInPlace(varCont,res)
        end if
    end if

end function calculate
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule jagged_array_generator_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
