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
submodule (param_wrapper_1i1_derivation_class) param_wrapper_1i1_derivation_procedures
!! author: Stefan Mijin 
!! 
!! Contains module procedures associated with the FunWrapper1I1 derivation class

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initWrapper1I1(this,fun,param,multConst,multConstInner)
    !! Initialize 1I1 function wrapper derivation

    class(FunWrapperDerivation1I1)   ,intent(inout) :: this
    procedure(realArrayFunctionIntParam)            :: fun
    integer(ik)                      ,intent(in)    :: param
    real(rk) ,optional               ,intent(in)    :: multConst
    real(rk) ,optional               ,intent(in)    :: multConstInner

    this%multConst = real(1,kind=rk)
    if (present(multConst)) this%multConst = multConst 

    this%multConstInner = real(1,kind=rk)
    if (present(multConstInner)) this%multConstInner = multConstInner

    this%param = param

    this%funPtr => fun

    call this%makeDefined()

end subroutine initWrapper1I1  
!-----------------------------------------------------------------------------------------------------------------------------------
module function calculateWrapper1I1(this,inputArray,indices) result(output)

    class(FunWrapperDerivation1I1)     ,intent(inout)    :: this 
    type(RealArray)       ,dimension(:) ,intent(in)    :: inputArray 
    integer(ik)           ,dimension(:) ,intent(in)    :: indices           
    real(rk) ,allocatable ,dimension(:)                :: output

    if (assertions) then 

        call assertPure(this%isDefined(),"calculateWrapper1I1 called on undefined derivation object")
        call assertPure(size(indices) == 1,"Number of indices passed to calculateWrapper1I1 must be 1")
        call assertPure(all(indices>0),"indices passed to calculateWrapper1I1 out of bounds - lower")
        call assertPure(all(indices<=size(inputArray)),"indices passed to calculateWrapper1I1 out of bounds - upper")

    end if

    output = this%multConst * this%funPtr(this%multConstInner * inputArray(indices(1))%entry,this%param)

end function calculateWrapper1I1
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeWrapper1I1(this) 
    !! Deallocate pointer component

    type(FunWrapperDerivation1I1) ,intent(inout) :: this

    nullify(this%funPtr)

end subroutine finalizeWrapper1I1 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule param_wrapper_1i1_derivation_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
