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
module inelastic_mapping
    !! author: Stefan Mijin 
    !! 
    !! Contains procedures associated with the construction of velocity space mapping for inelastic collisional processes

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use v_space_class               ,only: VSpace 
    use sparse_row_data_class       ,only: SparseRowData
    use support_types               ,only: IntArray ,RealArray
    use support_functions           ,only: findIndices
 
    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function inelWeights(space,E) result(wMat)
        !! Calculate velocity space inelastic mapping weight matrix for transition with energy E on given velocity space. Assumes that 
        !! velocity and energy are suitably normalized (i.e. v^2 is an energy)

        type(VSpace)        ,intent(in) :: space !! Velocity space object to calculate the mapping on
        real(rk)            ,intent(in) :: E    !! Mapping transition energy
        type(SparseRowData)             :: wMat

    end function inelWeights
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module subroutine fillEmissionVector(wMat,emit) 
        !! Fills passed velocity space vector emit with zeros or ones, depending on whether there are any entries in the corresponding 
        !! column of passed mapping weight matrix wMat

        type(SparseRowData)          ,intent(in) :: wMat
        real(rk)    ,dimension(:) ,intent(inout) :: emit

    end subroutine fillEmissionVector
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module inelastic_mapping
!-----------------------------------------------------------------------------------------------------------------------------------