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
module janev_fits
    !! author: Stefan Mijin 
    !! 
    !! Contains interfaces to fit functions from Janev's EIRENE report: http://www.eirene.de/report_4105.pdf

    use data_kinds                  ,only: rk, ik
    use physical_constants
    use support_functions           ,only: expInt1
 
    implicit none 

!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function radRecombRateHydrogen(temp,n) result(rate)
    !! Radiative recombination rate given by formulas (21) and (22) from Janev Section 2.1.4. The recombination for n=2 is obtained
    !! by a simple sum over the two l states. Result in 10^-14 cm^3/s

        real(rk) ,dimension(:) ,intent(in)   :: temp !! Temperature array in eV
        integer(ik)            ,intent(in)   :: n !! Final atomic state after recombination
        real(rk) ,allocatable ,dimension(:)  :: rate

    end function radRecombRateHydrogen
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function excitationCrossSectionHydrogen(eGrid,n,m) result(cs)
    !! Calculates hydrogen electron impact excitation cross-section from state n to state m  on given energy grid.
    !! NOTE: The returned value is in 10^-16 cm^2!

        real(rk) ,dimension(:) ,intent(in)   :: eGrid !! Electron energy grid in eV
        integer(ik)            ,intent(in)   :: n !! Initial atomic state
        integer(ik)            ,intent(in)   :: m !! Final atomic state
        real(rk) ,allocatable ,dimension(:)  :: cs

    end function excitationCrossSectionHydrogen
!-----------------------------------------------------------------------------------------------------------------------------------
    pure module function ionizationCrossSectionHydrogen(eGrid,n) result(cs)
    !! Calculates hydrogen electron impact ionization cross-section from state n on given energy grid.
    !! NOTE: The returned value is in 10^-16 cm^2!

        real(rk) ,dimension(:) ,intent(in)   :: eGrid !! Electron energy grid in eV
        integer(ik)            ,intent(in)   :: n !! Initial atomic state
        real(rk) ,allocatable ,dimension(:)  :: cs

    end function ionizationCrossSectionHydrogen
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
end module janev_fits
!-----------------------------------------------------------------------------------------------------------------------------------