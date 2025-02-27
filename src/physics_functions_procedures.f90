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
submodule (physics_functions) physics_functions_procedures
    !! author: Stefan Mijin 
    !! 
    !!  Contains the implementations of various physics routines

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function logLee(Te,ne) result(res)
    !! Calculate Coulomb logarithm for electron-electron collisions (NRL Formulary 2013 page 34 equation a)

    real(rk) ,intent(in) :: Te !! Electron temperature in eV
    real(rk) ,intent(in) :: ne !! Electron density in  m^{-3} 
    real(rk)             :: res

    res = real(23.5d0 - log(sqrt(ne * 1.00D-6) * Te**(-5.0d0/4.0d0)) - sqrt(1.0d-5 + (log(Te)-2.0d0)**2/16.0d0),kind=rk)

end function logLee
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function logLei(Te,ne,Z,removeDisc) result(res)
    !! Calculate Coulomb logarithm for electron-ion collisions (NRL Formulary 2013 page 34 equation b). Assumes Te > TiZm_e/m_i 

    real(rk) ,intent(in)           :: Te !! Electron temperature in eV
    real(rk) ,intent(in)           :: ne !! Electron density in m^{-3} 
    real(rk) ,intent(in)           :: Z !! Ion charge
    logical  ,intent(in), optional :: removeDisc !! Remove the discontinuity at 10eV by moving the branch switch to Te = Z **2 * e**2 eV. Defaults to .false.
    real(rk)             :: res

    real(rk) :: cutoff
    
    cutoff = 10*Z**2
    if (present(removeDisc)) then 
        if (removeDisc) cutoff = exp(real(2,kind=rk))*Z**2
    end if

    if (Te < cutoff) then 

        res = real(23.0d0 - log(sqrt(ne * 1.00D-6) * Z * Te ** (-3.0d0/2.0d0)),kind=rk)

    else

        res = real(24.0d0 - log(sqrt(ne * 1.00D-6) / Te ),kind=rk)

    end if

end function logLei
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function logLii(Z1,Z2,mRatio,n1,n2,T1,T2) result(res)
    !! Calculate Coulomb logarithm for ion-ion collisions (NRL Formulary 2013 page 34 equation c).

    real(rk) ,intent(in) :: Z1 !! Ion charge for first species
    real(rk) ,intent(in) :: Z2 !! Ion charge for second species
    real(rk) ,intent(in) :: mRatio !! Mass ratio m2/m1
    real(rk) ,intent(in) :: n1 !! First ion species density in m^{-3} 
    real(rk) ,intent(in) :: n2 !! Second ion species density in m^{-3} 
    real(rk) ,intent(in) :: T1 !! First ion temperature in eV
    real(rk) ,intent(in) :: T2 !! Second ion temperature in eV
    real(rk)             :: res

    res = real(23.0d0 - log(Z1 * Z2 * (1.0d0 + mRatio)*sqrt(n1 * 1.00D-6 * Z1**2/T1 + n2 * 1.00D-6 * Z2**2/T2)/(T2 + mRatio * T1)),&
               kind=rk)

end function logLii
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function shPotDrop(massRatio,tempRatio) result(res)
    !! Calculate sheath potential drop for 2 component plasma in units of kTe/e

    real(rk) ,intent(in) :: massRatio !! Mass ratio me/mi
    real(rk) ,intent(in) :: tempRatio !! Ion to electron temperature ratio
    real(rk)             :: res

    res = real(0.5d0*log(2*pi*massRatio*(1.0d0+tempRatio)),kind=rk)

end function shPotDrop
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function collTimeei(Te,ne,Z,logL) result(res)
    !! Calculate electron-ion collision time in seconds. This is not the Braginskii time, which can be calculated as this value 
    !! multiplied by 3 sqrt(pi)/4 .

    real(rk) ,intent(in) :: Te !! Electron temperature in eV
    real(rk) ,intent(in) :: ne !! Electron density in  m^{-3} 
    real(rk) ,intent(in) :: Z !! Ion charge
    real(rk) ,intent(in) :: logL !! Coulomb logarithm to be used
    real(rk)             :: res

    res = real(4*pi * epsilon0 ** 2 * sqrt(elMass/elCharge) * (2.0d0*Te) ** (3.0d0/2.0d0) /(Z * logL * ne * elCharge ** 2),kind=rk)

end function collTimeei
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function elPlasmaOscFreq(ne) result(res)
    !! Calculate cold electron plasma oscillation frequency in s^{-1} 

    real(rk) ,intent(in) :: ne !! Electron density in  m^{-3} 
    real(rk)             :: res

    res = real(sqrt(ne * elCharge**2/(elMass * epsilon0)),kind=rk)

end function elPlasmaOscFreq
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function elVthermal(Te) result(res)
    !! Mean thermal electron speed 

    real(rk) ,intent(in) :: Te !! Electron temperature in eV
    real(rk)             :: res

    res = real(sqrt(2*elCharge*Te/elMass),kind=rk)

end function elVthermal
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function vSonic(Te,Ti,mi,ePolytropic,iPolytropic) result(res)
    !! Sonic speed caclulated as sqrt((ePoly*Te + iPoly*Ti)/mi)

    real(rk) ,intent(in) :: Te !! Electron temperature in eV
    real(rk) ,intent(in) :: Ti !! Ion temperature in eV
    real(rk) ,intent(in) :: mi !! Ion mass 
    real(rk) ,intent(in) :: ePolytropic !! Electron polytropic coefficient
    real(rk) ,intent(in) :: iPolytropic !! Ion polytropic coefficient
    real(rk)             :: res

    res = real(sqrt(elCharge * (ePolytropic*Te + iPolytropic*Ti) / mi),kind=rk)

end function vSonic
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function normMaxwellian(n,T,vPoints) result(res)
    !! Returns a normalized stationary Maxwellian evaluated at vPoints. Assumes that velocity and temperature are normalized in 
    !! such way that  m v_0^2/2 = kT_0 , and the returned value is normalized to  1/v_0^(3/2) .

    real(rk)               ,intent(in)  :: n !! Density in arbitrary units
    real(rk)               ,intent(in)  :: T !! Temperature in units compatible with velocity grid
    real(rk) ,dimension(:) ,intent(in)  :: vPoints !! Velocity grid

    real(rk) ,allocatable ,dimension(:) :: res

    res = (pi*T) ** (-3.0d0/2.0d0) * n * exp(-vPoints**2 / T) 

end function normMaxwellian
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function simpleMoment(m,f,vPoints,vWidths,g) result(res)
    !! Returns the m-th moment of passed distribution function harmonic -  4*pi*int(v^(2+m)*f*dv) . Optionally takes moment of f*g
    !! where g conforms to f. The integration is a simple Riemann sum.

    integer(ik)                      ,intent(in)  :: m !! Moment degree
    real(rk)           ,dimension(:) ,intent(in)  :: f !! Distribution function to take the moment of
    real(rk)           ,dimension(:) ,intent(in)  :: vPoints !! Coordinates of cell centres in velocity space
    real(rk)           ,dimension(:) ,intent(in)  :: vWidths !! Widths of velocity space cells
    real(rk) ,optional ,dimension(:) ,intent(in)  :: g !! Optional velocity space array to include in moment integral

    real(rk)                                      :: res

    if (present(g)) then 
        res = 4*pi*sum(f*g*vWidths*vPoints**(2+m))
    else
        res = 4*pi*sum(f*vWidths*vPoints**(2+m))
    end if

end function simpleMoment
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule physics_functions_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
