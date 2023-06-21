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
submodule(janev_fits) janev_fits_procedures
    !! author: Stefan Mijin 
    !! 
    !! Contains implementation of functions for calculating rates and cross-sections based on the Janev EIRENE report: http://www.eirene.de/report_4105.pdf

    implicit none
    
contains
!-----------------------------------------------------------------------------------------------------------------------------------   
pure function johnsonY(n,m) result(ynm)
    !! y_nm factor showing up in Johnson's formula for hydrogen excitation from state n to state m
    integer(ik) ,intent(in) :: n ,m

    real(rk) :: ynm 

    ynm = real(1,kind=rk) - real(n,kind=rk) ** 2 / real(m,kind=rk) ** 2

end function johnsonY
!-----------------------------------------------------------------------------------------------------------------------------------   
pure function gFactor(n,ynm) result(g)
    !! g factor showing up in oscillator strength in Johnson's formula for excitation
    integer(ik) ,intent(in) :: n !! Starting state 

    real(rk)    ,intent(in) :: ynm !! Johnson y factor for transition from state n to state m

    real(rk) :: g 

    real(rk) :: g0 ,g1, g2

    if (n == 1) then 

        g0 = real(1.1330d00,kind=rk)
        g1 = real(-0.4059d00,kind=rk)
        g2 = real(0.0714d00,kind=rk)

    else if (n == 2) then 

        g0 = real(1.0785d00,kind=rk)
        g1 = real(-0.2319d00,kind=rk)
        g2 = real(0.0295d00,kind=rk)

    else

        g0 = real(0.9935d00,kind=rk) + real(0.2328d00,kind=rk) / real(n,kind=rk) - real(0.1296d00,kind=rk)  / real(n ** 2,kind=rk)
        g1 = -(real(0.6282d00,kind=rk) - real(0.5598d00,kind=rk) / real(n,kind=rk) &
        + real(0.5299d00,kind=rk) / real(n ** 2,kind=rk)) / real(n,kind=rk)
        g2 = (real(0.3887d00,kind=rk) -real(1.181d00,kind=rk) / real(n,kind=rk) &
        +real(1.1470d00,kind=rk) / real(n ** 2,kind=rk))/real(n ** 2,kind=rk)


    end if

    g = g0 + g1/ynm + g2/ynm

end function gFactor
!-----------------------------------------------------------------------------------------------------------------------------------   
pure function fOsc(n,m) result(f)
    !! Oscillator strength for hydrogen transition from state n to state m
    integer(ik) ,intent(in) :: n ,m

    real(rk) :: f 

    f = real(32.00d00,kind=rk) * real(n,kind=rk) * gFactor(n,johnsonY(n,m)) &
    / (johnsonY(n,m) ** 3 * real(m,kind=rk) ** 3 * real(3.00d00 * sqrt(3.00d00),kind=rk) * pi)

end function fOsc
!-----------------------------------------------------------------------------------------------------------------------------------
pure function csH1To2Ex(eGrid) result(crossSection)
    !! Integral cross section for ground to first excited state excitation of hydrogen - Janev equation (4) in section  2.1.1
    !! NOTE: The returned value is in 10^-16 cm^2!

    real(rk)  ,dimension(:) ,intent(in) :: eGrid !! Energy grid on which to calculate the cross section - energy in eV!
    real(rk) ,allocatable ,dimension(:) :: crossSection

    !Transition energy and Janev fit coefficients
    real(rk) ,parameter :: deltaE = real(10.20d00,kind=rk), & 
                           sigma0 = real(5.984d00,kind=rk), & 
                           a = real(0.228d00,kind=rk), &
                           b = real(0.1865d00,kind=rk), &
                           c = real(0.5025d00,kind=rk), &
                           A0 = real(4.4979d00,kind=rk), &
                           Av(1:5) = real([1.4182d00, -20.877d00, 49.735d00, -46.249d00, 17.442d00 ],kind=rk)

    real(rk) ,allocatable ,dimension(:) :: x ! Energy grid scaled to deltaE 

    integer(ik) :: i , j

    x = eGrid/deltaE 

    allocate(crossSection(size(eGrid)))

    crossSection = 0

    do i = 1, size(eGrid) 

        if (eGrid(i) > deltaE) then 

            if (eGrid(i) <= 11.56d00) then

                crossSection(i) = a + b * (eGrid(i) - deltaE)

            else if (eGrid(i) <= 12.23d00) then 

                crossSection(i) = c

            end if

                crossSection(i) = sigma0 * A0 * log(x(i)) / eGrid(i)

                do j = 1, 5
                    crossSection(i) = crossSection(i) + sigma0 * Av(j) / (deltaE * x(i)**j)
                end do
        end if 

    end do

end function csH1To2Ex
!-----------------------------------------------------------------------------------------------------------------------------------
pure function csH1To345Ex(eGrid,m) result(crossSection)
    !! Integral cross section for ground to m=3,4,5 excitation of hydrogen - Janev equation (5) in section  2.1.1
    !! NOTE: The returned value is in 10^-16 cm^2!

    real(rk)  ,dimension(:) ,intent(in) :: eGrid !! Energy grid on which to calculate the cross section - energy in eV!
    integer(ik)             ,intent(in) :: m !! Excited state (must be 3,4, or 5)
    real(rk) ,allocatable ,dimension(:) :: crossSection


    real(rk) ,allocatable ,dimension(:) :: x ! Energy grid scaled to deltaE 

    !Janev fit parameters
    real(rk) ,parameter :: sigma0 = real(5.984d00,kind=rk)
    real(rk)            :: deltaE ,alpha ,a0 ,A(4)

    integer(ik) :: i , j

    allocate(crossSection(size(eGrid)))

    crossSection = 0

    select case (m)
        case (3)

            deltaE = real(12.09d00,kind=rk)
            alpha = real(0.38277d00,kind=rk)
            A0 = real(0.75448d00,kind=rk)
            A = real([0.42956d00 ,-0.58288d00 ,1.0693d00 , 0.00d00],kind=rk)

        case (4)

            deltaE = real(12.75d00,kind=rk)
            alpha = real(0.41844d00,kind=rk)
            A0 = real(0.24300d00,kind=rk)
            A = real([0.24846d00 ,0.19701d00 ,0.00d00 , 0.00d00],kind=rk)

        case (5)

            deltaE = real(13.06d00,kind=rk)
            alpha = real(0.45929d00,kind=rk)
            A0 = real(0.11508d00,kind=rk)
            A = real([0.13092d00 ,0.23581d00 ,0.00d00 , 0.00d00],kind=rk)

        case default 
            error stop "Janev fit csH1To345Ex called for invalid m"
    end select

    x = eGrid/deltaE 

    do i = 1, size(eGrid) 

        if (eGrid(i) > deltaE) then 

                crossSection(i) = sigma0 * ((real(1,kind=rk) - x(i)**(-1)) ** alpha) * A0 * log(x(i)) / eGrid(i)

                do j = 1, 4
                    crossSection(i) = crossSection(i) &
                    + sigma0 * ((real(1,kind=rk) - x(i)**(-1)) ** alpha) * A(j) / (deltaE * x(i)**j)
                end do
        end if 

    end do

end function csH1To345Ex
!-----------------------------------------------------------------------------------------------------------------------------------
pure function csH1To6PlusEx(eGrid,m) result(crossSection)
    !! Integral cross section for ground to m>5 excitation of hydrogen - Janev equation (6) in section  2.1.1
    !! NOTE: The returned value is in 10^-16 cm^2!

    real(rk)  ,dimension(:) ,intent(in) :: eGrid !! Energy grid on which to calculate the cross section - energy in eV!
    integer(ik)             ,intent(in) :: m !! Excited state (should be >5)
    real(rk) ,allocatable ,dimension(:) :: crossSection


    real(rk) ,allocatable ,dimension(:) :: x ! Energy grid scaled to deltaE 

    !Janev fit parameters
    real(rk)            :: deltaE ,A,r,B ,ym

    integer(ik) :: i 

    allocate(crossSection(size(eGrid)))

    crossSection = 0

    deltaE = real(13.6d00,kind=rk) * (real(1,kind=rk) - real(m,kind=rk) ** (-2))
    x = eGrid/deltaE 

    ym = johnsonY(1,m)
    A = 2 * fOsc(1, m) / ym
    r = real(0.45d00,kind=rk)
    B = 4 * (real(1,kind=rk) + real(4,kind=rk) /(3 * ym) &
        - real(0.603d00,kind=rk)/ ym ** 2) / (real(m ** 3,kind=rk) * ym ** 2)
    do i = 1, size(eGrid) 

        if (eGrid(i) > deltaE) then 

            crossSection(i) = real(1.76d00,kind=rk) * (real(1,kind=rk) - exp(-r * ym * x(i))) &
            * (A * (log(x(i)) + real(1,kind=rk) / (2 * x(i))) + &
            (B - A * log(real(2,kind=rk)/ym)) * (real(1,kind=rk) - real(1,kind=rk)/x(i))) / (ym * x(i))
        end if 

    end do

end function csH1To6PlusEx
!-----------------------------------------------------------------------------------------------------------------------------------
pure function csH2To3Ex(eGrid) result(crossSection)
    !! Integral cross section for  n=2 -> n=3 transition of hydrogen - Janev equation (5) in section  2.1.1 with fitting parameters
    !! from page 10
    !! NOTE: The returned value is in 10^-16 cm^2!

    real(rk)  ,dimension(:) ,intent(in) :: eGrid !! Energy grid on which to calculate the cross section - energy in eV!
    real(rk) ,allocatable ,dimension(:) :: crossSection

    !Transition energy and Janev fit coefficients
    real(rk) ,parameter :: sigma0 = real(5.984d00,kind=rk), &
                           alpha = real(1.3196d00,kind=rk), &
                           A0 = real(38.906d00,kind=rk), &
                           deltaE = real(13.60d0,kind=rk) * (real(1,kind=rk)/real(4,kind=rk) - real(1,kind=rk)/real(9,kind=rk)), &
                           A(1:4) = real([ 5.2373d00, 119.25d00, -595.39d00, 816.71d00 ],kind=rk)

    real(rk) ,allocatable ,dimension(:) :: x ! Energy grid scaled to deltaE 

    integer(ik) :: i  ,j

    x = eGrid/deltaE 

    allocate(crossSection(size(eGrid)))

    crossSection = 0

    do i = 1, size(eGrid) 

        if (eGrid(i) > deltaE) then 

            crossSection(i) = sigma0 * ((real(1,kind=rk) - x(i)**(-1)) ** alpha) * A0 * log(x(i)) / eGrid(i)

            do j = 1, 4
                crossSection(i) = crossSection(i) &
                + sigma0 * ((real(1,kind=rk) - x(i)**(-1)) ** alpha) * A(j) / (deltaE * x(i)**j)
            end do

        end if 

    end do

end function csH2To3Ex
!-----------------------------------------------------------------------------------------------------------------------------------
pure function csH2PlusTo4PlusEx(eGrid,n,m) result(crossSection)
    !! Integral cross section for n>1 to m>3 excitation of hydrogen - Janev equation (9) in section  2.1.1
    !! NOTE: The returned value is in 10^-16 cm^2!

    real(rk)  ,dimension(:) ,intent(in) :: eGrid !! Energy grid on which to calculate the cross section - energy in eV!
    integer(ik)             ,intent(in) :: n !! Initial state (should be >1)
    integer(ik)             ,intent(in) :: m !! Final state (should be >3)
    real(rk) ,allocatable ,dimension(:) :: crossSection

    real(rk) ,allocatable ,dimension(:) :: x ! Energy grid scaled to deltaE 

    !Janev fit parameters
    real(rk)            :: deltaE ,A,r,B ,ynm ,bb

    integer(ik) :: i 

    allocate(crossSection(size(eGrid)))

    crossSection = 0

    deltaE = real(13.6d00,kind=rk) * (real(n,kind=rk)**(-2) - real(m,kind=rk)** (-2))

    x = eGrid/deltaE 

    ynm = johnsonY(n,m)
    A = 2 * real(n ** 2,kind=rk) * fOsc(n,m) / ynm
    r = real(1.94d00,kind=rk)/real(n,kind=rk) ** real(1.57d00,kind=rk)
    bb = (real(4,kind=rk) - real(18.63d00,kind=rk)/real(n,kind=rk) + real(36.24d00)/real(n ** 2,kind=rk)&
     - real(28.09d00,kind=rk)/real(n ** 3,kind=rk)) / real(n,kind=rk)
    B = 4 * real(n ** 4,kind=rk)&
     * (real(1,kind=rk) + real(4,kind=rk)/(3 * ynm) + bb/ynm ** 2) &
     / (real(m ** 3,kind=rk) * ynm ** 2)

    do i = 1, size(eGrid) 

        if (eGrid(i) > deltaE) then 

            crossSection(i) = real(1.76d00,kind=rk) * real(n ** 2,kind=rk) * (real(1,kind=rk) - exp(- r * ynm * x(i))) * &
               (A * (log(x(i)) + real(1,kind=rk)/(2 * x(i))) &
               + (B - A * log(2 * real(n ** 2,kind=rk) /ynm))&
                * (real(1,kind=rk) - real(1,kind=rk)/x(i))) / (ynm * x(i))
        end if 

    end do

end function csH2PlusTo4PlusEx
!-----------------------------------------------------------------------------------------------------------------------------------
pure function csH123Ion(eGrid,n) result(crossSection)
    !! Integral cross section for ionization of hydrogen states n=1,2,3 - Janev equation (14) in section  2.1.2
    !! NOTE: The returned value is in 10^-16 cm^2!

    real(rk)  ,dimension(:) ,intent(in) :: eGrid !! Energy grid on which to calculate the cross section - energy in eV!
    integer(ik)             ,intent(in) :: n !! Initial state (must be 1,2,or 3)
    real(rk) ,allocatable ,dimension(:) :: crossSection

    !Janev fit parameters
    real(rk)            :: ionPot ,A0 ,A(5)

    integer(ik) :: i , j

    allocate(crossSection(size(eGrid)))

    crossSection = 0

    ionPot = real(13.6d00,kind=rk)/real(n**2,kind=rk)

    select case (n)
        case (1)

            A0 = real(0.18450d00,kind=rk)
            A = real([ -0.032226d00, -0.034539d00, 1.4003d00, -2.8115d00, 2.2986d00 ],kind=rk)

        case (2)

            A0 = real(0.14784d00,kind=rk)
            A = real([ 0.0080871d00, -0.062270d00, 1.9414d00, - 2.1980d00, 0.95894d00 ],kind=rk)

        case (3)

            A0 = real(0.058463d00,kind=rk)
            A = real([ -0.051272d00, 0.85310d00, -0.57014d00, 0.76684d00, 0.00d00 ],kind=rk)

        case default 
            error stop "Janev fit csH123Ion called for invalid n"
    end select

    do i = 1, size(eGrid) 

        if (eGrid(i) > ionPot) then 

                crossSection(i) = 1.0d3 * A0 * log(eGrid(i)/ionPot)/(ionPot*eGrid(i))

                do j = 1, 5
                    crossSection(i) = crossSection(i) &
                    + 1.0d3 * A(j) * (real(1,kind=rk) - ionPot/eGrid(i)) ** j / (ionPot*eGrid(i))
                end do
        end if 

    end do

end function csH123Ion
!-----------------------------------------------------------------------------------------------------------------------------------
pure function csH4PlusIon(eGrid,n) result(crossSection)
    !! Integral cross section for ionization of hydrogen states n>3 - Janev equation (15) in section  2.1.2
    !! NOTE: The returned value is in 10^-16 cm^2!

    real(rk)  ,dimension(:) ,intent(in) :: eGrid !! Energy grid on which to calculate the cross section - energy in eV!
    integer(ik)             ,intent(in) :: n !! Initial state (should be n>3)
    real(rk) ,allocatable ,dimension(:) :: crossSection

    !Janev fit parameters
    real(rk) ,parameter :: sigma0 = real(5.984d00,kind=rk)
    real(rk)            :: ionPot ,g0,g1,g2,A,r,bb,B

    real(rk) ,allocatable ,dimension(:) :: x ! Energy grid scaled to ionPot

    integer(ik) :: i 

    allocate(crossSection(size(eGrid)))

    crossSection = 0

    ionPot = real(13.6d00,kind=rk)/real(n**2,kind=rk)

    g0 = real(0.9935d00,kind=rk) + real(0.2328d00,kind=rk) / real(n,kind=rk) - real(0.1296d00,kind=rk)  / real(n ** 2,kind=rk)

    g1 = -(real(0.6282d00,kind=rk) - real(0.5598d00,kind=rk) / real(n,kind=rk) &
    + real(0.5299d00,kind=rk) / real(n ** 2,kind=rk)) / real(n,kind=rk)

    g2 = (real(0.3887d00,kind=rk) -real(1.181d00,kind=rk) / real(n,kind=rk) &
    +real(1.1470d00,kind=rk) / real(n ** 2,kind=rk))/real(n ** 2,kind=rk)

    A = 32 * n * (g0/real(3,kind=rk) + g1/real(4,kind=rk) + g2/real(5,kind=rk)) / (3 * PI * sqrt(real(3,kind=rk)))
    r = real(1.94d00,kind=rk)/ real(n,kind=rk) ** real(1.57d00,kind=rk)
    bb = (real(4,kind=rk) - real(18.63d00,kind=rk)/real(n,kind=rk) + real(36.24d00,kind=rk)/real(n ** 2 ,kind=rk)&
    - real(28.09d00,kind=rk)/real(n ** 3,kind=rk)) &
    / real(n,kind=rk)
    B = 2 *real(n ** 2,kind=rk) * (real(5,kind=rk) + bb)/real(3,kind=rk)

    x = eGrid/ionPot

    do i = 1, size(eGrid) 

        if (eGrid(i) > ionPot) then 

            crossSection(i) = real(1.76d00,kind=rk) * real(n ** 2,kind=rk) * (real(1,kind=rk) - exp(- r * x(i))) * &
            (A * log(x(i)) + (B - A * log(real(2*n ** 2,kind=rk))) * (real(1,kind=rk) - real(1,kind=rk)/x(i)) ** 2) / x(i)

        end if 

    end do

end function csH4PlusIon
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function radRecombRateHydrogen(temp,n) result(rate)
    !! Radiative recombination rate given by formulas (21) and (22) from Janev Section 2.1.4. The recombination for n=2 is obtained
    !! by a simple sum over the two l states. Result in 10^-14 cm^3/s

    real(rk) ,dimension(:) ,intent(in)   :: temp !! Temperature array in eV
    integer(ik)            ,intent(in)   :: n !! Final atomic state after recombination
    real(rk) ,allocatable ,dimension(:)  :: rate
    
    real(rk) :: ionPot
    real(rk) ,allocatable ,dimension(:) :: beta

    ionPot = real(13.6d00,kind=rk)/real(n**2,kind=rk)
    beta = ionPot/temp

    select case (n)
    case (1)
        rate = real(3.92d00,kind=rk) * beta ** (real(1.5d0,kind=rk))/(beta + real(0.35d0,kind=rk)) 
    case (2)
        rate = real(2.47d00,kind=rk) * real(0.5d0,kind=rk) * beta ** (real(1.5d0,kind=rk))/(beta + real(0.12d0,kind=rk)) + &
               real(6.22d00,kind=rk) * real(0.5d0,kind=rk) * beta ** (real(1.5d0,kind=rk))/(beta + real(0.61d0,kind=rk))
    case default
        rate = real(5.201d00,kind=rk) * beta ** (real(1.5d0,kind=rk)) * expInt1(beta) * exp(beta)
    end select

end function radRecombRateHydrogen
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function excitationCrossSectionHydrogen(eGrid,n,m) result(cs)
    !! Calculates hydrogen electron impact excitation cross-section from state n to state m  on given energy grid.

    real(rk) ,dimension(:) ,intent(in)   :: eGrid !! Electron energy grid in eV
    integer(ik)            ,intent(in)   :: n !! Initial atomic state
    integer(ik)            ,intent(in)   :: m !! Final atomic state
    real(rk) ,allocatable ,dimension(:)  :: cs

    if (n == 1) then 

        if (m == 2) cs = csH1To2Ex(eGrid)
        if (any([3,4,5] == m)) cs = csH1To345Ex(eGrid,m)
        if (m > 5) cs = csH1To6PlusEx(eGrid,m)

    else if ((n == 2) .and. (m == 3)) then 

        cs = csH2To3Ex(eGrid)

    else 

        cs = csH2PlusTo4PlusEx(eGrid,n,m)

    end if

end function excitationCrossSectionHydrogen
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function ionizationCrossSectionHydrogen(eGrid,n) result(cs)
    !! Calculates hydrogen electron impact ionization cross-section from state n on given energy grid.

    real(rk) ,dimension(:) ,intent(in)   :: eGrid !! Electron energy grid in eV
    integer(ik)            ,intent(in)   :: n !! Initial atomic state
    real(rk) ,allocatable ,dimension(:)  :: cs

    if (any([1,2,3] == n)) then 
        cs = csH123Ion(eGrid,n)
    else
        cs = csH4PlusIon(eGrid,n)
    end if

end function ionizationCrossSectionHydrogen
!-----------------------------------------------------------------------------------------------------------------------------------   
end submodule janev_fits_procedures
!-----------------------------------------------------------------------------------------------------------------------------------