!================================================================================
!
!    BHAC (The Black Hole Accretion Code) solves the equations of
!    general relativistic magnetohydrodynamics and other hyperbolic systems
!    in curved spacetimes.
!
!    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
!    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
!
!    This file is part of BHAC.
!
!    BHAC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BHAC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

!=============================================================================
module constants
! module file for physical constants.
! 01.09.2012 Oliver Porth
! Sets constants in cgs units.
! For usage, the user has to provide the three scaling parameters
! UNIT_LENGTH, UNIT_DENSITY, UNIT_VELOCITY
! such that l[cgs]=UNIT_LENGTH*l[code] and so on.  
! best to set those in initglobaldata_usr.
! to include in a subroutine, just write
! use constants
! as with any other module.  
!============================================================================
implicit none
save

double precision, PARAMETER:: CONST_c   = 2.99792458D10    ! cm s^-1           ; Speed of light
double precision, PARAMETER:: CONST_me  = 9.1093897D-28    ! g                 ; Electron mass
double precision, PARAMETER:: CONST_mp  = 1.672621777D-24  ! g                 ; Proton mass
double precision, PARAMETER:: CONST_e   = 4.8032068D-10    ! g^1/2 cm^3/2 s^-1 ; Electron charge
double precision, PARAMETER:: CONST_MSun= 1.98892D33       ! g                 ; Solar mass
double precision, PARAMETER:: CONST_kB  = 1.3806488d-16    ! erg K^-1          ; Boltzmann constant
double precision, PARAMETER:: CONST_h   = 6.6260755d-27    ! erg s             ; Planck constant
! Conversion factors:
double precision, PARAMETER:: CONST_eV  = 1.6021772d-12    ! erg/eV            ; Electron volt
double precision, PARAMETER:: CONST_Tera= 1.d12            ! -                 ; Tera
double precision, PARAMETER:: CONST_Peta= 1.d15            ! -                 ; Peta
double precision, PARAMETER:: CONST_years= 3600*24*365     ! s year^-1         ; seconds in a year

end module constants
!=============================================================================
