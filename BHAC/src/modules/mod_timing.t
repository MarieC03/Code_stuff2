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
module mod_timing
! module file for counters of wallclock time.
! 14.05.2012 Oliver Porth
! to use the timers, just 
! use mod_timing
!============================================================================
implicit none
save
double precision       :: time_in, timeio0, timeio_tot=0.0d0
double precision       :: timegr0, timegr_tot=0.0d0, timeloop, timeloop0, timefirstbc
double precision       :: tpartc=0.0d0, tpartc_io=0.0d0, tpartc_int=0.0d0, tpartc_com=0.0d0, tpartc_grid=0.0d0
double precision       :: tpartc0, tpartc_int_0, tpartc_com0, tpartc_io_0, tpartc_grid_0

integer                :: itTimeLast
double precision       :: timeLast
end module mod_timing
!=============================================================================
