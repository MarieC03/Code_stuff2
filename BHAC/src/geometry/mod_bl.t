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
  ! This wraps a module around mod_coord_bl.t to be able to use
  ! the Boyer-Lindquist metric elements for initialization of a problem
  ! to be run in different coordinates.
  !
  ! Usage:
  ! use mod_bl, get_alpha_BL => get_alpha, get_beta_BL => get_beta, &
  ! get_g_component_BL => get_g_component, u4BLtoCoord_BL => u4BLtoCoord
  ! 
  ! Oliver Porth
  ! 2016-01-19
  !=============================================================================

!=============================================================================
module mod_bl

  use mod_metric_aux
  implicit none
  
  INCLUDE:geometry/mod_coord_bl.t


  !=============================================================================
  subroutine get_g4_BL(x^D,g)
    ! Return the four-metric at point x^D
    double precision, intent(in)                           :: x^D
    double precision, dimension(0:^NC,0:^NC),intent(out)   :: g
    ! .. local ..
    integer                                                :: i, j
    double precision, dimension(1:^NC)                     :: betaU, betaD
    double precision                                       :: beta2
    !-----------------------------------------------------------------------------


    do i=1,^NC
       do j=1,^NC
          call get_g_component(i,j,x^D,g(i,j))
       end do
    end do
    call get_alpha(x^D,g(0,0))
    do i=1,^NC
       call get_beta(i,x^D,betaU(i))
    end do
    ! lower the beta:
    betaD(:) = 0.0d0
    do i=1,^NC
       do j=1,^NC
          betaD(i) = betaD(i) + g(i,j)*betaU(j)
       end do
    end do

    
    beta2 = {^C& betaU(^C)*betaD(^C) |+}
    
    g(0,0) = -g(0,0)**2 + beta2

    do i=1,^NC
       g(i,0) = betaD(i)
       g(0,i) = betaD(i)
    end do


  end subroutine get_g4_BL
  !=============================================================================

  
end module mod_bl
