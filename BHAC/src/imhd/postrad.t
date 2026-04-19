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
!
! Assemble quantities used for radiation postprocessing
!
!=============================================================================
subroutine convert_postrad(ixI^L,ixO^L,igrid,nwpost,w,x,wpost)

!  use mod_metric, only: lower4
  use mod_physaux
  include 'amrvacdef.f'
  
  integer, intent(in)              :: ixI^L, ixO^L, igrid, nwpost
  double precision, intent(in)     :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
  double precision, intent(out)    :: wpost(ixI^S,1:nwpost)
  ! .. local ..
  double precision                 :: wprim(ixI^S,1:nw)
  integer, parameter               :: ir_=1, itheta_=2, iphi_=3, &
       irho_=iphi_+1, &
       ivr_=irho_+1, &
       ivtheta_=ivr_+1, &
       ivphi_ = ivtheta_+1, &
       ip_=ivphi_+1, &
       iBr_=ip_+1, &
       iBtheta_=iBr_+1, &
       iBphi_ = iBtheta_+1
  !-----------------------------------------------------------------------------
call mpistop('dont come to postrad')
  if (ibphi_ .gt. nwpost) &
       call mpistop('convert_postrad: nwpost array too small?!')
  
  wpost = zero
  
  ! First convert to primitive variables:
  wprim(ixO^S,1:nw) = w(ixO^S,1:nw)
  call primitive(ixI^L,ixO^L,wprim,x)


  ! Coordinates:
  wpost(ixO^S,ir_)     = x(ixO^S,r_)
  {^IFZIN
  wpost(ixO^S,itheta_) = x(ixO^S,z_)
  }{^IFZOUT
  wpost(ixO^S,itheta_) = dpi/2.0d0
  }{^IFPHIIN
  wpost(ixO^S,iphi_)   = x(ixO^S,phi_)
  }{^IFPHIOUT
  wpost(ixO^S,iphi_)   = zero
  }

  ! density and pressure:
  wpost(ixO^S,irho_) = wprim(ixO^S,rho_)
  wpost(ixO^S,ip_) = wprim(ixO^S,pp_)
  
  ! Fluid three-velocity:
  wpost(ixO^S,ivr_)     = wprim(ixO^S,u1_)/wprim(ixO^S,lfac_)
  {^IFZ
  wpost(ixO^S,ivtheta_) = wprim(ixO^S,u^Z_)/wprim(ixO^S,lfac_)
  }{^IFPHI
  wpost(ixO^S,ivphi_)   = wprim(ixO^S,u^PHI_)/wprim(ixO^S,lfac_)
  }

  ! Magnetic three-field in the Eulerian frame:
  wpost(ixO^S,iBr_)     = wprim(ixO^S,b1_)
  {^IFZ
  wpost(ixO^S,iBtheta_) = wprim(ixO^S,b^Z_)
  }{^IFPHI
  wpost(ixO^S,iBphi_)   = wprim(ixO^S,b^PHI_)
  }

end subroutine convert_postrad
!=============================================================================
