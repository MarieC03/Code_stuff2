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

{#IFDEF PARTICLES
!=============================================================================
subroutine integrate_particles()

!-----------------------------------------------------------------------------
{#IFDEF PARTICLES_LORENTZ
{^IFCOORDCART
call integrate_particles_lorentz
}
}
{#IFDEF PARTICLES_ADVECT
call integrate_particles_advect
}
{#IFDEF PARTICLES_GEODESIC
call integrate_particles_geo
}
{#IFDEF PARTICLES_GRLORENTZ
call integrate_particles_grlorentz
}

end subroutine integrate_particles
!=============================================================================
subroutine set_particles_dt()

!-----------------------------------------------------------------------------
{#IFDEF PARTICLES_LORENTZ
call set_particles_dt_lorentz
}
{#IFDEF PARTICLES_ADVECT
call set_particles_dt_advect
}
{#IFDEF PARTICLES_GEODESIC
call set_particles_dt_geo
}
{#IFDEF PARTICLES_GRLORENTZ
call set_particles_dt_grlorentz
}
end subroutine set_particles_dt
!=============================================================================
subroutine set_tmax_particles

  use mod_particles, only: tmax_particles
  include 'amrvacdef.f'
  !-----------------------------------------------------------------------------

  if (time_advance) tmax_particles = (t + dt)

  {#IFDEF PARTICLES_LORENTZ
  ! uses scaled fields:
  tmax_particles = tmax_particles * (UNIT_LENGTH/UNIT_VELOCITY)
  }

  
end subroutine set_tmax_particles
!=============================================================================
subroutine set_tsync_particles

  use mod_particles, only: tsync_particles, tmax_particles, dtsync_particles, t_particles
  include 'amrvacdef.f'
  logical, save               :: first=.true.
  !-----------------------------------------------------------------------------
  
  if (first) then
     tsync_particles = t
     first = .false.
  end if
  
  tsync_particles = min(tsync_particles + dtsync_particles, tmax_particles)

  {#IFDEF PARTICLES_LORENTZ
  ! uses scaled fields:
  call mpistop('set_tsync_particles: fixme for Lorentz pusher!')
  }

   end subroutine set_tsync_particles
!=============================================================================
}




{#IFDEF PARTICLES
{#IFDEF PARTICLES_ADVECT
!=============================================================================
subroutine integrate_particles_advect()
! this solves dx/dt=v for particles

use mod_odeint
use mod_particles
use mod_gridvars, only: get_vec, vp^C_, rhop_, pressp_, bb2p_, igrid_working, interpolate_var
include 'amrvacdef.f'
integer                             :: ipart, iipart
double precision                    :: dt_p
double precision, dimension(1:ndir) :: v, x
integer                             :: igrid
double precision                    :: rho, rho1, rho2, td, tloc
double precision, dimension(ixG^T,1:nw)   :: w
! for odeint:
integer                          :: nok, nbad
double precision                 :: h1
double precision,parameter       :: eps=1.0d-6, hmin=1.0d-8
integer, parameter               :: nvar = ^NC
double precision, dimension(1:bb2p_-ndir)  :: tmp_payload
external derivs_advect
!-----------------------------------------------------------------------------

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);
   
   dt_p = particle(ipart)%self%dt
   igrid = particle(ipart)%igrid
   igrid_working = igrid
   tloc = particle(ipart)%self%t
   x(1:ndir) = particle(ipart)%self%x(1:ndir)
   call set_tmpGlobals(igrid)

   ! **************************************************
   ! Position update
   ! **************************************************

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Simple forward Euler:
   call get_vec(igrid,x,tloc,v,vp1_,vp^NC_)
   particle(ipart)%self%u(1:ndir) = v(1:ndir)

   particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
        + dt_p * v(1:ndir)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Adaptive stepwidth RK4:
!   h1 = dt_p/2.0d0
!   call odeint(x,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_advect,rkqs)
!   particle(ipart)%self%x(1:ndir) = x(1:ndir)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




   ! **************************************************
   ! Velocity update
   ! **************************************************

   call get_vec(igrid,x,tloc+dt_p,v,vp1_,vp^NC_)
   particle(ipart)%self%u(1:ndir) = v(1:ndir)





   ! **************************************************
   ! Payload update
   ! **************************************************
   
   ! To give an example, we set the payload to the interpolated density at 
   ! the particles position.  
   ! In general, it would be better to add the auxilary variable to mod_gridvars 
   ! since then you can use the convenient subroutine get_vec() for this.  
   call get_vec(igrid,x,tloc+dt_p,tmp_payload,rhop_,bb2p_)
   
  
   particle(ipart)%self%payload(1) = tmp_payload(1)
   particle(ipart)%self%payload(2) = tmp_payload(2)
   particle(ipart)%self%payload(3) = tmp_payload(3)


   ! **************************************************
   ! Time update
   ! **************************************************
   particle(ipart)%self%t = particle(ipart)%self%t + dt_p


end do

!=============================================================================
end subroutine integrate_particles_advect
!=============================================================================
subroutine derivs_advect(t_s,x,dxdt)

use mod_gridvars, only: get_vec, vp^C_, igrid_working
include 'amrvacdef.f'

double precision                :: t_s, x(*)
double precision                :: dxdt(*)
! .. local ..
double precision                :: v(1:^NC)
!-----------------------------------------------------------------------------

call get_vec(igrid_working,x(1:^NC),t_s,v,vp1_,vp^NC_)
{^C& dxdt(^C) = v(^C);}

end subroutine derivs_advect
!=============================================================================
subroutine set_particles_dt_advect()

use mod_particles
include 'amrvacdef.f'

integer                         :: ipart, iipart
double precision                :: t_min_mype, tout, dt_particles_mype, dt_cfl
double precision                :: v(1:ndir)
!-----------------------------------------------------------------------------
t_min_mype        = tsync_particles
dt_particles_mype = bigdouble


do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

! make sure we step only one cell at a time:
{v(^C)   = abs(particle(ipart)%self%u(^C))\}

{#IFNDEF D1
   dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)})
}{#IFDEF D1
   dt_cfl = rnode(rpdx1_,particle(ipart)%igrid)/v(1)
}

   particle(ipart)%self%dt = limit_dt(particle(ipart)%self%t,dt_cfl)

   dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)

   t_min_mype = min(t_min_mype,particle(ipart)%self%t)
   
end do !ipart

call save_tmin(t_min_mype,dt_particles_mype)
 
end subroutine set_particles_dt_advect
}
}
!=============================================================================
{#IFDEF PARTICLES
{#IFDEF PARTICLES_LORENTZ
subroutine integrate_particles_lorentz()
! this is the relativistic Boris scheme, a leapfrog integrator
use constants
use mod_particles
use mod_metric
include 'amrvacdef.f'
integer                             :: ipart, iipart
double precision                    :: lfac, q, m, dt_p, cosphi, sinphi, phi1, phi2, r, re
double precision, dimension(1:ndir) :: b, e, emom, uminus, t_geom, s, udash, tmp, uplus, xcart1, xcart2, ucart2, radmom
double precision, dimension(1:ndir) :: uhalf, tau, uprime, ustar, tovergamma
double precision                    :: lfacprime, sscal, sigma
!-----------------------------------------------------------------------------

! Boris only implemented in Minkowski space (coord="cart")
if (coord .ne. "cart") then   
call mpistop('Boris scheme only implemented in Minkowski space')
end if

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

   q  = particle(ipart)%self%q
   m  = particle(ipart)%self%m
   dt_p = particle(ipart)%self%dt
   
   ! Push particle over half time step
   call get_lfac(particle(ipart)%self%u,lfac)
   particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * CONST_C / UNIT_LENGTH

   ! Get E, B at new position
   call get_b(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,b)
   call get_e(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,e)

select case (typeaxial)
! **************************************************
! CARTESIAN COORDINATES
! **************************************************
case ('slab')

! **************************************************
! Momentum update
! **************************************************
{#IFDEF PARTICLES_VAY

!The Vay mover
if (losses) then
call mpistop('radiative losses not implemented in vay, use Boris')
end if
        call cross(particle(ipart)%self%u,b,tmp)
        uhalf = particle(ipart)%self%u + &
             q * dt_p /(2.0d0 * m * CONST_C) * (e &
             + tmp/(lfac))
       
        tau = q * dt_p / (2.d0 * m * CONST_C) * b
        uprime = uhalf + q * dt_p / (2.d0 * m * CONST_C) * e
        call get_lfac(uprime,lfacprime)
        sigma = lfacprime**2 - ({tau(^C)*tau(^C)|+})
       
        ustar = ({uprime(^C)*tau(^C)|+})    
        lfac = sqrt((sigma + sqrt(sigma**2 + 4.d0 * &
             (({tau(^C)*tau(^C)|+}) + ({ustar(^C)*ustar(^C)|+})))) / 2.d0)

        tovergamma = tau / lfac
        sscal = 1.d0 / (1.d0 + ({tovergamma(^C)*tovergamma(^C)|+}))
        call cross(uprime,tovergamma,tmp)
        particle(ipart)%self%u = sscal * (uprime +
({uprime(^C)*tovergamma(^C)|+}) & 
        * tovergamma + tmp)      
}

!The Boris mover
{#IFNDEF PARTICLES_VAY
    emom = q * e * dt_p /(2.0d0 * m * CONST_C)
      if (losses) then
         call get_lfac(particle(ipart)%self%u,lfac)
         re = abs(q)**2 / (m * CONST_C**2)
         call cross(particle(ipart)%self%u,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*particle(ipart)%self%u(^C)|+})/lfac)**2) &
              * particle(ipart)%self%u / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if
      
      uminus = particle(ipart)%self%u + emom + radmom
      call get_lfac(uminus,lfac)
      call get_t(b,lfac,dt_p,q,m,t_geom)
      call get_s(t_geom,s)
      call cross(uminus,t_geom,tmp)
      udash = uminus + tmp
      call cross(udash,s,tmp)
      uplus = uminus + tmp
      if (losses) then
         call cross(uplus,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*uplus(^C)|+})/lfac)**2) &
              * uplus / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if
      
      particle(ipart)%self%u = uplus + emom + radmom
}
! **************************************************
! CYLINDRICAL COORDINATES
! **************************************************
case ('cylindrical')
      
{#IFDEF PARTICLES_VAY
call mpistop('Vay scheme not implemented in cylindrical coordinates')
}
! **************************************************
! Momentum update
! **************************************************
      emom = q * e * dt_p /(2.0d0 * m * CONST_C)

      if (losses) then
         call get_lfac(particle(ipart)%self%u,lfac)
         re = abs(q)**2 / (m * CONST_C**2)
         call cross(particle(ipart)%self%u,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*particle(ipart)%self%u(^C)|+})/lfac)**2) &
              * particle(ipart)%self%u / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if
      
      uminus = particle(ipart)%self%u + emom + radmom
      
      call get_lfac(uminus,lfac)
      call get_t(b,lfac,dt_p,q,m,t_geom)
      call get_s(t_geom,s)
      
      call cross(uminus,t_geom,tmp)
      udash = uminus + tmp
      
      call cross(udash,s,tmp)
      uplus = uminus + tmp

      if (losses) then
         call cross(uplus,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*uplus(^C)|+})/lfac)**2) &
              * uplus / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if

      
      particle(ipart)%self%u = uplus + emom + radmom

! **************************************************
! Position update
! **************************************************
      
      ! Get cartesian coordinates:
      phi1        = particle(ipart)%self%x(pphi_)
      cosphi     = cos(phi1)
      sinphi     = sin(phi1)

      xcart1(1)  = particle(ipart)%self%x(r_) * cosphi
      xcart1(2)  = particle(ipart)%self%x(r_) * sinphi
      xcart1(3)  = particle(ipart)%self%x(zz_)

      ucart2(1)   = cosphi * particle(ipart)%self%u(r_) - sinphi * particle(ipart)%self%u(pphi_)
      ucart2(2)   = cosphi * particle(ipart)%self%u(pphi_) + sinphi * particle(ipart)%self%u(r_)
      ucart2(3)   = particle(ipart)%self%u(zz_)

      ! update position
      xcart2(1:ndir) = xcart1(1:ndir) &
           + dt_p * ucart2(1:ndir)/lfac &
           * CONST_C/UNIT_LENGTH


      ! back to cylindrical coordinates
      phi2     = atan2(xcart2(2),xcart2(1))
      if (phi2 .lt. 0.0d0) phi2 = 2.0d0*dpi + phi2
      r       = sqrt(xcart2(1)**2 + xcart2(2)**2)
      particle(ipart)%self%x(r_)   = r
      particle(ipart)%self%x(pphi_) = phi2
      particle(ipart)%self%x(zz_)   = xcart2(3)

! **************************************************
! Rotate the momentum to the new cooridnates
! **************************************************

      ! rotate velocities
      cosphi     = cos(phi2-phi1)
      sinphi     = sin(phi2-phi1)

      tmp = particle(ipart)%self%u
      particle(ipart)%self%u(r_)   = cosphi * tmp(r_)   + sinphi * tmp(pphi_)
      particle(ipart)%self%u(pphi_) = cosphi * tmp(pphi_) - sinphi * tmp(r_)
      particle(ipart)%self%u(zz_)   = tmp(zz_)

end select

! **************************************************
!> Position update over half timestep at the end
! **************************************************
      ! update position
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * CONST_C / UNIT_LENGTH   

! **************************************************
! Time update
! **************************************************
   particle(ipart)%self%t = particle(ipart)%self%t + dt_p

! **************************************************
!> Payload update
! **************************************************
      if (npayload > 0) then
      !> current gyroradius
      call cross(particle(ipart)%self%u,b,tmp)
      tmp = tmp / sqrt({b(^C)**2|+})
      particle(ipart)%self%payload(1) = sqrt({tmp(^C)**2|+}) /sqrt({b(^C)**2|+}) * m / abs(q) * 8.9875d+20
      end if

      ! e.b payload
      if (npayload>1) then
      !> e.b (set npayload=2 first):
      particle(ipart)%self%payload(2) = ({e(^C)*b(^C)|+})/ (sqrt({b(^C)**2|+})*sqrt({e(^C)**2|+}))
      end if

end do ! ipart
!=============================================================================
contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine get_t(b,lfac,dt,q,m,t)

use constants
implicit none
double precision, dimension(^NC), intent(in)      :: b
double precision, intent(in)                      :: lfac, dt, q, m
double precision, dimension(^NC), intent(out)     :: t
!-----------------------------------------------------------------------------

t = q * b * dt / (2.0d0 * lfac * m * CONST_C)

end subroutine get_t
!=============================================================================
subroutine get_s(t,s)

implicit none
double precision, dimension(^NC), intent(in)   :: t
double precision, dimension(^NC), intent(out)  :: s
!-----------------------------------------------------------------------------

s = 2.0d0 * t / (1.0d0+{t(^C)**2|+})

end subroutine get_s
!=============================================================================
end subroutine integrate_particles_lorentz
!=============================================================================
subroutine set_particles_dt_lorentz()

use constants
use mod_particles
include 'amrvacdef.f'

integer                         :: ipart, iipart
double precision,dimension(^NC) :: b
double precision                :: lfac,absb,dt_particles_mype,dt_cfl,v(ndir)
double precision                :: t_min_mype, tout
double precision, parameter     :: cfl=0.8d0
!-----------------------------------------------------------------------------
t_min_mype        = tsync_particles
dt_particles_mype = bigdouble

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);
   
   call get_b(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,b)
   absb = sqrt({b(^C)**2|+})
   call get_lfac(particle(ipart)%self%u,lfac)

    !**************************************************
    ! CFL timestep
    !**************************************************
    ! make sure we step only one cell at a time:
    {v(^C)   = abs(CONST_C * particle(ipart)%self%u(^C) / lfac)\}

    {^IFPHI 
    ! convert to angular velocity:
    if (typeaxial =='cylindrical') v(phi_) = abs(v(phi_)/particle(ipart)%self%x(r_))
    }

    {#IFNDEF D1
       dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)})
    }
    {#IFDEF D1
       dt_cfl = rnode(rpdx1_,particle(ipart)%igrid)/v(1)
    }

    {^IFPHI
    if (typeaxial =='cylindrical') then 
    ! phi-momentum leads to radial velocity:
    if (phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
         sqrt(rnode(rpdx1_,particle(ipart)%igrid)/particle(ipart)%self%x(r_)) &
         / v(phi_))
    ! limit the delta phi of the orbit (just for aesthetic reasons):
       dt_cfl = min(dt_cfl,0.1d0/v(phi_))
    ! take some care at the axis:
       dt_cfl = min(dt_cfl,(particle(ipart)%self%x(r_)+smalldouble)/v(r_))
    end if
    }

    dt_cfl = dt_cfl * cfl
    !**************************************************
    !**************************************************

    ! bound by gyro-rotation:
    particle(ipart)%self%dt = &
         abs( dtheta * CONST_c/UNIT_LENGTH * particle(ipart)%self%m * lfac &
         / (particle(ipart)%self%q * absb) )

    particle(ipart)%self%dt = min(particle(ipart)%self%dt,dt_cfl)*UNIT_LENGTH

    particle(ipart)%self%dt = limit_dt(particle(ipart)%self%t,particle(ipart)%self%dt)

    dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)
 
    t_min_mype = min(t_min_mype,particle(ipart)%self%t)

end do !ipart

call save_tmin(t_min_mype,dt_particles_mype)

end subroutine set_particles_dt_lorentz
!=============================================================================
subroutine get_e(igrid,x,tloc,e)

! Get the electric field in the grid at postion x.
! For ideal SRMHD, we first interpolate b and u=lfac*v/c
! The electric field then follows from e = b x beta, where beta=u/lfac.  
! This ensures for the resulting e that e<b and e.b=0. Interpolating on u 
! avoids interpolation-errors leading to v>c.  
!
! For (non-ideal) MHD, we directly interpolate the electric field as 
! there is no such constraint.

use mod_particles
use mod_gridvars
include 'amrvacdef.f'

integer,intent(in)                                 :: igrid
double precision,dimension(^NC), intent(in)        :: x
double precision, intent(in)                       :: tloc
double precision,dimension(^NC), intent(out)       :: e
double precision,dimension(^NC)                    :: e1, e2
double precision,dimension(^NC)                    :: u, u1, u2, beta, b
double precision                                   :: lfac
integer                                            :: ic^D
double precision                                   :: td
!-----------------------------------------------------------------------------

{^IFPHYSRMHD

if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e2(^C))\}

td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt

{^C& e(^C) = e1(^C) * (1.0d0 - td) + e2(^C) * td\}

end if

}

{^IFPHYSRRMHD
if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e2(^C))\}

td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt

{^C& e(^C) = e1(^C) * (1.0d0 - td) + e2(^C) * td\}

end if
}

end subroutine get_e
!=============================================================================
subroutine get_b(igrid,x,tloc,b)

use mod_particles
use mod_gridvars
include 'amrvacdef.f'

integer,intent(in)                                 :: igrid
double precision,dimension(^NC), intent(in)        :: x
double precision, intent(in)                       :: tloc
double precision,dimension(^NC), intent(out)       :: b
integer                                            :: ic^D
double precision,dimension(^NC)                    :: b1, b2
double precision                                   :: td
!-----------------------------------------------------------------------------

if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b2(^C))\}

td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt

{^C& b(^C) = b1(^C) * (1.0d0 - td) + b2(^C) * td\}

end if

end subroutine get_b
!=============================================================================
subroutine get_lfac(u,lfac)

use constants

double precision,dimension(^NC), intent(in)        :: u
double precision, intent(out)                      :: lfac
!-----------------------------------------------------------------------------

lfac = sqrt(1.0d0 + ({u(^C)**2|+}))

end subroutine get_lfac
!=============================================================================
subroutine cross(a,b,c)

include 'amrvacdef.f'
double precision, dimension(^NC), intent(in)   :: a,b
double precision, dimension(^NC), intent(out)  :: c
!-----------------------------------------------------------------------------

! ndir needs to be three for this to work!!!
{#IFDEF C3
select case (typeaxial)
case ('slab')
   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)
case ('cylindrical')
   c(r_) = a(pphi_)*b(zz_) - a(zz_)*b(pphi_)
   c(pphi_) = a(zz_)*b(r_) - a(r_)*b(zz_)
   c(zz_) = a(r_)*b(pphi_) - a(pphi_)*b(r_)
case default
   call mpistop('geometry not implemented in cross(a,b,c)')
end select
}{#IFNDEF C3
call mpistop("Needs to be run with three components!")
}

end subroutine cross
!=============================================================================
}
}
!=============================================================================
{#IFDEF PARTICLES
{#IFDEF PARTICLES_GEODESIC
!=============================================================================
subroutine integrate_particles_geo()
! this integrates the geodesic equation for particles

use mod_metric
use mod_particles
use mod_odeint
use mod_gridvars, only: ipart_working
include 'amrvacdef.f'
integer                             :: ipart, iipart
double precision                    :: dt_p,alpha_p
double precision, dimension(1:^NC)  :: u, x
double precision, dimension(1:2*^NC):: k1,k2,k3,k4
double precision                    :: tloc
integer                             :: integrator=3 ! 1=RK4, 2=odeint, 3=IMR
! IMR variables
double precision                    :: IMRtol,er
integer                             :: IMRnmax,nk
double precision, dimension(1:^NC)  :: xk,uk
! odeint variables
integer                             :: nok, nbad
double precision                    :: h1, hmin
!!!!! Precision of time-integration: !!!!!!!!!!!!!
double precision,parameter          :: eps=1.0d-3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, parameter                  :: nvar = 2*^NC
double precision,dimension(1:nvar)  :: y
external derivs_geodesic
!-----------------------------------------------------------------------------

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

  dt_p = particle(ipart)%self%dt
  x(1:^NC) = particle(ipart)%self%x(1:^NC)
  u(1:^NC) = particle(ipart)%self%u(1:^NC)
!  call set_tmpGlobals(igrid)
   

  ! **************************************************
  ! Position and velocity update
  ! **************************************************

  select case(integrator)
  case(1) ! RK4
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Runge-Kutta 4th order:
   
    ! First sweep
    call get_geo_rk4_funcs(particle(nparticles)%self%m,x,u,k1(1:^NC),k1(1+^NC:2*^NC))
    y(1:^NC) = x(1:^NC) + dt_p/2.d0 * k1(1:^NC)
    y(1+^NC:2*^NC) = u(1:^NC) + dt_p/2.d0 * k1(1+^NC:2*^NC)
    ! Second sweep
    call get_geo_rk4_funcs(particle(nparticles)%self%m,y(1:^NC),y(1+^NC:2*^NC),k2(1:^NC),k2(1+^NC:2*^NC))
    y(1:^NC) = x(1:^NC) + dt_p/2.d0 * k2(1:^NC)
    y(1+^NC:2*^NC) = u(1:^NC) + dt_p/2.d0 * k2(1+^NC:2*^NC)
    ! Third sweep
    call get_geo_rk4_funcs(particle(nparticles)%self%m,y(1:^NC),y(1+^NC:2*^NC),k3(1:^NC),k3(1+^NC:2*^NC))
    y(1:^NC) = x(1:^NC) + dt_p * k3(1:^NC)
    y(1+^NC:2*^NC) = u(1:^NC) + dt_p * k3(1+^NC:2*^NC) 
    ! Fourth sweep
    call get_geo_rk4_funcs(particle(nparticles)%self%m,y(1:^NC),y(1+^NC:2*^NC),k4(1:^NC),k4(1+^NC:2*^NC))

    ! Advance
    particle(ipart)%self%x(1:^NC) = x(1:^NC) + dt_p/6.d0 * (k1(1:^NC)+2.d0*k2(1:^NC)+2.d0*k3(1:^NC)+k4(1:^NC))
    particle(ipart)%self%u(1:^NC) = u(1:^NC) + dt_p/6.d0 * (k1(1+^NC:2*^NC)+2.d0*k2(1+^NC:2*^NC)+2.d0*k3(1+^NC:2*^NC)+k4(1+^NC:2*^NC))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case(2) ! Adaptive RK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ODEint from Oliver's routines
    tloc = particle(ipart)%self%t
    ipart_working = particle(ipart)%self%index

    ! initial solution vector:
    y(1:^NC) = x(1:^NC) 
    y(^NC+1:2*^NC) = u(1:^NC)

    h1 = dt_p/2.0d0; hmin=h1/128.0d0
    call odeint(y,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_geodesic,rkqs)

    ! final solution vector:
    particle(ipart)%self%x(1:^NC) = y(1:^NC)
    particle(ipart)%self%u(1:^NC) = y(^NC+1:2*^NC)

   case(3) ! IMR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2nd-order Implicit Midpoint Rule (symplectic)
    ! Uses fixed-point iteration (full Jacobian requires some restructuring)

    IMRtol = 1.d-14 ! Absolute iteration tolerance
    IMRnmax = 100 ! Maximum iteration number

    ! Iteration variables
    xk(1:^NC) = x(1:^NC)
    uk(1:^NC) = u(1:^NC)

    ! Preconditioning with 1st order Euler
    call get_geo_rk4_funcs(particle(nparticles)%self%m,x,u,k1(1:^NC),k1(1+^NC:2*^NC))
    xk(1:^NC) = x(1:^NC) + dt_p/2.d0 * k1(1:^NC)
    uk(1:^NC) = u(1:^NC) + dt_p/2.d0 * k1(1+^NC:2*^NC)

    ! Start nonlinear iteration
    er = 1.d0 ! Initialise error
    do nk = 1,IMRnmax
      ! Compute residuals
      call get_geo_IMR_funcs(particle(nparticles)%self%m,xk,uk,k1(1:^NC),k1(1+^NC:2*^NC)) ! Get the RHS as a function of xk,uk
      k1(1:^NC) = xk(1:^NC) - x(1:^NC) - dt_p/2.d0*k1(1:^NC) ! F = xk - x - dt/2*RHS(xk,uk)
      k1(1+^NC:2*^NC) = uk(1:^NC) - u(1:^NC) - dt_p/2.d0*k1(1+^NC:2*^NC) ! H = uk - u - dt/2*RHS(xk,uk)
      ! Update and error check
      xk(1:^NC) = xk(1:^NC) - k1(1:^NC)
      uk(1:^NC) = uk(1:^NC) - k1(1+^NC:2*^NC)
      er = norm2(k1)

        if (er .lt. IMRtol) then
          exit
        end if
    end do

    ! Update variables
    x(1:^NC) = 2.d0*xk(1:^NC) - x(1:^NC)
    u(1:^NC) = 2.d0*uk(1:^NC) - u(1:^NC)

    ! Advance
    particle(ipart)%self%x(1:^NC) = x(1:^NC)
    particle(ipart)%self%u(1:^NC) = u(1:^NC)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end select


  ! **************************************************
  ! Payload update
  ! **************************************************
  
  ! To give an example, we set the payload to the interpolated density at 
  ! the particles position.  
  ! In general, it would be better to add the auxilary variable to mod_gridvars 
  ! since then you can use the convenient subroutine get_vec() for this.  
  
!   particle(ipart)%self%payload(1) = tmp_payload(1)

  ! **************************************************
  ! Time update
  ! **************************************************
  particle(ipart)%self%t = particle(ipart)%self%t + dt_p
  if (ipart .eq. 2 .and. (integrator .eq. 2 .or. integrator .eq. 3)) write(*,*) "Particle:",ipart, "Time:",particle(ipart)%self%t, "Final er:",er, "dt:",dt_p

end do

!=============================================================================
end subroutine integrate_particles_geo
!=============================================================================
subroutine get_geo_RK4_funcs(m,x,u,F,H)
! this provides the geodesic functions evaluated at x,u

use mod_metric
use mod_particles
use mod_gridvars, only: igrid_working
include 'amrvacdef.f'
double precision, intent(in), dimension(1:^NC)  :: u, x
double precision, intent(out), dimension(1:^NC) :: F, H
double precision, dimension(1:^NC)              :: tmp
double precision                                :: u0, uu, dummy
integer                                         :: igrid
double precision                                :: tloc
double precision, intent(in)                    :: m
double precision                                :: alpha_p
double precision, dimension(1:^NC)              :: d_alpha_p                
double precision, dimension(1:^NC)              :: betaup_p
double precision, dimension(1:^NC,1:^NC)        :: d_betaup_p
double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
double precision, dimension(1:^NC,1:^NC,1:^NC)  :: d_gammainv_p
integer                                         :: i, j, k
!-----------------------------------------------------------------------------

! Distinguish particle/photon
if (m==0.d0) then ! photon
  uu = 0.d0
else ! particle
  uu = 1.d0
end if

! Get metric functions at position x^i
call get_alpha({x(^D)},alpha_p)
!write(*,*) alpha_p

do i=1,^NC
  call get_beta(i,{x(^D)},betaup_p(i))
end do


do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{x(^D)},gammainv_p(i,j))
  end do
end do


! Compute u0
call get_u0(u,uu,gammainv_p,alpha_p,u0)

! RHS of equations for x
tmp = matmul(gammainv_p,u)
F(1:^NC) = tmp(1:^NC)/u0 - betaup_p(1:^NC)

! Get metric derivatives at position x^i
do i=1,^NC
  call get_alpha({x(^D)},dummy,dalphadj=d_alpha_p(i),jdir=i)
end do

do j=1,^NC
  do i=1,^NC	
    call get_beta(i,{x(^D)},dummy,dbetaidj=d_betaup_p(j,i),jdir=j)
  end do
end do

do k=1,^NC
  do j=1,^NC
    do i=1,^NC
      call get_gammainv_component_analytic(i,j,{x(^D)},dummy,dginvdk=d_gammainv_p(i,j,k),kdir=k)
    end do
  end do
end do

! RHS of equations for u
tmp = matmul(d_betaup_p,u)
H(1:^NC) = -alpha_p*u0*d_alpha_p(1:^NC) + tmp(1:^NC)
do k=1,^NC
  tmp(1:^NC) = matmul(d_gammainv_p(:,:,k),u)
  H(k) = H(k) - sum(tmp*u)/2.d0/u0
end do

!=============================================================================
end subroutine get_geo_RK4_funcs
!=============================================================================
subroutine get_geo_IMR_funcs(m,x,u,F,H)
! this provides the RHS of the geodesic functions evaluated at x,u

use mod_metric
use mod_particles
use mod_gridvars, only: igrid_working
include 'amrvacdef.f'
double precision, intent(in), dimension(1:^NC)  :: u, x
double precision, intent(out), dimension(1:^NC) :: F, H
double precision, dimension(1:^NC)              :: tmp
double precision                                :: u0, uu, dummy
integer                                         :: igrid
double precision                                :: tloc
double precision, intent(in)                    :: m
double precision                                :: alpha_p
double precision, dimension(1:^NC)              :: d_alpha_p                
double precision, dimension(1:^NC)              :: betaup_p
double precision, dimension(1:^NC,1:^NC)        :: d_betaup_p
double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
double precision, dimension(1:^NC,1:^NC,1:^NC)  :: d_gammainv_p
integer                                         :: i, j, k
!-----------------------------------------------------------------------------

! Distinguish particle/photon
if (m==0.d0) then ! photon
  uu = 0.d0
else ! particle
  uu = 1.d0
end if

! Get metric functions at position x^i
call get_alpha({x(^D)},alpha_p)
!write(*,*) alpha_p

do i=1,^NC
  call get_beta(i,{x(^D)},betaup_p(i))
end do


do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{x(^D)},gammainv_p(i,j))
  end do
end do


! Compute u0
call get_u0(u,uu,gammainv_p,alpha_p,u0)

! Rresidual equations for x
tmp = matmul(gammainv_p,u)
F(1:^NC) = tmp(1:^NC)/u0 - betaup_p(1:^NC)

! Get metric derivatives at position x^i
do i=1,^NC
  call get_alpha({x(^D)},dummy,dalphadj=d_alpha_p(i),jdir=i)
end do

do j=1,^NC
  do i=1,^NC	
    call get_beta(i,{x(^D)},dummy,dbetaidj=d_betaup_p(j,i),jdir=j)
  end do
end do

do k=1,^NC
  do j=1,^NC
    do i=1,^NC
      call get_gammainv_component_analytic(i,j,{x(^D)},dummy,dginvdk=d_gammainv_p(i,j,k),kdir=k)
    end do
  end do
end do

! RHS of equations for u
tmp = matmul(d_betaup_p,u)
H(1:^NC) = -alpha_p*u0*d_alpha_p(1:^NC) + tmp(1:^NC)
do k=1,^NC
  tmp(1:^NC) = matmul(d_gammainv_p(:,:,k),u)
  H(k) = H(k) - sum(tmp*u)/2.d0/u0
end do

!=============================================================================
end subroutine get_geo_IMR_funcs
!=============================================================================
subroutine derivs_geodesic(t_s,y,dydt)

use mod_metric
use mod_gridvars
use mod_particles, only: particle
use constants
include 'amrvacdef.f'

double precision                :: t_s, y(*)
double precision                :: dydt(*)
! .. local ..
double precision, dimension(1:^NC)              :: u, x
double precision, dimension(1:^NC)              :: tmp
double precision                                :: u0, uu, dummy
integer                                         :: igrid
double precision                                :: tloc
double precision                                :: m
double precision                                :: alpha_p
double precision, dimension(1:^NC)              :: d_alpha_p                
double precision, dimension(1:^NC)              :: betaup_p
double precision, dimension(1:^NC,1:^NC)        :: d_betaup_p
double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
double precision, dimension(1:^NC,1:^NC,1:^NC)  :: d_gammainv_p
integer                                         :: i, j, k
!-----------------------------------------------------------------------------

m = particle(ipart_working)%self%m
x = y(1:^NC)
u = y(^NC+1:2*^NC)

! Distinguish particle/photon
if (m .eq. 0.d0) then ! photon
  uu = 0.d0
else ! particle
  uu = 1.d0
end if

! Get metric functions at position x^i
call get_alpha({x(^D)},alpha_p)
!write(*,*) alpha_p

do i=1,^NC
  call get_beta(i,{x(^D)},betaup_p(i))
end do


do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{x(^D)},gammainv_p(i,j))
  end do
end do


! Compute u0
call get_u0(u,uu,gammainv_p,alpha_p,u0)

! RHS of equations for x
tmp = matmul(gammainv_p,u)
dydt(1:^NC) = tmp(1:^NC)/u0 - betaup_p(1:^NC)

! Get metric derivatives at position x^i
do i=1,^NC
  call get_alpha({x(^D)},dummy,dalphadj=d_alpha_p(i),jdir=i)
end do

do j=1,^NC
  do i=1,^NC	
    call get_beta(i,{x(^D)},dummy,dbetaidj=d_betaup_p(j,i),jdir=j)
  end do
end do

do k=1,^NC
  do j=1,^NC
    do i=1,^NC
      call get_gammainv_component_analytic(i,j,{x(^D)},dummy,dginvdk=d_gammainv_p(i,j,k),kdir=k)
    end do
  end do
end do

! RHS of equations for u
tmp = matmul(d_betaup_p,u)
dydt(^NC+1:2*^NC) = -alpha_p*u0*d_alpha_p(1:^NC) + tmp(1:^NC)
do k=1,^NC
  tmp(1:^NC) = matmul(d_gammainv_p(:,:,k),u)
  dydt(^NC+k) = dydt(^NC+k) - sum(tmp*u)/2.d0/u0
end do

end subroutine derivs_geodesic
!=============================================================================
subroutine set_particles_dt_geo()

use mod_metric
use mod_particles
include 'amrvacdef.f'

integer                                  :: ipart, iipart, i, j
double precision                         :: t_min_mype, tout, dt_particles_mype, dt_cfl
double precision                         :: v(1:^NC)
double precision                         :: alpha_p, u0, uu
double precision, dimension(1:^NC)       :: betaup_p, tmp
double precision, dimension(1:^NC,1:^NC) :: gammainv_p
double precision                         :: cfl=.1d0, dt0=.01d0
!-----------------------------------------------------------------------------
t_min_mype        = tsync_particles
dt_particles_mype = bigdouble


do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

! Reinitialise dt in case it became too small due to matching with tout
 particle(ipart)%self%dt = dt0*UNIT_LENGTH

! Distinguish particle/photon
if (particle(ipart)%self%m==0.d0) then ! photon
  uu = 0.d0
else ! particle
  uu = 1.d0
end if

! Get metric functions at position x^i
call get_alpha({particle(ipart)%self%x(^D)},alpha_p)

do i=1,^NC
  call get_beta(i,{particle(ipart)%self%x(^D)},betaup_p(i))
end do

do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{particle(ipart)%self%x(^D)},gammainv_p(i,j))
  end do
end do

! Compute u0
call get_u0(particle(ipart)%self%u,uu,gammainv_p,alpha_p,u0)

! make sure we step only one cell at a time:
tmp = matmul(gammainv_p,particle(ipart)%self%u)
{v(^C)   = abs(tmp(^C)-betaup_p(^C)*u0)/u0\}

{#IFNDEF D1
   dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)})
}
{#IFDEF D1
   dt_cfl = rnode(rpdx1_,particle(ipart)%igrid)/v(1)
}

   particle(ipart)%self%dt = min(dt_cfl*cfl,particle(ipart)%self%dt)*UNIT_LENGTH

   particle(ipart)%self%dt = limit_dt(particle(ipart)%self%t,particle(ipart)%self%dt)

   dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)
   t_min_mype = min(t_min_mype,particle(ipart)%self%t)

end do !ipart

call save_tmin(t_min_mype,dt_particles_mype)

end subroutine set_particles_dt_geo
!=============================================================================
subroutine get_u0(u,uu,gammainv_p,alpha_p,u0)

double precision, dimension(1:^NC), intent(in)         :: u
double precision, intent(in)                           :: alpha_p,uu
double precision, dimension(1:^NC,1:^NC), intent(in)   :: gammainv_p
double precision, intent(out)                          :: u0
double precision, dimension(1:^NC)                     :: tmp
!-----------------------------------------------------------------------------

tmp = matmul(gammainv_p,u)
u0 = sqrt(sum(u*tmp) + uu) / alpha_p

end subroutine get_u0
}
}
!=============================================================================
{#IFDEF PARTICLES
{#IFDEF PARTICLES_GRLORENTZ
!=============================================================================
subroutine integrate_particles_grlorentz()
! this integrates the GR-Lorentz equation for charged particles

use mod_metric
use mod_particles
use mod_odeint
use mod_gridvars, only: ipart_working
use mod_ngmres

include 'amrvacdef.f'
integer                             :: ipart, iipart, igridp, i, j
double precision                    :: dt_p,alpha_p, qp, mp, tp, lfac, gg
double precision, dimension(1:^NC)  :: u, x, e, b
double precision, dimension(1:2*^NC):: k1,k2,k3,k4
double precision                    :: tloc
integer                             :: integrator=5 ! 1=RK4, 2=IMR, 3=HAM
! IMR variables
double precision                    :: IMRtol=1.d-14,er=1.d0 ! Absolute iteration tolerance
integer                             :: IMRnmax=100,nk ! Maximum number of iterations
double precision, dimension(1:^NC)  :: xk,uk
double precision, dimension(1:2*^NC) :: xk_IMR,sol_IMR
! HAM variables
double precision                    :: HAMtol=1.d-14, HAMlim=2*sqrt(1.d-16) ! Absolute iteration tolerance
integer                             :: HAMnmax=100 ! Maximum number of iterations
integer,dimension(1:^NC)            :: di
double precision, dimension(1:2*^NC) :: xk_HAM,sol_HAM
double precision, dimension(1:^NC)  :: xm,um
! nsolgm variables
double precision                    :: abs_tol,rel_tol,etamax
integer                             :: ierr
!!!!! Precision of time-integration: !!!!!!!!!!!!!
double precision,parameter          :: eps=1.0d-3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, parameter                  :: nvar = 2*^NC
double precision,dimension(1:nvar)  :: y

common /resgrl/ x,u,mp,qp,tp,dt_p,igridp

external derivs_geodesic, res_grl_HAM, res_grl_IMR, res_geo_IMR, res_geo_HAM
!-----------------------------------------------------------------------------

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

  dt_p = particle(ipart)%self%dt
  x(1:^NC) = particle(ipart)%self%x(1:^NC)
  u(1:^NC) = particle(ipart)%self%u(1:^NC)
  qp = particle(ipart)%self%q
  mp = particle(ipart)%self%m
  tp = particle(ipart)%self%t
  igridp = particle(ipart)%igrid
!  call set_tmpGlobals(igrid)
   

  ! **************************************************
  ! Position and velocity update
  ! **************************************************

  select case(integrator)
  case(1) ! RK4
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Runge-Kutta 4th order:
   
    ! First sweep
    call get_grl_RHS_funcs(mp,qp,igridp,tp,x,u,k1(1:^NC),k1(1+^NC:2*^NC))
    y(1:^NC) = x(1:^NC) + dt_p/2.d0 * k1(1:^NC)
    y(1+^NC:2*^NC) = u(1:^NC) + dt_p/2.d0 * k1(1+^NC:2*^NC)
    ! Second sweep
    call get_grl_RHS_funcs(mp,qp,igridp,tp,y(1:^NC),y(1+^NC:2*^NC),k2(1:^NC),k2(1+^NC:2*^NC))
    y(1:^NC) = x(1:^NC) + dt_p/2.d0 * k2(1:^NC)
    y(1+^NC:2*^NC) = u(1:^NC) + dt_p/2.d0 * k2(1+^NC:2*^NC)
    ! Third sweep
    call get_grl_RHS_funcs(mp,qp,igridp,tp,y(1:^NC),y(1+^NC:2*^NC),k3(1:^NC),k3(1+^NC:2*^NC))
    y(1:^NC) = x(1:^NC) + dt_p * k3(1:^NC)
    y(1+^NC:2*^NC) = u(1:^NC) + dt_p * k3(1+^NC:2*^NC) 
    ! Fourth sweep
    call get_grl_RHS_funcs(mp,qp,igridp,tp,y(1:^NC),y(1+^NC:2*^NC),k4(1:^NC),k4(1+^NC:2*^NC))

    ! Advance
    particle(ipart)%self%x(1:^NC) = x(1:^NC) + dt_p/6.d0 * (k1(1:^NC)+2.d0*k2(1:^NC)+2.d0*k3(1:^NC)+k4(1:^NC))
    particle(ipart)%self%u(1:^NC) = u(1:^NC) + dt_p/6.d0 * (k1(1+^NC:2*^NC)+2.d0*k2(1+^NC:2*^NC)+2.d0*k3(1+^NC:2*^NC)+k4(1+^NC:2*^NC))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case(2) ! IMR (fixed-point iteration)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2nd-order Implicit Midpoint Rule (symplectic)
    ! Uses fixed-point iteration (full Jacobian requires some restructuring)

    ! Iteration variables
    xk(1:^NC) = x(1:^NC)
    uk(1:^NC) = u(1:^NC)

    ! Preconditioning with 1st order Euler
    call get_grl_RHS_funcs(mp,qp,igridp,tp,x,u,k1(1:^NC),k1(1+^NC:2*^NC))
    xk(1:^NC) = x(1:^NC) + dt_p/2.d0 * k1(1:^NC)
    uk(1:^NC) = u(1:^NC) + dt_p/2.d0 * k1(1+^NC:2*^NC)

    ! Start nonlinear iteration
    do nk = 1,IMRnmax
      ! Compute residuals
      call get_grl_RHS_funcs(mp,qp,igridp,tp+dt_p/2.d0,xk,uk,k1(1:^NC),k1(1+^NC:2*^NC)) ! Get the RHS as a function of xk,uk
      k1(1:^NC) = xk(1:^NC) - x(1:^NC) - dt_p/2.d0*k1(1:^NC) ! F = xk - x - dt/2*RHS(xk,uk)
      k1(1+^NC:2*^NC) = uk(1:^NC) - u(1:^NC) - dt_p/2.d0*k1(1+^NC:2*^NC) ! H = uk - u - dt/2*RHS(xk,uk)

      ! Update and error check
      xk(1:^NC) = xk(1:^NC) - k1(1:^NC)
      uk(1:^NC) = uk(1:^NC) - k1(1+^NC:2*^NC)
      er = norm2(k1)
      if (er .lt. IMRtol) then
        exit
      end if
    end do

    ! Update variables
    x(1:^NC) = 2.d0*xk(1:^NC) - x(1:^NC)
    u(1:^NC) = 2.d0*uk(1:^NC) - u(1:^NC)

    ! Advance
    particle(ipart)%self%x(1:^NC) = x(1:^NC)
    particle(ipart)%self%u(1:^NC) = u(1:^NC)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case(3) ! HAM (fixed-point iteration)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2nd-order implicit energy-conserving Hamiltonian scheme
    ! Uses fixed-point iteration (full Jacobian requires some restructuring)

    ! Iteration variables
    xk(1:^NC) = x(1:^NC)
    uk(1:^NC) = u(1:^NC)

    ! Preconditioning with 1st order Euler
    call get_grl_RHS_funcs(mp,qp,igridp,tp,x,u,k1(1:^NC),k1(1+^NC:2*^NC))
    xk(1:^NC) = x(1:^NC) + dt_p * k1(1:^NC)
    uk(1:^NC) = u(1:^NC) + dt_p * k1(1+^NC:2*^NC)

    ! Start nonlinear iteration
    do nk = 1,HAMnmax
      ! Compute residuals
      call get_grl_HAM_funcs(mp,qp,igridp,tp,dt_p,x,u,xk,uk,k1(1:^NC),k1(1+^NC:2*^NC)) ! Get the RHS as a function of xk,uk
      k1(1:^NC) = xk(1:^NC) - x(1:^NC) - dt_p*k1(1:^NC) ! F = xk - x - dt*RHS(xk,uk)
      k1(1+^NC:2*^NC) = uk(1:^NC) - u(1:^NC) - dt_p*k1(1+^NC:2*^NC) ! H = uk - u - dt*RHS(xk,uk)

      ! Update and error check
      xk(1:^NC) = xk(1:^NC) - k1(1:^NC)
      uk(1:^NC) = uk(1:^NC) - k1(1+^NC:2*^NC)
      er = norm2(k1)
      if (er .lt. HAMtol) then
        exit
      end if
    end do

    ! Update variables
    x(1:^NC) = xk(1:^NC)
    u(1:^NC) = uk(1:^NC)

    ! Advance
    particle(ipart)%self%x(1:^NC) = x(1:^NC)
    particle(ipart)%self%u(1:^NC) = u(1:^NC)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case(4) ! IMR (NK iteration with nsolgm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2nd-order Implicit Midpoint Rule (symplectic)

    ! Preconditioning with 1st order Euler
    call get_grl_RHS_funcs(mp,qp,igridp,tp,x,u,k1(1:^NC),k1(1+^NC:2*^NC))
    xk_IMR(1:^NC) = x(1:^NC) + dt_p/2.d0 * k1(1:^NC)
    xk_IMR(1+^NC:2*^NC) = u(1:^NC) + dt_p/2.d0 * k1(1+^NC:2*^NC)

    er = 1.d0
    ierr = 1
    abs_tol = IMRtol
    rel_tol = IMRtol
    etamax = 0.9d0

    ! Nonlinear iteration with nsolgm
    call nsolgm(xk_IMR, sol_IMR, 6, ierr, er, abs_tol, rel_tol, etamax, IMRnmax, 40, res_grl_IMR)

    ! Update variables
    x(1:^NC) = 2.d0*sol_IMR(1:^NC) - x(1:^NC)
    u(1:^NC) = 2.d0*sol_IMR(1+^NC:2*^NC) - u(1:^NC)

    ! Advance
    particle(ipart)%self%x(1:^NC) = x(1:^NC)
    particle(ipart)%self%u(1:^NC) = u(1:^NC)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  case(5) ! HAM (NK iteration with nsolgm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2nd-order implicit energy-conserving Hamiltonian scheme

    ! Preconditioning with 1st order Euler
    call get_grl_RHS_funcs(mp,qp,igridp,tp,x,u,k1(1:^NC),k1(1+^NC:2*^NC))
    xk_HAM(1:^NC) = x(1:^NC) + dt_p * k1(1:^NC)
    xk_HAM(1+^NC:2*^NC) = u(1:^NC) + dt_p * k1(1+^NC:2*^NC)

    er = 1.d0
    ierr = 1
    abs_tol = HAMtol
    rel_tol = HAMtol
    etamax = 0.9d0

    ! Nonlinear iteration with nsolgm
    call nsolgm(xk_HAM, sol_HAM, 6, ierr, er, abs_tol, rel_tol, etamax, HAMnmax, 40, res_grl_HAM)

    ! Update variables
    x(1:^NC) = sol_HAM(1:^NC)
    u(1:^NC) = sol_HAM(1+^NC:2*^NC)

    ! Advance
    particle(ipart)%self%x(1:^NC) = x(1:^NC)
    particle(ipart)%self%u(1:^NC) = u(1:^NC)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end select

  ! **************************************************
  ! Payload update: store lfac
  ! **************************************************
  
  lfac = 0.d0
  do i=1,^NC
    do j=1,^NC
      call get_gammainv_component_analytic(i,j,{particle(ipart)%self%x(^D)},gg)
      lfac = lfac + gg*particle(ipart)%self%u(j)*particle(ipart)%self%u(i)
    end do
  end do
  lfac = dsqrt(lfac + 1.d0)
  particle(ipart)%self%payload(1) = lfac

  ! **************************************************
  ! Time update
  ! **************************************************
  particle(ipart)%self%t = particle(ipart)%self%t + dt_p
!  if (ipart .eq. 2 .and. (integrator .eq. 2 .or. integrator .eq. 3)) &
!    write(*,*) "Particle:",ipart, "Time:",particle(ipart)%self%t, "Final er:",er, "Iterations:",nk, "dt:",dt_p

end do

!=============================================================================
end subroutine integrate_particles_grlorentz
!=============================================================================
subroutine get_grl_RHS_funcs(m,q,igridp,tp,x,u,F,H)
! this provides the right-hand side of the GR-Lorentz functions evaluated at x,u

use mod_metric
use mod_particles
use mod_gridvars, only: igrid_working
include 'amrvacdef.f'
double precision, intent(in), dimension(1:^NC)  :: u, x
double precision, intent(out), dimension(1:^NC) :: F, H
double precision, dimension(1:^NC)              :: tmp
double precision                                :: u0, uu, dummy
double precision                                :: tloc
double precision, intent(in)                    :: m, q, tp
integer, intent(in)                             :: igridp
! Variables for the geodesic part
double precision                                :: alpha_p
double precision, dimension(1:^NC)              :: d_alpha_p                
double precision, dimension(1:^NC)              :: betaup_p
double precision, dimension(1:^NC,1:^NC)        :: d_betaup_p
double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
double precision, dimension(1:^NC,1:^NC,1:^NC)  :: d_gammainv_p
integer                                         :: i, j, k
! Variables for the Lorentz force part
double precision                                :: qom, cross_D
double precision, dimension(1:^NC)              :: e, e_D, b, uxB_D
!-----------------------------------------------------------------------------

! Massive particles
uu = 1.d0 ! norm of the 4-velocity
qom = q/m ! Charge-to-mass ratio

! Get metric functions at position x^i
call get_alpha({x(^D)},alpha_p)
!write(*,*) alpha_p

do i=1,^NC
  call get_beta(i,{x(^D)},betaup_p(i))
end do


do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{x(^D)},gammainv_p(i,j))
  end do
end do

! Compute u0
call get_u0(u,uu,gammainv_p,alpha_p,u0)

! RHS of equations for x
tmp = matmul(gammainv_p,u)
F(1:^NC) = tmp(1:^NC)/u0 - betaup_p(1:^NC)

! Get metric derivatives at position x^i
do i=1,^NC
  call get_alpha({x(^D)},dummy,dalphadj=d_alpha_p(i),jdir=i)
end do

do j=1,^NC
  do i=1,^NC	
    call get_beta(i,{x(^D)},dummy,dbetaidj=d_betaup_p(j,i),jdir=j)
  end do
end do

do k=1,^NC
  do j=1,^NC
    do i=1,^NC
      call get_gammainv_component_analytic(i,j,{x(^D)},dummy,dginvdk=d_gammainv_p(i,j,k),kdir=k)
    end do
  end do
end do

! RHS of equations for u
tmp = matmul(d_betaup_p,u)
H(1:^NC) = -alpha_p*u0*d_alpha_p(1:^NC) + tmp(1:^NC)
do k=1,^NC
  tmp(1:^NC) = matmul(d_gammainv_p(:,:,k),u)
  H(k) = H(k) - sum(tmp*u)/2.d0/u0
end do

! Add Lorentz force
call get_e(igridp,x,tp,e) ! Interpolate e field
call lower3_point({x(^D)},e,e_D) !! gamma_ij D^j

call get_b(igridp,x,tp,b) ! Interpolate b field
do i=1,^NC
  uxB_D(i) = cross_D(matmul(gammainv_p,u)/u0,b,x,i) !! e_ijk gamma^jl u_l/u^0 b^k
end do

H = H + qom * (alpha_p*e_D + uxB_D) !! q/m (alpha D_i + e_ijk gamma^jl u_l/u^0 b^k)

!=============================================================================
end subroutine get_grl_RHS_funcs
!=============================================================================
subroutine get_grl_HAM_funcs(m,q,igridp,tp,dtp,xn,un,x,u,F,H)
! this provides the right-hand side of the GR-Lorentz functions evaluated at x,u

use mod_metric
use mod_particles
use mod_gridvars, only: igrid_working
include 'amrvacdef.f'
double precision, intent(in), dimension(1:^NC)  :: u, x, un, xn
double precision, intent(out), dimension(1:^NC) :: F, H
double precision, dimension(1:^NC)              :: tmp
double precision                                :: u0, uu, dummy
double precision                                :: tloc
double precision, intent(in)                    :: m, q, tp, dtp
integer, intent(in)                             :: igridp
! Variables for the geodesic part
double precision                                :: alpha_p
double precision, dimension(1:^NC)              :: d_alpha_p                
double precision, dimension(1:^NC)              :: betaup_p
double precision, dimension(1:^NC,1:^NC)        :: d_betaup_p
double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
double precision, dimension(1:^NC,1:^NC,1:^NC)  :: d_gammainv_p
integer                                         :: i, j, k
! Variables for the Lorentz force part
double precision                                :: qom, cross_D, HF, HH
double precision, dimension(1:^NC)              :: e, e_D, b, uxB_D, xm, bxB_D
!-----------------------------------------------------------------------------

! Massive particles
uu = 1.d0 ! norm of the 4-velocity
qom = q/m ! Charge-to-mass ratio

!!!!! GEODESIC PART
! Build terms of the Hamiltonian integrator
! RHS of equations for x
F(1) = (HF(x(1),xn(2),xn(3),u(1),un(2),un(3),un(1),1) &
       + HF(x(1),x(2),x(3),u(1),u(2),u(3),un(1),1) &
       + HF(xn(1),x(2),x(3),u(1),u(2),u(3),un(1),1) &
       + HF(x(1),xn(2),x(3),u(1),un(2),u(3),un(1),1) &
       + HF(xn(1),xn(2),xn(3),u(1),un(2),un(3),un(1),1) &
       + HF(xn(1),xn(2),x(3),u(1),un(2),u(3),un(1),1))/6.d0
F(2) = (HF(x(1),x(2),xn(3),u(1),u(2),un(3),un(2),2) &
       + HF(xn(1),x(2),xn(3),un(1),u(2),un(3),un(2),2) &
       + HF(xn(1),xn(2),xn(3),un(1),u(2),un(3),un(2),2) &
       + HF(x(1),x(2),x(3),u(1),u(2),u(3),un(2),2) &
       + HF(x(1),xn(2),xn(3),u(1),u(2),un(3),un(2),2) &
       + HF(x(1),xn(2),x(3),u(1),u(2),u(3),un(2),2))/6.d0
F(3) = (HF(x(1),x(2),x(3),u(1),u(2),u(3),un(3),3) &
       + HF(xn(1),x(2),x(3),un(1),u(2),u(3),un(3),3) &
       + HF(xn(1),x(2),xn(3),un(1),u(2),u(3),un(3),3) &
       + HF(xn(1),xn(2),x(3),un(1),un(2),u(3),un(3),3) &
       + HF(x(1),x(2),xn(3),u(1),u(2),u(3),un(3),3) &
       + HF(xn(1),xn(2),xn(3),un(1),un(2),u(3),un(3),3))/6.d0
! RHS of equations for u
H(1) = (HH(x(1),xn(2),xn(3),un(1),un(2),un(3),xn(1),1) &
       + HH(x(1),x(2),x(3),un(1),u(2),u(3),xn(1),1) &
       + HH(x(1),x(2),x(3),u(1),u(2),u(3),xn(1),1) &
       + HH(x(1),xn(2),x(3),un(1),un(2),u(3),xn(1),1) &
       + HH(x(1),xn(2),xn(3),u(1),un(2),un(3),xn(1),1) &
       + HH(x(1),xn(2),x(3),u(1),un(2),u(3),xn(1),1))/6.d0
H(2) = (HH(x(1),x(2),xn(3),u(1),un(2),un(3),xn(2),2) &
       + HH(xn(1),x(2),xn(3),un(1),un(2),un(3),xn(2),2) &
       + HH(xn(1),x(2),xn(3),un(1),u(2),un(3),xn(2),2) &
       + HH(x(1),x(2),x(3),u(1),un(2),u(3),xn(2),2) &
       + HH(x(1),x(2),xn(3),u(1),u(2),un(3),xn(2),2) &
       + HH(x(1),x(2),x(3),u(1),u(2),u(3),xn(2),2))/6.d0
H(3) = (HH(x(1),x(2),x(3),u(1),u(2),un(3),xn(3),3) &
       + HH(xn(1),x(2),x(3),un(1),u(2),un(3),xn(3),3) &
       + HH(xn(1),x(2),x(3),un(1),u(2),u(3),xn(3),3) &
       + HH(xn(1),xn(2),x(3),un(1),un(2),un(3),xn(3),3) &
       + HH(x(1),x(2),x(3),u(1),u(2),u(3),xn(3),3) &
       + HH(xn(1),xn(2),x(3),un(1),un(2),u(3),xn(3),3))/6.d0

!!!!! LORENTZ FORCE PART
! Get metric functions at position xbar
xm = .5d0 * (x + xn)
! alpha
call get_alpha({xm(^D)},alpha_p)
! beta^i
do i=1,^NC
  call get_beta(i,{xm(^D)},betaup_p(i))
end do
! gamma^ij
do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{xm(^D)},gammainv_p(i,j))
  end do
end do

! Get fields at xbar
call get_e(igridp,xm,tp+dtp/2.d0,e) ! Interpolate e field
call get_b(igridp,x,tp+dtp/2.d0,b) ! Interpolate b field

! Assemble
call lower3_point({xm(^D)},e,e_D) !! gamma_ij D^j
do i=1,^NC
  bxB_D(i) = cross_D(betaup_p,b,xm,i) !! e_ijk beta^j b^k
end do
do i=1,^NC
  uxB_D(i) = cross_D(F,b,xm,i) !! e_ijk dx^j/dt b^k
end do

H = H + qom * (alpha_p*e_D + bxB_D + uxB_D) !! q/m (alpha D_i + e_ijk beta^j b^k + e_ijk dx^j/dt b^k)

!=============================================================================
end subroutine get_grl_HAM_funcs
!=============================================================================
double precision function HF(x1,x2,x3,u1,u2,u3,uni,ii)
  ! Computes one term of the ii-th position equation for the Hamiltonian scheme

  use mod_metric
  include 'amrvacdef.f'

  double precision, intent(in)                    :: x1, x2, x3, u1, u2, u3, uni
  integer, intent(in)                             :: ii
  double precision, dimension(1:^NC)              :: x, u, un
  double precision                                :: alpha_p, betaup_p
  double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
  double precision                                :: tmp
  integer                                         :: i, j, k
  !-----------------------------------------------------------------------------

  ! Get metric functions at position x=(x1,x2,x3)
  x = (/ x1, x2, x3 /)
  ! alpha
  call get_alpha({x(^D)},alpha_p)
  ! beta^ii
  call get_beta(ii,{x(^D)},betaup_p)
  ! gamma^ij
  do i = 1,ndir
    do j = 1,ndir
      call get_gammainv_component_analytic(i,j,{x(^D)},gammainv_p(i,j))
    end do
  end do

  ! Momentum vectors
  u = (/ u1, u2, u3 /)
  un = u
  un(ii) = uni

  ! Build term
  HF = dot_product(gammainv_p(ii,:),u+un)
  tmp = sqrt(dot_product(matmul(gammainv_p,u),u) + 1.d0) &
        + sqrt(dot_product(matmul(gammainv_p,un),un) + 1.d0)
  HF = alpha_p * HF / tmp - betaup_p

end function HF
!=============================================================================
double precision function HH(x1,x2,x3,u1,u2,u3,xni,ii)
  ! Computes one term of the ii-th momentum equation for the Hamiltonian scheme
  ! In case di = 1, returns the analytic derivative-version of that term

  use mod_metric
  include 'amrvacdef.f'

  double precision, intent(in)                    :: x1, x2, x3, u1, u2, u3, xni
  integer, intent(in)                             :: ii
  double precision, dimension(1:^NC)              :: x, u, xn
  double precision                                :: alpha_p, alpha_pn, d_alpha_p, dummy
  double precision, dimension(1:^NC)              :: betaup_p, betaup_pn, d_betaup_p
  double precision, dimension(1:^NC,1:^NC)        :: gammainv_p, gammainv_pn, d_gammainv_p
  double precision                                :: tmp, dilimit=2*sqrt(1.d-16) ! Limit for the precision of finite differences
  integer                                         :: i, j, k, di=0
  !-----------------------------------------------------------------------------
  
  ! Position vectors
  x = (/ x1, x2, x3 /)
  xn = x
  xn(ii) = xni

  ! Detect vanishing position increments
  if (abs(x(ii) - xni) .le. dilimit) then
    di = 1
  end if

  ! Get metric functions at position x=(x1,x2,x3), xn
  ! alpha
  call get_alpha({x(^D)},alpha_p)
  call get_alpha({xn(^D)},alpha_pn)
  ! beta^i
  do i=1,^NC
    call get_beta(i,{x(^D)},betaup_p(i))
    call get_beta(i,{xn(^D)},betaup_pn(i))
  end do
  ! gamma^ij
  do i = 1,ndir
    do j = 1,ndir
      call get_gammainv_component_analytic(i,j,{x(^D)},gammainv_p(i,j))
      call get_gammainv_component_analytic(i,j,{xn(^D)},gammainv_pn(i,j))
    end do
  end do

  ! Momentum vector
  u = (/ u1, u2, u3 /)

  ! Build term
  if (di .eq. 0) then ! Finite-difference case
    tmp = sqrt(dot_product(matmul(gammainv_p,u),u) + 1.d0) &
          + sqrt(dot_product(matmul(gammainv_pn,u),u) + 1.d0)
    HH = - .5d0 * tmp * (alpha_p - alpha_pn) / (x(ii) - xni) &
         + dot_product(u,betaup_p-betaup_pn) / (x(ii) - xni) &
         - .5d0 * (alpha_p + alpha_pn) / tmp &
                * dot_product(matmul(gammainv_p-gammainv_pn,u),u) / (x(ii) - xni)
  else ! Analytic derivative case
    ! Get metric derivatives at position x^i
    call get_alpha({xn(^D)},dummy,dalphadj=d_alpha_p,jdir=ii) ! dalpha/dx^ii
    do i=1,^NC	
      call get_beta(i,{xn(^D)},dummy,dbetaidj=d_betaup_p(i),jdir=ii) ! dbeta^i/dx^ii
    end do
    do j=1,^NC
      do i=1,^NC
        call get_gammainv_component_analytic(i,j,{xn(^D)},dummy,dginvdk=d_gammainv_p(i,j),kdir=ii) ! dgamma^ij/dx^ii
      end do
    end do

    tmp = sqrt(dot_product(matmul(gammainv_p,u),u) + 1.d0) &
          + sqrt(dot_product(matmul(gammainv_pn,u),u) + 1.d0)
    HH = - .5d0 * tmp * d_alpha_p &
         + dot_product(u,d_betaup_p) &
         - .5d0 * (alpha_p + alpha_pn) / tmp &
                * dot_product(matmul(d_gammainv_p,u),u)
  end if
 
end function HH
!=============================================================================
subroutine get_geo_RHS_funcs(x,u,F,H)
! this provides the RHS of the geodesic functions evaluated at x,u

use mod_metric
use mod_particles
use mod_gridvars, only: igrid_working
include 'amrvacdef.f'
double precision, intent(in), dimension(1:^NC)  :: u, x
double precision, intent(out), dimension(1:^NC) :: F, H
double precision, dimension(1:^NC)              :: tmp
double precision                                :: u0, uu, dummy
integer                                         :: igrid
double precision                                :: tloc
double precision                                :: alpha_p
double precision, dimension(1:^NC)              :: d_alpha_p                
double precision, dimension(1:^NC)              :: betaup_p
double precision, dimension(1:^NC,1:^NC)        :: d_betaup_p
double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
double precision, dimension(1:^NC,1:^NC,1:^NC)  :: d_gammainv_p
integer                                         :: i, j, k
!-----------------------------------------------------------------------------

! Massive particles
uu = 1.d0 ! norm of the 4-velocity

! Get metric functions at position x^i
call get_alpha({x(^D)},alpha_p)
!write(*,*) alpha_p

do i=1,^NC
  call get_beta(i,{x(^D)},betaup_p(i))
end do


do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{x(^D)},gammainv_p(i,j))
  end do
end do


! Compute u0
call get_u0(u,uu,gammainv_p,alpha_p,u0)

! Rresidual equations for x
tmp = matmul(gammainv_p,u)
F(1:^NC) = tmp(1:^NC)/u0 - betaup_p(1:^NC)

! Get metric derivatives at position x^i
do i=1,^NC
  call get_alpha({x(^D)},dummy,dalphadj=d_alpha_p(i),jdir=i)
end do

do j=1,^NC
  do i=1,^NC	
    call get_beta(i,{x(^D)},dummy,dbetaidj=d_betaup_p(j,i),jdir=j)
  end do
end do

do k=1,^NC
  do j=1,^NC
    do i=1,^NC
      call get_gammainv_component_analytic(i,j,{x(^D)},dummy,dginvdk=d_gammainv_p(i,j,k),kdir=k)
    end do
  end do
end do

! RHS of equations for u
tmp = matmul(d_betaup_p,u)
H(1:^NC) = -alpha_p*u0*d_alpha_p(1:^NC) + tmp(1:^NC)
do k=1,^NC
  tmp(1:^NC) = matmul(d_gammainv_p(:,:,k),u)
  H(k) = H(k) - sum(tmp*u)/2.d0/u0
end do

end subroutine get_geo_RHS_funcs
!=============================================================================
subroutine get_geo_HAM_funcs(xn,un,x,u,F,H)
! this provides the right-hand side of the GR-geodesic functions evaluated at x,u

use mod_metric
use mod_particles
use mod_gridvars, only: igrid_working
include 'amrvacdef.f'
double precision, intent(in), dimension(1:^NC)  :: u, x, un, xn
double precision, intent(out), dimension(1:^NC) :: F, H
double precision, dimension(1:^NC)              :: tmp
double precision                                :: u0, uu, dummy
double precision                                :: tloc, HF, HH
! Variables for the geodesic part
double precision                                :: alpha_p
double precision, dimension(1:^NC)              :: d_alpha_p                
double precision, dimension(1:^NC)              :: betaup_p
double precision, dimension(1:^NC,1:^NC)        :: d_betaup_p
double precision, dimension(1:^NC,1:^NC)        :: gammainv_p
double precision, dimension(1:^NC,1:^NC,1:^NC)  :: d_gammainv_p
integer                                         :: i, j, k
!-----------------------------------------------------------------------------

! Massive particles
uu = 1.d0 ! norm of the 4-velocity

! Build terms of the Hamiltonian integrator
! RHS of equations for x
F(1) = (HF(x(1),xn(2),xn(3),u(1),un(2),un(3),un(1),1) &
       + HF(x(1),x(2),x(3),u(1),u(2),u(3),un(1),1) &
       + HF(xn(1),x(2),x(3),u(1),u(2),u(3),un(1),1) &
       + HF(x(1),xn(2),x(3),u(1),un(2),u(3),un(1),1) &
       + HF(xn(1),xn(2),xn(3),u(1),un(2),un(3),un(1),1) &
       + HF(xn(1),xn(2),x(3),u(1),un(2),u(3),un(1),1))/6.d0
F(2) = (HF(x(1),x(2),xn(3),u(1),u(2),un(3),un(2),2) &
       + HF(xn(1),x(2),xn(3),un(1),u(2),un(3),un(2),2) &
       + HF(xn(1),xn(2),xn(3),un(1),u(2),un(3),un(2),2) &
       + HF(x(1),x(2),x(3),u(1),u(2),u(3),un(2),2) &
       + HF(x(1),xn(2),xn(3),u(1),u(2),un(3),un(2),2) &
       + HF(x(1),xn(2),x(3),u(1),u(2),u(3),un(2),2))/6.d0
F(3) = (HF(x(1),x(2),x(3),u(1),u(2),u(3),un(3),3) &
       + HF(xn(1),x(2),x(3),un(1),u(2),u(3),un(3),3) &
       + HF(xn(1),x(2),xn(3),un(1),u(2),u(3),un(3),3) &
       + HF(xn(1),xn(2),x(3),un(1),un(2),u(3),un(3),3) &
       + HF(x(1),x(2),xn(3),u(1),u(2),u(3),un(3),3) &
       + HF(xn(1),xn(2),xn(3),un(1),un(2),u(3),un(3),3))/6.d0
! RHS of equations for u
H(1) = (HH(x(1),xn(2),xn(3),un(1),un(2),un(3),xn(1),1) &
       + HH(x(1),x(2),x(3),un(1),u(2),u(3),xn(1),1) &
       + HH(x(1),x(2),x(3),u(1),u(2),u(3),xn(1),1) &
       + HH(x(1),xn(2),x(3),un(1),un(2),u(3),xn(1),1) &
       + HH(x(1),xn(2),xn(3),u(1),un(2),un(3),xn(1),1) &
       + HH(x(1),xn(2),x(3),u(1),un(2),u(3),xn(1),1))/6.d0
H(2) = (HH(x(1),x(2),xn(3),u(1),un(2),un(3),xn(2),2) &
       + HH(xn(1),x(2),xn(3),un(1),un(2),un(3),xn(2),2) &
       + HH(xn(1),x(2),xn(3),un(1),u(2),un(3),xn(2),2) &
       + HH(x(1),x(2),x(3),u(1),un(2),u(3),xn(2),2) &
       + HH(x(1),x(2),xn(3),u(1),u(2),un(3),xn(2),2) &
       + HH(x(1),x(2),x(3),u(1),u(2),u(3),xn(2),2))/6.d0
H(3) = (HH(x(1),x(2),x(3),u(1),u(2),un(3),xn(3),3) &
       + HH(xn(1),x(2),x(3),un(1),u(2),un(3),xn(3),3) &
       + HH(xn(1),x(2),x(3),un(1),u(2),u(3),xn(3),3) &
       + HH(xn(1),xn(2),x(3),un(1),un(2),un(3),xn(3),3) &
       + HH(x(1),x(2),x(3),u(1),u(2),u(3),xn(3),3) &
       + HH(xn(1),xn(2),x(3),un(1),un(2),u(3),xn(3),3))/6.d0

end subroutine get_geo_HAM_funcs
!=============================================================================
subroutine res_grl_HAM(xk, res, n)
! this provides the residual functions for the mHAM scheme to be passed to nsolgm

include 'amrvacdef.f'

integer, intent(in)                         :: n
double precision, intent(in), dimension(n)  :: xk
double precision, intent(out), dimension(n) :: res
! Common variables
double precision, dimension(1:^NC)          :: x,u
double precision                            :: mp,qp,tp,dt_p
integer                                     :: igridp

common /resgrl/ x,u,mp,qp,tp,dt_p,igridp
!-----------------------------------------------------------------------------

! Get the RHS and store it in res
call get_grl_HAM_funcs(mp,qp,igridp,tp,dt_p,x,u,xk(1:^NC),xk(1+^NC:2*^NC),res(1:^NC),res(1+^NC:2*^NC))

! Finalise
res(1:^NC) = xk(1:^NC) - x(1:^NC) - dt_p*res(1:^NC) ! F = xk - x - dt*RHS(xk,uk)
res(1+^NC:2*^NC) = xk(1+^NC:2*^NC) - u(1:^NC) - dt_p*res(1+^NC:2*^NC) ! H = uk - u - dt*RHS(xk,uk)

end subroutine res_grl_HAM
!=============================================================================
subroutine res_grl_IMR(xk, res, n)
! this provides the residual functions for the IMR scheme to be passed to nsolgm

include 'amrvacdef.f'

integer, intent(in)                         :: n
double precision, intent(in), dimension(n)  :: xk
double precision, intent(out), dimension(n) :: res
! Common variables
double precision, dimension(1:^NC)          :: x,u
double precision                            :: mp,qp,tp,dt_p
integer                                     :: igridp

common /resgrl/ x,u,mp,qp,tp,dt_p,igridp
!-----------------------------------------------------------------------------

! Get the RHS and store it in res
call get_grl_RHS_funcs(mp,qp,igridp,tp+dt_p/2.d0,xk(1:^NC),xk(1+^NC:2*^NC),res(1:^NC),res(1+^NC:2*^NC))

! Finalise
res(1:^NC) = xk(1:^NC) - x(1:^NC) - dt_p/2.d0*res(1:^NC) ! F = xk - x - dt/2*RHS(xk,uk)
res(1+^NC:2*^NC) = xk(1+^NC:2*^NC) - u(1:^NC) - dt_p/2.d0*res(1+^NC:2*^NC) ! H = uk - u - dt/2*RHS(xk,uk)

end subroutine res_grl_IMR
!=============================================================================
subroutine res_geo_IMR(xk, res, n)
! this provides the residual functions for the mHAM scheme to be passed to nsolgm

include 'amrvacdef.f'

integer, intent(in)                         :: n
double precision, intent(in), dimension(n)  :: xk
double precision, intent(out), dimension(n) :: res
! Common variables
double precision, dimension(1:^NC)          :: x,u
double precision                            :: mp,qp,tp,dt_p
integer                                     :: igridp

common /resgrl/ x,u,mp,qp,tp,dt_p,igridp
!-----------------------------------------------------------------------------

! Get the RHS and store it in res
call get_geo_RHS_funcs(xk(1:^NC),xk(1+^NC:2*^NC),res(1:^NC),res(1+^NC:2*^NC))

! Finalise
res(1:^NC) = xk(1:^NC) - x(1:^NC) - dt_p/2.d0*res(1:^NC) ! F = xk - x - dt*RHS(xk,uk)
res(1+^NC:2*^NC) = xk(1+^NC:2*^NC) - u(1:^NC) - dt_p/2.d0*res(1+^NC:2*^NC) ! H = uk - u - dt*RHS(xk,uk)

end subroutine res_geo_IMR
!=============================================================================
subroutine res_geo_HAM(xk, res, n)
! this provides the residual functions for the geodesic Ham scheme to be passed to nsolgm

include 'amrvacdef.f'

integer, intent(in)                         :: n
double precision, intent(in), dimension(n)  :: xk
double precision, intent(out), dimension(n) :: res
! Common variables
double precision, dimension(1:^NC)          :: x,u
double precision                            :: mp,qp,tp,dt_p
integer                                     :: igridp

common /resgrl/ x,u,mp,qp,tp,dt_p,igridp
!-----------------------------------------------------------------------------

! Get the RHS and store it in res
call get_geo_HAM_funcs(x,u,xk(1:^NC),xk(1+^NC:2*^NC),res(1:^NC),res(1+^NC:2*^NC))

! Finalise
res(1:^NC) = xk(1:^NC) - x(1:^NC) - dt_p*res(1:^NC) ! F = xk - x - dt*RHS(xk,uk)
res(1+^NC:2*^NC) = xk(1+^NC:2*^NC) - u(1:^NC) - dt_p*res(1+^NC:2*^NC) ! H = uk - u - dt*RHS(xk,uk)

end subroutine res_geo_HAM
!=============================================================================
double precision function cross_D(a,b,x,idir)
  ! Takes two contravariant vectors and produces the single component idir of the covariant cross product

  use mod_metric, only: get_sqrtgamma
  include 'amrvacdef.f'

  double precision, intent(in)                 :: a(1:^NC), b(1:^NC)
  double precision, intent(in)                 :: x(1:^NC)
  integer, intent(in)                          :: idir

  ! .. local ..
  integer                                      :: i, j, k
  double precision                             :: sqrtgamma
  !-----------------------------------------------------------------------------

  call get_sqrtgamma({x(^D)},sqrtgamma)
  
  cross_D = zero
  do j=1,^NC
     do k=1,^NC
        if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
           cross_D = cross_D + &
                sqrtgamma*lvc(idir,j,k)*a(j)*b(k)
     end do
  end do

end function cross_D
!=============================================================================
subroutine set_particles_dt_grlorentz()

use mod_metric
use mod_particles
include 'amrvacdef.f'

integer                                  :: ipart, iipart, i, j
double precision                         :: t_min_mype, tout, dt_particles_mype, dt_cfl
double precision                         :: v(1:^NC)
double precision                         :: alpha_p, u0, uu
double precision, dimension(1:^NC)       :: betaup_p, tmp
double precision, dimension(1:^NC,1:^NC) :: gammainv_p
double precision                         :: cfl=.1d0, dt0=.01d0
!-----------------------------------------------------------------------------
t_min_mype        = tsync_particles
dt_particles_mype = bigdouble

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

! Reinitialise dt in case it becomes too small due to matching with tout
 particle(ipart)%self%dt = dt0

! Distinguish particle/photon
if (particle(ipart)%self%m==0.d0) then ! photon
  uu = 0.d0
else ! particle
  uu = 1.d0
end if

! Get metric functions at position x^i
call get_alpha({particle(ipart)%self%x(^D)},alpha_p)

do i=1,^NC
  call get_beta(i,{particle(ipart)%self%x(^D)},betaup_p(i))
end do

do i = 1,ndir
  do j = 1,ndir
    call get_gammainv_component_analytic(i,j,{particle(ipart)%self%x(^D)},gammainv_p(i,j))
  end do
end do

! Compute u0
call get_u0(particle(ipart)%self%u,uu,gammainv_p,alpha_p,u0)

! make sure we step only one cell at a time:
tmp = matmul(gammainv_p,particle(ipart)%self%u)
{v(^C)   = abs(tmp(^C)-betaup_p(^C)*u0)/u0\}

{#IFNDEF D1
   dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)})
}
{#IFDEF D1
   dt_cfl = rnode(rpdx1_,particle(ipart)%igrid)/v(1)
}

   particle(ipart)%self%dt = min(dt_cfl*cfl,particle(ipart)%self%dt)
   particle(ipart)%self%dt = limit_dt(particle(ipart)%self%t,particle(ipart)%self%dt)

   dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)
   t_min_mype = min(t_min_mype,particle(ipart)%self%t)

end do !ipart

call save_tmin(t_min_mype,dt_particles_mype)

end subroutine set_particles_dt_grlorentz
!=============================================================================
subroutine get_u0(u,uu,gammainv_p,alpha_p,u0)

double precision, dimension(1:^NC), intent(in)         :: u
double precision, intent(in)                           :: alpha_p,uu
double precision, dimension(1:^NC,1:^NC), intent(in)   :: gammainv_p
double precision, intent(out)                          :: u0
double precision, dimension(1:^NC)                     :: tmp
!-----------------------------------------------------------------------------

tmp = matmul(gammainv_p,u)
u0 = sqrt(sum(u*tmp) + uu) / alpha_p
!=============================================================================
end subroutine get_u0
!=============================================================================
subroutine get_e(igrid,x,tloc,e)

! Get the electric field in the grid at postion x.

use mod_particles
use mod_gridvars
include 'amrvacdef.f'

integer,intent(in)                                 :: igrid
double precision,dimension(^NC), intent(in)        :: x
double precision, intent(in)                       :: tloc
double precision,dimension(^NC), intent(out)       :: e
double precision,dimension(^NC)                    :: e1, e2
double precision,dimension(^NC)                    :: u, u1, u2, beta, b
double precision                                   :: lfac
integer                                            :: ic^D
double precision                                   :: td
!-----------------------------------------------------------------------------

if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e2(^C))\}

td = (tloc - t) / dt

{^C& e(^C) = e1(^C) * (1.0d0 - td) + e2(^C) * td\}

end if

end subroutine get_e
!=============================================================================
subroutine get_b(igrid,x,tloc,b)

use mod_particles
use mod_gridvars
include 'amrvacdef.f'

integer,intent(in)                                 :: igrid
double precision,dimension(^NC), intent(in)        :: x
double precision, intent(in)                       :: tloc
double precision,dimension(^NC), intent(out)       :: b
integer                                            :: ic^D
double precision,dimension(^NC)                    :: b1, b2
double precision                                   :: td
!-----------------------------------------------------------------------------

if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b2(^C))\}

td = (tloc - t) / dt

{^C& b(^C) = b1(^C) * (1.0d0 - td) + b2(^C) * td\}

end if

end subroutine get_b
!=============================================================================
}
}
