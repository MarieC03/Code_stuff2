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
subroutine advance(iit)

{#IFDEF PARTICLES
use mod_particles, only: tmax_particles
use mod_timing, only: tpartc, tpartc0
}
{#IFDEF DY_SP
use mod_cfc
use mod_cfc_parameters
use mod_gw_br
}
{#IFDEF M1
use mod_m1
use mod_m1_eas
}
include 'amrvacdef.f'

integer, intent(in) :: iit

integer :: iigrid, igrid, idimsplit, ix^D, i
{#IFDEF TCRKL2
integer :: s
}
logical :: firstsweep, lastsweep
common/first/firstsweep,lastsweep
double precision :: time_in, time_metric_in, time_metric_end, time_end
double precision :: max_local(10), max_global(10)
double precision :: min_local(1), min_global(1)
double precision :: mass_total_local, mass_total

!-----------------------------------------------------------------------------
! when computing to steady state, store old solution
if(residmin>smalldouble)then
!$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call copy_state(ps(igrid),psres(igrid))
  end do
!$OMP END PARALLEL DO
endif

{#IFDEF RAY
call update_rays
}

! split source addition
if(ssplitdivb .or. ssplitresis .or. ssplituser) &
  call addsource_all(.true.)

{#IFDEF TCRKL2
! split thermal conduction source addition 1/2
if(conduction) then
  if(dt/dtimpl < 0.5d0) then
    s=1
    else if(dt/dtimpl < 2.d0) then
    s=2
    else
    s=ceiling((dsqrt(9.d0+8.d0*dt/dtimpl)-1.d0)/2.d0)
    ! only use odd s number
    s=s/2*2+1
  endif
  if(mype==0 .and. .false.) print*,'supertime steps:', s, &
              ' normal subcycles:',ceiling(dt/dtimpl/2.d0), &
              'time step ratio:',dt/dtimpl
  call thermalconduct_RKL2(s,dt/2.d0,t) 
endif
}

if (use_t_dependent_output) then
!if (use_t_dependent_output .or. gw_br_use_I3_ij) then
  ! old solution values at t_n-1 no longer needed: make copy of w(t_n)
  !$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call copy_state(ps(igrid),psold(igrid))
  end do
  !$OMP END PARALLEL DO
endif

! split implicit source addition
if (sourceimpl) then
   if(mype==0.and..false.) print *,'implicit source addition'
   call addsource_impl
endif

firstsweep=.true.
if (dimsplit) then
   if ((iit/2)*2==iit .or. typedimsplit=='xy') then
      ! do the sweeps in order of increasing idim,
      do idimsplit=1,ndim
         lastsweep= idimsplit==ndim
         call advect(idimsplit,idimsplit)
      end do
   else
      ! If the parity of "iit" is odd and typedimsplit=xyyx,
      ! do sweeps backwards
      do idimsplit=ndim,1,-1
         lastsweep= idimsplit==1
         call advect(idimsplit,idimsplit)
      end do
   end if
else
   ! Add fluxes from all directions at once
   lastsweep= .true.
   call advect(1,ndim)
end if

! split source addition
if(ssplitdivb .or. ssplitresis .or. ssplituser) &
  call addsource_all(.false.)

{#IFDEF TCRKL2
! split thermal conduction source addition 2/2
if(conduction) then
  call thermalconduct_RKL2(s,dt/2.d0,t+dt) 
endif
}

{#IFDEF PARTICLES
tpartc0 = MPI_WTIME()
call handle_particles()
tpartc = tpartc + (MPI_WTIME() - tpartc0)
}

if(residmin>smalldouble)then
!$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     pwres(igrid)%w(ixG^T,1:nwflux)= &
        pw(igrid)%w(ixG^T,1:nwflux)-pwres(igrid)%w(ixG^T,1:nwflux)
     {#IFDEF STAGGERED
     pwsres(igrid)%w = pws(igrid)%w-pwsres(igrid)%w
     }
  end do
!$OMP END PARALLEL DO
endif

!time_metric_in = MPI_WTIME()
{#IFDEF DY_SP
! check if it's time to update metric
if (cfc_update_flag()) then
   ! primitive will be updated when solving cfc
   call cfc_solve_metric()
else
   ! transform the variables to primitive
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      !if(covariant)myM => mygeo%m
      !call set_tmpGlobals(igrid)
      call primitive(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
   end do
   !$OMP END PARALLEL DO
end if
!Extra c2p may include inside gw_br
if (use_gw_br) then
   if (gw_br_update_flag()) then
      call solve_gw_br()
   endif
endif
}
{#IFNDEF DY_SP
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call set_tmpGlobals(igrid)
      call primitive(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
   end do
   !$OMP END PARALLEL DO
}

!time_metric_end = MPI_WTIME()

! after everything, update cons for next timestep (especially with new metric) and then get ghostcells as well
call getbc(t,ps,psCoarse)

end subroutine advance
!=============================================================================
subroutine advect(idim^LIM)
 
!  integrate all grids by one step of its delta(t)

! This subroutine is in VAC terminology equivalent to
! `advect' (with the difference that it will `advect' all grids)

include 'amrvacdef.f'

integer, intent(in) :: idim^LIM

integer :: iigrid, igrid
logical :: firstsweep, lastsweep
common/first/firstsweep,lastsweep

double precision, parameter :: gamma_imex122 = 1.0d0 
double precision, parameter :: gamma_imex222 = 1.0d0 - 1.0d0/sqrt(2.0d0)
{#IFDEF UNITS_TESTS
double precision :: TESTqtC = 1.0d+10
}
!-----------------------------------------------------------------------------

! copy w instead of wold because of potential use of dimsplit or sourcesplit
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call copy_state(ps(igrid),ps1(igrid))
end do
!$OMP END PARALLEL DO

!-----------------------
! Preallocate psm1
{#IFDEF M1
  !$OMP PARALLEL DO PRIVATE(igrid)
  !do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
  do iigrid=1,igridstail; igrid=igrids(iigrid); ! with ghostzones?
    call alloc_m1_state(psm1(igrid),ixG^LL)
  end do
  !$OMP END PARALLEL DO

  {#IFDEF M1_TESTS
  ! this initialization happens here as psm1 is allocated here to save memory
  ! needed for testing m1 when microphysical rates are switched off
  !$OMP PARALLEL DO PRIVATE(igrid)
  !do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
  do iigrid=1,igridstail; igrid=igrids(iigrid);  
    call initonegrid_M1_usr(ixG^LL,ixM^LL,ps(igrid),psm1(igrid)) 
  end do
  !$OMP END PARALLEL DO
  }
}
!-----------------------------

istep=0

select case (typeadvance)
case ("onestep")
   call advect1(typefull1,one,    idim^LIM,t,          t, ps1,ps,psold,psm1)
case ("twostep")
   ! predictor step
   call advect1(typepred1,half,   idim^LIM,t,          t, ps,ps1,psold,psm1)
   ! corrector step
   call advect1(typefull1,one,    idim^LIM,t+half*dt,  t, ps1,ps,psold,psm1)
case ("ImEx12")
   ! ImEx method, e.g. Palenzuela et al. 2009 for srrmhd and Palenzuela et al,
   ! 2013
   ! and Palenzuela 2013 for grmhd.
   ! This is the first-second pseudo-ImEx scheme from Bucciantini and Del Zanna
   ! 2013

   ! predictor step
   call advect1(typepred1,half,   idim^LIM,t,          t, ps,ps1,psold,psm1)
   ! corrector step
   call advect1(typefull1,one,    idim^LIM,t+half*dt,  t, ps1,ps,psold,psm1)

case ("FullyImplicit")
   !This is a first order fully implicit (backward Euler) scheme
   !Simple and robust, uses only ONE implicit solve per timestep
   !
   ! The scheme is:
   ! U^{n+1} = U^{n} + dt * G(U^{n+1})
   !
   ! Where G contains the implicit collisional source terms.
   ! This is solved implicitly via Newton-Raphson rootfinding.
   !
   ! This scheme requires only ONE implicit solve and is unconditionally stable.
   ! It is first-order accurate in time.

   ! M1 EAS update (compute opacities and emissivities at current state)
   {#IFDEF M1
     {#IFNDEF M1_EXPLICIT
     {#IFDEF M1_RATES
      {#IFDEF M1_TAU_OPT 
        call m1_wradimpl_tau(dt, t, ps, psm1)
        psm1%calc_eqrates = .false.
      }
      psm1%calc_rates = .true.
      call m1_wradimpl_rates(dt, t, ps, psm1)
     }
     }
   }

   ! ============================================================
   ! SINGLE IMPLICIT UPDATE: Solve ps = ps_old + dt*G(ps)
   ! This updates ps in-place using backward Euler
   ! dtfactor = 1.0 means we use the full timestep dt
   ! ============================================================
   call global_implicit_update(one, dt, t, ps, ps, psm1)

   ! Apply boundary conditions at t^{n+1}
   call getbc(t+dt, ps, psCoarse)
case ("ImEx112")
   ! first m1_eas update
   {#IFDEF M1
     {#IFNDEF M1_EXPLICIT
     {#IFDEF M1_RATES
       {#IFDEF M1_TAU_OPT 
         call m1_wradimpl_tau(dt,t,ps,psm1)
         psm1%calc_eqrates = .false. !witout tau_opt is true
         call m1_wradimpl_rates(dt,t,ps,psm1)
       }
       {#IFNDEF M1_TAU_OPT
         psm1%calc_rates = .true.
         call m1_wradimpl_rates(dt,t,ps,psm1)
       }
     }
     }
   }
   ! ps1 = U^n (initial state before implicit solve)
   ! Implicit step: ps1 = U^n + γ Δt G(ps1)
   call global_implicit_update(gamma_imex122, dt, t, ps1, ps, psm1)

   ! Overwrite ps with explicit step: ps = ps + dt F(ps1)
   call advect1(typefull1, one, idim^LIM, t, t, ps1, ps, psold, psm1)

   ! Final averaging: ps = 0.5 * (ps + ps1)
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid = 1, igridstail; igrid = igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux) = 0.5d0 * (pw(igrid)%w(ixG^T,1:nwflux) + pw1(igrid)%w(ixG^T,1:nwflux))
      {#IFDEF STAGGERED
      pws(igrid)%w(ixG^T,1:nwflux) = 0.5d0 * (pws(igrid)%w(ixG^T,1:nwflux) + pws1(igrid)%w(ixG^T,1:nwflux))
      }
   end do
   !$OMP END PARALLEL DO
   ! Apply boundary conditions for the new state
   call getbc(t + dt, ps, psCoarse)


case ("ImEx122")
   ! IMEX(1,2,2): U^(1) = U^n + γ Δt G(U^(1))
   !             U^*  = U^n + Δt F(U^(1))         (explicit step -> overwrite ps)
   !             U^(n+1) = 0.5*(U^(1) + U^*)      (final average)
   ! explicit tableau
   !  0  |  0     0
   !  1  |  1     0
   ! ----|------------
   !     | 1/2   1/2
   ! implicit tablea
   ! γ  |  γ     0
   ! ----|------------
   !     | 1.0    0.0

   ! Preallocate ps2 as y^n for the implicit update (xi storage)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)
      ! store y^n into ps2 (including staggered)
      pw2(igrid)%w  = pw(igrid)%w
      {#IFDEF STAGGERED
      pws2(igrid)%w = pws(igrid)%w
      }
   end do

   ! first m1_eas update
   {#IFDEF M1
     {#IFNDEF M1_EXPLICIT
     {#IFDEF M1_RATES
       {#IFDEF M1_TAU_OPT 
         call m1_wradimpl_tau(dt,t,ps,psm1)
         psm1%calc_eqrates = .false. !witout tau_opt is true
         call m1_wradimpl_rates(dt,t,ps,psm1)
       }
       {#IFNDEF M1_TAU_OPT
         psm1%calc_rates = .true.
         call m1_wradimpl_rates(dt,t,ps,psm1)
       }
     }
     }
   }
   ! (1) Implicit solve: xi = y^n + gamma*dt*G(xi)
   !     ps2 plays role of xi (updated in-place)
   call global_implicit_update(gamma_imex122, dt, t, ps2, ps, psm1)

   ! (2) Explicit update 1:
   !     ps1 = y^n + dt * F_ex(xi)
   !     Make sure ps1 starts from y^n (it does: ps->ps1 copy at top)
   ! ------------------------------------------------------------
   call advect1(typefull1, one, idim^LIM, t, t, ps2, ps1, psold, psm1)

   ! Enforce BC before using ps1 as qtC in the second explicit update
   call getbc(t+dt, ps1, psCoarse)

   ! (3) Explicit update 2:
   !     ps = y^n + dt * F_ex(ps1)
   !     IMPORTANT: ps must be y^n as base state here (still unchanged).
   ! ------------------------------------------------------------
   call advect1(typefull1, one, idim^LIM, t, t, ps1, ps, psold, psm1)

   ! (4) Final combination:
   !     ps = 0.5*(xi + yE2) = 0.5*(ps2 + ps)
   ! ------------------------------------------------------------
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux) = half*( pw(igrid)%w(ixG^T,1:nwflux) + pw2(igrid)%w(ixG^T,1:nwflux) )
      {#IFDEF STAGGERED
      pws(igrid)%w = half*( pws(igrid)%w + pws2(igrid)%w )
      }
   end do
   !$OMP END PARALLEL DO

   ! Final BC at t^{n+1}
   call getbc(t+dt, ps, psCoarse)

 
case ("ImEx222")
   !This is the second order ImEx scheme IMEX-SSP2 (2,2,2)
   !(Palenzuela et al. 2009)
   !(Pareschi, Russo (https://arxiv.org/pdf/1009.2757.pdf))

   ! The method has two substeps, with the following 
   ! Butcher tableau for the explicit sector:
   ! 
   !    0 | 0   0
   !    1 | 1   0        := a_{ij}
   !    __|________
   !      | 1/2 1/2      := w_{j}
   !
   ! And the following for the implicit sector:
   !
   !    
   !    g | g      0
   !  1-g | 1-2g   g     := z_{ij}
   !    __|___________
   !      | 1/2    1/2   := k_{j}
   !
   ! Where g = 1 - 1/sqrt(2). 
   ! The convention here is as follows: 
   ! 
   ! U^{(0)} = U^{(n)}
   ! U^({i}) = U^{(0)} + \Delta t \sum_{j=1}^{i} a_{ij} X(U^{(j}))
   !                   + \Delta t \sum_{j=1}^{\nu} z_{ij} G(U^{(j}))  
   ! for i=1,..\nu 
   ! 
   ! U^{(n+1)} = U^{(0)} + \Delta t \sum_{j=1}^{\nu} w_j X(U^{(j}))
   !                     + \Delta t \sum_{j=1}^{\nu} k_{j} G(U^{(j}))
   ! 
   ! Where X contains the explicit source terms and G contains the 
   ! implicit ones. 
   ! 
   ! All in all, this scheme requires the storage of 1 explicit source 
   ! and 2 implicit ones.  

   !write(888,*)"------ Time ----",t
   !write(*,*)"------ Time ----",t

   ! Preallocate ps2 as y^n for the implicit update 
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)  
      {#IFDEF STAGGERED 
      !pws2(igrid)%w(ixG^T,1:nwflux) = pws(igrid)%w(ixG^T,1:nwflux)
      pws2(igrid)%w= pws(igrid)%w
      }    
   end do

    ! first m1_eas update
   {#IFDEF M1
     {#IFNDEF M1_EXPLICIT
     {#IFDEF M1_RATES
        {#IFDEF M1_TAU_OPT 
          call m1_wradimpl_tau(dt,t,ps,psm1)
          psm1%calc_eqrates = .false. !witout tau_opt is true
          call m1_wradimpl_rates(dt,t,ps,psm1)
        }
        {#IFNDEF M1_TAU_OPT
          psm1%calc_rates = .true.
          call m1_wradimpl_rates(dt,t,ps,psm1)
        }
     }
     }
   }
   ! Solve ps2 = ps + gamma_imex222*dt*F_im(ps2) 
   ! Solve xi1 = y^n + gamma.dt.F_im(xi1)
   call global_implicit_update(gamma_imex222, dt, t, ps2, ps,psm1)

   ! Set ps1 = ps1 + dt.F_ex(ps2)
   ! i.e Set ps1 = y^n + dt.F_ex(xi1)
   ! ... qtC,sCT,qt,s,sold,sm1)
   call advect1(typefull1, one, idim^LIM, t, t, ps2, ps1,psold,psm1)

   ! Set ps2 = dt.F_im(xi1)        (is at t^n)
   ! i.e Set ps2 = (ps2-ps)/gamma_imex222 
   ! Set ps  = y^n + dt/2 . F(xi1) (is at t^n+dt/2)
   ! Set ps1 = y^n + dt.F_ex(xi1) + (1-2.gamma).dt.F_im(xi1) and enforce BC (at t^n+dt)
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixG^T,1:nwflux) = (pw2(igrid)%w(ixG^T,1:nwflux) - pw(igrid)%w(ixG^T,1:nwflux)) / gamma_imex222
      {#IFDEF STAGGERED 
      pws2(igrid)%w = (pws2(igrid)%w - pws(igrid)%w) / gamma_imex222
      }
   end do
   !$OMP END PARALLEL DO

   ! ps2 = dt.F_im(xi1)
   ! ps1 = y^n + dt.F_ex(xi1)
   ! ps  = y^n
   ! Set ps  = y^n + dt/2 . F(xi1) (is at t^n+dt/2)
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux) = half*(pw(igrid)%w(ixG^T,1:nwflux) + pw1(igrid)%w(ixG^T,1:nwflux) + pw2(igrid)%w(ixG^T,1:nwflux))
      {#IFDEF STAGGERED 
      pws(igrid)%w = half*(pws(igrid)%w + pws1(igrid)%w + pws2(igrid)%w)
      }
   end do
   !$OMP END PARALLEL DO

   ! Set ps1 = y^n + dt.F_ex(xi1) + (1-2.gamma).dt.F_im(xi1) and enforce BC (at t^n+dt)
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux) = pw1(igrid)%w(ixG^T,1:nwflux) + (1.0d0 - 2.0d0*gamma_imex222)*pw2(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED 
      pws1(igrid)%w = pws1(igrid)%w + (1.0d0 - 2.0d0*gamma_imex222)*pws2(igrid)%w
      }
   end do
   !$OMP END PARALLEL DO

   !M1_TODO: comment this getbc
   !call getbc(t+dt,ps1,psCoarse)
   call getbc(t+dt,ps1,psCoarse) !KEN this was unnecessary

   ! ps2 = dt.F_im(xi1)
   ! ps1 = y^n + dt.F_ex(xi1) + (1-2.gamma).dt.F_im(xi1) 
   ! Preallocate ps2 as xi1 for the implicit update (is at t^n)
   ! xi1 = y^n + gamma.dt.F_im(xi1) <- This what we want
   ! (pws2 = 2ps -ps1 -gamma ps2)
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      !pw2(igrid)%w(ixG^T,1:nwflux) = 2.0d0*pw2(igrid)%w(ixG^T,1:nwflux) - pw1(igrid)%w(ixG^T,1:nwflux) - gamma_imex222*pw2(igrid)%w(ixG^T,1:nwflux)
      !{#IFDEF STAGGERED 
      !pws2(igrid)%w(ixG^T,1:nwflux) = 2.0d0*pws2(igrid)%w(ixG^T,1:nwflux) - pws1(igrid)%w(ixG^T,1:nwflux) - gamma_imex222*pws2(igrid)%w(ixG^T,1:nwflux)
      !}
      pw2(igrid)%w(ixG^T,1:nwflux) = 2.0d0*pw(igrid)%w(ixG^T,1:nwflux) - pw1(igrid)%w(ixG^T,1:nwflux) - gamma_imex222*pw2(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED 
      pws2(igrid)%w = 2.0d0*pws(igrid)%w - pws1(igrid)%w - gamma_imex222*pws2(igrid)%w
      }
   end do
   !$OMP END PARALLEL DO

   ! second m1_eas update
   {#IFDEF M1
     {#IFNDEF M1_EXPLICIT
     {#IFDEF M1_RATES
      !true if rates should be updated every substep
      if(m1_2_eas_updates)then
         {#IFDEF M1_TAU_OPT 
           call m1_wradimpl_tau(dt,t,ps1,psm1)
           psm1%calc_eqrates = .false. !without tau_opt is true
           call m1_wradimpl_rates(dt,t,ps1,psm1) 
         }
         {#IFNDEF M1_TAU_OPT
           psm1%calc_rates = .true.
           call m1_wradimpl_rates(dt,t,ps1,psm1) 
         }
      end if 
     }
     }
   }

   ! Solve ps2 = (ps1) + gamma.dt.F_im(ps2) 
   ! i.e Solve xi2 = (ps1) + gamma.dt.F_im(xi2)
   call global_implicit_update(gamma_imex222, dt, t, ps2, ps1,psm1)

   ! Add dt/2.F_im(xi2) to ps
   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux) = pw(igrid)%w(ixG^T,1:nwflux) + (pw2(igrid)%w(ixG^T,1:nwflux) - pw1(igrid)%w(ixG^T,1:nwflux)) / (2.0d0 * gamma_imex222)
      {#IFDEF STAGGERED 
      pws(igrid)%w = pws(igrid)%w + (pws2(igrid)%w - pws1(igrid)%w) / (2.0d0 * gamma_imex222)
      }
   end do
   !$OMP END PARALLEL DO
   
   ! ps  = y^n + dt/2 . (F(xi1)+F_im(xi2)) (is at t^n+dt/2)
   ! ps2 = xi2 = dt.F_im(xi1)
   ! ps2 = (ps1) + gamma.dt.F_im(ps2)
   ! i.e Set ps = y^n + dt/2.(F(xi1)+F(xi2)) = y^(n+1)
   call advect1(typefull1, half, idim^LIM, t+dt, t+half*dt, ps2, ps,psold,psm1)

case("ImEx3")
   !This is the third order ImEx scheme IMEX-SSP3 (4,3,3)
   !(Palenzuela et al. 2009)

   call mpistop('This ImEx scheme is not implemented yet')
case ("threestep")
  
   ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
   call advect1(typefull1,one, idim^LIM,t,             t, ps,ps1,psold,psm1)

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)      
   end do
   
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixG^T,1:nwflux)=0.75d0*pw(igrid)%w(ixG^T,1:nwflux)+0.25d0*&
        pw1(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
      {#IFDEF STAGGERED
      pws2(igrid)%w = 0.75d0*pws(igrid)%w + 0.25d0*pws1(igrid)%w
      }
   end do
!$OMP END PARALLEL DO

   call advect1(typefull1,0.25d0, idim^LIM,t+dt, t+dt*0.25d0, ps1,ps2,psold,psm1)

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)=1.0d0/3.0d0*pw(igrid)%w(ixG^T,1:nwflux)+&
        2.0d0/3.0d0*pw2(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws(igrid)%w=1.0d0/3.0d0*pws(igrid)%w + 2.0d0/3.0d0*pws2(igrid)%w
      }
   end do   
!$OMP END PARALLEL DO
   
   call advect1(typefull1,2.0d0/3.0d0, idim^LIM,t+dt/2.0d0, t+dt/3.0d0,  ps2,ps,psold,psm1)

case ("ssprk43")
! Strong stability preserving 4 stage RK 3rd order method by Ruuth and Spiteri
!
! Ruuth & Spiteri
! J. S C, 17 (2002) p. 211 - 220
!
! supposed to be stable up to CFL=2.
! don't use time-dependent sources since i did not bother to set those intermediate times.
! oliver.

! === First step ===
   call advect1(typefull1,0.5d0, idim^LIM,t,           t, ps,ps1,psold,psm1)

! === Second step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixG^T,1:nwflux)=pw1(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixG^T,nwflux+1:nw) = &
           pw(igrid)%w(ixG^T,nwflux+1:nw)
      {#IFDEF STAGGERED
      pws2(igrid)%w = pws1(igrid)%w
      }
   end do
   !$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0, idim^LIM,t,           t+dt, ps1,ps2,psold,psm1)

! === Third step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps3(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw3(igrid)%w(ixG^T,1:nwflux)=2.0d0/3.0d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 1.0d0/3.0d0 * pw2(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw3(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
      {#IFDEF STAGGERED
      pws3(igrid)%w=2.0d0/3.0d0 * pws(igrid)%w + 1.0d0/3.0d0 * pws2(igrid)%w
      }
   end do
   !$OMP END PARALLEL DO
   
   call advect1(typefull1,1.0d0/6.0d0, idim^LIM,t,       t+dt, ps2,ps3,psold,psm1)

! === Fourth step ===
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)=pw3(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws(igrid)%w=pws3(igrid)%w
      }
   end do   
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0, idim^LIM, t,            t+dt, ps3,ps,psold,psm1)



case ("ssprk54")
! Strong stability preserving 5 stage RK 4th order method by Ruuth and Spiteri
!
! SIAM J. NUMER. ANAL.
! Vol. 40, No. 2, pp. 469–491
! c 2002 Society for Industrial and Applied Mathematics
! A NEW CLASS OF OPTIMAL
! HIGH-ORDER STRONG-STABILITY-PRESERVING
! TIME DISCRETIZATION METHODS
!
! E.g. Table A.2
!
! I have the actual coefficients however from the overview article by Gottlieb, JoSC 25 (2005)
! ON HIGH ORDER STRONG STABILITY PRESERVING RUNGE-KUTTA AND MULTI STEP TIME DISCRETIZATIONS
!
! there are slight differences in the coefficients (~8th digit after the .)
! This is SSP till CFL number 1.508 which makes the effective CFL per step ceff=0.377, 
! in contrast to the ceff=0.3333 of the classical RK3.  
!
! coded by oliver on 11/05/2013.  Enjoy!

! === First step ===
   call advect1(typefull1,0.391752226571890d0, idim^LIM, t ,t, ps,ps1,psold,psm1)

! === Second step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixG^T,1:nwflux)=0.444370493651235d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.555629506348765d0 * pw1(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw2(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
      {#IFDEF STAGGERED
      pws2(igrid)%w=0.444370493651235d0 * pws(igrid)%w + 0.555629506348765d0 * pws1(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.368410593050371d0, idim^LIM, t+0.391752226571890d0*dt, t+0.2176690962611688d0*dt, ps1,ps2,psold,psm1)

! === Third step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps3(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw3(igrid)%w(ixG^T,1:nwflux)=0.620101851488403d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.379898148511597d0 * pw2(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw3(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
      {#IFDEF STAGGERED
      pws3(igrid)%w=0.620101851488403d0 * pws(igrid)%w + 0.379898148511597d0 * pws2(igrid)%w
      }

   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.251891774271694d0, idim^LIM, t+0.5860796893115398d0*dt, t+0.222650588849706d0*dt, ps2,ps3,psold,psm1)

! === Fourth step ===
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps4(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw4(igrid)%w(ixG^T,1:nwflux)=0.178079954393132d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.821920045606868d0 * pw3(igrid)%w(ixG^T,1:nwflux)
      if (nw>nwflux) pw4(igrid)%w(ixG^T,nwflux+1:nw) = &
                       pw(igrid)%w(ixG^T,nwflux+1:nw)
      {#IFDEF STAGGERED
      pws4(igrid)%w=0.178079954393132d0 * pws(igrid)%w + 0.821920045606868d0 * pws3(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.544974750228521d0, idim^LIM, t+0.4745423631214d0*dt, t+0.390035880739132d0*dt, ps3,ps4,psold,psm1)

! Now recover back the dt*L(u3), store in pw1:
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux) = ( pw4(igrid)%w(ixG^T,1:nwflux) &
           - (0.178079954393132d0 * pw(igrid)%w(ixG^T,1:nwflux) &
           + 0.821920045606868d0 * pw3(igrid)%w(ixG^T,1:nwflux)) ) / 0.544974750228521d0
      {#IFDEF STAGGERED
      pws1(igrid)%w = ( pws4(igrid)%w &
           - (0.178079954393132d0 * pws(igrid)%w &
           + 0.821920045606868d0 * pws3(igrid)%w) ) / 0.544974750228521d0
      }
   end do
!$OMP END PARALLEL DO

! === Fifth step ===
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)= 0.517231671970585d0 * pw2(igrid)%w(ixG^T,1:nwflux) &
           + 0.096059710526147d0 * pw3(igrid)%w(ixG^T,1:nwflux) &
           + 0.063692468666290d0 * pw1(igrid)%w(ixG^T,1:nwflux) & 
           + 0.386708617503269d0 * pw4(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws(igrid)%w= 0.517231671970585d0 * pws2(igrid)%w &
           + 0.096059710526147d0 * pws3(igrid)%w &
           + 0.063692468666290d0 * pws1(igrid)%w &
           + 0.386708617503269d0 * pws4(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   call advect1(typefull1,0.226007483236906d0, idim^LIM, t+0.935010630967653d0*dt, t+0.710300048096804d0*dt, ps4,ps,psold,psm1)



case ("rk4")
   ! classical RK4 four step scheme
   call advect1(typefull1,0.5d0, idim^LIM,t,         t, ps,ps1,psold,psm1)
   
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws2(igrid)%w=pws(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0, idim^LIM,t+dt/2d0,   t, ps1,ps2,psold,psm1)


   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps3(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw3(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws3(igrid)%w=pws(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,one,   idim^LIM,t+dt/2d0,   t, ps2,ps3,psold,psm1)


!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux)=(pw1(igrid)%w(ixG^T,1:nwflux) &
                               +two*pw2(igrid)%w(ixG^T,1:nwflux) &
                                   +pw3(igrid)%w(ixG^T,1:nwflux) &
                              -4.0d0*pw(igrid)%w(ixG^T,1:nwflux))/3.0d0
      {#IFDEF STAGGERED
      pws1(igrid)%w = ( pws1(igrid)%w + two*pws2(igrid)%w &
           + pws3(igrid)%w -4.0d0*pws(igrid)%w )/3.0d0
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,1.0d0/6.0d0, idim^LIM,t+dt,  t, ps3,ps,psold,psm1)

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)=pw1(igrid)%w(ixG^T,1:nwflux)+&
        pw(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws(igrid)%w = pws1(igrid)%w + pws(igrid)%w
      }
   end do
!$OMP END PARALLEL DO


case ("fourstep")
   ! four step scheme, variant Hans De Sterck
   call advect1(typefull1,0.12d0, idim^LIM,t,                t, ps,ps1,psold,psm1)
   
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws2(igrid)%w = pws(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.25d0, idim^LIM,t+dt*0.12d0,      t, ps1, ps2,psold,psm1)
   
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws1(igrid)%w = pws(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0,  idim^LIM,t+dt/4d0   ,t, ps2,ps1,psold,psm1)
   call advect1(typefull1,one,    idim^LIM,t+dt/2d0   ,t, ps1,ps,psold,psm1)


case ("jameson")
   ! four step scheme, variant jameson
   call advect1(typefull1,0.25d0, idim^LIM,t,                  t, ps,ps1,psold,psm1)
   
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call alloc_state(ps2(igrid),ixG^LL)
   end do

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw2(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws2(igrid)%w = pws(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,(1.0d0/3.0d0), idim^LIM,t+dt*0.25d0, t, ps1,ps2,psold,psm1)

!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw1(igrid)%w(ixG^T,1:nwflux)=pw(igrid)%w(ixG^T,1:nwflux)
      {#IFDEF STAGGERED
      pws1(igrid)%w = pws(igrid)%w
      }
   end do
!$OMP END PARALLEL DO
   
   call advect1(typefull1,0.5d0,  idim^LIM,t+dt/3d0,            t, ps2,ps1,psold,psm1)
   call advect1(typefull1,one,    idim^LIM,t+dt/2d0,            t, ps1,ps,psold,psm1)


case default
   write(unitterm,*) "typeadvance=",typeadvance
   write(unitterm,*) "Error in advect: Unknown time integration method"
   call mpistop("Correct typeadvance")
end select

!-----------------------
! Deallocate psm1
{#IFDEF M1
{#IFNDEF M1_EXPLICIT
!!$OMP PARALLEL DO PRIVATE(igrid)
!do iigrid=1,igridstail; igrid=igrids(iigrid);
!   call dealloc_m1_state(psm1(igrid)) !,ixG^LL)
!end do
!!$OMP END PARALLEL DO
}
}
!-----------------------------

firstsweep=.false.
end subroutine advect
!=============================================================================
subroutine advect1(method,dtfactor,idim^LIM,qtC,qt,psa,psb,psc,psm1grid)

!  integrate all grids by one partial step

! This subroutine is equivalent to VAC's `advect1', but does
! the advection for all grids
include 'amrvacdef.f'

integer, intent(in)          :: idim^LIM
double precision, intent(in) :: dtfactor, qtC, qt
character(len=*), intent(in) :: method(nlevelshi)
type(state)                  :: psa(ngridshi), psb(ngridshi), psc(ngridshi)
type(m1_impl_state)          :: psm1grid(ngridshi) ! can't do IFDEF M1 since its argument
! .. local ..
double precision             :: qdt
{#IFDEF UNIT_TESTS
TESTqtC
}
integer                      :: iigrid, igrid, level
logical                      :: do_fluxfix
integer                      :: ix^D
!-----------------------------------------------------------------------------
istep=istep+1

do_fluxfix = ((time_advance.and.levmax>levmin) .and. (istep==nstep.or.nstep>2))

if (istep > 1) then  !last time already has the final getbc
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      call set_tmpGlobals(igrid)
      call primitive(ixG^LL,ixM^LL,psa(igrid)%w%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
   end do
   call getbc(qt,psa,psCoarse)
endif

if (do_fluxfix) call init_comm_fix_conserve(idim^LIM)

!$OMP PARALLEL DO PRIVATE(igrid,level,qdt)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   level=node(plevel_,igrid)
   qdt=dtfactor*dt_grid(igrid)

   call set_tmpGlobals(igrid)

   ! psa(igrid) is the actually evolving one, and then become psb
   call process1_grid(method(level),igrid,qdt,ixG^LL,idim^LIM,qtC,&
        psa(igrid),qt,psb(igrid),psc(igrid),psm1grid(igrid))
end do
!$OMP END PARALLEL DO

if (do_fluxfix) then
   do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      call sendflux(igrid,idim^LIM) ! OMP: not threadsafe!
   end do
   call fix_conserve(psb,idim^LIM)
   {#IFDEF STAGGERED
     call fix_edges(psb,idim^LIM)
     {#IFDEF DY_SP
     !$OMP PARALLEL DO PRIVATE(igrid)
     do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        call faces2centers(ixM^LL,psb(igrid))
     end do
     !$OMP END PARALLEL DO
     }
   }
end if

!! Follow Old bhac, updating faces2centers every substep again
!! FIXME: rememebr, if you update faces2centers even not fix_conserve --> lfac>lfac_max a lot happens
{#IFNDEF DY_SP
  {#IFDEF STAGGERED
  ! Now we fill the centers for the staggered variables
  !$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
     call faces2centers(ixG^LL,psb(igrid))
  end do
  !$OMP END PARALLEL DO
  }
}

! Fill ghost cells
qdt=dtfactor*dt
call getbc(qt+qdt,psb,psCoarse)

end subroutine advect1
!=============================================================================
subroutine process1_grid(method,igrid,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold,sm1)

! This subroutine is equivalent to VAC's `advect1' for one grid

include 'amrvacdef.f'

character(len=*), intent(in) :: method
integer, intent(in)          :: igrid, ixI^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt
type(state)                  :: sCT, s, sold
type(m1_impl_state)          :: sm1
! .. local ..
double precision             :: dx^D
double precision             :: fC(ixI^S,1:nwflux,1:ndim)
logical                      :: do_fluxfix
double precision             :: fE(ixI^S,1:^NC)
!-----------------------------------------------------------------------------
do_fluxfix = ((time_advance.and.levmax>levmin) .and. (istep==nstep.or.nstep>2))

dx^D=rnode(rpdx^D_,igrid);
call set_tmpGlobals(igrid)

call advect1_grid(method,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D, &
                  px(igrid)%x,sm1)


if (do_fluxfix) then
   call storeflux(igrid,fC,idim^LIM)
{#IFDEF STAGGERED
   call storeedge(igrid,ixI^L,fE,idim^LIM) }
end if

end subroutine process1_grid
!=============================================================================
subroutine advect1_grid(method,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D,x,sm1)

!  integrate one grid by one partial step 
include 'amrvacdef.f'

character(len=*), intent(in) :: method
integer, intent(in) :: ixI^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
type(state)      :: sCT, s, sold
type(m1_impl_state)  :: sm1
double precision :: fC(ixI^S,1:nwflux,1:ndim)
double precision, dimension(ixI^S,1:nw) :: wprim
double precision :: fE(ixI^S,1:^NC)
integer          :: ixO^L,ix^D

!integer            :: ixtile^D, ntiles^D ! number of tiles
!integer, parameter :: ntile1=32, ntile2=16, ntile3=4 ! number of cells in tile
!-----------------------------------------------------------------------------
ixO^L=ixI^L^LSUBdixB;
{#IFNDEF DY_SP
if (covariant) myM => mygeo%m
}
!call primitive(ixI^L,ixI^L,sCT%w%w,x) ! also applies floors to wCT
wprim = sCT%w%w
!call conserve(ixI^L,ixI^L,sCT%w%w,x,patchfalse)

! Now tile it:
!ntiles^D=(ixImax^D-ixImin^D-2*dixB+1)/ntile^D;

!{^D& do ixtile^DB = 1, ntiles^DB\}

!{^D& ixOmin^D = ixImin^D + dixB + (ixtile^D-1)*ntile^D \}
!{^D& ixOmax^D = ixOmin^D + ntile^D - 1 \}

select case (method)
case ('tvdlf','tvdlf1','hll','hll1')
   call finite_volume(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,sold,wprim,fC,fE,dx^D,x,sm1)
case ('fd_tvdlf','fd_hll')
   call finite_difference(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,sold,wprim,fC,fE,dx^D,x,sm1)
case ('source')
   call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),&
                                      ixI^L,ixO^L,1,nw,qtC,sCT%w%w,qt,s%w%w,x,.false.)
case ('nul')
   ! There is nothing to do
case ('tvdlfpos')
   call mpistop('if you want positivity preserving, &
                 choose tvdlf/hll and then PP_limiter to be true')
case default
   write(unitterm,*)'Error in advect1_grid:',method,' is unknown!'
   call mpistop("")
end select

!{end do\}

end subroutine advect1_grid
!=============================================================================
subroutine addsource_all(prior)

include 'amrvacdef.f'

logical, intent(in) :: prior

integer :: iigrid, igrid
double precision :: qdt, qt
!-----------------------------------------------------------------------------

if ((.not.prior).and.&
    (typesourcesplit=='sf' .or. typesourcesplit=='ssf')) return

if (prior) then
   qt=t
else
   qt=t+dt
end if

!$OMP PARALLEL DO PRIVATE(igrid,qdt)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
         qdt=dt_grid(igrid)
         if (B0field) then
            myB0_cell => pB0_cell(igrid)
            {^D&myB0_face^D => pB0_face^D(igrid)\}
         end if
         call addsource1_grid(igrid,qdt,qt,pw(igrid)%w)
      end do
!$OMP END PARALLEL DO
      
      call getbc(qt,ps,psCoarse)

end subroutine addsource_all
!=============================================================================
subroutine addsource1_grid(igrid,qdt,qt,w)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: qdt, qt
double precision, intent(inout) :: w(ixG^T,nw)

double precision :: w1(ixG^T,nw)
!-----------------------------------------------------------------------------

call set_tmpGlobals(igrid)

w1(ixG^T,1:nw)=w(ixG^T,1:nw)

select case (typesourcesplit)
case ('sf')
   call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('sfs')
   call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('ssf')
   call addsource2(qdt/2,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
   call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case ('ssfss')
   call addsource2(qdt/4,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
   call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
case default
   write(unitterm,*)'No such typesourcesplit=',typesourcesplit
   call mpistop("Error: Unknown typesourcesplit!")
end select

end subroutine addsource1_grid
!=============================================================================
subroutine addsource2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! Add source within ixO for iws: w=w+qdt*S[wCT]

include 'amrvacdef.f'

! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
logical, intent(in) :: qsourcesplit
!-----------------------------------------------------------------------------

! user defined sources, typically explicitly added
if(qsourcesplit .eqv. ssplituser) &
   call specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! physics defined sources, typically explicitly added,
! along with geometrical source additions
call addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

end subroutine addsource2
!=============================================================================
subroutine process(iit,qt)

include 'amrvacdef.f'

integer,intent(in)            :: iit
double precision, intent(in)  :: qt
! .. local ..
integer:: iigrid, igrid, level
!-----------------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(igrid,level)
do iigrid=1,igridstail; igrid=igrids(iigrid);

   call set_tmpGlobals(igrid)
   level = node(plevel_,igrid)
 
   call process_grid_usr(igrid,level,ixG^LL,ixM^LL,qt,pw(igrid)%w,px(igrid)%x)
   
end do
!$OMP END PARALLEL DO

end subroutine process
!=============================================================================
subroutine addsource_impl

include 'amrvacdef.f'

integer :: iigrid, igrid, icycle, ncycle
double precision :: qdt, qt, sumqdt
!-----------------------------------------------------------------------------

energyonly=.true.
if(sourceimplcycle)then
   ncycle=ceiling(dt/dtimpl)
   if(ncycle<1) then
     ncycle=1
     dtimpl=dt
   endif
else
   ncycle=1
endif

if(mype==0.and..false.) then
  print *,'implicit source addition will subcycle with ',ncycle,' subtimesteps'
  print *,'dt and dtimpl= ',dt,dtimpl,' versus ncycle*dtimpl=',ncycle*dtimpl
endif

qt=t
sumqdt=zero
do icycle=1,ncycle
  if(sourceparasts)then
    qdt=dtimpl/(one+parastsnu+(parastsnu-one)* &
         dcos(half*dpi*(two*dble(icycle)-one)/dble(ncycle)))
  else
    qdt=dt/dble(ncycle)
  endif
!$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call addsourceimpl1_grid(igrid,qdt,qt,pw(igrid)%w,px(igrid)%x)
  end do
!$OMP END PARALLEL DO

  sumqdt=sumqdt+qdt
  qt=qt+qdt
  call getbc(qt,ps,psCoarse)
enddo

if(mype==0.and..false.) then
  if(sourceparasts)then
     print *,'qt-t=',qt-t,'versus ncycle2*dtimpl=',ncycle*ncycle*dtimpl,&
       ' and sumqdt=',sumqdt
  else
     print *,'qt-t=',qt-t,'versus ncycle*dtimpl=',ncycle*dtimpl,&
       ' and sumqdt=',sumqdt
  endif
endif
energyonly=.false.

end subroutine addsource_impl
!=============================================================================
subroutine addsourceimpl1_grid(igrid,qdt,qt,w,x)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: qdt, qt, x(ixG^T,1:ndim)
double precision, intent(inout) :: w(ixG^T,nw)

double precision :: w1(ixG^T,nw)
!-----------------------------------------------------------------------------

call set_tmpGlobals(igrid)

w1(ixG^T,1:nwflux)=w(ixG^T,1:nwflux)

call specialsource_impl(qdt,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,x)

end subroutine addsourceimpl1_grid
!=============================================================================
!< Globally apply implicit update:
!<    psa = psb + dtfactor * qdt * G(psa)
!< Calls bcs for psa and implicit_update (to be provided)
!< inside the phys module.
subroutine global_implicit_update(dtfactor,qdt,qtC,psa,psb,psm1impl)
   
    include 'amrvacdef.f'

    type(state)          :: psa(ngridshi)   !< Compute implicit part from this state and update it
    type(state)          :: psb(ngridshi)   !< Will be unchanged, as on entry
    !type(state), target :: psa(ngridshi)   !< Compute implicit part from this state and update it
    !type(state), target :: psb(ngridshi)   !< Will be unchanged, as on entry
    type(m1_impl_state)  :: psm1impl(ngridshi)   !< Compute implicit part
    double precision, intent(in) :: qdt      !< overall time step dt
    double precision, intent(in) :: qtC      !< Both states psa and psb at this time level
    double precision, intent(in) :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)

    integer                        :: iigrid, igrid

    !> First copy all variables from a to b, this is necessary to account for
    ! quantities is w with no implicit sourceterm
    !do iigrid=1,igridstail; igrid=igrids(iigrid);
    !   psa(igrid)%w = psb(igrid)%w
    !end do
    !MC: psa and psb on input are of type state,
    ! to do psa%w (like in pw%w) the variable psa should 
    ! be of type : type(walloc), target
    ! copy states via funtion or as psa%w%w=psb%w%w for centered vars
   !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      call copy_state(psb(igrid),psa(igrid)) !psa out
    end do
   !$OMP END PARALLEL DO

    ! would be nice to have this as abstract interface 
    ! but in bhac all physics functions are just global
    ! routines in amrvacphys so we stick to the (very 
    ! bad) convention.
    ! The following function does the following: 
    ! psa = psb + dtfactor*qdt*G(psa)
    call implicit_update(dtfactor,qdt,qtC,psa,psb,psm1impl)

    !do iigrid=1,igridstail; igrid=igrids(iigrid);    
    !  call implicit_update(dtfactor,igrid,qdt,qtC,psa(igrid),psb(igrid),psm1impl(igrid),ixG^LL)
    !end do 

    ! enforce boundary conditions for psa
    ! call getbc(qtC,psa,psCoarse)

end subroutine global_implicit_update


  !********************************************************************
{#IFNDEF M1_EXPLICIT

{#IFDEF M1_TAU_OPT
subroutine m1_wradimpl_tau(qdt,qt,psb,psm1impl)
   {#IFDEF UNIT_TESTS
     use mod_m1_tests
     use amrvacdef
   }
   {#IFDEF M1
   use mod_m1_tau_optical
   }
   {#IFNDEF UNIT_TESTS
     include 'amrvacdef.f' 
   }
   type(state), dimension(ngridshi)          :: psb !< Left unchanged 
   type(m1_impl_state), dimension(ngridshi)  :: psm1impl
   double precision, intent(in) :: qdt        !< Timestep 
   double precision, intent(in) :: qt        !< Current time  
   double precision :: dx^D !> separation
   integer                      :: iigrid, igrid
  !"""""""""""""""""""""""""""""""""""
  {#IFNDEF M1_EXPLICIT
   ! unfortunately we need an extra c2p to get the fluid vars for the rates & tau

   do iigrid=1,igridstail; igrid=igrids(iigrid);
    call set_tmpGlobals(igrid)
    call primitive(ixG^LL,ixM^LL,psb(igrid)%w%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
  end do
  call getbc(qt,psb,psCoarse)
  !-----------------------------
    {#IFDEF M1
    !if(psm1impl%calc_rates .eqv. .true.) then
     ! grid loop 
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        dx^D=rnode(rpdx^D_,igrid);
       !> fill in psm1 with equilibrium rates an tau_path:
       call m1_add_tau(qdt,qt,dx^D,igrid,px(igrid)%x(ixG^T,1:ndim),psb(igrid)%w%w,psm1impl(igrid)%pwrad%wradimpl,ixG^LL)
      
     end do
    !end if 
    }
   }!"""""""""""""""""""""""""""""""""
   !TODO: boundary conditions for tau_path!!

end subroutine m1_wradimpl_tau
}
!=====================================================

  subroutine m1_wradimpl_rates(qdt,qt,psb,psm1impl)
     {#IFDEF UNIT_TESTS
       use mod_m1_tests
       use amrvacdef
     }
     {#IFDEF M1
     use mod_m1_eas, only: m1_add_eas
     }
     {#IFNDEF UNIT_TESTS
       include 'amrvacdef.f' 
     }
     type(state), dimension(ngridshi)          :: psb !< Left unchanged 
     type(m1_impl_state), dimension(ngridshi)  :: psm1impl
     double precision, intent(in) :: qdt        !< Timestep 
     double precision, intent(in) :: qt        !< Current time  
     integer                      :: iigrid, igrid
     logical :: calc_equilib_rates
    !"""""""""""""""""""""""""""""""""""
    {#IFNDEF M1_EXPLICIT
     ! unfortunately we need an extra c2p to get the fluid vars for the rates

     do iigrid=1,igridstail; igrid=igrids(iigrid);
      call set_tmpGlobals(igrid)
      call primitive(ixG^LL,ixM^LL,psb(igrid)%w%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
    end do
    call getbc(qt,psb,psCoarse) 
    ! We probably need, since the c2p
    !-----------------------------
      {#IFDEF M1
      calc_equilib_rates = psm1impl(1)%calc_eqrates !psm1impl(igrid)%calc_eqrates
      !if(psm1impl%calc_rates .eqv. .true.) then
       ! grid loop 
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         !> fill in psm1 with rates:
         call m1_add_eas(qdt,qt,igrid,px(igrid)%x(ixG^T,1:ndim),psb(igrid)%w%w,psm1impl(igrid)%pwrad%wradimpl,ixG^LL,calc_equilib_rates)
        
       end do
      !end if 
      }
     }!"""""""""""""""""""""""""""""""""

  end subroutine m1_wradimpl_rates
}


