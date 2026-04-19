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
program amrvac
use mod_forest
use mod_multigrid_coupling
use mod_cfc_parameters
use mod_coord_cfc
{^IFTHREED
   use mod_healpix_det, only: initialize_healpix_detectors
}
!use mod_small_values

! AMRVAC solves a set of hyperbolic equations:
!              u  +  Div [f(u)]    =  s
!               t       x    
!  using adaptive mesh refinement.

! following line may avoid mystery problems in the combination
!  where we use ifort compiler and MPT (SGIs MPI implementation)
!DEC$ ATTRIBUTES NOINLINE :: read_snapshot

{#IFDEF PARTICLES
!use mod_gridvars, only: init_gridvars, finish_gridvars
}
{#IFDEF M1
use mod_m1, only: m1_startup
}
use mod_variables
include 'amrvacdef.f'
integer :: iigrid, igrid, total_cell, total_cell_local
integer :: itin
double precision :: time0, time_in
double precision :: max_local(10), max_global(10), u_max_local(1:ndir), u_max(1:ndir)
double precision :: min_local(1), min_global(1)
integer :: n_coarsen, n_refine, ix^D, i
{#IFDEF STAGGERED
integer                                       :: ixGs^L
}

!-----------------------------------------------------------------------------
call comm_start

time_advance=.false.

time0=MPI_WTIME()
time_bc=zero

! changed
call setup_variables

! get input parameters
call readcommandline
call readparameters
call initialize_vars
call init_comm_types

{#IFDEF DY_SP
   select case (coordinate)
   case(cartesian)
     call switch_coord_to_cart
   case(spherical)
     call switch_coord_to_sp
   case(cylindrical)
     call switch_coord_to_cylind
   case default
     call mpistop('no this coordinate for cfc framework')
   end select
}

{#IFDEF M1
call m1_startup
}

! Begin of the code
! -----------------
if (snapshotini/=-1) then
   ! restart from previous file or dat file conversion
   ! get input data from previous VAC/AMRVAC run
   itin=it
   tin=t

   ! order: first physics dependent part then user for defaults.
   call initglobal

{#IFDEF RAY
   call init_rays
}
! Here will reallocate grids for reading snapshot
{#IFDEF HDF5
   if (hdf5_ini) then
     call read_snapshot_hdf5
   else
     select case (typeparIO)
     case (1,-2)
        call read_snapshot
     case (0)
        call read_snapshot_nopar
     case(-1)
        call read_snapshot_noparf
     end select
   end if
\}

{#IFNDEF HDF5
   select case (typeparIO)
   case (1,-2)
      call read_snapshot
   case (0)
      call read_snapshot_nopar
   case(-1)
      call read_snapshot_noparf
   end select

\}

   ! modify globals
   if (changeglobals) call initglobal

{#IFDEF BOUNDARYDRIVER
   call read_boundary
}

{#IFDEF PARTICLES
   call init_tracerparticles
   call getbc(t,ps,psCoarse)
   call init_gridvars
   call read_particles_snapshot
   call finish_gridvars
}

   if (itreset) it=itin
   if (treset) t=tin
   ! modify initial condition
   if (firstprocess) then 
     call modify_IC
     call getbc(t,ps,psCoarse)
     {#IFDEF DY_SP
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         call set_tmpGlobals(igrid)
         call get_B_field(ixG^LL,ixM^LL,ps(igrid))
       enddo
       call getbc(t,ps,psCoarse)
       call fixer_after_usr2()
       call set_conserve()
     }
   endif

   ! Select active grids
   call selectgrids

   call getbc(t,ps,psCoarse)

   ! reset AMR grid
   if (resetgrid) call settree

   call getbc(t,ps,psCoarse)

   ! set up boundary flux conservation arrays
   if (levmax>levmin) call allocateBflux
{#IFDEF DY_SP
   if (use_multigrid) call mg_setup_multigrid()
   if (restart_init_metric) then
      write(11,*)"in restart init"
     call metric_initialize() 
     call getbc(t,ps,psCoarse)
     call set_conserve()
     call getbc(t,ps,psCoarse)
   endif
}

   if (convert) then
       {#IFDEF HDF5
       if (convert_type == 'hdf5') then
           call generate_plotfile
{#IFDEF PARTICLES
           call finish_tracerparticles
}
           call comm_finalize
           stop
       end if
       }
       if (npe/=1.and.(.not.(index(convert_type,'mpi')>=1)) &
            .and. convert_type .ne. 'user')  &
       call mpistop("non-mpi conversion only uses 1 cpu")
       call generate_plotfile
{#IFDEF PARTICLES
       call finish_tracerparticles
}
       call comm_finalize
       call mpistop('mpi stop after comm_finalize durign restart')
   end if
else  ! not from restart -----------------------------------------------------

   ! order: first physics dependent part then user for defaults.
   call initglobal
{#IFDEF RAY
   call init_rays
}
   ! form and initialize all grids at level one
   call initlevelone
   ! set up and initialize finer level grids, if needed
   call settree
   ! set up boundary flux conservation arrays
   if (levmax>levmin) call allocateBflux
{#IFDEF DY_SP

   if (use_multigrid) call mg_setup_multigrid()

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call set_tmpGlobals(igrid)
      call initonegrid_usr2(ixG^LL, ixM^LL, ps(igrid))  !call p2c
   enddo

   call getbc(t,ps,psCoarse)

   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call set_tmpGlobals(igrid)
     call get_B_field(ixG^LL,ixM^LL,ps(igrid))
   enddo

   call metric_initialize() 
   ! Set conserve variables after getting actual metric vars
   call set_conserve()
   call getbc(t,ps,psCoarse)

   {#IFDEF STAGGERED
   select case(clean_init_divB)
   case ('vecpot')
      ! Now that the grid has been set and the flux-conservation arrays are allocated,
      ! re-calculate the magnetic field from the vector potential in a completely
      ! divergence free way.
      if (ndim < 3) call mpistop('It is not 3D simulation, not working vecpot&
                                initial cleaning for divB')

!fixme: not sure the Bs still has psi6 inside?
      if (levmax>levmin .and. ndim .eq. 3) call recalculateB
   {^NOONED
   case ('mg')
      ! Project out the divB using poisson solver.
      ! Due to Jannis Teunissen. Thanks!   
      if (.not. use_multigrid) call mg_setup_multigrid()  ! "not" means you didnt turn it with CFC
      call clean_divb_multigrid()
   }
   case ('none')
   case default
      call mpistop('Unknown method selected for clean_init_divB!')
   end select
   }
   call getbc(t,ps,psCoarse)
   ! Warning: You should put normalize_initial B-field inside fixer_after_usr2, 
   ! If you did normalization, pls metric initialize again!

   !! TODO
   !! We need to do closure in initonegridusr2 to get Gamma to make nrad cons !!
   !do iigrid=1,igridstail; igrid=igrids(iigrid);
   !   call set_tmpGlobals(igrid)
   !   ! initonegrid_usr3
   !   call initonegrid_usr2(ixG^LL, ixM^LL, ps(igrid))  !call p2c
   !enddo

   !KEN wrote this to fix the inconsistency of sqrtg after metric initialize
   {#IFDEF M1
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call set_tmpGlobals(igrid)
      call fix_M1_after_usr2(ixG^LL, ixM^LL, ps(igrid))
   enddo
   }

   call fixer_after_usr2()
   ! Set conserve variables after fixing divB and B^i after clean_init_divB
   ! not using set_conserve --> I want the diagonal ghostzones as well
   call set_conserve()
   do iigrid=1,igridstail; igrid=igrids(iigrid);
     call set_tmpGlobals(igrid)
     !call conserve(ixG^LL, ixG^LL, pw(igrid)%w, px(igrid)%x, patchfalse)
     ! one more c2p to ensure more consistent with new B^i
     call primitive(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x)
   enddo
   call getbc(t,ps,psCoarse)
!TODO: later you need to add normalize_B-field for DY_SP!!
}

{#IFNDEF DY_SP
   {#IFDEF STAGGERED
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      ! fixme: staggered grid cannot use usr2
      call initonegrid_usr2(ixG^LL, ixM^LL, ps(igrid))  
   enddo

   select case(clean_init_divB)
   case ('vecpot')
      ! Now that the grid has been set and the flux-conservation arrays are allocated,
      ! re-calculate the magnetic field from the vector potential in a completely
      ! divergence free way.
      ! IF NOT TURN THIS IN 2D FMtorus/BHaccretion, divB ~ 0.1 for FMtorus, never jetting 
      !call mpistop('Clean_init_divB, please choose none or mg')
      if (levmax>levmin .and. ndim .eq. 3) call recalculateB
   {^NOONED
   case ('mg')
      ! Project out the divB using poisson solver.
      ! Due to Jannis Teunissen. Thanks!   
      call mg_setup_multigrid()
      call clean_divb_multigrid()
   }
   case ('none')
   case default
      call mpistop('Unknown method selected for clean_init_divB!')
   end select
   }
  ! update the boundaries
  call getbc(t,ps,psCoarse)
{#IFDEF NORMINIT
   ! normalize the B-field strength
   call normalize_initial
}   
}


{#IFDEF PARTICLES
   call init_tracerparticles
   call getbc(t,ps,psCoarse)
   call init_gridvars
   call init_particles
   call finish_gridvars
}

   call selectgrids

end if !end of restart or not

{^IFTHREED
   if (healpix_det_open) call initialize_healpix_detectors()
}

{#IFDEF DY_SP
    cfc_tol(1)      = cfc_tol_evolve(1)
    cfc_tol(2)      = cfc_tol_evolve(2)
    cfc_tol(3)      = cfc_tol_evolve(3)
    gw_br_tol(1)    = gw_br_tol_evolve(1)
    gw_br_tol(2)    = gw_br_tol_evolve(2)
    gw_br_tol(3)    = gw_br_tol_evolve(3)
}

if (mype==0) then
   print*,'-----------------------------------------------------------------------------'
   write(*,'(a,f17.3,a)')' Startup phase took : ',MPI_WTIME()-time0,' sec'
   print*,'-----------------------------------------------------------------------------'
end if


time_advance=.true.
! do time integration of all grids on all levels

call timeintegration

if (mype==0) then
   print*,'-----------------------------------------------------------------------------'
   write(*,'(a,f17.3,a)')' Finished BHAC in : ',MPI_WTIME()-time0,' sec'
   print*,'-----------------------------------------------------------------------------'
end if

{#IFDEF PARTICLES
call finish_tracerparticles
}
call comm_finalize

end program amrvac
!=============================================================================
subroutine timeintegration
{#IFDEF DY_SP
use mod_cfc_parameters
}
use mod_timing
{^IFTHREED
   use mod_healpix_det, only: destroy_healpix_detectors, initialize_healpix_detectors
}
include 'amrvacdef.f'

integer :: level, ifile, fixcount
double precision :: time_all_in, time_all_end
logical, external :: timetosave, fixgrid
integer :: iigrid, igrid, total_cell
double precision :: max_local(10), max_global(10)
double precision :: min_local(1), min_global(1)
double precision :: time_last_print
logical :: allow_output
{#IFDEF SAVENOW
logical :: alive
}
!-----------------------------------------------------------------------------
time_last_print = -bigdouble

it_start = 0
if (snapshotini/=-1) then
  ! When restart, decide to output t=0 (restart moment)
  if (need_restart_output) then
    allow_output = .true.
  else
    allow_output = .false.
  endif
else
  allow_output = .true.
endif

{#IFDEF DY_SP
if (use_gw_br) then
  !$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     pw_t_dot_old2(igrid)%w(ixG^T,:) = 0.0d0
     pw_t_dot_old(igrid)%w(ixG^T,:)  = 0.0d0
  end do
  !$OMP END PARALLEL DO
endif
}

time_in=MPI_WTIME()
timeio0=MPI_WTIME()
fixcount=1

n_saves(filelog_:fileout_) = snapshotini
n_saves_rt(filelog_:fileout_) = time_in

do ifile=nfile,1,-1
   if (time_accurate) tsavelast(ifile)=t
   tsavelast_rt(ifile) = timeio0 - time_in
   itsavelast(ifile)=it
end do

itmin=it
! the next two are used to keep track of the performance during runtime:
itTimeLast=it
timeLast=MPI_WTIME()

call getbc(t,ps,psCoarse)


!  ------ start of integration loop. ------------------
timeloop0=MPI_WTIME()
timefirstbc=timeloop0-time_in
time_bc=0.d0
if (mype==0) then
   write(*,'(a,f12.3,a)')&
      ' BCs before Advance took : ',timefirstbc,' sec'
end if
if (mype == 0) then
  write(*, *) "#         it,        t,        dt,       reality time"
endif

time_evol : do

   ! Select active grids
   call selectgrids

   ! Save previous dt
   dt_old3 = dt_old2
   dt_old2 = dt_old
   dt_old  = dt

   call setdt

   if(fixprocess) call process(it,t)

   timeio0=MPI_WTIME()

   if (allow_output) then
     do ifile=nfile,1,-1
        if(timetosave(ifile,timeio0-time_in)) call saveamrfile(ifile)
     end do
     ! later output must be allowed
   endif
   allow_output = .true.

   if (timeio0 - time_last_print > time_between_print) then
     time_last_print = timeio0
     if (mype == 0) then
       write(*, '(A4,I10,ES12.4,ES12.4,ES12.4)') " #", &
            it, t, dt, timeio0 - time_in
     end if
   end if


{#IFDEF SAVENOW   inquire(file='savenow',exist=alive)
   if(alive) then
     if(mype==0) write(*,'(a,i7,a,i7,a,es12.4)') ' save a snapshot No.',&
                      snapshot,' at it=',it,' t=',t
     call saveamrfile(1)
     call saveamrfile(2)
     call MPI_FILE_DELETE('savenow',MPI_INFO_NULL,ierrmpi)
   endif}
   timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)
   ! exit time loop criteria
   if (it>=itmax .or. (it-itmin)>=continue_it) exit time_evol
   if ((time_accurate .and. t>=tmax) .or. (t-tin)>=continue_t) exit time_evol

    !time_all_in = MPI_WTIME()  
   call advance(it)

   if((.not.time_accurate).or.(residmin>smalldouble)) then
      call getresidual(it)
   endif 
   if (residual<residmin .or. residual>residmax) exit time_evol

   ! resetting of tree BEFORE IO and setdt
   timegr0=MPI_WTIME()
   if(ditregrid>1) then
     if(fixcount<ditregrid) then
       fixcount=fixcount+1
     else
       if (mxnest>1 .and. .not.(fixgrid(0))) then
          call resettree
          {^IFTHREED
             if (healpix_det_open) call initialize_healpix_detectors()
          }
       endif
       fixcount=1
     endif
   else
     if (mxnest>1 .and. .not.(fixgrid(0))) then
        call resettree
        {^IFTHREED
           if (healpix_det_open) call initialize_healpix_detectors()
        }
     endif
   endif
   timegr_tot=timegr_tot+(MPI_WTIME()-timegr0)

   ! next timestep
   it = it + 1
   it_start = it_start + 1

   if (time_accurate) t = t + dt
   if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
   
   if(it>9000000)then
     it = slowsteps+10
     itmin=0
     itsavelast(:)=0
   end if

    !time_all_end = MPI_WTIME()  

    !write(*,*)  'time_all_end - time_all_in'
    !write(*,*)  time_all_end - time_all_in
end do time_evol

timeloop=MPI_WTIME()-timeloop0

if (mype==0) then
   write(*,'(a,f12.3,a)')' Total timeloop took        : ',timeloop,' sec'
   write(*,'(a,f12.3,a)')' Time spent on Regrid+Update: ',timegr_tot,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timegr_tot/timeloop,' %'
   write(*,'(a,f12.3,a)')' Time spent on IO in loop   : ',timeio_tot,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timeio_tot/timeloop,' %'
   write(*,'(a,f12.3,a)')' Time spent on BC           : ',time_bc,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*time_bc/timeloop,' %'
end if

timeio0=MPI_WTIME()
do ifile=nfile,1,-1
   if(itsavelast(ifile)<it)call saveamrfile(ifile)
enddo
if (mype==0) call MPI_FILE_CLOSE(log_fh,ierrmpi)
timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

{^IFTHREED
   if (healpix_det_open) call destroy_healpix_detectors
}

{#IFDEF RAY
call time_spent_on_rays
}
{#IFDEF PARTICLES
call time_spent_on_particles
}

if (mype==0) then
   write(*,'(a,f12.3,a)')' Total time spent on IO     : ',timeio_tot,' sec'
   write(*,'(a,f12.3,a)')' Total timeintegration took : ',MPI_WTIME()-time_in,' sec'
end if

end subroutine timeintegration
!=============================================================================
logical function timetosave(ifile, real_time)

! Save times are defined by either tsave(isavet(ifile),ifile) or
! itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile) or dtsave_rt(ifile)
! Other conditions may be included.

include 'amrvacdef.f'

integer:: ifile
logical:: oksave
double precision :: real_time
!-----------------------------------------------------------------------------
oksave=.false.
if (it==itsave(isaveit(ifile),ifile)) then
   oksave=.true.
   isaveit(ifile)=isaveit(ifile)+1
end if
if (it==itsavelast(ifile)+ditsave(ifile)) oksave=.true.
if (time_accurate) then
   if (t>=tsave(isavet(ifile),ifile)) then
      oksave=.true.
      isavet(ifile)=isavet(ifile)+1
   end if
   if (t>=tsavelast(ifile)+dtsave(ifile)-smalldouble)then
     oksave=.true.
     n_saves(ifile) = n_saves(ifile) + 1
   endif
end if
if (real_time>=tsavelast_rt(ifile)+dtsave_rt(ifile)-smalldouble)then
  oksave=.true.
  n_saves_rt(ifile) = n_saves_rt(ifile) + 1
endif

if (oksave) then
   if (time_accurate) tsavelast(ifile) =t
   tsavelast_rt(ifile) = real_time
   itsavelast(ifile)=it
end if
timetosave=oksave

return
end function timetosave
!=============================================================================
logical function fixgrid(dummy)

! fixing the grid at given times or given iteration is defined by either 
! tfixgrid or itfixgrid
! Other conditions may be included (dummy input integer unused).

include 'amrvacdef.f'

integer:: dummy
!-----------------------------------------------------------------------------

fixgrid= (t>=tfixgrid .or. it>=itfixgrid)

return
end function fixgrid
!=============================================================================
subroutine initglobal

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

call initglobaldata
call initglobaldata_usr
call checkglobaldata

end subroutine initglobal
!=============================================================================
