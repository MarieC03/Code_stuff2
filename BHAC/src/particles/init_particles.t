!=============================================================================
subroutine init_particle_integrator()

  use mod_particles
  use constants
  include 'amrvacdef.f'
  !-----------------------------------------------------------------------------

  itmax_particles = 10000000
  tmax_particles  = 10.0d0! * CONST_years
  dtsync_particles    = 0.1d0
  ditsave_particles = 1
  dtsave_ensemble   = dtsave(2)! * UNIT_LENGTH/UNIT_VELOCITY
  dtheta            = 2.0d0 * dpi / 60.0d0

  losses = .false.

end subroutine init_particle_integrator
!=============================================================================
subroutine init_particles()
  ! initialise the particles

  use constants
  use mod_particles
  use mod_gridvars
  use mod_random
  include 'amrvacdef.f'

  double precision, dimension(ndir)       :: x
  integer                                 :: igrid_particle, ipe_particle
  integer, parameter                      :: Npart=100
  integer(i8)                             :: seed(2)
  type(rng_t)                             :: myrand
  double precision                        :: r^D(1:Npart)
  double precision                        :: v(1:ndir)
  logical, dimension(1:Npart)             :: follow=.false.
  type(particle_node), dimension(1:Npart) :: particles_to_inject
  integer                                 :: ninject, nparticles_added
  double precision, dimension(1:3)        :: tmp_payload
  !-----------------------------------------------------------------------------

  ! initialise the random number generator
  seed = [310952_i8,24378_i8]
  call myrand%set_seed(seed)
  {^D& call myrand%unif_01_vec(r^D)\}

  ! flags to follow particles
  follow(31)   = .true.
  follow(81)  = .true.
  follow(49)  = .true.


  ninject = 0
  nparticles_added = 0
  x=0
  do while (nparticles_added .lt. Npart)
     nparticles_added = nparticles_added + 1

     {^D&x(^D) = xprobmin^D + r^D(nparticles_added) * (xprobmax^D - xprobmin^D)\}

     call find_particle_ipe(x,igrid_particle,ipe_particle)

     if (ipe_particle == mype) then 
        ninject = ninject + 1

        particles_to_inject(ninject)%follow = follow(nparticles_added)
        particles_to_inject(ninject)%q      = - CONST_e
        particles_to_inject(ninject)%m      =   CONST_me

        particles_to_inject(ninject)%t      = 0.0d0
        particles_to_inject(ninject)%dt     = bigdouble !important, or the particle will not be activated

        particles_to_inject(ninject)%x(:)   = x

        {#IFDEF PARTICLES_ADVECT
        call get_vec(igrid_particle,x,0.0d0,v,vp1_,vp^NC_)
        particles_to_inject(ninject)%u(:) = v(:)

        call get_vec(igrid_particle,x,0.0d0,tmp_payload,rhop_,bb2p_)
        particles_to_inject(ninject)%payload(:) = tmp_payload(:)
        }

     end if

  end do

  call inject_particles(particles_to_inject(1:ninject))

end subroutine init_particles
!=============================================================================
subroutine add_particles()
  ! add some more particles

  use constants
  use mod_particles
  use mod_gridvars
  use mod_random
  include 'amrvacdef.f'

  double precision, dimension(ndir)       :: x
  integer                                 :: igrid_particle, ipe_particle
  integer, parameter                      :: Npart=10
  type(rng_t)                             :: myrand
  integer(i8)                             :: seed(2)
  double precision                        :: r^D(1:Npart)
  double precision                        :: v(1:ndir)
  type(particle_node), dimension(1:Npart) :: particles_to_inject
  integer                                 :: ninject, nparticles_added
  double precision, dimension(1:3)        :: tmp_payload
  double precision, parameter             :: dt_inject=0.1d0
  double precision, save                  :: tinject_last=0.0d0
  !-----------------------------------------------------------------------------

  ! for static snapshot: inject every dtsync_particles interval
  if (.not. time_advance .or. time_advance .and. t_particles - tinject_last .gt. dt_inject) then
     
     seed = [310952_i8,24378_i8] + it_particles
     call myrand%set_seed(seed)
     {^D& call myrand%unif_01_vec(r^D)\}

     ninject = 0
     nparticles_added = 0
     x=0
     do while (nparticles_added .lt. Npart)
        nparticles_added = nparticles_added + 1

        {^D&x(^D) = xprobmin^D + r^D(nparticles_added) * (xprobmax^D - xprobmin^D)\}

        call find_particle_ipe(x,igrid_particle,ipe_particle)

        if (ipe_particle == mype) then 
           ninject = ninject + 1

           particles_to_inject(ninject)%follow = .false.
           particles_to_inject(ninject)%q      = - CONST_e
           particles_to_inject(ninject)%m      =   CONST_me

           particles_to_inject(ninject)%t      = t_particles
           particles_to_inject(ninject)%dt     = bigdouble ! important, or the particle will not be activated

           particles_to_inject(ninject)%x(:)   = x

           {#IFDEF PARTICLES_ADVECT
           call get_vec(igrid_particle,x,0.0d0,v,vp1_,vp^NC_)
           particles_to_inject(ninject)%u(:) = v(:)

           call get_vec(igrid_particle,x,0.0d0,tmp_payload,rhop_,bb2p_)
           particles_to_inject(ninject)%payload(:) = tmp_payload(:)
           }

        end if

     end do

     call inject_particles(particles_to_inject(1:ninject))

     tinject_last = t_particles

  end if

end subroutine add_particles
!=============================================================================
logical function user_destroy(myparticle)

  use mod_particles, only: particle_node

  type(particle_node), intent(in)          :: myparticle
  !-----------------------------------------------------------------------------
  user_destroy = .false.

end function user_destroy
!=============================================================================

