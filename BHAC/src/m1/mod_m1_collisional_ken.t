module mod_m1_collisional
  use mod_m1_internal
  use mod_m1_eas
  use mod_m1_closure
  {#IFDEF UNIT_TESTS
  use mod_m1_metric_interface, only: m1_metric_helper
  }
  implicit none

  integer :: ix_loc^D 
  integer, parameter :: m1_system_ndim = m1_numvars_internal
  double precision, dimension(1:m1_system_ndim) :: Ut_explicit ! Conservatives (E, Fi)
  double precision :: qdt_impl
  double precision :: alp, sqrtg
  double precision :: Fmag, fac
  double precision :: fmaxfact = 0.999d0 !0.99999d0 !0.999d0
  double precision, dimension(m1_num_eas) :: eas
  integer :: linearized_case
  
  type(m1_closure_helpers) :: stateM1 

  {#IFDEF UNIT_TESTS
  type(m1_metric_helper) :: metricM1_unit 
  }
  
  
contains
  !> calculate and add the collisional sources bar{G}(U) to the evolution eq.
  !> implicit update and solving for M1-variables U via multidimensional rootfinder
  !**********************************************************************************
  ! When this routine is called we assume the conservative variables to be already U_explicit  
  subroutine m1_get_implicit_collisional_sources(wrad, metricM1, vel_, eas_,ixI^L,ix^D, N, sources, qdt, speciesKSP, dtfactor)
    use mod_rootfinding, only: rootfinding_global_multid_newton_raphson, multistart_rootfinding!, rootfinding_multid_newton_raphson
    use mod_m1_metric_interface
    use mod_m1_internal
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    integer, intent(in) :: ix^D,ixI^L
    integer, intent(in) :: speciesKSP
    !> Note: input: wrad=wrad/sqrt, output: wrad=wrad*sqrt
    double precision, intent(inout) :: wrad(1:m1_numvars_internal) !E,Fi
    double precision, intent(in)    :: eas_(1:m1_num_eas)    ! m1_eas only on given ixD
    double precision, intent(in)    :: vel_(1:^NC)  
    double precision, intent(inout) :: N                     ! nrad = Gamma*nprim
    double precision, intent(out)   :: sources(m1_system_ndim+1) ! dE dS_i dYE
    double precision, intent(in)    :: qdt       !> time step
    double precision, intent(in)    :: dtfactor   !< Timestep factor 
	  type(m1_metric_helper), intent(in)  :: metricM1 
    ! internal
    double precision :: dYE
    double precision :: T_eq, Ye_eq
    integer :: iter, error
    integer :: i
    !------ for testing:
    logical :: linearized_eqs = .false. !true.
    
     linearized_case = 0   ! no linearization
     !linearized_case = 1   ! linearization around v_i = 0
     !linearized_case = 2   ! linearization around U_0 = 0
     !linearized_case = 3   ! linearization around v_i = 0 and k_s active
     !linearized_case = 4   ! set coll_sources to some value
    !-------

    ! Store variables needed for the iteration
    !Ut_explicit(1) = wrad(m1_energy_)
    !{C Ut_explicit(1+^C) = wrad(m1_flux^C_) }
     ! TEST
   ! Ut_explicit(1) = wrad(m1_energy_)/ metricM1%sqrtg(ix^D)
    !{C Ut_explicit(1+^C) = wrad(m1_flux^C_)/ metricM1%sqrtg(ix^D) }
   
    wrad(m1_energy_) = wrad(m1_energy_) / metricM1%sqrtg(ix^D)
    {^C& wrad(m1_flux^C_) = wrad(m1_flux^C_) / metricM1%sqrtg(ix^D) \}    

    !------------------------for singularieties
    if(metricM1%sqrtg(ix^D).eq.0.0d0) then
      wrad(m1_energy_) = 0.0d0
      {^C& wrad(m1_flux^C_) = 0.0d0 \}
    end if 
	  !--------------------------------------
    
    !> fill in stateM1, 
	  stateM1%E = wrad(m1_energy_)
    {^C& stateM1%F_low(^C) = wrad(m1_flux^C_) \} 
    {^C& stateM1%vel(^C)   = vel_(^C) \} 

    !! KEN: comment this!!
    !!do i=1,^NC
    !!  !if (dabs(wrad(m1_flux1_+(i-1))) < m1_E_atmo) then
    !!  if ((dabs(wrad(m1_flux1_+(i-1))) < m1_E_atmo) .and. (wrad(m1_flux1_+(i-1)) >= 0)) then
    !!     !stateM1%F_low(i) = max(m1_E_atmo, 1d-10*stateM1%E) !m1_E_atmo !0.0d0
    !!     stateM1%F_low(i) = 0.1d0 * m1_E_atmo
    !!  else if ((dabs(wrad(m1_flux1_+(i-1))) < m1_E_atmo) .and. (wrad(m1_flux1_+(i-1)) < 0)) then
    !!      !stateM1%F_low(i) = -1.0d0 * max(m1_E_atmo, 1d-10*stateM1%E) !m1_E_atmo !0.0d0
    !!      stateM1%F_low(i) = -0.1d0 * m1_E_atmo
    !!  else
    !!     stateM1%F_low(i) = wrad(m1_flux1_+(i-1))
    !!  end if
    !!end do    
    !!! and enforce Fmag <= fmax * E (typical fmax=1)
    !!Fmag = dsqrt(sum(stateM1%F_low(:)**2))
    !!if (Fmag > fmaxfact * stateM1%E) then
    !!  fac = (fmaxfact * stateM1%E) / Fmag
    !!  stateM1%F_low(:) = fac * stateM1%F_low(:)
    !!end if    


    !> update closure
	  call m1_update_closure_ixD(stateM1,metricM1,ix^D,.true.,get_vel_impl=.true.)
    !> Now J, H_i are up to date, fill wrad with update values
    wrad(m1_energy_) = stateM1%E
    {^C& wrad(m1_flux^C_) = stateM1%F_low(^C) \}

    Ut_explicit(1) = stateM1%E
    {^C& Ut_explicit(1+^C) = stateM1%F_low(^C) \}

    ! Store ix^D since it's needed during the iteration
    ! for raising / lowering indices
    {^D& ix_loc^D = ix^D \}
    qdt_impl = qdt*dtfactor 
    alp = metricM1%alp(ix^D)
    sqrtg = metricM1%sqrtg(ix^D)
    eas = eas_
    
    ! Compute initial guess for iteration
    call get_iter_initguess(wrad)  
   
    !--------------------------------------
    ! multidimensional newton raphson rootfinding
    if(linearized_eqs) then
      !linearized_eqs for testing coll_sources
      call rootfinding_global_multid_newton_raphson(m1_system_ndim, wrad, 1.0d-13, 150, &
       error, m1_sources_func_Lin,m1_sources_jacobian_Lin, m1_sources_min_Lin )
    else
       call rootfinding_global_multid_newton_raphson(m1_system_ndim, wrad, 1.0d-15, 600, &
          error, m1_sources_func,m1_sources_jacobian, m1_sources_min )
    end if 

    if ( error .ne. 0 ) then
       if( (error == 2) .or. (error == 4) .or.(error == 5) ) then
        call multistart_rootfinding(m1_system_ndim, wrad, 1.0d-15, 150, &
        error, m1_sources_func,m1_sources_jacobian, m1_sources_min, get_random_start )
         ! Use initial guess if iterative solution did not converge
        if( error .ne. 0 ) then
         call get_iter_initguess(wrad)
        end if
       end if
       
    end if
    !--------------------------------------

    ! fill state with updated wrad values
! KEN changed i to to C
    stateM1%E = max(wrad(m1_energy_),m1_E_atmo) 
    {^C& stateM1%F_low(^C) = wrad(m1_flux^C_) \}
    
    !!! KEN: comment this!!
    !!do i=1,^NC
    !!  !if (dabs(wrad(m1_flux1_+(i-1))) < m1_E_atmo) then
    !!  if ((dabs(wrad(m1_flux1_+(i-1))) < m1_E_atmo) .and. (wrad(m1_flux1_+(i-1)) >= 0)) then
    !!     !stateM1%F_low(i) = max(m1_E_atmo, 1d-10*stateM1%E) !m1_E_atmo !0.0d0
    !!     stateM1%F_low(i) = 0.1d0 * m1_E_atmo 
    !!  else if ((dabs(wrad(m1_flux1_+(i-1))) < m1_E_atmo) .and. (wrad(m1_flux1_+(i-1)) < 0)) then
    !!      !stateM1%F_low(i) = -1.0d0 * max(m1_E_atmo, 1d-10*stateM1%E) !m1_E_atmo !0.0d0
    !!      stateM1%F_low(i) = -0.1d0 * m1_E_atmo 
    !!  else
    !!     stateM1%F_low(i) = wrad(m1_flux1_+(i-1))
    !!  end if
    !!end do    
    !!! and enforce Fmag <= fmax * E (typical fmax=1)
    !!Fmag = dsqrt(sum(stateM1%F_low(:)**2))
    !!if (Fmag > fmaxfact * stateM1%E) then
    !!  fac = (fmaxfact * stateM1%E) / Fmag
    !!  stateM1%F_low(:) = fac * stateM1%F_low(:)
    !!end if    

    ! update the closure to make sure everything is consistent
	  !call m1_update_closure_ixD(stateM1,metricM1,ix^D,.false.,get_vel_impl=.false.)
	  call m1_update_closure_ixD(stateM1,metricM1,ix^D,.true.,get_vel_impl=.false.)
    
    !!!!!!! TEST
    if((stateM1%E - stateM1%Fv) < 0.0d0)then
      write(888,*)"Aiaiaiai",ix^D, stateM1%E - stateM1%Fv
      write(888,*)"E",stateM1%E 
      write(888,*)"E",stateM1%Fv
      write(*,*)"Aiaiaiai",ix^D, stateM1%E - stateM1%Fv
      write(*,*)"E",stateM1%E 
      write(*,*)"E",stateM1%Fv
    end if 

    ! fill wrad with updated state values 
    wrad(m1_energy_) = stateM1%E
    {^C& wrad(m1_flux^C_) = stateM1%F_low(^C) \}
    !! TEST:   
    !{C wrad(m1_flux^C_) = 0.0d0 }   

    ! fill source terms -- updates closure too! 
    if(linearized_eqs) then
      sources(1:m1_system_ndim) = m1_collisional_sources_Lin(wrad)
    else
      sources(1:m1_system_ndim) = m1_collisional_sources(wrad)
    end if 
    ! status: moved alpha out of collisional_sources-function
    sources = sources*sqrtg*alp

    wrad(m1_energy_) = wrad(m1_energy_)*sqrtg
    {^C& wrad(m1_flux^C_) = wrad(m1_flux^C_)*sqrtg \}
    
    !------------------------ ignore this 
    if(metricM1%sqrtg(ix^D).eq.0.0d0) then
      wrad(m1_energy_) = 0.0d0
      {^C& wrad(m1_flux^C_) = 0.0d0 \}
    end if 
    !-----------------------------------

    ! Finally compute N and YE update
    call m1_get_implicit_n_source(N,dYE)
  
    ! for adding dYE to sources 
    sources(m1_system_ndim+1) = dYE * sqrtg * alp   

   contains

  function get_random_start(z_vec)
    double precision, dimension(:), intent(inout) :: z_vec
    double precision, dimension(1:size(z_vec)) :: get_random_start
    integer :: i_z,n_z, clock
    integer, allocatable :: seed(:)
    double precision :: rand_num
    double precision :: z_low, z_high
       z_low = 0.0d0 + 1.d-40
       z_high = 1.0d0  

      ! Initialize the random seed
       call system_clock(count=clock)
       call random_seed(size=n_z)
       allocate(seed(n_z))
       seed = clock + 37 * (/ (i_z - 1, i_z = 1, n_z) /)
       call random_seed(put=seed)
       deallocate(seed)

      ! Generate and print 4 random numbers between 0 and 1
       do i_z = 1,m1_numvars_internal !size(z_vec)
          call random_number(rand_num)
          z_vec(i_z) = z_low + (z_high - z_low) * rand_num
          get_random_start(i_z) = 1.0d0
       end do
      
  end function get_random_start

  function m1_collisional_jacobian(z_vec)
    use mod_m1_eas, only: Q_ems, k_a, k_s
    use mod_m1_internal
    use mod_m1_closure
    double precision, dimension(1:m1_system_ndim,1:m1_system_ndim) :: m1_collisional_jacobian
    double precision, intent(in) :: z_vec(1:m1_system_ndim)
    !internal
    integer :: i,j
    double precision :: W2,W3
    double precision :: JvF, JfF, HvE, HfE, HdF, HvvF, HffF, HfvF, HvfF
    double precision :: dJdE
    double precision, dimension(^NC) :: dJdF, dHdE
    double precision, dimension(^NC,^NC) :: dHdF
    integer :: is_paper_musolino = .false.
    integer :: is_paper_radice = .true.

    stateM1%E = z_vec(1)
    {^C& stateM1%F_low(^C) = z_vec(1+^C) \} 
    ! update the closure 
    call m1_update_closure_ixD(stateM1,metricM1,ix_loc^D,.true.,get_vel_impl=.false.)
    ! Ken removed this closure in the high hopes that the module remembers it

    W2 = stateM1%W * stateM1%W
    W3 = W2 * stateM1%W
    
    !--- JvF -------------
    if(is_paper_musolino) then
      ! paper musolino:
      JvF = 2*stateM1%W**2*(-1+stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden - 2.0d0*stateM1%dthick*(W2-1.0d0)/(1.0d0 + 2.0d0*W2) )
    else if(is_paper_radice) then
      ! paper radice:
      JvF = 2*stateM1%W**2*(-1 + stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden + 2.0d0*stateM1%dthick*(W2-1.0d0)/(1.0d0 + 2.0d0*W2) )
    end if 

    JfF = -2.0d0*stateM1%dthin*W2*stateM1%E*stateM1%fhatv**2/stateM1%Fden

    HvE = W3*(-1-stateM1%dthin*stateM1%fhatv**2 + stateM1%dthick*(2.0d0*W2-3.0d0)/(1.0d0+2.0d0*W2) )
    HfE = -stateM1%dthin*stateM1%W*stateM1%fhatv

    HdF = stateM1%W*(1-stateM1%dthick*stateM1%v2-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden)
    !--- HvvF -------------
    if(is_paper_musolino) then
      ! paper musolino:
      HvvF = 2.0d0*W3*(1.0d0-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden-stateM1%dthick*(1.0d0 -1.0d0/(2.0d0*W2*(1.0d0+2.0d0*W2))))
    else if(is_paper_radice) then
      ! paper radice:
      HvvF = 2.0d0*W3*(1.0d0-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden-stateM1%dthick*(stateM1%v2+1.0d0/(2.0d0*W2*(1.0d0+2.0d0*W2))))
    end if 
    !--- HffF ------------
    HffF = 2.0d0*stateM1%dthin*stateM1%W*stateM1%E*stateM1%fhatv/stateM1%Fden
    !--- HvfF -------------
    if(is_paper_musolino) then
      ! paper musolino:
      HvfF = 2.0d0*stateM1%dthin*W3*stateM1%E*stateM1%fhatv/stateM1%Fden
    else if(is_paper_radice) then
      ! paper radice:
      HvfF = 2.0d0*stateM1%dthin*W3*stateM1%E*stateM1%fhatv**2/stateM1%Fden
    end if
    !--- HfvF -----------
    if(is_paper_musolino) then
      ! paper musolino:
      HfvF = -stateM1%dthin*stateM1%W*stateM1%E*stateM1%fhatv/stateM1%Fden
    else if(is_paper_radice) then
      ! paper radice:
      HfvF = -stateM1%dthin*stateM1%W*stateM1%E/stateM1%Fden
    end if

    dJdE = W2+stateM1%dthin*stateM1%fhatv**2*W2+stateM1%dthick*((3.0d0-2.0d0*W2)*(W2-1.0d0))/(1.0d0+2.0d0*W2)
    dJdF(:) = JvF*stateM1%vel(:)  + JfF*stateM1%fhat_up(:)
    dHdE(:) = HvE*stateM1%vlow(:) + HfE*stateM1%fhat_low(:)

    dJdE = max(0.d0, min(1.d0 - 1.d-12, dJdE))

    if((eas(k_a)+eas(k_s) - eas(k_s)*dJdE) .le. 0.0d0) then
      write(100,*)"jac: ks*dJdE too large --------------------------"
      write(100,*)"k_a+k_s-ks*dJdE ",eas(k_a)+eas(k_s) - eas(k_s)*dJdE
      write(100,*)"dJdE",dJdE
      write(100,*)"k_a+k_s",eas(k_a)+eas(k_s)
      write(100,*)"k_s*dJdE",eas(k_s)*dJdE
      write(100,*) stateM1%fhatv
    end if 

    do i=1, ^NC
       do j=1, ^NC
          dHdF(i,j) = HvvF*stateM1%vlow(i)*stateM1%vel(j) &
               + HffF * stateM1%fhat_low(i)*stateM1%fhat_up(j)  &
               + HvfF * stateM1%vlow(i) *stateM1%fhat_up(j)  &
               + HfvF * stateM1%fhat_low(i)*stateM1%vel(j)
          if( i .eq. j) dHdF(i,j) = dHdF(i,j) + HdF
          
         !!! TEST
         ! replace dHdF with purely diagonal, isotropic damping
          !dHdF(i,j) = 0.0d0
          !if (i == j) dHdF(i,j) = stateM1%W   ! just a minimal positive damping scale
     
       enddo
    enddo
    ! ------ fill jacobian --------
    !! OG:
    m1_collisional_jacobian(1,1) = -1.d0* stateM1%W*( (eas(k_a)+eas(k_s)) - eas(k_s)*dJdE )
    ! TEST 1:
    !m1_collisional_jacobian(1,1) = -1.d0* stateM1%W*( (eas(k_a)) - eas(k_s)*dJdE )
    !TEST 2:
    !m1_collisional_jacobian(1,1) = 1.d0* stateM1%W*(eas(k_s)*dJdE )
    ! TEST 3: kill k_s parts:
    !m1_collisional_jacobian(1,1) = -1.d0* stateM1%W*( (eas(k_a)+eas(k_s)))
    do i=1, ^NC
       !! OG:
       m1_collisional_jacobian(1,i+1) = stateM1%W* ( eas(k_s)*dJdF(i) + (eas(k_a)+eas(k_s))*stateM1%vel(i) )
       ! TEST 1:
       !m1_collisional_jacobian(1,i+1) = stateM1%W* ( eas(k_s)*dJdF(i) + (eas(k_a))*stateM1%vel(i) )
       ! TEST 2:
       !!m1_collisional_jacobian(1,i+1) = stateM1%W* ( eas(k_s)*dJdF(i) )
       ! TEST 3:
       !m1_collisional_jacobian(1,i+1) = stateM1%W* ( (eas(k_a)+eas(k_s))*stateM1%vel(i) )
       !! OG:
       m1_collisional_jacobian(i+1,1) = -1.d0*( (eas(k_a)+eas(k_s))*dHdE(i) + stateM1%W*eas(k_a)*dJdE*stateM1%vlow(i) )
       ! TEST 3:
       !m1_collisional_jacobian(i+1,1) = -1.d0*( (eas(k_a))*dHdE(i) + stateM1%W*eas(k_a)*dJdE*stateM1%vlow(i) )
    enddo
    do i=1, ^NC
       do j=1, ^NC
          m1_collisional_jacobian(i+1,j+1) = -1.d0*( (eas(k_a)+eas(k_s))*dHdF(i,j) + stateM1%W*eas(k_a)*stateM1%vlow(i)*dJdF(j) )
          ! TEST3:
          !m1_collisional_jacobian(i+1,j+1) = -1.d0*( (eas(k_a))*dHdF(i,j) + stateM1%W*eas(k_a)*stateM1%vlow(i)*dJdF(j) )
       enddo
    enddo
    ! done ! 
  end function m1_collisional_jacobian


  function m1_collisional_sources(z_vec)
    ! bar{G}
    use mod_m1_eas, only: Q_ems, k_a, k_s
    use mod_m1_internal
    use mod_m1_closure 
    double precision, dimension(1:m1_system_ndim) :: m1_collisional_sources
    double precision, intent(in) :: z_vec(1:m1_system_ndim)
    double precision :: Hn
    double precision :: deltaFrame, rel
    stateM1%E = z_vec(1)
    {^C& stateM1%F_low(^C) = z_vec(1+^C) \}
    ! update the closure  
    call m1_update_closure_ixD(stateM1,metricM1,ix_loc^D,.true.,get_vel_impl=.false.,use_init_zeta=.false.)

    Hn = 0.d0
    {^C& Hn = Hn + stateM1%vel(^C)*stateM1%H_low(^C) \}
    !!
   !! deltaFrame = (stateM1%E-stateM1%Fv) - stateM1%W*(stateM1%J + Hn)
   !! if(dabs(deltaFrame) > 1.0d-12) then
   !!   write(72,*)"---deltaFrame large",deltaFrame
   !!   write(72,*)"E-Fv",(stateM1%E-stateM1%Fv)
   !!   write(72,*)"W(J+Hn)",stateM1%W*(stateM1%J + Hn)
   !!   write(72,*)"E,W,Fv",stateM1%E,stateM1%W,stateM1%Fv
   !!   write(72,*)"J,Hn",stateM1%J,Hn
   !! end if 
   !!  rel = abs((stateM1%E - stateM1%Fv) - stateM1%W*(stateM1%J + Hn))/ max(1.d-30, stateM1%E + abs(stateM1%Fv))
   !! ! expect rel ≲ 1e-12 in smooth zones
   !!   if(rel .gt. 1.0d-12)then
   !!     write(*,*)"EF-JH not match: rel", rel
   !!     write(*,*)"E-Fv", (stateM1%E - stateM1%Fv)
   !!     write(*,*)"E,Fv",stateM1%E,stateM1%Fv
   !!     write(*,*)"W(J+Hn)", stateM1%W*(stateM1%J + Hn)
   !!     write(*,*)"J,Hn",stateM1%J,Hn
   !!   end if 

    ! If this does not work try setting k_s = 0 as gpt says to isolate issue e.g. F*v*kappa_s large
    !!!!m1_collisional_sources(1) = stateM1%W*( eas(Q_ems) + eas(k_s)*stateM1%J) - (eas(k_s)+eas(k_a))*Hn 
    ! Current Sources
    !m1_collisional_sources(1) = stateM1%W*( eas(Q_ems) + eas(k_s)*stateM1%J - (eas(k_s)+eas(k_a))*(stateM1%E-stateM1%Fv) )
    !{C m1_collisional_sources(1+^C) = stateM1%W*(eas(Q_ems)-eas(k_a)*stateM1%J)*stateM1%vlow(^C)&
    !    -(eas(k_a)+eas(k_s))*stateM1%H_low(^C) }

    ! OG:
    m1_collisional_sources(1) = stateM1%W*( eas(Q_ems) + eas(k_s)*stateM1%J - (eas(k_s)+eas(k_a))*(stateM1%E-stateM1%Fv) )
    !! TEST1 - kappa_a *(E-Fv) no kappa_S in energy source
    !m1_collisional_sources(1) = stateM1%W*( eas(Q_ems) + eas(k_s)*stateM1%J - (eas(k_a))*(stateM1%E-stateM1%Fv) )
    !! TEST2 cancel -kappa*(E-Fv) in energy source
    !m1_collisional_sources(1) = stateM1%W*( eas(Q_ems) + eas(k_s)*stateM1%J )
    {^C& m1_collisional_sources(1+^C) = stateM1%W*(eas(Q_ems)-eas(k_a)*stateM1%J)*stateM1%vlow(^C)&
         -(eas(k_a)+eas(k_s))*stateM1%H_low(^C) \}

  end function m1_collisional_sources


  function m1_sources_func(z_vec)
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec)) :: m1_sources_func
    ! internal
    double precision, dimension(1:m1_system_ndim) :: coll_sources

    coll_sources = m1_collisional_sources(z_vec) *alp * qdt_impl
   ! m1_sources_func = zero for root-findiing
   ! -> trying to find z_vec which is: zevc*sqrtg = (Ut_explicit=wrad*sqrtg)+ coll_source*sqrtg
   !  G = coll_sources = alp*sqrtg*G_bar 
    m1_sources_func(:) =  z_vec(:) -  coll_sources(:) - Ut_explicit(:) 
  end function m1_sources_func

  function m1_sources_jacobian(z_vec)
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec),1:size(z_vec)) :: m1_sources_jacobian
    ! internal
    integer :: i,j

    m1_sources_jacobian(:,:) = - m1_collisional_jacobian(z_vec) * qdt_impl *alp
    m1_sources_jacobian(1,1) = m1_sources_jacobian(1,1) + 1  
    {^C& m1_sources_jacobian(1+^C,1+^C) = m1_sources_jacobian(1+^C,1+^C) + 1  \}

    
  end function m1_sources_jacobian
  
  function m1_sources_min(z_vec)
    !use mod_metric
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    double precision :: m1_sources_min 
    double precision, dimension(:), intent(in) :: z_vec
    ! internal
    double precision, dimension(m1_system_ndim) :: Fsource_low, Fsource_hi

    !this sets the limits of the function for rootfinding
    Fsource_low = m1_sources_func(z_vec)
    !call raise_ixD(metricM1,ix_loc^D,Fsource_low(2:m1_system_ndim),Fsource_hi(2:m1_system_ndim))
    ! KEN I inlined this because raise_ixD is just slow
    Fsource_hi(2) = metricM1%gammaUPij(ix_loc^D,1)*Fsource_low(2) &
              + metricM1%gammaUPij(ix_loc^D,2)*Fsource_low(3) &
              + metricM1%gammaUPij(ix_loc^D,3)*Fsource_low(4)
    Fsource_hi(3) = metricM1%gammaUPij(ix_loc^D,2)*Fsource_low(2) &
              + metricM1%gammaUPij(ix_loc^D,4)*Fsource_low(3) &
              + metricM1%gammaUPij(ix_loc^D,5)*Fsource_low(4)
    Fsource_hi(4) = metricM1%gammaUPij(ix_loc^D,3)*Fsource_low(2) &
              + metricM1%gammaUPij(ix_loc^D,5)*Fsource_low(3) &
              + metricM1%gammaUPij(ix_loc^D,6)*Fsource_low(4)

    m1_sources_min = Fsource_low(1)*Fsource_low(1)
    {^C& m1_sources_min = m1_sources_min + Fsource_low(1+^C)*Fsource_hi(1+^C) \}
    !m1_sources_min = 0.5d0*m1_sources_min
    m1_sources_min = m1_sources_min**0.5
    
  end function m1_sources_min

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Linearized collisional sources:
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function m1_collisional_jacobian_Lin(z_vec)
    use mod_m1_eas, only: Q_ems, k_a, k_s
    use mod_m1_internal
    use mod_m1_closure
    double precision, dimension(1:m1_system_ndim,1:m1_system_ndim) :: m1_collisional_jacobian_Lin
    double precision, intent(in) :: z_vec(1:m1_system_ndim)
    !internal
    integer :: i,j
    double precision :: W2,W3
    double precision :: JvF, JfF, HvE, HfE, HdF, HvvF, HffF, HfvF, HvfF
    double precision :: dJdE
    double precision, dimension(^NC) :: dJdF, dHdE
    double precision, dimension(^NC,^NC) :: dHdF

    ! update the closure 
    call m1_update_closure_ixD(stateM1,metricM1,ix_loc^D,.false.,get_vel_impl=.false.)

    W2 = stateM1%W * stateM1%W
    W3 = W2 * stateM1%W
    
    if(linearized_case .eq. 1) then 
      m1_collisional_jacobian_Lin(:,:) = 0.0d0 
      do i = 1, m1_system_ndim
        do j = 1,m1_system_ndim
            if(i .eq. j) then
              m1_collisional_jacobian_Lin(i,j) = -1.0d0 * (eas(k_a) + eas(k_s))
              if(i .eq. 1 .and. j .eq. 1) then
                m1_collisional_jacobian_Lin(i,j) = -1.0d0 * (eas(k_a))
              end if 
            end if 
        end do          
      end do
    else if(linearized_case .eq. 2) then
        call mpistop("TODO jacobian algebra")
    else if(linearized_case .eq. 3) then
      JvF = 2*stateM1%W**2*(-1+stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden - 2.0d0*stateM1%dthick*(W2-1.0d0)/(1.0d0 + 2.0d0*W2) )
      !JvF = 2*stateM1%W**2*(-1-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden + 2.0d0*stateM1%dthick*(W2-1.0d0)/(1.0d0 + 2.0d0*W2) )
      JfF = -2.0d0*stateM1%dthin*W2*stateM1%E*stateM1%fhatv**2/stateM1%Fden
  
      HvE = W3*(-1-stateM1%dthin*stateM1%fhatv**2 + stateM1%dthick*(2.0d0*W2-3.0d0)/(1.0d0+2.0d0*W2) )
      HfE = -stateM1%dthin*stateM1%W*stateM1%fhatv
  
      HdF = stateM1%W*(1-stateM1%dthick*stateM1%v2-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden)
      HvvF = 2.0d0*W3*(1.0d0-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden-stateM1%dthick*(1.0d0 -1.0d0/(2.0d0*W2*(1.0d0+2.0d0*W2))))
      !HvvF = 2.0d0*W3*(1.0d0-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden-stateM1%dthick*(stateM1%v2+1.0d0/(2.0d0*W2*(1.0d0+2.0d0*W2))))
      HffF = 2.0d0*stateM1%dthin*stateM1%W*stateM1%E*stateM1%fhatv/stateM1%Fden
      HvfF = 2.0d0*stateM1%dthin*W3*stateM1%E*stateM1%fhatv/stateM1%Fden
      !HvfF = 2.0d0*stateM1%dthin*W3*stateM1%E*stateM1%fhatv**2/stateM1%Fden
      HfvF = -stateM1%dthin*stateM1%W*stateM1%E*stateM1%fhatv/stateM1%Fden
      !HfvF = -stateM1%dthin*stateM1%W*stateM1%E/stateM1%Fden
  
      dJdE = W2+stateM1%dthin*stateM1%fhatv**2*W2+stateM1%dthick*((3.0d0-2.0d0*W2)*(W2-1.0d0))/(1.0d0+2.0d0*W2)
      dJdF(:) = JvF*stateM1%vel(:)  + JfF*stateM1%fhat_up(:)
      dHdE(:) = HvE*stateM1%vlow(:) + HfE*stateM1%fhat_low(:)
  
      do i=1, ^NC
         do j=1, ^NC
            dHdF(i,j) = HvvF*stateM1%vlow(i)*stateM1%vel(j) &
                 + HffF * stateM1%fhat_low(i)*stateM1%fhat_up(j)  &
                 + HvfF * stateM1%vlow(i) *stateM1%fhat_up(j)  &
                 + HfvF * stateM1%fhat_low(i)*stateM1%vel(j)
            if( i .eq. j) dHdF(i,j) = dHdF(i,j) + HdF
         enddo
      enddo
      ! ------ fill jacobian --------
      m1_collisional_jacobian_Lin(1,1) = -1.d0* stateM1%W*( (eas(k_s)) - eas(k_s)*dJdE )
      do i=1, ^NC
         m1_collisional_jacobian_Lin(1,i+1) = stateM1%W*(eas(k_s))*( dJdF(i))
         m1_collisional_jacobian_Lin(i+1,1) = -1.d0*( (eas(k_s))*dHdE(i))
      enddo
      do i=1, ^NC
         do j=1, ^NC
            m1_collisional_jacobian_Lin(i+1,j+1) = -1.d0*( (eas(k_s))*dHdF(i,j))
         enddo
      enddo
    else if(linearized_case .eq. 4) then
      m1_collisional_jacobian_Lin(:,:) = 0.0d0 
    end if ! end if linearized

    ! done ! 
  end function m1_collisional_jacobian_Lin


  function m1_collisional_sources_Lin(z_vec)
    ! bar{G}
    use mod_m1_eas, only: Q_ems, k_a, k_s
    use mod_m1_internal
    use mod_m1_closure 
    double precision, dimension(1:m1_system_ndim) :: m1_collisional_sources_Lin
    double precision, intent(in) :: z_vec(1:m1_system_ndim)
    ! update the closure 
    call m1_update_closure_ixD(stateM1,metricM1,ix_loc^D,.false.,get_vel_impl=.false.)

    if(linearized_case .eq. 1) then
      m1_collisional_sources_Lin(1) = eas(Q_ems) -  (eas(k_s)+eas(k_a))*stateM1%E
      {^C& m1_collisional_sources_Lin(1+^C) = - (eas(k_s)+eas(k_a))*stateM1%H_low(^C) \}
    else if(linearized_case .eq. 2) then
      m1_collisional_sources_Lin(1) = stateM1%W*( eas(Q_ems) )
      {^C& m1_collisional_sources_Lin(1+^C) = stateM1%W*eas(Q_ems)*stateM1%vlow(^C) \}
    else if(linearized_case .eq. 3) then
      m1_collisional_sources_Lin(1) = stateM1%W*eas(k_s) *(stateM1%J - stateM1%E )
      {^C& m1_collisional_sources_Lin(1+^C) = - (eas(k_s))*stateM1%H_low(^C) \}
    else if(linearized_case .eq. 4) then
      m1_collisional_sources_Lin(1) = 0.0d0
      m1_collisional_sources_Lin(2) = 10.0d0 
      m1_collisional_sources_Lin(3) = 0.0d0 
      m1_collisional_sources_Lin(4) = 0.0d0 
    end if 
  end function m1_collisional_sources_Lin

  function m1_sources_func_Lin(z_vec)
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec)) :: m1_sources_func_Lin
    ! internal
    double precision, dimension(1:m1_system_ndim) :: coll_sources

    coll_sources = m1_collisional_sources_Lin(z_vec) *alp * qdt_impl
   ! m1_sources_func = zero for root-findiing
   ! -> trying to find z_vec which is: zevc*sqrtg = (Ut_explicit=wrad*sqrtg)+ coll_source*sqrtg
   !  G = coll_sources = alp*sqrtg*G_bar 
    m1_sources_func_Lin(:) = Ut_explicit(:) + sqrtg*( coll_sources(:) - z_vec(:)) 
  end function m1_sources_func_Lin

  function m1_sources_jacobian_Lin(z_vec)
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec),1:size(z_vec)) :: m1_sources_jacobian_Lin
    ! internal
    integer :: i,j

    m1_sources_jacobian_Lin(:,:) = m1_collisional_jacobian_Lin(z_vec) * qdt_impl *alp

    m1_sources_jacobian_Lin(1,1) = m1_sources_jacobian_Lin(1,1) - 1  
    {^C& m1_sources_jacobian_Lin(1+^C,1+^C) = m1_sources_jacobian_Lin(1+^C,1+^C) -1  \}

    m1_sources_jacobian_Lin = m1_sources_jacobian_Lin * sqrtg
    
  end function m1_sources_jacobian_Lin
  
  function m1_sources_min_Lin(z_vec)
    !use mod_metric
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    double precision :: m1_sources_min_Lin 
    double precision, dimension(:), intent(in) :: z_vec
    ! internal
    double precision, dimension(m1_system_ndim) :: Fsource_low, Fsource_hi

    Fsource_low = m1_sources_func(z_vec)
    call raise_ixD(metricM1,ix_loc^D,Fsource_low(2:m1_system_ndim),Fsource_hi(2:m1_system_ndim))

    m1_sources_min_Lin = Fsource_low(1)*Fsource_low(1)
    {^C& m1_sources_min_Lin = m1_sources_min_Lin + Fsource_low(1+^C)*Fsource_hi(1+^C) \}
    m1_sources_min_Lin = 0.5d0*m1_sources_min_Lin
    
  end function m1_sources_min_Lin
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end subroutine m1_get_implicit_collisional_sources

  ! Important! This routine      !
  ! Assumes that the closure is  ! 
  ! already up-to-date           ! 
  ! compute initial guess for rootfinder
  subroutine get_iter_initguess(z_vec)
    use mod_m1_closure
    use mod_m1_eas, only: k_a, k_s, Q_ems
    double precision, intent(inout) :: z_vec(1:m1_system_ndim)
    ! internal
    double precision :: Hn
    double precision :: W2

    W2 = stateM1%W * stateM1%W
        
    stateM1%J = (stateM1%J + qdt_impl/stateM1%W*eas(Q_ems) ) / ( 1.0d0 + eas(k_a)*qdt_impl/stateM1%W )
    {^C& stateM1%H_low(^C) = stateM1%H_low(^C)/(1.0d0 + qdt_impl/stateM1%W*(eas(k_a)+eas(k_s)) ) \}

    Hn = {^C& stateM1%vel(^C)*stateM1%H_low(^C)+}
    Hn = -1.0d0 * Hn 

    stateM1%E = stateM1%J/3.0d0*(4.0d0*W2-1.0d0) - 2.0d0*Hn*stateM1%W

    {^C& stateM1%F_low(^C) = stateM1%W*stateM1%H_low(^C) + ( 4.0d0/3.0d0*W2*stateM1%J - stateM1%W*Hn ) * stateM1%vlow(^C) \}
  
    z_vec(m1_energy_) = stateM1%E
    {^C&
    z_vec(m1_flux^C_) = stateM1%F_low(^C)
    \}
  end subroutine get_iter_initguess

  ! implicit update of number density!
  subroutine m1_get_implicit_n_source(N, dYE)
    use mod_m1_eas, only: k_n, Q_ems_n
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    double precision, intent(inout) :: N
    double precision, intent(out) :: dYE
    !internal
    double precision :: n_prim
    !Note! Q_ems_n : Q_ems for number density
    !Note! k_n : kappa for number density

    N = ( N + qdt_impl*alp*sqrtg*eas(Q_ems_n) )/ ( 1.0d0 + alp*qdt_impl*eas(k_n)/stateM1%Gamma) 
    N = max(N,m1_E_atmo) 
    n_prim = N / stateM1%Gamma / sqrtg

    dYE  = (eas(Q_ems_n) - eas(k_n)*n_prim )

  end subroutine m1_get_implicit_n_source
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! -- all below only for UNIT_TESTS ------------------------
  {#IFDEF UNIT_TESTS
  subroutine fill_globals_metricM1_unit(metricM1_unit2,wrad,vel_,qdt,ix^D,eas_)
    use mod_m1_metric_interface
    use mod_m1_tests
    use mod_m1_closure
    use mod_m1_internal
	  type(m1_metric_helper), intent(inout)  :: metricM1_unit2 
    integer, intent(in) :: ix^D
    double precision, intent(inout) :: wrad(1:m1_numvars_internal)
    double precision, intent(in)    :: eas_(1:m1_num_eas)
    double precision, intent(in)    :: vel_(1:^NC)  !actually not used
    double precision, intent(in)    :: qdt

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1_unit)
    call fill_metric(metricM1_unit2)
    }

    {^D& ix_loc^D = ix^D \}
    qdt_impl = qdt
    alp = metricM1_unit%alp(ix^D)
    sqrtg = metricM1_unit%sqrtg(ix^D)
    eas = eas_

    Ut_explicit(1) = wrad(m1_energy_)
    {^C& Ut_explicit(1+^C) = wrad(m1_flux^C_) \}

    wrad(m1_energy_) = wrad(m1_energy_) / metricM1_unit%sqrtg(ix^D)
    {^C& wrad(m1_flux^C_) = wrad(m1_flux^C_) / metricM1_unit%sqrtg(ix^D) \}  

    stateM1%E = wrad(m1_energy_)
    {^C& stateM1%F_low(^C) = wrad(m1_flux^C_) \} 
    {^C& stateM1%vel(^C)   = vel_(^C) \} 

    call m1_update_closure_ixD(stateM1,metricM1_unit,ix^D,.false.)
    wrad(m1_energy_) = stateM1%E
    {^C& wrad(m1_flux^C_) = stateM1%F_low(^C) \}

    wrad(m1_energy_) = wrad(m1_energy_)*sqrtg
    {^C& wrad(m1_flux^C_) = wrad(m1_flux^C_)*sqrtg \}

  end subroutine fill_globals_metricM1_unit

  function get_random_start_unit(z_vec)
    double precision, dimension(:), intent(inout) :: z_vec
    double precision, dimension(1:size(z_vec)) :: get_random_start_unit
    integer :: i_z,n_z, clock
    integer, allocatable :: seed(:)
    double precision :: rand_num
    double precision :: z_low, z_high
       z_low = 0.0d0 + 1.d-40
       z_high = 1.0d0  

      ! Initialize the random seed
       call system_clock(count=clock)
       call random_seed(size=n_z)
       allocate(seed(n_z))
       seed = clock + 37 * (/ (i_z - 1, i_z = 1, n_z) /)
       call random_seed(put=seed)
       deallocate(seed)

      ! Generate and print 4 random numbers between 0 and 1
       do i_z = 1,m1_numvars_internal !size(z_vec)
          call random_number(rand_num)
          z_vec(i_z) = z_low + (z_high - z_low) * rand_num
          get_random_start_unit(i_z) = 1.0d0
       end do
      
  end function get_random_start_unit

  function m1_collisional_jacobian_unit(z_vec)
    use mod_m1_eas, only: Q_ems, k_a, k_s
    use mod_m1_internal
    use mod_m1_closure
    double precision, dimension(1:m1_system_ndim,1:m1_system_ndim) :: m1_collisional_jacobian_unit
    double precision, intent(in) :: z_vec(1:m1_system_ndim)
    !internal
    integer :: i,j
    double precision :: W2,W3
    double precision :: JvF, JfF, HvE, HfE, HdF, HvvF, HffF, HfvF, HvfF
    double precision :: dJdE
    double precision, dimension(^NC) :: dJdF, dHdE
    double precision, dimension(^NC,^NC) :: dHdF

    ! update the closure 
    call m1_update_closure_ixD(stateM1,metricM1_unit,ix_loc^D,.false.,get_vel_impl=.false.)

    W2 = stateM1%W * stateM1%W
    W3 = W2 * stateM1%W
      
    JvF = 2*stateM1%W**2*(-1+stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden - 2.0d0*stateM1%dthick*(W2-1.0d0)/(1.0d0 + 2.0d0*W2) )
    !JvF = 2*stateM1%W**2*(-1-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden + 2.0d0*stateM1%dthick*(W2-1.0d0)/(1.0d0 + 2.0d0*W2) )
    JfF = -2.0d0*stateM1%dthin*W2*stateM1%E*stateM1%fhatv**2/stateM1%Fden

    HvE = W3*(-1-stateM1%dthin*stateM1%fhatv**2 + stateM1%dthick*(2.0d0*W2-3.0d0)/(1.0d0+2.0d0*W2) )
    HfE = -stateM1%dthin*stateM1%W*stateM1%fhatv

    HdF = stateM1%W*(1-stateM1%dthick*stateM1%v2-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden)
    HvvF = 2.0d0*W3*(1.0d0-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden-stateM1%dthick*(1.0d0 -1.0d0/(2.0d0*W2*(1.0d0+2.0d0*W2))))
    !HvvF = 2.0d0*W3*(1.0d0-stateM1%dthin*stateM1%E*stateM1%fhatv/stateM1%Fden-stateM1%dthick*(stateM1%v2+1.0d0/(2.0d0*W2*(1.0d0+2.0d0*W2))))
    HffF = 2.0d0*stateM1%dthin*stateM1%W*stateM1%E*stateM1%fhatv/stateM1%Fden
    HvfF = 2.0d0*stateM1%dthin*W3*stateM1%E*stateM1%fhatv/stateM1%Fden
    !HvfF = 2.0d0*stateM1%dthin*W3*stateM1%E*stateM1%fhatv**2/stateM1%Fden
    HfvF = -stateM1%dthin*stateM1%W*stateM1%E*stateM1%fhatv/stateM1%Fden
    !HfvF = -stateM1%dthin*stateM1%W*stateM1%E/stateM1%Fden


    dJdE = W2+stateM1%dthin*stateM1%fhatv**2*W2+stateM1%dthick*((3.0d0-2.0d0*W2)*(W2-1.0d0))/(1.0d0+2.0d0*W2)
    dJdF(:) = JvF*stateM1%vel(:)  + JfF*stateM1%fhat_up(:)
    dHdE(:) = HvE*stateM1%vlow(:) + HfE*stateM1%fhat_low(:)

    do i=1, ^NC
       do j=1, ^NC
          dHdF(i,j) = HvvF*stateM1%vlow(i)*stateM1%vel(j) &
               + HffF * stateM1%fhat_low(i)*stateM1%fhat_up(j)  &
               + HvfF * stateM1%vlow(i) *stateM1%fhat_up(j)  &
               + HfvF * stateM1%fhat_low(i)*stateM1%vel(j)
          if( i .eq. j) dHdF(i,j) = dHdF(i,j) + HdF
       enddo
    enddo
    m1_collisional_jacobian_unit(1,1) = -1.d0* stateM1%W*( (eas(k_a)+eas(k_s)) - eas(k_s)*dJdE )
    do i=1, ^NC
       m1_collisional_jacobian_unit(1,i+1) = stateM1%W* ( dJdF(i) + (eas(k_a)+eas(k_s))*stateM1%vel(i) )
       m1_collisional_jacobian_unit(i+1,1) = -1.d0*( (eas(k_a)+eas(k_s))*dHdE(i) + stateM1%W*eas(k_a)*dJdE*stateM1%vlow(i) )
    enddo
    do i=1, ^NC
       do j=1, ^NC
          m1_collisional_jacobian_unit(i+1,j+1) = -1.d0*( (eas(k_a)+eas(k_s))*dHdF(i,j) + stateM1%W*eas(k_a)*stateM1%vlow(i)*dJdF(j) )
       enddo
    enddo
  end function m1_collisional_jacobian_unit


  function m1_collisional_sources_unit(z_vec)
    ! bar{G}
    use mod_m1_eas, only: Q_ems, k_a, k_s
    use mod_m1_internal
    use mod_m1_closure 
    double precision, dimension(1:m1_system_ndim) :: m1_collisional_sources_unit
    double precision, intent(in) :: z_vec(1:m1_system_ndim)
    ! update the closure 
    call m1_update_closure_ixD(stateM1,metricM1_unit,ix_loc^D,.false.,get_vel_impl=.false.)
    m1_collisional_sources_unit(1) = stateM1%W*( eas(Q_ems) + eas(k_s)*stateM1%J - (eas(k_s)+eas(k_a))*(stateM1%E-stateM1%Fv) )
    {^C& m1_collisional_sources_unit(1+^C) = stateM1%W*(eas(Q_ems)-eas(k_a)*stateM1%J)*stateM1%vlow(^C)&
         -(eas(k_a)+eas(k_s))*stateM1%H_low(^C) \}
  end function m1_collisional_sources_unit

  function m1_sources_func_unit(z_vec)
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec)) :: m1_sources_func_unit
    ! internal
    double precision, dimension(1:m1_system_ndim) :: coll_sources

    coll_sources = m1_collisional_sources_unit(z_vec) *alp * qdt_impl
    m1_sources_func_unit(:) = Ut_explicit(:) + sqrtg*( coll_sources(:) - z_vec(:)) 
  end function m1_sources_func_unit

  function m1_sources_jacobian_unit(z_vec)
    double precision, dimension(:), intent(in) :: z_vec
    double precision, dimension(1:size(z_vec),1:size(z_vec)) :: m1_sources_jacobian_unit
    ! internal
    integer :: i,j

    m1_sources_jacobian_unit(:,:) = m1_collisional_jacobian_unit(z_vec) * qdt_impl *alp
    m1_sources_jacobian_unit(1,1) = m1_sources_jacobian_unit(1,1) - 1
    {^C& m1_sources_jacobian_unit(1+^C,1+^C) = m1_sources_jacobian_unit(1+^C,1+^C) -1  \}
    m1_sources_jacobian_unit = m1_sources_jacobian_unit * sqrtg
    
  end function m1_sources_jacobian_unit
  
  function m1_sources_min_unit(z_vec)
    use mod_m1_metric_interface
    use mod_m1_tests
    double precision :: m1_sources_min_unit  
    double precision, dimension(:), intent(in) :: z_vec
    ! internal
    double precision, dimension(m1_system_ndim) :: Fsource_low, Fsource_hi

    Fsource_low = m1_sources_func_unit(z_vec)
    call raise_ixD(metricM1_unit,ix_loc^D,Fsource_low(2:m1_system_ndim),Fsource_hi(2:m1_system_ndim))

    m1_sources_min_unit = Fsource_low(1)*Fsource_low(1)
    {^C& m1_sources_min_unit = m1_sources_min_unit + Fsource_low(1+^C)*Fsource_hi(1+^C) \}
    m1_sources_min_unit = 0.5d0*m1_sources_min_unit
    
  end function m1_sources_min_unit
  }
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
end module mod_m1_collisional
