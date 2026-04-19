module mod_m1

  use mod_m1_closure 
  !use mod_m1_closure, disabled => m1_update_closure 

  implicit none

  public 

  contains
  {#IFNDEF UNIT_TESTS   
  !< Startup m1 modules, make sure routines are associated.
  subroutine m1_startup()
    use mod_m1_closure, only: m1_closure_activate
    !{#IFNDEF UNIT_TESTS_EXPLICIT
    {#IFNDEF M1_EXPLICIT
    use mod_m1_eas, only: m1_eas_activate, m1_eas_gray_activate
    use mod_m1_eas_test, only: m1_eas_test_activate

    use mod_Weakhub_reader
    }
    use mod_m1_eas_param
    !}
    {#IFDEF M1_EXPLICIT  
    include "amrvacdef.f"
    }
    integer :: spec_counter

    ! Initialize all indices to -10 (i.e. inactive)
    m1_i_nue = -10
    m1_i_nuebar = -10
    m1_i_nux = -10
    m1_i_mu = -10
    m1_i_mubar = -10
    m1_i_photon = -10
    
    spec_counter = 1
    
    if (m1_use_neutrinos) then
        if (mype==0) write(*,*) "M1: You will use M1 with neutrinos!"
        m1_i_nue    = spec_counter
        spec_counter = spec_counter + 1
        m1_i_nuebar = spec_counter
        spec_counter = spec_counter + 1
        m1_i_nux   = spec_counter
        spec_counter = spec_counter + 1
    else
        if (mype==0) write(*,*) "M1: You will NOT use M1 with neutrinos!"
    end if
    
    if (m1_use_muons) then
        if (mype==0) write(*,*) "M1: You will use M1 with muons!"
        m1_i_mu    = spec_counter
        spec_counter = spec_counter + 1
        m1_i_mubar = spec_counter
        spec_counter = spec_counter + 1
    else
        if (mype==0) write(*,*) "M1: You will NOT use M1 with muons!"
    end if
    
    if (m1_use_photons) then
        if (mype==0) write(*,*) "M1: You will use M1 with photons!"
        m1_i_photon = spec_counter
    else
        if (mype==0) write(*,*) "M1: You will NOT use M1 with photons!"
    end if


    if(mype==0) write(*,*) "M1: The indices of radiation in species loop are:"   
    if(mype==0) write(*,*) "M1: m1_i_nue    ",m1_i_nue
    if(mype==0) write(*,*) "M1: m1_i_nuebar ",m1_i_nuebar
    if(mype==0) write(*,*) "M1: m1_i_nux   ",m1_i_nux
    if(mype==0) write(*,*) "M1: m1_i_mu     ",m1_i_mu
    if(mype==0) write(*,*) "M1: m1_i_mubar  ",m1_i_mubar
    if(mype==0) write(*,*) "M1: m1_i_photon ",m1_i_photon
    
    call m1_closure_activate()
    ! m1_eas_tests:
    {#IFDEF M1_RATES_TEST1
    call m1_eas_tests_activate(1)
    ! activates the eas if m1_get_eas is associated to test
    call m1_eas_activate()
    }
    {#IFDEF M1_RATES_TEST2
    call m1_eas_tests_activate(2)
    ! activates the eas if m1_get_eas is associated to test
    call m1_eas_activate()
    }
    {#IFDEF M1_RATES 
    ! standard gray microphysics
      call m1_eas_gray_activate()
      call m1_read_Weakhub()  
      if(mype==0)  write(*,*)"M1: Actived tabulated rates from Weakhub"
    ! activates the eas if m1_get_eas is associated to gray rates
      call m1_eas_activate()
    }
    if(mype==0) write(*,*) "Rise and shine, M1 is ready!"

  end subroutine m1_startup
  }

  !-----------------------------------------------
  subroutine m1_read_Weakhub()
   use mod_Weakhub_reader
   {#IFNDEF UNIT_TESTS
   !include "amrvacdef.f"
   }
   {#IFDEF UNIT_TESTS
   character(len=512) :: fileWeakhub
   fileWeakhub = "/mnt/rafast/mcassing/Weakhub/dURCA_convenforbetaproc_DD2_grey_npe_20240723_Weakhub_grey.h5"
   }

   call read_h5_Weakhub_table_gray(fileWeakhub)

  end subroutine m1_read_Weakhub
   !-------------------------------------------

  !< Compute lambda_p/m for M1 variables. This routine is a level 1 
  !< M1 routine and is called inside of get_lambda (amrvacphys.t)
  subroutine m1_get_wavespeeds(ixI^L, ixO^L, idim, wprim, x, lambda, qtC)
    ! This routine assumes the  !
    ! primitives are up to date !
    ! ====================================================================== 
    use mod_m1_closure
	  use mod_m1_metric_interface
    ! ====================================================================== 
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    integer, intent(in) :: ixI^L, ixO^L, idim
    double precision, intent(inout) :: wprim(ixI^S,1:nw) 
	  double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: qtC
	  double precision, intent(out)   :: lambda(ixI^S, 1:ncons, 1:2)
    ! ====================================================================== 
    !internal variables
    integer :: ix^D,i,idir
    double precision, dimension(ixI^S) :: dthick, dthin
    double precision, dimension(ixI^S) :: lmthin, lmthick, lpthin, lpthick
    double precision, dimension(ixI^S) :: p_tmp, r_tmp, chi, lfact
    double precision, dimension(ixI^S) :: F_sq
    double precision, dimension(ixI^S,^NC) :: F_hi, vel
    integer :: metric_comp 
	  type(m1_metric_helper) :: metricM1  
    double precision, dimension(ixI^S,^NC) :: testa, testb
    double precision :: M1_mini_pl, M1_mini_min
    {#IFDEF UNITS_TESTS
    double precision :: TESTqtC = 1.0d+10
    }

    ! ====================================================================== 
    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wprim,x,ixI^L,ixO^L) 
    }

    select case(idim)
    case (1)
      metric_comp = 1 
    case (2)
      metric_comp = 4 
    case (3)
      metric_comp = 6 
    end select 

    
    !###############################################################################################
    !###############################################################################################
	  {^KSP& ! Species loop 
    ! Remove normalisation from radiation moments (temporarily)
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) / metricM1%sqrtg(ixO^S)
    wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) / metricM1%sqrtg(ixO^S)
    {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) / metricM1%sqrtg(ixO^S) \}  

    if(m1_actual_speeds) then
        ! in case singularity and sqrtg = 0:
        {^D& do ix^D=ixOmin^D,ixOmax^D \}
        if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
          wprim(ix^D,nrad^KSP_) = 0.0d0
          wprim(ix^D,erad^KSP_) = 0.0d0
          {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
        end if 
        {^D& end do \}

        ! update closure 
        ! This routine in general will have reconstructed 
        ! cell face values as input, so the closure always needs to 
        ! be updated.
        call m1_update_closure(metricM1, wprim, x, ixI^L, ixO^L, ^KSP, .true., chi=chi, W=lfact, vel=vel)
    
        ! Compute F_sq 
        call metricM1%raise(ixI^L,ixO^L,wprim(ixI^S,frad^KSP1_:frad^KSP^NC_), F_hi)

        F_sq(ixO^S)=0.0d0
        {^C& 
        F_sq(ixO^S) =F_sq(ixO^S) + wprim(ixO^S,frad^KSP^C_)*F_hi(ixO^S,^C)
        \}
    
        ! ====================================================================== 
        ! lambda_m^(thin)
        lmthin(ixO^S) = -metricM1%beta(ixO^S,idim) &
        - metricM1%alp(ixO^S)*dabs(F_hi(ixO^S,idim))/dsqrt(F_sq(ixO^S)+1.d-30)
        ! lambda_p^(thin)
        lpthin(ixO^S) = -metricM1%beta(ixO^S,idim) &
        + metricM1%alp(ixO^S)*dabs(F_hi(ixO^S,idim))/dsqrt(F_sq(ixO^S)+1.d-30)
        ! Auxiliaries for thick lambdas 
        p_tmp(ixO^S) = metricM1%alp(ixO^S) * vel(ixO^S,idim)/lfact(ixO^S)
        r_tmp(ixO^S) = dsqrt( metricM1%alp(ixO^S)*metricM1%alp(ixO^S) * metricM1%gammaUPij(ixO^S,metric_comp) * &
             (2.0d0*lfact(ixO^S)*lfact(ixO^S) + 1.0d0) - 2*lfact(ixO^S)*lfact(ixO^S)*p_tmp(ixO^S)*p_tmp(ixO^S))  
        ! lambda_m^(thick)
        lmthick(ixO^S) = MIN(-metricM1%beta(ixO^S,idim)+p_tmp(ixO^S), &
             -metricM1%beta(ixO^S,idim)+ (2.d0*p_tmp(ixO^S)*lfact(ixO^S)*lfact(ixO^S) - r_tmp(ixO^S))/(2.0d0*lfact(ixO^S)*lfact(ixO^S)+1.0d0) )
        ! lambda_p^(thick) 
        lpthick(ixO^S) = MAX(-metricM1%beta(ixO^S,idim)+p_tmp(ixO^S), &
             -metricM1%beta(ixO^S,idim)+ (2.0d0*p_tmp(ixO^S)*lfact(ixO^S)*lfact(ixO^S) + r_tmp(ixO^S))/(2.0d0*lfact(ixO^S)*lfact(ixO^S)+1.0d0) )
        ! ====================================================================== 
        dthin(ixO^S) = 1.5d0*chi(ixO^S) - 0.5d0
        dthick(ixO^S) = 1.0d0-dthin(ixO^S)

        M1_mini_pl = 1.0d-30
        M1_mini_min = 1.0d-30
    
        ! Lambda for n, F are identical! 
        ! _____________ Store lambda + ____________________
        lambda(ixO^S,erad^KSP_,1) = dthick(ixO^S)*lpthick(ixO^S) + dthin(ixO^S)*lpthin(ixO^S) + M1_mini_pl
        lambda(ixO^S,nrad^KSP_,1) = dthick(ixO^S)*lpthick(ixO^S) + dthin(ixO^S)*lpthin(ixO^S)+ M1_mini_pl
        {^C& lambda(ixO^S,frad^KSP^C_,1) = dthick(ixO^S)*lpthick(ixO^S) + dthin(ixO^S)*lpthin(ixO^S) + M1_mini_pl \} 
        ! ______________ Store lambda - ______________________
        lambda(ixO^S,erad^KSP_,2) = dthick(ixO^S)*lmthick(ixO^S) + dthin(ixO^S)*lmthin(ixO^S) + M1_mini_min
        lambda(ixO^S,nrad^KSP_,2) = dthick(ixO^S)*lmthick(ixO^S) + dthin(ixO^S)*lmthin(ixO^S) + M1_mini_min
        {^C& lambda(ixO^S,frad^KSP^C_,2) = dthick(ixO^S)*lmthick(ixO^S) + dthin(ixO^S)*lmthin(ixO^S) + M1_mini_min \} 
        ! ======================================================================
    !++++++++++++++ Speeds as in Radice et al +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if(m1_radice_speeds) then
        ! _____________ Store lambda + ____________________
        lambda(ixO^S,erad^KSP_,1) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) + metricM1%beta(ixO^S,1)) + M1_mini_pl
        lambda(ixO^S,nrad^KSP_,1) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) + metricM1%beta(ixO^S,1)) + M1_mini_pl
        {^C& lambda(ixO^S,frad^KSP^C_,1) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) + metricM1%beta(ixO^S,1)) + M1_mini_pl \} 
        ! ______________ Store lambda - ______________________
        lambda(ixO^S,erad^KSP_,2) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) - metricM1%beta(ixO^S,1)) + M1_mini_min
        lambda(ixO^S,nrad^KSP_,2) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) - metricM1%beta(ixO^S,1)) + M1_mini_min
        {^C& lambda(ixO^S,frad^KSP^C_,2) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) - metricM1%beta(ixO^S,1)) + M1_mini_min\} 
        ! ====================================================================== 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else 
      lambda(ix^D,nrad^KSP_,1) = 1.0d0
      lambda(ix^D,erad^KSP_,1) = 1.0d0
      {^C& lambda(ix^D,frad^KSP^C_,1) = 1.0d0 \}
      lambda(ix^D,nrad^KSP_,2) = - 1.0d0
      lambda(ix^D,erad^KSP_,2) = - 1.0d0
      {^C& lambda(ix^D,frad^KSP^C_,2) = - 1.0d0 \}
    end if 
    !++++++++++++++ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    ! rescale back rad. variables to evolved ones
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) * metricM1%sqrtg(ixO^S) 
    wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) * metricM1%sqrtg(ixO^S) 
    {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) *metricM1%sqrtg(ixO^S) \}

    ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
     if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
       wprim(ix^D,nrad^KSP_) = 0.0d0
       wprim(ix^D,erad^KSP_) = 0.0d0
       {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
       lambda(ix^D,nrad^KSP_,1) = 0.577350269189625d0
       lambda(ix^D,erad^KSP_,1) = 0.577350269189625d0
       {^C& lambda(ix^D,frad^KSP^C_,1) = 0.577350269189625d0 \}
       lambda(ix^D,nrad^KSP_,2) = - 0.577350269189625d0
       lambda(ix^D,erad^KSP_,2) = - 0.577350269189625d0
       {^C& lambda(ix^D,frad^KSP^C_,2) = - 0.577350269189625d0 \}
     end if 
    {^D& end do \}

    ! in case wavespeeds too close to zero:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
     if(lambda(ix^D,nrad^KSP_,1) .le. 1.d-4 .and. lambda(ix^D,nrad^KSP_,2) .ge. -1.d-4) then
       lambda(ix^D,nrad^KSP_,1) = 1.0d0
       lambda(ix^D,erad^KSP_,1) = 1.0d0
       {^C& lambda(ix^D,frad^KSP^C_,1) = 1.0d0 \}
       lambda(ix^D,nrad^KSP_,2) = - 1.0d0
       lambda(ix^D,erad^KSP_,2) = - 1.0d0
       {^C& lambda(ix^D,frad^KSP^C_,2) = - 1.0d0 \}
     end if 
    {^D& end do \}

    \} ! end ^ KSP
    ! ====================================================================== 
    !###############################################################################################
    !###############################################################################################
      
    call metricM1%destroy()

  end subroutine m1_get_wavespeeds
  
  subroutine m1_get_fluxes(ixI^L, ixO^L, idim, x, wprim, flux, qtC)
    use mod_m1_closure
	  use mod_m1_metric_interface
    ! ====================================================================== 
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    ! ====================================================================== 
    integer, intent(in) :: ixI^L, ixO^L, idim
	  double precision, intent(in)  :: qtC
	  double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout)  :: wprim(ixI^S,1:nw) 
    double precision, intent(out) :: flux(ixI^S,1:nwflux)
    ! ====================================================================== 
    ! internal variables
    ! ====================================================================== 
    integer iw, ix^D, i, idir
    double precision, dimension(ixI^S,^NC) :: F, Hup, vel
    double precision, dimension(ixI^S,m1_npress) :: P
    double precision, dimension(ixI^S)    :: Gamma, Jradi, E, Wlor
    double precision, dimension(ixI^S,^NC) :: Pimu
    double precision, dimension(ixI^S,^NC) ::  Pud
	  double precision :: zetaL
    {#IFDEF UNITS_TESTS
    double precision :: TESTqtC = 1.0d+10
    }
	  logical :: get_zeta
    logical :: check_atmo_in_getflux = .false.
    !integer, dimension(3,3) :: idimtoimu  
    ! ====================================================================== 
	  type(m1_metric_helper)   :: metricM1 
    ! ====================================================================== 

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wprim,x,ixI^L,ixO^L)
    }

    {^KSP& ! Species loop 
    ! Remove normalisation 
    ! IMPORTANT: always needs to be added back at the end of 
    ! the routine. 
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) / metricM1%sqrtg(ixO^S)
    wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) / metricM1%sqrtg(ixO^S)
    {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) / metricM1%sqrtg(ixO^S) \}   
    
     ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}  
    ! update closure to get E,Fi,zeta
    call m1_update_closure(metricM1, wprim, x, ixI^L, ixO^L,^KSP,.true.,Gamma,Hup,Jradi,P,W=Wlor,vel=vel)

    ! Remove Gamma from n 
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_)/Gamma(ixO^S) 
    ! get energy
    E(ixO^S)  = wprim(ixO^S,erad^KSP_)
    ! Pimu = ( PUPxx, PUPxy, PUPxz ) if idim = 1 
    ! Pimu = ( PUPxy, PUPyy, PUPyz ) if idim = 2
    ! Pimu = ( PUPxz, PUPyz, PUPzz ) if idim = 3 
    {^C& Pimu(ixO^S,^C) = P(ixO^S,metricM1%idimtoimu(idim,^C) ) \} 

    ! Shuffle some indices up and down cause why not 
    call metricM1%lower(ixI^L,ixO^L,Pimu, Pud) !ixO
		call metricM1%raise(ixI^L,ixO^L,wprim(ixI^S,frad^KSP1_:frad^KSP^NC_), F) 

    ! Get fluxes 
    ! ====================================================================== 
    ! Nflux^idim = sqrtg * alp * n * (( W*v^idim-beta^idim/alp) + H^idim/Jrad )
    ! See Radice, Bernuzzi et al. Eqs. (23) and (28)
    flux(ixO^S,nrad^KSP_) = metricM1%sqrtg(ixO^S) * metricM1%alp(ixO^S) * wprim(ixO^S,nrad^KSP_) &
         * (Wlor(ixO^S)*( vel(ixO^S,idim) - metricM1%beta(ixO^S,idim)/metricM1%alp(ixO^S)) &
         + Hup(ixO^S,idim)/(Jradi(ixO^S) + M1_TINY))  
    ! ====================================================================== 
    ! Eflux^idim  = sqrtg * ( alp * F^idim - beta^idim* E )
     flux(ixO^S,erad^KSP_) = metricM1%sqrtg(ixO^S) * ( metricM1%alp(ixO^S) * F(ixO^S,idim) &
    - metricM1%beta(ixO^S,idim) * E(ixO^S) )
    !flux(ixO^S,erad^KSP_) = metricM1%sqrtg(ixO^S) * metricM1%alp(ixO^S) *( F(ixO^S,idim) &
    !- metricM1%beta(ixO^S,idim) * E(ixO^S)/metricM1%alp(ixO^S)  )
    !!!TEST  
    ! ====================================================================== 
    ! Fflux^idim_i  = sqrtg * ( alp * P^idim_i - beta^idim * F_i )
    {^C&flux(ixO^S,frad^KSP^C_) = metricM1%sqrtg(ixO^S) * ( metricM1%alp(ixO^S) * Pud(ixO^S,^C) &
    - metricM1%beta(ixO^S,idim) * wprim(ixO^S,frad^KSP^C_) ) \}
    ! ====================================================================== 
    ! ====================================================================== 
    ! ====================================================================== 
    ! NB here is where we re-densitize the radiation vars 
 		 wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_)*metricM1%sqrtg(ixO^S)*Gamma(ixO^S)
		 wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_)*metricM1%sqrtg(ixO^S)
		 {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_)*metricM1%sqrtg(ixO^S) \}     
     
     ! in case singularity and sqrtg = 0:
     {^D& do ix^D=ixOmin^D,ixOmax^D \}
     if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
     flux(ix^D,nrad^KSP_) = 0.0d0
     flux(ix^D,erad^KSP_) = 0.0d0
     {^C& flux(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}

    !---------------------------------------
    ! is set to false:
    ! check atmosphere
    if(check_atmo_in_getflux) then
       {^D& do ix^D=ixImin^D,ixImax^D \}
       ! check wprim
       if(wprim(ix^D,erad^KSP_) .lt. m1_E_atmo) then
         wprim(ix^D,nrad^KSP_) = m1_E_atmo
         wprim(ix^D,erad^KSP_) = m1_E_atmo
       end if 
       {^C&
       if((wprim(ix^D,frad^KSP^C_) .lt. m1_E_atmo) .and. (wprim(ix^D,frad^KSP^C_) .gt. -1.0d0*m1_E_atmo)) then
          wprim(ix^D,frad^KSP^C_) = m1_E_atmo 
          if(wprim(ix^D,frad^KSP^C_) .gt. -1.0d0*m1_E_atmo) wprim(ix^D,frad^KSP^C_) = -1.0d0*m1_E_atmo
       end if 
       \}
       ! check flux
       if(flux(ix^D,erad^KSP_) .lt. m1_E_atmo) then
         flux(ix^D,nrad^KSP_) = m1_E_atmo
         flux(ix^D,erad^KSP_) = m1_E_atmo
       end if 
       {^C&
       if((flux(ix^D,frad^KSP^C_) .lt. m1_E_atmo) .and. (flux(ix^D,frad^KSP^C_) .gt. -1.0d0*m1_E_atmo)) then
          flux(ix^D,frad^KSP^C_) = m1_E_atmo 
       end if 
       \}
       {^D& end do \}
    endif
    !---------------------------------------

    \} ! End of species loop 
    call metricM1%destroy()
    ! ====================================================================== 
    ! ====================================================================== 
    !                               All done! 
    ! ====================================================================== 
    ! ====================================================================== 
  end subroutine m1_get_fluxes

  !< Apply correction to hll flux 
  !< to preserve asymptotic (diffusion)
  !< limit of RTE
  subroutine m1_correct_asymptotic_fluxes(ix^D, ixI^L, tvdlfeps1, dxdim, idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)
  {#IFDEF UNIT_TESTS
  use mod_m1_tests
  }
  {#IFNDEF UNIT_TESTS
  include "amrvacdef.f"
  }
  integer, intent(in) :: ix^D, ixI^L, idims  
  double precision, intent(in) :: tvdlfeps1, qtC
  double precision, intent(in), dimension(1:ndim)                   :: dxdim        !< grid spacing 
  double precision, intent(in), dimension(ixI^S,1:nwflux)           :: fLC, fRC     !< left/right fluxes 
  double precision, intent(inout), dimension(ixI^S,1:nwflux,1:ndim) :: fC           !< flux 
  double precision, intent(in), dimension(ixI^S,1:nw)               :: wLC, wRC     !< left/right state 
  double precision, intent(in), dimension(ixI^S,1:ncons)            :: cminC, cmaxC !< min / max wavespeeds 
  double precision, intent(in), dimension(ixI^S,1:nm1rad_eas)           :: wradimpl     !old state 

  ! internals 
  double precision :: A 
  {#IFDEF UNITS_TESTS
  double precision :: TESTqtC = 1.0d+10
  }
  double precision :: fc_loc
  integer          :: iw, ixR^D, i,j

  {#IFDEF UNIT_TESTS
  !Kronecker delta and Levi-Civita tensors
  INTEGER:: kr(0:3,0:3),lvc(1:3,1:3,1:3)
  }

  {#IFDEF UNIT_TESTS    
    ! Kronecker delta defined in mod_m1_tests    
   do i= 0,3
    do j=0,3
       kr(i,j)=0            
       if(i.eq.j) then
        kr(i,j)=1
       end if
    end do 
   end do 
  }
  {^D& ixR^D = ix^D + kr(idims,^D)\}
  
  {^KSP& 
    ! Compute A factor 
    A = MIN(1.0d0, dabs(4.0d0/ dabs(dxdim(idims)) /(wradimpl(ixR^D,kappa_a^KSP_) + wradimpl(ix^D,kappa_a^KSP_) &
                        + wradimpl(ixR^D,kappa_s^KSP_) + wradimpl(ix^D,kappa_s^KSP_) + 1.0d-40)))

      if(A .ne. A) then
        A = 1.0d0
      end if 
      if(A .gt. 1.0d0) then 
        A = 1.0d0
      end if
                        
    ! Correct energy and number fluxes
    ! F(erad) ==> F(erad) + (A-1)/(cmax-cmin) * cmax*cmin * (UR-UL)

    !-------------------- erad
    ! overwrite HLL like this else in fc already is cmax*cmin/(cmax-xmin)*(WR-wL) without A
    !------------------------                    
    iw = erad^KSP_    
    if(m1_erad_LLF) then
      call mpistop("m1 LLF not implemented yet")
    else if(m1_parabolic) then
    ! to test diffusion limit and equations get parabolic
    fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw))&
                  /(cmaxC(ix^D,iw)-cminC(ix^D,iw))    
    else
      fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw)&
      +A*tvdlfeps1*cminC(ix^D,iw)&
      *cmaxC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
      /(cmaxC(ix^D,iw)-cminC(ix^D,iw)) 
    end if 
    !-------------------- nrad
    ! overwrite HLL like this else in fc already is cmax*cmin/(cmax-xmin)*(WR-wL) without A
    !------------------------ 
    iw = nrad^KSP_ 
    if(m1_parabolic) then
      ! test diffusion limit:
      fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw))&
                    /(cmaxC(ix^D,iw)-cminC(ix^D,iw))
    else 
      fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw)&
                  +A*tvdlfeps1*cminC(ix^D,iw)&
                  *cmaxC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
                  /(cmaxC(ix^D,iw)-cminC(ix^D,iw))
    end if 

     !-------------------- frad-C
    ! overwrite HLL like this else in fc already is cmax*cmin/(cmax-xmin)*(WR-wL) without A
    !------------------------   
    {^C&
    iw = frad^KSP^C_ 
    if(m1_frad_A2_A3) then
       fC(ix^D, iw,idims) = A**2*(cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw)&
       +A*tvdlfeps1*cminC(ix^D,iw)&
       *cmaxC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
       /(cmaxC(ix^D,iw)-cminC(ix^D,iw))
   
       fC(ix^D,iw,idims) = fC(ix^D, iw,idims) +0.5d0 *(1.0d0-A**2)* ( fLC(ix^D,iw) + fRC(ix^D,iw) )
    else if(m1_frad_A2_A) then
      fC(ix^D, iw,idims) = (A**2*(cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw))&
                          +A*tvdlfeps1*cmaxC(ix^D,iw)*cminC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))/&
                          (cmaxC(ix^D,iw)-cminC(ix^D,iw)) + 0.5d0*(1.0d0-A**2)*(fLC(ix^D,iw)+fRC(ix^D,iw))
    else if(m1_frad_LLF) then
      call mpistop("LLF for frad not implemented yet")
    else if(m1_parabolic) then
       !diffusive lim
       fC(ix^D, iw,idims) =  0.5d0 * ( fLC(ix^D,iw) + fRC(ix^D,iw) )
    else 
      fC(ix^D, iw,idims) = A**2 * fC(ix^D, iw,idims) + 0.5d0 * (1.0d0 - A**2) * ( fLC(ix^D,iw) + fRC(ix^D,iw) )
    end if 
    \}

    !---------------------------------------
    ! check atmosphere
    ! check flux
    !!if(fC(ix^D,erad^KSP_,idims) .lt. m1_E_atmo) then
    !!  fC(ix^D,nrad^KSP_,idims) = 0.0d0 !m1_E_atmo
    !!  fC(ix^D,erad^KSP_,idims) = 0.0d0 !m1_E_atmo
    !!end if 
    !!{C
    !!if((fC(ix^D,frad^KSP^C_,idims) .lt. m1_E_atmo) .and. (fC(ix^D,frad^KSP^C_,idims) .gt. -1.0d0*m1_E_atmo)) then
    !!  fC(ix^D,frad^KSP^C_,idims) = 0.0d0 !m1_E_atmo 
    !!end if 
    !!}
    !------------------
    ! check for NAN :
    !{#IFDEF NAN_CHECK
    if(fC(ix^D,nrad^KSP_,idims) .ne. fC(ix^D,nrad^KSP_,idims)) then
      write(777,*) "---hll: NAN detected: FC(nrad)", fC(ix^D,nrad^KSP_,idims)
      fC(ix^D,nrad^KSP_,idims) = 0.0d0
    end if 
    if(fC(ix^D,erad^KSP_,idims) .ne. fC(ix^D,erad^KSP_,idims)) then
      write(777,*) "---hll: NAN detected: FC(erad)", fC(ix^D,erad^KSP_,idims)
      fC(ix^D,erad^KSP_,idims) = 0.0d0
    end if 
    {^C&
    if(fC(ix^D,frad^KSP^C_,idims) .ne. fC(ix^D,frad^KSP^C_,idims))then
      write(777,*) "---hll: NAN detected: FC(frad)", fC(ix^D,frad^KSP^C_,idims)
      fC(ix^D,frad^KSP^C_,idims) = 0.0d0
    end if
    \}
    !}
    !---------------------------------------

  \}   ! end KSP loop

  end subroutine m1_correct_asymptotic_fluxes

  ! ====================================================================== 
  ! ====================================================================== 
  !< compute geometrical sources for M1
  subroutine m1_add_geometrical_sources(ixI^L,ixO^L,x,wprim,wcons, qdt,d_k_gamma_ij,Kij_local,qtC)
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    use mod_m1_closure 
	  use mod_m1_metric_interface
    {#IFNDEF UNIT_TESTS
    use mod_metric
    include "amrvacdef.f"
    }
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, qtC
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
	  double precision, intent(inout)    :: Kij_local(ixI^S,^NC,^NC)
	  double precision, intent(inout)    :: d_k_gamma_ij(ixI^S,^NC,^NC,^NC)
    double precision, intent(inout) :: wcons(ixI^S, 1:nw) 
	  double precision, intent(inout) :: wprim(ixI^S, 1:nw)
    ! internal
    integer :: inonzero, i, j, pressind, ix^D
    double precision :: symfact
    double precision, dimension(ixI^S) :: rhs
    double precision, dimension(ixI^S,m1_npress) :: P
    double precision, dimension(ixI^S,^NC) :: Pimu
    double precision, dimension(ixI^S,^NC) :: FU
	  type(m1_metric_helper) :: metricM1 
    {#IFDEF UNITS_TESTS
    double precision :: TESTqtC = 1.0d+10
    }
    
    !===============================

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wprim,x,ixI^L,ixO^L)
    }
    
	{^KSP&
     rhs(ixO^S) = 0.0d0
     
     ! multiply with sqrtg at end of add_geometrical_sources
 		 wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) &
     / metricM1%sqrtg(ixO^S)
		 wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) / metricM1%sqrtg(ixO^S)
		 {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) &
     / metricM1%sqrtg(ixO^S) \}   

    ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}


    ! update closure to get E,Fi,zeta and pressure
    call m1_update_closure(metricM1, wprim, x, ixI^L, ixO^L,^KSP,.true.,Press=P) 

    !------------  N-density -------------
	  ! geom-source is zero:
	  ! wcons(ixO^S,nrad^KSP_) =  wcons(ixO^S,nrad^KSP_) +qdt*0.0d0*metricM1%sqrtg(ixO^S)

    ! ------------ Energy ------------------------------
    ! -- 
    ! Source term is:
    ! sqrtg [ alp PK - F^i partial_i alp ]
    ! --
      {#IFDEF M1_FLAT_KIJ
      Kij_local(ixO^S,:,:) = 0.0d0
      }
      {#IFDEF M1_FLAT_DKDGIJ
       d_k_gamma_ij(ixO^S,:,:,:)=0.0d0
      }

    rhs(ixO^S) = Kij_local(ixO^S,1,1)*P(ixO^S,1) &
         + Kij_local(ixO^S,1,2)*P(ixO^S,2) &
         + Kij_local(ixO^S,1,3)*P(ixO^S,3) &
         + Kij_local(ixO^S,2,2)*P(ixO^S,4) &
         + Kij_local(ixO^S,2,3)*P(ixO^S,5) &
         + Kij_local(ixO^S,3,3)*P(ixO^S,6) & 
         + Kij_local(ixO^S,2,1)*P(ixO^S,2) & 
         + Kij_local(ixO^S,3,1)*P(ixO^S,3) & 
         + Kij_local(ixO^S,3,2)*P(ixO^S,5) 
   
    rhs(ixO^S) = rhs(ixO^S) * metricM1%alp(ixO^S)

    call metricM1%raise(ixI^L,ixO^L,wprim(ixI^S,frad^KSP1_:frad^KSP^NC_),FU)   

	  do i = 1,^NC 
       rhs(ixO^S) = rhs(ixO^S) - FU(ixO^S,i)* metricM1%dalp(ixO^S,i)  
    end do 
	
    wcons(ixO^S,erad^KSP_) = wcons(ixO^S,erad^KSP_) + qdt*rhs(ixO^S)*metricM1%sqrtg(ixO^S)

    ! ------------ Fluxes ----------------------------
    ! --
    ! Source term is:
    ! sqrtgamma [ F_i partial_k beta^i -
    ! E partial_k alp + alp/2*P^{ij} partial_k gamma_{ij} ]
    ! --
    
    {^C&
    ! first set all to zero if all metric derivs are zero 
    rhs(ixO^S) = 0.0d0
    ! - E partial alp
       rhs(ixO^S) = -1.0d0 * wprim(ixO^S,erad^KSP_) * metricM1%dalp(ixO^S,^C)

    ! F_i partial_k beta^i
    do i=1, ^NC
       rhs(ixO^S) = rhs(ixO^S) + wprim(ixO^S,frad^KSP1_+i-1) * metricM1%dbeta(ixO^S,i,^C) 
    end do

    pressind = 1
    symfact  = 1.0d0
    do i= 1, ^NC
       do j=i, ^NC
          if ( i .ne. j) then
             symfact = 2.0d0
          else
             symfact = 1.0d0
          end if
          rhs(ixO^S) = rhs(ixO^S) + 0.5d0 * metricM1%alp(ixO^S) * symfact &
          * P(ixO^S,pressind) * d_k_gamma_ij(ixO^S,i,j,^C)
          pressind = pressind + 1 
       end do
    end do
    wcons(ixO^S,frad^KSP^C_) = wcons(ixO^S,frad^KSP^C_) + qdt*rhs(ixO^S)*metricM1%sqrtg(ixO^S)
    \}   
    
  !-------------------------
     ! rescale back wprim with sqrtg 
 		 wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_)*metricM1%sqrtg(ixO^S)
		 wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_)*metricM1%sqrtg(ixO^S)
		 {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) *metricM1%sqrtg(ixO^S) \}   
     
    ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}

    \} ! end ^KSP

    ! All done
    call metricM1%destroy()

    !===============================

  end subroutine m1_add_geometrical_sources

  !****************************************************************************

  {#IFNDEF UNIT_TESTS_EXPLICIT
  {#IFNDEF M1_EXPLICIT
  subroutine m1_add_collisional_sources(dtfactor,qdt,qtC,x,wcons,wold,wradimpl,ixI^L)
    ! wcons = psa%w%w out, at return it will be =U+F+S+G
    ! wold = psb%w%w old state without implicit part, at return = U+F+S
    ! wradimpl = psm1%%pwrad%wradimpl

    use mod_m1_internal
    use mod_m1_collisional, only: m1_get_implicit_collisional_sources
    use mod_m1_eas
    use mod_eos, only: small_rho, small_temp, big_ye, eos_yemin 
    use mod_Weakhub_reader, only: logtemp_min_IV, logtemp_max_IV,ye_min_IV, ye_max_IV,logrho_min_IV,logrho_max_IV
    use mod_m1_eas_param
    {#IFNDEF UNIT_TESTS
    use mod_m1_metric_interface
    include "amrvacdef.f"
    }
    integer, intent(in) :: ixI^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(in)    :: dtfactor   !< Timestep factor 
    double precision, intent(inout) :: wcons(ixI^S,1:nw)
    double precision, intent(inout) :: wold(ixI^S,1:nw)
    double precision, intent(in) :: wradimpl(ixI^S,1:nm1rad_eas)
    double precision, intent(in) :: qdt,qtC
    !internal
    integer :: ix^D, ixO^L
    double precision :: N
    double precision :: eas_ixD(1:m1_num_eas)
    double precision, dimension(1:m1_numvars_internal) :: wrad
    double precision, dimension(1:m1_numvars_internal) :: wradold
    double precision, dimension(1:m1_numvars_internal+1,^NS) :: sources
    double precision, dimension(1:^NC) :: vel_
    double precision :: t_backreact
    !double precision :: couple_fluid = 0.02d0
    double precision :: couple_tau = 1.0d-3
    double precision :: couple_momenta = 1.0d-3
    double precision :: couple_ye = 1.0d-4
    double precision :: couple_ymu = 1.0d-4
    double precision :: dE_source, dS_source1, dS_source2, dS_source3, dS_norm, dYeC_source, dYmuC_source
    double precision:: tiny, Snorm, fracE, fracS, rel_dye_den, fracYe_rel, absYe_cap, fracYe_abs
    double precision :: rel_dymu_den,fracYmu_rel, absYmu_cap, fracYmu_abs
    double precision:: capE, capS, capYe, capYmu, factorA
    double precision :: deltaErad, deltaN
    double precision :: Er,Fmag,facFrad
    double precision :: fmaxfact = 0.999d0 !0.99999d0 !0.999d0
    double precision :: deltaFrad(1:^NC)
    double precision :: tau_tmp, ye_tmp, dye_tmp, error_Ye, Delta_Ye
    double precision :: rel_error_Ye = 1.0d-3
    double precision :: damping
    double precision :: stepsBetweenDamp = 10.d0 !5.0d0 
    logical :: m1_apply_limits = .true. !.false. !whether to limit source changes and limit wrad
    logical :: m1_Iterative_Damping = .false. !.true.
    logical :: energy_good
    type(m1_metric_helper) :: metricM1 	
    integer, dimension(1:6) :: signum
    integer :: i

    ixO^L=ixI^L^LSUBdixB;

    t_backreact = m1_tset_backreact !2.0d0 * m1_tset

	  !need wprim here for metric!
    {#IFDEF UNITT_TESTS2
     call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wold,x,ixI^L,ixO^L)  
    }
    
    if(m1_use_neutrinos) then
      signum(m1_i_nue) = 1
      signum(m1_i_nuebar) = -1
      signum(m1_i_nux) = 0
    end if 
    if(m1_use_muons) then
      signum(m1_i_mu) = 1 !M1_muons_TODO
      signum(m1_i_mubar) = -1  !M1_muons_TODO
    end if 
    if(m1_use_photons)then
      signum(m1_i_photon) = 0  !M1_photons_TODO
    end if     

    {^KSP&
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    energy_good = .true.
    !! ps2 = wcons, ps1 = wold
    !! Ut_explicit has to be wold
    N = wold(ix^D,nrad^KSP_)  
    wrad(m1_energy_) = wold(ix^D,erad^KSP_)
    {^C& wrad(m1_flux^C_) = wold(ix^D,frad^KSP^C_) \}
    {^C& vel_(^C) = wold(ix^D,u0_+^C) \}

    eas_ixD(Q_ems) = wradimpl(ix^D,Q_er^KSP_)
    eas_ixD(Q_ems) = wradimpl(ix^D,Q_er^KSP_)
    eas_ixD(k_a) = wradimpl(ix^D,kappa_a^KSP_)
    eas_ixD(k_s) = wradimpl(ix^D,kappa_s^KSP_)
    ! TEST::
    !eas_ixD(k_s) = 0.0d0
    eas_ixD(Q_ems_n) = wradimpl(ix^D,Q_nr^KSP_)
    eas_ixD(k_n) = wradimpl(ix^D,kappa_nr^KSP_)

    {#IFDEF M1_TESTS
    call m1_get_implicit_collisional_sources(wrad,metricM1,vel_,eas_ixD,&
    ixI^L,ix^D,N,sources(:,^KSP),qdt,^KSP,dtfactor)
    }
    {#IFNDEF M1_TESTS
    if((wcons(ix^D,rho_) .le. small_rho) .or. (wcons(ix^D,T_eps_) .le. small_temp) .or. (wcons(ix^D,ye_) .ge. big_ye) .or. (wcons(ix^D,ye_) .le. eos_yemin)) then
      !M1_TODO: once testing is over remove setting eas to zero here
      eas_ixD(:)  = 1.0d-30 !0.0d0
      !N = m1_E_atmo
      sources(:,:) = 0.0d0
      {#IFDEF DEBUG_BACKREACT
      if(qtc .gt. t_backreact) then
      write(88,*)"DYSP ixD",ix^D
      write(88,*)"in correct kappa: at t", qtC
      write(88,*)"set eas to zero in fluid lims"
      end if 
      }
    else if((wcons(ix^D,rho_) .le. exp(logrho_min_IV)) .or. (wcons(ix^D,rho_) .ge. exp(logrho_max_IV)) .or. &
          (wcons(ix^D,T_eps_) .le. exp(logtemp_min_IV)) .or. (wcons(ix^D,T_eps_) .ge. exp(logtemp_max_IV)) .or. &
          (wcons(ix^D,ye_) .le. ye_min_IV) .or. (wcons(ix^D,ye_) .ge. ye_max_IV)) then
      !M1_TODO: once testing is over remove setting eas to zero here
      eas_ixD(:)  = 1.0d-30 !0.0d0
      !N = m1_E_atmo
      sources(:,:) = 0.0d0
      {#IFDEF DEBUG_BACKREACT
      if(qtc .gt. t_backreact) then
      write(88,*)"DYSP ixD",ix^D
      write(88,*)"in correct kappa: at t", qtC
      write(88,*)"set eas to zero in weakhub lims"
      endif
      }
    else
       call m1_get_implicit_collisional_sources(wrad,metricM1,vel_,eas_ixD,&
         ixI^L,ix^D,N,sources(:,^KSP),qdt,^KSP,dtfactor)
    end if 
    }

    if(M1_FLUID_BACKREACT) then

      !----------------- DYSP
      {#IFDEF DY_SP
      !if((wcons(ix^D,rho_) .ge. small_rho) .and. (wcons(ix^D,T_eps_) .ge. small_temp) .and. (wcons(ix^D,ye_) .le. big_ye) ) then     if((wcons(ix^D,rho_) .ge. small_rho) .and. (wcons(ix^D,T_eps_) .ge. small_temp) .and. (wcons(ix^D,ye_) .le. big_ye) .and. (wcons(ix^D,ye_) .ge. eos_yemin)) then
      !if((wcons(ix^D,rho_) .ge. small_rho) .and. (wcons(ix^D,T_eps_) .ge. small_temp) .and. (wcons(ix^D,ye_) .le. big_ye) .and. (wcons(ix^D,ye_) .ge. eos_yemin)) then 
      if((qtC .gt. t_backreact) .and. (wcons(ix^D,rho_) .gt. small_rho) .and. (wcons(ix^D,T_eps_) .gt. small_temp) .and. (wcons(ix^D,ye_) .lt. big_ye) .and. (wcons(ix^D,ye_) .gt. eos_yemin)) then
      !if((wcons(ix^D,rho_) .gt. small_rho) .and. (wcons(ix^D,T_eps_) .gt. exp(logtemp_min_IV)) .and. (wcons(ix^D,ye_) .lt. ye_max_IV) .and. (wcons(ix^D,ye_) .gt. ye_min_IV)) then
       ! NOTE: qdt*dtfactor
       ! update fluid, m1 backreaction on fluid
        !-------------------------

        if(wcons(ix^D,T_eps_) .gt. 1000.d0 ) then
          write(88,*)"------------------------"
          write(88,*)"DYSP ixD",ix^D
          write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
          write(88,*)"in correct kappa: at t", qtC
          write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
          write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
          write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
          write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
          write(88,*) "wimpl: Qems Qn",wradimpl(ix^D,Q_er^KSP_),wradimpl(ix^D,Q_nr^KSP_)
          write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
          write(88,*) "wimpl: ka ks kn",wradimpl(ix^D,kappa_a^KSP_),wradimpl(ix^D,kappa_s^KSP_),wradimpl(ix^D,kappa_nr^KSP_)
          do i=1,m1_numvars_internal+1
            write(88,*) "i sources", i, sources(i,^KSP)
          end do 
          if(wcons(ix^D,T_eps_) .gt. 1000.d0 ) then
            write(88,*)"DANGER: fluid T too large!",wcons(ix^D,T_eps_)
          end if 
          {^C& write(88,*)"vel_i",vel_(^C) \}
          write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
          write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
          {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
          write(88,*)"wprim: T",wcons(ix^D,T_eps_)
          write(88,*)"wprim: Ye",wcons(ix^D,ye_)
          write(88,*)"wprim: Rho",wcons(ix^D,rho_)
          write(88,*)"wcons: tau",wcons(ix^D,tau_),wcons(ix^D,tau_) - qdt*dtfactor*sources(1,^KSP)
          {^C& write(88,*)"wcons: sC",wcons(ix^D,s^C_) ,wcons(ix^D,s^C_) - qdt*dtfactor*sources(1+^C,^KSP) \}
          write(88,*)"wcons: dye",wcons(ix^D,tau_),wcons(ix^D,dye_) - qdt*dtfactor * signum(^KSP) * sources(1+^NC+1,^KSP)   
        end if 


        do i=1,m1_numvars_internal+1
        if(sources(i,^KSP) .ne. sources(i,^KSP)) then
        !!  write(88,*)"DYSP ixD",ix^D
        !!  write(88,*)"in correct kappa: at t", qtC
        !!  write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
        !!  write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
        !!  write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
        !!  write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
        !!  write(88,*) "wimpl: Qems Qn",wradimpl(ix^D,Q_er^KSP_),wradimpl(ix^D,Q_nr^KSP_)
        !!  write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
        !!  write(88,*) "wimpl: ka ks kn",wradimpl(ix^D,kappa_a^KSP_),wradimpl(ix^D,kappa_s^KSP_),wradimpl(ix^D,kappa_nr^KSP_)
        !!  write(88,*) "i sources", i, sources(i,^KSP)
          write(88,*) "SOURCES ARE NAN"
          write(88,*) "i sources", i, sources(i,^KSP)
          sources(i,^KSP) = 0.0d0
          wrad(m1_energy_) = wold(ix^D,erad^KSP_)
          {^C& wrad(m1_flux^C_) = wold(ix^D,frad^KSP^C_) \}
           N = wold(ix^D,nrad^KSP_)
        end if 
        end do 


        ! propose raw updates (fluid gets -sources)
        dE_source   = - qdt*dtfactor * sources(1,^KSP)
        dS_source1  = - qdt*dtfactor * sources(1+1,^KSP)
        dS_source2  = - qdt*dtfactor * sources(1+2,^KSP)
        dS_source3  = - qdt*dtfactor * sources(1+3,^KSP)   ! if 3D
        if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
          dYeC_source = - qdt*dtfactor * signum(^KSP) * sources(1+^NC+1,^KSP)
        else
          dYeC_source = 0.0d0
        endif

        !if(m1_use_muons) then
        !     if((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mubar)) then
        !      dYmuC_source = - qdt*dtfactor * signum(^KSP) * sources(1+^NC+2, ^KSP)
        !     else
        !       dYmuC_source = 0.0d0
        !     end if 
        !end if 

        damping = 1.0d0
        if(M1_Iterative_Damping) then
            if((qtC > t_backreact) .and. (qtC < t_backreact + stepsBetweenDamp*qdt)) then
              damping = 0.005d0
            else if((qtC > t_backreact+ stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 2.0 * stepsBetweenDamp*qdt)) then
              damping = 0.01d0
            else if((qtC > t_backreact+ 2.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 3.0 * stepsBetweenDamp*qdt)) then
              damping = 0.02d0
            else if((qtC > t_backreact+ 3.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 4.0 * stepsBetweenDamp*qdt)) then
              damping = 0.05d0
            else if((qtC > t_backreact+ 4.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 5.0*stepsBetweenDamp*qdt)) then
              damping = 0.1d0
            else if((qtC > t_backreact+ 5.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 6.0*stepsBetweenDamp*qdt)) then
              damping = 0.2d0
            else if((qtC > t_backreact+ 6.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 7.0*stepsBetweenDamp*qdt)) then
              damping = 0.3d0
            else if((qtC > t_backreact+ 7.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 8.0*stepsBetweenDamp*qdt)) then
              damping = 0.4d0
            else if((qtC > t_backreact+ 8.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 9.0*stepsBetweenDamp*qdt)) then
              damping = 0.5d0
            else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
              damping = 0.6d0
            else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
              damping = 0.7d0
            else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
              damping = 0.8d0
            else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
              damping = 0.9d0
            else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
              damping = 1.0d0
            end if 
        end if 

        if(m1_apply_limits) then
           ! fractional changes (use vector norm for momentum)
           tiny = 1.0d-60
           Snorm    = sqrt(wcons(ix^D,s1_)**2 + wcons(ix^D,s2_)**2 + wcons(ix^D,s3_)**2)
           dS_norm  = sqrt(dS_source1**2 + dS_source2**2 + dS_source3**2)
           
           fracE  = dabs(dE_source) / max(tiny, dabs(wcons(ix^D,tau_)))
           fracS  = dS_norm / max(tiny, Snorm)
           fracS = 0.d0
           fracS = max(fracS, dabs(dS_source1) / max(tiny, dabs(wcons(ix^D,s1_))))
           fracS = max(fracS, dabs(dS_source2) / max(tiny, dabs(wcons(ix^D,s2_))))
           fracS = max(fracS, dabs(dS_source3) / max(tiny, dabs(wcons(ix^D,s3_)))) 
   
           
           ! for lepton: combine a relative cap (if dye is a conserved density)
           ! and an absolute cap on deltaYe per step to protect EOS
           {#IFDEF TABEOS
           rel_dye_den = max(tiny, dabs(wcons(ix^D,dye_)))
           fracYe_rel  = dabs(dYeC_source) / rel_dye_den
           absYe_cap   = 2.0d-3     ! e.g. limit absolute deltaYe per substep
           fracYe_abs  = 0.0d0
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             fracYe_abs = dabs(dYeC_source) / max(tiny, absYe_cap)  ! treat like a fraction vs cap
           endif      
           !if(m1_use_muons) then
           !   ! for lepton: combine a relative cap (if dye is a conserved density)
           !   ! and an absolute cap on deltaYmu per step to protect EOS
           !   rel_dymu_den = max(tiny, dabs(wcons(ix^D,dymu_)))
           !   fracYmu_rel  = dabs(dYmuC_source) / rel_dymu_den
           !   absYmu_cap   = 2.0d-3     ! e.g. limit absolute deltaYmu per substep
           !   fracYmu_abs  = 0.0d0
           !   if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !     fracYmu_abs = dabs(dYmuC_source) / max(tiny, absYu_cap)  ! treat like a fraction vs cap
           !   endif
           !end if 
           }
           
           ! user-set caps
           capE  = couple_tau      ! e.g. 0.02
           capS  = couple_momenta    ! e.g. 0.02
           capYe = couple_ye         ! e.g. 0.02 (relative), plus abs cap above
           !capYmu = couple_ymu      ! e.g. 0.02 (relative), plus abs cap above
           
           factorA = 1.0d0
           factorA = min(factorA, capE  / max(tiny, fracE))
           factorA = min(factorA, capS  / max(tiny, fracS))
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             factorA = min(factorA, capYe / max(tiny, fracYe_rel))
             factorA = min(factorA, 1.0d0 / max(tiny, fracYe_abs))
           endif
           !if(m1_use_muons) then
           !  if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !    factorA = min(factorA, capYmu / max(tiny, fracYmu_rel))
           !    factorA = min(factorA, 1.0d0 / max(tiny, fracYmu_abs))
           !  endif
           !end if 

           factorA = factorA * damping
           ! rescale the upadted radiation by same factorA as fluid coupling
           deltaErad  = (wrad(m1_energy_) - wold(ix^D,erad^KSP_)) * factorA
           {^C& deltaFrad(^C)  = (wrad(m1_flux^C_) - wold(ix^D,frad^KSP^C_)) * factorA \}
           deltaN = (N - wold(ix^D,nrad^KSP_)) * factorA
           ! add to old wrad
           !N = wold(ix^D,nrad^KSP_) + deltaN
           !wrad(m1_energy_) = wold(ix^D,erad^KSP_) + deltaErad
           !{C wrad(m1_flux^C_)  = wold(ix^D,frad^KSP^C_)  + deltaFrad(^C)}
           !!!
           !wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * deltaErad
           !{C wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * deltaFrad(^C) }
        
           ! now apply the SAME factorA to all fluid channels
           wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * dE_source
           {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * dS_source^C \}
           {#IFDEF TABEOS
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             dye_tmp = wcons(ix^D,dye_) + factorA * dYeC_source
             ye_tmp = dye_tmp / wcons(ix^D,d_)
             Delta_Ye = dabs( ye_tmp - wcons(ix^D,dye_)/ wcons(ix^D,d_))
             error_Ye = Delta_Ye / rel_error_Ye
             if(ye_tmp < eos_yemin .or. ye_tmp > big_ye) then
                write(88,*) "m1-backreaction: ye is out of eos-bounds" 
                write(88,*) "ye yemin yemax",ye_tmp,eos_yemin,big_ye
             else if(error_Ye > 10.0d0) then
                write(88,*)"-----------------------------"
                write(88,*) "m1-backreaction: the relative ye-error is > 10"               
                write(88,*)"DYSP ixD",ix^D
                write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
                write(88,*)"in correct kappa: at t", qtC
                write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
                write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
                write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
                write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
                write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
                write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
                write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
                write(88,*) "dYe-source", dYeC_source
                {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
             else 
                wcons(ix^D,dye_) = wcons(ix^D,dye_) + factorA * dYeC_source
             end if 
             ! TEST !!!
             !wcons(ix^D,dye_) = wcons(ix^D,dye_) - factorA * dYeC_source
           endif
           !if(m1_use_muons) then
           !  if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !    wcons(ix^D,dymu_) = wcons(ix^D,dymu_) + factorA * dYmuC_source
           !  endif
           !end if 
           }
        else 
           factorA = 1.0d0
           !factorA = 0.95d0
           !damping = 1.0d0
           factorA = factorA * damping
           ! now apply the SAME factorA to all fluid channels
           tau_tmp =  wcons(ix^D,tau_) + factorA * dE_source
           if(tau_tmp<0.0d0) then
              energy_good = .false.
              write(88,*) "m1-backreaction: energy tau is < 0"
           else
              wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * dE_source
           endif 
           ! TEST !!!
           !wcons(ix^D,tau_) = wcons(ix^D,tau_) - factorA * dE_source
           {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * dS_source^C \}
           ! TEST: !!
           !{C wcons(ix^D,s^C_) = wcons(ix^D,s^C_) - factorA * dS_source^C}
           {#IFDEF TABEOS
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             dye_tmp = wcons(ix^D,dye_) + factorA * dYeC_source
             ye_tmp = dye_tmp / wcons(ix^D,d_)
             Delta_Ye = dabs( ye_tmp - wcons(ix^D,dye_)/ wcons(ix^D,d_))
             error_Ye = Delta_Ye / rel_error_Ye
             if(ye_tmp < eos_yemin .or. ye_tmp > big_ye) then
                write(88,*) "m1-backreaction: ye is out of eos-bounds" 
                write(88,*) "ye yemin yemax",ye_tmp,eos_yemin,big_ye
             else if(error_Ye > 10.0d0) then
                write(88,*)"-----------------------------"
                write(88,*) "m1-backreaction: the relative ye-error is > 10"               
                write(88,*)"DYSP ixD",ix^D
                write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
                write(88,*)"in correct kappa: at t", qtC
                write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
                write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
                write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
                write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
                write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
                write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
                write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
                write(88,*) "dYe-source", dYeC_source
                {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
             else 
                wcons(ix^D,dye_) = wcons(ix^D,dye_) + factorA * dYeC_source
             end if 
             ! TEST !!!
             !wcons(ix^D,dye_) = wcons(ix^D,dye_) - factorA * dYeC_source
           endif
           !if(m1_use_muons) then
           !  if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !    wcons(ix^D,dymu_) = wcons(ix^D,dymu_) + factorA * dYmuC_source
           !  endif
           !end if 
           }
        end if 
        ! now apply the SAME factorA to all fluid channels
         !! !! TEST
         !! !factorA = 1.0d0
         !! !!
         !! wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * dE_source
         !! {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * dS_source^C \}
         !! {#IFDEF TABEOS
         !! if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
         !!   wcons(ix^D,dye_) = wcons(ix^D,dye_) + factorA * dYeC_source
         !! endif
         !! !if(m1_use_muons) then
         !! !  if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
         !! !    wcons(ix^D,dymu_) = wcons(ix^D,dymu_) + factorA * dYmuC_source
         !! !  endif
         !! !end if 
         !! }
      else
          !!write(88,*)"could not add sources to ixd",ix^D
          !!write(88,*)"rho  small_rho", wcons(ix^D,rho_),small_rho
          !!write(88,*)"T    small_temp", wcons(ix^D,T_eps_),small_temp
          !!write(88,*)"ye     yemax", wcons(ix^D,ye_),big_ye
      end if 
      } ! end if DYSP

      !*****************************************************
      !----------------- nonDYSP
      {#IFNDEF DY_SP
      if((qtC .gt. t_backreact) .and. (wcons(ix^D,rho_) .ge. small_rho) .and. (wcons(ix^D,ye_) .le. big_ye) .and. (wcons(ix^D,ye_) .ge. eos_yemin)) then
        ! NOTE: qdt*dtfactor
        ! update fluid, m1 backreaction on fluid
         !-------------------------
         do i=1,m1_numvars_internal+1
         if(sources(i,^KSP) .ne. sources(i,^KSP)) then
           write(88,*)"nonDYSP ixD",ix^D
           write(88,*) "Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
           write(88,*) "ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
           write(88,*)"i sources", i, sources(i,^KSP)
         end if 
         end do 
         !---------------------         
        ! propose raw updates (fluid gets -sources)
        dE_source   = - qdt*dtfactor * sources(1,^KSP)
        dS_source1  = - qdt*dtfactor * sources(1+1,^KSP)
        dS_source2  = - qdt*dtfactor * sources(1+2,^KSP)
        dS_source3  = - qdt*dtfactor * sources(1+3,^KSP)   ! if 3D
        if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
          dYeC_source = - qdt*dtfactor * signum(^KSP) * sources(1+^NC+1,^KSP)
        else
          dYeC_source = 0.0d0
        endif
        !if(m1_use_muons) then
        !     if((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mubar)) then
        !      dYmuC_source = - qdt*dtfactor * signum(^KSP) * sources(1+^NC+2, ^KSP)
        !     else
        !       dYmuC_source = 0.0d0
        !     end if 
        !end if 

        damping = 1.0d0
        if(M1_Iterative_Damping) then
          if((qtC > t_backreact) .and. (qtC < t_backreact + stepsBetweenDamp*qdt)) then
            damping = 0.005d0
          else if((qtC > t_backreact+ stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 2.0 * stepsBetweenDamp*qdt)) then
            damping = 0.01d0
          else if((qtC > t_backreact+ 2.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 3.0 * stepsBetweenDamp*qdt)) then
            damping = 0.02d0
          else if((qtC > t_backreact+ 3.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 4.0 * stepsBetweenDamp*qdt)) then
            damping = 0.05d0
          else if((qtC > t_backreact+ 4.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 5.0*stepsBetweenDamp*qdt)) then
            damping = 0.1d0
          else if((qtC > t_backreact+ 5.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 6.0*stepsBetweenDamp*qdt)) then
            damping = 0.2d0
          else if((qtC > t_backreact+ 6.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 7.0*stepsBetweenDamp*qdt)) then
            damping = 0.3d0
          else if((qtC > t_backreact+ 7.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 8.0*stepsBetweenDamp*qdt)) then
            damping = 0.4d0
          else if((qtC > t_backreact+ 8.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 9.0*stepsBetweenDamp*qdt)) then
            damping = 0.5d0
          else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
            damping = 0.6d0
          else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
            damping = 0.7d0
          else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
            damping = 0.8d0
          else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
            damping = 0.9d0
          else if((qtC > t_backreact+ 9.0*stepsBetweenDamp*qdt) .and. (qtC < t_backreact + 10.0*stepsBetweenDamp*qdt)) then
            damping = 1.0d0
          end if 
      end if 
        
        if(m1_apply_limits) then
           ! fractional changes (use vector norm for momentum)
           tiny = 1.0d-60
           Snorm    = sqrt(wcons(ix^D,s1_)**2 + wcons(ix^D,s2_)**2 + wcons(ix^D,s3_)**2)
           dS_norm  = sqrt(dS_source1**2 + dS_source2**2 + dS_source3**2)
           
           fracE  = dabs(dE_source) / max(tiny, dabs(wcons(ix^D,tau_)))
           !fracS  = dS_norm / max(tiny, Snorm)
           fracS = 0.d0
           fracS = max(fracS, dabs(dS_source1) / max(tiny, dabs(wcons(ix^D,s1_))))
           fracS = max(fracS, dabs(dS_source2) / max(tiny, dabs(wcons(ix^D,s2_))))
           fracS = max(fracS, dabs(dS_source3) / max(tiny, dabs(wcons(ix^D,s3_))))       
           ! for lepton: combine a relative cap (if dye is a conserved density)
           ! and an absolute cap on deltaYe per step to protect EOS
           {#IFDEF TABEOS
           rel_dye_den = max(tiny, dabs(wcons(ix^D,dye_)))
           fracYe_rel  = dabs(dYeC_source) / rel_dye_den
           absYe_cap   = 2.0d-3     ! e.g. limit absolute deltaYe per substep
           fracYe_abs  = 0.0d0
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             fracYe_abs = dabs(dYeC_source) / max(tiny, absYe_cap)  ! treat like a fraction vs cap
           endif      
           !if(m1_use_muons) then
           !   ! for lepton: combine a relative cap (if dye is a conserved density)
           !   ! and an absolute cap on deltaYmu per step to protect EOS
           !   rel_dymu_den = max(tiny, dabs(wcons(ix^D,dymu_)))
           !   fracYmu_rel  = dabs(dYmuC_source) / rel_dymu_den
           !   absYmu_cap   = 2.0d-3     ! e.g. limit absolute deltaYmu per substep
           !   fracYmu_abs  = 0.0d0
           !   if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !     fracYmu_abs = dabs(dYmuC_source) / max(tiny, absYu_cap)  ! treat like a fraction vs cap
           !   endif
           !end if 
           }
           
           ! user-set caps
           capE  = couple_tau      ! e.g. 0.02
           capS  = couple_momenta    ! e.g. 0.02
           capYe = couple_ye         ! e.g. 0.02 (relative), plus abs cap above
           !capYmu = couple_ymu      ! e.g. 0.02 (relative), plus abs cap above
           
           factorA = 1.0d0
           factorA = min(factorA, capE  / max(tiny, fracE))
           factorA = min(factorA, capS  / max(tiny, fracS))
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             factorA = min(factorA, capYe / max(tiny, fracYe_rel))
             factorA = min(factorA, 1.0d0 / max(tiny, fracYe_abs))
           endif
           !if(m1_use_muons) then
           !  if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !    factorA = min(factorA, capYmu / max(tiny, fracYmu_rel))
           !    factorA = min(factorA, 1.0d0 / max(tiny, fracYmu_abs))
           !  endif
           !end if 
           
           factorA = factorA * damping
           ! rescale the upadted radiation by same factorA as fluid coupling
           deltaErad  = (wrad(m1_energy_) - wold(ix^D,erad^KSP_)) * factorA
           {^C& deltaFrad(^C)  = (wrad(m1_flux^C_) - wold(ix^D,frad^KSP^C_)) * factorA \}
           deltaN = (N - wold(ix^D,nrad^KSP_)) * factorA
           ! add to old wrad
           !N = wold(ix^D,nrad^KSP_) + deltaN
           !wrad(m1_energy_) = wold(ix^D,erad^KSP_) + deltaErad
           !{C wrad(m1_flux^C_)  = wold(ix^D,frad^KSP^C_)  + deltaFrad(^C)}
           !!!
           !wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * deltaErad
           !{C wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * deltaFrad(^C) }
        
           ! now apply the SAME factorA to all fluid channels
           wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * dE_source
           {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * dS_source^C \}
           {#IFDEF TABEOS
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             dye_tmp = wcons(ix^D,dye_) + factorA * dYeC_source
             ye_tmp = dye_tmp / wcons(ix^D,d_)
             Delta_Ye = dabs( ye_tmp - wcons(ix^D,dye_)/ wcons(ix^D,d_))
             error_Ye = Delta_Ye / rel_error_Ye
             if(ye_tmp < eos_yemin .or. ye_tmp > big_ye) then
                write(88,*) "m1-backreaction: ye is out of eos-bounds" 
                write(88,*) "ye yemin yemax",ye_tmp,eos_yemin,big_ye
             else if(error_Ye > 10.0d0) then
                write(88,*)"-----------------------------"
                write(88,*) "m1-backreaction: the relative ye-error is > 10"               
                write(88,*)"DYSP ixD",ix^D
                write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
                write(88,*)"in correct kappa: at t", qtC
                write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
                write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
                write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
                write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
                write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
                write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
                write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
                write(88,*) "dYe-source", dYeC_source
                {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
             else 
                wcons(ix^D,dye_) = wcons(ix^D,dye_) + factorA * dYeC_source
             end if 
             ! TEST !!!
             !wcons(ix^D,dye_) = wcons(ix^D,dye_) - factorA * dYeC_source
           endif
           !if(m1_use_muons) then
           !  if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !    wcons(ix^D,dymu_) = wcons(ix^D,dymu_) + factorA * dYmuC_source
           !  endif
           !end if 
           }
        else 
           factorA = 1.0d0
           !factorA = 0.95d0
           !damping = 1.0d0
           factorA = factorA * damping
           ! now apply the SAME factorA to all fluid channels
           tau_tmp =  wcons(ix^D,tau_) + factorA * dE_source
           if(tau_tmp<0.0d0) then
              energy_good = .false.
              write(88,*) "m1-backreaction: energy tau is < 0"
           else
              wcons(ix^D,tau_) = wcons(ix^D,tau_) + factorA * dE_source
           endif 
           ! TEST !!!
           !wcons(ix^D,tau_) = wcons(ix^D,tau_) - factorA * dE_source
           {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) + factorA * dS_source^C \}
           ! TEST: !!
           !{C wcons(ix^D,s^C_) = wcons(ix^D,s^C_) - factorA * dS_source^C}
           {#IFDEF TABEOS
           if ((^KSP .eq. m1_i_nue) .or. (^KSP .eq. m1_i_nuebar)) then
             dye_tmp = wcons(ix^D,dye_) + factorA * dYeC_source
             ye_tmp = dye_tmp / wcons(ix^D,d_)
             Delta_Ye = dabs( ye_tmp - wcons(ix^D,dye_)/ wcons(ix^D,d_))
             error_Ye = Delta_Ye / rel_error_Ye
             if(ye_tmp < eos_yemin .or. ye_tmp > big_ye) then
                write(88,*) "m1-backreaction: ye is out of eos-bounds" 
                write(88,*) "ye yemin yemax",ye_tmp,eos_yemin,big_ye
             else if(error_Ye > 10.0d0) then
                write(88,*)"-----------------------------"
                write(88,*) "m1-backreaction: the relative ye-error is > 10"               
                write(88,*)"DYSP ixD",ix^D
                write(88,*)"x",x(ix^D,1),x(ix^D,2),x(ix^D,3)
                write(88,*)"in correct kappa: at t", qtC
                write(88,*)"fluid: rho, smallRho",wcons(ix^D,rho_),small_rho
                write(88,*)"fluid: T, smallT",wcons(ix^D,T_eps_),small_temp
                write(88,*)"fluid: ye, yemin,bigye",wcons(ix^D,ye_),eos_yemin,big_ye
                write(88,*) "eas:   Qems,Qn",eas_ixD(Q_ems),eas_ixD(Q_ems_n)
                write(88,*) "eas:   ka ks kn",eas_ixD(k_a),eas_ixD(k_s),eas_ixD(k_n)
                write(88,*)"wrad: N",N/metricM1%sqrtg(ix^D)
                write(88,*)"wrad: E",wrad(m1_energy_)/metricM1%sqrtg(ix^D)
                write(88,*) "dYe-source", dYeC_source
                {^C& write(88,*)"wrad: Fi",wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) \}
             else 
                wcons(ix^D,dye_) = wcons(ix^D,dye_) + factorA * dYeC_source
             end if 
             ! TEST !!!
             !wcons(ix^D,dye_) = wcons(ix^D,dye_) - factorA * dYeC_source
           endif
           !if(m1_use_muons) then
           !  if ((^KSP .eq. m1_i_mu) .or. (^KSP .eq. m1_i_mu)) then
           !    wcons(ix^D,dymu_) = wcons(ix^D,dymu_) + factorA * dYmuC_source
           !  endif
           !end if 
           }
        end if ! endif apply_limits
      else
          !! write(88,*)"could not add sources to ixd",ix^D
          !! write(88,*)"rho  small_rho", wcons(ix^D,rho_),small_rho
          !! write(88,*)"T    small_temp", wcons(ix^D,T_eps_),small_temp
          !! write(88,*)"ye     yemax", wcons(ix^D,ye_),big_ye
      end if 
      }  !-------- endif nonDYSP
      !

    end if  ! end M1_BAKCREACT

    !--------------------------
    !!if(m1_apply_limits) then
    !!   {^C& wrad(m1_flux^C_) = dsign(max(m1_E_atmo, 1d-10*wrad(m1_energy_)), wrad(m1_flux^C_)) \}
    !!   
    !!   ! Ensure positivity and Fmag <=E
    !!   Er = max(wrad(m1_energy_)/metricM1%sqrtg(ix^D), m1_E_atmo)
    !!   Fmag = {^C&wrad(m1_flux^C_)/metricM1%sqrtg(ix^D)*wrad(m1_flux^C_)/metricM1%sqrtg(ix^D) +}
    !!   Fmag = dsqrt( Fmag)
    !!   if (Fmag > fmaxfact * wrad(m1_energy_)/metricM1%sqrtg(ix^D)) then
    !!     facFrad = (fmaxfact * wrad(m1_energy_)/metricM1%sqrtg(ix^D)) / max(Fmag,1.d-300)
    !!     {^C& wrad(m1_flux^C_) = facFrad * wrad(m1_flux^C_) \}
    !!   endif
    !!   !____
    !!   !Er = max(wrad(m1_energy_), m1_E_atmo)
    !!   !Fmag = {^C&wrad(m1_flux^C_)*wrad(m1_flux^C_) +}
    !!   !Fmag = dsqrt( Fmag)
    !!   !if (Fmag > fmaxfact * wrad(m1_energy_)) then
    !!   !  facFrad = (fmaxfact * wrad(m1_energy_)) / max(Fmag,1.d-300)
    !!  ! {C wrad(m1_flux^C_) = facFrad * wrad(m1_flux^C_) }
    !!  ! endif
    !!end if 
    !!!------------------------

    ! update rad
    if(energy_good) then
    wcons(ix^D, nrad^KSP_) = max(N,m1_E_atmo) !0.0d0) !N
    wcons(ix^D, erad^KSP_) = max(wrad(m1_energy_),m1_E_atmo) !0.0d0)  !wrad(m1_energy_)
    {^C& wcons(ix^D, frad^KSP^C_) = wrad(m1_flux^C_) \}
    end if

    {^D& end do \}

    \} ! end KSP

    call metricM1%destroy()

    
  end subroutine m1_add_collisional_sources
   } ! end IFNDEF M1_EXPLICIT
   } ! end IFNDEF UNIT_TESTS 
  
end module mod_m1
