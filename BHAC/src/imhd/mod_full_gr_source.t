module mod_full_gr_source
  implicit none
  public


  contains
  ! it is only for spatial diagonal metric
  subroutine get_g_up_munu(g, alp, beta, gamma, ixI^L, ixO^L)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: alp(ixI^S), beta(ixI^S, 1:3)
    double precision, intent(in)    :: gamma(ixI^S, 1:3, 1:3)
    double precision, intent(out)   :: g(ixI^S, 0:3, 0:3)

    integer                         :: mu, nu

    g(ixO^S, 0:3,0:3) = 0.0d0

    g(ixO^S, 0,0) = - 1.0d0 / alp(ixO^S)**2
    do mu = 1, 3
       g(ixO^S, 0, mu) = beta(ixO^S, mu) / alp(ixO^S)**2
       g(ixO^S, mu, 0) = g(ixO^S, 0, mu)
    end do

    ! this is actually gamma^{ij}
    do mu = 1, 3
       g(ixO^S, mu, mu) = 1.0d0 / gamma(ixO^S, mu, mu)
    end do

    ! g^{ij} = gamma^{ij} - beta^i beta^j / alpha^2
    do mu = 1, 3
       do nu = 1, 3
          g(ixO^S, mu, nu) = g(ixO^S, mu, nu) - &
                             beta(ixO^S, mu) * beta(ixO^S, nu) / alp(ixO^S)**2
       end do
    end do

  end subroutine get_g_up_munu

  ! it is only for spatial diagonal metric
  subroutine get_g_down_munu(g, alp, beta, gamma, ixI^L, ixO^L)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: alp(ixI^S), beta(ixI^S, 1:3)
    double precision, intent(in)    :: gamma(ixI^S, 1:3, 1:3)
    double precision, intent(out)   :: g(ixI^S, 0:3, 0:3)

    integer                         :: mu, nu

    g(ixO^S, 0:3,0:3) = 0.0d0

    g(ixO^S, 0,0) = - alp(ixO^S)**2 + {^C& beta(ixO^S, ^C)**2 * gamma(ixO^S, ^C, ^C) +}
    do mu = 1, 3
       g(ixO^S, 0, mu) = beta(ixO^S, mu) * gamma(ixO^S, mu, mu)
       g(ixO^S, mu, 0) = g(ixO^S, 0, mu)
    end do

    ! this is actually gamma_{ij}
    do mu = 1, 3
       g(ixO^S, mu, mu) = gamma(ixO^S, mu, mu)
    end do

  end subroutine get_g_down_munu
  
  subroutine addgeometry_source_full_gr(qdt,ixI^L,ixO^L,wCT,w,wold,x,qtC, metric_M1)
    use mod_metric
    use mod_eos
    use mod_imhd_intermediate	
    use mod_cfc_parameters
    use mod_m1_metric_interface
	 {#IFDEF M1
    use mod_m1 
    }
    include 'amrvacdef.f'

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim), qtC
    double precision, intent(inout) :: wCT(ixI^S, 1:nw)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: wold(ixI^S, 1:nw)
	 type(m1_metric_helper), intent(in) :: metric_M1 

    {^NOONED
    double precision                :: cot_theta(ixI^S)
    }
    double precision                :: add_source_tau(ixI^S)
    double precision                :: add_source_mom(ixI^S,1:ndir)
    double precision,dimension(ixI^S,1:ndir)  :: vU, bsU, bsD, sU, wikdjgammaik
    double precision,dimension(ixI^S)         :: VdotB, sqrB, tmp, tmp2, wikbetajdjgammaik, sijkij
    double precision                :: lrho(ixI^S), lye(ixI^S), ptot_mom(ixI^S)

    double precision                :: htot(ixI^S) ! enthalpy
    double precision                :: Ptot(ixI^S) 
    double precision                :: prs_tmp(ixI^S), eps_tmp(ixI^S), temp_local, rho_local
    double precision                :: lfac(ixI^S) ! Lorentz factor
    ! metric variables
    double precision                :: alp(ixI^S), alp_prime(ixI^S)
    double precision                :: psi(ixI^S)
    double precision                :: psi4(ixI^S)
    double precision                :: betai(ixI^S,1:3)
    ! dervitaves of the metric variables
    double precision                :: dalp(ixI^S,1:3), dpsi(ixI^S,1:3), dalp_prime(ixI^S,1:3)
    double precision                :: dbeta(ixI^S,1:3,1:3)

    ! 3-metric gamma_ij
    double precision                :: gamma(ixI^S,1:3,1:3)
    double precision                :: d_k_gamma_ij(ixI^S,1:3,1:3,1:3)
    double precision                :: d_gamma_ij_hat(ixI^S,1:3,1:3,1:3)
    ! 3-metric gamma^ij
    double precision                :: gammainv(ixI^S,1:3,1:3)
	! sqrt gamma
	double precision                :: gammasqrt(ixI^S) 

    ! extrinsic curvature K_ij
    double precision                :: K_ij(ixI^S,1:3,1:3)
    ! 4-metric g^{munu}
    double precision                :: g(ixI^S,0:3,0:3)
    ! dervitaves of the 3-metric gamma_{jk},i
    double precision                :: Dgamma(ixI^S,1:3,1:3,1:3)
    ! Christoffel symbols of the reference 3-metric gamma_hat
    double precision                :: christoffel(ixI^S,1:3,1:3,1:3)
    ! covariant dervitaves of beta D_i beta^k
    double precision                :: D_beta(ixI^S,1:3,1:3)
    ! energy-momentum tensor T^{munu}
    double precision                :: Tmunu(ixI^S,0:3,0:3)
    ! 4-velocity u^{mu}
    double precision                :: u(ixI^S,0:3), h_th(ixI^S), xi(ixI^S)
    double precision                :: bmu(ixI^S,0:3)
    integer                         :: iw,idir,jdir,kdir,ix^D,i,j,k,idims,inonzero,jnonzero,hxO^L

    call imhd_get_intermediate_variables(ixI^L, ixO^L, wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                gammainv=gammainv(ixI^S,1:3,1:3),&
                bmu=bmu(ixI^S,0:ndir), &
                lfac=lfac(ixI^S), &
                h_th=h_th(ixI^S), &
                htot=htot(ixI^S), &
                Ptot=Ptot(ixI^S), &
                P_th=prs_tmp(ixI^S), &
                eps =eps_tmp(ixI^S) &
                 )

    xi(ixO^S) = lfac(ixO^S)**2 * wCT(ixO^S, rho_) * h_th(ixO^S)

    ! source term for mom in bhac version
    VdotB(ixO^S)= ({^C& wCT(ixO^S,s^C_)*wCT(ixO^S,b^C_)+})/xi(ixO^S)

    sqrB(ixO^S) = {^C& gamma(ixO^S,^C,^C) * wCT(ixO^S,b^C_)**2 +}
    {^C&  sU(ixO^S,^C) = wCT(ixO^S, s^C_) * gammainv(ixO^S,^C,^C) \}

    {^C&
         vU(ixO^S,^C)=(sU(ixO^S,^C)+VdotB(ixO^S)*wCT(ixO^S,b0_+^C))/ &
         (xi(ixO^S)+sqrB(ixO^S))
    bsU(ixO^S,^C) = wCT(ixO^S,b0_+^C)/lfac(ixO^S) + lfac(ixO^S)*VdotB(ixO^S)*vU(ixO^S,^C)
    \}

    {^C&  bsD(ixO^S,^C) = gamma(ixO^S,^C,^C) * bsU(ixO^S,^C) \}


    ptot_mom(ixO^S) = prs_tmp(ixO^S)
    ptot_mom(ixO^S) = half*(VdotB(ixO^S)**2  &
         + sqrB(ixO^S)/lfac(ixO^S)**2)  &
         + ptot_mom(ixO^S)
    ! end of source term for mom in bhac version

    add_source_tau(ixO^S)      = 0.0d0
    {^C& add_source_mom(ixO^S, ^C) = 0.0d0  \}

    if (use_gw_br .and. use_h00_coupling) then
       alp_prime(ixI^S) = dsqrt( wCT(ixI^S,alp_metric_)**2 - wCT(ixI^S,h_00_) )
       !alp_prime(ixI^S) = wCT(ixI^S,alp_metric_) - wCT(ixI^S,h_00_)/2.0d0/wCT(ixI^S,alp_metric_)
    else
       alp_prime(ixI^S) = wCT(ixI^S,alp_metric_)
    endif

    ! initialize the metric
    psi(ixI^S) = wCT(ixI^S,psi_metric_)
    if (gw_br_couple_weakly) then
      alp(ixI^S) = wCT(ixI^S,alp_metric_)
    else
      alp(ixI^S) = alp_prime(ixI^S)
    endif
    {^C&  betai(ixI^S,^C) = wCT(ixI^S, beta_metric^C_)  \}
    psi4(ixI^S)= psi(ixI^S)**4

    ! volume averaged Christoffel symbols
    !call get_christoffel(ixI^L,ixI^L,x(ixI^S,1:ndim),christoffel(ixI^S,1:3,1:3,1:3))
    !christoffel = 0.0d0

    !-----------------------------------------------------------------------
    ! Part 2: gravitational source terms, only needed when use_GR = .true.
    !-----------------------------------------------------------------------
!    if ( .not. use_GR ) return

    ! calculate derivitives of the metric variables
    dalp(ixI^S,1:3) = 0.0d0
    dbeta(ixI^S,1:3,1:3) = 0.0d0
    dpsi(ixI^S,1:3) = 0.0d0

    do idir = 1, ndim
       call partial_d( alp(ixI^S) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dalp(ixI^S,idir) )
       !call partial_d( alp_prime(ixI^S) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dalp_prime(ixI^S,idir) )
       call partial_d( psi(ixI^S) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dpsi(ixI^S,idir) )
       dpsi(ixO^S,idir) = dpsi(ixO^S,idir) / psi(ixO^S)
       do jdir = 1, ndir
          call partial_d( betai(ixI^S,jdir) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dbeta(ixI^S,jdir,idir) )
       end do
    end do

    ! covariant derivative of beta: partial_i beta^k + Gamma^k_{ij} beta^j
    !D_beta(ixO^S, 1:3, 1:3) = dbeta(ixO^S, 1:3, 1:3)
    !if ( coordinate /= cartesian ) then
    !   do idir = 1, 3
    !      do kdir = 1, 3
    !         do jdir = 1, 3
    !            D_beta(ixO^S, kdir, idir) = D_beta(ixO^S, kdir, idir) &
    !                                      + christoffel(ixO^S, kdir,idir,jdir) * betai(ixO^S,jdir)
    !         end do
    !      end do
    !   end do
    !end if

    ! partial dervitaves of metric D_i gamma_jk
    !Dgamma(ixI^S,1:3,1:3,1:3) = 0.0d0
    !do kdir = 1, ndim
    !   do idir = 1, 3
    !      Dgamma(ixO^S,idir,idir,kdir) = 4.0d0 * gamma(ixO^S,idir,idir) * dpsi(ixO^S,kdir)
    !   end do
    !end do

    ! get the d_gamma_ij_hat in particular coordinate
       call partial_d_gamma_ij_hat( x(ixI^S, 1:ndim), ixI^L, ixO^L, d_gamma_ij_hat)


    ! partial dervitaves of metric d_i gamma_jk
    d_k_gamma_ij(ixI^S,1:3,1:3,1:3) = 0.0d0
!    do kdir = 1, ndim
!     do jdir = 1, ndim
!      do idir = 1, ndim
!          d_k_gamma_ij(ixO^S,idir,jdir,kdir) = 4.0d0 * gamma(ixO^S,idir,jdir) * dpsi(ixO^S,kdir) +&
!                                          psi4(ixO^S) * d_gamma_ij_hat(ixO^S,idir,jdir,kdir)
!      end do
!     end do
!    end do

    do kdir = 1, ndim
     do jdir = 1, 3
      do idir = 1, 3
          d_k_gamma_ij(ixO^S,idir,jdir,kdir) = 4.0d0 * gamma(ixO^S,idir,jdir) * dpsi(ixO^S,kdir) +&
                                          psi4(ixO^S) * d_gamma_ij_hat(ixO^S,idir,jdir,kdir)
      end do
     end do
    end do

    ! Note that g(mu,nu) here is g^{munu}
    !call get_g_up_munu(g(ixI^S,0:3,0:3), alp(ixI^S), betai(ixI^S,1:3), gamma(ixI^S,1:3,1:3), ixI^L, ixO^L)

    ! Calculate the 4-velocity
    !u(ixO^S,0:3) = 0.0d0
    !u(ixO^S,0) = lfac(ixO^S) / alp(ixO^S)
    ! {^C&  u(ixO^S, ^C) = wCT(ixO^S, u^C_) - lfac(ixO^S) * betai(ixO^S,^C) / alp(ixO^S) \}

    !! energy-momentum tensor T^{munu}
    !do idir = 0,3
    !   do jdir = 0,3
    !      Tmunu(ixO^S,idir,jdir) = wCT(ixO^S, rho_) * htot(ixO^S) * u(ixO^S, idir) * u(ixO^S, jdir) &
    !                           + Ptot(ixO^S) * g(ixO^S,idir,jdir) &
    !                           - bmu(ixO^S, idir) * bmu(ixO^S, jdir)
    !   enddo
    !enddo

    !! Now calculate the source terms for HD (mom, tau) in the compact form
    !{^C&
    !   do jdir = 1, 3
    !      do kdir = 1, 3
    !         add_source_mom(ixO^S, ^C) = add_source_mom(ixO^S, ^C) &
    !                                   + 0.5d0 * Dgamma(ixO^S, jdir,kdir,^C) * ( &
    !                                    Tmunu(ixO^S,0,0)*betai(ixO^S,kdir)*betai(ixO^S,jdir) &
    !                                     + 2.0d0*Tmunu(ixO^S,0,jdir)*betai(ixO^S,kdir) + Tmunu(ixO^S, jdir,kdir) )
    !      end do
    !   end do
    !\}

    !{^C&
    !   add_source_mom(ixO^S, ^C) = add_source_mom(ixO^S, ^C) &
    !                             - Tmunu(ixO^S,0,0)*alp(ixO^S)*dalp(ixO^S, ^C)
    !   ! T^{0j} * gamma_jk * D_i beta^k
    !   do kdir = 1, 3
    !      do jdir = 1, 3
    !         add_source_mom(ixO^S, ^C) = add_source_mom(ixO^S, ^C) &
    !                                + Tmunu(ixO^S,0,jdir) * gamma(ixO^S,jdir,kdir) * D_beta(ixO^S,kdir,^C)
    !      end do
    !   end do
    !\}

    ! To compute the tau source term ( add_source(:,:,:,tau_) ), we need to work out the extrinsic curvature K_{ij}
    ! fixme: make get_K_ij as a stand alone subroutine
    K_ij(ixO^S,1:3,1:3) = 0.0d0

    select case (coordinate)
    case (cartesian)
       K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*alp_prime(ixO^S)) &
       !K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*alp(ixO^S)) &
                    * ( 2.0d0*dbeta(ixO^S,1,1) &
                      {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

       K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*alp_prime(ixO^S)) &
       !K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) &
                      {^NOONED + 2.0d0*dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

       K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*alp_prime(ixO^S)) &
       !K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) &
                      {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
    case (cylindrical)
       K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*alp(ixO^S)) &
                    * ( 2.0d0*dbeta(ixO^S,1,1) - betai(ixO^S,1)/x(ixO^S,r_) &
                      {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

       K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) - betai(ixO^S,1)/x(ixO^S,r_) &
                      {^NOONED + 2.0d0 * dbeta(ixO^S,2,2) } {^IFTHREED - dbeta(ixO^S,3,3) } )

       K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*alp(ixO^S)) &
                    * ( - dbeta(ixO^S,1,1) + 2.0d0 * betai(ixO^S,1)/x(ixO^S,r_) &
                      {^NOONED - dbeta(ixO^S,2,2) } {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
    case (spherical)
       {^NOONED
       cot_theta(ixO^S) = dcos(x(ixO^S,^Z))/dsin(x(ixO^S,^Z))
       !cot_theta(ixO^S) = dcos(x(ixO^S,theta_))/dsin(x(ixO^S,theta_))
       }

       K_ij(ixO^S,1,1) = gamma(ixO^S,1,1)/(3.0d0*alp(ixO^S)) &
                    * ( 2.0d0*dbeta(ixO^S,1,1) - 2.0d0*betai(ixO^S,1)/x(ixO^S,1) &
                      {^NOONED - dbeta(ixO^S,2,2) - betai(ixO^S,2)*cot_theta(ixO^S) } {^IFTHREED - dbeta(ixO^S,3,3) } )

       K_ij(ixO^S,2,2) = gamma(ixO^S,2,2)/(3.0d0*alp(ixO^S)) &
                    * ( -dbeta(ixO^S,1,1) + betai(ixO^S,1)/x(ixO^S,1) &
                      {^NOONED + 2.0d0*dbeta(ixO^S,2,2) - betai(ixO^S,2)*cot_theta(ixO^S) } {^IFTHREED - dbeta(ixO^S,3,3) } )

       K_ij(ixO^S,3,3) = gamma(ixO^S,3,3)/(3.0d0*alp(ixO^S)) &
                    * ( -dbeta(ixO^S,1,1) + betai(ixO^S,1)/x(ixO^S,1) &
                      {^NOONED - dbeta(ixO^S,2,2) + 2.0d0*betai(ixO^S,2)*cot_theta(ixO^S) } {^IFTHREED + 2.0d0 * dbeta(ixO^S,3,3) } )
    end select
    {^NOONED
    K_ij(ixO^S,1,2) = ( dbeta(ixO^S,1,2)*gamma(ixO^S,1,1) + dbeta(ixO^S,2,1)*gamma(ixO^S,2,2) ) / (2.0d0*alp_prime(ixO^S))
    K_ij(ixO^S,1,3) = ( dbeta(ixO^S,3,1)*gamma(ixO^S,3,3) + dbeta(ixO^S,1,3)*gamma(ixO^S,1,1) ) / (2.0d0*alp_prime(ixO^S))
    K_ij(ixO^S,2,3) = ( dbeta(ixO^S,3,2)*gamma(ixO^S,3,3) + dbeta(ixO^S,2,3)*gamma(ixO^S,2,2) ) / (2.0d0*alp_prime(ixO^S))
    !K_ij(ixO^S,1,2) = ( dbeta(ixO^S,1,2)*gamma(ixO^S,1,1) + dbeta(ixO^S,2,1)*gamma(ixO^S,2,2) ) / (2.0d0*alp(ixO^S))
    !K_ij(ixO^S,1,3) = ( dbeta(ixO^S,3,1)*gamma(ixO^S,3,3) + dbeta(ixO^S,1,3)*gamma(ixO^S,1,1) ) / (2.0d0*alp(ixO^S))
    !K_ij(ixO^S,2,3) = ( dbeta(ixO^S,3,2)*gamma(ixO^S,3,3) + dbeta(ixO^S,2,3)*gamma(ixO^S,2,2) ) / (2.0d0*alp(ixO^S))
    ! K_ji=K_ij
    do idir=1,2
       do jdir=idir+1,3
          K_ij(ixO^S,jdir,idir) = K_ij(ixO^S,idir,jdir)
       end do
    end do
    }

    ! Old
!    do idir=1,3
!       add_source_tau(ixO^S) = add_source_tau(ixO^S) &
!                  + Tmunu(ixO^S,0,0) * ( -betai(ixO^S,idir)*dalp(ixO^S,idir) ) &
!                  + Tmunu(ixO^S,0,idir) * ( -dalp(ixO^S,idir) )
!    enddo
!
!    do idir=1,3
!       do jdir=1,3
!          add_source_tau(ixO^S) = add_source_tau(ixO^S) &
!                    + Tmunu(ixO^S,0,idir)*( 2.0d0*betai(ixO^S,jdir)*K_ij(ixO^S,idir,jdir) ) &
!                    + Tmunu(ixO^S,0,0) * ( K_ij(ixO^S,idir,jdir)*betai(ixO^S,idir)*betai(ixO^S,jdir) ) &
!                    + Tmunu(ixO^S,idir,jdir) * K_ij(ixO^S,idir,jdir)
!       enddo
!    enddo
!
    ! New
    !do idir = 1, ndim
    !   add_source_tau(ixO^S) = add_source_tau(ixO^S) &
    !              + Tmunu(ixO^S,0,0) * ( -betai(ixO^S,idir)*dalp(ixO^S,idir) ) &
    !              + Tmunu(ixO^S,0,idir) * ( -dalp(ixO^S,idir) )
    !end do

    !do idir = 1, ndir; do jdir = 1, ndir
    !   add_source_tau(ixO^S) = add_source_tau(ixO^S) &
    !             + Tmunu(ixO^S,0,idir) * ( 2.0d0*betai(ixO^S,jdir)*K_ij(ixO^S,idir,jdir) ) &
    !             + Tmunu(ixO^S,0,0) * ( K_ij(ixO^S,idir,jdir)*betai(ixO^S,idir)*betai(ixO^S,jdir) )
    !end do; end do

    !do idir = 1, 3; do jdir = 1, 3
    !   add_source_tau(ixO^S) = add_source_tau(ixO^S) &
    !             + Tmunu(ixO^S,idir,jdir) * K_ij(ixO^S,idir,jdir)
    !end do; end do

  if (evolve_hydro) then

    ! wikdjgammaik is without the total pressure term
    wikdjgammaik(ixO^S,:) = zero
    do k = 1, ndir
    do j = 1, ndir
    do i = 1, ndir
       tmp(ixO^S) = sU(ixO^S,i)*vU(ixO^S,k) &
            !        + ptot(ixO^S)*myM%gammainv(i,k)%elem(ixO^S) &
            - bsU(ixO^S,i)*wCT(ixO^S,b0_+k)/lfac(ixO^S)
       wikdjgammaik(ixO^S,j) = wikdjgammaik(ixO^S,j) + tmp(ixO^S) * d_k_gamma_ij(ixO^S,i,k,j)
       !wikdjgammaik(ixO^S,j) = wikdjgammaik(ixO^S,j) + tmp(ixO^S) * myM%DgDk(i,k,j)%elem(ixO^S)
    end do
    end do
    end do

    sijkij = zero
    tmp(ixO^S) = zero
    do j = 1, ndir
    do i = 1, ndir
       tmp(ixO^S) = sU(ixO^S,i)*vU(ixO^S,j)  &
                     + ptot_mom(ixO^S) * gammainv(ixO^S,i,j) & 
                    -bsU(ixO^S,i) * wCT(ixO^S,b0_+j)/lfac(ixO^S)
       sijkij(ixO^S) = sijkij(ixO^S) + tmp(ixO^S) * K_ij(ixO^S,i,j)
    end do
    end do

!     tau_
            ! alp SijK_ij  - Sj dalpj
!code test
            tmp(ixO^S) = alp_prime(ixO^S) * sijkij(ixO^S) 
            !tmp(ixO^S) = alp(ixO^S) * sijkij(ixO^S) 
             do j = 1, ndir
             tmp(ixO^S) = tmp(ixO^S) - sU(ixO^S,j)* dalp(ixO^S,j)
             !tmp(ixO^S) = tmp(ixO^S) - sU(ixO^S,j)*myM%nonzeroDalphaDj(inonzero)%elem(ixO^S)
             enddo
          w(ixO^S,tau_) = w(ixO^S,tau_) + qdt*tmp(ixO^S)


!    s^C
      {^C& 
          ! s[s^C_] = 1/2 alpha W**ik d^Cdgammaik + S_i d^Cbeta**i -U d^Calpha

          tmp(ixO^S) = half * alp_prime(ixO^S) * wikdjgammaik(ixO^S,^C)
          !tmp(ixO^S) = half * alp(ixO^S) * wikdjgammaik(ixO^S,^C)
          !tmp(ixO^S) = half * myM%alpha(ixO^S) * wikdjgammaik(ixO^S,^C)
          ! Treat the total pressure separately to make discretized gradient:
          hxOmin^D=ixOmin^D-kr(^D,^C);hxOmax^D=ixOmax^D-kr(^D,^C);
          idims = ^C
          select case(idims)
             {case(^D)
             tmp(ixO^S) = tmp(ixO^S) + alp_prime(ixO^S) * ptot_mom(ixO^S) &
             !tmp(ixO^S) = tmp(ixO^S) + alp(ixO^S) * ptot_mom(ixO^S) &
             !tmp(ixO^S) = tmp(ixO^S) + myM%alpha(ixO^S)*ptot(ixO^S) &
                  *(mygeo%surfaceC^D(ixO^S)  - &
                    mygeo%surfaceC^D(hxO^S)  ) &
                  /mygeo%dvolume(ixO^S) \}
          end select

          do i = 1, ^NC
             tmp(ixO^S) = tmp(ixO^S) + wCT(ixO^S,s0_+i) * dbeta(ixO^S,i,^C)
             !tmp(ixO^S) = tmp(ixO^S) + wCT(ixO^S,s0_+i) * myM%dbetaiDj(i,^C)%elem(ixO^S)
          end do

             tmp(ixO^S) = tmp(ixO^S) - (wCT(ixO^S,tau_)+wCT(ixO^S,d_))* dalp(ixO^S,^C)
             !tmp(ixO^S) = tmp(ixO^S) - (wCT(ixO^S,tau_)+wCT(ixO^S,d_))*myM%dalphaDj(^C)%elem(ixO^S)


          w(ixO^S,s^C_) = w(ixO^S,s^C_) + qdt*tmp(ixO^S)



!          w(ixO^S,s^C_) = w(ixO^S,s^C_) + qdt * alp(ixO^S) * add_source_mom(ixO^S,^C)   ! from christoffel
          \}

    !      {#IFDEF GLM
    !      {^C&
    !                            ! s[B^C_] = -B**i/alpha dbetaid^C + beta**^C/alpha**2 * B**i * dalphadi
    !           case(b^C_)

    !      tmp(ixO^S) = zero
    !      do i = 1, ^NC
    !         if (dbetaidj_is_zero(i,^C)) cycle
    !         tmp(ixO^S) = tmp(ixO^S) - wCT(ixO^S,b0_+i) * myM%dbetaidj(i,^C)%elem(ixO^S)
    !      end do
    !      tmp(ixO^S) = tmp(ixO^S)/myM%alpha(ixO^S)

    !      if (.not. beta_is_zero(^C)) then
    !         tmp2(ixO^S) = zero
    !         do inonzero = 1, myM%nnonzeroDalphaDj
    !            i = myM%nonzeroDalphaDj(inonzero)%j
    !            tmp2(ixO^S) = tmp2(ixO^S) + wCT(ixO^S,b0_+i)*myM%nonzeroDalphaDj(inonzero)%elem(ixO^S)
    !         end do
    !         tmp(ixO^S) = tmp(ixO^S) + myM%beta(^C)%elem(ixO^S)/myM%alpha(ixO^S)**2 * tmp2(ixO^S)
    !      end if

    !      w(ixO^S,iw) = w(ixO^S,iw) + qdt*tmp(ixO^S)
    !      }

    !      ! s[psi_] = - phi DbetaiDi - 1/2 phi gammainv**ij * beta**k * DgijDk
    !   case(psi_)

    !      tmp(ixO^S) = zero
    !      do i = 1, ^NC
    !         if (dbetaidj_is_zero(i,i)) cycle
    !         tmp(ixO^S) = tmp(ixO^S) - myM%dbetaidj(i,i)%elem(ixO^S)
    !      end do
    !      tmp(ixO^S) = tmp(ixO^S) * wCT(ixO^S,psi_)

    !      tmp2(ixO^S) = zero
    !      do inonzero = 1, myM%nnonzeroDgDk
    !         i = myM%nonzeroDgDk(inonzero)%i
    !         j = myM%nonzeroDgDk(inonzero)%j
    !         k = myM%nonzeroDgDk(inonzero)%k
    !         if ( beta_is_zero(k) ) cycle

    !         tmp2(ixO^S) = tmp2(ixO^S) &
    !              + myM%beta(k)%elem(ixO^S) * myM%gammainv(i,j)%elem(ixO^S) * myM%DgDk(i,j,k)%elem(ixO^S)

    !      end do
    !      tmp(ixO^S) = tmp(ixO^S) - 0.5d0 * wCT(ixO^S,psi_) * tmp2(ixO^S)

    !      w(ixO^S,iw) = w(ixO^S,iw) + qdt*tmp(ixO^S)
    !      }
   endif
   
   ! Source terms for m1-radiation equations
   {#IFDEF M1
      call m1_add_geometrical_sources(ixI^L,ixO^L,x,wCT(ixI^S,1:nw),w, wold, qdt,d_k_gamma_ij,K_ij,qtC, metric_M1)
   } ! end DEF M1
	

  end subroutine addgeometry_source_full_gr


  subroutine addgeometry_source_cowling_cfc(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add geometrical source terms to w
    use mod_eos
    use mod_metric
    use mod_imhd_intermediate
    include 'amrvacdef.f'

    integer, intent(in)                       :: ixI^L, ixO^L
    double precision, intent(in)              :: qdt
    double precision, intent(in)              :: x(ixI^S,1:ndim)
    double precision, intent(inout)           :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    ! .. local ..
    integer                            :: iw, i,j,k, inonzero, jnonzero, hxO^L, idims, idir, jdir, kdir
    double precision,dimension(ixI^S,1:ndir)  :: vU, bsU, bsD, sU, wikdjgammaik
    double precision,dimension(ixI^S)         :: VdotB, sqrB, tmp, tmp2, wikbetajdjgammaik
    double precision,dimension(ixI^S)         :: ptot, lrho, lye, prs_tmp, eps_tmp
    integer :: ix^D
    double precision                          :: alp(ixI^S), temp_local, rho_local
    double precision                          :: betai(ixI^S,1:3)
    double precision                          :: psi(ixI^S), psi4(ixI^S), tryfunction(ixI^S), dtryfunction(ixI^S,1:ndim)

    ! dervitaves of the metric variables
    double precision                          :: dalp(ixI^S,1:3), dpsi(ixI^S,1:3)
    double precision                          :: dbeta(ixI^S,1:3,1:3)
    double precision                          :: d_k_gamma_ij(ixI^S,1:3,1:3,1:3)
    double precision                          :: d_gamma_ij_hat(ixI^S,1:3,1:3,1:3)
    ! 3-metric gamma^ij
    double precision                          :: gammainv(ixI^S,1:3,1:3)
    ! 3-metric gamma_ij
    double precision                          :: gamma(ixI^S,1:3,1:3)


    {#IFDEF DY_SP 
          call mpistop('why call cowling in DY_SP framework?')
    }


    VdotB(ixO^S)= ({^C& wCT(ixO^S,s^C_)*wCT(ixO^S,b^C_)+})/wCT(ixO^S,xi_)
    call square3u(ixI^L,ixO^L,myM,wCT(ixI^S,b1_:b^NC_),sqrB)
    call raise3(ixI^L,ixO^L,myM,wCT(ixI^S,s1_:s^NC_),sU)

    {^C&
         vU(ixO^S,^C)=(sU(ixO^S,^C)+VdotB(ixO^S)*wCT(ixO^S,b0_+^C))/ &
         (wCT(ixO^S,xi_)+sqrB(ixO^S))
    bsU(ixO^S,^C) = wCT(ixO^S,b0_+^C)/wCT(ixO^S,lfac_) + wCT(ixO^S,lfac_)*VdotB(ixO^S)*vU(ixO^S,^C)
    \}
    call lower3(ixI^L,ixO^L,myM,bsU,bsD)

    {do ix^D = ixO^LIM^D \}
        if (eos_type == tabulated) then
          temp_local = wCT(ix^D,T_eps_)
          rho_local  = wCT(ix^D,rho_)
          call eos_temp_get_all_one_grid(rho_local,temp_local,wCT(ix^D,ye_),&
                                         eps_tmp(ix^D),prs=prs_tmp(ix^D))
        else
          eps_tmp(ix^D) = wCT(ix^D, T_eps_)
          call eos_get_pressure_one_grid(prs_tmp(ix^D),wCT(ix^D,rho_),eps_tmp(ix^D))
        endif
    {enddo^D&\}

    ptot(ixO^S) = prs_tmp(ixO^S)
    ptot(ixO^S) = half*(VdotB(ixO^S)**2  &
         + sqrB(ixO^S)/wCT(ixO^S,lfac_)**2)  &
         + ptot(ixO^S)


   call imhd_get_intermediate_variables(ixI^L, ixI^L, wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                gammainv=gammainv(ixI^S,1:3,1:3))


    psi(ixI^S)       = wCT(ixI^S,psi_metric_)
    alp(ixI^S)       = wCT(ixI^S,alp_metric_)
{^C& betai(ixI^S,^C)  = wCT(ixI^S, beta_metric^C_)    \}
    psi4(ixI^S)      = psi(ixI^S)**4

    ! calculate derivitives of the metric variables
    dalp(ixI^S,1:3) = 0.0d0
    dbeta(ixI^S,1:3,1:3) = 0.0d0
    dpsi(ixI^S,1:3) = 0.0d0
    do idir = 1, ndim
       call partial_d( alp(ixI^S) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dalp(ixI^S,idir) )
       call partial_d( psi(ixI^S) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dpsi(ixI^S,idir) )
       dpsi(ixO^S,idir) = dpsi(ixO^S,idir) / psi(ixO^S)
       do jdir = 1, ndir
          call partial_d( betai(ixI^S,jdir) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dbeta(ixI^S,jdir,idir) )
       end do
    end do

!   write(*,*) dalp(ixO^S,1)
!   write(*,*) dalp(ixO^S,2)
!stop 'dalp12'

!    write(*,*)  dalp(:,5,1)
!  write(*,*) '------------------------'
!    write(*,*)  dalp(:,5,2)
!  write(*,*) '------------------------'
!    write(*,*)  dalp(:,5,3)
!stop

!code test
!   tryfunction(ixI^S) = x(ixI^S,1) **2
!    do idir = 1, ndim
!          call partial_d( tryfunction(ixI^S) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,dtryfunction(ixI^S,idir) )
!    end do
!
!
!  write(*,*)  tryfunction(:,5)
!
!  write(*,*) '------------------------'
!   write(*,*) dtryfunction(:,5,1) 
!  write(*,*) '------------------------'
!   write(*,*) dtryfunction(:,5,2) 
!  write(*,*) '------------------------'
!   write(*,*) dtryfunction(:,5,3) 
!  write(*,*) '------------------------'


    ! get the d_gamma_ij_hat in particular coordinate
       call partial_d_gamma_ij_hat( x(ixI^S, 1:ndim), ixI^L, ixI^L, d_gamma_ij_hat)

    ! partial dervitaves of metric d_i gamma_jk
    d_k_gamma_ij(ixI^S,1:3,1:3,1:3) = 0.0d0
!    do kdir = 1, ndim
!     do jdir = 1, ndim
!      do idir = 1, ndim
!          d_k_gamma_ij(ixO^S,idir,jdir,kdir) = 4.0d0 * gamma(ixO^S,idir,jdir) * dpsi(ixO^S,kdir) +&
!                                          psi4(ixO^S) * d_gamma_ij_hat(ixO^S,idir,jdir,kdir)
!      end do
!     end do
!    end do
    do kdir = 1, ndim
     do jdir = 1, ndir
      do idir = 1, ndir
          d_k_gamma_ij(ixO^S,idir,jdir,kdir) = 4.0d0 * gamma(ixO^S,idir,jdir) * dpsi(ixO^S,kdir) +&
                                          psi4(ixO^S) * d_gamma_ij_hat(ixO^S,idir,jdir,kdir)
      end do
     end do
    end do



  if (evolve_hydro) then
    ! wikdjgammaik is without the total pressure term
    wikdjgammaik(ixO^S,:) = zero
    do k = 1, ndir
    do j = 1, ndir
    do i = 1, ndir
       tmp(ixO^S) = sU(ixO^S,i)*vU(ixO^S,k) &
            !        + ptot(ixO^S)*myM%gammainv(i,k)%elem(ixO^S) &
            - bsU(ixO^S,i)*wCT(ixO^S,b0_+k)/wCT(ixO^S,lfac_)
       wikdjgammaik(ixO^S,j) = wikdjgammaik(ixO^S,j) + tmp(ixO^S) * d_k_gamma_ij(ixO^S,i,k,j)
       !wikdjgammaik(ixO^S,j) = wikdjgammaik(ixO^S,j) + tmp(ixO^S) * myM%DgDk(i,k,j)%elem(ixO^S)
    end do
    end do
    end do


    wikbetajdjgammaik(ixO^S) = zero
    do k = 1, ndir
    do j = 1, ndir
    do i = 1, ndir
       tmp(ixO^S) = sU(ixO^S,i)*vU(ixO^S,k) &
            + ptot(ixO^S)* gammainv(ixO^S,i,k) &
            !+ ptot(ixO^S)*myM%gammainv(i,k)%elem(ixO^S) &
            - bsU(ixO^S,i)*wCT(ixO^S,b0_+k)/wCT(ixO^S,lfac_)

       wikbetajdjgammaik(ixO^S) = wikbetajdjgammaik(ixO^S) &
            + tmp(ixO^S) * betai(ixO^S,j) * d_k_gamma_ij(ixO^S,i,k,j)
            !+ tmp(ixO^S) * myM%beta(j)%elem(ixO^S) * myM%DgDk(i,k,j)%elem(ixO^S)
    end do
    end do
    end do


!          s^C
      {^C& 
                                ! s[s^C_] = 1/2 alpha W**ik d^Cdgammaik + S_i d^Cbeta**i -U d^Calpha
          tmp(ixO^S) = half * alp(ixO^S) * wikdjgammaik(ixO^S,^C)
          !tmp(ixO^S) = half * myM%alpha(ixO^S) * wikdjgammaik(ixO^S,^C)
          ! Treat the total pressure separately to make discretized gradient:
          hxOmin^D=ixOmin^D-kr(^D,^C);hxOmax^D=ixOmax^D-kr(^D,^C);
          idims = ^C
          select case(idims)
             {case(^D)
             tmp(ixO^S) = tmp(ixO^S) + alp(ixO^S) *ptot(ixO^S) &
             !tmp(ixO^S) = tmp(ixO^S) + myM%alpha(ixO^S)*ptot(ixO^S) &
                  *(mygeo%surfaceC^D(ixO^S)-mygeo%surfaceC^D(hxO^S)) &
                  /mygeo%dvolume(ixO^S)\}
          end select

          do i = 1, ^NC
             tmp(ixO^S) = tmp(ixO^S) + wCT(ixO^S,s0_+i) * dbeta(ixO^S,i,^C)
             !tmp(ixO^S) = tmp(ixO^S) + wCT(ixO^S,s0_+i) * myM%dbetaiDj(i,^C)%elem(ixO^S)
          end do

             tmp(ixO^S) = tmp(ixO^S) - (wCT(ixO^S,tau_)+wCT(ixO^S,d_))* dalp(ixO^S,^C)
             !tmp(ixO^S) = tmp(ixO^S) - (wCT(ixO^S,tau_)+wCT(ixO^S,d_))*myM%dalphaDj(^C)%elem(ixO^S)


          w(ixO^S,s^C_) = w(ixO^S,s^C_) + qdt*tmp(ixO^S)

          \}
     

 
!           tau

 
          ! s[tau_] = 1/2 W**ik beta**j dgammaikdj + W**j_i*dbetaidj - S**j dalphadj
          tmp(ixO^S) = half*wikbetajdjgammaik(ixO^S)

!not sure  fixme
          do i = 1, ndir
          do j = 1, ndir
             tmp2(ixO^S) = vU(ixO^S,j)*wCT(ixO^S,s0_+i)
             if (i .eq. j) &
                  tmp2(ixO^S) = tmp2(ixO^S) + ptot(ixO^S)
             tmp2(ixO^S) = tmp2(ixO^S) - bsD(ixO^S,i)*wCT(ixO^S,b0_+j)/wCT(ixO^S,lfac_)
             !            print*, 'Source:',i,j,maxval(abs(tmp2(ixO^S)))
             tmp(ixO^S) = tmp(ixO^S) + tmp2(ixO^S)  * dbeta(ixO^S,i,j)
             !tmp(ixO^S) = tmp(ixO^S) + tmp2(ixO^S)  * myM%nonzeroDbetaiDj(inonzero)%elem(ixO^S)
          end do
          end do

          do j = 1, ndir
             tmp(ixO^S) = tmp(ixO^S) - sU(ixO^S,j)* dalp(ixO^S,j)
             !tmp(ixO^S) = tmp(ixO^S) - sU(ixO^S,j)*myM%nonzeroDalphaDj(inonzero)%elem(ixO^S)
          enddo
          w(ixO^S,tau_) = w(ixO^S,tau_) + qdt*tmp(ixO^S)


!          {#IFDEF GLM
!          {^C&
!                                ! s[B^C_] = -B**i/alpha dbetaid^C + beta**^C/alpha**2 * B**i * dalphadi
!               case(b^C_)
!
!          tmp(ixO^S) = zero
!          do i = 1, ^NC
!             if (dbetaidj_is_zero(i,^C)) cycle
!             tmp(ixO^S) = tmp(ixO^S) - wCT(ixO^S,b0_+i) * myM%dbetaidj(i,^C)%elem(ixO^S)
!          end do
!          tmp(ixO^S) = tmp(ixO^S)/myM%alpha(ixO^S)
!
!          if (.not. beta_is_zero(^C)) then
!             tmp2(ixO^S) = zero
!             do inonzero = 1, myM%nnonzeroDalphaDj
!                i = myM%nonzeroDalphaDj(inonzero)%j
!                tmp2(ixO^S) = tmp2(ixO^S) + wCT(ixO^S,b0_+i)*myM%nonzeroDalphaDj(inonzero)%elem(ixO^S)
!             end do
!             tmp(ixO^S) = tmp(ixO^S) + myM%beta(^C)%elem(ixO^S)/myM%alpha(ixO^S)**2 * tmp2(ixO^S)
!          end if
!
!          w(ixO^S,iw) = w(ixO^S,iw) + qdt*tmp(ixO^S)
!          \}
!
!          ! s[psi_] = - phi DbetaiDi - 1/2 phi gammainv**ij * beta**k * DgijDk
!       case(psi_)
!
!          tmp(ixO^S) = zero
!          do i = 1, ^NC
!             if (dbetaidj_is_zero(i,i)) cycle
!             tmp(ixO^S) = tmp(ixO^S) - myM%dbetaidj(i,i)%elem(ixO^S)
!          end do
!          tmp(ixO^S) = tmp(ixO^S) * wCT(ixO^S,psi_)
!
!          tmp2(ixO^S) = zero
!          do inonzero = 1, myM%nnonzeroDgDk
!             i = myM%nonzeroDgDk(inonzero)%i
!             j = myM%nonzeroDgDk(inonzero)%j
!             k = myM%nonzeroDgDk(inonzero)%k
!             if ( beta_is_zero(k) ) cycle
!
!             tmp2(ixO^S) = tmp2(ixO^S) &
!                  + myM%beta(k)%elem(ixO^S) * myM%gammainv(i,j)%elem(ixO^S) * myM%DgDk(i,j,k)%elem(ixO^S)
!
!          end do
!          tmp(ixO^S) = tmp(ixO^S) - 0.5d0 * wCT(ixO^S,psi_) * tmp2(ixO^S)
!
!          w(ixO^S,iw) = w(ixO^S,iw) + qdt*tmp(ixO^S)
!          }
   endif

  end subroutine addgeometry_source_cowling_cfc

   

end module mod_full_gr_source
