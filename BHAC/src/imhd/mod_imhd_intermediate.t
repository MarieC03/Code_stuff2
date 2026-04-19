module mod_imhd_intermediate
  implicit none
  public

  ! entire module is for DY_SP dynamical spacetime (CFC, BSSN)

  contains


  subroutine imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, &
             gamma, gammainv, lfac2, lfac, v2, v_hat, b2, B_dot_v, Bvec2, bmu, Ptot, h_th, htot, P_th, eps, u4, ub_bernoulli)
    use mod_eos
    use mod_metric
    use mod_cfc_parameters
    include 'amrvacdef.f'
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(in)            :: w(ixI^S, 1:nw)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)

    double precision, intent(out), optional :: B_dot_v(ixI^S)
    double precision, intent(out), optional :: v_hat(ixI^S,1:ndir)
    double precision, intent(out), optional :: bmu(ixI^S,0:ndir)   ! projection of B^mu along fluid four-velocity u^nu
    double precision, intent(out), optional :: b2(ixI^S)        
    double precision, intent(out), optional :: Bvec2(ixI^S)        
    double precision, intent(out), optional :: gamma(ixI^S,1:3,1:3)
    double precision, intent(out), optional :: gammainv(ixI^S,1:3,1:3)
    double precision, intent(out), optional :: v2(ixI^S)        
    double precision, intent(out), optional :: lfac2(ixI^S) ! Lorentz factor square
    double precision, intent(out), optional :: lfac(ixI^S)  ! Lorentz factor
    double precision, intent(out), optional :: h_th(ixI^S)     ! thermal enthalpy: h = 1 + eps + p/rho
    double precision, intent(out), optional :: htot(ixI^S)  ! modified enthalpy:(h + b2/rho)
    double precision, intent(out), optional :: Ptot(ixI^S)  ! total pressure
    double precision, intent(out), optional :: P_th(ixI^S)  ! thermal pressure
    double precision, intent(out), optional :: eps(ixI^S)  ! specific internal energy
    double precision, intent(out), optional :: u4(ixI^S,0:ndir)
    double precision, intent(out), optional :: ub_bernoulli(ixI^S) ! Bernoulli criteria for unbound mass

    integer                                 :: idir, ix^D
    double precision                        :: B_dot_v_tmp(ixI^S)
    double precision                        :: bmu_tmp(ixI^S,0:ndir)   ! projection of B^mu along fluid four-velocity u^nu
    double precision                        :: b2_tmp(ixI^S), htot_tmp(ixI^S)        
    double precision                        :: prs_tmp(ixI^S), eps_tmp(ixI^S)
    double precision                        :: W2v2(ixI^S), u_t(ixI^S)        
    double precision                        :: vU(ixI^S,1:ndir) 
    double precision                        :: bD(ixI^S,1:ndir) 
    double precision                        :: v_hat_tmp(ixI^S,1:ndir)
    double precision                        :: lfac2_tmp(ixI^S) ! Lorentz factor square
    double precision                        :: lfac_tmp(ixI^S) ! Lorentz factor
    double precision                        :: gamma_tmp(ixI^S,1:3,1:3), alp_prime(ixI^S)
    double precision                        :: gammainv_tmp(ixI^S,1:3,1:3)
    double precision                        :: temp_local, rho_local

  

    {#IFNDEF DY_SP
       call mpistop('non cfc framework should not come here')
    }

    if ( .not. useprimitiveRel ) then
      call mpistop('current get intermediate only support useprimitiveRel')
    endif

    ! get the metric
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma_tmp(ixO^S,idir,idir) = gamma_tmp(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
    end do
    if ( present(gamma) ) then
       gamma(ixO^S,1:3,1:3) = gamma_tmp(ixO^S,1:3,1:3)
    end if

    ! get the metric
    call get_gammainvij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gammainv_tmp(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gammainv_tmp(ixO^S,idir,idir) = gammainv_tmp(ixO^S,idir,idir) / w(ixO^S, psi_metric_)**4
    end do
    if ( present(gammainv) ) then
       gammainv(ixO^S,1:3,1:3) = gammainv_tmp(ixO^S,1:3,1:3)
    end if

    W2v2 = 0.0d0
    ! calculate W^2 * v^2 first
    W2v2(ixO^S) = {^C& gamma_tmp(ixO^S,^C,^C)*w(ixO^S, u^C_)**2 +}

    ! Calculate the Lorentz factor from velocity
    lfac2_tmp(ixO^S) = W2v2(ixO^S) + 1.0d0
    lfac_tmp(ixO^S) = dsqrt( lfac2_tmp(ixO^S) )

    if ( present(lfac2) ) lfac2(ixO^S) = lfac2_tmp(ixO^S)
    if ( present(lfac) ) lfac(ixO^S) = lfac_tmp(ixO^S)
    if ( present(v2) ) v2(ixO^S) = W2v2(ixO^S) / lfac2_tmp(ixO^S)

    if (use_gw_br .and. use_h00_coupling) then
       alp_prime(ixO^S) = dsqrt( w(ixO^S,alp_metric_)**2 - w(ixO^S,h_00_) )
       !alp_prime(ixO^S) = w(ixO^S,alp_metric_) - w(ixO^S,h_00_)/2.0d0/w(ixO^S,alp_metric_)
    else
       alp_prime(ixO^S) = w(ixO^S,alp_metric_)
    endif

    ! v_hat^i = v^i * alp - beta
    {^C&   v_hat_tmp(ixO^S, ^C) = alp_prime(ixO^S) * w(ixO^S, u^C_) / lfac_tmp(ixO^S) - w(ixO^S, beta_metric^C_)  \}
    if ( present(v_hat) ) v_hat(ixO^S,1:ndir) = v_hat_tmp(ixO^S,1:ndir)

    B_dot_v_tmp(ixO^S) = 0.0d0

    B_dot_v_tmp(ixO^S) = {^C& w(ixO^S, b^C_) * gamma_tmp(ixO^S,^C,^C) &
               * w(ixO^S, u^C_) / lfac_tmp(ixO^S)  +}

    if ( present(B_dot_v) ) B_dot_v(ixO^S) = B_dot_v_tmp(ixO^S)

    b2_tmp(ixO^S) = 0.0d0

    b2_tmp(ixO^S) = {^C& w(ixO^S, B^C_)**2 * gamma_tmp(ixO^S,^C,^C) +}

    
    if (present(Bvec2))  Bvec2(ixO^S) = b2_tmp(ixO^S)

    b2_tmp(ixO^S) = b2_tmp(ixO^S) / lfac2_tmp(ixO^S) + B_dot_v_tmp(ixO^S)**2
    if (present(b2)) b2(ixO^S) = b2_tmp(ixO^S)

    if ( present(bmu) ) then
       ! Calculate the projection of B^mu along fluid four-velocity u^nu
       bmu(ixO^S,0)  = B_dot_v_tmp(ixO^S) * lfac_tmp(ixO^S) / alp_prime(ixO^S)

  {^C& bmu(ixO^S,^C) = w(ixO^S, b^C_) / lfac_tmp(ixO^S) + bmu(ixO^S,0) * v_hat_tmp(ixO^S, ^C) \}
    end if
 
    if ( present(Ptot) .or. present(h_th) .or. present(htot) .or. present(P_th) .or. present(eps) .or. present(ub_bernoulli) ) then
       {do ix^D = ixO^LIM^D \}
         if (eos_type == tabulated) then
           temp_local = w(ix^D, T_eps_)
           rho_local  = w(ix^D, rho_)
           call eos_temp_get_all_one_grid(rho_local,temp_local,w(ix^D,ye_),&
                                          eps_tmp(ix^D),prs=prs_tmp(ix^D))
         else
           eps_tmp(ix^D) = w(ix^D, T_eps_)
           call eos_get_pressure_one_grid(prs_tmp(ix^D),w(ix^D,rho_),eps_tmp(ix^D))
         endif
       {enddo^D&\}
    endif

    if ( present(Ptot) ) then
       ! Calculate the total pressure
       Ptot(ixO^S) = prs_tmp(ixO^S) + 0.5d0 * b2_tmp(ixO^S)
    end if

    if ( present(h_th) ) then
       ! Calculate the thermal specific enthalpy
       h_th(ixO^S) = 1.0d0 + eps_tmp(ixO^S) &
                + ( prs_tmp(ixO^S) ) / w(ixO^S, rho_)
    end if

    if ( present(htot)) then
       ! Calculate the magnetically modified specific enthalpy
       htot(ixO^S) = 1.0d0 + eps_tmp(ixO^S) &
                + ( prs_tmp(ixO^S) + b2_tmp(ixO^S) ) / w(ixO^S, rho_)
    end if

    if ( present(P_th) ) P_th(ixO^S) = prs_tmp(ixO^S)
    if ( present(eps) ) eps(ixO^S) = eps_tmp(ixO^S)

    ! it is vector not 1-form
    if ( present(u4) ) then
      u4(ixO^S, 0) = lfac_tmp(ixO^S)/ alp_prime(ixO^S)
      {^C& u4(ixO^S, ^C) = (w(ixO^S, u^C_)/lfac_tmp(ixO^S) - w(ixO^S, beta_metric^C_)/alp_prime(ixO^S)) * lfac_tmp(ixO^S) \}
    endif
    if ( present(ub_bernoulli) ) then
        htot_tmp(ixO^S) = 1.0d0+eps_tmp(ixO^S)+(prs_tmp(ixO^S))/w(ixO^S, rho_)
        u_t(ixO^S) = -lfac_tmp(ixO^S)*w(ixO^S, alp_metric_)
        {^C& u_t(ixO^S) = u_t(ixO^S)+w(ixO^S, psi_metric_)**4*w(ixO^S, beta_metric^C_)*w(ixO^S, u^C_) \}
        ub_bernoulli(ixO^S) = -htot_tmp(ixO^S)*u_t(ixO^S)-eos_hmin
    endif

  end subroutine imhd_get_intermediate_variables

  ! you should have updated u^C
  subroutine dysp_get_lfac(ixI^L, ixO^L, w, x, lfac)
    include 'amrvacdef.f'
    integer, intent(in)                     :: ixI^L, ixO^L
    double precision, intent(in)            :: w(ixI^S, 1:nw)
    double precision, intent(in)            :: x(ixI^S, 1:ndim)
    double precision, intent(out)           :: lfac(ixI^S)
  
    double precision                        :: W2v2(ixI^S), gamma_tmp(ixI^S,1:3,1:3)
    integer                                 :: idir

    {#IFNDEF DY_SP
       call mpistop('non cfc framework should not come here')
    }

    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma_tmp(ixO^S,idir,idir) = gamma_tmp(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
    end do

    W2v2(ixO^S) = {^C& gamma_tmp(ixO^S,^C,^C)*w(ixO^S, u^C_)**2 +}
    lfac(ixO^S) = dsqrt( W2v2(ixO^S) + 1.0d0 )

  end subroutine dysp_get_lfac

  subroutine imhd_get_Efield_Eulerian(ixI^L, ixO^L, w, x, Ei)
      include 'amrvacdef.f'
      integer, intent(in)                         :: ixI^L, ixO^L
      double precision, intent(in)                :: w(ixI^S, 1:nw)
      double precision, intent(in)                :: x(ixI^S, 1:ndim)
      double precision, intent(out)               :: Ei(ixI^S, 1:ndir)        ! contravariant 1 vector in 3 space
      double precision, dimension(ixI^S, 1:ndim)  :: E_i(ixI^S, 1:ndir)    ! covariant 1 form in 3 space
      integer                                     :: idir1, idir2, idir
      double precision                            :: gamma_tmp(ixI^S,1:3,1:3), lfac(ixI^S), sqrt_gamma(ixI^S), W2v2(ixI^S)

      call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_tmp(ixI^S,1:3,1:3))
      do idir = 1, ndir
         gamma_tmp(ixO^S,idir,idir) = gamma_tmp(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
      end do
  
      W2v2(ixO^S) = {^C& gamma_tmp(ixO^S,^C,^C)*w(ixO^S, u^C_)**2 +}
      lfac(ixO^S) = dsqrt( W2v2(ixO^S) + 1.0d0 )

      call get_sqrt_gamma_hat(x, ixI^L, ixO^L, sqrt_gamma)

      sqrt_gamma(ixO^S) = sqrt_gamma(ixO^S) * w(ixO^S, psi_metric_)**6

      E_i(ixI^S, 1:ndim) = zero

      do idir1=1, ndir
          do idir2=1, ndir
              do idir=1, ndir
                  E_i(ixO^S, idir) = E_i(ixO^S, idir)+ sqrt_gamma(ixO^S) * lvc(idir, idir1, idir2)* &
                      w(ixO^S, b0_+idir1)*w(ixO^S, u0_+idir2)/lfac(ixO^S)
              enddo
          enddo
      enddo

      do idir = 1, ^NC
        Ei(ixO^S, idir) = E_i(ixO^S, idir) / gamma_tmp(ixO^S, idir, idir)
      enddo
  end subroutine imhd_get_Efield_Eulerian

  !> This subroutine fix the abnormal values in primitive variables !
  subroutine imhd_handle_small_values(w, x, ixI^L, ixO^L, update_eos)
    use mod_eos
    use mod_cfc_parameters
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    logical, intent(in)             :: update_eos

    integer                         :: idir
    integer                         :: ix^D
    double precision                :: eps_min, eps_max, eps_tmp
    logical                         :: need_to_fix_eos

    ! avoid coordinate singularities
    ! fixme: this might depends on different bc, but in general this should work
{#IFDEF DY_SP
    if ( coordinate /= cartesian ) then
       where ( dabs(x(ixO^S,1)) < smalldouble )
          w(ixO^S, u1_)           = 0.0d0
          w(ixO^S, s1_)           = 0.0d0
          !w(ixO^S, beta_metric1_) = 0.0d0
       end where
    end if
}
    select case (small_values_method)
    case ("replace")
       ! check the prim variables one by one
       {do ix^D = ixO^LIM^D \}
          need_to_fix_eos = .false.
          if ( w(ix^D, rho_) < small_rho_thr ) then
             ! atmosphere handling
             w(ix^D, rho_)   = small_rho 
             if (eos_type == tabulated) then
               w(ix^D, T_eps_) = small_temp
               w(ix^D, ye_)    = big_ye    ! depends on input data
             else
               w(ix^D, T_eps_) = small_eps
             endif
         {^C& w(ix^D, u^C_)   = 0.0d0 \}
             {#IFNDEF DY_SP
               w(ix^D, pp_)    = small_press
               w(ix^D, lfac_)  = 1.0d0
               ! fixme: add small_xi
             }
          else
             ! this is not atmosphere
             if (eos_type == tabulated) then
               if ( ( w(ix^D, ye_) < eos_yemin ) .or. ( w(ix^D, ye_) > eos_yemax ) ) then
                  w(ix^D, ye_) = max(eos_yemin, min(eos_yemax, w(ix^D, ye_)))
                  need_to_fix_eos = .True.
               end if

               if ( ( w(ix^D, T_eps_) < eos_tempmin ) .or. ( w(ix^D, T_eps_) > eos_tempmax ) ) then
                  w(ix^D, T_eps_) = max(eos_tempmin, min(eos_tempmax, w(ix^D, T_eps_)))
                  need_to_fix_eos = .True.
               end if 
             endif

             if ( need_to_fix_eos .and. update_eos) then
                call Eos_update_one_point(w(ix^D,1:nw))
             end if
          end if
       {enddo^D&\}
    case default
       ! nothing to do here
       !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       return
    end select
  end subroutine imhd_handle_small_values

  ! it will be called in amrvacusr.t and reconstruction
  subroutine Eos_update_one_grid(ixI^L, ixO^L, w)
    include 'amrvacdef.f'
    integer, intent(in)                       :: ixI^L, ixO^L
    double precision, intent(inout)           :: w(ixI^S, 1:nw)
    integer                                   :: ix^D
    {do ix^D = ixO^LIM^D \}
       call Eos_update_one_point(w(ix^D,1:nw))
    {end do^D&\}
  end subroutine Eos_update_one_grid

  subroutine Eos_update_one_point(w)
    use mod_eos
    include 'amrvacdef.f'
    double precision, intent(inout)           :: w(1:nw)
    double precision                          :: eps_tmp, dummy

    {#IFDEF DY_SP
    ! nothing, becoz there is no pp_, eps_, cs2_ anymore
    }
    {#IFNDEF DY_SP
    ! it has pp_ only
    if (eos_type == tabulated) then
         call eos_temp_get_all_one_grid(w(rho_),w(T_eps_),w(ye_),&
                                        dummy, prs=w(pp_))
    else
       eps_tmp = w(T_eps_)
       call eos_get_pressure_one_grid(w(pp_),w(rho_),eps_tmp)
    end if
    }

  end subroutine Eos_update_one_point

  subroutine imhd_modify_wLR(ixI^L,ixO^L,qt,wL,wR,s,idir)
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wL(ixI^S,1:nw), wR(ixI^S,1:nw)
    type(state)                     :: s

{#IFDEF STAGGERED
    associate(w=>s%w%w, ws=>s%ws%w, mygeo=>s%geo)
    {^C&
      if (^C == idir) then
        wL(ixO^S,b^C_) = ws(ixO^S, idir)
        wR(ixO^S,b^C_) = ws(ixO^S, idir)
      endif
    }
    end associate
}
  end subroutine imhd_modify_wLR

  ! get psi^6 U, which is part of the source terms in cfc solver
  subroutine imhd_get_tilde_U(ixI^L,ixO^L,w,x,tilde_U)
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(out)   :: tilde_U(ixI^S)
    double precision :: Tmunu_m1 = 0.0d0
    ! cons are * psi**6
    {#IFNDEF M1
    tilde_U(ixO^S) = w(ixO^S, tau_) + w(ixO^S, D_)
    }
    {#IFDEF M1
    tilde_U(ixO^S) = w(ixO^S, tau_) + w(ixO^S, D_) !+ Tmunu_m1
    }
  end subroutine imhd_get_tilde_U

  ! get psi^6 S_i, which is part of the source terms in cfc solver
  subroutine imhd_get_tilde_S_i(ixI^L,ixO^L,w,x,tilde_S_i)
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(out)   :: tilde_S_i(ixI^S, 1:3)

    integer                         :: idir
    double precision :: Tmunu_m1 = 0.0d0

    ! cons are * psi**6
    tilde_S_i = 0.0d0
    {#IFNDEF M1
    {^C&   tilde_S_i(ixO^S, ^C) = w(ixO^S, s^C_)  \}
    }
    {#IFDEF M1
    {^C&   tilde_S_i(ixO^S, ^C) = w(ixO^S, s^C_) \} ! +Tmunu_m1
    }
  end subroutine imhd_get_tilde_S_i

  ! get psi^6 (gamma_{ij}S^{ij}), which is part of the source terms in cfc_alp
  subroutine imhd_get_tilde_S(ixI^L,ixO^L,w,x,tilde_S)
    include 'amrvacdef.f'

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(out)   :: tilde_S(ixI^S)

    double precision                :: B_dot_v(ixI^S)
    double precision                :: b2(ixI^S)
    double precision                :: htot(ixI^S)    ! modified enthalpy:(h + b2/rho)
    double precision                :: Ptot(ixI^S) ! total pressure
    double precision                :: lfac(ixI^S) ! Lorentz factor
    double precision                :: gamma(ixI^S,1:3,1:3)
    integer                         :: idir
  
    ! cons are * psi**6

    call imhd_get_intermediate_variables(ixI^L, ixO^L, w(ixI^S, 1:nw), x(ixI^S, 1:ndim), &
                gamma=gamma(ixI^S,1:3,1:3), &
                lfac=lfac(ixI^S), &
                B_dot_v=B_dot_v(ixI^S), b2=b2(ixI^S), &
                Ptot=Ptot(ixI^S), &
                htot=htot(ixI^S) )

    ! nota that S = S_i v^i + 3 p* - b^2
    tilde_S(ixO^S) = 3.0d0 * Ptot(ixO^S) - b2(ixO^S)
    tilde_S(ixO^S) = tilde_S(ixO^S) * w(ixO^S, psi_metric_)**6
     {^C&
       tilde_S(ixO^S) = tilde_S(ixO^S) &
        + w(ixO^S, s^C_) * w(ixO^S, u^C_) / lfac(ixO^S)
     \}
     !{#IFDEF M1
      !... M1_TODO
      !}
  end subroutine imhd_get_tilde_S

! all including ghostzones
  subroutine conformal_transformation(ixI^L,ixO^L,w,x)
    use mod_eos
    include 'amrvacdef.f'
    integer, intent(in)               :: ixI^L, ixO^L
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    double precision, intent(inout)   :: w(ixI^S,1:nw)

{^C& w(ixO^S, b^C_) = w(ixO^S, b^C_) * w(ixO^S, psi_metric_)**6    \}
     w(ixO^S, D_)   = w(ixO^S, D_) * w(ixO^S, psi_metric_)**6
     w(ixO^S, tau_) = w(ixO^S, tau_) * w(ixO^S, psi_metric_)**6
{^C& w(ixO^S, s^C_) = w(ixO^S, s^C_) * w(ixO^S, psi_metric_)**6    \}
    if (eos_type == tabulated) w(ixO^S, Dye_) = w(ixO^S, Dye_) * w(ixO^S, psi_metric_)**6

  end subroutine conformal_transformation

  subroutine inverse_conformal_trans(ixI^L,ixO^L,w,x)
    use mod_eos
    include 'amrvacdef.f'
    integer, intent(in)               :: ixI^L, ixO^L
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    double precision, intent(inout)   :: w(ixI^S,1:nw)

{^C& w(ixO^S, b^C_) = w(ixO^S, b^C_) / w(ixO^S, psi_metric_)**6    \}
     w(ixO^S, D_)   = w(ixO^S, D_) / w(ixO^S, psi_metric_)**6
     w(ixO^S, tau_) = w(ixO^S, tau_) / w(ixO^S, psi_metric_)**6
{^C& w(ixO^S, s^C_) = w(ixO^S, s^C_) / w(ixO^S, psi_metric_)**6    \}
    if (eos_type == tabulated) w(ixO^S, Dye_) = w(ixO^S, Dye_) / w(ixO^S, psi_metric_)**6

  end subroutine inverse_conformal_trans

!ppm related
  subroutine imhd_ppm_flatcd(ixI^L, ixO^L, ixL^L, ixR^L, w, d2w, drho, dp)
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, 1:nw),d2w(ixI^S,1:nw)
    double precision, intent(inout) :: drho(ixI^S), dp(ixI^S)

    drho(ixO^S) = abs(d2w(ixO^S,rho_)) &
                  / min( w(ixL^S, rho_), w(ixR^S, rho_) )

    !dp(ixO^S) = abs(d2w(ixO^S,pp_)) &
    !              / min( w(ixL^S, pp_), w(ixR^S, pp_) )
  end subroutine imhd_ppm_flatcd

  subroutine imhd_ppm_flatsh(ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L, idims, w, beta, z, dv)
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: beta(ixI^S), z(ixI^S), dv(ixI^S)

    ! beta here is the shock width
    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    !where ( abs( w(ixRR^S, pp_) - w(ixLL^S, pp_) )> smalldouble )
    !   beta(ixO^S) = abs( w(ixR^S, pp_) - w(ixL^S, pp_) ) &
    !                 /abs( w(ixRR^S, pp_) - w(ixLL^S, pp_) )
    !else where
    !   beta(ixO^S) = 0.0d0
    !end where

    ! Z here is the shock strength
    ! in eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
    ! where the denominator is rho * cs2, which is gamma * P in ideal gas.
    ! However, in eq. B16, page 218, Mignone and Bodo 2005, ApJS (beta1),
    ! this should be min(p(i-1), p(i+1)).
    ! Here we follow the later approach, as usually we don't have gamma
    ! in general
    !z(ixO^S) = abs( w(ixR^S, pp_) - w(ixL^S, pp_) ) &
    !                 /min( w(ixR^S, pp_), w(ixL^S, pp_) )

    {^C&
      if (^C == idims) then
        dv(ixO^S)=w(ixR^S,u^C_) - w(ixL^S,u^C_)
      endif
    \}
    

  end subroutine imhd_ppm_flatsh

end module mod_imhd_intermediate
