!0. the whole con2prim is using eps instead of xi, with eps --> calculate the xi/prs/etc.
!1. enthlpay output is specific thermo enthplay, not the total enthplay including magnetic pressure
module mod_imhd_con2prim
  implicit none
  save
  double precision                :: v_max
  double precision                :: h0           ! lower bound for the relativistic enthalpy over entire validity region of eos
  double precision                :: one_over_h0  ! 1 / h0
  integer                         :: iter_max = 50000
  double precision                :: Kc2p_tolerance = 1.0d-15
  logical                         :: usr_W_limit_scheme = .false.

contains

  subroutine imhd_con2prim_setup
!###############################################################################################################!
     ! nothing need to set the typeinversion. As this con2prim can do all the EOSs with all hydro variables
!###############################################################################################################!
    use mod_eos, only: eos_hmin
    include 'amrvacdef.f'
    if (lfac_max .le. 1) call mpistop('lfac_max has not been setup, set it in para file')
    ! initialize 1/h0
    h0 = eos_hmin 
    one_over_h0 = 1.0d0 / h0
    ! initialize the k_max and v_max
    v_max = dsqrt(1.0d0 - 1.0d0 / lfac_max**2)
  end subroutine imhd_con2prim_setup

  subroutine imhd_con2prim(ixI^L, ixO^L, x, w, patchw, need_aux_only)
    use mod_metric, only: raise3_ixD, lower3_ixD, square3u_ixD
    use mod_eos
    use mod_rootfinding
    use mod_small_values
    use mod_imhd_intermediate
    include 'amrvacdef.f'
    
    logical, intent(in)            :: need_aux_only
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(inout):: w(ixI^S,1:nw)
    double precision, intent(in)   :: x(ixI^S,1:ndim)
    logical, intent(in),dimension(ixG^T)   :: patchw

!    ! .. local ..
    double precision                :: w_nan(ixI^S,1:nw)
    integer                         :: idir, flag(ixI^S)
    double precision                :: psi6(ixI^S)
    double precision                :: cons_tmp(ixI^S,1:nw)
    double precision                :: cons_tmp_new(ixI^S,1:nw)
    double precision                :: cons_tmp_adjust(ixI^S,1:nw)
    double precision                :: cons_sC_tmp_U_adjust(ixI^S,1:ndir)
    double precision                :: cons_sC_tmp_D_adjust(ixI^S,1:ndir)
    double precision                :: cons_sC_tmp_new_U(ixI^S,1:ndir)
    double precision                :: cons_sC_tmp_new_D(ixI^S,1:ndir)
    integer                         :: ix^D, i, j 
    double precision                :: gamma(ixI^S,1:3,1:3) ! metric gamma_ij
    ! rescaled variables
    double precision                :: q
    double precision                :: r_i(1:ndir) !r_i
    double precision                :: ri(1:ndir) !r^i
    double precision                :: bi(1:ndir) ! b^i, note that this is a rescaled Bfield, do not mix with the one in prim2con
    double precision                :: r_dot_b, b2, r2
    double precision                :: b_sqr_r_norm_sqr
    double precision                :: vi_hat(1:ndir)   ! v_hat^i
    double precision                :: v_i_hat(1:ndir)   ! v_hat_i
    double precision                :: rescale_factor
    ! variables needed during the root founding
    logical                         :: adjustment
    integer                         :: error_code
    double precision                :: eps_min, eps_max
    double precision                :: W_hat, eps_hat, rho_hat, prs_hat, dummy
    double precision                :: v_hat_sqr, v0_sqr
    double precision                :: chi
    double precision                :: r_bar_sqr
    double precision                :: q_bar
    double precision                :: mu_bounds(2) ! bounds of mu_plus
    double precision                :: mu_plus ! upper bound of mu
    double precision                :: mu ! this is the root of the equation
    ! extra variables needed after the root founding
    double precision                :: B_dot_v
    double precision                :: h_hat
    double precision                :: b2_tmp
    ! parameters should be initialized
    double precision                :: ye_hat 
    double precision                :: temp_hat 
    double precision                :: prs_tmp
    double precision                :: cs2_tmp
    double precision                :: lfac_old, lfac2

    logical                         :: fail_recovery

    cons_tmp = 0.0d0
{#IFDEF DY_SP
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
    end do
}

    flag(ixO^S) = 0

    cons_tmp(ixO^S, d_)    = w(ixO^S,d_)
    cons_tmp(ixO^S, tau_)  = w(ixO^S,tau_)

{^C& cons_tmp(ixO^S, s^C_)  = w(ixO^S,s^C_) \}
! Nothing to do with B field actually as well as Bphi(if glm)
{^C& cons_tmp(ixO^S, b^C_)  = w(ixO^S,b^C_) \}

    if (eos_type == tabulated) cons_tmp(ixO^S, Dye_)  = w(ixO^S, Dye_)

    !Note that v is contravariant
    {do ix^D = ixO^LIM^D \}

      if (eos_type == tabulated) then
        ye_hat = cons_tmp(ix^D, Dye_)/cons_tmp(ix^D, D_)  
        ye_hat = max(eos_yemin, min(eos_yemax,ye_hat))
      else
        ye_hat = 0.0d0
      endif

      ! code test
      !do i = 1, nwflux
      !  if (cons_tmp(ix^D, i) .ne. cons_tmp(ix^D, i) ) then
      !     write(*,*) cons_tmp(ix^D, i), i
      !     write(*,*) 'cons_tmp(ix^D, s1_), cons_tmp(ix^D, s2_), cons_tmp(ix^D, s3_), cons_tmp(ix^D, b1_), &
      !                 cons_tmp(ix^D, b2_), cons_tmp(ix^D, b3_)'
      !     write(*,*) cons_tmp(ix^D, s1_), cons_tmp(ix^D, s2_), cons_tmp(ix^D, s3_), cons_tmp(ix^D, b1_), &
      !                cons_tmp(ix^D, b2_), cons_tmp(ix^D, b3_)
      !     write(*,*) x(ix^D, :)
      !     call mpistop('NaN cons_tmp in C2P')
      !  endif
      !enddo

      adjustment = .False.
      if ( (cons_tmp(ix^D, D_) <= small_D)) then
         ! atmosphere handling 
         ! skip the primitive recovery
         if (old_bhac_safety) then
             call usr_atmo_pt(ix^D,w(ix^D,1:nw),x(ix^D,1:ndim))
            cycle
         endif
         W_hat          = 1.0d0
         vi_hat(1:ndir) = 0.0d0
         rho_hat        = small_rho
         ye_hat         = big_ye       ! depends on input data
         eps_hat        = small_eps
         adjustment     = .True.    ! adjust all conserved variables
         ! fixme: add small_xi
      else 
         ! Step 0: initialize all the variables
         ! rescale the conserved variables
         q = cons_tmp(ix^D, tau_) / cons_tmp(ix^D, D_)

         {^C& r_i(^C) = cons_tmp(ix^D, s^C_) / cons_tmp(ix^D, D_) \}

         {^C& bi(^C) = w(ix^D, b^C_)/dsqrt( cons_tmp(ix^D, D_) ) \}
   
         ! some useful relations
         r_dot_b = 0.0d0
         b2 = 0.0d0 
         r2 = 0.0d0

         do idir = 1, ndir
            r_dot_b = r_dot_b + r_i(idir) * bi(idir)
         end do
    
{#IFNDEF DY_SP    
         call square3u_ixD(ix^D,myM,bi(1:ndir),b2)
         call raise3_ixD(ix^D,myM,r_i(1:ndir),ri(1:ndir))
         call square3u_ixD(ix^D,myM,ri(1:ndir),r2)
}
{#IFDEF DY_SP
         do idir = 1, ndir
            b2 = b2 + bi(idir)**2 * gamma(ix^D,idir,idir)
            r2 = r2 + r_i(idir)**2 / gamma(ix^D,idir,idir)
         enddo
}

         v0_sqr = r2 / ( h0**2 + r2 )
   
         ! decompose the momentum in parts parallel and normal to the magnetic
         ! field, and we need b**2 * r_norm**2 only:
         b_sqr_r_norm_sqr = b2 * r2 - r_dot_b**2
   
         ! Step 1: solve mu+
         ! note that func_for_mu_plus(mu) <= 0 is always true.
         if ( ( dabs(func_for_mu_plus(one_over_h0)) <= Kc2p_tolerance )) then
            mu_plus = one_over_h0
         else
            ! get the bounds for the root
            mu_bounds(1) = 1.0d0 / dsqrt( h0**2 + r2 )
            call get_vars(one_over_h0, r_bar_sqr_out=mu_bounds(2))
            mu_bounds(2) = 1.0d0 / dsqrt( h0**2 + mu_bounds(2) ) 

            ! mathematically, mu_max >= mu_min as long as r^2,b^2, r^l b_l
            ! respect Schwarz inequality. However, it is possible that 
            ! mu_min=mu_max, which is why we have to account for roundoff 
            ! errors
            mu_bounds(1) = mu_bounds(1) * ( 1.0d0 - smalldouble )
            mu_bounds(2) = mu_bounds(2) * ( 1.0d0 + smalldouble )

            ! If this ad-hoc margin was not enough, we just use a much wider 
            ! bracket (the original one given in the article).
            if ( mu_bounds(2) < mu_bounds(1) ) then
               mu_bounds(1) = 0.0d0
               mu_bounds(2) = one_over_h0 * ( 1.0d0 + smalldouble ) 
            end if

            ! try constrained newton raphson first
            call rootfinding_constrained_newton_raphson(mu_plus, mu_bounds(1), mu_bounds(2), &
                   Kc2p_tolerance, 50, error_code, func_for_mu_plus, d_mu_plus)
            !call rootfinding_brent(mu_plus, mu_bounds(1), mu_bounds(2), Kc2p_tolerance, 50, error_code, func_for_mu_plus)

            if ( error_code /= 0 ) then
               !if newton raphson failed, try illinois
               write(*,*) "Warning: Newton raphson failed when finding mu+, now try brent"
               if (any(mu_bounds(:) /= mu_bounds(:))) then
                  mu_bounds(1) = 0.0d0
                  mu_bounds(2) = one_over_h0 * (1.0d0 + smalldouble) 
               end if
               call rootfinding_brent(mu_plus, mu_bounds(1), mu_bounds(2), Kc2p_tolerance, iter_max, error_code, func_for_mu_plus)
            end if
            ! check mu+
            select case (error_code)
            !case (0) ! root is found
            !   write(*,*) "z= ", z
            case (-1) ! nothing happened
               call mpistop("have you ever attemp to find mu+ in con2prim?")
            case (1) ! failed to find the root
               call mpistop("Fail to find the root for mu_plus in con2prim")
               !write(*,*) "Warning: Fail to find the root for mu_plus, now use the maximum value"
               !mu_plus = one_over_h0
            case (2) ! root is NaN
               call mpistop("NaN mu_plus in con2prim")
               !write(*,*) "Warning: mu_plus is NaN, now use the maximum value"
               !mu_plus = one_over_h0
            case (3) ! z is not bracketed
               write(*,*) "cons_tmp(ix^D,D_), cons_tmp(ix^D,tau_),&
                           cons_tmp(ix^D,s1_), cons_tmp(ix^D,s2_), cons_tmp(ix^D,s3_)"
               write(*,*) cons_tmp(ix^D,D_), cons_tmp(ix^D,tau_),  &
                         cons_tmp(ix^D,s1_), cons_tmp(ix^D,s2_), cons_tmp(ix^D,s3_)
               call mpistop("the root is not bracketed for mu_plus in con2prim")
            end select
         end if
   
         ! Step 2: solve mu
         ! check if the bound make sense, use mu_bounds(1) as a tmp variable
         mu_bounds(1) = func_for_mu(mu_plus)
         if ( dabs(mu_bounds(1)) <= Kc2p_tolerance ) then
            mu = mu_plus
         else
            if ( mu_bounds(1) < 0.0d0 ) mu_plus = one_over_h0
            call rootfinding_brent(mu, 0.0d0, mu_plus, Kc2p_tolerance, iter_max, error_code, func_for_mu)
            !call rootfinding_brent(mu, 0.0d0, mu_plus*(1.0d0+smalldouble), Kc2p_tolerance, iter_max, error_code, func_for_mu)
            ! illinois is slighty slower
            !call rootfinding_illinois(mu, 0.0d0, mu_plus*(1.0d0+smalldouble), Kc2p_tolerance, iter_max, error_code, func_for_mu)
            ! check mu
            select case (error_code)
            !case (0) ! root is found
            !   write(*,*) "z= ", z
            case (-1) ! nothing happened
               call mpistop("have you ever attemp to find mu in con2prim?")
            case (1) ! failed to find the root
               !if (usr_W_limit_scheme) then
               !  call fixp_usr_pt(ix^D,w(ix^D,1:nw),x(ix^D,1:ndim))
               !  write(*,*) 'Fail to find the root for mu'
               !else
                 call mpistop("Fail to find the root for mu")
               !endif
            case (2) ! root is NaN
               ! This special flag is used to indicate NaN error. 
               ! Normally it wont be negative
               flag(ix^D) = -1 
            case (3) ! z is not bracketed
               write(*,*) "rho_hat, eps_hat, ye_hat"
               write(*,*) rho_hat, eps_hat, ye_hat
               write(*,*) "cons_tmp(ix^D,D_), cons_tmp(ix^D,tau_), &
                           cons_tmp(ix^D,s1_), cons_tmp(ix^D,s2_), cons_tmp(ix^D,s3_)"
               write(*,*) cons_tmp(ix^D,D_), cons_tmp(ix^D,tau_), &
                          cons_tmp(ix^D,s1_), cons_tmp(ix^D,s2_), cons_tmp(ix^D,s3_)
               call mpistop("the root is not bracketed for mu")
            end select
         end if
   
         ! Step 3: calculate all the primitive variables from the solution mu
         call get_vars(mu, chi_out=chi, W_out=W_hat, v_hat_sqr_out=v_hat_sqr, rho_out=rho_hat, &
                       eps_out=eps_hat, temp_out=temp_hat, ye_out = ye_hat)
         ! calculate velocities
{#IFNDEF DY_SP
         call raise3_ixD(ix^D,myM,r_i(1:ndir),ri(1:ndir))
         do idir = 1, ndir
            vi_hat(idir) = mu * chi * ( ri(idir) + mu * r_dot_b * bi(idir) )
         end do
}
{#IFDEF DY_SP
         do idir = 1, ndir
            vi_hat(idir) = mu * chi * ( r_i(idir) / gamma(ix^D,idir,idir) + mu * r_dot_b * bi(idir) )
         end do
}
!         do idir = 1, ndir
!            vi_hat(idir) = mu * chi * ( r_i(idir) / gamma(ix^D,idir,idir) + mu * r_dot_b * bi(idir) )
!         end do
         ! adjust the results if it is invalid
         if ( rho_hat <= small_rho_thr ) then
               if (old_bhac_safety) then
                  call usr_atmo_pt(ix^D,w(ix^D,1:nw),x(ix^D,1:ndim))
                  cycle
               endif
               ! reset the primitive variables
               W_hat          = 1.0d0
               vi_hat(1:ndir) = 0.0d0
               rho_hat        = small_rho
               eps_hat        = small_eps  
               ye_hat         = big_ye    
               adjustment     = .True. 
         else
            ! limit the velocities
            if ( W_hat > lfac_max ) then
               if (usr_W_limit_scheme) then
               else ! default one
                 ! rescale the velocity such that v = v_max and keeping D constant
                 rescale_factor = v_max / dsqrt( v_hat_sqr )
                 vi_hat(1:ndir) = vi_hat(1:ndir) * rescale_factor
                 !v_hat_sqr = v_max**2
                 W_hat = lfac_max
                 ! although D is kept constant, density is changed a bit
                 rho_hat = cons_tmp(ix^D, D_) / W_hat
               end if

               !! caseI: should I this?
               !! Since the rho_hat or sth else is changed
               !! check if eps fails into the validity range
               !call eos_get_eps_range(rho_hat, eps_min, eps_max, ye=ye_hat)
   
               !if ( eps_hat < eps_min ) then
               !   eps_hat = eps_min
               !else if ( eps_hat > eps_max ) then
               !   eps_hat = eps_max
               !end if

               ! case II: maybe I can just bound with eos_epsmin/max?
               eps_hat = max( min( eos_epsmax, eps_hat), eos_epsmin)
            end if !endif of lfac_max

         end if !endif of small_rho_thr
      end if ! end of con2prim or small_D checker

      ! update the primitive variables,
      ! note that B-field was updated at the begining of the loop
    if (.not. need_aux_only) then
          ! Update prim
          w(ix^D, rho_)   = rho_hat
          if (eos_type == tabulated) then
            w(ix^D, ye_)  = ye_hat  
            call eos_get_temp_one_grid(w(ix^D, rho_), eps_hat, w(ix^D, T_eps_), w(ix^D, ye_))
          else
            w(ix^D, T_eps_) = eps_hat
          endif

          {#IFDEF DY_SP
          {^C& w(ix^D, u^C_) = vi_hat(^C) * W_hat\} 
          }

          {#IFNDEF DY_SP
          if (eos_type == tabulated) then
            call eos_temp_get_all_one_grid(w(ix^D,rho_), w(ix^D,T_eps_), w(ix^D,ye_), dummy, &
                                           prs = prs_tmp)
          else 
            call eos_get_pressure_one_grid(prs_tmp,w(ix^D,rho_),eps_hat)
          endif

          w(ix^D, pp_)   = prs_tmp
          ! Update aux
          w(ix^D, lfac_) = W_hat      
          w(ix^D, xi_)   = (1.0d0 + eps_hat + prs_tmp/rho_hat) * W_hat**2 * rho_hat   !xi = lfac**2*rho*h

          !update vel
          if (useprimitiveRel) then
             if (w(ix^D,xi_) <smallxi) then
               {^C& w(ix^D, u^C_) = zero \}
             else
               {^C& w(ix^D, u^C_) = vi_hat(^C) * W_hat\} 
             endif
          else
             if (w(ix^D,xi_) <smallxi) then 
               {^C& w(ix^D, u^C_) = zero \}
             else
               {^C& w(ix^D, v^C_) = vi_hat(^C) \} 
             endif
          endif
          }
    else  ! if need_aux_only
          call mpistop('never come to getaux c2p')
    endif   ! endif for need_xi_lfac_only
 
    {enddo^D&\}

    ! Only apply it after reconstruction ! If apply this, it means follow the safety as original BHAC
    if (old_bhac_safety) then
      if (tlow>zero) call fixp_usr(ixI^L,ixO^L,w,x)
    endif

    {do ix^D = ixO^LIM^D \}
      call Simple_check_data_correctness_pt(w(ix^D,1:nw), x(ix^D,1:ndim), cons_tmp(ix^D,1:nw), 'after c2p' )
    {enddo^D&\}

    call conserven(ixI^L,ixO^L, w, x, patchfalse)
    contains

   subroutine get_vars(mu, f_mu_plus, df_mu_plus, f_mu, chi_out, r_bar_sqr_out, W_out, v_hat_sqr_out, rho_out, eps_out, temp_out, ye_out, h_out)
      implicit none
      double precision, intent(in)              :: mu
      double precision, intent(out), optional   :: f_mu_plus
      double precision, intent(out), optional   :: df_mu_plus
      double precision, intent(out), optional   :: f_mu
      double precision, intent(out), optional   :: chi_out
      double precision, intent(out), optional   :: W_out, v_hat_sqr_out
      double precision, intent(out), optional   :: r_bar_sqr_out
      double precision, intent(out), optional   :: rho_out
      double precision, intent(out), optional   :: eps_out
      double precision, intent(out), optional   :: temp_out
      double precision, intent(out), optional   :: ye_out
      double precision, intent(out), optional   :: h_out

      double precision                          :: chi, d_chi
      double precision                          :: r_bar_sqr, d_r_bar_sqr
      double precision                          :: q_bar
      double precision                          :: W_hat, eps_hat, a_hat, rho_hat, prs_hat, ye_hat
      double precision                          :: v_hat_sqr
      double precision                          :: mu_hat
      double precision                          :: nu_A, nu_B, nu_hat
   
      ye_hat = 0.0d0

      chi = 1.0d0 / (1.0d0 + mu * b2)
      if (present(chi_out)) chi_out = chi
      r_bar_sqr = r2 * chi**2 + mu * chi * ( 1.0d0 + chi ) * r_dot_b**2

      if (present(r_bar_sqr_out)) then
             r_bar_sqr_out = r_bar_sqr
             return
      end if

      if (present(f_mu_plus)) then
         f_mu_plus = mu * dsqrt( h0**2 + r_bar_sqr ) - 1.0d0
         return
      end if
      if (present(df_mu_plus)) then
         d_chi = - b2 / (1.0d0 + mu * b2)**2
         r_bar_sqr = r2 * chi**2 + mu * chi * ( 1.0d0 + chi ) * r_dot_b**2 
         d_r_bar_sqr = 2.0d0 * r2 * chi * d_chi &
               + ( chi + mu * d_chi + chi**2 + 2.0d0 * mu * chi * d_chi ) * r_dot_b**2 
         df_mu_plus = dsqrt( h0**2 + r_bar_sqr )
         df_mu_plus = df_mu_plus + mu * d_r_bar_sqr / df_mu_plus
         return
      end if

      q_bar = q - 0.5d0 * b2 &
              - 0.5d0 * mu**2 * chi**2 * b_sqr_r_norm_sqr
      v_hat_sqr = min( mu**2 * r_bar_sqr , v0_sqr)
      if (present(v_hat_sqr_out)) v_hat_sqr_out = v_hat_sqr
      W_hat = 1.0d0 / dsqrt( 1.0d0 - v_hat_sqr )
      if (present(W_out)) W_out = W_hat
      rho_hat = cons_tmp(ix^D, D_) / W_hat
      rho_hat = max( min( eos_rhomax, rho_hat ), eos_rhomin )
      if (present(rho_out)) rho_out = rho_hat
      eps_hat = W_hat * ( q_bar - mu * r_bar_sqr ) &
                + v_hat_sqr * W_hat**2 / ( 1.0d0 + W_hat )

      if (eos_type == tabulated) then 
         ye_hat = cons_tmp(ix^D, Dye_)/ cons_tmp(ix^D, D_)
         ye_hat = max( min(eos_yemax, ye_hat), eos_yemin)
      endif

      if (present(ye_out)) ye_out = ye_hat
      call eos_get_eps_range(rho_hat, eps_min, eps_max, ye=ye_hat)

      eps_hat = max( min( eps_max, eps_hat ), eps_min )

      if (present(eps_out)) eps_out = eps_hat

      ! Never input temp to find pressure here
      call eos_get_pressure_one_grid(prs_hat,rho_hat,eps_hat,ye=ye_hat)

      if (present(f_mu)) then
         a_hat = prs_hat / (rho_hat*(1.0d0+eps_hat))
         nu_A = (1.0d0 + a_hat) * (1.0d0+eps_hat) / W_hat ! h/W
         nu_B = (1.0d0 + a_hat) * ( 1.0d0 + q_bar - mu * r_bar_sqr ) 
         nu_hat = max(nu_A,nu_B)
         mu_hat = 1.0d0 / ( nu_hat + r_bar_sqr * mu )
         f_mu = mu - mu_hat
      end if
   end subroutine get_vars          

   double precision function func_for_mu_plus(mu)
      double precision, intent(in)    :: mu
      call get_vars(mu, f_mu_plus=func_for_mu_plus)
   end function func_for_mu_plus

   double precision function d_mu_plus(mu)
      double precision, intent(in)    :: mu
      call get_vars(mu, df_mu_plus=d_mu_plus)
   end function d_mu_plus

   !> master function f(mu) for finding the root mu
   double precision function func_for_mu(mu)
      double precision, intent(in)    :: mu
      call get_vars(mu, f_mu=func_for_mu)
   end function func_for_mu

  end subroutine imhd_con2prim
end module mod_imhd_con2prim
