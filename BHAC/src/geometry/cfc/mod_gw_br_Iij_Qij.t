!> Reference: L. Blanchet+ 1990; G. Faye and G. Schafer 2003; R. Oechslin 2007
!> The integration of Iij is volume integral: \int_{cell} dV (sqrtgamma/sqrtgamma_hat) (xxx)
!> Order of index = 1: I11   2: I12  3: I13  4: I22  5: I23  6: I33  7: \sqrtgamma D  8: eps

module mod_gw_br_Iij_Qij
  use mod_cfc_parameters
  use mod_imhd_intermediate
  use mod_eos
  implicit none
  public

  contains

  ! fixme: CFC metric init will call here if it needs
  !> compute Q^ij[3], but it is not Q^ij_3dot, see references. Remember using primitive to get here

  ! only for Q3_ij at a cell, without integration
  subroutine compute_Q3_ij(ixI^L, ixO^L, w, x, Q3_ij, sigma, gamma, lfac, w_i)
    include 'amrvacdef.f'

    double precision, intent(in)    :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: Q3_ij(ixI^S,1:6)
    double precision, intent(in)    :: gamma(ixI^S,1:3,1:3), lfac(ixI^S), w_i(ixI^S,1:3), sigma(ixI^S)
    integer, intent(in)             :: ixI^L, ixO^L

    double precision                :: Temp(ixI^S), U(ixI^S), U_i(ixI^S, 1:3), htot(ixI^S), &
                                       d_i_htot(ixI^S,1:3), S(ixI^S), d_i_S(ixI^S,1:3), d_i_prs(ixI^S,1:3)
    integer                         :: ix^D, idir, jdir, pdir, kdir
    double precision                :: u4_D(ixI^S, 1:3), b2(ixI^S), prs_st(ixI^S)
    double precision                :: d_U(ixI^S,1:3), d_U_i(ixI^S,1:3,1:3), dU_i(ixI^S,1:3,1:3), d_w_i(ixI^S,1:3,1:3)
    double precision                :: d_d_U(ixI^S,1:3,1:3), dpU_p(ixI^S), d_dpU_p(ixI^S, 1:3)
    double precision                :: d_d_U_i(ixI^S,1:3,1:3,1:3)! new
    double precision                :: rho_st, rho, p_st_tot, dV, wp(1:3), wpd_d_U(1:3), eps_tmp, temp_tmp, dummy
    double precision                :: Q3_ij_local(1:3,1:3)
    double precision                :: Q3_ij_STF(1:3,1:3)

     U(ixI^S) = w(ixI^S, U_br_)
     if (gw_br_include_dU_i) then
        {^C& U_i(ixI^S,^C) = w(ixI^S, U_br^C_) \}
     else
        U_i = 0.0d0
     endif

    ! take also the ghostcell
    call imhd_get_intermediate_variables(ixI^L, ixI^L, w, x, b2=b2(ixI^S), htot=htot(ixI^S))

    {do ix^D = ixI^LIM^D \}
        if (eos_type == tabulated) then
          rho_st = sigma(ix^D)
          !rho_st = w(ix^D, psi_metric_)**6 * w(ix^D, D_)
          rho = w(ix^D, rho_)
          !rho_st = w(ix^D,rho_)
          Temp(ix^D) = w(ix^D, T_eps_)
          eps_tmp = 0.0d0

          ! fixme: HOTFIX bound the rho_max for khi proj.
          rho_st = max( min( eos_rhomax, rho_st), eos_rhomin)

          ! using d_h - T d_S
          !! compose manual script:  entropy per baryon is dimensionless
          !!call eos_get_pressure_one_grid(prs_st, rho, dummy, temp=Temp(ix^D), ye=w(ix^D,ye_))
          !! then the contribution of T d_i_S is super small
          if (gw_br_couple_weakly) then
            ! Blanchet 1990, h is enthaply - rest mass
            htot(ix^D) = htot(ix^D) - 1.0d0
          else
            call eos_temp_get_all_one_grid(rho_st,Temp(ix^D),w(ix^D,ye_),eps_tmp,prs=prs_st(ix^D))
            !Temp(ix^D) = Temp(ix^D) * mev_gf
            !! htot_st:
            htot(ix^D) = eps_tmp + (prs_st(ix^D) + b2(ix^D)) / rho_st 
          endif
        else
          call mpistop('No eos call for entropy when not using tabulated')
        endif
    {enddo^D&\}


    d_U   = 0.0d0
    d_U_i = 0.0d0
    d_w_i = 0.0d0
    ! 1st derivative \partial_p U  and \partial_p U_i
    do idir = 1, ndim
      call partial_d_5pts( U(ixI^S),ixI^L,ixO^L^LADD2,x(ixI^S,1:ndim),idir,d_U(ixI^S,idir) )
      !call partial_d( S(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_i_S(ixI^S,idir) )
      call partial_d( htot(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_i_htot(ixI^S,idir) )
      !call partial_d( prs_st(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_i_prs(ixI^S,idir) )
        if (gw_br_include_dU_i) then
          do jdir = 1, 3
            call partial_d_5pts( U_i(ixI^S,jdir),ixI^L,ixO^L^LADD2,x(ixI^S,1:ndim),idir,d_U_i(ixI^S,jdir,idir) )
            !call partial_d( U_i(ixI^S,jdir),ixI^L,ixO^L^LADD3,x(ixI^S,1:ndim),idir,d_U_i(ixI^S,jdir,idir) )
            !call partial_d( w_i(ixI^S,jdir),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_w_i(ixI^S,jdir,idir) )
          enddo
        endif
    enddo

    ! make d_p U_i to  dp U_i first
    if (gw_br_include_dU_i) then
      if (use_index_contract) then
        do jdir = 1, ndir
           do idir = 1, ndir
             dU_i(ixI^S,jdir,idir) = d_U_i(ixI^S, jdir, idir) / gamma(ixI^S, idir, idir)
           enddo
        enddo
      else
        dU_i(ixI^S,1:3,1:3) = d_U_i(ixI^S,1:3,1:3)
      endif
    endif

    ! Old method for dpUp got very similar to new method

    ! then, find dpU_p   old method
    dpU_p = 0.0d0
    !do idir = 1, 3
    !   dpU_p(ixI^S) = dpU_p(ixI^S) + dU_i(ixI^S,idir,idir) 
    !enddo

    ! new method: 
    if (gw_br_include_dU_i) then
      d_d_U_i = 0.0d0
      do idir = 1, ndim ! direction of the derivative
         do jdir = 1, 3
           do kdir = 1, 3
             call partial_d( dU_i(ixI^S,jdir,kdir),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_d_U_i(ixI^S,jdir,kdir,idir) )
           enddo
         enddo
      enddo
    endif
    
    d_dpU_p = 0.0d0
    d_d_U = 0.0d0
    ! 2nd derivative 
    do idir = 1, ndim
       do jdir = 1, 3
         call partial_d( d_U(ixI^S,jdir),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_d_U(ixI^S,jdir,idir) )
       enddo

       if (gw_br_include_dU_i) then
         !call partial_d( dpU_p(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_dpU_p(ixI^S,idir) )
         do jdir = 1, 3
           d_dpU_p(ixO^S,idir) = d_dpU_p(ixO^S,idir) + d_d_U_i(ixO^S, jdir, jdir, idir)
         enddo
       else
         d_dpU_p = 0.0d0 
       endif
    enddo

    {do ix^D = ixO^LIM^D \}
      ! Excluding atmosphere contribution
      if (w(ix^D,rho_) .gt. 1e-12) then
        !rho_st = w(ix^D, rho_)
        rho_st = sigma(ix^D)
        !rho_st = w(ix^D, psi_metric_)**6 * w(ix^D, D_)

        !Assumed Minkowsi dV, should I?
        dV = mygeo%dvolume(ix^D) ! it includes sqrtgamma_hat already

        !! get prs_st from rho_st
        !prs_st = 0.0d0
        !if (eos_type == tabulated) then
        !  eps_tmp = 0.0d0
        !  temp_tmp = w(ix^D, T_eps_)
        !  call eos_get_pressure_one_grid(prs_st,w(ix^D,rho_),eps_tmp,temp=temp_tmp, ye=w(ix^D,ye_))
        !  !call eos_get_pressure_one_grid(prs_st,rho_st,eps_tmp,temp=temp_tmp, ye=w(ix^D,ye_))
        !else
        !  eps_tmp = w(ix^D, T_eps_)
        !  call eos_get_pressure_one_grid(prs_st,w(ix^D,rho_),eps_tmp)
        !  !call eos_get_pressure_one_grid(prs_st,rho_st,eps_tmp)
        !endif
        !
        !! adding magnetic contribution
        !prs_st = prs_st + b2(ix^D)/2.0d0

        do idir = 1, 3
          if (use_index_contract) then
            wp(idir) = w_i(ix^D, idir) / gamma(ix^D, idir, idir)
          else
            wp(idir) = w_i(ix^D,idir)
          endif
        enddo

        wpd_d_U = 0.0d0
        ! w^p d_j d_p U
        do jdir = 1, 3
          do pdir = 1, 3
            wpd_d_U(jdir) = wpd_d_U(jdir) + wp(pdir) * d_d_U(ix^D, pdir, jdir)
          enddo
        enddo

        do idir = 1, 3
           do jdir = 1, 3
              !Q3_ij_local(idir,jdir) = dV * (4.0d0 * rho_st * w_i(ix^D, idir) * d_i_prs(ix^D,jdir)/rho_st + &
              Q3_ij_local(idir,jdir) = dV * (4.0d0 * rho_st * w_i(ix^D, idir) * (d_i_htot(ix^D,jdir) ) + &
              !Q3_ij_local(idir,jdir) = dV * (4.0d0 * rho_st * w_i(ix^D, idir) * (d_i_htot(ix^D,jdir) - Temp(ix^D) * d_i_S(ix^D,jdir)) + &
              !Q3_ij_local(idir,jdir) = dV * (4.0d0 * prs_st * d_w_i(ix^D,jdir,idir) + &
                                             6.0d0 * rho_st * w_i(ix^D,jdir) * d_U(ix^D,idir) + &
                                             2.0d0 * rho_st * x(ix^D,idir) *&
                                               (wpd_d_U(jdir) - d_dpU_p(ix^D, jdir))  )
           enddo
        enddo
    
        ! fixme: Aij should be same as A_ij for STF?
        call STF(Q3_ij_local(1:3,1:3), Q3_ij_STF(1:3,1:3))
   
        ! NaN checker, pls comment after code test
        !do idir = 1, ndir
        !   do jdir = 1, ndir
        !      if (Q3_ij_STF(idir,jdir) .ne. Q3_ij_STF(idir,jdir)) then
        !         write(*,*) Q3_ij_STF(idir,jdir), idir, jdir
        !         call mpistop('NaN Q3_ij at GW_BR solver')
        !      endif
        !   enddo
        !enddo
   
        Q3_ij(ix^D, 1) = Q3_ij_STF(1,1)
        Q3_ij(ix^D, 2) = Q3_ij_STF(1,2)
        Q3_ij(ix^D, 3) = Q3_ij_STF(1,3)
        Q3_ij(ix^D, 4) = Q3_ij_STF(2,2)
        Q3_ij(ix^D, 5) = Q3_ij_STF(2,3)
        Q3_ij(ix^D, 6) = Q3_ij_STF(3,3)
      else
        Q3_ij(ix^D,:) = 0.0d0
      endif
    {enddo^D&\}

  end subroutine compute_Q3_ij

  ! Eq.~(6.5) in Blanchet 1990
 
  ! only for Q2_ij at a cell, without integration
  subroutine compute_Q2_ij(ixI^L, ixO^L, w, x, Q2_ij, lfac, dV_star)
    include 'amrvacdef.f'

    double precision, intent(in)    :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: Q2_ij(ixI^S,1:6)
    double precision, intent(in)    :: lfac(ixI^S), dV_star(ixI^S)
    integer, intent(in)             :: ixI^L, ixO^L

    integer                         :: ix^D, idir, jdir
    double precision                :: Q2_ij_local(1:3,1:3), Q2_ij_STF(1:3,1:3), vi(1:3), alp_prime
    double precision                :: U(ixI^S), d_U(ixI^S,1:3)


    Q2_ij = 0.0d0

    U(ixI^S) = w(ixI^S, U_br_)
 
    do idir = 1, ndim
      call partial_d( U(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_U(ixI^S,idir) )
    enddo

    {do ix^D = ixO^LIM^D \}
        ! using alp_prime for postprocessing --> lower freq. --> not that match full GR results
        !if (use_h00_coupling) then
        !  alp_prime = dsqrt( w(ix^D, alp_metric_)**2 - w(ix^D,h_00_) ) 
        !else
          alp_prime = w(ix^D, alp_metric_)
        !endif

        {^C& vi(^C) = alp_prime * w(ix^D, u^C_)/lfac(ix^D) - w(ix^D, beta_metric^C_) \}

        do idir = 1, 3
           do jdir = 1, 3
              Q2_ij_local(idir,jdir) = 2.0d0 * dV_star(ix^D) * w(ix^D, D_) * &
                                         ( vi(idir) * vi(jdir) + x(ix^D, idir) * d_U(ix^D, jdir) )
           enddo
        enddo

        call STF(Q2_ij_local(1:3,1:3), Q2_ij_STF(1:3,1:3))

        Q2_ij(ix^D, 1) = Q2_ij_STF(1,1)
        Q2_ij(ix^D, 2) = Q2_ij_STF(1,2)
        Q2_ij(ix^D, 3) = Q2_ij_STF(1,3)
        Q2_ij(ix^D, 4) = Q2_ij_STF(2,2)
        Q2_ij(ix^D, 5) = Q2_ij_STF(2,3)
        Q2_ij(ix^D, 6) = Q2_ij_STF(3,3)
    {enddo^D&\}
  end subroutine compute_Q2_ij

!> old : t-1 slice;  old2 : t-2 slice
!> I3dot(t-1) = [ Idot(t) - 2* Idot(t-1) + Idot(t-2) ] / delta_t^2

  subroutine compute_Iij_dot_grid(ixI^L, ixO^L, s, Iij_dot, w_old, delta_t)
    use mod_full_gr_source

    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L

    double precision, intent(inout) :: Iij_dot(ixI^S,1:6)
    double precision, intent(inout) :: w_old(ixI^S,1:6)
    double precision, intent(in)    :: delta_t
    type(state), intent(in)         :: s
    integer                         :: ix^D, i, j, idir, k
    double precision                :: r, A, B, xkdkU, xkvk, xdotkdkU, dV, dummy, v, v_dot, xkvkdot, xdotvk
    double precision                :: rho, D_st, D_st_old, eps_old, eps, eps_dot, D_st_dot, r_dot 
    double precision                :: vi_old(1:3), vi(1:3), vi_dot(1:3), U_old, U_dot, wi(1:3), wi_old(1:3), wi_dot(1:3), x_dot(1:3)
    double precision                :: Iij_dot_STF(1:3,1:3), Iij_dot_local(1:3,1:3)

    double precision                :: d_U(ixI^S,1:3), U(ixI^S), v2(ixI^S), lfac(ixI^S), u4_D(ixI^S,0:3), w_i(ixI^S,0:3)
    double precision                :: g(ixI^S,0:3,0:3), alp(ixI^S), beta(ixI^S,1:3), gamma(ixI^S,1:3,1:3), htot(ixI^S)

    associate(x=>s%x%x,w=>s%w%w,mygeo=>s%geo)
 
    if (delta_t == 0.0d0) call mpistop('delta_t = 0 in compute_Iij_dot')

     if (gw_br_use_wi) then
      call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, gamma=gamma(ixI^S,1:3,1:3), &
                                           lfac=lfac(ixI^S))
           alp(ixI^S)     = w(ixI^S, alp_metric_)
      {^C& beta(ixI^S,^C) = w(ixI^S, beta_metric^C_) \}
      call get_g_up_munu(g(ixI^S,0:3,0:3), alp(ixI^S), beta(ixI^S,1:3), gamma(ixI^S,1:3,1:3), ixI^L, ixO^L)

      ! Rezzolla Book Eq. 7.22
      {^C& u4_D(ixI^S,^C) = w(ixI^S,u^C_) * gamma(ixI^S,^C,^C) \}

      ! Rezzolla Book Eq. 7.21
      u4_D(ixI^S,0) = lfac(ixI^S) * (-alp(ixI^S) + {^C& gamma(ixI^S,^C,^C) * beta(ixI^S,^C) * w(ixI^S,u^C_)/lfac(ixI^S) +} )

      ! w_i, but it is not purely spatial just like four velocity
      do idir = 0, 3
        w_i(ixI^S,idir) = (1.0d0 + htot(ixI^S) ) * u4_D(ixI^S,idir)
      enddo

      v2(ixI^S) = 0.0d0
      do i = 1, 3
        v2(ixI^S) = v2(ixI^S) + g(ixI^S,i,i) * w_i(ixI^S,i)**2
      enddo

    else
      call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, &
                                           v2=v2(ixI^S), lfac=lfac(ixI^S))
    endif

    U(ixI^S) = w(ixI^S, U_br_)

    do idir = 1, ndim
      call partial_d( U(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_U(ixI^S,idir) )
    enddo
     
    {do ix^D = ixO^LIM^D \}
        ! we take the upper instead of lower indices, since it doesnt care up or down
        rho      = w(ix^D, rho_)
        D_st     = w(ix^D, D_) * w(ix^D, psi_metric_)**6
        D_st_old = w_old(ix^D,1)

        ! eps should be same definition when B-field is presented
        if (eos_type == tabulated) then
          call eos_get_eps_one_grid(dummy,rho,eps,temp=w(ix^D,T_eps_), ye=w(ix^D,ye_))
        else
          eps = w(ix^D, T_eps_)
        endif
        eps_old  = w_old(ix^D,2)

 
        if (gw_br_use_wi) then
          wi(:) = 0.0d0
          do i = 1,3
             do j = 0,3
                wi(i) = wi(i) + g(ix^D, i,j) * w_i(ix^D,j)
             enddo
          enddo
          {^C& vi(^C) = wi(^C) \} 
!code test
!vi = 0.0d0
        else
          {^C& vi(^C)   = w(ix^D, u^C_)/lfac(ix^D) \}
        end if

        {^C& vi_old(^C) = w_old(ix^D,2+^C) \}
        {^C& x_dot(^C)  = w(ix^D, u^C_)/lfac(ix^D) \} ! it must be the actual v^i = dx^i/dt
!code test
! without x_dot, no weird values
!x_dot = 0.0d0
!v2 = 0.0d0
  
        U_old = w_old(ix^D,6)

        ! time deviatives
        D_st_dot   = (D_st - D_st_old)/delta_t
        eps_dot    = (eps - eps_old)/delta_t
        ! higher order terms
    {^C& vi_dot(^C) = (vi(^C) - vi_old(^C))/delta_t \}
        U_dot      = (U(ix^D) - U_old)/delta_t
        v = dsqrt( vi(1)**2 + vi(2)**2 + vi(3)**2 )
        if (v .ne. 0.0d0) then
          v_dot = (vi(1) * vi_dot(1) + vi(2) * vi_dot(2) + vi(3) * vi_dot(3)) / v
        else
          v_dot = 0.0d0
        endif
!code test
!eps_dot = 0.0d0
!d_U(ix^D,:) = 0.0d0



        !Assumed Minkowsi dV, should I?
        dV = mygeo%dvolume(ix^D) ! it includes sqrtgamma_hat already
        r = dsqrt( x(ix^D,1)**2 + x(ix^D,2)**2 + x(ix^D,3)**2 ) 

        xkvk = 0.0d0
        xkvk = {^C& x(ix^D,^C) * x_dot(^C) +}
        r_dot = xkvk/r
   

        xkvk = 0.0d0
        do k = 1, 3
          xkvk = xkvk + x(ix^D,k) * vi(k)
        enddo

        xkdkU = 0.0d0
        do k = 1, 3
          xkdkU = xkdkU + x(ix^D,k) * d_U(ix^D,k)
        enddo

        xdotkdkU = 0.0d0
        do k = 1, 3
          xdotkdkU = xdotkdkU + x_dot(k) * d_U(ix^D,k)
        enddo

        xkvkdot = 0.0d0
        do k = 1, 3
          xkvkdot = xkvkdot + x(ix^D,k) * vi_dot(k)
        enddo

        xdotvk = 0.0d0
        do k = 1, 3
          xdotvk = xdotvk + x_dot(k) * vi(k)
        enddo

        do i = 1, 3
           do j = 1, 3


              A = x(ix^D,i) * x(ix^D,j) * (1.0d0 + v2(ix^D)/2.0d0 - U(ix^D) + eps) +&
                  11.0d0/21.0d0 * r**2 * vi(i) * vi(j) -&
                  4.0d0/7.0d0 * x(ix^D,i) * xkvk * vi(j) +&
                  4.0d0/21.0d0 * v2(ix^D) * x(ix^D,i) * x(ix^D,j) +&
                  11.0d0/21.0d0 * r**2 * x(ix^D,i) * d_U(ix^D,j) -&
                  17.0d0/21.0d0 * x(ix^D,i) * x(ix^D,j) * xkdkU

   

              B = x_dot(i) * x(ix^D,j) * (1.0d0 + v2(ix^D)/2.0d0 - U(ix^D) + eps) +&
                  x(ix^D,i) * x_dot(j) * (1.0d0 + v2(ix^D)/2.0d0 - U(ix^D) + eps) +&
                  x(ix^D, i) * x(ix^D,j) * eps_dot +&  ! 1st
                  11.0d0/21.0d0 * 2.0d0 * r * r_dot * vi(i) * vi(j) -& ! 2nd 
                  4.0d0/7.0d0 * (x_dot(i) * xkvk * vi(j) + x(ix^D,i) * vi(j) * xdotvk) +& !3rd
                  4.0d0/21.0d0 * (v2(ix^D) * x_dot(i) * x(ix^D,j) + v2(ix^D) * x(ix^D,i) * x_dot(j)) +& !4th
                  11.0d0/21.0d0 * (2.0d0 * r * r_dot * x(ix^D,i) * d_U(ix^D,j) + r**2 * x_dot(i) * d_U(ix^D,j)) -& ! 5th
                  17.0d0/21.0d0 * (x_dot(i) * x(ix^D,j) * xkdkU + x(ix^D,i) * x_dot(j) * xkdkU + x(ix^D,i) * x(ix^D,j) * xdotkdkU)

              ! adding higher order terms: vi_dot, U_dot
!code test
! U_dot, has huge contributions during inspiral
!U_dot = 0.0d0
!v_dot, vi_dot have 1e-4 - 1e-5 contributions on h_00
!v_dot = 0.0d0
!vi_dot = 0.0d0

!code test, D_st_dot --> no changes
D_st_dot = 0.0d0

              B = B + x(ix^D,i) * x(ix^D,j) * ( v * v_dot - U_dot) +&
                  11.0d0/21.0d0 * (r**2 * vi_dot(i) * vi(j) + r**2 * vi(i) * vi_dot(j)) -&
                  4.0d0/7.0d0 * (x(ix^D,i) * xkvk * vi_dot(j) + x(ix^D,i) * xkvkdot * vi(j)) +&
                  4.0d0/21.0d0 * (2.0d0 * v * v_dot * x(ix^D,i) * x(ix^D,j))
            
              Iij_dot_local(i,j) = dV * ( D_st_dot * A + D_st * B)

              ! NaN checker
              if (Iij_dot_local(i,j) .ne. Iij_dot_local(i,j) ) then
                 write(*,*) A, B, D_st_dot, vi_dot(:), U_dot, v_dot
                 call mpistop('NaN in Iij')
              endif

           enddo ! enddo of j
        enddo ! enddo of i

        call STF(Iij_dot_local(1:3,1:3), Iij_dot_STF(1:3,1:3))

        Iij_dot(ix^D, 1) = Iij_dot_STF(1,1)
        Iij_dot(ix^D, 2) = Iij_dot_STF(1,2)
        Iij_dot(ix^D, 3) = Iij_dot_STF(1,3)
        Iij_dot(ix^D, 4) = Iij_dot_STF(2,2)
        Iij_dot(ix^D, 5) = Iij_dot_STF(2,3)
        Iij_dot(ix^D, 6) = Iij_dot_STF(3,3) 

!code test no NaN in c2p
!Iij_dot(ix^D,:) = 0.0d0
       
         w_old(ix^D, 1)    = w(ix^D, D_) * w(ix^D, psi_metric_)**6
         w_old(ix^D, 2)    = eps
    {^C& w_old(ix^D, 2+^C) = vi(^C) \}
         w_old(ix^D, 6)    = U(ix^D)
    {enddo^D&\}

    ! at the end, store the newest w for w_old

    end associate
  end subroutine compute_Iij_dot_grid

  subroutine compute_Iij_dot_grid_Oechslin2007(ixI^L, ixO^L, s, Iij_dot, w_old, delta_t, delta_t_old)
    use mod_full_gr_source

    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L

    double precision, intent(inout) :: Iij_dot(ixI^S,1:6)
    double precision, intent(inout) :: w_old(ixI^S,1:6)
    double precision, intent(in)    :: delta_t, delta_t_old
    type(state), intent(in)         :: s
    integer                         :: ix^D, i, j, idir, k
    double precision                :: A, B, xkdkU, xkvk, xdotkdkU, dV, dummy, v, v_dot, xkvkdot, vkdkU
    double precision                :: rho, D_st, D_st_old, eps_old, eps, eps_dot, D_st_dot, r_dot, r, U
    double precision                :: bracket, v2, bracket_dot, bracket_old
    double precision                :: vi_old(1:3), vi(1:3), vi_dot(1:3), U_old, U_dot, wi(1:3), wi_old(1:3), wi_dot(1:3), x_dot(1:3)
    double precision                :: Iij_dot_STF(1:3,1:3), Iij_dot_local(1:3,1:3)

    double precision                :: d_psi(ixI^S,1:3), psi(ixI^S), lfac(ixI^S), u4_D(ixI^S,0:3), w_i(ixI^S,0:3)
    double precision                :: U_st(ixI^S), d_U_st(ixI^S,1:3)
    double precision                :: g(ixI^S,0:3,0:3), alp(ixI^S), beta(ixI^S,1:3), gamma(ixI^S,1:3,1:3)
    double precision, parameter     :: scale_fac = 0.0d0

    associate(x=>s%x%x,w=>s%w%w,mygeo=>s%geo)

    if (delta_t == 0.0d0) call mpistop('delta_t = 0 in compute_Iij_dot')

    call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, gamma(ixI^S,1:3,1:3),&
                                         lfac=lfac(ixI^S))

    psi(ixI^S) = w(ixI^S,psi_metric_)
    U_st(ixI^S)   = w(ixI^S, U_br_)

    do idir = 1, 3
      call partial_d( U_st(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_U_st(ixI^S,idir) )
      !call partial_d( psi(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_psi(ixI^S,idir) )
    enddo

    {do ix^D = ixO^LIM^D \}
        ! old quantities
        D_st_old    = w_old(ix^D,1)
        bracket_old = w_old(ix^D,2)

        D_st     = w(ix^D, D_) * w(ix^D, psi_metric_)**6
        rho      = w(ix^D, rho_)
        ! Newtonian potential, which is > U_br
        !U = 2.0d0 * (psi(ix^D) - 1.0d0)

        U = U_st(ix^D)
        ! eps should be same definition when B-field is presented
        if (eos_type == tabulated) then
          ! using rho* to get eps
          call eos_get_eps_one_grid(dummy,D_st,eps,temp=w(ix^D,T_eps_), ye=w(ix^D,ye_))
          !call eos_get_eps_one_grid(dummy,rho,eps,temp=w(ix^D,T_eps_), ye=w(ix^D,ye_))
        else
          !call eos_get_eps_one_grid(dummy,D_st,eps)
          call mpistop('Fixme for the eps(D_st)')
          eps = w(ix^D, T_eps_)
        endif

        ! vi defintion, case I: vi = vi in bhac+ code
        {^C& vi(^C) = w(ix^D, u^C_)/lfac(ix^D) \}
        ! case II: vi = u^i/u_0 following Lioutas2024
        !{^C& vi(^C) = w(ix^D,alp_metric_)*w(ix^D, u^C_)/lfac(ix^D) - w(ix^D, beta_metric^C_)\}

        ! No gamma_ij
        v2 = 0.0d0
        if (use_index_contract) then
          v2 = {^C& vi(^C) * vi(^C) * gamma(ix^D,^C,^C) +}
        else
          v2 = {^C& vi(^C) * vi(^C) +}
        endif
        !Assumed Minkowsi dV, should I?
        dV = mygeo%dvolume(ix^D) ! it includes sqrtgamma_hat already

        r = dsqrt( x(ix^D,1)**2 + x(ix^D,2)**2 + x(ix^D,3)**2 )

        xkvk = 0.0d0
! no need contract indices here?
        xkvk = {^C& x(ix^D,^C) * vi(^C) +}
        r_dot = xkvk/r

        bracket = 2.0d0 + v2/2.0d0 - U + eps

        if ( abs(delta_t - delta_t_old)/delta_t .le. tol_rel_dt ) then
          ! Time derivatives:
            ! restrict the changes of time d to avoid time shock
            !if  ( abs(D_st - D_st_old)/D_st_old < 0.2d0) then
            !  D_st_dot = (D_st - D_st_old)/delta_t
            !else
            !  D_st_dot = 0.0d0
            !endif
          D_st_dot = (D_st - D_st_old)/delta_t
          bracket_dot = (bracket - bracket_old)/delta_t
        else
          D_st_dot = 0.0d0
        endif

!code test less dotted quantities --> more stable 
!1. both D_st_dot and bracket_dot non zeros --> very big h00 --> BH ; update every
!2. bracket_dot = 0, still very big h00 --> BH ; update every
!3.  ADD ANY time derivative inside time-D --> intabilities
!D_st_dot = 0.0d0  ! range: [-1e-5  to 1e-5, mostly 1e-10]
bracket_dot = 0.0d0

        xkdkU = 0.0d0
        !xkdkU = 2.0d0 * ({^C& x(ix^D,^C)*d_psi(ix^D,^C) +})
        xkdkU = {^C& x(ix^D,^C)*d_U_st(ix^D,^C) +}
        vkdkU = 0.0d0
        !vkdkU = 2.0d0 * ({^C& vi(^C)*d_psi(ix^D,^C) +})
        vkdkU = {^C& vi(^C)*d_U_st(ix^D,^C) +}
        ! No gamma_ij
        xkvk = 0.0d0
        if (use_index_contract) then
          xkvk = {^C& x(ix^D,^C)*vi(^C)*gamma(ix^D,^C,^C) +}
        else
          xkvk = {^C& x(ix^D,^C)*vi(^C) +}
        endif

        do i = 1, 3
           do j = 1, 3
              A = x(ix^D,i) * x(ix^D,j) * bracket + &
                  11.0d0/21.0d0 * r**2 * vi(i) * vi(j) -&
                  4.0d0/7.0d0 * x(ix^D,i) * xkvk * vi(j) +&
                  4.0d0/21.0d0 * v2 * x(ix^D,i) * x(ix^D,j) +&
                  11.0d0/21.0d0 * r**2 * x(ix^D,i) * d_U_st(ix^D,j)  -&
                  17.0d0/21.0d0 * x(ix^D,i) * x(ix^D,j) * xkdkU

              B =  (bracket - 17.0d0/21.0d0 * xkdkU + 4.0d0/21.0d0*v2)*(x(ix^D,i)*vi(j) + vi(i)*x(ix^D,j)) &
              !B =  (bracket - 17.0d0/21.0d0 * xkdkU )*(x(ix^D,i)*vi(j) + vi(i)*x(ix^D,j)) &
                   - 17.0d0/21.0d0 * x(ix^D,i) * x(ix^D,j) * vkdkU &
                   + x(ix^D,i) * x(ix^D,j) * bracket_dot &
                   + 11.0d0/21.0d0 * d_U_st(ix^D,j) * (2.0d0*xkvk*x(ix^D,i) + r**2 *vi(i)) &
                   + 11.0d0/21.0d0 * 2.0d0 * r * r_dot * vi(i) * vi(j) &
                   - 4.0d0/7.0d0 * ( vi(i) * xkvk * vi(j) + x(ix^D,i) * v2 * vi(j))

              Iij_dot_local(i,j) = dV * ( D_st_dot * scale_fac * A + D_st * B)

              !if (Iij_dot_local(i,j) .ne. Iij_dot_local(i,j) ) then
              !   write(*,*) vi(:), psi(ix^D), lfac(ix^D)
              !   call mpistop('NaN in Iij')
              !endif
           enddo
        enddo


        call STF(Iij_dot_local(1:3,1:3), Iij_dot_STF(1:3,1:3))

        Iij_dot(ix^D, 1) = Iij_dot_STF(1,1)
        Iij_dot(ix^D, 2) = Iij_dot_STF(1,2)
        Iij_dot(ix^D, 3) = Iij_dot_STF(1,3)
        Iij_dot(ix^D, 4) = Iij_dot_STF(2,2)
        Iij_dot(ix^D, 5) = Iij_dot_STF(2,3)
        Iij_dot(ix^D, 6) = Iij_dot_STF(3,3)

        w_old(ix^D, 1)    = D_st
        w_old(ix^D, 2)    = bracket
  !  {^C& w_old(ix^D, 2+^C) = vi(^C) \}
  !       w_old(ix^D, 6)    = U
    {enddo^D&\}

    end associate
  end subroutine compute_Iij_dot_grid_Oechslin2007

  subroutine compute_Iij_grid(ixI^L, ixO^L, s, Iij_grid)
    use mod_full_gr_source

    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L

    double precision, intent(inout) :: Iij_grid(ixI^S,1:6)
    type(state), intent(in)         :: s
    integer                         :: ix^D, i, j, idir
    double precision                :: d_psi(ixI^S,1:3), psi(ixI^S), lfac(ixI^S), u4_D(ixI^S,0:3), w_i(ixI^S,0:3)
    double precision                :: gamma(ixI^S,1:3,1:3)
    double precision                :: Iij_local(1:3,1:3), Iij_STF(1:3,1:3)
    double precision                :: A, D_st, xkdkU, xkvk, v2, r_dot, bracket, vi(1:3), r, U, dV
    double precision                :: rho, eps, dummy

    associate(x=>s%x%x,w=>s%w%w,mygeo=>s%geo)

    call imhd_get_intermediate_variables(ixI^L, ixO^L, w, x, gamma(ixI^S,1:3,1:3),&
                                         lfac=lfac(ixI^S))

    psi(ixI^S) = w(ixI^S,psi_metric_)

    do idir = 1, 3
      call partial_d( psi(ixI^S),ixI^L,ixO^L,x(ixI^S,1:ndim),idir,d_psi(ixI^S,idir) )
    enddo

    {do ix^D = ixO^LIM^D \}
        D_st     = w(ix^D, D_) * w(ix^D, psi_metric_)**6
        ! Newtonian potential, which is > U_br
        U = 2.0d0 * (psi(ix^D) - 1.0d0)
        rho      = w(ix^D, rho_)
        ! eps should be same definition when B-field is presented
        if (eos_type == tabulated) then
          call eos_get_eps_one_grid(dummy,rho,eps,temp=w(ix^D,T_eps_), ye=w(ix^D,ye_))
        else
          eps = w(ix^D, T_eps_)
        endif

        ! vi defintion, case I: vi = vi in bhac+ code
        !{^C& vi(^C) = w(ix^D, u^C_)/lfac(ix^D) \}
        ! case II: vi = u^i/u_0 following Lioutas2024
        {^C& vi(^C) = w(ix^D,alp_metric_)*w(ix^D, u^C_)/lfac(ix^D) - w(ix^D, beta_metric^C_)\}

        ! No gamma_ij
        v2 = 0.0d0
        if (use_index_contract) then
          v2 = {^C& vi(^C) * vi(^C) * gamma(ix^D,^C,^C) +}
        else
          v2 = {^C& vi(^C) * vi(^C) +}
        endif
        !Assumed Minkowsi dV, should I?
        dV = mygeo%dvolume(ix^D) ! it includes sqrtgamma_hat already

        r = dsqrt( x(ix^D,1)**2 + x(ix^D,2)**2 + x(ix^D,3)**2 )

        xkvk = 0.0d0
! no need contract indices here?
        xkvk = {^C& x(ix^D,^C) * vi(^C) +}
        r_dot = xkvk/r

        bracket = 1.0d0 + v2/2.0d0 - U + eps

        xkdkU = 0.0d0
        xkdkU = 2.0d0 * ({^C& x(ix^D,^C)*d_psi(ix^D,^C) +})
        ! No gamma_ij
        xkvk = 0.0d0
        if (use_index_contract) then
          xkvk = {^C& x(ix^D,^C)*vi(^C)*gamma(ix^D,^C,^C) +}
        else
          xkvk = {^C& x(ix^D,^C)*vi(^C) +}
        endif

        do i = 1, 3
           do j = 1, 3
              A = x(ix^D,i) * x(ix^D,j) * bracket + &
                  11.0d0/21.0d0 * r**2 * vi(i) * vi(j) -&
                  4.0d0/7.0d0 * x(ix^D,i) * xkvk * vi(j) +&
                  4.0d0/21.0d0 * v2 * x(ix^D,i) * x(ix^D,j) +&
                  11.0d0/21.0d0 * r**2 * x(ix^D,i) * 2.0d0 * d_psi(ix^D,j) -&
                  17.0d0/21.0d0 * x(ix^D,i) * x(ix^D,j) * xkdkU
              Iij_local(i,j) = dV * D_st * A
           enddo
        enddo

        call STF(Iij_local(1:3,1:3), Iij_STF(1:3,1:3))

        Iij_grid(ix^D, 1) = Iij_STF(1,1)
        Iij_grid(ix^D, 2) = Iij_STF(1,2)
        Iij_grid(ix^D, 3) = Iij_STF(1,3)
        Iij_grid(ix^D, 4) = Iij_STF(2,2)
        Iij_grid(ix^D, 5) = Iij_STF(2,3)
        Iij_grid(ix^D, 6) = Iij_STF(3,3)
    {enddo^D&\}
    end associate
  end subroutine compute_Iij_grid

  ! Taking the Symmetric trace-free; same for covariant and contravariant
  subroutine STF(Aij, STF_Aij)
     double precision, intent(in) :: Aij(1:3, 1:3)
     double precision, intent(out):: STF_Aij(1:3, 1:3)
     integer :: i, j

     STF_Aij = 0.0d0
     do i = 1, 3
       do j = 1, 3
         if (i == j) then
           STF_Aij(i,j) = 0.5d0 * Aij(i,j) + 0.5d0 * Aij(j,i) - 1.0d0/3.0d0 * Aij(i,j)
         else
           STF_Aij(i,j) = 0.5d0 * Aij(i,j) + 0.5d0 * Aij(j,i) 
         endif
       enddo
     enddo
  end subroutine STF

end module mod_gw_br_Iij_Qij


