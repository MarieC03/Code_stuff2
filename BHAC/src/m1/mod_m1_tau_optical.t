module mod_m1_tau_optical
    use mod_m1_eas
    use mod_m1_eas_param
    use mod_m1_metric_interface

   contains
!****************************************************************************
  !> this rouitne calculates the optical depth via path integral, eikonal eq.
  !> this routines also calculates the equilibrium rates for tau
  subroutine m1_add_tau(qdt,qtC,dx^D,igrid,x,wprim,wradimpl,ixI^L)
    !Note: have to make sure that wprim is prim since we need fluid vars
    use mod_m1_internal
    use mod_m1_eas_microphysical_gray
    use mod_m1_constants 
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    use amrvacdef
    }
    {#IFNDEF UNIT_TESTS
      include 'amrvacdef.f' 
    }
    integer, intent(in) :: igrid,ixI^L
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout)    :: wprim(ixI^S,1:nw) ! wprim is not output
    double precision, intent(inout) :: wradimpl(ixI^S,1:nm1rad_eas)
    double precision, intent(in) :: qdt,qtC
    !internal
    integer :: ix^D, ixO^L
    double precision :: dummy = 1.0d0
    double precision :: eas_eq(1:m1_num_eas)        !>eq. rates
    double precision, dimension(1:3) :: fluid_Prim  !> fluid vars
    double precision, dimension(1:m1_numvars_internal) :: wrad
    type(m1_metric_helper) :: metricM1
    ! internal for tau
    integer :: ixs1,ixs2,ixs3, icomp
    double precision :: tau_empirical(ixI^S) !> empirical optical depth
    double precision :: metric_av(1:6)       !> average metric components
    double precision :: delta(1:3)           !> coord difference
    double precision :: ds2_line, ds_line    !> line element
    double precision :: tau_opt_path         !> optical depth
    double precision :: kappa_ave             !> average opacities
    integer :: ni,nj,nk

    ixO^L=ixI^L^LSUBdixB;

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wprim,x,ixI^L,ixO^L)
    }

    !********************************************************************* 
    !> get equilibrium opacities kappa
    {do ix^D=ixOmin^D,ixOmax^D \}

      fluid_Prim(idx_rho) = wprim(ix^D,rho_) 
      fluid_Prim(idx_T) = wprim(ix^D,T_eps_) 
      fluid_Prim(idx_Ye) = wprim(ix^D,Ye_) 

      {^KSP&  !+++++++++++

        wrad(m1_energy_) = dummy
        {^C& wrad(m1_flux^C_) = dummy \}
        call m1_get_eas(wrad,^KSP,ix^D,eas_eq,fluid_Prim)

        !> convert Kappa_eq from cgs (1/cm) to code units
        eas_eq(:) = eas_eq/LENGTHGF

        wradimpl(ix^D,Q_er^KSP_) = eas_eq(Q_ems) 
        wradimpl(ix^D,kappa_a^KSP_) = eas_eq(k_a)
        wradimpl(ix^D,kappa_s^KSP_) = eas_eq(k_s) 
        wradimpl(ix^D,Q_nr^KSP_) = eas_eq(Q_ems_n)
        wradimpl(ix^D,kappa_nr^KSP_) = eas_eq(k_n) 
    
      \} ! end KSP !+++++++
       
   {enddo ^D&\}
   !********************************************************************* 

    !> empirical formula for tau   
   tau_empirical(ixO^S) = exp(DLOG10(10.d0) &
           *( 0.96d0 * ( DLOG10(wprim(ixO^S,rho_)*INVRHOGF ) / DLOG10(10.0d0) - 11.7d0)))

   !> calculate tau
   {do ix^D=ixOmin^D,ixOmax^D \}

     {^KSP&

     {^IFTWOD
       do ni = -1,1
        do nj = -1,1            
            if((ni .eq. 0) .and. (nj .eq. 0)) then
               continue
            end if
  
            !> finding neighbourin cells
            ixs1 = ix1 + ni
            ixs2 = ix2 + nj
  
            do icomp = 1,6
            metric_av(icomp) = 0.5 *(metricM1%gammaij(ix^D,icomp) + metricM1%gammaij(ixs^D,icomp))
            end do
  
            delta(1) = dx1 *ni
            delta(2) = dx2 *nj
            delta(3) = 0.0d0
  
            !> line element 
            ds2_line = metric_av(1) * delta(1) * delta(1) &
               + metric_av(4) * delta(2) * delta(2) &
               + metric_av(6) * delta(3) * delta(3) &
               + 2.d0 * ( metric_av(2) * delta(1) * delta(2) &
                 + metric_av(3) * delta(1) * delta(3) &
                 + metric_av(5) * delta(2) * delta(3) ) 
  
            ds_line = Dsqrt(ds2_line)
  
            !> average kappa
            kappa_ave = 0.5*(wradimpl(ix^D,kappa_a^KSP_) + wradimpl(ixs^D,kappa_a^KSP_)&
                        + wradimpl(ix^D,kappa_s^KSP_) + wradimpl(ixs^D,kappa_s^KSP_))
  
            !> path integral wie sweeping method
            tau_opt_path = kappa_ave * ds_line + tau_empirical(ix^D)
  
        end do
      end do
      wradimpl(ix^D,tau_path^KSP_) = max(0.0d0,tau_opt_path)
     }

     {^IFTHREED
      do ni = -1,1
        do nj = -1,1
           do nk = -1,1
            
            if((ni .eq. 0) .and. (nj .eq. 0) .and. (nk .eq. 0)) then
               continue
            end if

            !> finding neighbourin cells
            ixs1 = ix1 + ni
            if(^ND .ge. 2) then
            ixs2 = ix2 + nj
              if(^ND .eq. 3) then
                ixs3 = ix3 + nk
              end if 
            end if 

            do icomp = 1,6
            metric_av(icomp) = 0.5 *(metricM1%gammaij(ix^D,icomp) + metricM1%gammaij(ixs^D,icomp))
            end do

            delta(1) = dx1 *ni
            delta(2) = dx2 *nj
            delta(3) = dx3 *nk

            !> line element 
            ds2_line = metric_av(1) * delta(1) * delta(1) &
               + metric_av(4) * delta(2) * delta(2) &
               + metric_av(6) * delta(3) * delta(3) &
               + 2.d0 * ( metric_av(2) * delta(1) * delta(2) &
                 + metric_av(3) * delta(1) * delta(3) &
                 + metric_av(5) * delta(2) * delta(3) ) 

            ds_line = Dsqrt(ds2_line)

            !> average kappa
            kappa_ave = 0.5*(wradimpl(ix^D,kappa_a^KSP_) + wradimpl(ixs^D,kappa_a^KSP_)&
                        + wradimpl(ix^D,kappa_s^KSP_) + wradimpl(ixs^D,kappa_s^KSP_))

            !> path integral wie sweeping method
            tau_opt_path = kappa_ave * ds_line + tau_empirical(ix^D)

           end do
        end do
      end do
      wradimpl(ix^D,tau_path^KSP_) = max(0.0d0,tau_opt_path)
      } !end ifdef threed

    \} !end KSP loop

   {enddo ^D&\} 

   !
  !
  end subroutine m1_add_tau

end module mod_m1_tau_optical

