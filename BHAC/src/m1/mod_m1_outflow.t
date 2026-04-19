module mod_m1_outflow
    implicit none
    public
contains

    {#IFDEF M1
    ! calculate the flux for the integrand of the luminosity for healpix
    subroutine calculate_m1_outflow_flux(ixI^L, ixO^L, w, x, Outflux, species)
        use mod_m1_metric_interface
        use mod_m1_eas_param
        use mod_m1_constants
        {#IFDEF UNIT_TESTS
        use mod_m1_tests
        }
        {#IFNDEF UNIT_TESTS
        include "amrvacdef.f"
        }
        integer, intent(in) :: ixI^L, ixO^L
        double precision, intent(in)     :: x(ixI^S, 1:ndim)
        double precision, intent(inout)  :: w(ixI^S, 1:nw) 
        double precision, intent(inout)  :: Outflux(ixI^S, 1:ndim)
        integer, intent(in)     :: species
        ! local variables
        double precision  :: betaD(ixI^S,1:^NC)
        double precision  :: Enu(ixI^S),Fnu(ixI^S,1:^NC)
        double precision  :: betaDown_ixD(1:^NC),betaUp_ixD(1:^NC)
        type(m1_metric_helper)  :: metricM1 
        integer  :: m1_nvars = 1+1+(^NC)
        integer  :: i, ix^D
        integer  :: eidx, fidx^C

        {#IFDEF UNIT_TESTS2
        call fill_metric(metricM1)
        }
        {#IFNDEF UNIT_TESTS2
          call metricM1%fill_metric(w,x,ixI^L,ixO^L)
        }    
        {^D& do ix^D=ixOmin^D,ixOmax^D \}
          {^C& betaUp_ixD(^C) = metricM1%beta(ix^D,^C) \}
          call metricM1%lower_ixD(ix^D,betaUp_ixD, betaDown_ixD)
          {^C& betaD(ix^D,^C)  = betaDown_ixD(^C) \}
        {^D& enddo \}

        !finding indices for given species
        eidx = erad1_ + (species-1) * m1_nvars 
        {^C& fidx^C = frad1^C_ + (species-1) * m1_nvars \} 
   
        Enu(ixO^S) = w(ixO^S,eidx) / metricM1%sqrtg(ixO^S)
        {^C& Fnu(ixO^S,^C) = w(ixO^S,fidx^C) / metricM1%sqrtg(ixO^S) \}  

        Outflux = 0.0d0 
        !> Fout_i = alpha*F_i - beta_i*E
        do i=1,ndim
            Outflux(ixO^S,i) = metricM1%alp(ixO^S) * Fnu(ixO^S,i) &
                - betaD(ixO^S,i) * Enu(ixO^S) 
        end do    

        !TODO: unit conversion:
        !Q: Does outflow alredy have to be in correct units?
        !If yes: convert code-units to erg/cm^3  : /RHOGF/EPSGF
        !If yes: convert code-units to erg/s  : * TIMEGF / (LENGTHGF**3) 
        !If yes: convert code-units to erg/s  : * TIMEGF * MSUN_GRAMM * C2_CGS
        Outflux(ixO^S,:)  = Outflux(ixO^S,:) * TIMEGF * C2_CGS * MSUN_GRAMM

    end subroutine calculate_m1_outflow_flux
    !------------------------------------------------
    }
end module mod_m1_outflow


!blabl
