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

  !#############################################################################
  ! Covariant imhd module, based on the public srmhd version of MPI-AMRVAC
  ! 20.10.2015, Oliver Porth
  
  !=============================================================================
  subroutine getdt(w,ixI^L,ixO^L,dtnew,dx^D,x)

    include 'amrvacdef.f'

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
    !-----------------------------------------------------------------------------

    dtnew=bigdouble

  end subroutine getdt
  !=============================================================================
  subroutine checkglobaldata
    use mod_eos
    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------
    if (eos_type == tabulated) then
      ! nothing
      smallxi = -1.0d90
    else if (eos_type == idealgas) then
      minrho = max(zero,smallrho)
      govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
      minp  = max(zero,smallp)
      smallxi=minrho+minp*govergminone
      smalltau = minp/(eqpar(gamma_)-one)
    else
      minrho = max(zero,smallrho)
      {#IFDEF ISO
      minp=eqpar(adiab_)*minrho**eqpar(gamma_)
      govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
      }
      {#IFDEF SYNGE
      !call smallvaluesEOS
      }
    endif 
  end subroutine checkglobaldata
  !=============================================================================
  subroutine initglobaldata

    ! place to initialize some globals
    !use mod_con2prim, only: setup_con2prim
    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------
    eqpar(adiab_)=1.0d0

    {#IFNDEF SYNGE
    eqpar(gamma_)=4.0d0/3.0d0
    }
    {#IFDEF SYNGE
    eqpar(gamma_)=5.d0/3.d0
    }

    {#IFDEF GLM
    eqpar(Cr_)= 0.18d0
    }

    eqpar(a_) = 0.0d0
    eqpar(m_) = 1.0d0

    {#IFDEF ELECTRONS
    eqpar(gammae_)=4.d0/3.d0
    }
    
    limitvalue=smalldouble**2

    !call setup_con2prim

  end subroutine initglobaldata
  !=============================================================================
  subroutine conserve(ixI^L,ixO^L,w,x,patchw)

    ! Transform primitive variables into conservative ones
    ! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
    ! v is contravariant, S is covarariant  
    ! call to smallvalues
    ! --> latter only used for correcting procedure in correctaux
    ! --> input array patchw for spatially selective transformation

    include 'amrvacdef.f'

    integer, intent(in)               :: ixI^L,ixO^L
    double precision, intent(inout)   :: w(ixI^S,nw)
    logical, intent(in)               :: patchw(ixG^T)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    !-----------------------------------------------------------------------------

    call conserven(ixI^L,ixO^L,w,x,patchw)

  end subroutine conserve

  !=============================================================================
{#IFDEF DY_SP
  subroutine conserven(ixI^L,ixO^L,w,x,patchw)
    ! Use RelPrimitive => use  4-veloc as primitive, otherwise uses v^i
    ! Transform primitive variables into conservative ones
    ! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
    ! no call to smallvalues
    ! --> latter only used for correcting procedure in correctaux
    ! --> input array patchw for spatially selective transformation
    use mod_imhd_intermediate
    use mod_imhd_con2prim
    use mod_metric
    use mod_small_values
    use mod_eos
    include 'amrvacdef.f'

    integer, intent(in)                       :: ixI^L,ixO^L
    double precision, intent(inout)           :: w(ixI^S,1:nw)
    double precision, intent(in)              :: x(ixI^S,1:ndim)
    logical, intent(in)                       :: patchw(ixG^T)
    double precision, dimension(ixI^S)        :: sqrV,sqrU,sqrB,VdotB,rhoh,prs_tmp,eps_tmp,lfac,xi
    ! 3-metric gamma_ij
    double precision                          :: gamma(ixI^S,1:3,1:3), temp_local, rho_local
    double precision, dimension(ixI^S,1:^NC)  :: vD, bD, bD_ixD, d, u 
    integer                                   :: inonzero, i, j, ix^D, idir

    !-----------------------------------------------------------------------------
 
    !if (old_bhac_safety) then
    !  if (tlow>zero) call fixp_usr(ixI^L,ixO^L,w,x)
    !else
    !  call imhd_handle_small_values(w(ixI^S,1:nw),x(ixI^S,1:ndim),ixI^L,ixO^L,.true.)
    !endif

    d(ixO^S,:) = 0.0d0
    prs_tmp(ixO^S) = 0.0d0

    call get_gamma_ij_hat(x(ixI^S,1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
    end do
	
      ! diagonal metric
      {^C&  vD(ixO^S, ^C)  = w(ixO^S, u^C_) * gamma(ixO^S,^C,^C)  \}
      {^C&  bD(ixO^S, ^C)  = w(ixO^S, b^C_) * gamma(ixO^S,^C,^C)  \}

    {do ix^D = ixO^LIM^D \}
        if (eos_type == tabulated) then
          temp_local = w(ix^D,T_eps_)
          rho_local  = w(ix^D,rho_)
          call eos_temp_get_all_one_grid(rho_local,temp_local,w(ix^D,ye_),&
                                         eps_tmp(ix^D),prs=prs_tmp(ix^D))
        else
          eps_tmp(ix^D) = w(ix^D, T_eps_)
          call eos_get_pressure_one_grid(prs_tmp(ix^D),w(ix^D,rho_),eps_tmp(ix^D))
        endif
    {enddo^D&\}

    if(useprimitiveRel)then
       ! assumes four velocity computed in primitive (rho u p B) with u=lfac*v
       where(.not.patchw(ixO^S))
          sqrU(ixO^S) ={^C& w(ixO^S,u^C_)*vD(ixO^S,^C)+} 
          sqrV(ixO^S) =sqrU(ixO^S)/(one+sqrU(ixO^S))
          sqrB(ixO^S) ={^C& w(ixO^S,b^C_)*bD(ixO^S,^C)+}
       endwhere
    else
       ! assumes velocity in primitive (rho v p B) 
       where(.not.patchw(ixO^S))
          sqrV(ixO^S) ={^C& w(ixO^S,v^C_)*vD(ixO^S,^C)+}
          sqrB(ixO^S) ={^C& w(ixO^S,b^C_)*bD(ixO^S,^C)+} 
       endwhere
    endif

    ! fill the auxilary variable lfac (lorentz factor)
    if(useprimitiveRel)then
       where(.not.patchw(ixO^S))
          lfac(ixO^S)=dsqrt(one+sqrU(ixO^S)) 
          VdotB(ixO^S)  =({^C& bD(ixO^S,^C)*w(ixO^S,u^C_)+})/lfac(ixO^S)
       endwhere
    else
       call mpistop('DY_SP needs useprimitiveRel')
    endif

    rhoh(ixO^S) = (1.0d0 + eps_tmp(ixO^S) + prs_tmp(ixO^S)/w(ixO^S,rho_)) * w(ixO^S,rho_)

    where(.not.patchw(ixO^S))
       xi(ixO^S) = lfac(ixO^S)**2 * rhoh(ixO^S)
       w(ixO^S,d_) =w(ixO^S,rho_)*lfac(ixO^S)
       ! recycle sqrU array for storing temporary positive 
       ! array for use in energy variable
       sqrU(ixO^S)=sqrB(ixO^S)*sqrV(ixO^S)-VdotB(ixO^S)**2
    endwhere
    where(.not.patchw(ixO^S).and.sqrU(ixO^S)<zero)
       sqrU(ixO^S)=zero
    endwhere

    {#IFDEF ENTROPY
    ! Get back the conserved entropy Gamma rho s = Gamma p/rho**(gamma)
    where(.not.patchw(ixO^S))
       w(ixO^S,Ds_) = w(ixO^S,d_)*w(ixO^S,s_)
    end where
    }

    {#IFDEF ELECTRONS
    ! Get back the conserved electron entropy Gamma rho s_e = Gamma p_e/rho**(gamma)
    where(.not.patchw(ixO^S))
       w(ixO^S,Dse_) = w(ixO^S,d_)*w(ixO^S,se_)
    end where
    }

    {#IFDEF TRACER
    ! We got D, now we can get the conserved tracers:
    where(.not.patchw(ixO^S))
       {^FL&w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)\}
    endwhere
    }

    {#IFDEF EPSINF
    where(.not.patchw(ixO^S))
       w(ixO^S,Drho1_) =w(ixO^S,rho1_)*lfac(ixO^S)
       w(ixO^S,Drho0_) =w(ixO^S,Drho1_)*w(ixO^S,rho0_)
       w(ixO^S,Dn_) =w(ixO^S,n_)*lfac(ixO^S)
       w(ixO^S,Dn0_) =w(ixO^S,Dn_)*w(ixO^S,n0_)
       w(ixO^S,Depsinf_) = w(ixO^S,epsinf_)*w(ixO^S,Drho1_)**(2.0d0/3.0d0) &
            *lfac(ixO^S)**(1.0d0/3.0d0)
    endwhere
    }
    ! fill the vector S
    ! s= (xi + B^2) * v - (v.B) * B
    if(useprimitiveRel)then
       where(.not.patchw(ixO^S))
          {^C& w(ixO^S,s^C_)=(xi(ixO^S)+sqrB(ixO^S))&
               *vD(ixO^S,^C)/lfac(ixO^S) - &
               VdotB(ixO^S)*bD(ixO^S,^C);}
       endwhere
    endif


    !{#IFDEF ENERGY
    ! E = xi - p +B^2/2 + (v^2 B^2 - (v.B)^2)/2 
    ! instead of E use tau= E - D
    where(.not.patchw(ixO^S))
    w(ixO^S,tau_)=xi(ixO^S) - prs_tmp(ixO^S) - w(ixO^S,d_) +& 
            half*(sqrB(ixO^S) + sqrU(ixO^S))
    endwhere
    !}

    if (eos_type == tabulated) then
      where(.not.patchw(ixO^S))
         w(ixO^S,Dye_) =  w(ixO^S, d_) * w(ixO^S, ye_)
      endwhere
    endif

  end subroutine conserven
}
! =========================================================================
{#IFNDEF DY_SP
  subroutine conserven(ixI^L,ixO^L,w,x,patchw)
    ! Use RelPrimitive => use  4-veloc as primitive, otherwise uses v^i
    ! Transform primitive variables into conservative ones
    ! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
    ! no call to smallvalues
    ! --> latter only used for correcting procedure in correctaux
    ! --> input array patchw for spatially selective transformation
    use mod_eos
    use mod_imhd_intermediate
    use mod_imhd_con2prim
    use mod_metric
    use mod_small_values

    include 'amrvacdef.f'

    integer, intent(in)                       :: ixI^L,ixO^L
    double precision, intent(inout)           :: w(ixI^S,1:nw)
    double precision, intent(in)              :: x(ixI^S,1:ndim)
    logical, intent(in)                       :: patchw(ixG^T)
    integer :: ix^D
    double precision, dimension(ixI^S)        :: sqrV,sqrU,sqrB,VdotB,rhoh,prs_tmp,eps_tmp
    ! 3-metric gamma_ij
    double precision                          :: gamma(ixI^S,1:3,1:3)
    double precision, dimension(ixI^S,1:^NC)  :: vD, bD, bD_ixD, d, u 
    integer                                   :: inonzero, i, j

    !-----------------------------------------------------------------------------
      ! code test
      !{do ix^D=ixOmin^D, ixOmax^D \}
      !  call Simple_check_data_correctness_pt(w(ix^D,1:nw), x(ix^D,1:ndim), w(ix^D,1:nw), 'b4 p2c' )
      !{enddo^D&\}
     
    !if (old_bhac_safety) then
    !  if (tlow>zero) call fixp_usr(ixI^L,ixO^L,w,x)
    !else
    !  call imhd_handle_small_values(w(ixI^S,1:nw),x(ixI^S,1:ndim),ixI^L,ixO^L,.true.)
    !endif

    d(ixO^S,:) = 0.0d0

    call lower3(ixI^L,ixO^L,myM,w(ixI^S,u1_:u^NC_),vD)
    call lower3(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),bD)

    {do ix^D = ixO^LIM^D \}
        prs_tmp(ix^D) = w(ix^D, pp_)
        if (eos_type == tabulated) then
          call eos_get_eps_one_grid(prs_tmp(ix^D), w(ix^D,rho_), eps_tmp(ix^D), w(ix^D,T_eps_),&
                                    ye=w(ix^D,ye_))
        else
          eps_tmp(ix^D) = w(ix^D, T_eps_)
        endif
    {enddo^D&\}

    if(useprimitiveRel)then
       ! assumes four velocity computed in primitive (rho u p B) with u=lfac*v
       where(.not.patchw(ixO^S))
          sqrU(ixO^S) ={^C& w(ixO^S,u^C_)*vD(ixO^S,^C)+} 
          sqrV(ixO^S) =sqrU(ixO^S)/(one+sqrU(ixO^S))
          sqrB(ixO^S) ={^C& w(ixO^S,b^C_)*bD(ixO^S,^C)+}
       endwhere
    else
       ! assumes velocity in primitive (rho v p B) 
       where(.not.patchw(ixO^S))
          sqrV(ixO^S) ={^C& w(ixO^S,v^C_)*vD(ixO^S,^C)+}
          sqrB(ixO^S) ={^C& w(ixO^S,b^C_)*bD(ixO^S,^C)+} 
       endwhere
    endif

    ! fill the auxilary variable lfac (lorentz factor)
    if(useprimitiveRel)then
       where(.not.patchw(ixO^S))
          w(ixO^S,lfac_)=dsqrt(one+sqrU(ixO^S)) 
          VdotB(ixO^S)  =({^C& bD(ixO^S,^C)*w(ixO^S,u^C_)+})/w(ixO^S,lfac_)
       endwhere
    else
       where(.not.patchw(ixO^S))
          w(ixO^S,lfac_)=one/dsqrt(one-sqrV(ixO^S))
          VdotB(ixO^S)  ={^C& bD(ixO^S,^C)*w(ixO^S,v^C_)+}
       endwhere
    endif

    rhoh(ixO^S) = (1.0d0 + eps_tmp(ixO^S) + prs_tmp(ixO^S)/w(ixO^S,rho_)) * w(ixO^S,rho_)

    where(.not.patchw(ixO^S))
       w(ixO^S,xi_)=w(ixO^S,lfac_)*w(ixO^S,lfac_)* rhoh(ixO^S)
       w(ixO^S,d_) =w(ixO^S,rho_)*w(ixO^S,lfac_)
       ! recycle sqrU array for storing temporary positive 
       ! array for use in energy variable
       sqrU(ixO^S)=sqrB(ixO^S)*sqrV(ixO^S)-VdotB(ixO^S)**2
    endwhere
    where(.not.patchw(ixO^S).and.sqrU(ixO^S)<zero)
       sqrU(ixO^S)=zero
    endwhere

    {#IFDEF ENTROPY
    ! Get back the conserved entropy Gamma rho s = Gamma p/rho**(gamma)
    where(.not.patchw(ixO^S))
       w(ixO^S,Ds_) = w(ixO^S,d_)*w(ixO^S,s_)
    end where
    }

    {#IFDEF ELECTRONS
    ! Get back the conserved electron entropy Gamma rho s_e = Gamma p_e/rho**(gamma)
    where(.not.patchw(ixO^S))
       w(ixO^S,Dse_) = w(ixO^S,d_)*w(ixO^S,se_)
    end where
    }

    {#IFDEF TRACER
    ! We got D, now we can get the conserved tracers:
    where(.not.patchw(ixO^S))
       {^FL&w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)\}
    endwhere
    }

    {#IFDEF EPSINF
    where(.not.patchw(ixO^S))
       w(ixO^S,Drho1_) =w(ixO^S,rho1_)*w(ixO^S,lfac_)
       w(ixO^S,Drho0_) =w(ixO^S,Drho1_)*w(ixO^S,rho0_)
       w(ixO^S,Dn_) =w(ixO^S,n_)*w(ixO^S,lfac_)
       w(ixO^S,Dn0_) =w(ixO^S,Dn_)*w(ixO^S,n0_)
       w(ixO^S,Depsinf_) = w(ixO^S,epsinf_)*w(ixO^S,Drho1_)**(2.0d0/3.0d0) &
            *w(ixO^S,lfac_)**(1.0d0/3.0d0)
    endwhere
    }
    ! fill the vector S
    ! s= (xi + B^2) * v - (v.B) * B
    if(useprimitiveRel)then
       where(.not.patchw(ixO^S))
          {^C& w(ixO^S,s^C_)=(w(ixO^S,xi_)+sqrB(ixO^S))&
               *vD(ixO^S,^C)/w(ixO^S,lfac_) - &
               VdotB(ixO^S)*bD(ixO^S,^C);}
       endwhere
    else
       where(.not.patchw(ixO^S))
          {^C& w(ixO^S,s^C_)=(w(ixO^S,xi_)+sqrB(ixO^S))*vD(ixO^S,^C) - &
               VdotB(ixO^S)*bD(ixO^S,^C);}
       endwhere
    endif


    !{#IFDEF ENERGY
    ! E = xi - p +B^2/2 + (v^2 B^2 - (v.B)^2)/2 
    ! instead of E use tau= E - D
    where(.not.patchw(ixO^S))
    w(ixO^S,tau_)=w(ixO^S,xi_) - prs_tmp(ixO^S) - w(ixO^S,d_) +& 
            half*(sqrB(ixO^S) + sqrU(ixO^S))
    endwhere
    !}

    if (eos_type == tabulated) then
      where(.not.patchw(ixO^S))
         w(ixO^S,Dye_) =  w(ixO^S, d_) * w(ixO^S, ye_)
      endwhere
    endif

  end subroutine conserven
}
  !=============================================================================
  subroutine primitive(ixI^L,ixO^L,w,x)

    ! Transform conservative variables into primitive ones
    ! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)
    use mod_imhd_con2prim
    use mod_small_values
    include 'amrvacdef.f'

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(inout)    :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)

    logical, dimension(ixG^T)          :: patchw

    !-----------------------------------------------------------------------------
    patchw(ixO^S) = .false.

    call imhd_con2prim(ixI^L, ixO^L, x, w, patchw, .false.)

  end subroutine primitive

  !==============================================================================
  subroutine getv(wprim,x,ixI^L,ixO^L,idims,v)

    ! Returns the contravariant (velocity*lapse - shift)
    ! Now takes input in primitive form for performance!
    use mod_imhd_intermediate
    use mod_cfc_parameters
    include 'amrvacdef.f'

    integer, intent(in)                              :: ixI^L,ixO^L,idims
    double precision, intent(in)                     :: wprim(ixI^S,1:nw)
    double precision, intent(in)                     :: x(ixI^S,1:ndim)
    double precision, dimension(ixI^S), intent(out)  :: v
    double precision                                 :: lfac(ixI^S), alp_prime(ixI^S), gamma_tmp(ixI^S,1:3,1:3), W2v2(ixI^S)
    !-----------------------------------------------------------------------------

{#IFNDEF DY_SP
   v(ixO^S) = myM%alpha(ixO^S)*wprim(ixO^S,u0_+idims)/wprim(ixO^S,lfac_) &
        - myM%beta(idims)%elem(ixO^S)
}
{#IFDEF DY_SP
    call dysp_get_lfac(ixI^L, ixO^L, wprim, x, lfac)
 
    if (use_gw_br .and. use_h00_coupling) then
       alp_prime(ixO^S) = dsqrt( wprim(ixO^S,alp_metric_)**2 - wprim(ixO^S,h_00_) )
       !alp_prime(ixO^S) = wprim(ixO^S,alp_metric_) - wprim(ixO^S,h_00_)/2.0d0/wprim(ixO^S,alp_metric_)
    else
       alp_prime(ixO^S) = wprim(ixO^S,alp_metric_) 
    endif

    v(ixO^S) = alp_prime(ixO^S)*wprim(ixO^S,u0_+idims)/lfac(ixO^S) &
         - wprim(ixO^S,alp_metric_+idims)
}
  end subroutine getv
  !=============================================================================
  subroutine get_cbounds(wL, wR, x, ixI^L, ixO^L, idims, cmax, cmin, method, qtC, metricM1_L, metricM1_R)
   use mod_m1_metric_interface
   include 'amrvacdef.f'
	
   integer, intent(in)             :: ixI^L, ixO^L, idims
   double precision, intent(in)    :: qtC
   double precision, intent(in)    :: wL(ixI^S, 1:nw), wR(ixI^S, 1:nw)
   double precision, intent(in)    :: x(ixI^S, 1:ndim)
   double precision, intent(inout) :: cmax(ixI^S, 1:ncons)
   double precision, intent(inout) :: cmin(ixI^S, 1:ncons)
   character(len=*), intent(in)    :: method
   type(m1_metric_helper), intent(in)   :: metricM1_L, metricM1_R
   integer :: ix^D,ii

   double precision, dimension(ixI^S, 1:ncons,1:2) :: tmp_c
   double precision, dimension(ixI^S, 1:ncons,1:2) :: lambdaL, lambdaR

   {#IFDEF UNITS_TESTS
   double precision :: TESTqtC = 1.0d+10
   }

   call get_lambda(ixI^L, ixO^L, idims, wL(ixI^S, 1:nw), x(ixI^S, 1:ndim), lambdaL(ixI^S,1:ncons,1:2), qtC, metricM1_L) !KEN
   call get_lambda(ixI^L, ixO^L, idims, wR(ixI^S, 1:nw), x(ixI^S, 1:ndim), lambdaR(ixI^S,1:ncons,1:2), qtC, metricM1_R) !KEN

   tmp_c(ixO^S,1:ncons,1)=max(0.0d0, lambdaL(ixO^S,1:ncons,1), lambdaR(ixO^S,1:ncons,1) )
   tmp_c(ixO^S,1:ncons,2)=min(0.0d0, lambdaL(ixO^S,1:ncons,2), lambdaR(ixO^S,1:ncons,2) )

   select case(method)
   case('tvdlf','tvdlf1','fd_tvdlf')
   cmax(ixO^S,1:ncons) = max(dabs(tmp_c(ixO^S,1:ncons,1)), dabs(tmp_c(ixO^S,1:ncons,2)))
   cmin(ixO^S,1:ncons) = 0.0d0 !useless
   case('hll','hll1','fd_hll')
   cmax(ixO^S,1:ncons) = tmp_c(ixO^S,1:ncons,1)
   cmin(ixO^S,1:ncons) = tmp_c(ixO^S,1:ncons,2)
   case default
      call mpistop('Pls specify FV method for cbounds')
   end select
  

   ! test------------
   ! M1_cbounds test, set to 1
   {#IFDEF M1_SPEED1
   {^KSP& 
   cmax(ixI^S,nrad^KSP_) = 1.0d0
   cmax(ixI^S,erad^KSP_) = 1.0d0
   {^C& cmax(ixI^S,frad^KSP^C_) = 1.0d0 \}
   
   cmin(ixI^S,nrad^KSP_) = -1.0d0
   cmin(ixI^S,erad^KSP_) = -1.0d0
   {^C& cmin(ixI^S,frad^KSP^C_) = -1.0d0 \}
   \}
   }
  !-----------------

  end subroutine get_cbounds

  subroutine get_lambda(ixI^L, ixO^L, idims, w, x, lambda, qtC, metricM1) !KEN
    use mod_cfc_parameters
    use mod_imhd_intermediate
    use mod_metric
    use mod_eos
    use mod_m1_metric_interface !KEN
    {#IFDEF M1 
    use mod_m1, only: m1_get_wavespeeds 
    }
    include 'amrvacdef.f'
    integer, intent(in)                       :: ixI^L, ixO^L, idims
    double precision, intent(inout)              :: w(ixI^S, 1:nw)
    double precision, intent(in)              :: qtC
    double precision, intent(in)              :: x(ixI^S, 1:ndim)
    double precision, intent(out)             :: lambda(ixI^S, 1:ncons, 1:2)
    type(m1_metric_helper), intent(in)   :: metricM1 !KEN

    !local quantities
    double precision                          :: VdotB(ixI^S), bD(ixI^S,1:^NC), v2(ixI^S), b2(ixI^S), htot(ixI^S)
    double precision                          :: ca2(ixI^S), cs2(ixI^S), vel(ixI^S), inv_gamma_ii(ixI^S), W2v2(ixI^S)
    double precision                          :: clight(ixI^S), a2(ixI^S), tmp1(ixI^S), tmp2(ixI^S), tmp3(ixI^S)
    double precision                          :: gamma(ixI^S,1:3,1:3), dummy, eps_tmp(ixI^S), prs_tmp(ixI^S)
    integer                                   :: idir, ix^D
    double precision                          :: cs2_local(ixI^S), lfac(ixI^S), temp_local, rho_local, alp_prime(ixI^S)
    {#IFDEF UNITS_TESTS
    double precision :: TESTqtC = 1.0d+10
    }

    if (.not. useprimitiveRel) call mpistop('Stop, the code needs u = Wv instead of v')

    {do ix^D = ixO^LIM^D \}
        {#IFNDEF DY_SP
          prs_tmp(ixO^S) = w(ixO^S, pp_)
        }
        if (eos_type == tabulated) then
          temp_local = w(ix^D,T_eps_)
          rho_local  = w(ix^D,rho_)
          call eos_temp_get_all_one_grid(rho_local,temp_local,w(ix^D,ye_),&
                                         eps_tmp(ix^D),prs=prs_tmp(ix^D),cs2=cs2_local(ix^D))
        else
          eps_tmp(ix^D) = w(ix^D, T_eps_)
          call eos_get_cs2_one_grid(cs2_local(ix^D),w(ix^D,rho_),eps_tmp(ix^D))
          {#IFDEF DY_SP
            call eos_get_pressure_one_grid(prs_tmp(ix^D),w(ix^D,rho_),eps_tmp(ix^D))
          }
        endif
    {enddo^D&\}
   
{#IFDEF DY_SP
     call imhd_get_intermediate_variables(ixI^L, ixO^L, w(ixI^S, 1:nw), x(ixI^S, 1:ndim), &
                   gamma=gamma(ixI^S,1:3,1:3), lfac=lfac(ixI^S), &
                   v2=v2(ixI^S), &
                   b2=b2(ixI^S))

     htot(ixO^S) = 1.0d0 + eps_tmp(ixO^S) &
              + ( prs_tmp(ixO^S) + b2(ixO^S) ) / w(ixO^S, rho_)

     ! beware of coordinate singularities
     where ( dabs(gamma(ixO^S, idims, idims)) > smalldouble )
        ! sound speed square
        cs2(ixO^S) = cs2_local(ixO^S)
        ! Calculate Alfven speed
        ca2(ixO^S) = b2(ixO^S) / ( w(ixO^S, rho_) * htot(ixO^S) )
        vel(ixO^S) = w(ixO^S, u0_+idims) / lfac(ixO^S)
        inv_gamma_ii(ixO^S) = 1.0d0 / gamma(ixO^S, idims, idims)
     else where
        cs2(ixO^S) = 0.0d0
        ca2(ixO^S) = 0.0d0
        vel(ixO^S) = 0.0d0
        inv_gamma_ii(ixO^S) = 0.0d0
        v2(ixO^S) = 0.0d0
        b2(ixO^S) = 0.0d0
     end where

}

{#IFNDEF DY_SP
     b2(ixO^S)      = 0.0d0
     VdotB(ixO^S)   = 0.0d0 
     W2v2(ixO^S)    = 0.0d0
     htot(ixO^S)    = 0.0d0
  
     call square3u(ixI^L,ixO^L,myM,w(ixI^S,u1_:u^NC_),W2v2)

     call square3u(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),b2)

     call lower3(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),bD)

     v2(ixO^S)      = W2v2(ixO^S) / w(ixO^S, lfac_)**2.0d0
     VdotB(ixO^S)   = ({^C& bD(ixO^S,^C)*w(ixO^S,u^C_)+})/w(ixO^S,lfac_)
     b2(ixO^S)      = b2(ixO^S)/w(ixO^S, lfac_)**2.0d0 + VdotB(ixO^S)**2.0d0

     ! Calculate the magnetically modified specific enthalpy
     htot(ixO^S)    = 1.0d0 + eps_tmp(ixO^S) &
                      + ( prs_tmp(ixO^S) + b2(ixO^S) ) / w(ixO^S, rho_)

     where ( dabs(myM%g(idims,idims)%elem(ixO^S)) > smalldouble )
        cs2(ixO^S) = cs2_local(ixO^S)
        ! Calculate Alfven speed
        ca2(ixO^S) = b2(ixO^S) / ( w(ixO^S, rho_) * htot(ixO^S) )
        vel(ixO^S) = w(ixO^S, u0_+idims) / w(ixO^S, lfac_)
        inv_gamma_ii(ixO^S) = myM%gammainv(idims,idims)%elem(ixO^S)
     else where
        cs2(ixO^S) = 0.0d0
        ca2(ixO^S) = 0.0d0
        vel(ixO^S) = 0.0d0
        inv_gamma_ii(ixO^S) = 0.0d0
        v2(ixO^S) = 0.0d0
        b2(ixO^S) = 0.0d0
     end where
}
     ! speed of light
     clight(ixO^S) = dsqrt( inv_gamma_ii(ixO^S) )

     ! when using GLM, we need to include two additional modes with speed of
     ! light, to without violating causality
{#IFDEF GLM
        lambda(ixO^S,D_,1) =  clight(ixO^S)
        lambda(ixO^S,D_,2) = -clight(ixO^S)
}
{#IFNDEF GLM
        ! upper bound for the fast wave speed
        a2(ixO^S) = cs2(ixO^S) + ca2(ixO^S) - ca2(ixO^S) * cs2(ixO^S)

        tmp1(ixO^S) = ( 1.0d0 - a2(ixO^S) ) * vel(ixO^S)
        tmp2(ixO^S) = ( 1.0d0 - v2(ixO^S) * a2(ixO^S) )
        tmp3(ixO^S) = ( a2(ixO^S) * ( 1.0d0 - v2(ixO^S) ) * &
                        ( tmp2(ixO^S) * inv_gamma_ii(ixO^S) &
                          - ( 1.0d0 - a2(ixO^S) ) * vel(ixO^S)**2 ) )

        where ( tmp2(ixO^S)==0.0d0 .or. tmp3(ixO^S)<0.0d0 )
           lambda(ixO^S,D_,1) = 0.0d0
           lambda(ixO^S,D_,2) = 0.0d0
        else where
           tmp3(ixO^S) = dsqrt( tmp3(ixO^S) )
           lambda(ixO^S,D_,1) = ( tmp1(ixO^S) + tmp3(ixO^S) ) / tmp2(ixO^S)
           lambda(ixO^S,D_,2) = ( tmp1(ixO^S) - tmp3(ixO^S) ) / tmp2(ixO^S)
        end where

        ! limit with speed of light
        clight(ixO^S) = dsqrt( inv_gamma_ii(ixO^S) )
        lambda(ixO^S,D_,1) = max( min( lambda(ixO^S,D_,1), clight(ixO^S) ), -clight(ixO^S) )
        lambda(ixO^S,D_,2) = max( min( lambda(ixO^S,D_,2), clight(ixO^S) ), -clight(ixO^S) )
}

{#IFDEF DY_SP
     if (use_gw_br  .and. use_h00_coupling) then
        alp_prime(ixO^S) = dsqrt( w(ixO^S,alp_metric_)**2 - w(ixO^S,h_00_) )
        !alp_prime(ixO^S) = w(ixO^S,alp_metric_) - w(ixO^S,h_00_)/2.0d0/w(ixO^S,alp_metric_)
     else
        alp_prime(ixO^S) = w(ixO^S,alp_metric_)
     endif

     lambda(ixO^S,D_,1) = alp_prime(ixO^S) * lambda(ixO^S,D_,1) - w(ixO^S, alp_metric_+idims)
     lambda(ixO^S,D_,2) = alp_prime(ixO^S) * lambda(ixO^S,D_,2) - w(ixO^S, alp_metric_+idims)
}
{#IFNDEF DY_SP
     lambda(ixO^S,D_,1) = myM%alpha(ixO^S) * lambda(ixO^S,D_,1) - myM%beta(idims)%elem(ixO^S) 
     lambda(ixO^S,D_,2) = myM%alpha(ixO^S) * lambda(ixO^S,D_,2) - myM%beta(idims)%elem(ixO^S) 
}

     ! for magnetised fluid, use the same lambda
     do idir = d_+1, ncons ! start from d_+1
        lambda(ixO^S,idir,1) = lambda(ixO^S,D_,1)
        lambda(ixO^S,idir,2) = lambda(ixO^S,D_,2)
     end do

  {#IFDEF M1 
   call m1_get_wavespeeds(ixI^L,ixO^L, idims, w, x, lambda, qtC, metricM1) !KEN added metricM1
  } 
  end subroutine get_lambda

  subroutine getcmax(w, x, ixI^L, ixO^L, idims, cmax)
    use mod_m1_metric_interface !KEN
    include 'amrvacdef.f'
    integer, intent(in)                       :: ixI^L, ixO^L, idims
    double precision, intent(in)              :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    integer                                   :: iw, ix^D
    double precision                          :: lambda(ixI^S,1:ncons,1:2)
    type(m1_metric_helper)   :: metricM1 !KEN

    double precision :: qtC
    qtC = -2.0
    {#IFDEF M1
      call metricM1%fill_metric(w,x,ixI^L, ixO^L)!KEN
    }

    call get_lambda(ixI^L, ixO^L, idims, w(ixI^S, 1:nw), x(ixI^S, 1:ndim), lambda(ixI^S,1:ncons,1:2), qtC, metricM1)
    cmax(ixO^S) = 0.0d0
    do iw = 1, ncons
       cmax(ixO^S) = max(cmax(ixO^S), abs(lambda(ixO^S,iw,1)), abs(lambda(ixO^S,iw,2)))
    end do

   !! test------------
   ! M1_cbounds test
   !cmax(ixO^S) = 1.0d0
   {#IFDEF M1_SPEED1
   cmax(ixI^S) = 1.0d0
   }
   !!-----------------
   {#IFDEF M1
      call metricM1%destroy()!KEN
    }

  end subroutine getcmax


  subroutine getcmax_old(w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)

    ! Calculate cmax_idim within ixO^L
    use mod_imhd_intermediate
    use mod_eos
    use mod_metric, only: square3u
    use mod_cfc_parameters
    include 'amrvacdef.f'

    logical, intent(in)               :: needcmin
    integer, intent(in)               :: ixI^L,ixO^L,idims
    double precision, intent(inout)   :: w(ixI^S,1:nw)
    double precision                  :: wprim_tmp(ixI^S,1:nw), lrho(ixI^S), lye(ixI^S)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    double precision, intent(out)     :: cmax(ixI^S),cmin(ixI^S)

    ! .. local ..
    double precision, dimension(ixI^S)  :: sqrB,csound2
    double precision, dimension(ixI^S)  :: A,B,C
    double precision, dimension(ixI^S)  :: B2
    double precision, dimension(ixI^S)  :: sUidims
    ! 3-metric gamma_ij and gammaij
    double precision                    :: gamma(ixI^S,1:3,1:3), lfac(ixI^S), h_th(ixI^S), xi(ixI^S), eps_tmp(ixI^S), prs_tmp(ixI^S)
    double precision                    :: gammainv(ixI^S,1:3,1:3), cs2_local(ixI^S), temp_local, rho_local, alp_prime(ixI^S), W2v2(ixI^S)
    integer                             :: j, ix^D, idir
    !-----------------------------------------------------------------------------

    {do ix^D = ixO^LIM^D \}
        if (eos_type == tabulated) then
          temp_local = w(ix^D,T_eps_)
          rho_local  = w(ix^D,rho_)
          call eos_temp_get_all_one_grid(rho_local,temp_local,w(ix^D,ye_),&
                                         eps_tmp(ix^D),prs=prs_tmp(ix^D),cs2=cs2_local(ix^D))
        else
          eps_tmp(ix^D) = w(ix^D, T_eps_)
          call eos_get_cs2_one_grid(cs2_local(ix^D),w(ix^D,rho_),eps_tmp(ix^D))
          {#IFNDEF DY_SP
             prs_tmp(ix^D) = w(ix^D, pp_)
          }
          {#IFDEF DY_SP
            call eos_get_pressure_one_grid(prs_tmp(ix^D),w(ix^D,rho_),eps_tmp(ix^D))
          }
        endif
    {enddo^D&\}


    ! squared Eulerian field strength (big B^2)
{#IFNDEF DY_SP
    call square3u(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),B2)
    lfac(ixO^S) = w(ixO^S, lfac_)
    ! C = VdotB
    C(ixO^S)= ({^C& w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/w(ixO^S,xi_)
}
{#IFDEF DY_SP
    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
    end do

    W2v2(ixO^S) = {^C& gamma(ixO^S,^C,^C)*w(ixO^S, u^C_)**2 +}
    lfac(ixO^S) = dsqrt( W2v2(ixO^S) + 1.0d0 )


    call get_gammainvij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gammainv(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gammainv(ixO^S,idir,idir) = gammainv(ixO^S,idir,idir) / w(ixO^S, psi_metric_)**4
    end do

    h_th(ixO^S) = 1.0d0 + eps_tmp(ixO^S) + prs_tmp(ixO^S) / w(ixO^S, rho_)

    B2(ixO^S) = {^C& gamma(ixO^S,^C,^C) * w(ixO^S,b^C_)**2 +}
    xi(ixO^S) = lfac(ixO^S)**2 * w(ixO^S, rho_) * h_th(ixO^S)
    ! C = VdotB
    C(ixO^S) = ({^C& w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/xi(ixO^S)
}

    ! squared comoving field strength (little b^2)
    sqrB(ixO^S) = B2(ixO^S)/lfac(ixO^S)**2 + C(ixO^S)**2

    csound2(ixO^S) = cs2_local(ixO^S)

    ! contravariant momentum used for velocity
    ! get only the required component:
    sUidims(ixO^S) = zero
{#IFNDEF DY_SP 
    do j=1,^NC
       sUidims(ixO^S) = sUidims(ixO^S) &
            + myM%gammainv(idims,j)%elem(ixO^S)*w(ixO^S,s0_+j)
    end do
}
{#IFDEF DY_SP
    do j=1,^NC
       sUidims(ixO^S) = sUidims(ixO^S) &
            + w(ixO^S,s0_+j) * gammainv(ixO^S,idims,j)
    end do
}

    select case(typepoly)

    case('gammie')

       A(ixO^S)   = xi(ixO^S)/lfac(ixO^S)**2 ! rhoh
       ! squared Alfven speed
       B(ixO^S)  = sqrB(ixO^S)/(sqrB(ixO^S)+A(ixO^S))
       ! equation 72 Del zanna et al. 2007A&A...473...11D
       A(ixO^S)  = csound2(ixO^S) + B(ixO^S) - csound2(ixO^S)*B(ixO^S)

       sqrB(ixO^S)     = one - one/lfac(ixO^S)**2 !v2
       csound2(ixO^S)  = (sUidims(ixO^S) + C(ixO^S) * w(ixO^S,b0_+idims)) &
            / ( xi(ixO^S)+B2(ixO^S) ) ! vidim

       ! reuse B for square root:
{#IFNDEF DY_SP
      B(ixO^S) = dsqrt( A(ixO^S) * (one-sqrB(ixO^S)) &
           * ( (one-sqrB(ixO^S)*A(ixO^S)) * myM%gammainv(idims,idims)%elem(ixO^S) &
           - csound2(ixO^S)**2 * (one-A(ixO^S)) ) )
}
{#IFDEF DY_SP
       B(ixO^S) = dsqrt( A(ixO^S) * (one-sqrB(ixO^S)) &
            * ( (one-sqrB(ixO^S)*A(ixO^S)) * gammainv(ixO^S,idims,idims) &
            - csound2(ixO^S)**2 * (one-A(ixO^S)) ) )
}

       ! reuse C for first term:
       C(ixO^S) = csound2(ixO^S)*(one-A(ixO^S))

       ! characteristic speed Eq. 76 Del zanna et al. 2007A&A...473...11D
       cmax(ixO^S) = ( C(ixO^S) + B(ixO^S) ) &
            / ( one - sqrB(ixO^S)*A(ixO^S) )
       cmin(ixO^S) = ( C(ixO^S) - B(ixO^S) ) &
            / ( one - sqrB(ixO^S)*A(ixO^S) )

    case default

       call mpistop("getcmax old: Unknown polynomial type")

    end select

{#IFNDEF DY_SP    
    ! Limit by speed of light:
    cmin(ixO^S) = max(cmin(ixO^S), - 1.0d0/sqrt(myM%g(idims,idims)%elem(ixO^S)))
    cmin(ixO^S) = min(cmin(ixO^S),   1.0d0/sqrt(myM%g(idims,idims)%elem(ixO^S)))
    cmax(ixO^S) = max(cmax(ixO^S), - 1.0d0/sqrt(myM%g(idims,idims)%elem(ixO^S)))
    cmax(ixO^S) = min(cmax(ixO^S),   1.0d0/sqrt(myM%g(idims,idims)%elem(ixO^S)))
}
{#IFDEF DY_SP
    ! Limit by speed of light:
    cmin(ixO^S) = max(cmin(ixO^S), - 1.0d0/sqrt(gamma(ixO^S,idims,idims) ))
    cmin(ixO^S) = min(cmin(ixO^S),   1.0d0/sqrt(gamma(ixO^S,idims,idims) ))
    cmax(ixO^S) = max(cmax(ixO^S), - 1.0d0/sqrt(gamma(ixO^S,idims,idims) ))
    cmax(ixO^S) = min(cmax(ixO^S),   1.0d0/sqrt(gamma(ixO^S,idims,idims) ))
}

{#IFNDEF DY_SP 
   ! add shift-part and mutliply by lapse
   if (.not. needcmin) then
      cmax(ixO^S) = max( abs(myM%alpha(ixO^S)*cmax(ixO^S) - myM%beta(idims)%elem(ixO^S)), &
           abs(myM%alpha(ixO^S)*cmin(ixO^S) - myM%beta(idims)%elem(ixO^S)) )
   else
      cmax(ixO^S) = myM%alpha(ixO^S)*cmax(ixO^S) - myM%beta(idims)%elem(ixO^S)
      cmin(ixO^S) = myM%alpha(ixO^S)*cmin(ixO^S) - myM%beta(idims)%elem(ixO^S)
   end if
}
{#IFDEF DY_SP
    if (use_gw_br  .and. use_h00_coupling) then
       alp_prime(ixO^S) = dsqrt( w(ixO^S,alp_metric_)**2 - w(ixO^S,h_00_) )
       !alp_prime(ixO^S) = w(ixO^S,alp_metric_) - w(ixO^S,h_00_)/2.0d0/w(ixO^S,alp_metric_)
    else
       alp_prime(ixO^S) = w(ixO^S,alp_metric_)
    endif

    ! add shift-part and mutliply by lapse
    if (.not. needcmin) then
       cmax(ixO^S) = max( abs(alp_prime(ixO^S)*cmax(ixO^S) - w(ixO^S,alp_metric_+idims) ), &
            abs(alp_prime(ixO^S)*cmin(ixO^S) - w(ixO^S,alp_metric_+idims) ) )
    else
       cmax(ixO^S) = alp_prime(ixO^S)*cmax(ixO^S) - w(ixO^S,alp_metric_+idims) ! alp_metric+idims is beta(idims)
       cmin(ixO^S) = alp_prime(ixO^S)*cmin(ixO^S) - w(ixO^S,alp_metric_+idims)
    end if
}



  end subroutine getcmax_old
  !=============================================================================
  subroutine getflux(w,wprim,x,ixI^L,ixO^L,idims,f,transport, qtC, metric_M1) !KEN

    ! Calculate non-transport flux f_idim[iw] within ixO^L.
    use mod_imhd_intermediate
    use mod_metric
    use mod_eos
    use mod_cfc_parameters
    use mod_m1_metric_interface !KEN
    {#IFDEF M1
    use mod_m1
    }
    include 'amrvacdef.f'
    integer, intent(in)                :: ixI^L,ixO^L,idims
    double precision, intent(in)       :: qtC
    double precision, intent(in)       :: w(ixI^S,1:nw)     ! conserved variables
    double precision, intent(inout)       :: wprim(ixI^S,1:nw) ! primitive variables
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: f(ixI^S,1:nwflux)
    logical, intent(out)               :: transport(1:nwflux)
    type(m1_metric_helper), intent(in)   :: metric_M1 !KEN

    ! .. local ..
    double precision, dimension(ixI^S) :: VdotB, sqrB, ptot, v, xi, h_th, prs_tmp, eps_tmp, alp_prime
    double precision, dimension(ixI^S,1:ndir)  :: alphabjU, alphabjD, vU
    ! 3-metric gamma_ij
    double precision                           :: gamma(ixI^S,1:3,1:3), lfac(ixI^S), temp_local, rho_local, W2v2(ixI^S)
    integer                                    :: i, iw, idir, ix^D
    !-----------------------------------------------------------------------------
    if (evolve_hydro) then
       transport=.true.
    else
       return
    endif

    {do ix^D = ixO^LIM^D \}
        if (eos_type == tabulated) then
          temp_local = wprim(ix^D,T_eps_)
          rho_local  = wprim(ix^D,rho_)
          call eos_temp_get_all_one_grid(rho_local,temp_local,wprim(ix^D,ye_),&
                                         eps_tmp(ix^D),prs=prs_tmp(ix^D))
        else
          eps_tmp(ix^D) = wprim(ix^D, T_eps_)
          call eos_get_pressure_one_grid(prs_tmp(ix^D),wprim(ix^D,rho_),eps_tmp(ix^D))
        endif
    {enddo^D&\}

{#IFDEF DY_SP
    if (use_gw_br  .and. use_h00_coupling) then
       alp_prime(ixO^S) = dsqrt( wprim(ixO^S,alp_metric_)**2 - wprim(ixO^S,h_00_) )
       !alp_prime(ixO^S) = wprim(ixO^S,alp_metric_) - wprim(ixO^S,h_00_)/2.0d0/wprim(ixO^S,alp_metric_)
    else
       alp_prime(ixO^S) = wprim(ixO^S,alp_metric_)
    endif

    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
    do idir = 1, ndir
       gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * wprim(ixO^S, psi_metric_)**4
    end do

    W2v2(ixO^S) = {^C& gamma(ixO^S,^C,^C)*wprim(ixO^S, u^C_)**2 +}
    lfac(ixO^S) = dsqrt( W2v2(ixO^S) + 1.0d0 )

    h_th(ixO^S) = 1.0d0 + eps_tmp(ixO^S) + prs_tmp(ixO^S) / wprim(ixO^S, rho_)


    xi(ixO^S)    = lfac(ixO^S)**2 * wprim(ixO^S,rho_) * h_th(ixO^S)
    VdotB(ixO^S) = ({^C& w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/xi(ixO^S)
    sqrB(ixO^S)  = {^C& gamma(ixO^S,^C,^C)*w(ixO^S,b^C_)**2 +}

    do i=1, ^NC
       vU(ixO^S,i) = wprim(ixO^S,u0_+i)/lfac(ixO^S)
       alphabjU(ixO^S,i) = alp_prime(ixO^S) * (w(ixO^S,b0_+i)/lfac(ixO^S)**2 &
            + VdotB(ixO^S)*vU(ixO^S,i) )
    end do

    do i=1, ^NC
       alphabjD(ixO^S,i) = alphabjU(ixO^S,i) * gamma(ixO^S,i,i)
    enddo
}

{#IFNDEF DY_SP
    xi(ixO^S) = w(ixO^S,xi_)
    VdotB(ixO^S) = ({^C& w(ixO^S,s^C_)*w(ixO^S,b^C_)+})/xi(ixO^S)
    call square3u(ixI^L,ixO^L,myM,w(ixI^S,b1_:b^NC_),sqrB(ixI^S))
    lfac(ixO^S) = wprim(ixO^S, lfac_)

    do i=1, ^NC
       vU(ixO^S,i) = wprim(ixO^S,u0_+i)/lfac(ixO^S)
       alphabjU(ixO^S,i) = myM%alpha(ixO^S) * (w(ixO^S,b0_+i)/lfac(ixO^S)**2 &
            + VdotB(ixO^S)*vU(ixO^S,i) )
    end do
    call lower3(ixI^L,ixO^L,myM,alphabjU,alphabjD)
}

    ptot(ixO^S) = half*(VdotB(ixO^S)**2  & 
         + sqrB(ixO^S)/(lfac(ixO^S)**2)) + prs_tmp(ixO^S)


! ------  prims  ------- !
        {#IFNDEF DY_SP
        f(ixO^S, pp_)     =  0.0d0
        transport(pp_)    = .false.
        }

        f(ixO^S, T_eps_)  =  0.0d0
        transport(T_eps_) = .false.
        f(ixO^S, rho_)    =  0.0d0
        transport(rho_)   = .false.
{^C&    f(ixO^S, u^C_)    =  0.0d0
        transport(u^C_)   = .false.  \}

        if (eos_type == tabulated) then
          f(ixO^S, ye_)     =  0.0d0
          transport(ye_)    = .false.
! ------  cons  ------- !
          f(ixO^S, Dye_)    = 0.0d0
        endif

        f(ixO^S, d_)      = 0.0d0

{#IFDEF TRACER
{^FL&       f(ixO^S, tr^FL_)  = 0.0d0 \}
} 
       
{#IFDEF EPSINF
        f(ixO^S,epsinf_)  = 0.0d0
        f(ixO^S,rho0_)    = 0.0d0
        f(ixO^S,rho1_)    = 0.0d0
        f(ixO^S,n0_)      = 0.0d0
        f(ixO^S,n_)       = 0.0d0
}

{#IFDEF ENTROPY
        f(ixO^S,Ds_)      = 0.0d0
}

{#IFDEF ELECTRONS
        f(ixO^S,Dse_)     = 0.0d0
}
{#IFDEF M1
  {^KSP&
      f(ixO^S,nrad^KSP_)  = 0.0d0
	   f(ixO^S,erad^KSP_)  = 0.0d0
	   {^C& f(ixO^S,frad^KSP^C_)  = 0.0d0 \}
      transport(nrad^KSP_)    = .false.
      transport(erad^KSP_)    = .false.
      {^C& transport(frad^KSP^C_)    = .false. \} 
  \}
}          

         

      !  psi
      {#IFDEF GLM
      !!!f_i[psi]=B**i - psi beta**i
      ! Eq. 24e Dedner et al 2002 JCP, 175, 645
          {#IFNDEF DY_SP
          f(ixO^S,psi_)=w(ixO^S,b0_+idims) - w(ixO^S,psi_)*myM%beta(idims)%elem(ixO^S)
          }
          {#IFDEF DY_SP
            f(ixO^S,psi_)=w(ixO^S,b0_+idims) - w(ixO^S,psi_)*wprim(ixO^S,alp_metric_+idims) ! this is beta_metric^C                                                                    instead of alp_metric
          }

      }


     !  b^C

     {^C&    !!! F^i[B_c] = B_c*v_i - v_c*B_i 
          if (idims==^C) then
             ! f_i[b_i] should be exactly 0, so we do not use the transport flux
             f(ixO^S,b^C_)   =zero
             transport(b^C_) =.false.
          else
             ! transport velocity in direction ^C
             call getv(wprim,x,ixI^L,ixO^L,^C,v)
             ! multiply by B in direction idims
             f(ixO^S,b^C_) = - v(ixO^S) * w(ixO^S,b0_+idims)
          end if
          {#IFDEF GLM
             {#IFNDEF DY_SP
!!! - B**i beta**j / alpha
          f(ixO^S,b^C_) = f(ixO^S,b^C_) - w(ixO^S,b0_+idims)*myM%beta(^C)%elem(ixO^S)/myM%alpha(ixO^S)
             }
             {#IFDEF DY_SP
!!! - B**i beta**j / alpha
          f(ixO^S,b^C_) = f(ixO^S,b^C_) - w(ixO^S,b0_+idims)*wprim(ixO^S,beta_metric^C_)/alp_prime(ixO^S)
             }
          }

          \}

     !  tau

          where(xi(ixO^S)<smallxi)
             f(ixO^S,e_)=zero
          elsewhere
             ! f=ptot*v^i
             f(ixO^S,e_) = ptot(ixO^S)*vU(ixO^S,idims)
             ! f=ptot*v^i - v.B*B^i
             f(ixO^S,e_) = f(ixO^S,e_) - VdotB(ixO^S)*w(ixO^S,b0_+idims)
             ! f=(ptot*v^i - v.B*B^i)*alpha
             {#IFNDEF DY_SP
             f(ixO^S,e_) = f(ixO^S,e_) * myM%alpha(ixO^S)
             }
             {#IFDEF DY_SP
             f(ixO^S,e_) = f(ixO^S,e_) * alp_prime(ixO^S)
             }

          end where


    !   s^C
       {^C& ! momentum is covariant, transport flux done.
          ! f=-bj Bi /lfac

          f(ixO^S,s^C_) = - alphabjD(ixO^S,^C) * w(ixO^S,b0_+idims)

          if (idims==^C) then !!! add ptot
             {#IFNDEF DY_SP
             f(ixO^S,s^C_)=f(ixO^S,s^C_) + ptot(ixO^S) * myM%alpha(ixO^S)
             }
             {#IFDEF DY_SP
             f(ixO^S,s^C_)=f(ixO^S,s^C_) + ptot(ixO^S) * alp_prime(ixO^S)
             }
          end if
          \}


    !> M1-fluxes:
    {#IFDEF M1
	  call m1_get_fluxes(ixI^L, ixO^L, idims, x, wprim, f, qtC, metric_M1) !KEN
	 }
   
  end subroutine getflux
  !=============================================================================
  subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,wold,x,qtC,metric_M1)

    ! Add geometrical source terms to w
    use mod_metric, only: lower3, raise3, dalphadj_is_zero, beta_is_zero, dbetaidj_is_zero, dgdk_is_zero, square3u
    use mod_full_gr_source
    use mod_m1_metric_interface
    {#IFDEF M1 
    use mod_m1
    } 
    include 'amrvacdef.f'

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qdt
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in)       :: wold(ixI^S,1:nw)
    double precision, intent(in)       :: qtC   !TestM1
	 type(m1_metric_helper), intent(in) :: metric_M1 
    ! .. local ..
    integer                            :: iw, i,j,k, inonzero, jnonzero, hxO^L, idims
    double precision,dimension(ixI^S,1:ndir)  :: vU, bsU, bsD, sU, wikdjgammaik
    double precision,dimension(ixI^S)         :: VdotB, sqrB, tmp, tmp2, wikbetajdjgammaik
    double precision,dimension(ixI^S)         :: ptot, lrho, lye
    integer :: ix^D
    logical :: cfc_geosource_terms
	{#IFDEF M1 ! m1:
	integer ::  pressind
    double precision :: symfact
     !double precision, dimension(ixI^S) :: rhs
    double precision, dimension(ixI^S,m1_npress,^NS) :: P
     !double precision, dimension(ixI^S,^NC)           :: FU
	double precision, dimension(ixI^S,^NC,^NC)       :: Kij_local 
	double precision, dimension(ixI^S,^NC,^NC,^NC)       :: d_k_gamma_ij
    double precision, dimension(ixI^S)       :: m1_tmp1
	double precision, dimension(ixI^S)       :: m1_tmp2
	double precision, dimension(ixI^S)       :: m1_tmp3
    !double precision, dimension(ixI^S) ::  Gamma_updt
	} 
    !-----------------------------------------------------------------------------

    {#IFDEF DY_SP
       !cfc_geosource_terms = .false.    
       cfc_geosource_terms = .true.    
    }
    {#IFNDEF DY_SP
       cfc_geosource_terms = .false.    
    }


! for DY_SP
    if (cfc_geosource_terms) then
       !call addgeometry_source_cowling_cfc(qdt,ixI^L,ixO^L,wCT,w,x)
       call addgeometry_source_full_gr(qdt,ixI^L,ixO^L,wCT,w,wold,x,qtC, metric_M1)
       return
    endif

!  below is non DY_SP
    VdotB(ixO^S)= ({^C& wCT(ixO^S,s^C_)*wCT(ixO^S,b^C_)+})/wCT(ixO^S,xi_)
    call square3u(ixI^L,ixO^L,myM,wCT(ixI^S,b1_:b^NC_),sqrB)
    call raise3(ixI^L,ixO^L,myM,wCT(ixI^S,s1_:s^NC_),sU)

    {^C&
         vU(ixO^S,^C)=(sU(ixO^S,^C)+VdotB(ixO^S)*wCT(ixO^S,b0_+^C))/ &
         (wCT(ixO^S,xi_)+sqrB(ixO^S))
    bsU(ixO^S,^C) = wCT(ixO^S,b0_+^C)/wCT(ixO^S,lfac_) + wCT(ixO^S,lfac_)*VdotB(ixO^S)*vU(ixO^S,^C)
    \}
    call lower3(ixI^L,ixO^L,myM,bsU,bsD)

    ptot(ixO^S) = wCT(ixO^S,pp_)
    ptot(ixO^S) = half*(VdotB(ixO^S)**2  & 
         + sqrB(ixO^S)/wCT(ixO^S,lfac_)**2)  &
         + ptot(ixO^S)


  if (evolve_hydro) then

    ! wikdjgammaik is without the total pressure term
    wikdjgammaik(ixO^S,:) = zero
    do inonzero = 1, myM%nnonzeroDgDk
       i = myM%nonzeroDgDk(inonzero)%i
       k = myM%nonzeroDgDk(inonzero)%j
       j = myM%nonzeroDgDk(inonzero)%k

       tmp(ixO^S) = sU(ixO^S,i)*vU(ixO^S,k) &
            !        + ptot(ixO^S)*myM%gammainv(i,k)%elem(ixO^S) &
            - bsU(ixO^S,i)*wCT(ixO^S,b0_+k)/wCT(ixO^S,lfac_)
       wikdjgammaik(ixO^S,j) = wikdjgammaik(ixO^S,j) + tmp(ixO^S) * myM%DgDk(i,k,j)%elem(ixO^S)
    end do


    wikbetajdjgammaik(ixO^S) = zero
    do inonzero = 1, myM%nnonzeroDgDk
       i = myM%nonzeroDgDk(inonzero)%i
       k = myM%nonzeroDgDk(inonzero)%j
       j = myM%nonzeroDgDk(inonzero)%k
       
       if  (beta_is_zero(j)) cycle

       tmp(ixO^S) = sU(ixO^S,i)*vU(ixO^S,k) &
            + ptot(ixO^S)*myM%gammainv(i,k)%elem(ixO^S) &
            - bsU(ixO^S,i)*wCT(ixO^S,b0_+k)/wCT(ixO^S,lfac_)

       wikbetajdjgammaik(ixO^S) = wikbetajdjgammaik(ixO^S) &
            + tmp(ixO^S) * myM%beta(j)%elem(ixO^S) * myM%DgDk(i,k,j)%elem(ixO^S)
    end do
    

!   s^C
          {^C&
                                ! s[s^C_] = 1/2 alpha W**ik d^Cdgammaik + S_i d^Cbeta**i -U d^Calpha
          tmp(ixO^S) = half * myM%alpha(ixO^S) * wikdjgammaik(ixO^S,^C)
          ! Treat the total pressure separately to make discretized gradient:
          hxOmin^D=ixOmin^D-kr(^D,^C);hxOmax^D=ixOmax^D-kr(^D,^C);
          idims = ^C
          select case(idims)
             {case(^D)
             tmp(ixO^S) = tmp(ixO^S) + myM%alpha(ixO^S)*ptot(ixO^S) &
                  *(mygeo%surfaceC^D(ixO^S)-mygeo%surfaceC^D(hxO^S)) &
                  /mygeo%dvolume(ixO^S)\}
          end select

          do i = 1, ^NC
             if (dbetaidj_is_zero(i,^C)) cycle
             tmp(ixO^S) = tmp(ixO^S) + wCT(ixO^S,s0_+i) * myM%dbetaiDj(i,^C)%elem(ixO^S)
          end do

          if (.not. dalphadj_is_zero(^C)) then
             tmp(ixO^S) = tmp(ixO^S) - (wCT(ixO^S,tau_)+wCT(ixO^S,d_))*myM%dalphaDj(^C)%elem(ixO^S)
          end if


          w(ixO^S,s^C_) = w(ixO^S,s^C_) + qdt*tmp(ixO^S)

          \}


!  tau
          ! s[tau_] = 1/2 W**ik beta**j dgammaikdj + W**j_i*dbetaidj - S**j dalphadj
          tmp(ixO^S) = half*wikbetajdjgammaik(ixO^S)
          do inonzero = 1, myM%nnonzeroDbetaiDj
             i = myM%nonzeroDbetaiDj(inonzero)%i
             j = myM%nonzeroDbetaiDj(inonzero)%j
             tmp2(ixO^S) = vU(ixO^S,j)*wCT(ixO^S,s0_+i)
             if (i .eq. j) &
                  tmp2(ixO^S) = tmp2(ixO^S) + ptot(ixO^S)
             tmp2(ixO^S) = tmp2(ixO^S) - bsD(ixO^S,i)*wCT(ixO^S,b0_+j)/wCT(ixO^S,lfac_)
             !            print*, 'Source:',i,j,maxval(abs(tmp2(ixO^S)))
             tmp(ixO^S) = tmp(ixO^S) + tmp2(ixO^S)  * myM%nonzeroDbetaiDj(inonzero)%elem(ixO^S)
          end do

          do inonzero = 1, myM%nnonzeroDalphaDj
             j = myM%nonzeroDalphaDj(inonzero)%j
             tmp(ixO^S) = tmp(ixO^S) - sU(ixO^S,j)*myM%nonzeroDalphaDj(inonzero)%elem(ixO^S)
          end do

          w(ixO^S,tau_) = w(ixO^S,tau_) + qdt*tmp(ixO^S)

          {#IFDEF GLM
!        b^C
          {^C&
                                ! s[B^C_] = -B**i/alpha dbetaid^C + beta**^C/alpha**2 * B**i * dalphadi
          tmp(ixO^S) = zero
          do i = 1, ^NC
             if (dbetaidj_is_zero(i,^C)) cycle
             tmp(ixO^S) = tmp(ixO^S) - wCT(ixO^S,b0_+i) * myM%dbetaidj(i,^C)%elem(ixO^S)            
          end do
          tmp(ixO^S) = tmp(ixO^S)/myM%alpha(ixO^S)

          if (.not. beta_is_zero(^C)) then
             tmp2(ixO^S) = zero
             do inonzero = 1, myM%nnonzeroDalphaDj
                i = myM%nonzeroDalphaDj(inonzero)%j
                tmp2(ixO^S) = tmp2(ixO^S) + wCT(ixO^S,b0_+i)*myM%nonzeroDalphaDj(inonzero)%elem(ixO^S)
             end do
             tmp(ixO^S) = tmp(ixO^S) + myM%beta(^C)%elem(ixO^S)/myM%alpha(ixO^S)**2 * tmp2(ixO^S)
          end if

          w(ixO^S,b^C_) = w(ixO^S,b^C_) + qdt*tmp(ixO^S)
          \}


!   psi
          ! s[psi_] = - phi DbetaiDi - 1/2 phi gammainv**ij * beta**k * DgijDk

          tmp(ixO^S) = zero
          do i = 1, ^NC
             if (dbetaidj_is_zero(i,i)) cycle
             tmp(ixO^S) = tmp(ixO^S) - myM%dbetaidj(i,i)%elem(ixO^S)
          end do
          tmp(ixO^S) = tmp(ixO^S) * wCT(ixO^S,psi_)

          tmp2(ixO^S) = zero
          do inonzero = 1, myM%nnonzeroDgDk
             i = myM%nonzeroDgDk(inonzero)%i
             j = myM%nonzeroDgDk(inonzero)%j
             k = myM%nonzeroDgDk(inonzero)%k
             if ( beta_is_zero(k) ) cycle

             tmp2(ixO^S) = tmp2(ixO^S) &
                  + myM%beta(k)%elem(ixO^S) * myM%gammainv(i,j)%elem(ixO^S) * myM%DgDk(i,j,k)%elem(ixO^S)

          end do
          tmp(ixO^S) = tmp(ixO^S) - 0.5d0 * wCT(ixO^S,psi_) * tmp2(ixO^S)

          w(ixO^S,psi_) = w(ixO^S,psi_) + qdt*tmp(ixO^S)
          }
  endif
  
   ! Source terms for m1-radiation equations
   ! dkdg and Kij calculated for M1 in non-DY_SP
   {#IFDEF M1
   	 do i= 1, ^NC
       do j=i, ^NC
	      ! this order of indices to match DY_SP case
          {^C& d_k_gamma_ij(ixO^S,i,j,^C) = myM%DgDk(i,j,^C)%elem(ixO^S) \}
	   end do
	   end do
	   
	   ! need K(1,1),K(1,2),K(1,3),K(2,2),K(2,3),K(3,3)
      !> K_ij = 1/(2*alpha) *(gamma_ik*dbetakdj + gamma_jk*dbetakdi + beta^k *dgamma_ij_dk)
      do i=1,^NC
         do j=1,^NC
            m1_tmp1(ixO^S)=0.0d0
            m1_tmp2(ixO^S)=0.0d0
            m1_tmp3(ixO^S)=0.0d0
            do k = 1, ^NC
               m1_tmp1(ixO^S)=m1_tmp1(ixO^S)+myM%g(i,k)%elem(ixO^S)*myM%dbetaidj(k,j)%elem(ixO^S)
               m1_tmp2(ixO^S)=m1_tmp2(ixO^S)+myM%g(j,k)%elem(ixO^S)*myM%dbetaidj(k,i)%elem(ixO^S)
               m1_tmp3(ixO^S)=m1_tmp3(ixO^S)+myM%beta(k)%elem(ixO^S)*myM%DgDk(i,j,k)%elem(ixO^S)
            end do
            Kij_local(ixO^S,i,j) = (0.5d0/myM%alpha(ixO^S))* (m1_tmp1(ixO^S) + m1_tmp2(ixO^S) + m1_tmp3(ixO^S))
         end do
      end do

      call m1_add_geometrical_sources(ixI^L,ixO^L,x,wCT(ixI^S,1:nw),w, wold, qdt,d_k_gamma_ij,Kij_local,qtC, metric_M1)
    } ! end M1
	
! based on these conservative to find new aux
!    call getaux(.true.,w,x,ixI^L,ixO^L,'addgeometry')
  end subroutine addgeometry
  !=============================================================================
  subroutine addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

    ! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

    include 'amrvacdef.f'

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit

    double precision :: dx^D
    !-----------------------------------------------------------------------------

    dx^D=dxlevel(^D);

    ! Sources related to div B

    {#IFDEF GLM
    if(qsourcesplit) then
       call addsource_glmB(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
    else
       call addsource_glmA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)   
    end if
    }

    {#IFDEF ELECTRONS
     call addsource_elecheating(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
    }
    
  end subroutine addsource
  !=============================================================================
  subroutine getcurrent(ixI^L,ixO^L,w,current)

    use mod_metric, only: lower3
    include 'amrvacdef.f'

    integer, intent(in)           :: ixO^L, ixI^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(out) :: current(ixI^S,1:3)

    ! .. local ..
    double precision              :: bvec(ixI^S,1:ndir)
    !-----------------------------------------------------------------------------
    call mpistop('Dont call here, not fix yet: getcurrent')
    call lower3(ixI^L,ixO^L^LADD1,myM,w(ixI^S,b1_:b^NC_),bvec)

    call curl3(bvec,ixI^L,ixO^L,current)

  end subroutine getcurrent

  {#IFDEF ELECTRONS
  !=============================================================================
subroutine addsource_elecheating(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

  ! Add electron heating related sources to w (electron entropy) within ixO

  include 'amrvacdef.f'

  integer, intent(in) :: ixI^L, ixO^L, iw^LIM
  double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
  double precision, intent(in) :: dx^D
  double precision, intent(inout) :: wCT(ixI^S,1:nw)  
  double precision, intent(inout) :: w(ixI^S,1:nw)
  ! .. local ..

  double precision,dimension(ixI^S)          :: kappag, kappagh, kappae, kappaeh, fe, tmp1
  double precision :: tmp2
  !-----------------------------------------------------------------------------
call mpistop("elecheating stop, you need to modify it!")
  
  ! Calculate of advected gas entropy and advected electron entropy  

  kappagh(ixO^S)=w(ixO^S,Ds_)/w(ixO^S,d_)
  kappaeh(ixO^S)=w(ixO^S,Dse_)/w(ixO^S,d_)

  ! Transform w to primitive variables
  call primitive(ixI^L,ixO^L,w,x)

  ! Calculate of gas entropy (including dissipative heating)
  kappag(ixO^S)=w(ixO^S,s_)

  ! Transform wCT to primitive variables (need density at n+1/2)
  call primitive(ixI^L,ixO^L,wCT,x) 

  ! Calculate electron heating rate
  call calheatingrate(wCT,x,ixI^L,ixO^L,fe)
  
  ! Calculate electron entropy by adding electron heating
  ! Eq. (27) in Ressler et al. (2015)
  
  tmp2 = (eqpar(gammae_)-1.d0)/(eqpar(gamma_)-1.d0)

  tmp1(ixO^S) = wCT(ixO^S,rho_)**(eqpar(gamma_)-eqpar(gammae_))
  
  w(ixO^S,se_) = kappaeh(ixO^S) + tmp1(ixO^S)*tmp2*fe(ixO^S)*(kappag(ixO^S)-kappagh(ixO^S)) *qdt/dt

!  where(w(ixO^S,se_) .gt. kappag(ixO^S))
!    w(ixO^S,se_) = 0.99d0*kappag(ixO^S)
!  endwhere   
  
  ! Transform w back to conserved variables
  call conserve(ixI^L,ixO^L,wCT,x,patchfalse)
  call conserve(ixI^L,ixO^L,w,x,patchfalse)

end subroutine addsource_elecheating
  !=============================================================================
subroutine calheatingrate(wCT,x,ixI^L,ixO^L,fe)

    use mod_physaux, only: get_b2
    include 'amrvacdef.f'
  
    integer, intent(in)               :: ixI^L,ixO^L
    double precision, intent(inout)   :: wCT(ixI^S,1:nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    double precision, intent(out)     :: fe(ixI^S)
    ! .. local ..

    character*31                        :: typeelheat
    double precision,dimension(ixI^S)   :: Qpe, rhogas, Pe, Pgas, Pp, Te, Tgas, Tp, b2, beta, c1, c2, c3, pow1, tmp1, tmp2, tmp3, ent, sigmaw, betamax
    double precision, parameter         :: small=1.d-16, mpme=1836.15267389d0
    !-----------------------------------------------------------------------------
!  input prim
   call mpistop("dont touch calheatingrate first")


   typeelheat = 'kawazura'

    select case (typeelheat)

    case('constant')
    
      fe(ixO^S) = 0.5d0
!
!      
    case('turbulent')
      ! heating ratio of tuebulent model (Howes (2011))
 
      ! Calculation of pressures
      rhogas(ixO^S) = wCT(ixO^S,rho_)
      
      Pe(ixO^S) = wCT(ixO^S,se_) * rhogas(ixO^S)**eqpar(gammae_)
      Pgas(ixO^S) = wCT(ixO^S,pp_)

      ! Pe > Pgas => keep Pe < Pgas
      !where( Pe(ixO^S) .gt. Pgas(ixO^S) )
      !   print*, 'Pe > pgas'
      !   Pe(ixO^S) = 0.999d0*Pgas(ixO^S)  
      !endwhere
    
      ! Calculation of temperatures
      Tgas(ixO^S) = Pgas(ixO^S)/rhogas(ixO^S)
      Te(ixO^S) = Pe(ixO^S) /rhogas(ixO^S)
      Tp(ixO^S) = ((eqpar(gammap_)-1.d0)/(eqpar(gamma_)-1.d0))*Tgas(ixO^S) - ((eqpar(gammap_)-1.d0)/(eqpar(gammae_)-1.d0))*Te(ixO^S)
      where(Tp(ixO^S) .le. 0.d0)
        Tp(ixO^S) = 0.01d0*Te(ixO^S)
      endwhere

      ! calculation of plasma beta

      call get_b2(ixI^L,ixO^L,wCT(ixI^S,1:nw),x,b2,w_is_primitive=.true.)
    
      beta(ixO^S) =  2.d0*wCT(ixO^S,pp_)/(b2(ixO^S)+small)

      ! calculation of Qp/Qe
      ! Eq. (48) in Ressler et al. (2015)
      where (Tp(ixO^S) .ge. Te(ixO^S))
         c1(ixO^S) = 0.92d0
         c2(ixO^S) = 1.6d0*(Te(ixO^S)/Tp(ixO^S))
         c3(ixO^S) = 18.d0 + 5.d0*log10(Tp(ixO^S)/Te(ixO^S))
      elsewhere(Tp(ixO^S) .lt. Te(ixO^S))   
         c1(ixO^S) = 0.92d0
         c2(ixO^S) = 1.2d0*(Te(ixO^S)/Tp(ixO^S))
         c3(ixO^S) = 18.d0
      endwhere

      pow1(ixO^S) = 2.d0-0.2d0*log10(Tp(ixO^S)/Te(ixO^S))
      tmp1(ixO^S) = c2(ixO^S)**2 + beta(ixO^S)**pow1(ixO^S)
      tmp2(ixO^S) = c3(ixO^S)**2 + beta(ixO^S)**pow1(ixO^S)
      tmp3(ixO^S) = sqrt(mpme*(Tp(ixO^S)/Te(ixO^S)))
      Qpe(ixO^S)  = c1(ixO^S)*(tmp1(ixO^S)/tmp2(ixO^S))*tmp3(ixO^S)*exp(-1.d0/beta(ixO^S))

      fe(ixO^S) = 1.d0/(1.d0+Qpe(ixO^S))
!
!      
    case('kawazura')
      ! heating ratio of tuebulent model (Kawazura et al. (2018))
 
      ! Calculation of pressures
      rhogas(ixO^S) = wCT(ixO^S,rho_)
      
      Pe(ixO^S) = wCT(ixO^S,se_) * rhogas(ixO^S)**eqpar(gammae_)
      Pgas(ixO^S) = wCT(ixO^S,pp_)

      ! Pe > Pgas => keep Pe < Pgas
      !where( Pe(ixO^S) .gt. Pgas(ixO^S) )
      !   print*, 'Pe > pgas'
      !   Pe(ixO^S) = 0.999d0*Pgas(ixO^S)  
      !endwhere
    
      ! Calculation of temperatures
      Tgas(ixO^S) = Pgas(ixO^S)/rhogas(ixO^S)
      Te(ixO^S) = Pe(ixO^S) /rhogas(ixO^S)
      Tp(ixO^S) = ((eqpar(gammap_)-1.d0)/(eqpar(gamma_)-1.d0))*Tgas(ixO^S) - ((eqpar(gammap_)-1.d0)/(eqpar(gammae_)-1.d0))*Te(ixO^S)
      where(Tp(ixO^S) .le. 0.d0)
        Tp(ixO^S) = 0.01d0*Te(ixO^S)
      endwhere

      ! calculation of plasma beta
      call get_b2(ixI^L,ixO^L,wCT(ixI^S,1:nw),x,b2,w_is_primitive=.true.)
    
      beta(ixO^S) =  2.d0*wCT(ixO^S,pp_)/(b2(ixO^S)+small)

      ! calculation of Qp/Qe
      ! Eq (2) in Kawazura et al. (2018)
      tmp1(ixO^S) = (beta(ixO^S)/15.d0)**(-1.4d0)
      tmp2(ixO^S) = exp(-0.1d0*(Te(ixO^S)/Tp(ixO^S)))
      Qpe(ixO^S) = 35.d0/(1.d0+tmp1(ixO^S)*tmp2(ixO^S))
      fe(ixO^S) = 1.d0/(1.d0+Qpe(ixO^S))
!
!      
    case('magrec')
      ! heating ratio of magnetic reconnection model (Rowan et al. (2017))
      
      ! calculation of plasma beta
      call get_b2(ixI^L,ixO^L,wCT(ixI^S,1:nw),x,b2,w_is_primitive=.true.)
    
      beta(ixO^S) =  2.d0*wCT(ixO^S,pp_)/(b2(ixO^S)+small)

      ent(ixO^S) = wCT(ixO^S,rho_)+(eqpar(gamma_)/(eqpar(gamma_)-1.d0))*wCT(ixO^S,pp_)
      sigmaw(ixO^S) = b2(ixO^S)/ent(ixO^S)
      betamax(ixO^S) = 1.d0/(4.d0*sigmaw(ixO^S))

      ! calculation of Qp/Qe
      ! Eq (13) in Chael et al. (2018)
      tmp1(ixO^S) = -(1.d0-beta(ixO^S)/betamax(ixO^S))
      Qpe(ixO^S) = 0.5d0*exp(-tmp1(ixO^S)/(0.8d0+sqrt(sigmaw(ixO^S))))

      fe(ixO^S) = 1.d0/(1.d0+Qpe(ixO^S))


    
    case default
      call mpistop('Unknown typeelheating in calheatingrate')
    end select   
!
  end subroutine calheatingrate
}
  !=============================================================================
  ! just dummies:
  !=============================================================================
  subroutine e_to_rhos(ixI^L,ixO^L,w,x)

    include 'amrvacdef.f'
    integer:: ixI^L,ixO^L
    double precision:: w(ixI^S,nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    !-------------------------------------------------------------------------

    call mpistop("e to rhos unavailable")

  end subroutine e_to_rhos
  !=============================================================================
  subroutine rhos_to_e(ixI^L,ixO^L,w,x)

    include 'amrvacdef.f'

    integer:: ixI^L,ixO^L
    double precision:: w(ixI^S,nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    !-----------------------------------------------------------------------------

    call mpistop("e to rhos unavailable")

  end subroutine rhos_to_e

  !< This function performs all the required implicit 
  !< updates globally (on all grids). This means
  !<    psa = psb + dtfactor * qdt * F_im(psa)
  !< Note that psa is already set to psb, so nothing needs
  !< to be done for variables where there are no implicit 
  !< sources. This routine is called in advect (advance.t)
  !< when an implicit evolution scheme is selected 
  !< (e.g. IMEX222) 

  subroutine implicit_update(dtfactor,qdt,qtC,psa,psb, psm1impl)
  !subroutine implicit_update(dtfactor,igrid,qdt,qtC,sa,sb, sm1impl, ixI^L)
      {#IFNDEF M1_EXPLICIT
      {#IFDEF M1 
      use mod_m1, only: m1_add_collisional_sources
      }
      }
      include 'amrvacdef.f' 

      type(state), dimension(ngridshi)    :: psa !< Implicit update computed here 
      type(state), dimension(ngridshi)    :: psb !< Left unchanged 
      type(m1_impl_state), dimension(ngridshi)  :: psm1impl

      !type(state)      :: sa !< Implicit update computed here 
      !type(state)      :: sb !< Left unchanged 
      !type(m1_impl_state)   :: sm1impl
      
      double precision, intent(in) :: dtfactor   !< Timestep factor 
      double precision, intent(in) :: qdt        !< Timestep 
      double precision, intent(in) :: qtC        !< Current time  
      !integer, intent(in) :: igrid, ixI^L      
      integer :: iigrid, igrid 
      !internal
      !double precision, dimension(ixI^L, 1:nw) :: wold, wcons
      !double precision, dimension(ixI^L, 1:nm1rad_eas) :: wradimpl_r

      {#IFNDEF M1_EXPLICIT
      ! grid loop 
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        {#IFDEF M1 
        !write(*,*) "implicit update, igrid",igrid
        ! dummy routine for now psa%w%w for centered variables
        !call m1_add_collisional_sources(dtfactor,qdt,qtC,px(igrid)%x(ixG^T,1:ndim),sa%w%w(ixG^T,1:nw),sb%w%w(ixG^T,1:nw),sm1impl%pwrad%wradimpl(ixG^T,1:nm1rad_eas),ixG^LL)
        !call m1_add_collisional_sources(dtfactor,qdt,qtC,px(igrid)%x(ixG^T,1:ndim),sa%w%w(ixG^T,1:nw),sb%w%w(ixG^T,1:nw),sm1impl%pwrad%wradimpl(ixG^T,1:nm1rad_eas),ixI^L)
        call m1_add_collisional_sources(dtfactor,qdt,qtC,px(igrid)%x(ixG^T,1:ndim),psa(igrid)%w%w(ixG^T,1:nw),psb(igrid)%w%w(ixG^T,1:nw),psm1impl(igrid)%pwrad%wradimpl(ixG^T,1:nm1rad_eas),ixG^LL)
        !call m1_add_collisional_sources(dtfactor,qdt,qtC,psb%x%x,psa%w%w,psb%w%w,psm1impl%pwrad%wradimpl,ixG^LL)
        }
      end do
      }

  end subroutine implicit_update

  !=============================================================================
  ! end module amrvacphys
  !=============================================================================

