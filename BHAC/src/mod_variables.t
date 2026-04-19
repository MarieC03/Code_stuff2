module mod_variables

  implicit none
  public

  logical            :: var_reconstruct(1000)
  logical            :: var_fillbc(1000)

  ! types of variables
  integer, parameter :: cons_var           = 1
  integer, parameter :: prim_var           = 2
  integer, parameter :: tracer_var         = 3
  integer, parameter :: aux_var            = 4
  integer, parameter :: extra_var          = 5
  integer, parameter :: gwbr_var           = 6
  integer, parameter :: radiation_var      = 7 
  integer, parameter :: radiation_impl_var = 8 

  ! DY_SP:    D_,  s^C_, tau_, b1_, b2_, b3_, T_eps_, rho_, u^C_, metric
  ! nonDY_SP: D_,  s^C_, tau_, b1_, b2_, b3_, pp_, T_eps_, rho_, u^C_, lfac_, xi_, mygeo (myM)
  ! define TABEOS, have two extra: Dye_, ye_
  ! define GW_BR, have three extra: U_br_, U_br1_, U_br2_, U_br3_, R_br_, U5_br_, h_00_
  ! define M1, have 5species extra: nrad1_, erad1_, frad11_, frad12_, frad13_ ...

contains

  subroutine setup_variables
    include 'amrvacdef.f'

!still need to follow order
      nw = 0

!     Cons
      D_     = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.)  
{#IFDEF TABEOS
      Dye_   = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.)  
}
{^C&  s^C_   = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true., &
                            vector=.true., ix=^C)  \}
      s0_    = s1_ - 1
      tau_   = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      e_     = tau_
      rhos_  = tau_

    !  foo_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 

! staggered CT, no need to reconstruct B^C before passing CT routines, except ppm limiter
{#IFDEF STAGGERED
{^C&  b^C_   = var_set_wvar(cons_var, need_rec=.true., fill_bc=.true.,&
                            vector=.true., ix=^C)   \}
}
{#IFNDEF STAGGERED
{^C&  b^C_   = var_set_wvar(cons_var, need_rec=.true., fill_bc=.true.,&
                            vector=.true., ix=^C)   \}
}
      b0_    = b1_ - 1

      {#IFDEF FOO
      erad1_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      }
      {#IFDEF FOO4
      foo_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo1_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo2_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo3_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      }
      {#IFDEF FOO5
      foo_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo1_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo2_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo3_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo4_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      }

      {#IFDEF FOO5_RAD
      foo_ = var_set_wvar(radiation_var, need_rec=.false., fill_bc=.true.) 
      foo1_ = var_set_wvar(radiation_var, need_rec=.false., fill_bc=.true.) 
      foo2_ = var_set_wvar(radiation_var, need_rec=.false., fill_bc=.true.) 
      foo3_ = var_set_wvar(radiation_var, need_rec=.false., fill_bc=.true.) 
      foo4_ = var_set_wvar(radiation_var, need_rec=.false., fill_bc=.true.) 
      }

      {#IFDEF FOO6
      foo_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo1_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo2_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo3_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo4_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo5_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      }

      {#IFDEF FOO10
      foo_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo1_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo2_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo3_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo4_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo5_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo6_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo7_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo8_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      foo9_ = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
      }


{#IFDEF ENTROPY
      Ds_    = var_set_wvar(cons_var, need_rec=.false., fill_bc=.true.) 
}
      
{#IFDEF TRACER
{^FL& Dtr^FL_ = var_set_wvar(cons_var, need_rec=.true., fill_bc=.true.)  \}
}

!--- m1 has to be here because is ncons too
! rad 
{#IFDEF M1
{^KSP&
   nrad^KSP_ = var_set_wvar(radiation_var,need_rec=.true., fill_bc=.true.)
   erad^KSP_ = var_set_wvar(radiation_var,need_rec=.true., fill_bc=.true.)
   {^C& frad^KSP^C_ = var_set_wvar(radiation_var,need_rec=.true., fill_bc=.true.&
                                 , vector=.true., ix=^C) \}

    ! ----- rates variables: ----
     kappa_a^KSP_  = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.) 
     kappa_s^KSP_  = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.)  
     Q_er^KSP_     = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.)  
     !eta^KSP_     = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.)  
     !eta_n^KSP_   = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.) 
     Q_nr^KSP_     = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.)  
     kappa_nr^KSP_ = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.)  
     tau_path^KSP_ = var_set_wvar(radiation_impl_var,need_rec=.false., fill_bc=.false.)  
  \} 
} 

!     prim
{#IFNDEF DY_SP
!  save mem
  {#IFDEF RECONSTRUCT_PRESS
        pp_    = var_set_wvar(prim_var, need_rec=.true., fill_bc=.true.) 
  }
  {#IFNDEF RECONSTRUCT_PRESS
        pp_    = var_set_wvar(prim_var, need_rec=.false., fill_bc=.true.) 
  }
}
      T_eps_ = var_set_wvar(prim_var, need_rec=.true., fill_bc=.true.) 
{#IFDEF TABEOS
      ye_    = var_set_wvar(prim_var, need_rec=.true., fill_bc=.true.) 
}
      rho_   = var_set_wvar(prim_var, need_rec=.true., fill_bc=.true.) 
{^C&  u^C_   = var_set_wvar(prim_var, need_rec=.true., fill_bc=.true.,&
                            vector=.true., ix=^C)   \}
      u0_    = u1_ - 1
{^C&  v^C_   = u^C_  \}       
      v0_    = v1_ - 1
{#IFDEF ENTROPY
      s_    = var_set_wvar(prim_var, need_rec=.true., fill_bc=.true.) 
}
   

{#IFDEF GLM
   ! later on
}


{#IFNDEF DY_SP
! Do not allocate lfac, xi for saving memory
      lfac_   = var_set_wvar(aux_var, need_rec=.false., fill_bc=.true.) 
      ! xi = W**2 rho * h
      xi_     = var_set_wvar(aux_var, need_rec=.false., fill_bc=.true.) 
}
{#IFDEF C2PREL
      c2p_d_  = var_set_wvar(aux_var, need_rec=.false., fill_bc=.true.) 
{^C&  c2p_s^C_= var_set_wvar(aux_var, need_rec=.false., fill_bc=.true.)  \}
      c2p_tau_= var_set_wvar(aux_var, need_rec=.false., fill_bc=.true.) 
} 

! fixme: do we need to set beta, Xvec to be vector type?
{#IFDEF DY_SP
          alp_metric_    = var_set_wvar(aux_var, need_rec=.false., fill_bc=.true., is_metric=.true.) 
    {^C&  beta_metric^C_ = var_set_wvar(aux_var, need_rec=.false., fill_bc=.true., is_metric=.true.) \}
    {^C&  Xvec_metric^C_ = var_set_wvar(aux_var, need_rec=.false., fill_bc=.true., is_metric=.true.) \}
          psi_metric_    = var_set_wvar(aux_var, need_rec=.false., fill_bc=.true., is_metric=.true.)
}

{#IFDEF GW_BR
! do not casify to be metric, skip interpolation for the interface
          U_br_     = var_set_wvar(gwbr_var, need_rec=.false., fill_bc=.true., is_metric=.true.)
    ! U_bri is very neglectable --> save cost and memory, no need them
          ! Need this to be metric variable as well
          h_00_     = var_set_wvar(gwbr_var, need_rec=.false., fill_bc=.true., is_metric=.true.)
          {#IFDEF DELTA         
          delta_br_ = var_set_wvar(gwbr_var, need_rec=.false., fill_bc=.true., is_metric=.true.)
          }
}

          

{#IFDEF TRACER
{^FL&  tr^FL_ = Dtr^FL_  \}
}
  se_= Dse_
!  s_ = Ds_
  epsinf_=Depsinf_
  rho0_=Drho0_
  rho1_=Drho1_
  n0_=Dn0_
  n_=Dn_
 
 ! polar variables names
  sr_=s0_+r_
  sphi_=s0_+phi_
  sz_=s0_+z_
  vr_=v0_+r_
  vphi_=v0_+phi_
  vz_=v0_+z_
  uz_=v0_+z_
  ur_=v0_+r_
  uphi_=v0_+phi_
  br_=b0_+r_
  bphi_=b0_+phi_
  bz_=b0_+z_
  ee_=e_


  if (mype == 0) then

    write(*,*) '------------The number of the bounds of indices--------------'

    print*,'nw          :', nw
    print*,'nwflux      :', nwflux
    print*,'ncons       :', ncons

    print*,'nprim       :', nprim
    print*,'nwfluxtr    :', nwfluxtr
    print*,'nwaux       :', nwaux
    print*,'nwextra     :', nwextra
    print*,'nmetric     :', nmetric
    print*,'nwgwbr      :', nwgwbr  
    print*,'nm1rad      :', nm1rad
    print*,'nm1rad_eas  :', nm1rad_eas
    print*,'nwflux_lo   :', nwflux_lo
    print*,'nwflux_hi   :', nwflux_hi
    print*,'ncons_lo    :', ncons_lo
    print*,'ncons_hi    :', ncons_hi
    print*,'nprim_lo    :', nprim_lo
    print*,'nprim_hi    :', nprim_hi
    print*,'nwfluxtr_lo :', nwfluxtr_lo
    print*,'nwfluxtr_hi :', nwfluxtr_hi
    print*,'nwaux_lo    :', nwaux_lo
    print*,'nwaux_hi    :', nwaux_hi
    print*,'nwextra_lo  :', nwextra_lo
    print*,'nwextra_hi  :', nwextra_hi
    print*,'nmetric_lo  :', nmetric_lo
    print*,'nmetric_hi  :', nmetric_hi
    print*,'nwgwbr_lo   :', nwgwbr_lo
    print*,'nwgwbr_hi   :', nwgwbr_hi
    print*,'nm1rad_lo   :', nm1rad_lo
    print*,'nm1rad_hi   :', nm1rad_hi
    print*,'nm1rad_eas_lo   :', nm1rad_eas_lo
    print*,'nm1rad_eas_hi   :', nm1rad_eas_hi


!    write(*,*) '--------------The indices of Variables------------------'
!
!    print*,'d_            :', d_           , 'reconstruct?', var_reconstruct(d_)
!    print*,'dye_          :', dye_         , 'reconstruct?', var_reconstruct(dye_)
!    print*,'s1_           :', s1_          , 'reconstruct?', var_reconstruct(s1_) 
!    print*,'s2_           :', s2_          , 'reconstruct?', var_reconstruct(s2_)
!    print*,'s3_           :', s3_          , 'reconstruct?', var_reconstruct(s3_)
!    print*,'tau_          :', tau_         , 'reconstruct?', var_reconstruct(tau_)
!    print*,'Ds_           :', Ds_          , 'reconstruct?', var_reconstruct(Ds_)
!    print*,'Dtr1_         :', Dtr1_        , 'reconstruct?', var_reconstruct(Dtr1_)
!    print*,'pp_           :', pp_          , 'reconstruct?', var_reconstruct(pp_)
!    print*,'ye_           :', ye_          , 'reconstruct?', var_reconstruct(ye_)
!    print*,'rho_          :', rho_         , 'reconstruct?', var_reconstruct(rho_)
!    print*,'u1_           :', u1_          , 'reconstruct?', var_reconstruct(u1_)
!    print*,'u2_           :', u2_          , 'reconstruct?', var_reconstruct(u2_)
!    print*,'u3_           :', u3_          , 'reconstruct?', var_reconstruct(u3_)
!    print*,'b1_           :', b1_          , 'reconstruct?', var_reconstruct(b1_)
!    print*,'b2_           :', b2_          , 'reconstruct?', var_reconstruct(b2_)
!    print*,'b3_           :', b3_          , 'reconstruct?', var_reconstruct(b3_)
!    print*,'T_eps_        :', T_eps_       , 'reconstruct?', var_reconstruct(T_eps_)
!    print*,'cs2_          :', cs2_         , 'reconstruct?', var_reconstruct(cs2_)
!    print*,'tr1_          :', tr1_         , 'reconstruct?', var_reconstruct(tr1_)
!    print*,'lfac_         :', lfac_        , 'reconstruct?', var_reconstruct(_lfac)
!    print*,'xi_           :', xi_          , 'reconstruct?', var_reconstruct(xi_)
!    print*,'alp_metric_   :', alp_metric_  , 'reconstruct?', var_reconstruct(alp_metric_)
!    print*,'beta_metric1_ :', beta_metric1_, 'reconstruct?', var_reconstruct(beta_metric1_)
!    print*,'beta_metric2_ :', beta_metric2_, 'reconstruct?', var_reconstruct(beta_metric2_)
!    print*,'beta_metric3_ :', beta_metric3_, 'reconstruct?', var_reconstruct(beta_metric3_)
!    print*,'Xvec_metric1_ :', Xvec_metric1_, 'reconstruct?', var_reconstruct(Xvec_metric1_)
!    print*,'Xvec_metric2_ :', Xvec_metric2_, 'reconstruct?', var_reconstruct(Xvec_metric2_)
!    print*,'Xvec_metric3_ :', Xvec_metric3_, 'reconstruct?', var_reconstruct(Xvec_metric3_)
!    print*,'psi_metric_   :', psi_metric_  , 'reconstruct?', var_reconstruct(psi_metric_)
!    print*,'U_br_         :', U_br_        , 'reconstruct?', var_reconstruct(U_br_)
!    print*,'U_br1_        :', U_br1_       , 'reconstruct?', var_reconstruct(U_br1_)
!    print*,'U_br2_        :', U_br2_       , 'reconstruct?', var_reconstruct(U_br2_)
!    print*,'U_br3_        :', U_br3_       , 'reconstruct?', var_reconstruct(U_br3_)
!    print*,'h_00_         :', h_00_        , 'reconstruct?', var_reconstruct(h_00_)
!    print*,'delta_br_     :', delta_br_    , 'reconstruct?', var_reconstruct(delta_br_)
!    {^KSP&
!    print*,'nrad',^KSP,'_   :',nrad^KSP_    , 'reconstruct?', var_reconstruct(nrad^KSP_)
!    print*,'erad',^KSP,'_   :',erad^KSP_    , 'reconstruct?', var_reconstruct(erad^KSP_)
!    print*,'frad',^KSP,'1_  :',frad^KSP1_    , 'reconstruct?', var_reconstruct(frad^KSP1_)
!    print*,'frad',^KSP,'2_  :',frad^KSP2_    , 'reconstruct?', var_reconstruct(frad^KSP2_)
!    print*,'frad',^KSP,'3_  :',frad^KSP3_    , 'reconstruct?', var_reconstruct(frad^KSP3_)
!    \}    
{#IFDEF STAGGERED
    {^D& print*,'bs^D_          :', bs^D_ \}
}

    write(*,*)' iw_vector(1:2), s0_, b0_'
    write(*,*) iw_vector(1:2), s0_, b0_
    write(*,*) nvector, 'vector'
  endif
 
  call allocation_amrvacdef
  end subroutine setup_variables

  subroutine allocation_amrvacdef
    include 'amrvacdef.f'
     nflag_ = nw+1
     allocate(flags(nflag_))
     allocate(wflags(nflag_))

     allocate(iws(nwflux+nwaux))
{#IFDEF STAGGERED
     iws = [(itmp-b0_, itmp = 1, nwflux+nwaux)]
}
{#IFNDEF STAGGERED
     iws = [(0, itmp = 1, nwflux+nwaux)]
}

     ! from amrvacdef
     allocate(typeentropy(1:nw))
     allocate(entropycoef(1:nw))
     allocate(loglimit(1:nw))
     allocate(logflag(1:nw))
     allocate(writew(1:nw))
     allocate(typeB(nw,nhiB))
     allocate(normvar(0:nw))

!   write(*,*) shape(writew), 'shape of writew'

  end subroutine allocation_amrvacdef

   
  !> Set w variable, separate their types in to prim, flux, metric
  !> vector and ix are only for conservative vector
  function var_set_wvar(var_type, need_rec, fill_bc, vector, ix, is_metric) result(iw)
    include 'amrvacdef.f'
    integer, intent(in) :: var_type  !< variable type
    logical, intent(in) :: need_rec  !< Require reconstruction for hydro only (default: true)
    logical, optional, intent(in) :: vector    ! vector or not
    logical, optional, intent(in) :: is_metric ! metric vars?
    logical, intent(in)           :: fill_bc   ! need to fill bc?
    integer, intent(in), optional :: ix        !< Optional index for vector

    integer                       :: iw, var_t
    logical                       :: flag

    ! total number of w variables
    nw  = nw + 1
    ! the current number of the variable inside w
    iw  = nw

	
    if (present(vector)) then
      if (vector) then 
       if (present(ix)) then
         if (ix==1) then  ! only for the first index of vector
           nvector = nvector + 1
           iw_vector(nvector) = iw - 1
         endif
       endif
      end if
    endif

    var_reconstruct(iw) = need_rec

    var_fillbc(iw) = fill_bc

    var_t = -1
    var_t = var_type

    select case (var_t)
    case (cons_var)
       if ( ncons_lo == 0 ) ncons_lo = iw
       ncons = ncons + 1
       ncons_hi = iw

       ! as cons is also nwflux vars
       if ( nwflux_lo == 0 ) nwflux_lo = iw
       nwflux = nwflux + 1
       nwflux_hi = iw
    case (prim_var)
       if ( nprim_lo == 0 ) nprim_lo = iw
       nprim = nprim + 1
       nprim_hi = iw
    
       ! as prim is also nwflux vars 
       if ( nwflux_lo == 0 ) call mpistop('prim vars cannot be the first nwflux')
       nwflux = nwflux + 1
       nwflux_hi = iw

    case (tracer_var)
       if ( nwfluxtr_lo == 0 ) nwfluxtr_lo = iw
       nwfluxtr = nwfluxtr + 1
       nwfluxtr_hi = iw

       ! as nwfluxtr is also nwflux vars 
       if ( nwflux_lo == 0 ) call mpistop('tracer vars cannot be the first nwflux')
       nwflux = nwflux + 1
       nwflux_hi = iw

    case (aux_var, gwbr_var)
       if ( nwaux_lo == 0 ) nwaux_lo = iw
       nwaux = nwaux + 1
       nwaux_hi = iw


    case (extra_var)
       if ( nwextra_lo == 0 ) nwextra_lo = iw
       nwextra = nwextra + 1
       nwextra_hi = iw

    !case (gwbr_var)
    !   if ( nwgwbr_lo == 0 ) nwgwbr_lo = iw
    !   nwgwbr    = nwgwbr + 1
    !   nwgwbr_hi = iw

    case (radiation_var) 
       if ( ncons_lo == 0 ) stop 'Radiation vars cannot be the first nwcons'
       ncons = ncons + 1
       ncons_hi = iw
       if ( nwflux_lo == 0 ) stop 'Radiation vars cannot be the first nwflux'
       nwflux = nwflux + 1
       nwflux_hi = iw
       if ( nm1rad_lo == 0 ) nm1rad_lo = iw
       nm1rad = nm1rad + 1
       nm1rad_hi = iw

     case (radiation_impl_var) 
          {#IFDEF M1
            ! remove this var-call from nw as:
            nw = nw - 1
            ! add to implicit radiation variables:
            nm1rad_eas = nm1rad_eas + 1
            if ( nm1rad_eas_lo == 0 .and. nm1rad_eas == 1) nm1rad_eas_lo = 1
            iw = nm1rad_eas
            nm1rad_eas_hi = iw          
          }
    case default
        call mpistop('cannot find this variable type')
    end select

    ! extra type for gwbr
    if (var_t == gwbr_var) then
       if ( nwgwbr_lo == 0 ) nwgwbr_lo = iw
       nwgwbr    = nwgwbr + 1
       nwgwbr_hi = iw
    endif
	
    if (present(is_metric)) then
      if (is_metric) then
        if ( nmetric_lo == 0 ) nmetric_lo = iw
        nmetric = nmetric + 1
        nmetric_hi = iw
      endif  
    endif


  end function var_set_wvar

end module mod_variables
