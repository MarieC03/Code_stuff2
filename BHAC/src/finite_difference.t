!============================================================================
! sCT current w state as input for calculating snew;  snew is the output new w state
subroutine finite_difference(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,sCT,qt,snew,sold,wprim,fC,fE,dx^D,x,sm1)
use mod_interpolate
use mod_imhd_intermediate
{#IFDEF M1
use mod_m1, only: m1_correct_asymptotic_fluxes
}
use mod_m1_metric_interface ! We now need this even when no M1

use mod_limiter
include 'amrvacdef.f'

character(len=*), intent(in)                           :: method
double precision, intent(in)                           :: qdt, qtC, qt, dx^D
integer, intent(in)                                    :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in)  :: x
type(state)                                            :: sCT, snew, sold
type(m1_impl_state)                                    :: sm1
double precision, dimension(ixI^S,1:nw), intent(in)    :: wprim
double precision, dimension(ixI^S,1:nw)                :: metric_i 
double precision, dimension(ixI^S,1:nwflux,1:ndim)     :: fC
double precision, dimension(ixI^S,1:nwflux,1:ndim)     :: fC_local
double precision, dimension(ixI^S,1:^NC)               :: fE

double precision, dimension(ixI^S,1:ndim)              :: xi
double precision, dimension(ixI^S,1:nw)                :: wLC, wRC, wLp, wRp
double precision, dimension(ixI^S,1:nwflux)            :: fLC, fRC, fCT
double precision, dimension(ixI^S)                     :: fadd, fdiss
double precision, dimension(ixI^S)                     :: vLC, vRC, vCT
double precision, dimension(ixI^S,1:ncons)             :: cmaxC, cminC

double precision, dimension(ixI^S)                     :: psi6_surfC^D ! cell interface psi6
double precision, dimension(ixI^S)                     :: psi6         ! cell center psi6
{#IFDEF STAGGERED
! UCT related storage
double precision, dimension(ixI^S,1:ndir,2)            :: vbarC
double precision, dimension(ixI^S,1:ndir,2)            :: vbarLC,vbarRC
double precision, dimension(ixI^S,ndim)                :: cbarmin,cbarmax
integer                                                :: idimE, idimN
}
double precision                                       :: dxinv(1:ndim),dxdim(1:ndim)
integer                                                :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ixtest^L, ix^D
integer                                                :: ixIs^L
integer, dimension(ixI^S)                              :: patchf
logical                                                :: transport(1:nwflux), logiB
{#IFDEF M1
double precision, dimension(ixI^S,1:nm1rad_eas)            :: wradimpl
}
type(m1_metric_helper)   :: metricM1, metricM1_L, metricM1_R    ! MetricM1_L is equal to MetricM1
{#IFDEF UNITS_TESTS
double precision :: TESTqtC = 1.0d+10
}

! variables that related to positivity preserving limiter
! low order flux
double precision, dimension(:^D&,:,:), allocatable    :: fC_low
!double precision, parameter                           :: epsD = smalldouble, epstau = smalldouble
!double precision, dimension(ixI^S,1:nwflux) :: fLC_pp, fRC_pp

!-----------------------------------------------------------------------------
associate(wCT=>sCT%w%w, wnew=>snew%w%w, wold=>sold%w%w)

{#IFDEF M1
  wradimpl=sm1%pwrad%wradimpl !associate(wradimpl=>sm1%pwrad%wradimpl)
}

{#IFDEF DY_SP
  !psi6(ixI^S)         = 1.0d0
  psi6(ixI^S)         = wprim(ixI^S, psi_metric_)**6
  mygeo%dvolume(ixI^S)     = mygeo%dvolume(ixI^S) * psi6(ixI^S)
}

{#IFNDEF DY_SP 
  psi6 = 1.0d0
  {^D& psi6_surfC^D = 1.0d0 \} 
}

  
{#IFDEF STAGGERED
ixIsmin^D=sCT%ws%ixGmin^D;
ixIsmax^D=sCT%ws%ixGmax^D;

staggered : associate(wCTs=>sCT%ws%w, wnews=>snew%ws%w, wolds=>sold%ws%w)
}

logiB=(BnormLF.and.b0_>0)

if (PP_limiter) then
   allocate(fC_low(ixI^S,1:ncons,1:ndim))
   fC_low=0.0d0
end if


if (covariant .and. typelimited .ne. 'predictor') call mpistop('FD: Covariant formulation only implemented with typelimited=predictor')


!if (idimmax>idimmin .and. typelimited=='original' .and. &
!   method/='tvdlf1')&
!   call mpistop("Error in finite difference: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in finite difference: Nonconforming input limits")


^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;

call metricM1%fill_metric(wprim,x,ixI^L,ixI^L)

do idims= idim^LIM
   
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   {#IFNDEF DY_SP
   if (covariant) then
      select case (idims)
      {case (^D)
         myM => mygeo%mSurface^D\}
      end select
   end if
   }

   ! ================================
   ! ====== Assemble indices: =======  
   ! ================================
   hxO^L=ixO^L-kr(idims,^D);
   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
   kxR^L=kxC^L+kr(idims,^D);

   select case (typeemf)
   case default
      ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
      ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
   case ('average')
      ! Flux-interpolated constrained transport needs one more layer in the tangential components:
      ixCmax^D=ixOmax^D+1-kr(idims,^D); ixCmin^D=hxOmin^D-1+kr(idims,^D);
   case ('uct1','uct2','uct2+av')
      ! Upwind constrained transport needs dixB more layers in the tangential components
      ixCmax^D=ixOmax^D+dixB-dixB*kr(idims,^D); ixCmin^D=hxOmin^D-dixB+dixB*kr(idims,^D);
   end select

   ! code test, use del_zanna flux, need one more extra cell, try to extend one more
   if (use_del_zanna) then
      {^D& ixCmax^D = ixCmax^D + 1 \}; {^D& ixCmin^D = ixCmin^D - 1 \};
   endif

   !ixCmin^D=hxOmin^D-(dixB)*(1-kr(idims,^D));
   !ixCmax^D=ixOmax^D+(dixB)*(1-kr(idims,^D));

   !ixCmin^D=hxOmin^D-(dixB-phys_extra_ghostcells)*(1-kr(idims,^D));
   !ixCmax^D=ixOmax^D+(dixB-phys_extra_ghostcells)*(1-kr(idims,^D));

   {ixCRmin^D = max(ixCmin^D,ixGlo^D)\}
   {ixCRmax^D = min(ixCmax^D,ixGhi^D)\}

   ! ================================
   ! ====== Done with indices ======
   ! ================================
   
   ! Get interface positions:
   ! Use this form, the ending cell = 0
   !xi(kxC^S,1:ndim) = x(kxC^S,1:ndim)
   !xi(kxC^S,idims) = half* ( x(kxR^S,idims)+x(kxC^S,idims) )
   ! Instead we use this
   xi(ixI^S,1:ndim) = x(ixI^S,1:ndim)
   xi(ixI^S,idims)  = x(ixI^S,idims) + 0.5d0* mygeo%dx(ixI^S,idims)

   
   !wRp(kxC^S,1:nwflux+nwaux)=wprim(kxR^S,1:nwflux+nwaux)
   !wLp(kxC^S,1:nwflux+nwaux)=wprim(kxC^S,1:nwflux+nwaux)
   wRp(kxC^S,1:nw)=wprim(kxR^S,1:nw)
   wLp(kxC^S,1:nw)=wprim(kxC^S,1:nw)
   ! ================================
{#IFDEF DY_SP
   ! divB ~ 1e-18
   call metric_interpolation(ixI^L,ixI^L,idims,wprim,x,metric_i,xi)
   wLp(ixI^S,nmetric_lo:nmetric_hi) = metric_i(ixI^S,nmetric_lo:nmetric_hi)
   wRp(ixI^S,nmetric_lo:nmetric_hi) = metric_i(ixI^S,nmetric_lo:nmetric_hi)
   ! OLD way divB ~ 1e-7
   !call metric_interpolation(ixI^L,ixCR^L,idims,wprim,x,wLp,xi)
   !call metric_interpolation(ixI^L,ixCR^L,idims,wprim,x,wRp,xi)
   !call metric_interpolation(ixI^L,ixO^L^LADD3,idims,wprim,x,metric_i,xi)

   select case (idims)
   {case (^D)
   !psi6_surfC^D(ixI^S) = 1.0d0
   psi6_surfC^D(ixI^S) = metric_i(ixI^S, psi_metric_)**6
   mygeo%surfaceC^D(ixI^S)  = mygeo%surfaceC^D(ixI^S) * psi6_surfC^D(ixI^S)
   \}
   end select
}
   !KEN
   {#IFDEF M1
   call metricM1_L%fill_metric(wLp,x,ixI^L,ixI^L)
   call metricM1_R%fill_metric(wRp,x,ixI^L,ixI^L)
   }

   ! positivity preserving limiter always uses tvdlf and it can apply to tvdlf and hll cases for fC_low
   if (PP_limiter) then
     wRC(kxC^S,1:nw)=wCT(kxR^S,1:nw)
     wLC(kxC^S,1:nw)=wCT(kxC^S,1:nw)

     {#IFDEF STAGGERED
     ! Staggered grid, B-field do not need reconstruction, interface B_field = ws(B)
        {#IFNDEF DY_SP
          wLp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims); wRp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)
          wLC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims); wRC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)
        }
        {#IFDEF DY_SP
          ! Since ws has psi6 inside, we need pure B^i for calculating cbounds, getflux
          wLC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wLp(ixCR^S, psi_metric_)**6
          wRC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wRp(ixCR^S, psi_metric_)**6
          wLp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wLp(ixCR^S, psi_metric_)**6
          wRp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wRp(ixCR^S, psi_metric_)**6
        }
     }

     call get_cbounds(wLp,wRp,xi,ixI^L,ixC^L,idims,cmaxC,cminC,method,qtC, metricM1_L, metricM1_R)
     call getv(wLp,xi,ixI^L,ixC^L,idims,vLC)
     call getv(wRp,xi,ixI^L,ixC^L,idims,vRC)
     call getflux(wLC,wLp,xi,ixI^L,ixC^L,idims,fLC,transport,qtC, metricM1_L) !KEN
     call getflux(wRC,wRp,xi,ixI^L,ixC^L,idims,fRC,transport,qtC, metricM1_R) !KEN
     call get_Riemann_flux_tvdlf()
   
     do iw=1,ncons
           fC_low(ixC^S,iw,idims)=fC(ixC^S,iw,idims)
     end do ! Next iw
   end if ! end pp limiter


   if (method=='fd_tvdlf' .or. method=='fd_hll') then 
      select case (typelimited)
      case ('previous')
         call mpistop('never use previous typelimited')
         call primitive(ixI^L,ixI^L,wold,x)
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wold,wprim,wLC,wRC,wLp,wRp,x,xi,dxdim(idims), metricM1, metricM1_L, metricM1_R)
         call conserve(ixI^L,ixI^L,wold,x,patchfalse)
         {#IFDEF STAGGERED
         wLC(ixCR^S,b0_+idims)=wolds(ixCR^S,idims); wRC(ixCR^S,b0_+idims)=wolds(ixCR^S,idims)
         wLp(ixCR^S,b0_+idims)=wolds(ixCR^S,idims); wRp(ixCR^S,b0_+idims)=wolds(ixCR^S,idims)
         }
      case ('predictor') !default
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wprim,wprim,wLC,wRC,wLp,wRp,x,xi,dxdim(idims), metricM1, metricM1_L, metricM1_R)
         {#IFDEF STAGGERED
         ! Staggered grid, B-field do not need reconstruction, interface B_field = ws(B)
            {#IFNDEF DY_SP
              wLp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims); wRp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)
              wLC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims); wRC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)
            }
            {#IFDEF DY_SP
              ! Since ws has psi6 inside, we need pure B^i for calculating cbounds, getflux
              wLC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wLp(ixCR^S, psi_metric_)**6
              wRC(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wRp(ixCR^S, psi_metric_)**6
              wLp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wLp(ixCR^S, psi_metric_)**6
              wRp(ixCR^S,b0_+idims)=wCTs(ixCR^S,idims)/ wRp(ixCR^S, psi_metric_)**6
            }
         }

      case ('original')
         call mpistop('never use original typelimited')
         call primitive(ixI^L,ixI^L,wnew,x)
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wnew,wprim,wLC,wRC,wLp,wRp,x,xi,dxdim(idims), metricM1, metricM1_L, metricM1_R)
         call conserve(ixI^L,ixI^L,wnew,x,patchfalse)
         {#IFDEF STAGGERED
         wLC(ixCR^S,b0_+idims)=wnews(ixCR^S,idims); wRC(ixCR^S,b0_+idims)=wnews(ixCR^S,idims)
         wLp(ixCR^S,b0_+idims)=wnews(ixCR^S,idims); wRp(ixCR^S,b0_+idims)=wnews(ixCR^S,idims)
         }
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   else
      wRC(kxC^S,1:nwflux+nwaux) = wCT(kxR^S,1:nwflux+nwaux)
      wLC(kxC^S,1:nwflux+nwaux) = wCT(kxC^S,1:nwflux+nwaux)
   end if
   ! ==== Done with limiting ========
   !call imhd_modify_wLR(ixI^L,ixC^L,qt,wLp,wRp,sCT,idims)

   call get_cbounds(wLp,wRp,xi,ixI^L,ixC^L,idims,cmaxC,cminC,method, qtC, metricM1_L, metricM1_R)

   ! ================================
   ! Calculate velocities for transport fluxes
   ! ================================
   {#IFNDEF STAGGERED
   call getv(wLp,xi,ixI^L,ixC^L,idims,vLC)
   call getv(wRp,xi,ixI^L,ixC^L,idims,vRC)
   }
 
{#IFDEF STAGGERED
   ! ================================
   ! Calculate velocities for transport fluxes
   ! ================================
   call getv(wLp,xi,ixI^L,ixC^L,idims,vLC)
   call getv(wRp,xi,ixI^L,ixC^L,idims,vRC)
   ! ================================
   ! Storage for various ways of EMF averaging:
   ! ================================
   select case (typeemf)
   case default
   case ('average')
   case ('uct1')
      select case(method)
      case('fd_tvdlf')
        ! Store magnitude of characteristics
        cbarmax(ixC^S,idims)=max(cmaxC(ixC^S,d_),zero)
        cbarmin(ixC^S,idims)=cbarmax(ixC^S,idims)
      case('fd_hll')
        ! Store magnitude of characteristics
        cbarmax(ixC^S,idims)=max( cmaxC(ixC^S,d_),zero)
        cbarmin(ixC^S,idims)=max(-cminC(ixC^S,d_),zero)
      end select
      
      ! If the current direction will be useful for a component of the electric field
      ! (1 on 2D or 3 in 3D), then store relevant velocities.
      if (idims.le.2*ndim-3) then
         idimN=mod(idims,ndir)+1 ! 'Next' direction
         idimE=mod(idims+1,ndir)+1 ! Electric field direction

         vbarLC(ixC^S,idimE,1)=vLC(ixC^S)
         vbarRC(ixC^S,idimE,1)=vRC(ixC^S)
         ! Calculate additional component of velocity needed
         call getv(wLp,xi,ixI^L,ixC^L,idimN,vbarLC(ixI^S,idimE,2))
         call getv(wRp,xi,ixI^L,ixC^L,idimN,vbarRC(ixI^S,idimE,2))
      end if

   case ('uct2','uct2+av')
      select case(method)
      case('fd_tvdlf')
        ! Store magnitude of characteristics
        cbarmax(ixC^S,idims)=max(cmaxC(ixC^S,d_),zero)
        cbarmin(ixC^S,idims)=cbarmax(ixC^S,idims)
      case('fd_hll')
        ! Store magnitude of characteristics
        cbarmax(ixC^S,idims)=max( cmaxC(ixC^S,d_),zero)
        cbarmin(ixC^S,idims)=max(-cminC(ixC^S,d_),zero)
      end select

      idimN=mod(idims,ndir)+1 ! 'Next' direction
      idimE=mod(idims+1,ndir)+1 ! Electric field direction

      call getv(wLp,xi,ixI^L,ixC^L,idimN,vbarLC(ixI^S,idims,1))
      call getv(wRp,xi,ixI^L,ixC^L,idimN,vbarRC(ixI^S,idims,1))

      vbarC(ixC^S,idims,1)=half*(vbarLC(ixC^S,idims,1) &
           +vbarRC(ixC^S,idims,1))

      call getv(wLp,xi,ixI^L,ixC^L,idimE,vbarLC(ixI^S,idims,2))
      call getv(wRp,xi,ixI^L,ixC^L,idimE,vbarRC(ixI^S,idims,2))

      select case(method)
      case('fd_tvdlf')
        vbarC(ixC^S,idims,2)=half*(vbarLC(ixC^S,idims,2) &
             +vbarRC(ixC^S,idims,2))
      case('fd_hll')
        vbarC(ixC^S,idims,2)=(cbarmax(ixC^S,idims)*vbarLC(ixC^S,idims,2) &
             +cbarmin(ixC^S,idims)*vbarRC(ixC^S,idims,2))&
             /(cbarmax(ixC^S,idims) + cbarmin(ixC^S,idims))
      end select
   end select
   ! ================================
}
   
   ! ================================
   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   ! ================================
   call getflux(wLC,wLp,xi,ixI^L,ixC^L,idims,fLC,transport,qtC, metricM1_L) !KEN
   call getflux(wRC,wRp,xi,ixI^L,ixC^L,idims,fRC,transport,qtC, metricM1_R) !KEN
      
   select case(method)
   case('fd_tvdlf')
      call get_Riemann_flux_tvdlf()
   case('fd_hll')
      call get_Riemann_flux_hll()
   case default
      call mpistop('You are using new framework of BHAC, pls specific FD_HLL,FD_TVDLF')
   end select

   if (use_4th_order_fd .or. use_6th_order_fd) then
     if (.not. use_del_zanna) then
       ! get the cell center flux 
       call getv(wprim,x,ixI^L,ixI^L,idims,vCT)
       call getflux(wCT,wprim,x,ixI^L,ixI^L,idims,fCT,transport,qtC, metricM1_L) ! KEN i just use metricM1_L
       do iw=1, ncons
         if (transport(iw)) then
           fCT(ixI^S,iw) = fCT(ixI^S,iw)+vCT(ixI^S)*wCT(ixI^S,iw)
         endif
       enddo 
     endif
     call high_order_flux(ixI^L,ixC^L,idims,fC(ixI^S,1:nwflux,idims),fCT)
   endif
   ! ================================
   ! PP limiter with fC_low and fC to modify the flux
   if (PP_limiter) call positivity_preserving_limiter()

   do iw=1,ncons
   !do iw=1,nwflux
      if (slab) then
         !fC(ixC^S,iw,idims)=fC(ixC^S)
      else
         select case (idims)
            {case (^D)
            {#IFNDEF HARDBC
            !This where statement catches the axis where inverse metric becomes inf:
            where(mygeo%surfaceC^D(ixC^S) .gt. 1.0d-9*mygeo%dvolume(ixC^S))
             fC(ixC^S,iw,^D) = mygeo%surfaceC^D(ixC^S)*fC(ixC^S,iw,^D)
            elsewhere
               fC(ixC^S,iw,^D) = zero
            end where
            }{#IFDEF HARDBC
            if (typeaxial == 'cartesian') then
               fC(ixC^S,iw,^D)=mygeo%surfaceC^D(ixC^S)*fC(ixC^S,iw,^D)
            else
              if(idims.eq.^Z) then

                 where(abs(xi(ixC^S,idims)-xprobmin^Z).le.smalldouble .or. &
                      abs(xi(ixC^S,idims)-xprobmax^Z).le.smalldouble)
                    fC(ixC^S,iw,^D) = zero
                 elsewhere
                    fC(ixC^S,iw,^D)=mygeo%surfaceC^D(ixC^S)*fC(ixC^S,iw,^D)
                 end where

              else
                 fC(ixC^S,iw,^D)=mygeo%surfaceC^D(ixC^S)*fC(ixC^S,iw,^D)
              endif
            endif
            }
            \}
         end select
      end if

   end do ! Next iw

end do ! Next idims

{#IFNDEF STAGGERED
if (typeemf .eq. 'average') call fct_average(ixI^L,ixO^L,fC)
{#IFDEF HARDBC
if (typeaxial .ne. 'cartesian') then
  {^IFZIN
  do idims=1,ndir
     hxO^L=ixO^L-kr(^Z,^D);
     ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
     ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
  
     kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(^Z,^D);
     kxR^L=kxC^L+kr(^Z,^D);
     ! Get interface positions:
     xi(kxC^S,^Z) = half* ( x(kxR^S,^Z)+x(kxC^S,^Z) )
  
     where(abs(xi(ixC^S,^Z)-xprobmin^Z).le.smalldouble .or. &
          abs(xi(ixC^S,^Z)-xprobmax^Z).le.smalldouble)
        fC(ixC^S,b0_+idims,^Z) = zero
     end where
  end do
  }
endif
}
}{#IFDEF STAGGERED
  ! Before updatefaces, ws has psi6 inside already
  {#IFDEF DY_SP
    mygeo%dvolume(ixI^S)                   = mygeo%dvolume(ixI^S) / psi6(ixI^S)
    {^D& mygeo%surfaceC^D(ixI^S)           = mygeo%surfaceC^D(ixI^S) / psi6_surfC^D(ixI^S) \}
    ! The sqrtgamma inside should be Minkowski, since ws has psi6 inside
  }
  
  ! Update values at faces
  select case (typeemf)
  case default
     call mpistop('please specify typeemf for staggered field upwinding')
  case ('average')
     call mpistop('It will be very wrong when you use <average> for typeemf when turned on Staggered')
     call updatefaces(ixI^L,ixO^L,qdt,fC,fE,snew)
  case ('uct1')
     call updatefacesuct1(ixI^L,ixO^L,qdt,vbarRC,vbarLC,cbarmin,cbarmax,fE,snew)
  case ('uct2')
     call updatefacesuct2(ixI^L,ixO^L,qdt,vbarC,cbarmin,cbarmax,fE,snew)
  case ('uct2+av')
     call updatefacesuct2av(ixI^L,ixO^L,qdt,vbarC,cbarmin,cbarmax,fC,fE,snew)
  end select
  
  ! after updatefaces
  {#IFDEF DY_SP
    mygeo%dvolume(ixI^S)                   = mygeo%dvolume(ixI^S) * psi6(ixI^S)
    {^D& mygeo%surfaceC^D(ixI^S)           = mygeo%surfaceC^D(ixI^S) * psi6_surfC^D(ixI^S) \}
  }
}


! ================================
! Now update the state:
! ================================
do idims= idim^LIM

   if (.not. evolve_hydro) then
      fC(ixI^S, 1:ncons, :) = 0.0d0
   endif
   
   hxO^L=ixO^L-kr(idims,^D);
   
   do iw=1,ncons
   !do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)
         wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ixI^S,iw,^D)=-qdt*fC(ixI^S,iw,idims) !* mygeo%surfaceC^D(ixI^S)
           wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims



{#IFDEF STAGGERED
 {#IFDEF DY_SP
 ! Before updatefaces, ws has psi6 inside already
 mygeo%dvolume(ixI^S)                   = mygeo%dvolume(ixI^S) / psi6(ixI^S)
 {^D& mygeo%surfaceC^D(ixI^S)           = mygeo%surfaceC^D(ixI^S) / psi6_surfC^D(ixI^S) \}
 }
 call faces2centers(ixO^L,snew)
 ! after updatefaces
 {#IFDEF DY_SP
   mygeo%dvolume(ixI^S)                   = mygeo%dvolume(ixI^S) * psi6(ixI^S)
   {^D& mygeo%surfaceC^D(ixI^S)           = mygeo%surfaceC^D(ixI^S) * psi6_surfC^D(ixI^S) \}
 }
}
! ================================
! === Done updating state ========
! ================================


{#IFNDEF DY_SP
  if (covariant) myM => mygeo%m
}
if (.not.slab.and.idimmin==1) then 
   call addgeometry(qdt,ixI^L,ixO^L,wCT,wnew,wold,x, metricM1_L)
else
   call mpistop ('do not pass getaux in FD')
   if(nwaux>0) call getaux(.true.,wnew,x,ixI^L,ixO^L,'FD_new')
end if

call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixI^L,ixO^L,1,nw,qtC,&
                wCT,qt,wnew,x,.false.)

if (PP_limiter) then
   deallocate(fC_low)
end if

{#IFDEF DY_SP
  mygeo%dvolume(ixI^S)     = mygeo%dvolume(ixI^S) / psi6(ixI^S)
  {^D&  mygeo%surfaceC^D(ixI^S)  = mygeo%surfaceC^D(ixI^S) / psi6_surfC^D(ixI^S) \}
}

{#IFDEF STAGGERED
end associate staggered
}

call metricM1_L%destroy()
call metricM1_R%destroy()

end associate

contains 

   subroutine get_Riemann_flux_tvdlf()
   double precision  :: tvdlfeps_tmp =0.5
     
     do iw=1, ncons
       if (transport(iw)) then
          fadd(ixC^S)  = fLC(ixC^S,iw)+vLC(ixC^S)*wLC(ixC^S,iw)
          fdiss(ixC^S) = fRC(ixC^S,iw)+vRC(ixC^S)*wRC(ixC^S,iw)
       else
          fadd(ixC^S)  = fLC(ixC^S,iw)
          fdiss(ixC^S) = fRC(ixC^S,iw)
       end if

       ! Mean flux:
       fadd(ixC^S)=half*(fadd(ixC^S)+fdiss(ixC^S))
       ! Add TVDLF dissipation:
       fdiss(ixC^S)=-tvdlfeps*cmaxC(ixC^S,iw)*half*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
       ! fLC contains physical+dissipative fluxes
       fC(ixC^S,iw,idims)=fadd(ixC^S)+fdiss(ixC^S)
     enddo
      !KEN
     ! m1: correct fluxes form hll for diffusive lim
     !{#IFDEF M1 
     ! {do ix^DB=ixCmin^DB,ixCmax^DB\} 
     ! call m1_correct_asymptotic_fluxes(ix^D,ixI^L,tvdlfeps_tmp,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)    
     ! {end do\}
     !} 
     {#IFDEF M1
     call m1_correct_asymptotic_fluxes(ixI^L,ixC^L,tvdlfeps_tmp,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)   
     }
   end subroutine get_Riemann_flux_tvdlf

   subroutine get_Riemann_flux_hll()
       integer :: ix^D
     double precision  :: tvdlfeps_tmp =0.5
     do iw=1, ncons
       if (transport(iw)) then
          fLC(ixC^S,iw)  = fLC(ixC^S,iw)+vLC(ixC^S)*wLC(ixC^S,iw)
          fRC(ixC^S,iw)  = fRC(ixC^S,iw)+vRC(ixC^S)*wRC(ixC^S,iw)
       else
       end if

       {do ix^DB=ixCmin^DB,ixCmax^DB\}
         if(cminC(ix^D,iw) >= 0.0d0) then
           fC(ix^D,iw,idims)=fLC(ix^D,iw)
         else if(cmaxC(ix^D,iw) <= 0.0d0) then
           fC(ix^D,iw,idims)=fRC(ix^D,iw)
         else
           ! Add hll dissipation to the flux
           fC(ix^D,iw,idims)=(cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw)&
                 +tvdlfeps*cminC(ix^D,iw)&
                  *cmaxC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
                 /(cmaxC(ix^D,iw)-cminC(ix^D,iw))
         end if
       {end do\}

            ! m1: correct fluxes form hll for diffusive lim
       !{#IFDEF M1 
       !{do ix^DB=ixCmin^DB,ixCmax^DB\} 
       !call m1_correct_asymptotic_fluxes(ix^D,ixI^L,tvdlfeps,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)    
       !{end do\}
       !} 
       {#IFDEF M1
         call m1_correct_asymptotic_fluxes(ixI^L,ixC^L,tvdlfeps,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC) 
       }
      enddo

   end subroutine get_Riemann_flux_hll

   subroutine high_order_flux(ixI^L,ixC^L,idims,fC,fCT)
     integer, intent(in)                :: ixI^L, ixC^L, idims
     double precision, intent(in)       :: fCT(ixI^S,1:nwflux)
     double precision, intent(inout)    :: fC(ixI^S,1:nwflux)
     double precision, dimension(ixI^S) :: df2, df4
     integer                            :: ixCpp^L, ixCp^L, ixCm^L, ixCmm^L
     ixCp^L=ixC^L+kr(idims,^D);
     ixCpp^L=ixCp^L+kr(idims,^D);
     ixCm^L=ixC^L-kr(idims,^D);

     if (use_del_zanna) then
     ! have problem !fixme
     ! Del Zanna + 2007
       ixCmm^L=ixCm^L-kr(idims,^D);
       do iw = 1,nwflux
          df2(ixC^S) = ( fC(ixCm^S,iw) - 2.0d0 * fC(ixC^S,iw) + fC(ixCp^S,iw) ) / 2.4d1
          if (use_4th_order_fd) then
            fC(ixC^S,iw) = fC(ixC^S,iw) - df2(ixC^S)
          else if (use_6th_order_fd) then
            df4(ixC^S) = ( fC(ixCmm^S,iw) - 4.0d0 * fC(ixCm^S,iw) &
                         + 6.0d0 * fC(ixC^S,iw) &
                         + fC(ixCpp^S,iw) - 4.0d0 * fC(ixCp^S,iw) &
                             ) * 3.0d0 / 6.4d2
            fC(ixC^S,iw) = fC(ixC^S,iw) - df2(ixC^S) + df4(ixC^S)
          endif
       end do ! Next iw
     else
     ! Yuxi Chen+ 2016
       do iw = 1,nwflux
          df2(ixC^S) = ( fCT(ixC^S,iw) - 2.0d0 * fC(ixC^S,iw) + fCT(ixCp^S,iw) ) / 6.0d0
          if (use_4th_order_fd) then
            fC(ixC^S,iw) = fC(ixC^S,iw) - df2(ixC^S)
          else if (use_6th_order_fd) then
            df4(ixC^S) = ( fCT(ixCm^S,iw) - 9.0d0 * fCT(ixC^S,iw) &
                         + 1.6d1 * fC(ixC^S,iw) &
                         + fCT(ixCpp^S,iw) - 9.0d0 * fCT(ixCp^S,iw) &
                             ) / 1.8d2
            fC(ixC^S,iw) = fC(ixC^S,iw) - df2(ixC^S) + df4(ixC^S)
          endif
       end do ! Next iw
     endif
     
     ! m1: correct fluxes form hll for diffusive lim
   !  {#IFDEF M1 
   !  {do ix^DB=ixCmin^DB,ixCmax^DB\} 
   !  call m1_correct_asymptotic_fluxes(ix^D,ixI^L,tvdlfeps,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)    
   !  {end do\}
   !  } 

   end subroutine high_order_flux

   subroutine positivity_preserving_limiter()
     use mod_eos, only: small_rho
     implicit none
     double precision, dimension(ixI^S)   :: inv_lambda, epsD
     double precision, dimension(ixI^S)   :: theta, theta_tau
     double precision, dimension(ixI^S)   :: flow_dA, fhigh_dA
     double precision, parameter          :: eps_fac = 1.0d-6
     integer                              :: ix^D

     {do ix^D = ixI^LIM^D \}
        epsD(ix^D)   = eps_fac * small_rho
     {end do^D&\}

     if (slab) then
        inv_lambda(ixI^S) = dxdim(idims)/(qdt)
        flow_dA(ixI^S)  = fC_low(ixI^S,D_,idims)
        fhigh_dA(ixI^S) = fC(ixI^S,D_,idims)
     else
        inv_lambda(ixI^S) = mygeo%dvolume(ixI^S)/(qdt)
        select case (idims)
        {case (^D)
        flow_dA(ixI^S)  = fC_low(ixI^S,D_,idims)   * mygeo%surfaceC^D(ixI^S)  
        fhigh_dA(ixI^S) = fC(ixI^S,D_,idims) * mygeo%surfaceC^D(ixI^S)
        \}
        end select
     end if

     call get_theta(ixI^L,ixC^L,idims,epsD,inv_lambda(ixI^S),sCT%w%w(ixI^S,d_),flow_dA,fhigh_dA,theta)
     ! tau_ may be -ve for tabulate EOS
     !call get_theta(ixI^L,ixC^L,idims,epstau,inv_lambda(ixI^S),sCT%cons(ixI^S,tau_),fC_low(ixI^S,tau_,idims),fC(ixI^S,tau_,idims),theta_tau)

     !theta(ixC^S) = min(theta(ixC^S),theta_tau(ixC^S))

     do iw = 1,ncons
        ! note that pp limiter cannot act on stagger grid
        {#IFDEF STAGGERED
           if ( (iw>=b1_) .and. (iw<=b^NC_) ) cycle
        }
        fC(ixC^S,iw,idims) = theta(ixC^S)*fC(ixC^S,iw,idims) &
             + (1.0d0-theta(ixC^S))*fC_low(ixC^S,iw,idims)
     end do ! Next iw
   end subroutine positivity_preserving_limiter

   subroutine get_theta(ixI^L,ixC^L,idims,eps,inv_lambda,u,flow_dA,fhigh_dA,theta)
     integer, intent(in)                             :: ixI^L, ixC^L, idims
     double precision, intent(in)                    :: eps(ixI^S)
     double precision, dimension(ixI^S), intent(in)  :: inv_lambda, u
     double precision, dimension(ixI^S), intent(in)  :: flow_dA, fhigh_dA
     double precision, dimension(ixI^S), intent(out) :: theta
   
     integer                                         :: ixCp^L, ixOp^L
     double precision, dimension(ixI^S)              :: tmp, thp, thm
     double precision, dimension(ixI^S)              :: diff_fdA
   
     ! Note: here we assume that u( i=0 ) is given
     ixCp^L=ixC^L+kr(idims,^D);
     ixOpmin^D=ixCmin^D; ixOpmax^D=ixCpmax^D;
   
     thm(ixC^S) = 1.0d0
     thp(ixC^S) = 1.0d0
   
     tmp(ixOp^S) = 0.5d0*inv_lambda(ixOp^S)*(u(ixOp^S)-eps(ixOp^S))
   
     diff_fdA(ixC^S) = -flow_dA(ixC^S) + fhigh_dA(ixC^S)
     where (diff_fdA(ixC^S) == 0.0d0)
        diff_fdA(ixC^S) = smalldouble ! avoid flow = fhight case
     end where
   
     where (tmp(ixC^S) < fhigh_dA(ixC^S))
        thm(ixC^S) = tmp(ixC^S) - flow_dA(ixC^S)
        thm(ixC^S) = thm(ixC^S) / (diff_fdA(ixC^S))
     end where
   
     where (tmp(ixCp^S) < -fhigh_dA(ixC^S))
        thp(ixC^S) = - tmp(ixCp^S) - flow_dA(ixC^S)
        thp(ixC^S) = thp(ixC^S) / (diff_fdA(ixC^S))
     end where
   
     theta(ixC^S) = min(thm(ixC^S),thp(ixC^S))
     theta(ixC^S) = min(max(theta(ixC^S),0.0d0),1.0d0)
   end subroutine get_theta

end subroutine finite_difference
