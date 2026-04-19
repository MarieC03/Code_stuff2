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

!=============================================================================
subroutine upwindLR(ixI^L,ixL^L,ixR^L,idims,w,wCT,wLC,wRC,wLp,wRp,x,xi,dxdim, metricM1, metricM1_L, metricM1_R) !KEN

! Determine the upwinded wLC(ixL) and wRC(ixR) from w. 
! the wCT is only used when PPM is exploited.
use mod_eos
use mod_limiter
use mod_imhd_intermediate
use mod_variables
use mod_small_values
use mod_ppm, only:  PPMflattening
use mod_m1_metric_interface
use mod_m1_closure
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
double precision, intent(in) :: dxdim
double precision, dimension(ixI^S,1:nw) :: w, wCT
double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S,1:nw) :: wLp, wRp
double precision, dimension(ixI^S,1:ndim) :: x
double precision, dimension(ixI^S,1:ndim) :: xi
type(m1_metric_helper), intent(in)   :: metricM1, metricM1_L, metricM1_R !KEN
!type(m1_metric_helper)   :: metricM1 !KEN
! .. local ..
integer :: jxR^L, ixC^L, jxC^L, iw, ix^D, ix^L
double precision :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
double precision :: fmag(ixI^S), scaleF(ixI^S)
double precision :: fmaxfact = 0.999d0
logical          :: fattening = .False.

{#IFDEF M1
!type(m1_metric_helper)   :: metricM1_L, metricM1_R, metricM1
double precision :: M1_mini = 1.d-45
!integer :: ix^D
}
!-----------------------------------------------------------------------------
!    w(ixI^S,eps_)  = w(ixI^S,eps_)  * w(ixI^S, rho_)
!    wLp(ixI^S,eps_) = wLp(ixI^S,eps_) * wLp(ixI^S, rho_)
!    wRp(ixI^S,eps_) = wRp(ixI^S,eps_) * wRp(ixI^S, rho_)
!{do ix^D = ixL^LIM^D \}
!call Simple_check_data_correctness_pt(w(ix^D,1:nw), xi(ix^D,1:ndim), w(ix^D,1:nw), 'w b4 reconstruction small values')
!{end do\}
!write(*,*) 'b4 reconstruction'

 !-------- M1 for reconstruction ---------------------------------------------
{#IFDEF M1
   {#IFDEF UNIT_TESTS
   call fill_metric(metricM1_L)
   call fill_metric(metricM1_R)
   call fill_metric(metricM1)
   }
   {#IFNDEF UNIT_TESTS
    !call metricM1_L%fill_metric(wLp,x,ixI^L,ixL^L)
    !call metricM1_R%fill_metric(wRp,x,ixI^L,ixR^L)
    !call metricM1%fill_metric(w,x,ixI^L,ixI^L)! KEN
    !KEN Additionaly metricM1=metricM1_L
   }
   {^KSP&
     w(ixI^S,erad^KSP_) = (w(ixI^S,erad^KSP_)) / metricM1%sqrtg(ixI^S) + M1_mini !KEN KEN wrote it behind, and then I removed the M1_mini
     w(ixI^S,erad^KSP_) = max(w(ixI^S,erad^KSP_),m1_E_atmo)
     w(ixI^S,nrad^KSP_) = max(w(ixI^S,nrad^KSP_), m1_E_atmo) /(w(ixI^S,erad^KSP_)) / metricM1%sqrtg(ixI^S)
     !w(ixI^S,nrad^KSP_) = max(w(ixI^S,nrad^KSP_) / metricM1%sqrtg(ixI^S), m1_E_atmo) /(w(ixI^S,erad^KSP_))
     !!w(ixI^S,nrad^KSP_) = max(w(ixI^S,nrad^KSP_),m1_E_atmo)
     {^C& w(ixI^S,frad^KSP^C_) = w(ixI^S,frad^KSP^C_) / metricM1%sqrtg(ixI^S) / (w(ixI^S,erad^KSP_) )\}
     
    !call m1_update_closure(metricM1, w, x, ixI^L, ixO^L,^KSP,.true.,Gamma)

     wLp(ixL^S,erad^KSP_) = (wLp(ixL^S,erad^KSP_)) / metricM1_L%sqrtg(ixL^S) + M1_mini
     wLp(ixL^S,erad^KSP_) = max(wLp(ixL^S,erad^KSP_),m1_E_atmo)
     wLp(ixL^S,nrad^KSP_) = max(wLp(ixL^S,nrad^KSP_),m1_E_atmo) /(wLp(ixL^S,erad^KSP_) ) / metricM1_L%sqrtg(ixL^S)
     !wLp(ixL^S,nrad^KSP_) = max(wLp(ixL^S,nrad^KSP_) / metricM1_L%sqrtg(ixL^S),m1_E_atmo) /(wLp(ixL^S,erad^KSP_) )
     {^C& wLp(ixL^S,frad^KSP^C_) = wLp(ixL^S,frad^KSP^C_) / metricM1_L%sqrtg(ixL^S) / (wLp(ixL^S,erad^KSP_)) \}

     wRp(ixR^S,erad^KSP_) = (wRp(ixR^S,erad^KSP_)) / metricM1_R%sqrtg(ixR^S) + M1_mini
     wRp(ixR^S,erad^KSP_) = max(wRp(ixR^S,erad^KSP_),m1_E_atmo)
     wRp(ixR^S,nrad^KSP_) = max(wRp(ixR^S,nrad^KSP_), m1_E_atmo) /(wRp(ixR^S,erad^KSP_)) / metricM1_R%sqrtg(ixR^S)
     !wRp(ixR^S,nrad^KSP_) = max(wRp(ixR^S,nrad^KSP_) / metricM1_R%sqrtg(ixR^S), m1_E_atmo) /(wRp(ixR^S,erad^KSP_))
     {^C& wRp(ixR^S,frad^KSP^C_) = wRp(ixR^S,frad^KSP^C_) / metricM1_R%sqrtg(ixR^S) / (wRp(ixR^S,erad^KSP_))\}
   \}  
   ! now wLp,wRp is : nGamma/E, E, F_i/E   
}
!-----------------------------------------------------------------------------

    do iw=1, nwflux
       if (.not. var_reconstruct(iw) ) cycle

       !{#IFDEF DY_SP
       if (iw == b1_ .or. iw == b2_ .or. iw == b3_) then
          call PPMlimiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),PPM_extrema)
       else
       !}
         select case (typelimiter)
         case (limiter_mp5)
            call MP5limiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw))
         case (limiter_weno5)
            call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),1,.false.)
         case (limiter_wenoZ5)
            call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),2,.false.)
         case (limiter_wenoZP)
            call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),3,.false.)
         case (limiter_ppm)
            !fattening = .True. ! fixme: remove this
            ! Since PPM uses the ordinary grid-index:
            !ixCmin^D=ixLmin^D+kr(^D,idims);
            !ixCmax^D=ixLmax^D;
            ! our fattening is only available for hydro
            call PPMlimiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),PPM_extrema)
            !call PPMlimiter(ixI^L,ixC^L,idims,w(ixI^S,iw),wCT(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),fattening)
         case default

            jxR^L=ixR^L+kr(idims,^D);
            ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
            jxC^L=ixC^L+kr(idims,^D);
  
            if (loglimit(iw)) then
               w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
               wLp(ixL^S,iw)=dlog10(wLp(ixL^S,iw))
               wRp(ixR^S,iw)=dlog10(wRp(ixR^S,iw))
            end if
  
            ! original version
            dwC(ixC^S)=( w(jxC^S,iw)-w(ixC^S,iw) )
            ! limit flux from left and/or right
            call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter,ldw,rdw)
            wLp(ixL^S,iw)=wLp(ixL^S,iw)+0.5d0*ldw(ixL^S)
            wRp(ixR^S,iw)=wRp(ixR^S,iw)-0.5d0*rdw(jxR^S)
  
            ! new version
            !dxi(ixC^S) = dxdim!xi(jxC^S,idims) - xi(ixC^S,idims)
            !dwC(ixC^S)=( w(jxC^S,iw)-w(ixC^S,iw) ) &
            !            / (x(jxC^S,idims) - x(ixC^S,idims)) &
            !             * dxdim
  
            !! limit flux from left and/or right
            !call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter,ldw,rdw)
            !wLp(ixL^S,iw)=w(ixL^S,iw)+ldw(ixL^S)*(xi(ixL^S,idims)-x(ixL^S,idims))/dxdim
            !wRp(ixR^S,iw)=w(jxR^S,iw)+rdw(jxR^S)*(xi(ixR^S,idims)-x(jxR^S,idims))/dxdim
  
            if (loglimit(iw)) then
               w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
               wLp(ixL^S,iw)=10.0d0**wLp(ixL^S,iw)
               wRp(ixR^S,iw)=10.0d0**wRp(ixR^S,iw)
            end if
           !code test since minmod for b field always NaN
           {^C& wLp(ixL^S,b^C_) = 0.0d0\} 
           {^C& wRp(ixR^S,b^C_) = 0.0d0\}
         end select
         !{#IFDEF DY_SP
         endif ! endif of using PPM for b^C
         !}
  
      end do            

!{do ix^D = ixL^LIM^D \}
!call Simple_check_data_correctness_pt(wLp(ix^D,1:nw), x(ix^D,1:ndim), wLp(ix^D,1:nw), 'wLp after handle small values')
!{end do\}
!{do ix^D = ixR^LIM^D \}
!call Simple_check_data_correctness_pt(wRp(ix^D,1:nw), x(ix^D,1:ndim), wRp(ix^D,1:nw), 'wRp after handle small values')
!{end do\}

if (fixp_after_reconstruct) then
   ! fixp usr has eos calls with different EOS for the atmosphere
   if (tlow>zero) then
    !call fixp_usr(ixI^L,ixL^L,wLp,xi)
    !call fixp_usr(ixI^L,ixR^L,wRp,xi)
    call fixp_usr(ixI^L,ixL^L,wLp,x)
    call fixp_usr(ixI^L,ixR^L,wRp,x)
   endif
 else
  call imhd_handle_small_values(wLp,xi,ixI^L,ixL^L,.false.)
  call imhd_handle_small_values(wRp,xi,ixI^L,ixR^L,.false.)
  !code test
  !{do ix^D = ixL^LIM^D \}
  !call Simple_check_data_correctness_pt(wLp(ix^D,1:nw), xi(ix^D,1:ndim), wLp(ix^D,1:nw), 'wLp after handle small values')
  !{end do\}
  !{do ix^D = ixR^LIM^D \}
  !call Simple_check_data_correctness_pt(wRp(ix^D,1:nw), xi(ix^D,1:ndim), wRp(ix^D,1:nw), 'wRp after handle small values')
  !{end do\}
  call Eos_update_one_grid(ixI^L, ixL^L, wLp(ixI^S,1:nw))
  call Eos_update_one_grid(ixI^L, ixR^L, wRp(ixI^S,1:nw))
 end if
 
 
 !call Eos_update_one_grid(ixI^L, ixL^L, wLp(ixI^S,1:nw))
 !call Eos_update_one_grid(ixI^L, ixR^L, wRp(ixI^S,1:nw))
       

if (typelimiter == limiter_ppm .and. PPM_flatten) &
     call PPMflattening(ixI^L,ixL^L,idims,w,wLp,wRp,PPM_flatcd,PPM_flatsh)

     !----------- M1 after reconstruction  ------------------------------------
!------------ here multiply by E sqrtg again ------------
! conserven() does not do anything to rad-vars so we have to 
! multiply sqrtg by ourselves, since in getcbouds wLp and wRp is used as sqrtg*prim
{#IFDEF M1
    {^KSP&
       w(ixI^S,erad^KSP_) = max(w(ixI^S,erad^KSP_),m1_E_atmo)
       w(ixI^S,nrad^KSP_) = max(w(ixI^S,nrad^KSP_) * metricM1%sqrtg(ixI^S) * (w(ixI^S,erad^KSP_)) ,m1_E_atmo)!KEN KEN removed M1Mini. No safeguard is needed for mult.
       !w(ixI^S,nrad^KSP_) = max(w(ixI^S,nrad^KSP_),m1_E_atmo) * metricM1%sqrtg(ixI^S) * (w(ixI^S,erad^KSP_)) !KEN KEN removed M1Mini. No safeguard is needed for mult.
       {^C& w(ixI^S,frad^KSP^C_) = w(ixI^S,frad^KSP^C_) * metricM1%sqrtg(ixI^S) * (w(ixI^S,erad^KSP_)) \}
       w(ixI^S,erad^KSP_) = (w(ixI^S,erad^KSP_)) * metricM1%sqrtg(ixI^S)

      wLp(ixL^S,erad^KSP_) = max(wLp(ixL^S,erad^KSP_),m1_E_atmo)
      wLp(ixL^S,nrad^KSP_) = max(wLp(ixL^S,nrad^KSP_)* metricM1_L%sqrtg(ixL^S) * (wLp(ixL^S,erad^KSP_)),m1_E_atmo) 
      !wLp(ixL^S,nrad^KSP_) = max(wLp(ixL^S,nrad^KSP_),m1_E_atmo) * metricM1_L%sqrtg(ixL^S) * (wLp(ixL^S,erad^KSP_))
      {^C& wLp(ixL^S,frad^KSP^C_) = wLp(ixL^S,frad^KSP^C_) * metricM1_L%sqrtg(ixL^S) * (wLp(ixL^S,erad^KSP_)) \}
      wLp(ixL^S,erad^KSP_) = (wLp(ixL^S,erad^KSP_)) * metricM1_L%sqrtg(ixL^S)
 
      wRp(ixL^S,erad^KSP_) = max(wRp(ixL^S,erad^KSP_),m1_E_atmo)
      wRp(ixR^S,nrad^KSP_) = max(wRp(ixR^S,nrad^KSP_) * metricM1_R%sqrtg(ixR^S) * (wRp(ixR^S,erad^KSP_)),m1_E_atmo)
      !wRp(ixR^S,nrad^KSP_) = max(wRp(ixR^S,nrad^KSP_),m1_E_atmo) * metricM1_R%sqrtg(ixR^S) * (wRp(ixR^S,erad^KSP_))
      {^C& wRp(ixR^S,frad^KSP^C_) = wRp(ixR^S,frad^KSP^C_) * metricM1_R%sqrtg(ixR^S) * (wRp(ixR^S,erad^KSP_)) \}
      wRp(ixR^S,erad^KSP_) = (wRp(ixR^S,erad^KSP_)) * metricM1_R%sqrtg(ixR^S)
 
    \}  
    ! now wLp,wRp is : n*Gamma*sqrtg, E*sqrtg, F_i*sqrtg  
 }
 !-----------------------------------------------------------------------------
 

  wLC(ixL^S,1:nw) = wLp(ixL^S,1:nw)
  wRC(ixR^S,1:nw) = wRp(ixR^S,1:nw)
  call conserve(ixI^L,ixL^L,wLC,xi,patchfalse)   ! in DYSP, p2c depends on xi
  call conserve(ixI^L,ixR^L,wRC,xi,patchfalse)   
  wLp(ixL^S,1:nw) = wLC(ixL^S,1:nw)
  wRp(ixR^S,1:nw) = wRC(ixR^S,1:nw)
  ! wRC wLC have all newest metric, prim cons aux 

end subroutine upwindLR

!============================================================================
! sCT current w state as input for calculating snew;  snew is the output new w state
subroutine finite_volume(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,sCT,qt,snew,sold,wprim,fC,fE,dx^D,x,sm1)
use mod_interpolate
use mod_imhd_intermediate
{#IFDEF M1
use mod_m1, only: m1_correct_asymptotic_fluxes
!use mod_m1_metric_interface
}
use mod_limiter
use mod_m1_metric_interface

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
double precision, dimension(ixI^S,1:^NC)               :: fE

double precision, dimension(ixI^S,1:ndim)              :: xi
double precision, dimension(ixI^S,1:nw)                :: wLC, wRC, wLp, wRp, wi
double precision, dimension(ixI^S,1:nwflux)            :: fLC, fRC
double precision, dimension(ixI^S)                     :: fadd, fdiss
double precision, dimension(ixI^S)                     :: vLC, vRC
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
integer                                                :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L, ixtest^L
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


if (covariant .and. typelimited .ne. 'predictor') call mpistop('FV: Covariant formulation only implemented with typelimited=predictor')


!if (idimmax>idimmin .and. typelimited=='original' .and. &
!   method/='tvdlf1')&
!   call mpistop("Error in finite volume: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in finite volume: Nonconforming input limits")

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
   !code test
   !ixCR^L=ixC^L; ! indices to reconstruct to
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
   
   ! ================================
   ! === Start with flat (Godunov): =
   ! ================================
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
   call metricM1_L%fill_metric(wLp,x,ixI^L,ixI^L, get_deriv=.true.)
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
     
     call get_cbounds(wLp,wRp,xi,ixI^L,ixC^L,idims,cmaxC,cminC,method, qtC, metricM1_L, metricM1_R)

     call getv(wLp,xi,ixI^L,ixC^L,idims,vLC)
     call getv(wRp,xi,ixI^L,ixC^L,idims,vRC)
     call getflux(wLC,wLp,xi,ixI^L,ixC^L,idims,fLC,transport, qtC, metricM1_L) !KEN
     call getflux(wRC,wRp,xi,ixI^L,ixC^L,idims,fRC,transport, qtC, metricM1_R) !KEN
     call get_Riemann_flux_tvdlf()
   
     do iw=1,ncons
           fC_low(ixC^S,iw,idims)=fC(ixC^S,iw,idims)
     end do ! Next iw
   end if ! end pp limiter

   if (method=='tvdlf' .or. method=='hll') then  ! if tvdlf1 or hll1 means no reconstruction (first order FV)
      select case (typelimited)
      case ('previous')
         call mpistop('never use previous typelimited')
         call primitive(ixI^L,ixI^L,wold,x)
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wold,wprim,wLC,wRC,wLp,wRp,x,xi,dxdim(idims), metricM1,metricM1_L,metricM1_R)
         call conserve(ixI^L,ixI^L,wold,x,patchfalse)
         {#IFDEF STAGGERED
         wLC(ixCR^S,b0_+idims)=wolds(ixCR^S,idims); wRC(ixCR^S,b0_+idims)=wolds(ixCR^S,idims)
         wLp(ixCR^S,b0_+idims)=wolds(ixCR^S,idims); wRp(ixCR^S,b0_+idims)=wolds(ixCR^S,idims)
         }
      case ('predictor') !default
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wprim,wprim,wLC,wRC,wLp,wRp,x,xi,dxdim(idims),metricM1,metricM1_L,metricM1_R)
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
         call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wnew,wprim,wLC,wRC,wLp,wRp,x,xi,dxdim(idims),metricM1,metricM1_L,metricM1_R)
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

   select case(method)
   case('tvdlf','tvdlf1')
      ! ================================
      ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
      ! the maximum eigenvalue, it is calculated in advance.
      ! ================================
      call get_cbounds(wLp,wRp,xi,ixI^L,ixC^L,idims,cmaxC,cminC,method,qtC, metricM1_L, metricM1_R)

      ! ================================
      ! Calculate velocities for transport fluxes
      ! ================================
      {#IFNDEF STAGGERED
      call getv(wLp,xi,ixI^L,ixC^L,idims,vLC)
      call getv(wRp,xi,ixI^L,ixC^L,idims,vRC)
      }
   case('hll','hll1')
      ! ================================
      ! For the high order hll scheme the limiter is based on
      ! the maximum eigenvalue, it is calculated in advance.
      ! ================================
      call get_cbounds(wLp,wRp,xi,ixI^L,ixC^L,idims,cmaxC,cminC,method,qtC, metricM1_L, metricM1_R)

      ! ================================
      ! Calculate velocities for transport fluxes
      ! ================================
      {#IFNDEF STAGGERED
      !if(any(patchf(ixC^S)/= 2).or.(logiB)) call getv(wLp,xi,ixI^L,ixC^L,idims,vLC)
      !if(any(patchf(ixC^S)/=-2).or.(logiB)) call getv(wRp,xi,ixI^L,ixC^L,idims,vRC)
      call getv(wLp,xi,ixI^L,ixC^L,idims,vLC)
      call getv(wRp,xi,ixI^L,ixC^L,idims,vRC)
      }
   case default
         call mpistop('You are using new framework of BHAC, pls specific HLL, or TVDLF')
   end select
 
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
      case('tvdlf','tvdlf1')
        ! Store magnitude of characteristics
        cbarmax(ixC^S,idims)=max(cmaxC(ixC^S,d_),zero)
        cbarmin(ixC^S,idims)=cbarmax(ixC^S,idims)
      case('hll','hll1')
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
      case('tvdlf','tvdlf1')
        ! Store magnitude of characteristics
        cbarmax(ixC^S,idims)=max(cmaxC(ixC^S,d_),zero)
        cbarmin(ixC^S,idims)=cbarmax(ixC^S,idims)
      case('hll','hll1')
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
      case('tvdlf','tvdlf1')
        vbarC(ixC^S,idims,2)=half*(vbarLC(ixC^S,idims,2) &
             +vbarRC(ixC^S,idims,2))
      case('hll','hll1')
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
   call getflux(wLC,wLp,xi,ixI^L,ixC^L,idims,fLC,transport, qtC, metricM1_L) !KEN
   call getflux(wRC,wRp,xi,ixI^L,ixC^L,idims,fRC,transport, qtC, metricM1_R) !KEN

   select case(method)
   case('tvdlf','tvdlf1')
      call get_Riemann_flux_tvdlf()
   case('hll','hll1')
      call get_Riemann_flux_hll()
   case default
      call mpistop('You are using new framework of BHAC, pls specific HLL,TVDLF')
   end select

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

if (any(wprim(ixO^S,20) .lt. 0.0)) then
   print *, "nrad is negative before finite volume geometry"
endif

if (any(wnew(ixO^S,20) .lt. 0.0)) then
   print *, "nrad is negative before add geometry"
endif

{#IFNDEF DY_SP
  if (covariant) myM => mygeo%m
}
if (.not.slab.and.idimmin==1) then 
   call addgeometry(qdt,ixI^L,ixO^L,wCT,wnew,wold,x,qtC, metricM1_L)
else
   call mpistop ('do not pass getaux in FV')
   if(nwaux>0) call getaux(.true.,wnew,x,ixI^L,ixO^L,'FV_new')
end if

if (any(wnew(ixO^S,20) .le. 0)) then
   print *, "nrad is negative after add geometry"
endif

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
     integer :: ix^D
     double precision  :: tvdlfeps_tmp =0.5
     
     do iw=1, ncons
       !if (iw .eq. nrad1_) exit !KEN. We have our correct asymptotic flux
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
     !call m1_correct_asymptotic_fluxes(ixI^L,ixC^L,tvdlfeps_tmp,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)    

   end subroutine get_Riemann_flux_tvdlf


   subroutine get_Riemann_flux_hll()
       integer :: ix^D

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
     enddo ! 1:nwcons 

     ! m1: correct fluxes form hll for diffusive lim
     !{#IFDEF M1 
     ! {do ix^DB=ixCmin^DB,ixCmax^DB\} 
     ! call m1_correct_asymptotic_fluxes(ix^D,ixI^L,tvdlfeps,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)    
     ! {end do\}
     !} 
     call m1_correct_asymptotic_fluxes(ixI^L,ixC^L,tvdlfeps,dxdim,idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)    

   end subroutine get_Riemann_flux_hll


   subroutine positivity_preserving_limiter()
     use mod_eos, only: small_rho
     implicit none
   
     double precision, parameter :: eps_fac = 1.0d-6
     double precision, dimension(ixI^S) :: inv_lambda
     double precision, dimension(ixI^S) :: epsD, theta_hydro
     double precision, dimension(ixI^S) :: flow_dA, fhigh_dA
   
     {#IFDEF M1
     double precision, dimension(ixI^S) :: epsR, theta_n, theta_e, theta_rad
     }
   
     integer :: iw
     !-----------------------------------------------------------------------------
   
     !--------------------------
     ! eps for hydro density D_
     !--------------------------
     epsD(ixI^S) = eps_fac * small_rho
   
     !-----------------------------------------
     ! inv_lambda and fluxes must be in *update*
     ! form consistent with get_theta()
     !-----------------------------------------
     if (slab) then
        inv_lambda(ixI^S) = dxdim(idims)/qdt
        flow_dA(ixI^S)    = fC_low(ixI^S,D_,idims)
        fhigh_dA(ixI^S)   = fC(ixI^S,D_,idims)
     else
        inv_lambda(ixI^S) = mygeo%dvolume(ixI^S)/qdt
        select case (idims)
        {case (^D)
           flow_dA(ixI^S)  = fC_low(ixI^S,D_,idims) * mygeo%surfaceC^D(ixI^S)
           fhigh_dA(ixI^S) = fC(ixI^S,D_,idims)     * mygeo%surfaceC^D(ixI^S)
        \}
        end select
     end if
   
     call get_theta(ixI^L,ixC^L,idims,epsD,inv_lambda(ixI^S),sCT%w%w(ixI^S,D_), &
                    flow_dA,fhigh_dA,theta_hydro)
   
     !---------------------------------------------------------
     ! Blend high-order flux -> low-order flux for hydro cons.
     ! Note: do NOT "exit" early; skip staggered B as before.
     !---------------------------------------------------------
     do iw = 1, ncons
        {#IFDEF STAGGERED
           if ((iw>=b1_) .and. (iw<=b^NC_)) cycle
        }
        {#IFDEF M1
           !if ( (iw>=nrad1_) .and. (iw<=frad^NS3_) ) cycle
           if (iw>=nrad1_) exit
        }
        fC(ixC^S,iw,idims) = theta_hydro(ixC^S)*fC(ixC^S,iw,idims) &
                           + (1.0d0-theta_hydro(ixC^S))*fC_low(ixC^S,iw,idims)
     end do
   
     !======================================================================
     !========================  M1 radiation part  ==========================
     !======================================================================
     {#IFDEF M1
   
     ! Choose epsR in the SAME UNITS as the CONSERVED variables.
     ! In BHAC-M1 your evolved nrad/erad are densitized already, so m1_E_atmo
     ! is a reasonable minimal admissible value in conserved form.
     epsR(ixI^S) = m1_E_atmo
   
     {^KSP&
       !-------------------------------------------------------------
       ! Build theta_rad from BOTH nrad and erad positivity constraints
       !-------------------------------------------------------------
   
       ! --- nrad constraint ---
       if (slab) then
          flow_dA(ixI^S)    = fC_low(ixI^S,nrad^KSP_,idims)
          fhigh_dA(ixI^S)   = fC(ixI^S,nrad^KSP_,idims)
       else
          select case (idims)
          {case (^D)
             flow_dA(ixI^S)  = fC_low(ixI^S,nrad^KSP_,idims) * mygeo%surfaceC^D(ixI^S)
             fhigh_dA(ixI^S) = fC(ixI^S,nrad^KSP_,idims)     * mygeo%surfaceC^D(ixI^S)
          \}
          end select
       end if
   
       call get_theta(ixI^L,ixC^L,idims,epsR,inv_lambda(ixI^S),sCT%w%w(ixI^S,nrad^KSP_), &
                      flow_dA,fhigh_dA,theta_n)
   
       ! --- erad constraint ---
       if (slab) then
          flow_dA(ixI^S)    = fC_low(ixI^S,erad^KSP_,idims)
          fhigh_dA(ixI^S)   = fC(ixI^S,erad^KSP_,idims)
       else
          select case (idims)
          {case (^D)
             flow_dA(ixI^S)  = fC_low(ixI^S,erad^KSP_,idims) * mygeo%surfaceC^D(ixI^S)
             fhigh_dA(ixI^S) = fC(ixI^S,erad^KSP_,idims)     * mygeo%surfaceC^D(ixI^S)
          \}
          end select
       end if
   
       call get_theta(ixI^L,ixC^L,idims,epsR,inv_lambda(ixI^S),sCT%w%w(ixI^S,erad^KSP_), &
                      flow_dA,fhigh_dA,theta_e)
   
       ! Combine: enforce BOTH constraints
       theta_rad(ixC^S) = min(theta_n(ixC^S), theta_e(ixC^S))
       theta_rad(ixC^S) = min(max(theta_rad(ixC^S), 0.0d0), 1.0d0)
   
       !-------------------------------------------------------------
       ! Apply SAME theta_rad to all moments of this radiation species
       ! (keeps the system consistent; avoids breaking realizability
       !  by mixing moments with different thetas).
       !-------------------------------------------------------------
       fC(ixC^S,nrad^KSP_,idims) = theta_rad(ixC^S)*fC(ixC^S,nrad^KSP_,idims) &
                                + (1.0d0-theta_rad(ixC^S))*fC_low(ixC^S,nrad^KSP_,idims)
   
       fC(ixC^S,erad^KSP_,idims) = theta_rad(ixC^S)*fC(ixC^S,erad^KSP_,idims) &
                                + (1.0d0-theta_rad(ixC^S))*fC_low(ixC^S,erad^KSP_,idims)
   
       {^C&
       fC(ixC^S,frad^KSP^C_,idims) = theta_rad(ixC^S)*fC(ixC^S,frad^KSP^C_,idims) &
                                  + (1.0d0-theta_rad(ixC^S))*fC_low(ixC^S,frad^KSP^C_,idims)
       \}
   
     \} ! end KSP loop
   
     } ! end IFDEF M1

   end subroutine positivity_preserving_limiter

   subroutine positivity_preserving_limiter_KEEEEN()
     use mod_eos, only: small_rho
     implicit none
     double precision, dimension(ixI^S)   :: inv_lambda, epsD
     double precision, dimension(ixI^S)   :: theta, theta_tau
     double precision, dimension(ixI^S)   :: flow_dA, fhigh_dA
     double precision, parameter          :: eps_fac = 1.0d-6
     integer                              :: ix^D, ixR^L
     double precision, dimension(ixI^S) :: A, phi
     double precision :: kappa_bar


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
         {#IFDEF M1 
         if (iw .eq. nrad1_) exit
         } 
        {#IFDEF STAGGERED
           if ( (iw>=b1_) .and. (iw<=b^NC_) ) cycle
        }
        fC(ixC^S,iw,idims) = theta(ixC^S)*fC(ixC^S,iw,idims) &
             + (1.0d0-theta(ixC^S))*fC_low(ixC^S,iw,idims)
     end do ! Next iw
   
   {#IFDEF M1
   ! theta(ixI^S) = 1.0d0
   !theta_tmp(ixI^S) = 1.0d0
   {^D& ixRmin^D = ixOmin^D + kr(idims,^D); ixRmax^D = ixOmax^D + kr(idims,^D) \}
   
   {do ix^D = ixI^LIM^D \}
        epsD(ix^D)   = m1_E_atmo
   {end do^D&\}

   {^KSP& 
   !theta_tmp(ixI^S) = 1.0d0 !KEN KEN wird sowieso initialisiert
   !A(ixI^S) = 1.0d0

   if (slab) then
        !inv_lambda(ixI^S) = dxdim(idims)/(qdt)
        flow_dA(ixI^S)  = fC_low(ixI^S,nrad^KSP_,idims)
        fhigh_dA(ixI^S) = fC(ixI^S,nrad^KSP_,idims)
     else
        !inv_lambda(ixI^S) = mygeo%dvolume(ixI^S)/(qdt)
        select case (idims)
        {case (^D)
        flow_dA(ixI^S)  = fC_low(ixI^S,nrad^KSP_,idims)   * mygeo%surfaceC^D(ixI^S)  
        fhigh_dA(ixI^S) = fC(ixI^S,nrad^KSP_,idims) * mygeo%surfaceC^D(ixI^S)
        \}
        end select
     end if

     call get_theta(ixI^L,ixC^L,idims,epsD,inv_lambda(ixI^S),sCT%w%w(ixI^S,nrad^KSP_),&
                    flow_dA,fhigh_dA,theta)
   
          ! Compute A per interface
     
    !{do ix^D = ixC^LIM^D }
    !  ! Get kappa_bar at interface (geometric mean of neighbors)
    !  kappa_bar = (wradimpl(ix^D, kappa_a^KSP_) + wradimpl(ix^D, kappa_s^KSP_)) + &
    !               (wradimpl(ix^D+kr(idims,^D), kappa_a^KSP_) + &
    !                wradimpl(ix^D+kr(idims,^D), kappa_s^KSP_))       
    !   ! A = min(1, 4/(kappa_sum*dx))  -- matches correct_asymptotic_flux
    !                !KEN thinks 4 may be wrong
    !   A(ix^D) = min(1.0d0, 4.0d0 / (kappa_bar * dxdim(idims) + 1.0d-45))
    !   A(ix^D) = max(0.0d0, min(1.0d0, A(ix^D)))
!
    !  if(A(ix^D) .ne. A(ix^D)) then
    !    A(ix^D) = 1.0d0
    !  end if 
    !{end do}
    ! Compute A factor (vectorized over all cells in ixO^S)
   A(ixC^S) = MIN(1.0d0, dabs(4.0d0 / dabs(dxdim(idims)) / &
               (wradimpl(ixR^S,kappa_a^KSP_) + wradimpl(ixC^S,kappa_a^KSP_) + &
                wradimpl(ixR^S,kappa_s^KSP_) + wradimpl(ixC^S,kappa_s^KSP_) + 1.0d-40)))

   ! Handle NaN (branchless check: A /= A is true only for NaN)
   where(A(ixC^S) /= A(ixC^S)) A(ixC^S) = 1.0d0
    
   ! Handle overflow (already handled by MIN above, but explicit check)
   ! where(A(ixC^S) > 1.0d0) A(ixC^S) = 1.0d0
    
    ! Compute phi = (1 - theta) * A
    phi(ixC^S) = (1.0d0 - theta(ixC^S)) * A(ixC^S)
    
    ! Apply blending to ALL radiation variables
    fC(ixC^S,erad^KSP_,idims) = phi(ixC^S) * fC_low(ixC^S,erad^KSP_,idims) + &
                                 (1.0d0 - phi(ixC^S)) * fC(ixC^S,erad^KSP_,idims)
    
    fC(ixC^S,nrad^KSP_,idims) = phi(ixC^S) * fC_low(ixC^S,nrad^KSP_,idims) + &
                                 (1.0d0 - phi(ixC^S)) * fC(ixC^S,nrad^KSP_,idims)
    
    {^C&
    fC(ixC^S,frad^KSP^C_,idims) = phi(ixC^S) * fC_low(ixC^S,frad^KSP^C_,idims) + &
                                   (1.0d0 - phi(ixC^S)) * fC(ixC^S,frad^KSP^C_,idims)
    \}
    
    ! Update global theta for this species
    ! theta(ixC^S) = min(theta(ixC^S), theta_tmp(ixC^S)) !KEN KEN
   \} ! end KSP loop


     ! m1: correct fluxes form hll for diffusive lim
     ! KEN
     } 


   end subroutine positivity_preserving_limiter_KEEEEN

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

end subroutine finite_volume
