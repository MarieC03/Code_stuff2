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

{#IFDEF GLM
!=============================================================================
subroutine addsource_glmb(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

  ! Add divB related sources to w within ixO
  ! Split part
  include 'amrvacdef.f'

  integer, intent(in) :: ixI^L, ixO^L, iw^LIM
  double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
  double precision, intent(in) :: dx^D
  double precision, intent(inout) :: w(ixI^S,1:nw)
  !-----------------------------------------------------------------------------

  ! implicit update of psi variable
  !w(ixO^S,psi_) = exp(-qdt*cmax_global/eqpar(Cr_))*w(ixO^S,psi_)
  w(ixO^S,psi_) = exp(-eqpar(Cr_))*w(ixO^S,psi_)

end subroutine addsource_glmb
!=============================================================================
subroutine addsource_glma(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

  ! Add divB related sources to w within ixO
  ! Unsplit part
  use mod_metric, only: raise3
  include 'amrvacdef.f'

  integer, intent(in) :: ixI^L, ixO^L, iw^LIM
  double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
  double precision, intent(in) :: dx^D
  double precision, intent(inout) :: w(ixI^S,1:nw)
  ! .. local ..
  integer                                    :: idims
  double precision, dimension(ixG^T,1:ndir)  :: GradPsi
  double precision, dimension(ixI^S,1:ndir)  :: GradPsiU
  double precision,dimension(ixG^T)          :: Psi
  !-----------------------------------------------------------------------------

  Psi(ixI^S)=wCT(ixI^S,psi_)
  do idims=1,ndir
     if (idims .gt. ndim) then
        gradPsi(ixO^S,idims) = zero
     else
        select case(typegrad)
        case("central")
           call gradient(Psi,ixO^L,idims,gradPsi(ixG^T,idims))
        case("limited")
           call gradientS(Psi,ixO^L,idims,gradPsi(ixG^T,idims))
        case("upwind")
           call upwindGradientS(ixI^L,ixO^L,idims,wCT,x,Psi(ixI^S),gradPsi(ixI^S,idims))
        end select
     end if
  end do
  call raise3(ixI^L,ixO^L,myM,gradPsi(ixI^S,1:ndir),GradPsiU)

  ! B**j = B**j - qdt * alpha * gamma**ji del_i psi
  do idims=1, ndim
     w(ixO^S,b0_+idims) = w(ixO^S,b0_+idims) - qdt*GradPsiU(ixO^S,idims)*myM%alpha(ixO^S)
  end do

  call getaux(.true.,w,x,ixI^L,ixO^L,'addsource_glma')

end subroutine addsource_glma
}
!=============================================================================
subroutine get_divb(ixI^L,ixO^L,w,divb)

  ! Calculate div B within ixO

  include 'amrvacdef.f'

  integer, intent(in)                :: ixI^L, ixO^L
  double precision, intent(in)       :: w(ixI^S,1:nw)
  double precision, intent(out)      :: divb(ixI^S)
  ! .. local ..
  double precision                   :: bvec(ixI^S,1:ndir)

  integer                            :: ixC^L, idir, idim, ixJp^L, ic^D, ix^L
  integer                            :: ixKp^L, ixJpKp^L, ixJm^L, ixJmKm^L
  double precision                   :: divb_corner(ixI^S), sign
  double precision                   :: aux_vol(ixI^S)
  !-----------------------------------------------------------------------------

  bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)

  if (typeemf .eq. 'none') then
     ! We don't enforce divB to machine precision:
     select case(typediv)
     case("central")
        call divvector(bvec,ixI^L,ixO^L,divb)
     case("limited")
        call divvectorS(bvec,ixI^L,ixO^L,divb)
     end select
     
  else
     ! Use the FCT (larger stencil) formula for divB:

     ! For fct, we calculate the divB on the corners according to Toth (2000), 
     ! eq. (27) and average to the cell centers for output.

     {ixCmax^D=ixOmax^D;}
     {ixCmin^D=ixOmin^D-1;} ! Extend range by one


     ! Get the corner-centered divb:

     divb_corner(ixC^S) = zero
     aux_vol(ixC^S) = zero
     do idir = 1, ndim ! idir is the component of the field to consider (j)
        {do ic^DB=0,1\}
        {ix^L=ixC^L+ic^D;}

        select case(idir)
           {^D&   
                case(^D)
           sign = dble(ic^D*2 - 1)
           \}
        end select

        if (slab) then
           divb_corner(ixC^S) = divb_corner(ixC^S) &
                + sign * bvec(ix^S,idir)/dxlevel(idir)
        else
           divb_corner(ixC^S) = divb_corner(ixC^S) &
                + sign * mygeo%dvolume(ix^S) * bvec(ix^S,idir)/dxlevel(idir)
           aux_vol(ixC^S) = aux_vol(ixC^S) + mygeo%dvolume(ix^S)
        end if
        {end do\}
     end do

     divb_corner(ixC^S) = 2.0d0 * ndim * divb_corner(ixC^S) / aux_vol(ixC^S)

     if (slab) divb_corner(ixC^S) = divb_corner(ixC^S) / 2.0d0**(ndim-1)



     ! Now average back to the cell centers:

     divb(ixO^S) = zero
     {do ic^DB=-1,0\}
     {ixC^L=ixO^L+ic^D;}

     divb(ixO^S) = divb(ixO^S) &
          + divb_corner(ixC^S)

     {end do\}
     divb(ixO^S) = divb(ixO^S) / 2.0d0**(ndim)

  end if

end subroutine get_divb
!=============================================================================
