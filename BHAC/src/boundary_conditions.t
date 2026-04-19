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
subroutine bc_phys(iside,idims,time,s,ixB^L)

include 'amrvacdef.f'

integer, intent(in)          :: iside, idims, ixB^L
double precision, intent(in) :: time
type(state), intent(inout)   :: s

integer :: iw, iB, ix^D, ixI^L, ixD^L, ixDB^L, ixIp^L, ixG^L
{#IFDEF STAGGERED
integer :: idir, is
integer :: ixIs^L,hxI^L,jxI^L
double precision :: Q(ixG^T),Qp(ixG^T) 
logical, save    :: zerodiv=.true.
double precision, dimension(:^D&), pointer :: surf
}
integer, parameter :: dixPolefix=2, dixHotaka=1
logical            :: is_coarse(ndim,2)
integer         :: iksp
!-----------------------------------------------------------------------------
associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})
ixG^L=s%w%ixG^L;

! Test if we have a coarse neighbor:
call is_neighbor_coarse(s,is_coarse)

select case (idims)
{case (^D)
if (iside==2) then
       
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! maximal boundary
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iB=ismax^D
      ixImin^DD=ixBmax^D+1-dixB^D%ixImin^DD=ixBmin^DD;
      ixImax^DD=ixBmax^DD;

      !!!!!!! Make tangential ranges greater for staggered components
      !!!!!!! Add fill the normal component

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         call mpistop('Turn off primitiveB for every boundary')
         ixDmin^DD=ixBmin^D+dixB^D%ixDmin^DD=ixBmin^DD;
         ixDmax^DD=ixBmax^D-dixB^D%ixDmax^DD=ixBmax^DD;
         call primitive(ixG^L,ixD^L,w,x)
      end if
      
      ! cont/symm/asymm types
      do iw=1,nw
{#IFDEF STAGGERED ! If the variable has no staggered counterpart
         if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then
}
         select case (typeB(iw,iB))
         case ("symm")
            w(ixI^S,iw) = w(ixImin^D-1:ixImin^D-dixB:-1^D%ixI^S,iw)
         case ("asymm")
            w(ixI^S,iw) =-w(ixImin^D-1:ixImin^D-dixB:-1^D%ixI^S,iw)
         case ("cont")
            do ix^D=ixImin^D,ixImax^D
               w(ix^D^D%ixI^S,iw) = w(ixImin^D-1^D%ixI^S,iw)
            end do
         case("noinflow")
            ! only need velocities and momentum to be noinflow,
            ! all the other is const actually
            if (iw == u^D_) then
            !if (iw == u^D_ .or. iw == s^D_) then
              do ix^D=ixImin^D,ixImax^D
                  w(ix^D^D%ixI^S,iw) = max(w(ixImin^D-1^D%ixI^S,iw),zero)
              end do
               {#IFDEF M1
               else if({^KSP& (iw == frad^KSP^D_) .or. &\} 
                  .false.) then
                 do ix^D=ixImin^D,ixImax^D
                   w(ix^D^D%ixI^S,iw) = max(w(ixImin^D-1^D%ixI^S,iw),zero)
                 end do
               }
               else
                  do ix^D=ixImin^D,ixImax^D
                        w(ix^D^D%ixI^S,iw) = w(ixImin^D-1^D%ixI^S,iw)
                  end do
               end if 
    
               

         case("limitinflow")
            if (iw==1+^D)then
              do ix^D=ixImin^D,ixImax^D
                  w(ix^D^D%ixI^S,iw) = max(w(ixImin^D-1^D%ixI^S,iw), &
                                           w(ixImin^D-1^D%ixI^S,iw)*ratebdflux)
              end do
            else
              do ix^D=ixImin^D,ixImax^D
                  w(ix^D^D%ixI^S,iw) = w(ixImin^D-1^D%ixI^S,iw)
              end do
           end if
        case ("polefix")
           ! First overwrite cells in domain
           ixIpmin^DD=ixImin^D-dixPolefix^D%ixIpmin^DD=ixImin^DD;
           ixIpmax^DD=ixImin^D-1^D%ixIpmax^DD=ixImax^DD;
           do ix^D=ixIpmin^D,ixIpmax^D
              w(ix^D^D%ixIp^S,iw) = w(ixIpmin^D-1^D%ixIp^S,iw)
           end do
           ! Now apply symm boundary
            w(ixI^S,iw) = w(ixImin^D-1:ixImin^D-dixB:-1^D%ixI^S,iw)
        case ("apolefix")
           ! First interpolate cells in domain down to zero at pole
           ixIpmin^DD=ixImin^D-dixPolefix^D%ixIpmin^DD=ixImin^DD;
           ixIpmax^DD=ixImax^D-1^D%ixIpmax^DD=ixImax^DD;
           do ix^D=ixIpmin^D,ixIpmax^D
              w(ix^D^D%ixIp^S,iw) = w(ixIpmin^D-1^D%ixIp^S,iw) &
                   * (xprobmax^D-x(ix^D^D%ixIp^S,^D)) &
                   / (xprobmax^D-x(ixIpmin^D-1^D%ixIp^S,^D))
           end do
        case ("hotaka")
           ! First overwrite cells in domain with zero
           ixIpmin^DD=ixImin^D-dixHotaka^D%ixIpmin^DD=ixImin^DD;
           ixIpmax^DD=ixImin^D-1^D%ixIpmax^DD=ixImax^DD;
           do ix^D=ixIpmin^D,ixIpmax^D
              w(ix^D^D%ixIp^S,iw) = zero
           end do
           ! Now apply symm boundary
            w(ixI^S,iw) = w(ixImin^D-1:ixImin^D-dixB:-1^D%ixI^S,iw)
        case ("ahotaka")
            ! First overwrite cells in domain with zero
           ixIpmin^DD=ixImin^D-dixHotaka^D%ixIpmin^DD=ixImin^DD;
           ixIpmax^DD=ixImax^D-1^D%ixIpmax^DD=ixImax^DD;
           do ix^D=ixIpmin^D,ixIpmax^D
              w(ix^D^D%ixIp^S,iw) = zero
           end do
           ! Now apply asymm boundary
            w(ixI^S,iw) =-w(ixImin^D-1:ixImin^D-dixB:-1^D%ixI^S,iw)
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
            !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixI^S,iw) = - w(ixI^S,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select

{#IFDEF STAGGERED
      else !If there is a staggered counterpart
      ! At this stage, extrapolation is applied only to the tangential components
        idir=mod(iws(iw)-1,ndir)+1
        
        if (idir.ne.^D) then
        ixIsmax^DD=ixImax^DD;
        ixIsmin^DD=ixImin^DD-kr(^DD,idir);
        select case(typeB(iw,iB))
        case ("symm")
          ws(ixIs^S,iws(iw)) = ws(ixIsmin^D-1:ixIsmin^D-dixB:-1^D%ixIs^S,iws(iw))
        case ("asymm")
          ws(ixIs^S,iws(iw)) =-ws(ixIsmin^D-1:ixIsmin^D-dixB:-1^D%ixIs^S,iws(iw))
        case ("cont")
          do ix^D=ixIsmin^D,ixIsmax^D
             ws(ix^D^D%ixIs^S,iws(iw)) = ws(ixIsmin^D-1^D%ixIs^S,iws(iw))
          end do
        case ("periodic")
        case ("zero")
           ws(ixIs^S,iws(iw)) = 0.0d0
        case ("specialT")
           call specialbound_usr(time,ixG^L,ixI^L,-iB,s)
        case ("special")
           ! skip it here, do AFTER all normal type boundaries are set
        case default
          write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
             "for variable iw=",iw," and side iB=",iB
        end select

{#IFNDEF D1              

      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin^DD=ixImin^DD-kr(^DD,idir);
            ixIsmax^DD=ixImax^DD;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  {^DD&
                  case (^DD)
                     ixIsmax^DD=min(ixGlo^DD+dixB-1,ixIsmax^DD);
                     surf=>s%geo%surfaceC^DD
                  \}! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  {^DD&
                  case (^DD)
                     ixIsmin^DD=max(ixGhi^DD-dixB,ixIsmin^DD);
                     surf=>s%geo%surfaceC^DD
                  \}! I am hating LASY just right now.
               end select

            end if

            select case(typeB(iw,iB))
            case ("cont")
               do ix^D=ixIsmin^D,ixIsmax^D
                  where ((surf(ixIsmin^D-1^D%ixIs^S) + surf(ixIsmin^D-2^D%ixIs^S)) .gt. smalldouble )
                     ws(ix^D^D%ixIs^S,iws(iw)) = (ws(ixIsmin^D-1^D%ixIs^S,iws(iw))*surf(ixIsmin^D-1^D%ixIs^S) &
                          + ws(ixIsmin^D-2^D%ixIs^S,iws(iw))*surf(ixIsmin^D-2^D%ixIs^S) ) / (surf(ixIsmin^D-1^D%ixIs^S)+surf(ixIsmin^D-2^D%ixIs^S))
                  elsewhere
                     ws(ix^D^D%ixIs^S,iws(iw)) = zero
                  end where
               end do

            case default
               ! Only implemented cont type, symm and asymm will also work if problem respects symmetry
            end select

         end if
      end do
}
   end if
end if
}
      end do
{#IFDEF STAGGERED
      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
      do iw=1,nws
         
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.^D) then
          ixIs^L=ixI^L;
          hxI^L=ixI^L-dixB*kr(^DD,^D);

          ! Calculate divergence and partial divergence
          Q=zero
          if (.not.zerodiv) call div_staggered(ixG^L,hxI^L,s,Q(ixG^T))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIs^S,iw)=zero
            do ix^D=0,dixB-1
              call div_staggered(ixG^L,ixI^L,s,Qp(ixG^T))
              ws(ixIsmin^D+ix^D^D%ixIs^S,iw)=&
              (Q(hxImax^D-ix^D^D%hxI^S)*s%geo%dvolume(hxImax^D-ix^D^D%hxI^S)&
             -Qp(ixImin^D+ix^D^D%ixI^S)*s%geo%dvolume(ixImin^D+ix^D^D%ixI^S))&
              /s%geo%surfaceC^D(ixIsmin^D+ix^D^D%ixIs^S)
            end do
          case("asymm")
            ws(ixIs^S,iw)=zero
            do ix^D=0,dixB-1
              call div_staggered(ixG^L,ixI^L,s,Qp(ixG^T))
              ws(ixIsmin^D+ix^D^D%ixIs^S,iw)=&
             (-Q(hxImax^D-ix^D^D%hxI^S)*s%geo%dvolume(hxImax^D-ix^D^D%hxI^S)&
             -Qp(ixImin^D+ix^D^D%ixI^S)*s%geo%dvolume(ixImin^D+ix^D^D%ixI^S))&
              /s%geo%surfaceC^D(ixIsmin^D+ix^D^D%ixIs^S)
            end do
          case("cont", "zero")
            ws(ixIs^S,iw)=zero
            do ix^D=0,dixB-1
              call div_staggered(ixG^L,ixI^L,s,Qp(ixG^T))
              ws(ixIsmin^D+ix^D^D%ixIs^S,iw)=&
              (Q(hxImax^D^D%hxI^S)*s%geo%dvolume(hxImax^D^D%hxI^S)&
             -Qp(ixImin^D+ix^D^D%ixI^S)*s%geo%dvolume(ixImin^D+ix^D^D%ixI^S))&
              /s%geo%surfaceC^D(ixIsmin^D+ix^D^D%ixIs^S)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixI^L,s)

}

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         call mpistop('Turn off primitiveB for every boundary')
         ixDBmin^DD=ixBmin^D+dixB^D%ixDBmin^DD=ixBmin^DD;
         ixDBmax^DD=ixBmax^DD;
         call conserve(ixG^L,ixDB^L,w,x,patchfalse)
      end if
      
      
   else
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! minimal boundary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      iB=ismin^D
      ixImin^DD=ixBmin^DD;
      ixImax^DD=ixBmin^D-1+dixB^D%ixImax^DD=ixBmax^DD;
      

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then         
         call mpistop('Turn off primitiveB for every boundary')
         ixDmin^DD=ixBmin^D+dixB^D%ixDmin^DD=ixBmin^DD;
         ixDmax^DD=ixBmax^D-dixB^D%ixDmax^DD=ixBmax^DD;
         call primitive(ixG^L,ixD^L,w,x)
      end if
 
      ! cont/symm/asymm types
      do iw=1,nw
         
{#IFDEF STAGGERED ! If the variable has no staggered counterpart
if ((iws(iw).lt.1).or.(iws(iw).gt.nws)) then
   }
   select case (typeB(iw,iB))
   case ("symm")
      w(ixI^S,iw) = w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
   case ("asymm")
      w(ixI^S,iw) =-w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
   case ("cont")
      do ix^D=ixImin^D,ixImax^D
         w(ix^D^D%ixI^S,iw) = w(ixImax^D+1^D%ixI^S,iw)
      end do
   case("noinflow")
      ! only need velocities and momentum to be noinflow,
      ! all the other is const actually
      if (iw == u^D_) then
      !if (iw == u^D_ .or. iw == s^D_) then
        do ix^D=ixImin^D,ixImax^D
            w(ix^D^D%ixI^S,iw) = min(w(ixImax^D+1^D%ixI^S,iw),zero)
        end do
      {#IFDEF M1
               else if({^KSP& (iw == frad^KSP^D_) .or. &\} 
                  .false.) then
                 do ix^D=ixImin^D,ixImax^D
                   w(ix^D^D%ixI^S,iw) = min(w(ixImax^D+1^D%ixI^S,iw),zero)
                 end do
      }
               else
                  do ix^D=ixImin^D,ixImax^D
                        w(ix^D^D%ixI^S,iw) = w(ixImax^D+1^D%ixI^S,iw)
                  end do
               end if 
   case("limitinflow")
      if (iw==1+^D)then
         do ix^D=ixImin^D,ixImax^D
            w(ix^D^D%ixI^S,iw) = min(w(ixImax^D+1^D%ixI^S,iw), &
                 w(ixImax^D+1^D%ixI^S,iw)*ratebdflux)
         end do
      else
         do ix^D=ixImin^D,ixImax^D
            w(ix^D^D%ixI^S,iw) = w(ixImax^D+1^D%ixI^S,iw)
         end do
      end if
   case ("polefix")
      ! First overwrite cells in domain
      ixIpmin^DD=ixImax^D+1^D%ixIpmin^DD=ixImin^DD;
      ixIpmax^DD=ixImax^D+dixPolefix^D%ixIpmax^DD=ixImax^DD;
      do ix^D=ixIpmin^D,ixIpmax^D
         w(ix^D^D%ixIp^S,iw) = w(ixIpmax^D+1^D%ixIp^S,iw)
      end do
      ! Now apply symm boundary
      w(ixI^S,iw) = w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
   case ("apolefix")
      ! First interpolate cells in domain down to zero at pole
      ixIpmin^DD=ixImax^D+1^D%ixIpmin^DD=ixImin^DD;
      ixIpmax^DD=ixImax^D+dixPolefix^D%ixIpmax^DD=ixImax^DD;
      do ix^D=ixIpmin^D,ixIpmax^D
         w(ix^D^D%ixIp^S,iw) = w(ixIpmax^D+1^D%ixIp^S,iw) &
              * (x(ix^D^D%ixIp^S,^D))/x(ixIpmax^D+1^D%ixIp^S,^D)
      end do
      ! Now apply asymm boundary
      w(ixI^S,iw) =-w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
   case ("hotaka")
      ! First overwrite cells in domain with zero
      ixIpmin^DD=ixImax^D+1^D%ixIpmin^DD=ixImin^DD;
      ixIpmax^DD=ixImax^D+dixHotaka^D%ixIpmax^DD=ixImax^DD;
      do ix^D=ixIpmin^D,ixIpmax^D
         w(ix^D^D%ixIp^S,iw) = zero
      end do
      ! Now apply symm boundary
      w(ixI^S,iw) = w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
   case ("ahotaka")
      ! First overwrite cells in domain with zero
      ixIpmin^DD=ixImax^D+1^D%ixIpmin^DD=ixImin^DD;
      ixIpmax^DD=ixImax^D+dixHotaka^D%ixIpmax^DD=ixImax^DD;
      do ix^D=ixIpmin^D,ixIpmax^D
         w(ix^D^D%ixIp^S,iw) = zero
      end do
      ! Now apply asymm boundary
      w(ixI^S,iw) =-w(ixImax^D+dixB:ixImax^D+1:-1^D%ixI^S,iw)
   case ("special")
      ! skip it here, do AFTER all normal type boundaries are set
   case ("aperiodic")
      !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
      w(ixI^S,iw) = - w(ixI^S,iw)
   case ("periodic")
      !            call mpistop("periodic bc info should come from neighbors")
   case default
      write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
           "for variable iw=",iw," and side iB=",iB
   end select

{#IFDEF STAGGERED
else !If there is a staggered counterpart
   ! At this stage, extrapolation is applied only to the tangential components
   idir=mod(iws(iw)-1,ndir)+1 ! vector direction
   
   if (idir.ne.^D) then
      ixIsmax^DD=ixImax^DD;
      ixIsmin^DD=ixImin^DD-kr(^DD,idir);
      select case(typeB(iw,iB))
      case ("symm")
         ws(ixIs^S,iws(iw)) = ws(ixIsmax^D+dixB:ixIsmax^D+1:-1^D%ixIs^S,iws(iw))
      case ("asymm")
         ws(ixIs^S,iws(iw)) =-ws(ixIsmax^D+dixB:ixIsmax^D+1:-1^D%ixIs^S,iws(iw))
      case ("cont")
         do ix^D=ixIsmin^D,ixIsmax^D
            ws(ix^D^D%ixIs^S,iws(iw)) = ws(ixIsmax^D+1^D%ixIs^S,iws(iw))
         end do
      case("periodic")
      case ("zero")
         ws(ixIs^S,iws(iw)) = 0.0d0
      case ("specialT")
      call specialbound_usr(time,ixG^L,ixI^L,-iB,s)
      case ("special")
         ! skip it here, do AFTER all normal type boundaries are set
      case default
         write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
              "for variable iw=",iw," and side iB=",iB
      end select

{#IFNDEF D1

      ! Special treatment when we have coarse neighbors in tangential dir:
      do is=1,2
         if (is_coarse(idir,is)) then

            ixIsmin^DD=ixImin^DD-kr(^DD,idir);
            ixIsmax^DD=ixImax^DD;

            if (is.eq.1) then ! coarse neighbor at beginning of block

               select case(idir)
                  {^DD&
                  case (^DD)
                     ixIsmax^DD=min(ixGlo^DD+dixB-1,ixIsmax^DD);
                     surf=>s%geo%surfaceC^DD
                  \}! I am hating LASY just right now.
               end select

            else ! coarse neighbor at end of block

               select case(idir)
                  {^DD&
                  case (^DD)
                     ixIsmin^DD=max(ixGhi^DD-dixB,ixIsmin^DD);
                     surf=>s%geo%surfaceC^DD
                  \}! I am hating LASY just right now.
               end select

            end if


            select case(typeB(iw,iB))
            case ("cont")
               do ix^D=ixIsmin^D,ixIsmax^D
                  where ( (surf(ixIsmax^D+1^D%ixIs^S) + surf(ixIsmax^D+2^D%ixIs^S)) .gt. smalldouble )
                     ws(ix^D^D%ixIs^S,iws(iw)) = (ws(ixIsmax^D+1^D%ixIs^S,iws(iw))*surf(ixIsmax^D+1^D%ixIs^S) &
                          + ws(ixIsmax^D+2^D%ixIs^S,iws(iw))*surf(ixIsmax^D+2^D%ixIs^S) ) / (surf(ixIsmax^D+1^D%ixIs^S)+surf(ixIsmax^D+2^D%ixIs^S))
                  elsewhere
                     ws(ix^D^D%ixIs^S,iws(iw)) = zero
                  end where
               end do

            case default
               ! Only implemented cont type, symm and asymm will also work if problem respects symmetry
            end select

         end if
      end do
}
   end if

end if
}
end do

{#IFDEF STAGGERED
      ! Now that the tangential components are set,
      ! fill the normal components using a prescription for the divergence.
      ! This prescription is given by the typeB for the normal component.
do iw=1,nws
   
        ! When using more than one staggered vector field
        idir=mod(iw-1,ndir)+1
        ! Consider only normal direction
        if (idir.eq.^D) then
          ixIs^L=ixI^L-kr(^DD,^D);
          jxI^L=ixI^L+dixB*kr(^DD,^D);

          ! Calculate divergence and partial divergence
          Q=0
          if (.not.zerodiv) call div_staggered(ixG^L,jxI^L,s,Q(ixG^T))
          select case(typeB(iw+b0_,iB))
          case("symm")
            ws(ixIs^S,iw)=zero
            do ix^D=0,dixB-1
              call div_staggered(ixG^L,ixI^L,s,Qp(ixG^T))
              ws(ixIsmax^D-ix^D^D%ixIs^S,iw)=&
             -(Q(jxImin^D+ix^D^D%jxI^S)*s%geo%dvolume(jxImin^D+ix^D^D%jxI^S)&
             -Qp(ixImax^D-ix^D^D%ixI^S)*s%geo%dvolume(ixImax^D-ix^D^D%ixI^S))&
             /s%geo%surfaceC^D(ixIsmax^D-ix^D^D%ixIs^S)
            end do
          case("asymm")
             ws(ixIs^S,iw)=zero
            do ix^D=0,dixB-1
              call div_staggered(ixG^L,ixI^L,s,Qp(ixG^T))
              ws(ixIsmax^D-ix^D^D%ixIs^S,iw)=&
            -(-Q(jxImin^D+ix^D^D%jxI^S)*s%geo%dvolume(jxImin^D+ix^D^D%jxI^S)&
             -Qp(ixImax^D-ix^D^D%ixI^S)*s%geo%dvolume(ixImax^D-ix^D^D%ixI^S))&
             /s%geo%surfaceC^D(ixIsmax^D-ix^D^D%ixIs^S)
            end do
          case("cont", "zero")
            ws(ixIs^S,iw)=zero
            do ix^D=0,dixB-1
              call div_staggered(ixG^L,ixI^L,s,Qp(ixG^T))
              ws(ixIsmax^D-ix^D^D%ixIs^S,iw)=&
             -(Q(jxImin^D^D%jxI^S)*s%geo%dvolume(jxImin^D^D%jxI^S)&
             -Qp(ixImax^D-ix^D^D%ixI^S)*s%geo%dvolume(ixImax^D-ix^D^D%ixI^S))&
              /s%geo%surfaceC^D(ixIsmax^D-ix^D^D%ixIs^S)
            end do
          case("periodic")
          end select
        end if
      end do

      ! Fill cell averages
      call faces2centers(ixI^L,s)
}
      
      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         call mpistop('Turn off primitiveB for every boundary')
         ixDBmin^DD=ixBmin^DD;
         ixDBmax^DD=ixBmax^D-dixB^D%ixDBmax^DD=ixBmax^DD;
         call conserve(ixG^L,ixDB^L,w,x,patchfalse)
      end if
      
   end if \}
end select

!{#IFDEF STAGGERED
!   call faces2centers(ixI^L,s)
!}

! do special case AFTER all normal cases are set
if (any(typeB(1:nwflux+nwaux,iB)=="special")) then
   call specialbound_usr(time,ixG^L,ixI^L,iB,s)
end if


   end associate
end subroutine bc_phys
!=============================================================================
subroutine getintbc(time,psuse)

include 'amrvacdef.f'

double precision, intent(in)              :: time
type(state), dimension(ngridshi)          :: psuse

! .. local ..
integer :: iigrid, igrid, ixO^L,level
!----------------------------------------------------------------------------
ixO^L=ixG^LL^LSUBdixB;

!$OMP PARALLEL DO PRIVATE(igrid,level)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call set_tmpGlobals(igrid)
   level = node(plevel_,igrid)
   call bc_int(level,time,ixG^LL,ixO^L,psuse(igrid)%w%w,px(igrid)%x)
end do
!$OMP END PARALLEL DO

      
end subroutine getintbc
!=============================================================================
subroutine is_neighbor_coarse(s,is_coarse)

include 'amrvacdef.f'

type(state), intent(in)          :: s
logical, intent(out)             :: is_coarse(ndim,2)
! .. local ..
integer                          :: i^D, i
!-----------------------------------------------------------------------------

is_coarse=.false.
! If we are the coarse version of a buffer
! dont consider neighbor as coarser
if (s%is_coarse .eqv. .true.) return 

   {do i^DB=-1,1\}
      if (({i^D.eq.0|.and.}).or.({abs(i^D)|+}).ne.1) cycle
      if (neighbor_type(i^D,s%igrid).eq.2) then
         {if (i^D.eq.-1) is_coarse(^D,1) = .true.\}
         {if (i^D.eq.+1) is_coarse(^D,2) = .true.\}
      end if
   {end do\}

end subroutine is_neighbor_coarse
!=============================================================================
