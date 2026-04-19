  !=============================================================================
  ! amrvacusr.t
  !=============================================================================
  INCLUDE:amrvacnul/speciallog.t
  INCLUDE:amrvacnul/specialbound.t
  INCLUDE:amrvacnul/specialsource.t
  INCLUDE:amrvacnul/specialimpl.t
  INCLUDE:amrvacnul/usrflags.t
  INCLUDE:amrvacnul/correctaux_usr.t
  !=============================================================================
  subroutine initglobaldata_usr

    include 'amrvacdef.f'
    !-----------------------------------------------------------------------------

!    Here you set global parameters, e.g. the adiabatic index:

!    eqpar(gamma_) = 5.0d0/3.0d0
    
  end subroutine initglobaldata_usr
  !=============================================================================
  subroutine initonegrid_usr(ixI^L,ixO^L,s)

    ! initialize one grid within ixO^L

    include 'amrvacdef.f'

    integer, intent(in) :: ixI^L, ixO^L
    type(state)         :: s
    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})


!      Set your initial conditions here
!      e.g.:
!      w(ixG^S,rho_) = 1.0d0
!      w(ixO^S,pp_)  = 1.0d0
      
!      w(ixO^S,u1_)  = 0.0d0
!      w(ixO^S,u2_)  = 0.0d0       
!      w(ixO^S,u3_)  = 0.0d0 
      
      
!      For divB-free fields, obtain B-fields like this:
!      {#IFNDEF STAGGERED
!      call b_from_vectorpotential(ixI^L,ixO^L,w,x)
!      }{#IFDEF STAGGERED
!      call b_from_vectorpotential(s%ws%ixG^L,ixI^L,ixO^L,ws,x)
!      call faces2centers(ixO^L,s)
!      }

!      Eventually convert to primitive variables:
!      call conserve(ixI^L,ixO^L,w,x,patchfalse)


    end associate
  end subroutine initonegrid_usr
  !=============================================================================
  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)

    ! initialize the vectorpotential on the corners
    ! used by b_from_vectorpotential()


    include 'amrvacdef.f'

    integer, intent(in)                :: ixI^L, ixC^L, idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    !-----------------------------------------------------------------------------

    ! Set your vectorpotential in direction idir here

    A(ixC^S) = zero


  end subroutine initvecpot_usr
  !=============================================================================
  ! amrvacusr.t
  !=============================================================================
