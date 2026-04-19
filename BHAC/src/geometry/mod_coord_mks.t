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
  ! Metric components for modified Kerr-Schild coordinates -coord=mks
  !
  ! Check modKSorig.nb for the definitions.
  ! Oliver Porth
  ! 04.02.2017
  !
  !=============================================================================
  
  !-----------------------------------------
  ! Define constants specific to coordinates
  !-----------------------------------------

  character*20, parameter              :: coord="mks"
  logical, save                        :: init_from_g4 = .false.  ! We don't provide the four-metric

  ! Migrating to have coordinate-specific constants here.
  ! a_ and m_ however are in the physics part
  ! so put the MKS parameters, e.g. R0_ and h_ here:
  integer, parameter                   :: R0_=1, h_=2
  integer, parameter                   :: ncoordpar = 2
  double precision, save               :: coordpar(ncoordpar)
  ! not yet implemented...

contains
!=============================================================================
  subroutine init_coord

    get_gammainv_component_analytic => mks_get_gammainv_component_analytic
    get_gammainv_component => get_gammainv_component_analytic
    LRNFu => MKS_LRNFu
    u4CoordToKS => MKS_u4CoordToKS
    u4CoordToBL => MKS_u4CoordToBL
    
  end subroutine init_coord
  !=============================================================================
  double precision function s(rks)

    include 'amrvacdef.f'
    double precision                               :: rks
    ! .. local ..
    !-----------------------------------------------------------------------------

    s = log( rks - coordpar(R0_) )

  end function s
  !=============================================================================
  double precision function vartheta(theta)

    use mod_rtsafe2D, only: MAXIT, rt1D
    include 'amrvacdef.f'
    double precision                               :: theta
    ! .. local ..
    double precision                               :: rtmp
    integer                                        :: ierror
    !-----------------------------------------------------------------------------

    vartheta = theta
    if (coordpar(h_) .eq. zero) then 
       return
    else

       MAXIT = 30
       call rt1D(func_vartheta,vartheta,absaccnr,ierror)

       if (ierror .ne. 0) call mpistop('vartheta: Does not converge')

    end if
  contains
    !=============================================================================
    subroutine func_vartheta(myvartheta,f,df)

      double precision, intent(in)      :: myvartheta
      double precision, intent(out)     :: f, df
      ! .. local ..
      !-----------------------------------------------------------------------------

      f  = myvartheta+(coordpar(h_)*Sin(2.0d0*myvartheta))/2.0d0 - theta
      df = 1.0d0 + coordpar(h_) * Cos(2.0d0*myvartheta)

    end subroutine func_vartheta
    !=============================================================================
  end function vartheta
  !=============================================================================
  double precision function rks(s)

    double precision                               :: s
    ! .. local ..
    !-----------------------------------------------------------------------------

    rks = coordpar(R0_) + exp(s)

  end function rks
  !=============================================================================
  double precision function thetaks(vartheta)

    double precision                               :: vartheta

    ! .. local ..
    !-----------------------------------------------------------------------------

    thetaks = vartheta+(coordpar(h_)*Sin(2.0d0*vartheta))/2.0d0

  end function thetaks
  !=============================================================================
  subroutine get_sqrtgamma_analytic(x^D,sqrtgamma,is_analytic,w_pt)

    include 'amrvacdef.f'

    ! You can specify sqrtgamma here.
    ! If you don't want to, set is_analytic = .false.
    ! In that case, get_sqrtgamma() will calculate from the metric,
    ! which can be a major slowdown.
    ! init_metric() will set the flag sqrtgamma_is_analytic from the return-value
    ! of this subroutine.  

    double precision, intent(in)                     :: x^D
    double precision, intent(out)                    :: sqrtgamma
    logical, optional                                :: is_analytic
    double precision, intent(in), optional           :: w_pt(1:nw)

    ! .. local ..
    {^IFZIN
    double precision                                 :: mythetaks}
    !-----------------------------------------------------------------------------

    if(present(is_analytic)) is_analytic = .true.

    {^IFZIN
    mythetaks = thetaks(x^Z)
    }
    
    {^IFZIN
    sqrtgamma = sqrt(Exp(2.0d0*x1)*((coordpar(R0_)+  &
         Exp(x1))*(coordpar(R0_)+2.0d0*eqpar(m_)+Exp(x1))+  &
         eqpar(a_)**2*Cos(mythetaks)**2)*(1.0d0+  &
         coordpar(h_)*Cos(2.0d0*x^Z))**2*Sin(mythetaks)**2*((coordpar(R0_)  &
         +Exp(x1))**2+eqpar(a_)**2-1.0d0*eqpar(a_)**2*Sin(mythetaks)**2))
    }{^IFZOUT
    sqrtgamma = sqrt((-1.0d0+  &
         coordpar(h_))**2*Exp(2.0d0*x1)*(coordpar(R0_)+  &
         Exp(x1))**3*(coordpar(R0_)+2.0d0*eqpar(m_)+Exp(x1)))
    }


  end subroutine get_sqrtgamma_analytic
  !=============================================================================
  subroutine mks_get_gammainv_component_analytic(iin,jin,x^D,ginv,iszero,dginvdk_iszero,dginvdk,kdir,w_pt)

    include 'amrvacdef.f'
    
    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: ginv
    logical, optional, intent(out)           :: iszero, dginvdk_iszero
    double precision, optional, intent(out)  :: dginvdk
    double precision, optional, intent(in)   :: w_pt(1:nw)

    ! .. local ..
    double precision                         :: myrks, mythetaks, rho2
    double precision                         :: d1rks, d2thks, d1rho2, d2rho2
    integer                                  :: i, j
    !-----------------------------------------------------------------------------
    if (iin .gt. jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    ! Auxiliary quantities
    myrks = rks(x1)
    {^IFZIN
    mythetaks = thetaks(x^Z)
    }
    {^IFZOUT
    mythetaks = dpi/2.0d0
    }
    rho2 = myrks**2 + (eqpar(a_) * cos(mythetaks))**2

    ! Depending on the direction of derivation,
    ! compute derivatives of the auxiliary quantities
    if (present(kdir)) then
       if (kdir .eq. 1) then
	  d1rks = exp(x1)
          d1rho2 = 2.d0 * d1rks * myrks
       else if (kdir .eq. 2) then
	  {^IFZIN
	  d2thks = 1.d0 + eqpar(h_) * cos(2.d0*x^Z)
	  }{^IFZOUT
	  d2thks = 0.d0
	  }
	  d2rho2 = -2.d0 * eqpar(a_)**2 * cos(mythetaks) * sin(mythetaks) * d2thks
       end if
    end if

    ! Get components of gamma^ij
    if (i .eq. 1 .and. j .eq. 1) then ! rr-component:
              
       ginv = 2.0d0*Exp(-2.0d0*x1)*(-((eqpar(m_)*myrks)/(2*eqpar(m_)*myrks  &
            +rho2))+(eqpar(a_)**2+myrks**2)/(eqpar(a_)**2+  &
            2*myrks**2+eqpar(a_)**2*Cos(2*mythetaks)))
       if (present(iszero)) iszero = .false.

       ! Derivatives of gamma^rr
       if (present(kdir)) then
          if (kdir .eq. 1) then
             if (present(dginvdk)) &
                dginvdk = - 2.d0*ginv &
                          + 2.d0*exp(-x1) * (eqpar(m_)*(2.d0*myrks**2-rho2)/(2.d0*eqpar(m_)*myrks+rho2)**2 &
                          + 2.d0*myrks*eqpar(a_)**2*(cos(2.d0*mythetaks)-1.d0) &
                            /(eqpar(a_)**2+2.d0*myrks**2+eqpar(a_)**2*cos(2.d0*mythetaks))**2)
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
          else if (kdir .eq. ^Z) then
	     {^IFZIN
             if (present(dginvdk)) &
                dginvdk = 2.d0*exp(-x1) * (eqpar(m_)*myrks*d2rho2/(2.d0*eqpar(m_)*myrks+rho2)**2 &
                          + 2.d0*(eqpar(a_)**2+myrks**2)*eqpar(a_)**2*sin(2.d0*mythetaks)*d2thks &
                            /(eqpar(a_)**2+2.d0*myrks**2+eqpar(a_)**2*cos(2.d0*mythetaks))**2)
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             }{^IFZOUT
             if (present(dginvdk)) dginvdk = 0.d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             }
          else
             if (present(dginvdk)) dginvdk = 0.0d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
          end if
       end if
       
    else if (i .eq. 1 .and. j .eq. ^PHI) then ! r-phi-component:
       
       ginv = (2.0d0*eqpar(a_)*Exp(-1.0d0*x1))/(eqpar(a_)**2+  &
            2*myrks**2+eqpar(a_)**2*Cos(2*mythetaks))
       if (present(iszero)) iszero = .false.

       ! Derivatives of gamma^rphi
       if (present(kdir)) then
          if (kdir .eq. 1) then
             if (present(dginvdk)) &
                dginvdk = - ginv * (1.d0 + 4.d0*myrks*exp(x1) &
                            /(eqpar(a_)**2+2.d0*myrks**2+eqpar(a_)**2*cos(2.d0*mythetaks)))
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
          else if (kdir .eq. ^Z) then
	     {^IFZIN
             if (present(dginvdk)) &
                dginvdk = ginv * 2.d0*eqpar(a_)**2*sin(2.d0*mythetaks)*d2thks/(eqpar(a_)**2 &
                            +2.d0*myrks**2+eqpar(a_)**2*cos(2.d0*mythetaks))
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             }{^IFZOUT
             if (present(dginvdk)) dginvdk = 0.d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             }
          else
             if (present(dginvdk)) dginvdk = 0.0d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
          end if
       end if
       
    else if (i .eq. ^Z .and. j .eq. ^Z) then ! theta-theta-component:
       
       {^IFZIN
       ginv = 1/(rho2*(1+coordpar(h_)*Cos(2.0d0*x^Z))**2)
       }{^IFZOUT
       ! equatorial plane
       ginv = 1/((1-1.0d0*coordpar(h_))**2*rho2)
       }
       if (present(iszero)) iszero = .false.

       ! Derivatives of gamma^thetatheta
       if (present(kdir)) then
          if (kdir .eq. 1) then
             {^IFZIN
             if (present(dginvdk)) &
                dginvdk = - ginv**2 * (1.d0+eqpar(h_)*cos(2.d0*x^Z))**2 * d1rho2
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             }{^IFZOUT
             if (present(dginvdk)) &
                dginvdk = - ginv**2 * (1.d0-eqpar(h_)**2) * d1rho2
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             }
          else if (kdir .eq. ^Z) then
	     {^IFZIN
             if (present(dginvdk)) &
                dginvdk = - ginv**2 * (d2rho2 * (1.d0+eqpar(h_)*cos(2.d0*x^Z))**2 &
                           - 2.d0*rho2*eqpar(h_)**2*sin(2.d0*x^Z))
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             }{^IFZOUT
             if (present(dginvdk)) dginvdk = 0.d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             }
          else
             if (present(dginvdk)) dginvdk = 0.0d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
          end if
       end if
       
    else if (i .eq. ^PHI .and. j .eq. ^PHI) then ! phi-phi component
       
       ginv = 1/sin(mythetaks)**2/(eqpar(a_)**2+myrks**2-  &
            1.0d0*eqpar(a_)**2*Sin(mythetaks)**2)
       if (present(iszero)) iszero = .false.

       ! Derivatives of gamma^phiphi
       if (present(kdir)) then
          if (kdir .eq. 1) then
             if (present(dginvdk)) &
                dginvdk = - ginv**2 * sin(mythetaks)**2 * 2.d0*exp(x1)*myrks
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
          else if (kdir .eq. ^Z) then
	     {^IFZIN
             if (present(dginvdk)) &
                dginvdk = - ginv**2 * d2thks * sin(2.d0*mythetaks) & 
                          * (eqpar(a_)**2+myrks**2-2.d0*eqpar(a_)**2*sin(mythetaks)**2)
             if (present(dginvdk_iszero)) dginvdk_iszero = .false.
             }{^IFZOUT
             if (present(dginvdk)) dginvdk = 0.d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
             }
          else
             if (present(dginvdk)) dginvdk = 0.0d0
             if (present(dginvdk_iszero)) dginvdk_iszero = .true.
          end if
       end if
       
    else
       
       ginv = 0.0d0
       if (present(iszero)) iszero = .true.

       if (present(dginvdk)) dginvdk = 0.0d0
       if (present(dginvdk_iszero)) dginvdk_iszero = .true.
       
    end if
       
!    if (present(dginvdk_iszero)) dginvdk_iszero = .false. ! At least we cant tell with certainty
!    if (present(dginvdk)) call mpistop('mks_get_gammainv_component_analytical: dginvdk not yet implemented.')

    
  end subroutine mks_get_gammainv_component_analytic
  !=============================================================================
  subroutine get_alpha(x^D,alpha,iszero,dalphadj_iszero,dalphadj,jdir,w_pt)

    include 'amrvacdef.f'

    ! get the lapse.  Optional parameter is true if lapse is
    ! identically zero (does not really make sense)
    ! Optional parameters jdir and dalphadj request derivatives
    ! \partial_j \alpha ; j=jdir
    double precision, intent(in)                     :: x^D
    integer, optional, intent(in)                    :: jdir
    double precision, intent(out)                    :: alpha
    logical, optional, intent(out)                   :: iszero, dalphadj_iszero
    double precision, optional, intent(out)          :: dalphadj
    double precision, intent(in), optional           :: w_pt(1:nw)

    ! .. local ..
    integer                                  :: i,j
    double precision                         :: myrks, mythetaks
    !-----------------------------------------------------------------------------
    if(present(dalphadj) .and. .not. present(jdir) .or. &
         present(dalphadj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_alpha: derivatives requested without direction or output-slot given.")

    if(present(iszero)) iszero = .false.

    myrks = rks(x1)
    {^IFZIN
    mythetaks = thetaks(x^Z)
    }
    {^IFZOUT
    mythetaks = dpi/2.0d0
    }  

    {^IFZIN
    alpha = sqrt((myrks**2+  &
         eqpar(a_)**2*Cos(mythetaks)**2)/(myrks*(2.0d0*eqpar(m_)  &
         +myrks)+eqpar(a_)**2*Cos(mythetaks)**2))
    }{^IFZOUT
    alpha =sqrt(myrks/(2.0d0*eqpar(m_)+myrks))
    }

    if (present(jdir)) then

       select case(jdir)

       case(1)
          ! Radial derivative:
          if (present(dalphadj)) then
             {^IFZIN
             dalphadj = (eqpar(m_)*Exp(x1)*(myrks**2-  &
                  1.0d0*eqpar(a_)**2*Cos(mythetaks)**2))/((myrks*(2.0d0*eqpar(m_)  &
                  +myrks)+eqpar(a_)**2*Cos(mythetaks)**2)**2*Sqrt(1-  &
                  (2.0d0*eqpar(m_)*myrks)/(myrks*(2*eqpar(m_)+myrks)  &
                  +eqpar(a_)**2*Cos(mythetaks)**2)))
             }{^IFZOUT
             dalphadj = (eqpar(m_)*exp(x1))/(Sqrt(myrks/(2.0d0*eqpar(m_)+  &
                  myrks))*(2.0d0*eqpar(m_)+myrks)**2)
             }
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .false.
          end if

       case(^Z)
          ! Theta-derivative:
          {^IFZIN
          if (present(dalphadj)) then
             dalphadj = (-1.0d0*eqpar(a_)**2*eqpar(m_)*myrks*(1  &
                  + coordpar(h_)*Cos(2*x^Z))*Sin(2*x^Z  &
                  + coordpar(h_)*Sin(2*x^Z)))/((myrks*(2*eqpar(m_)+  &
                  myrks)+eqpar(a_)**2*Cos(mythetaks)**2)**2*Sqrt(1-  &
                  (2*eqpar(m_)*myrks)/(myrks*(2*eqpar(m_)+myrks)+  &
                  eqpar(a_)**2*Cos(mythetaks)**2)))
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .false.
          end if
          }{^IFZOUT
          if (present(dalphadj)) then
             dalphadj = 0.0d0
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true.
          end if
          }

       case default
          ! Other cases (phi):
          if (present(dalphadj)) then
             dalphadj = 0.0d0
          end if
          if (present(dalphadj_iszero)) then
             dalphadj_iszero = .true.
          end if
       end select

    end if ! present(jdir)

  end subroutine get_alpha
  !=============================================================================
  subroutine get_beta(idir,x^D,beta,iszero,dbetaidj_iszero,dbetaidj,jdir,w_pt)

    include 'amrvacdef.f'

    ! get the (contravariant!!) shift vector.
    ! The optional argument iszero is true if shift-component is 
    ! identically zero.
    ! if requested, dbetaidj is the derivative of the contravariant shift.
    ! \partial_j \beta^i ; i=idir, j=jdir
    integer, intent(in)                      :: idir
    double precision, intent(in)             :: x^D
    integer, optional, intent(in)            :: jdir
    double precision, intent(out)            :: beta
    logical, optional, intent(out)           :: iszero, dbetaidj_iszero
    double precision, optional, intent(out)  :: dbetaidj
    double precision, intent(in), optional   :: w_pt(1:nw)

    ! .. local ..
    double precision                         :: myrks, mythetaks
    !-----------------------------------------------------------------------------
    if(present(dbetaidj) .and. .not. present(jdir) .or. &
         present(dbetaidj_iszero) .and. .not. present(jdir)) &
         call mpistop("get_beta: derivatives requested &
         &without direction or output-slot given.")

    myrks = rks(x1)
    {^IFZIN
    mythetaks = thetaks(x^Z)
    }
    {^IFZOUT
    mythetaks = dpi/2.0d0
    }  


    select case(idir)

    case(1)
       ! Betar
       {^IFZIN
       beta = (2.0d0*eqpar(m_)*myrks)/(exp(x1)*(myrks*(2*eqpar(m_)  &
            +myrks)+eqpar(a_)**2*Cos(mythetaks)**2))
       }{^IFZOUT
       beta = (2.0d0*eqpar(m_))/(exp(x1)*(2*eqpar(m_)+myrks))
       }

       if(present(iszero)) then
          iszero = .false.
       end if

       if (present(jdir)) then
          select case(jdir)

          case(1)
             ! dbetardr:
             if (present(dbetaidj)) then
                {^IFZIN
            dbetaidj = (-2.0d0*eqpar(m_)*((2.0d0*exp(x1)+2.0d0*eqpar(m_)  &
                 + coordpar(R0_))*myrks**2  &
                 + eqpar(a_)**2*coordpar(R0_)*Cos(mythetaks)**2))/(exp(x1)*(myrks &
                 * (2.0d0 * eqpar(m_)  &
                 + myrks)+eqpar(a_)**2*Cos(mythetaks)**2)**2)              
                }{^IFZOUT
                dbetaidj = (-2.0d0*eqpar(m_)*(2*exp(x1)+2*eqpar(m_)+  &
                     coordpar(R0_)))/(exp(x1)*(2.0d0*eqpar(m_)+myrks)**2)
                }
             end if ! present(dbetaidj)

             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .false.
             end if
          case(^Z)
             ! dbetardtheta:
             {^IFZIN
             if (present(dbetaidj)) then
                dbetaidj = (2.0d0*eqpar(a_)**2*eqpar(m_)*myrks*exp(-x1)*(1  &
                     + coordpar(h_)*Cos(2*x^Z))*Sin(2*x^Z  &
                     + coordpar(h_)*Sin(2*x^Z)))/(myrks*(2*eqpar(m_)+  &
                     myrks)+eqpar(a_)**2*Cos(mythetaks)**2)**2
             end if
             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .false.
             end if
             }{^IFZOUT
             if (present(dbetaidj)) then
                dbetaidj = 0.0d0
             end if
             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .true.
             end if
             }


          case default
             ! dbbetardphi is zero
             if (present(dbetaidj)) then
                dbetaidj = 0.0d0
             end if ! present(dbetaidj)           

             if (present(dbetaidj_iszero)) then
                dbetaidj_iszero = .true.
             end if

          end select
       end if ! present(jdir)


    case default
       ! All zeroes except r:
       beta = 0.0d0
       if(present(iszero))           iszero          = .true.
       if (present(dbetaidj))         dbetaidj       = 0.0d0
       if (present(dbetaidj_iszero)) dbetaidj_iszero = .true.

    end select

  end subroutine get_beta
  !=============================================================================
  subroutine get_g_component(iin,jin,x^D,g,iszero,dgdk_iszero,dgdk,kdir,w_pt)

    include 'amrvacdef.f'

    ! This is at the heart of the scheme: Set the (spatial) metric components here
    ! and only here...
    ! Indices of the metric are down (covariant) g_{ij}
    ! The optional argument iszero is true if the element is identically zero
    ! The optional arguments dgdk and kdir request derivatives of the metric
    ! \partial_k g_{ij} ; i=iin, j=jin, k=kdir
    ! The optional argument dgdk_iszero
    integer, intent(in)                      :: iin,jin
    integer, optional, intent(in)            :: kdir
    double precision, intent(in)             :: x^D
    double precision, intent(out)            :: g
    logical, optional, intent(out)           :: iszero, dgdk_iszero
    double precision, optional, intent(out)  :: dgdk
    double precision, intent(in), optional   :: w_pt(1:nw)

    ! .. local ..
    integer                                  :: i,j
    double precision                         :: myrks, mythetaks
    !-----------------------------------------------------------------------------
    if(present(dgdk) .and. .not. present(kdir) .or. &
         present(dgdk_iszero) .and. .not. present(kdir)) &
         call mpistop("get_g_component: derivatives requested without &
         &direction or output-slot given.")

    ! metric is symmetric: swap indices if needed:
    ! User needs only to provide values for i<=j (upper triangle).  
    if (iin .gt. jin) then
       i=jin; j=iin
    else
       i=iin; j=jin
    end if

    if(present(iszero)) iszero = .false.

    myrks = rks(x1)
    {^IFZIN
    mythetaks = thetaks(x^Z)
    }
    {^IFZOUT
    mythetaks = dpi/2.0d0
    }  
    ! Non-diagonal!
    if (i .eq. 1 .and. j .eq. 1) then
       !grr:
       {^IFZIN
       g = exp(2.0d0*x1)*(1.0d0  &
            + (2.0d0*eqpar(m_)*myrks)/(myrks**2  &
            + eqpar(a_)**2*Cos(mythetaks)**2))
       }{^IFZOUT
       g = (exp(2.0d0*x1)*(2.0d0*eqpar(m_)+myrks))/myrks
       }
       if (present(kdir)) then
          ! Derivative info requested

          select case(kdir)

          case (1)
             ! dgrrdr:
             if (present(dgdk)) then
                {^IFZIN
                dgdk = (2.0d0*exp(2*x1)*(eqpar(a_)**2*(2*coordpar(R0_)*(coordpar(R0_)  &
                     +eqpar(m_))+(4*coordpar(R0_)+3*eqpar(m_))*exp(x1)+  &
                     2*exp(2*x1))*Cos(mythetaks)**2+  &
                     eqpar(a_)**4*Cos(mythetaks)**4+  &
                     myrks**2*(coordpar(R0_)*(coordpar(R0_)+2*eqpar(m_))+  &
                     (2*coordpar(R0_)+eqpar(m_))*exp(x1)+  &
                     exp(2*x1))))/(myrks**2+  &
                     eqpar(a_)**2*Cos(mythetaks)**2)**2
                }{^IFZOUT
                dgdk = (2.0d0*exp(2*x1)*(coordpar(R0_)**2+  &
                     2*coordpar(R0_)*(eqpar(m_)+exp(x1))+exp(x1)*(eqpar(m_)  &
                     +exp(x1))))/myrks**2
                }
             end if
             if (present(dgdk_iszero)) dgdk_iszero = .false.

          case(^Z)
             !dgrrdtheta:
             {^IFZIN
             if (present(dgdk)) then
                dgdk = (2.0d0*eqpar(a_)**2*eqpar(m_)*myrks*exp(2*x1)*(1  &
                     + coordpar(h_)*Cos(2*x^Z))*Sin(2*x^Z  &
                     + coordpar(h_)*Sin(2*x^Z)))/(myrks**2+  &
                     eqpar(a_)**2*Cos(mythetaks)**2)**2
             end if
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) then
                dgdk = 0.0d0
             end if
             if (present(dgdk_iszero)) dgdk_iszero = .true.
             }

          case default
             ! phi derivatives zero:
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.

          end select
       end if

    else if (i .eq. ^PHI .and. j .eq. ^PHI) then
       ! gphiphi
       {^IFZIN
       g = Sin(mythetaks)**2*(eqpar(a_)**2+myrks**2  &
            + (2.0d0*eqpar(a_)**2*eqpar(m_)*myrks*Sin(mythetaks)**2)/(myrks**2  &
            + eqpar(a_)**2*Cos(mythetaks)**2))
       }{^IFZOUT
       ! equatorial plane
       g = myrks**2+(eqpar(a_)**2*(2.0d0*eqpar(m_)  &
            + myrks))/myrks
       }

       if (present(kdir)) then
          ! Derivative info requested
          select case(kdir)
          case(1)
             !DgphiphiDr
             {^IFZIN
             if (present(dgdk)) dgdk = 2.0d0*exp(x1)*Sin(mythetaks)**2*(myrks+  &
                  (eqpar(a_)**2*eqpar(m_)*(-myrks**2+  &
                  eqpar(a_)**2*Cos(mythetaks)**2)*Sin(mythetaks)**2)/(myrks**2  &
                  +eqpar(a_)**2*Cos(mythetaks)**2)**2)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 2.0d0*exp(x1)*(-((eqpar(a_)**2*eqpar(m_))/myrks**2)  &
                  + myrks)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }
          case(^Z)
             !DgphiphiDtheta
             {^IFZIN
             if (present(dgdk)) then
                dgdk = (2.0d0*Cos(mythetaks)*(1+  &
                     coordpar(h_)*Cos(2*x^Z))*Sin(mythetaks)*(2*eqpar(a_)**2*eqpar(m_) &
                     *myrks*(eqpar(a_)**2  &
                     +myrks**2)*Sin(mythetaks)**2+(myrks**2+  &
                     eqpar(a_)**2*Cos(mythetaks)**2)*(eqpar(a_)**2*(myrks**2  &
                     + eqpar(a_)**2*Cos(mythetaks)**2)+myrks**2*(myrks**2+  &
                     eqpar(a_)**2*Cos(mythetaks)**2)+  &
                     2*eqpar(a_)**2*eqpar(m_)*myrks*Sin(mythetaks)**2)))/(myrks**2  &
                     +eqpar(a_)**2*Cos(mythetaks)**2)**2
             end if 
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
             }
          case default
             !DphiphiDphi
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.
          end select
       end if

    else if (i .eq. ^Z .and. j .eq. ^Z) then
       ! gthetatheta
       {^IFZIN
       g = (myrks**2+eqpar(a_)**2*Cos(mythetaks)**2)*(1+  &
            coordpar(h_)*Cos(2.0d0*x^Z))**2
       }{^IFZOUT
       g = (-1+coordpar(h_))**2*myrks**2
       }

       if (present(kdir)) then
          ! Derivative info requested
          select case(kdir)

          case(1)
             ! DgthetathetaDs:
             {^IFZIN
             if (present(dgdk)) dgdk = 2.0d0*myrks*exp(x1)*(1+coordpar(h_)*Cos(2*x^Z))**2
             }{^IFZOUT
             if (present(dgdk)) dgdk = 2.0d0*(-1+coordpar(h_))**2*myrks*exp(x1)
             }
             if (present(dgdk_iszero)) dgdk_iszero = .false.

          case(^Z)
             ! DthetathetaDtheta:
             {^IFZIN
             if (present(dgdk)) dgdk = -4.0d0*coordpar(h_)*(myrks**2  &
                  + eqpar(a_)**2*Cos(mythetaks)**2)*(1  &
                  + coordpar(h_)*Cos(2*x^Z))*Sin(2*x^Z)  &
                  - 1.0d0*eqpar(a_)**2*(1  &
                  + coordpar(h_)*Cos(2*x^Z))**3*Sin(2*x^Z  &
                  + coordpar(h_)*Sin(2*x^Z))
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.              
             }

          case default
             ! DthetathetaDphi
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.

          end select
       end if

    else if (i .eq. 1 .and. j .eq. ^PHI) then
       ! grphi
       {^IFZIN
       g = -1.0d0*exp(x1)*eqpar(a_)*(1.0d0  &
            + (2.0d0*eqpar(m_)*myrks)/(myrks**2  &
            + eqpar(a_)**2*Cos(mythetaks)**2))*Sin(mythetaks)**2
       }{^IFZOUT
       g = (-1.0d0*exp(x1)*eqpar(a_)*(2.0d0*eqpar(m_)  &
            + myrks))/myrks
       }
 
       if (present(kdir)) then
          ! Derivative info requested
          select case(kdir)

          case(1)
             ! DgrphiDr:
             {^IFZIN
             if (present(dgdk)) dgdk = (-1.0d0*eqpar(a_)*exp(x1)*(2*eqpar(a_)**2 &
                  *(coordpar(R0_)*(coordpar(R0_)  &
                  +eqpar(m_))+2*(coordpar(R0_)+eqpar(m_))*exp(x1)+  &
                  exp(2*x1))*Cos(mythetaks)**2+  &
                  eqpar(a_)**4*Cos(mythetaks)**4+  &
                  myrks**2*(coordpar(R0_)*(coordpar(R0_)+2*eqpar(m_))+  &
                  2*coordpar(R0_)*exp(x1)+  &
                  exp(2*x1)))*Sin(mythetaks)**2)/(myrks**2+  &
                  eqpar(a_)**2*Cos(mythetaks)**2)**2
             }{^IFZOUT
             if (present(dgdk)) dgdk = (-1.0d0*eqpar(a_)*exp(x1)*(coordpar(R0_)**2+  &
                  2*coordpar(R0_)*(eqpar(m_)+exp(x1))+  &
                  exp(2*x1)))/myrks**2
             }
             if (present(dgdk_iszero)) dgdk_iszero = .false.

          case(^Z)
             ! DgrphiDtheta:
             {^IFZIN
             if (present(dgdk)) dgdk = 2.0d0*eqpar(a_)*exp(x1)*Cos(mythetaks)*(1+  &
                  coordpar(h_)*Cos(2*x^Z))*Sin(mythetaks)*(-1-  &
                  (2*eqpar(m_)*myrks)/(myrks**2+  &
                  eqpar(a_)**2*Cos(mythetaks)**2)-  &
                  (2*eqpar(a_)**2*eqpar(m_)*myrks*Sin(mythetaks)**2)/(myrks**2  &
                  +eqpar(a_)**2*Cos(mythetaks)**2)**2)
             if (present(dgdk_iszero)) dgdk_iszero = .false.
             }{^IFZOUT
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.              
             }

          case default
             ! DrphiDphi
             if (present(dgdk)) dgdk = 0.0d0
             if (present(dgdk_iszero)) dgdk_iszero = .true.

          end select
       end if

    else

       g = 0.0d0

       if(present(iszero))      iszero       = .true.
       if (present(dgdk))       dgdk         = 0.0d0
       if(present(dgdk_iszero)) dgdk_iszero  = .true.

    end if

  end subroutine get_g_component
  !=============================================================================
  double precision function outerhorizon()
    
    include 'amrvacdef.f'
    double precision                       :: rksout
    !-----------------------------------------------------------------------------

    rksout = eqpar(m_) + sqrt(eqpar(m_)**2 - eqpar(a_)**2)
    outerhorizon = s(rksout)
    
  end function outerhorizon
  !=============================================================================
  subroutine u4BLtoCoord(ixI^L,ixO^L,xBL,u4BL,u4Coord)

    ! Transforms the (contravariant) four-velocity u4BL from Boyer-Lindquist coordinates
    ! to the current (modified KS) coordinates u4Coord.  Often initial conditions are
    ! given in terms of BL coordinates and this routine comes in handy.
    !
    !
    ! utKS       = utBL + 2 M rBL / Delta * urBL
    ! urKS       = e**-s * urBL
    ! uthetaKS   = uthetaBL / (1+h*Cos(2.0d0*vartheta))
    ! uphiKS     = uphiBL + a/Delta * urBL
    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4BL
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4Coord
    ! .. local ..
    integer                                                :: ix^D
    double precision, dimension(ixI^S)                     :: Delta, twoMr, ems, varth
    !-----------------------------------------------------------------------------

    ems(ixO^S)   = one/(xBL(ixO^S,1) - coordpar(R0_))

    {^IFZIN
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    varth(ix^D) = vartheta(xBL(ix^D,^Z))
    {end do \}
    }{^IFZOUT
    varth(ixO^S) = dpi/2.0d0
    }

    twoMr(ixO^S) = 2.0d0*eqpar(m_)*abs(xBL(ixO^S,1))
    Delta(ixO^S) = xBL(ixO^S,1)**2 - twoMr(ixO^S) + eqpar(a_)**2

    u4Coord(ixO^S,0) = u4BL(ixO^S,0) + twoMr(ixO^S)/Delta(ixO^S) * u4BL(ixO^S,1)
    u4Coord(ixO^S,1) = ems(ixO^S) * u4BL(ixO^S,1)

    {^IFZ
    u4Coord(ixO^S,^Z) = 1/(1+coordpar(h_)*Cos(2.0d0*varth(ixO^S))) * u4BL(ixO^S,^Z)
    }
    {^IFPHI
    u4Coord(ixO^S,^PHI) = u4BL(ixO^S,^PHI) + eqpar(a_)/Delta(ixO^S) * u4BL(ixO^S,1)
    }

  end subroutine u4BLtoCoord
  !=============================================================================
  subroutine BLToCoord(ixI^L,ixO^L,xBL,xCoord)

    ! Transforms the (spatial) Boyer-Lindquist coordinates to the current coordinates
    ! We don't return the new time-component although it changes.
    !
    use mod_transform, only: BLToKS
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xBL
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    !-----------------------------------------------------------------------------

    call BLToKS(ixI^L,ixO^L,xBL,xKS)
    call KSToCoord(ixI^L,ixO^L,xKS,xCoord)

  end subroutine BLToCoord
  !=============================================================================
  subroutine CoordToBL(ixI^L,ixO^L,xCoord,xBL)

    ! Transforms the (spatial) modified KS coordinate to Boyer-Lindquist coordinates.
    ! We don't return the new time-component although it changes.
    !
    use mod_transform, only: KSToBL
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xBL
    ! .. local ..
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    !-----------------------------------------------------------------------------

    call CoordToKS(ixI^L,ixO^L,xCoord,xKS)
    call KSToBL(ixI^L,ixO^L,xKS,xBL)

  end subroutine CoordToBL

  subroutine CoordToBL_pt(xCoord,xBL)
    ! Transforms the (spatial) modified KS coordinate to Boyer-Lindquist coordinates.
    ! We don't return the new time-component although it changes.
    !
    use mod_transform, only: KSToBL_pt
    include 'amrvacdef.f'

    double precision, dimension(1:^ND), intent(in)   :: xCoord
    double precision, dimension(1:^ND), intent(out)  :: xBL
    ! .. local ..
    double precision, dimension(1:^ND)               :: xKS

    call CoordToKS_pt(xCoord,xKS)
    call KSToBL_pt(xKS,xBL)
  end subroutine CoordToBL_pt

  !=============================================================================
  subroutine KSToCoord(ixI^L,ixO^L,xKS,xCoord)

    ! Transforms the KS coordinate to the modified KS-coordinates.
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCoord
    ! .. local ..
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    xCoord(ix^D,1) = s(xKS(ix^D,1))
    {^IFZIN
    xCoord(ix^D,^Z) = vartheta(xKS(ix^D,^Z))
    }
    {^IFPHIIN
    xCoord(ix^D,^PHI) = xKS(ix^D,^PHI)
    }
    {end do\}

  end subroutine KSToCoord
  !=============================================================================
  subroutine CoordToKS(ixI^L,ixO^L,xCoord,xKS)

    ! Transforms the modified KS coordinate to KS-coordinates.
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xKS
    ! .. local ..
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    xKS(ix^D,1) = rks(xCoord(ix^D,1))
    {^IFZIN
    xKS(ix^D,^Z) = thetaks(xCoord(ix^D,^Z))
    }
    {^IFPHIIN
    xKS(ix^D,^PHI) = xCoord(ix^D,^PHI)
    }
    {end do\}

  end subroutine CoordToKS

  subroutine CoordToKS_pt(xCoord,xKS)

    ! Transforms the modified KS coordinate to KS-coordinates.
    !
    include 'amrvacdef.f'

    double precision, dimension(1:^ND), intent(in)   :: xCoord
    double precision, dimension(1:^ND), intent(out)  :: xKS
    !-----------------------------------------------------------------------------

    xKS(1) = rks(xCoord(1))
    {^IFZIN
    xKS(^Z) = thetaks(xCoord(^Z))
    }
    {^IFPHIIN
    xKS(^PHI) = xCoord(^PHI)
    }

  end subroutine CoordToKS_pt


  !=============================================================================
  subroutine MKS_u4CoordToKS(ixI^L,ixO^L,xCoord,u4Coord,u4KS,J)

    ! Transforms the modified KS four vector u4Coord to
    ! the ordinary KS four vector u4KS.
    !
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4Coord
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4KS
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    integer                                                :: ix^D
    !-----------------------------------------------------------------------------

    u4KS(ixO^S,0)    = u4Coord(ixO^S,0)
    u4KS(ixO^S,r_)   = exp(xCoord(ixO^S,r_)) * u4Coord(ixO^S,r_)
    {^IFZIN
    u4KS(ixO^S,z_)   = (1.0d0 + coordpar(h_)*cos(2.0d0*xCoord(ixO^S,z_))) * u4Coord(ixO^S,z_)
    }
    {^IFPHI
    u4KS(ixO^S,phi_) = u4Coord(ixO^S,phi_)
    }

    if (present(J)) then
       J(ixO^S,:,:) = 0.0d0
       J(ixO^S,r_,r_) = exp(xCoord(ixO^S,r_))
       {^IFZIN
       J(ixO^S,z_,z_) = (1.0d0 + coordpar(h_)*cos(2.0d0*xCoord(ixO^S,z_)))
       }
    end if

  end subroutine MKS_u4CoordToKS
  !=============================================================================
  subroutine MKS_u4CoordToBL(ixI^L,ixO^L,xCoord,u4Coord,u4BL,J)

    ! Transforms the modified KS four vector u4Coord to
    ! the BL four vector u4BL.
    !

    use mod_transform, only: u4KStoBL, compose
    include 'amrvacdef.f'

    integer,intent(in)                                     :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: xCoord
    double precision, dimension(ixI^S,0:^NC), intent(in)   :: u4Coord
    double precision, dimension(ixI^S,0:^NC), intent(out)  :: u4BL
    double precision, dimension(ixI^S,0:^NC,0:^NC), optional, intent(out)  :: J
    ! .. local ..
    integer                                                :: ix^D
    double precision, dimension(ixI^S,0:^NC)               :: u4KS
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    double precision, dimension(ixI^S,0:^NC,0:^NC)         :: Jac1, Jac2
    !-----------------------------------------------------------------------------

    u4KS(ixO^S,0)    = u4Coord(ixO^S,0)
    u4KS(ixO^S,r_)   = exp(xCoord(ixO^S,r_)) * u4Coord(ixO^S,r_)
    {^IFZIN
    u4KS(ixO^S,z_)   = (1.0d0 + coordpar(h_)*cos(2.0d0*xCoord(ixO^S,z_))) * u4Coord(ixO^S,z_)
    }
    {^IFPHI
    u4KS(ixO^S,phi_) = u4Coord(ixO^S,phi_)
    }

    Jac1(ixO^S,:,:) = 0.0d0
    Jac1(ixO^S,r_,r_) = exp(xCoord(ixO^S,r_))
    {^IFZIN
    Jac1(ixO^S,z_,z_) = (1.0d0 + coordpar(h_)*cos(2.0d0*xCoord(ixO^S,z_)))
    }
    
    call CoordToKS(ixI^L,ixO^L,xCoord,xKS)
    call u4KStoBL(ixI^L,ixO^L,xKS,u4KS,u4BL,J=Jac2)

    if (present(J)) then
       call compose(ixI^L,ixO^L,Jac2,Jac1,J)
    end if
    
  end subroutine MKS_u4CoordToBL
  !=============================================================================
  subroutine CoordToCart(ixI^L,ixO^L,x,xCart)

    ! Transforms the coordinates to "Kerr-Schild Cartesian coordinates"
    ! Two transformations are needed:  
    ! rKS     = R0 + exp(s)
    ! thetaKS = vartheta + h/2 Sin(2 vartheta)
    ! phiKS   = phiMKS
    ! 
    ! x = rKS Sin(thetaKS) Cos(phiKS)
    ! y = rKS Sin(thetaKS) Sin(phiKS)
    ! z = rKS Cos(thetaKS)
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'
    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCart
    ! .. local ..
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, ra
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    integer                                                :: zcart_, ycart_
    !-----------------------------------------------------------------------------

    if (ndim.eq.3) then
       ycart_=2
       zcart_=3
    else if (ndim.eq.2 .and. ^Z.eq.2) then
       ycart_=1 ! must catch ycart_=1 case and don't use
       zcart_=2
    else if (ndim.eq.2 .and. ^PHI.eq.2) then
       ycart_=2
       zcart_=1 ! must catch zcart_=1 case and don't use
    else if (ndim.eq.1) then
       ycart_=1
       zcart_=1
    else
       call mpistop("CoordToCart: unknown parameter combination of -d, -phi and -z!")
    end if
    !-----------------------------------------------------------------------------

    call CoordToKS(ixI^L,ixO^L,x,xKS)

    {^IFPHIIN
    sph(ixO^S) = sin(xKS(ixO^S,phi_))
    cph(ixO^S) = cos(xKS(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(xKS(ixO^S,z_))
    cth(ixO^S) = cos(xKS(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }

    ra(ixO^S) = xKS(ixO^S,1)!sqrt(xKS(ixO^S,1)**2+eqpar(a_)**2)

    xCart(ixO^S,1) = ra(ixO^S) * sth(ixO^S)*cph(ixO^S)
    if (ycart_.ne.1) &
         xCart(ixO^S,ycart_) = ra(ixO^S) * sth(ixO^S)*sph(ixO^S)
    if (zcart_.ne.1) &
         xCart(ixO^S,zcart_) = ra(ixO^S) * cth(ixO^S)

  end subroutine CoordToCart
  !=============================================================================
  subroutine u3CoordToCart(ixI^L,ixO^L,x,u3Coord,u3Cart,J)

    ! Transforms any contravariant three vector u3Coord in the current coordinates
    ! to a vector in the basis of "Kerr-Schild Cartesian coordinates".
    !
    ! Two transformations are needed:  
    ! rKS     = R0 + exp(s)
    ! thetaKS = vartheta + h/2 Sin(2 vartheta)
    ! phiKS   = phiMKS
    ! 
    ! x = rKS Sin(thetaKS) Cos(phiKS)
    ! y = rKS Sin(thetaKS) Sin(phiKS)
    ! z = rKS Cos(thetaKS)
    ! 
    ! The transformation matrix J=del(x,y,z)/del(s,vartheta,phiMKS)
    ! Since the coordinates don't depend on time, this is the same transformation
    ! as for the four-vector u4CoordToCart()
    !
    ! u3Cart = J u3Coord
    !
    ! see modKS.nb for details
    !-----------------------------------------------------------------------------

    include 'amrvacdef.f'

    integer, intent(in)                                    :: ixI^L, ixO^L
    double precision, dimension(ixI^S,1:^ND), intent(in)   :: x
    double precision, dimension(ixI^S,1:^NC), intent(in)   :: u3Coord
    double precision, dimension(ixI^S,1:^NC), intent(out)  :: u3Cart
    double precision, dimension(ixI^S,1:^NC,1:^NC), optional, intent(out)  :: J
    ! .. local ..
    double precision, dimension(ixI^S,1:^NC,1:^NC)         :: Jac
    double precision, dimension(ixI^S,1:^ND)               :: xKS
    double precision, dimension(ixI^S)                     :: sth, cth, sph, cph, tmp
    integer                                                :: zcart_, ycart_, ix, jx
    !-----------------------------------------------------------------------------

    if (ndim.eq.3.and.ndir.eq.3) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. ^PHI .eq. 2) then
       ycart_=2
       zcart_=3
    else if ((ndim.eq.2.or.ndim.eq.1).and.ndir.eq.3 .and. ^Z .eq. 2) then
       ycart_=3
       zcart_=2
    else if (ndir.eq.2 .and. ^Z.eq.2) then
       ycart_=1 ! must catch ycart_=1 case and don't use
       zcart_=2
    else if (ndir.eq.2 .and. ^PHI.eq.2) then
       ycart_=2
       zcart_=1 ! must catch zcart_=1 case and don't use
    else if (ndir.eq.1) then
       ycart_=1
       zcart_=1
    else
       call mpistop("u3CoordToCart: unknown parameter combination of -d, -phi and -z!")
    end if

    call CoordToKS(ixI^L,ixO^L,x,xKS)

    {^IFPHIIN
    sph(ixO^S) = sin(xKS(ixO^S,phi_))
    cph(ixO^S) = cos(xKS(ixO^S,phi_))
    }{^IFPHIOUT
    sph(ixO^S) = zero
    cph(ixO^S) = one
    }

    {^IFZIN
    sth(ixO^S) = sin(xKS(ixO^S,z_))
    cth(ixO^S) = cos(xKS(ixO^S,z_))
    }{^IFZOUT
    sth(ixO^S) = one
    cth(ixO^S) = zero
    }

    {^IFZIN
    tmp(ixO^S) = 1+coordpar(h_)*Cos(2.0d0*x(ixO^S,^Z))
    }{^IFZOUT
    tmp(ixO^S) = 1.0d0 - coordpar(h_)
    }

    Jac(ixO^S,1,1) = exp(x(ixO^S,1)) * sth(ixO^S)*cph(ixO^S)
    {^IFZ
    Jac(ixO^S,1,^Z) = (cph(ixO^S)*cth(ixO^S)*tmp(ixO^S)*xKS(ixO^S,1))
    }
    {^IFPHI
    Jac(ixO^S,1,^PHI) = - sph(ixO^S)*sth(ixO^S)*xKS(ixO^S,1)
    }

    if (ycart_.ne.1) then
       Jac(ixO^S,ycart_,1) = exp(x(ixO^S,1))*sph(ixO^S)*sth(ixO^S)
       {^IFZ
       Jac(ixO^S,ycart_,^ZZ) = cth(ixO^S)*sph(ixO^S)*tmp(ixO^S)*xKS(ixO^S,1)
       }
       {^IFPHI
       Jac(ixO^S,ycart_,^PPHI) = cph(ixO^S)*sth(ixO^S)*xKS(ixO^S,1)
       }
    end if

    if(zcart_.ne.1) then
       Jac(ixO^S,zcart_,1) = cth(ixO^S)*exp(x(ixO^S,1))
       {^IFZ
       Jac(ixO^S,zcart_,^ZZ) = - sth(ixO^S)*tmp(ixO^S)*xKS(ixO^S,1)
       }
       {^IFPHI
       Jac(ixO^S,zcart_,^PPHI) = 0.0d0
       }
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Jacobian fully assembled, now
    ! transform contravariant three-vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    u3Cart(ixO^S,:) = zero
    do ix=1,ndir
       do jx=1,ndir
          u3Cart(ixO^S,ix) = u3Cart(ixO^S,ix) + Jac(ixO^S,ix,jx) * u3Coord(ixO^S,jx)
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (present(J)) J = Jac


  end subroutine u3CoordToCart
  !=============================================================================
  subroutine MKS_LRNFu(ixI^L,ixO^L,m,ehatu)

    ! Obtains orthonormal tetrad vectors in matrix form.
    ! This is the locally non-rotating reference frame (LRNF)
    ! which to my understanding is just the Eulerian one, but its normalized.  
    ! Follows Rohta Takahashi, MNRAS, 383, 1155-116, 2008
    ! His equations (11) - (14)
    ! Should obey the same structure as KS case.
    ! This takes in already the MKS metric datastructure 'm', I might regret
    ! this choice later on. 
    
    include 'amrvacdef.f'
    integer, intent(in)                                         :: ixI^L, ixO^L
    type(metric)                                                :: m
    double precision, dimension(ixI^S,0:^NC,0:^NC), intent(out) :: ehatu
    ! .. local ..
    !-----------------------------------------------------------------------------

    ehatu = zero
    
    ehatu(ixO^S,0,0) = m%alpha(ixO^S)

    ehatu(ixO^S,r_,0)  = m%beta(r_)%elem(ixO^S)/sqrt(m%gammainv(r_,r_)%elem(ixO^S))
    ehatu(ixO^S,r_,r_) = one/sqrt(m%gammainv(r_,r_)%elem(ixO^S))

    {^IFZ
    ehatu(ixO^S,z_,z_) = sqrt(m%g(z_,z_)%elem(ixO^S))
    }

    {^IFPHI
    ehatu(ixO^S,phi_,0)    = m%beta(r_)%elem(ixO^S) * m%g(r_,phi_)%elem(ixO^S) / sqrt(m%g(phi_,phi_)%elem(ixO^S))
    ehatu(ixO^S,phi_,r_)   = m%g(r_,phi_)%elem(ixO^S) / sqrt(m%g(phi_,phi_)%elem(ixO^S))
    ehatu(ixO^S,phi_,phi_) = sqrt(m%g(phi_,phi_)%elem(ixO^S))
    }
    
  end subroutine MKS_LRNFu
  !=============================================================================
  ! Dummies
  !=============================================================================
  subroutine get_dgammainvdk(x^D,dgammainvdk,k)
    ! Obtains the derivatives of the inverse _spatial_ metric
    ! \partial_k gamma^{ij} for a given k. 
    ! This is not required when initialising gamma, lapse and shift
    ! e.g. when set init_from_g4 = .false.

    double precision, intent(in)           :: x^D
    double precision, intent(out)          :: dgammainvdk(1:^NC,1:^NC)
    integer, intent(in)                    :: k
    ! .. local ..
    !-----------------------------------------------------------------------------

    call mpistop("get_dgammainvdk: Not required and not implemented.")

    dgammainvdk = 0.0d0
    
  end subroutine get_dgammainvdk
  !=============================================================================
  ! End of coordinate-specific definitions.
  !=============================================================================
