! remember reconstruction now
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
! include file amrvacpar.t
module amrvacpar
{#IFDEF UNIT_TESTS
 use mod_m1_tests
}
character*5,parameter:: typephys='imhd'            ! VACPHYS module name

{#IFNDEF GLM
{#IFNDEF EPSINF
{#IFNDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab m a'    ! Equation parameter names
}{#IFDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab gammae gammap m a'    ! Equation parameter names
}}
{#IFDEF EPSINF
{#IFNDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab epsfloor rho0floor rho1e m a'    ! Equation parameter names
}{#IFDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab gammae gammap epsfloor rho0floor rho1e m a'    ! Equation parameter names
}}}

{#IFDEF GLM
{#IFNDEF EPSINF
{#IFNDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab Cr m a'    ! Equation parameter names
}{#IFDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab gammae gammap Cr m a'    ! Equation parameter names
}}
{#IFDEF EPSINF
{#IFNDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab Cr epsfloor rho0floor rho1e m a'    ! Equation parameter names
}{#IFDEF ELECTRONS
CHARACTER*9,PARAMETER:: eqparname='gamma adiab gammae gammap Cr epsfloor rho0floor rho1e m a'    ! Equation parameter names
}}}

!---------------
{#IFNDEF UNIT_TESTS
!> Number of w variables (including cons prim)
!integer:: nw
INTEGER :: nw

! nw = nwflux + nprim + ...... etc

!> Number of conservatives variables
integer :: ncons
!> Number of w variables which need evolve
integer :: nwflux
!> Number of prim variables
integer :: nprim
}!-------------------
!> Number of variables with tracer
integer :: nwfluxtr
!> Number of auxiliary variables
integer :: nwaux
!> Number of variables with extra
integer :: nwextra
!> Number of metric variables
integer :: nmetric
!> Number of GW backreaction variables
integer :: nwgwbr

!> index of min ncons variables which need evolve
integer :: ncons_lo
!> index of max ncons variables which need evolve
integer :: ncons_hi

!> index of min nwflux variables 
integer :: nwflux_lo
!> index of max nwflux variables 
integer :: nwflux_hi

!> index of min prim hydro variables, exclude metric
integer :: nprim_lo
!> index of max prim hydro variables, exclude metric
integer :: nprim_hi

!> index of min nwfluxtr variables
integer :: nwfluxtr_lo
!> index of max nwfluxtr variables
integer :: nwfluxtr_hi

!> index of min aux variables
integer :: nwaux_lo
!> index of max aux variables
integer :: nwaux_hi

!> index of min nwextra variables
integer :: nwextra_lo
!> index of max nwextra variables
integer :: nwextra_hi

!> index of min metric variables
integer :: nmetric_lo
!> index of max metric variables
integer :: nmetric_hi

!> index of min GW br variables
integer :: nwgwbr_lo
!> index of max GW br variables
integer :: nwgwbr_hi

integer :: nm1rad
!> index of min m1 variables
integer :: nm1rad_lo
!> index of max m1 variables
integer :: nm1rad_hi

integer :: nm1rad_eas
!> index of min m1_eas variables
integer :: nm1rad_eas_lo
!> index of max m1_eas variables
integer :: nm1rad_eas_hi

!  nws not count inside nw
integer :: itmp

{#IFDEF STAGGERED
! Staggered magnetic field
integer, parameter :: bs0_ = 0
integer, parameter :: bs^D_=bs0_+^D
integer, parameter :: nws=^ND
! Staggered counterparts for centered variables:
! set in mod_variables
integer, allocatable  :: iws(:) 
! The number of staggered components of the vector fields is
! the number of dimensions. (In 2D, the 3rd component is centred)
!integer, parameter :: nws=^ND
}{#IFNDEF STAGGERED
integer, parameter :: nws=0
integer, allocatable  :: iws(:) 
}




! flow variables
!=====Conserve variables=====!
integer:: d_           = 1000
integer:: Dye_         = 1000
integer:: s0_          = 1000
{#IFNDEF UNIT_TESTS
integer:: s1_          = 1000
integer:: s2_          = 1000
integer:: s3_          = 1000
}
integer:: e_           = 1000
{#IFNDEF UNIT_TESTS
integer:: tau_         = 1000
}
integer:: b0_          = 1000
integer:: b1_          = 1000
integer:: b2_          = 1000
integer:: b3_          = 1000
integer:: rhos_        = 1000
!=====Primitive variables=====!
{#IFNDEF UNIT_TESTS
integer:: ye_          = 1000
integer:: cs2_         = 1000
integer:: rho_         = 1000
integer:: u0_          = 1000
}
integer:: u1_          = 1000
integer:: u2_          = 1000
integer:: u3_          = 1000
integer:: v0_          = 1000
integer:: v1_          = 1000
integer:: v2_          = 1000
integer:: v3_          = 1000
integer:: pp_          = 1000
{#IFNDEF UNIT_TESTS
integer:: T_eps_       = 1000
integer:: temp_        = 1000
integer:: eps_         = 1000
}
integer:: psi_         = 1000
integer:: Ds_          = 1000
integer:: Dse_         = 1000

integer:: Dtr1_        = 1000
integer:: Dtr2_        = 1000
integer:: Dtr3_        = 1000
!======Cutoff energy ========!
integer:: Depsinf_     = 1000
!======Advected electron density ========!
integer:: Drho0_       = 1000
!======Conserved electron density ========!
integer:: Drho1_       = 1000
!======Advected electron density v2======!
integer:: Dn0_          = 1000
!======Conserved electron density v2=====!
integer:: Dn_           = 1000

   !---start of aux vars ---!
   integer:: lfac_         = 1000
   integer:: xi_           = 1000

   !---- metric vars ---- !
   integer:: alp_metric_  = 1000
   integer:: beta_metric1_= 1000
   integer:: beta_metric2_= 1000
   integer:: beta_metric3_= 1000
   integer:: Xvec_metric1_= 1000
   integer:: Xvec_metric2_= 1000
   integer:: Xvec_metric3_= 1000
   integer:: psi_metric_  = 1000

   !---- GW back reaction ---- !
   integer:: U_br_     = 1000
   integer:: U_br1_    = 1000
   integer:: U_br2_    = 1000
   integer:: U_br3_    = 1000
   integer:: h_00_     = 1000
   integer:: delta_br_ = 1000
 
   ! TEST 
    integer :: foo_ = 1000
    integer :: foo1_ = 1000
    integer :: foo2_ = 1000
    integer :: foo3_ = 1000
    integer :: foo4_ = 1000
    integer :: foo5_ = 1000
    integer :: foo6_ = 1000
    integer :: foo7_ = 1000
    integer :: foo8_ = 1000
    integer :: foo9_ = 1000

  !------ M1 radiation -------!
  {#IFNDEF UNIT_TESTS
   {^KSP&
     integer :: nrad^KSP_     = 1000 
     integer :: erad^KSP_     = 1000 
     {^C& 
      integer :: frad^KSP^C_  = 1000
     \}
     ! energy rates
     integer :: kappa_a^KSP_   = 1000
     integer :: kappa_s^KSP_   = 1000
     integer :: Q_er^KSP_      = 1000 
     ! number rates 
     integer :: Q_nr^KSP_      = 1000 
     integer :: kappa_nr^KSP_  = 1000 
     integer :: tau_path^KSP_   = 1000 
   \}
   }

   ! relative difference of c2p
   integer:: c2p_d_       = 1000
   integer:: c2p_s1_      = 1000
   integer:: c2p_s2_      = 1000
   integer:: c2p_s3_      = 1000
   integer:: c2p_tau_     = 1000
   !--- end of aux vars --- !
!=============================!

integer:: ierr_           = 1000
integer:: nerr_           = 1000

integer:: se_             = 1000
integer:: tr1_            = 1000
integer:: tr2_            = 1000
integer:: tr3_            = 1000
integer:: epsinf_         = 1000
integer:: rho0_           = 1000
integer:: rho1_           = 1000
integer:: n0_             = 1000
integer:: n_              = 1000

integer:: s_              = 1000
!=============================!
! polar variable names
integer:: sr_             = 1000
integer:: sphi_           = 1000
integer:: sz_             = 1000
integer:: vr_             = 1000
integer:: vphi_           = 1000
integer:: vz_             = 1000
integer:: uz_             = 1000
integer:: ur_             = 1000
integer:: uphi_           = 1000
integer:: br_             = 1000
integer:: bphi_           = 1000
integer:: bz_             = 1000
integer:: ee_             = 1000

! set inside imhd activate
integer:: nvector  = 0                             ! No. vector vars
integer, dimension(100) :: iw_vector

integer              :: nflag_
integer, allocatable :: flags(:)
double precision, allocatable :: wflags(:)



! if define in amrvacdef, cannot allocate
double precision, allocatable :: entropycoef(:)
double precision, allocatable :: normvar(:)

character*131, allocatable :: typeB(:,:)
character*131, allocatable :: typeentropy(:)

logical, allocatable :: loglimit(:)
logical, allocatable :: writew(:)
logical, allocatable :: logflag(:)

!  parameters
integer,parameter:: fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4 ! Characteristic
integer,parameter:: entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8 ! waves
integer,parameter:: nworkroe=15


{#IFNDEF GLM
{#IFNDEF EPSINF
{#IFNDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2, m_=3, a_=4, neqpar=4     ! equation params
}{#IFDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2, gammae_=3, gammap_=4, m_=5, a_=6, neqpar=6     ! equation params
}}
{#IFDEF EPSINF
{#IFNDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2, epsfloor_=3, rho0floor_=4, rho1e_=5, m_=6, a_=7, neqpar=7     ! equation params
}{IFDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2,gammae_=3, gammap_=4, epsfloor_=5, rho0floor_=6, rho1e_=7, m_=8, a_=9, neqpar=9     ! equation params
}}}

{#IFDEF GLM
{#IFNDEF EPSINF
{#IFNDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2, Cr_=3, m_=4, a_=5, neqpar=5     ! equation params
}{#IFDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2, gammae_=3, gammap_=4, Cr_=5, m_=6, a_=7, neqpar=7     ! equation params
}}
{#IFDEF EPSINF
{#IFNDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2, Cr_=3, epsfloor_=4, rho0floor_=5, rho1e_=6, m_=7, a_=8,  neqpar=8     ! equation params
}{#IFDEF ELECTRONS
INTEGER,PARAMETER:: gamma_=1, adiab_=2, gammae_=3, gammap_=4, Cr_=5, epsfloor_=6, rho0floor_=7, rho1e_=8, m_=9, a_=10,  neqpar=10     ! equation params
}}}

!COMMON, DOUBLE PRECISION::minp,minrho,smallxi,smalltau{#IFNDEF SYNGE , govergminone}
!COMMON, DOUBLE PRECISION::limitvalue
!
!! xprob: problem box; iprob: problem
!COMMON, INTEGER:: iprob
!COMMON, DOUBLE PRECISION:: xprob^L

! end include file amrvacpar.t
!=============================================================================
end module amrvacpar
