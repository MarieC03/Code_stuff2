module mod_m1_fermi
  use mod_m1_constants
  
  implicit none
  public
  integer, parameter :: NfermiDim = 20000 !10000
  double precision, parameter :: tol1 = 1.0d-16 !1.0d-8  8 digits accuracy
  double precision, parameter :: ETA_MAX = 2.0d0 !not necessary or really used

  !************************************************************************
  
contains
 !===================================================================
 !==================== new fermi ==================================== 
 
 ! calculates fermi integral upto N=5 for negative and positive eta
 ! uses fermi_calc(..) for eta < 2 and direct_integral(...) for eta >2
  function fermi_dirac(eta,N)
    double precision :: fermi_dirac
    double precision, intent(in) :: eta
    integer, intent(in) :: N
    double precision :: FII(0:6,NfermiDim)
	double precision :: fermi
   double precision :: eta_set
	
	    call prepare_fermi(FII)
		
      if( N .eq. 0 ) then
       if ( eta < ETA_MAX) then
	      fermi_dirac = fermi_calc(0,eta,FII)
       else
         call direct_integral(0,eta,fermi)
	      fermi_dirac = fermi
       end if
    else if ( N .eq. 1 ) then 
       if ( eta < ETA_MAX) then
	      fermi_dirac = fermi_calc(1,eta,FII)
       else
         call direct_integral(1,eta,fermi)
	      fermi_dirac = fermi
       end if
    else
      if ( eta < ETA_MAX) then
	      fermi_dirac = fermi_calc(N,eta,FII)
       else
         call direct_integral(N,eta,fermi)
	      fermi_dirac = fermi
       end if
	end if 
  end function fermi_dirac

  subroutine direct_integral(N,X,Fout)
      integer, intent(in) :: N  
      double precision, intent(in) :: X
      double precision, intent(out) :: Fout
	  ! internal
	  double precision :: Z,F,Sum,DT,T
	  integer :: i
	  Fout=0.0d0
      DT=0.005d0
      Sum=0.0d0
      do i=1,NfermiDim 
        T=float(i)*DT
        F=T**N/(exp(T-X)+1.d0)
        Sum=Sum+DT*F
      end do
      Fout=Sum 
  end subroutine direct_integral
  
  ! Integer order N; real argument Z of function Li_N(Z)
  function LiN(N,Z)
      double precision :: LiN
      double precision, intent(in) :: Z
      integer, intent(in) :: N
      ! internals
	  double precision :: sum1,arg,argu
	  ! absolute values to estimate convergence
	  double precision :: Bval,Bval2
	  integer :: k
	  
	  LiN = 0.0d0
      if(Z.EQ.0.0d0) return
      if(N.EQ.0) then
        LiN=Z/(1.d0-Z)
        return
      endif
      if(N.EQ.1) then
        LiN= -DLOG(1.d0-Z)
        return
      endif
      if(abs(Z).ge.1.0d0) then
        return
      endif
      if(N.GT.1) then
        sum1=Z
        arg=Z
        do k=2,NfermiDim 
          arg=arg*Z
          argu = arg/float(k)**N
          sum1 = sum1+argu 
          Bval = abs(argu)
          Bval2 = abs(sum1)
          if(Bval/Bval2 .lt. tol1) exit ! accuracy
        end do
        !write(2,*) Z,float(k),argu,sum1,float(N) ! test convergence
      endif
      LiN=sum1
  end function LiN

  ! Integer order N=0,..,5; real argument X of Fermi integral
  ! FF_N(X)=\int_0^\infty dt t^N/(exp(t-X)+1)
  ! Taylor expansion for |arg|<1; Integral else  
  function fermi_calc(N,X,FII)
    double precision :: fermi_calc
    double precision, intent(in) :: X
    integer, intent(in) :: N
    double precision :: FII(0:6,NfermiDim)
	!internal
	double precision ::FNX,arg,range1
	double precision :: GammaF(0:10)
	integer :: iarg,irange
	
	fermi_calc = 10.d0
	irange = 10
	! Gamma function:
	GammaF(0) = 1.d0
	GammaF(1) = 1.d0
    GammaF(2) = 2.d0
    GammaF(3) = GammaF(2)*3.d0
    GammaF(4) = GammaF(3)*4.d0
    GammaF(5) = GammaF(4)*5.d0
	! Compute first F_n(X) ! different normalization
    arg = -exp(X)   ! <0
	if(abs(arg) .lt. 1.0d0) then
       FNX=-LiN(N+1,arg)    ! relation between functions
       !if(N.eq.0) write(*,*) 'taylor used for',arg 
    else
	   iarg = NINT(-arg*float(NfermiDim)/irange)
       if(iarg.gt.NfermiDim) then
          !write(*,*)'Argument out of range',iarg
          {#IFNDEF UNIT_TESTS
          write(76,*) "iarg,irange,arg,x",iarg,irange,arg,X
          call mpistop("Argument for X in fermi integral out of range")
          }
          {#IFDEF UNIT_TESTS
             write(99,*)"Argument for X in fermi integral out of range"
          }
       end if 
       FNX=-FII(N+1,iarg)
	   !if(N.eq.0) write(*,*) 'Storage used for',arg
    end if
    ! Normalize
	fermi_calc = GammaF(N)*FNX
  end function fermi_calc

  ! Fi(N,X) for -range < X < 0 in NfermiDim steps
  subroutine prepare_fermi(FII)
      double precision, intent(inout) :: FII(0:6,NfermiDim)
	   ! internal 
	   double precision :: SumTot(10)
	   double precision :: range1
      double precision :: DZ,X
	   integer :: i
	   range1 = 10.0d0
	   DZ = -range1/float(NfermiDim)
	   SumTot(:) = 0.0d0
	   do i=1,NfermiDim
	     X=float(i)*DZ
		 FII(1,i) = -DLOG(1.0d0-X)
         SumTot(1)=SumTot(1)+FII(1,i)/X*DZ
         FII(2,i)=SumTot(1)
         SumTot(2)=SumTot(2)+FII(2,i)/X*DZ
         FII(3,i)=SumTot(2)
         SumTot(3)=SumTot(3)+FII(3,i)/X*DZ
         FII(4,i)=SumTot(3)
         SumTot(4)=SumTot(4)+FII(4,i)/X*DZ
         FII(5,i)=SumTot(4)
         SumTot(5)=SumTot(5)+FII(5,i)/X*DZ
         FII(6,i)=SumTot(5)
	   end do	   
  end subroutine prepare_fermi
 !===============================================================
 !=============================================================== 
!KEN fixed this entire function.
!KEN We forgot to change this in Maries code
subroutine get_neutrino_temp_ixD(wrad,N,ix^D,Gamma,velU,Wlor,Jrad, eta, speciesKSP, fluid_Prim, av_energy, T_nu, T_fluid, timeCurrent1,x)    
   !> Calculates the neutrino effective temperature in MeV using 
   !> Fermi integral and average energy of neutrinos
   use mod_m1_internal
   use mod_m1_eas_param 
   use mod_Weakhub_reader, only: logtemp_min_IV
   use ieee_arithmetic
   {#IFDEF UNIT_TESTS
   use mod_m1_tests
   }
   {#IFNDEF UNIT_TESTS
   include "amrvacdef.f"
   }
   integer, intent(in) :: ix^D
   integer, intent(in) :: speciesKSP
   double precision, intent(inout) :: wrad(1:m1_numvars_internal)
   double precision, intent(inout) :: fluid_Prim(1:3)    ! not used
   double precision, intent(in)    :: Gamma, Wlor, Jrad
   double precision, intent(in)    :: velU(1:^NC) 
   double precision, intent(in)    :: N
   double precision, intent(in)    :: eta
   double precision, intent(out)   :: av_energy !> average neutrino energy (code unit)
   double precision, intent(out)   :: T_nu   !> neutrino temperature (MeV)
   double precision, intent(in)    :: T_fluid   !> fluid temperature
   double precision, intent(in)    :: timeCurrent1
   double precision, intent(in)    :: x(1:ndim)
   ! internal:
   double precision :: E         !> neutrino energy
   double precision :: n_dens    !> number density
   double precision :: Fv        !> F_i*v^i
   double precision :: av_energy_MEV !> average energy in MeV
   double precision :: av_energy1  !> average energy code unit
   double precision :: F3 !> Need this for small eta
   !double precision :: E_atmo = 1.d-16 !1d-14

     Fv  = {^C&wrad(m1_flux^C_)*velU(^C)+}

     E = wrad(m1_energy_) 
     n_dens = N/Gamma    !theoretically not needed

      ! neutrinos atmosphere?
      if((n_dens .lt. m1_E_atmo*100) .or. (E .lt. m1_E_atmo*100) .or. (T_fluid .le. exp(logtemp_min_IV)*1.05)) then
           av_energy = m1_E_atmo
           T_nu = 1.d-2
           return
      endif

        !> av_energy needed in MeV thus code unit -> MeV

        av_energy1 = Wlor*(E-Fv)/N
        !av_energy = Jrad/n_dens ! equivalent to av_energy = Wlor*(E-Fv)/N

        !> units of av_energy in MeV:
        !TODO make this pretty:
        av_energy_MEV = av_energy1 !* MNUC_MEV
        av_energy = av_energy_MEV
        
        ! Thoughts: I don't need MNUC_MEV and the value is already close to MEV,
        ! I could also ty something like *MNUC /Srqtg

        ! calculate the neutrino temperature from Ynu not e_av
        ! set atmosphere 1d-14

        F3 = fermi_dirac(eta,3)
        if (F3<1.0d-15) then
          T_nu = 1.0/3.0 * av_energy_MEV
          return
        endif

        !> neutrino temperature: now in MeV, fermi-Dirac dimensionless
        T_nu = fermi_dirac(eta,2)/F3*av_energy_MEV

        if (.not. ieee_is_finite(T_nu)) then
          write(*,*) "T_nu is NaN/Inf!"
          write(*,*) "eta=", eta
          write(*,*) "fermi_dirac(eta,2)=", fermi_dirac(eta,2)
          write(*,*) "fermi_dirac(eta,3)=", fermi_dirac(eta,3)
          write(*,*) "av_energy_MEV=", av_energy_MEV
          write(*,*) "E=", E
          write(*,*) "Fv=", Fv
          write(*,*) "N=", N
          write(*,*) "n_dens=", n_dens
        end if
 end subroutine get_neutrino_temp_ixD 

 subroutine blackbody_ixD(wrad,ix^D,speciesKSP,fluid_Prim,T_fluid,eta,Black_er,Black_nr)
      use mod_m1_internal 
      {#IFDEF UNIT_TESTS
      use mod_m1_tests
      }
      use mod_m1_eas_param
      {#IFNDEF UNIT_TESTS
      include "amrvacdef.f"
      }
      integer, intent(in) :: ix^D, speciesKSP
      double precision, intent(in) :: eta
      double precision, intent(in) :: T_fluid 
      double precision, intent(inout) :: wrad(1:m1_numvars_internal)  ! not used
      double precision, intent(in) :: fluid_Prim(1:3)          ! not used
      double precision, intent(out) :: Black_er,Black_nr
      ! internal: 
      double precision, dimension(1:6) :: g_fact
      double precision :: Black_er_MEV, Black_nr_MEV
      double precision :: Units, Units2, UnitsCGS_er, UnitsCGS_nr
   
      Units2 = 4*PI/C2_CGS/HC_MEVCM  !not used 

      Units = 4*PI*CLITE_CM/(HC_MEVCM)**3
      !-----------
      if(m1_use_muons) then
        if(m1_use_neutrinos) then
            g_fact(m1_i_nue) = 1.0d0  ! nu_e
            g_fact(m1_i_nuebar) = 1.0d0  ! nu_bar_e
            g_fact(m1_i_nux) = 2.0d0  ! nu_x
            g_fact(m1_i_mu) = 1.0d0  ! mu
            g_fact(m1_i_mubar) = 1.0d0  ! mu_bar
        else
            g_fact(m1_i_mu) = 1.0d0  ! mu
            g_fact(m1_i_mubar) = 1.0d0  ! mu_bar
        end if 
      else
        if(m1_use_neutrinos) then
            g_fact(m1_i_nue) = 1.0d0  ! nu_e
            g_fact(m1_i_nuebar) = 1.0d0  ! nu_bar_e
            g_fact(m1_i_nux) = 4.0d0  ! nu_x
            ! for test:
            g_fact(4) = 1.0d0  ! nu_..
            g_fact(5) = 1.0d0  ! nu_..
        end if 
      end if 

       !> blackbody radiation in MeV:
       Black_nr_MEV = g_fact(speciesKSP) * Units * T_fluid**3 * fermi_dirac(eta,2)
       Black_er_MEV = g_fact(speciesKSP) * Units * T_fluid**4 * fermi_dirac(eta,3)

       !> blackbody radiation in code units:
       Black_nr = Black_nr_MEV * MNUC_CGS * RHOGF/TIMEGF * LENGTHGF 
       Black_er = Black_er_MEV * MEV_TO_ERG * EPSGF * RHOGF/TIMEGF * LENGTHGF 

   end subroutine blackbody_ixD

end module mod_m1_fermi
