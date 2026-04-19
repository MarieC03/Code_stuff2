module mod_fermi
  !use mod_m1_constants
  
  implicit none
  public
  integer, parameter :: NfermiDim = 10000
  double precision, parameter :: tol1 = 1.0d-16 !1.0d-8  8 digits accuracy
  double precision, parameter :: ETA_MIN = 1.0d-3
  ! eta = chem_mu_eq/Temp
    ! unit conversions 
  double precision, parameter :: LENGTHGF = 6.77269222552442d-6   ! cm to CU
  double precision, parameter :: TIMEGF   = 2.03040204956746d5    ! s to CU
  double precision, parameter :: RHOGF    = 1.61887093132742d-18  ! g/cm^3 to CU
  double precision, parameter :: PRESSGF  = 1.80123683248503d-39  ! dyn / cm^2 to CU
  double precision, parameter :: EPSGF    = 1.11265005605362d-21  ! 1. / c_cgs**2

  ! inverses
  double precision, parameter :: INVRHOGF   = 6.1771447040563d+17
  double precision, parameter :: INVEPSGF   = 8.98755178736818d20
  double precision, parameter :: INVPRESSGF = 5.55174079257738d38

  ! physical constants
  double precision, parameter :: C2_CGS        = 8.9875517873681764d+20
  double precision, parameter :: MNUC_MEV      = 931.494061
  double precision, parameter :: MNUC_CGS      = 1.660539040d-24
  double precision, parameter :: MNUC_CU       = MNUC_CGS * PRESSGF*(LENGTHGF*LENGTHGF*LENGTHGF)*C2_CGS
  double precision, parameter :: MEV_TO_ERG    = 1.60217733d-6
  double precision, parameter :: CM3_TO_FM3    = 1.0d39
  double precision, parameter :: AVOGADRO      = 6.0221367d23
  double precision, parameter :: M_NEUTRON_MEV = 939.565379
  double precision, parameter :: SIGMA_0       = 1.76d-44
  double precision, parameter :: ALPHA         = 1.23d0
  double precision, parameter :: QNP           = 1.293333
  double precision, parameter :: HC_MEVCM      = 1.23984172d-10
  double precision, parameter :: CV            = 0.5d0 + 2.0d0 * 0.23d0
  double precision, parameter :: CA            = 0.5d0
  double precision, parameter :: GAMMA_0       = 5.565d-2
  double precision, parameter :: PI            = 3.14159265358979323846
  double precision, parameter :: kB_MEV        = 8.61738568d-11

   

  !************************************************************************
  
contains
 !==================================================================
 {#IFDEF NOTUSED
 !========================= old fermi ============================== 
  function fermi_dirac_0(eta)
    double precision :: fermi_dirac_0
    double precision, intent(in) :: eta
    if ( eta > ETA_MIN) then
          fermi_dirac_0 = eta + LOG( 1.0d0 + EXP(-eta) )
       else
          fermi_dirac_0 = LOG(1.0d0 + EXP(eta) )
       end if
     end function fermi_dirac_0
  
  function fermi_dirac_old(eta,N)
    double precision :: fermi_dirac_old
    integer, intent(in) :: N
    double precision, intent(in) :: eta

    if( N .eq. 0 ) then
       if ( eta > ETA_MIN) then
          fermi_dirac_old = eta + LOG( 1.0d0 + EXP(-eta) )
       else
          fermi_dirac_old = LOG(1.0d0 + EXP(eta) )
       end if
    else if ( N .eq. 1 ) then 
       if ( eta > ETA_MIN ) then
          fermi_dirac_old = horner(eta,1.6449d0,0.0d0,0.5d0)/(1.0d0 + EXP(-1.6855d0*eta))
       else
          fermi_dirac_old = EXP(eta) / (1.0d0 + 0.2159d0 * EXP(0.8857d0*eta))
       end if

    else
       !M1_ToCheck TODO  
    end if 
  end function fermi_dirac_old
  }

 !===================================================================
 !==================== new fermi ==================================== 
 
 ! calculates fermi integral upto N=5 for negative and positive(<2) eta
  function fermi_dirac(eta,N)
    double precision :: fermi_dirac
    double precision, intent(in) :: eta
    integer, intent(in) :: N
    double precision :: FII(0:6,NfermiDim)
	double precision :: fermi
	
	    call prepare_fermi(FII)
		
      if( N .eq. 0 ) then
       if ( eta > ETA_MIN) then
	      fermi_dirac = fermi_calc(0,eta,FII)
          !fermi_dirac = eta + LOG( 1.0d0 + EXP(-eta) )
       else
	      fermi_dirac = fermi_calc(0,eta,FII) !ETA_MIN,FII)
          !fermi_dirac = LOG(1.0d0 + EXP(eta) )
       end if
    else if ( N .eq. 1 ) then 
       if ( eta > ETA_MIN ) then
	      fermi_dirac = fermi_calc(1,eta,FII)
          !fermi_dirac = horner(eta,1.6449d0,0.0d0,0.5d0)/(1.0d0 + EXP(-1.6855d0*eta))
       else
	      fermi_dirac = fermi_calc(1,eta,FII)
          !fermi_dirac = EXP(eta) / (1.0d0 + 0.2159d0 * EXP(0.8857d0*eta))
       end if
    else
	      fermi_dirac = fermi_calc(N,eta,FII)
		  ! for comparison with direct integration
		  ! call direct_integral(N,eta,fermi)
		  ! fermi_dirac = fermi
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
      do i=1,NfermiDim !5000
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
        !write(*,*) 'Argument > 1'
        return
      endif
      if(N.GT.1) then
        sum1=Z
        arg=Z
        do k=2,NfermiDim ! 10000
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
  subroutine get_neutrino_temp_ixD(wrad,N,ix^D,Gamma,velU,Wlor, eta, speciesKSP, fluid_Prim, av_energy, T_nu)    
   !> Calculates the neutrino effective temperature in MeV using 
   !> Fermi integral and average energy of neutrinos
   use mod_m1_internal
   use mod_m1_closure
   use mod_m1_eas_param !, only: idx_T
   {#IFDEF UNIT_TESTS
   use mod_m1_tests
   }
   {#IFNDEF UNIT_TESTS
   include "amrvacdef.f"
   }
   integer, intent(in) :: ix^D
   integer, intent(in) :: speciesKSP
   double precision, intent(inout) :: wrad(1:m1_numvars_internal)
   !double precision, intent(inout) :: wfluidprim_ixD(1:nw)    ! not used
   double precision, intent(inout) :: fluid_Prim(1:3)    ! not used
   !double precision, intent(inout) :: eas(1:m1_num_eas)       ! not used
   double precision, intent(in)    :: Gamma, Wlor
   double precision, intent(in)    :: velU(1:^NC) 
   double precision, intent(in)    :: N
   double precision, intent(in)    :: eta
   ! double precision, intent(in) :: wprim(ixI^S,1:nw)
   double precision, intent(out)   :: av_energy !> average neutrino energy
   double precision, intent(out)   :: T_nu   !> neutrino temperature
   ! internal:
   double precision :: E         !> neutrino energy
   double precision :: n_dens    !> number density
   double precision :: Fv        !> F_i*v^i
   double precision :: T_fluid   !> fluid temperature


     Fv  = {^C&wrad(m1_flux^C_)*velU(^C)+}

     E = wrad(m1_energy_) 
     n_dens = N/Gamma    
     if(n_dens .eq. 0.d0) then
        av_energy = 0.d0
        T_nu = 0.d0
     end if 
     if((E .eq. 0.d0) .and. (Fv .eq. 0.d0)) then
        av_energy = 0.d0
        T_nu = 0.d0
     else 
        if(n_dens < 1.d-16) then
           n_dens = 1.d-16
        endif
        !> av_energy needed in MeV thus transfer erg/cm^3 to MeV/cm^3
        E = E/MEV_TO_ERG
        Fv = Fv/MEV_TO_ERG
        !> now E[MeV/cm^3], Fv[MeV/cm^3], n_dens[1/cm^3]
        !> units of av_energy in MeV
        av_energy = Wlor*(E-Fv)/n_dens

        !> neutrino temperature: now in MeV, fermi-Dirac dimensionless
        T_nu = fermi_dirac(eta,2)/fermi_dirac(eta,3)*av_energy
     end if

     T_fluid = fluid_Prim(idx_T)
     if(T_nu > T_fluid) then
      {#IFDEF DEBUGMEERATES
        write(93,*)"Tnu>T_f",T_nu,T_fluid
        if(T_nu > 50.0d0) write(93,*)"Tnu>50Mev aiaiai"
      }
     end if 

     {#IFDEF DEBUGMEERATES
     write(93,*)"E",E 
     write(93,*)"n",n_dens
     write(93,*)"Fv",Fv
     write(93,*)"Wlor",Wlor
     write(93,*)"av_energy",av_energy  
     write(93,*)"fermi_dirac(eta,2)",fermi_dirac(eta,2) 
     write(93,*)"fermi_dirac(eta,3)",fermi_dirac(eta,3) 
     write(93,*)"T_nu",T_nu
     T_fluid = fluid_Prim(idx_T)
     if(T_nu > T_fluid) then
        write(93,*)"Tnu>T_f",T_nu,T_fluid
        if(T_nu > 50.0d0) write(93,*)"Tnu>50Mev aiaiai"
     end if 
     }

 end subroutine get_neutrino_temp_ixD 

 subroutine blackbody_ixD(wrad,ix^D,speciesKSP,fluid_Prim,T_fluid,eta,Black_er,Black_nr)
     ! use mod_m1_constants
      use mod_m1_internal 
      {#IFDEF UNIT_TESTS
      use mod_m1_tests
      }
      {#IFNDEF UNIT_TESTS
      include "amrvacdef.f"
      }
      integer, intent(in) :: ix^D, speciesKSP
      double precision, intent(in) :: eta
      double precision, intent(in) :: T_fluid 
      double precision, intent(inout) :: wrad(1:m1_numvars_internal)  ! not used
      !double precision, intent(inout) :: eas(1:m1_num_eas)      ! not used
      double precision, intent(in) :: fluid_Prim(1:3)          ! not used
      !double precision, intent(in) :: wfluidprim_ixD(1:nw)     ! not used
      double precision, intent(out) :: Black_er,Black_nr
      ! internal: 
      double precision, dimension(1:^NS) :: g_fact
      double precision :: Units, UnitsCGS_er, UnitsCGS_nr
   
      Units = 4*PI/C2_CGS/HC_MEVCM  !M1_ToCheck Units, MeV ->CGS fermi_dirc which units, T in MEV?
     
      UnitsCGS_nr = 6.593d+12 ! in 1/cm^3 (no kB) this is the prefactor 4pi/c^3/h^3
      !Q: what are the correct units for Blackbody?
      !UnitsCGS_er = 9.103d-4  ! in erg/cm^3 including kB
      !UnitsCGS_er = 568.2     ! in MeV/cm^3 including kB (kb would be in MeV/K)
      !TODO !!!!
      {#IFDEF DEBUGMEERATES
      write(93,*)"Units nr",Units,UnitsCGS_nr
      !write(93,*)"Units er",Units*kB_MEV,UnitsCGS_er
      }

       g_fact(1) = 1.0d0  ! nu_e
       if(^NS .ge. 2) then
       g_fact(2) = 1.0d0  ! nu_bar_e
       end if 
       if(^NS .ge. 3) then
       g_fact(3) = 4.0d0  ! nu_x
       end if 
       if(^NS .ge. 4) then
       g_fact(4) = 4.0d0  ! nu_x
       end if 

       !Check: use T_fluid here or T_nu neutrinos?
       Black_nr = g_fact(speciesKSP) * Units * T_fluid**3 * fermi_dirac(eta,2)
       Black_er = g_fact(speciesKSP) * Units * T_fluid**4 * fermi_dirac(eta,3)
       
       !wrong: dont need kb
       !Black_nr = g_fact(speciesKSP) * UnitsCGS_nr * T_fluid**3 * fermi_dirac(eta,2)
       !Black_er = g_fact(speciesKSP) * UnitsCGS_er * T_fluid**4 * fermi_dirac(eta,3)
      
   end subroutine blackbody_ixD

 !========================= NOT USED : ======================================
   {#IFDEF NOTUSED
 !===============================================================
  subroutine get_neutrino_temp(wprim,ixI^L,ixO^L,idim,eta,T_vu)
       use mod_m1_closure
       {#IFDEF UNIT_TESTS
       use mod_m1_tests
       }
       {#IFNDEF UNIT_TESTS
       include "amrvacdef.f"
       }
      integer, intent(in) :: ixI^L, ixO^L,idim
	   double precision, intent(in) :: eta
       double precision, intent(in) :: wprim(ixI^S,1:nw)
	   double precision, intent(out) :: T_vu(ixI^S,1:^SP) 
	   ! internal:
	   double precision, dimension(ixI^S,^NC) :: velU !M1_ToCheck
	   double precision :: av_energy,n_dens
	   
	   {^C& velU(ixO^S,^C) =  wprim(ixO^S,u0_+idim) \} !or wprim(ixO^S,u0_+^C)
	   
	   Fv   = (^C&wprim(ixO^S,frad^KSP^C_)*velU(ixO^S,^C)+)
	   
	   !{^SP&
	   !call get_av_energy(wprim,av_energy)
	     {^C& F(:,^C) = wprim(ixO^S,frad^KSP^C_)\}
         E = wprim(ixO^S,erad^KSP_) 

         !Check: devide by Gamma!
		 n_dens = wprim(ixO^S,nrad^KSP_) 
		 
	   av_energy = wprim(ixO^S,lfac_)*(E-Fv)/n_dens
	   
	   T_vu(ixO^S,^KSP) = fermi_dirac(eta,2)/fermi_dirac(eta,3)*av_energy
	   
	   !\}
	   
  end subroutine get_neutrino_temp
  
  subroutine blackbody(wprim,ixI^L,ixO^L,idim,eta,Temp,Black_er,Black_nr)
      {#IFDEF UNIT_TESTS
      use mod_m1_tests
      }
      {#IFNDEF UNIT_TESTS
      include "amrvacdef.f"
      }
      !use mod_m1_constants
       integer, intent(in) :: ixI^L, ixO^L,idim,iSP
	   double precision, intent(in) :: eta
	   double precision, intent(in) :: Temp(ixI^S,1:^SP) 
       double precision, intent(in) :: wprim(ixI^S,1:nw)
	   double precision, intent(out) :: Black_er(ixI^L,1:^SP), Black_nr(ixI^L,1:^SP)
	   ! internal: 
	   double precision, dimension(ixI^S,^NC) :: velU !M1_ToCheck
	   double precision, dimesion(1:^SP) :: g_fact
	   double precision :: av_energy,n_dens,Units
  
       Units = 4*PI/C2_CGS/HC_MEVCM  !M1_ToCheck Units, MeV ->CGS fermi_dirc which units, T in MEV?
	   
	   {^SP& g_fact(^SP) = 4.d0 \}
	   
	   !M1_ToCheck: SP indices for neutrinos and antineurinos? g_fact for nu_tau?
	   g_fact(1) = 1.d0  ! nu_e
	   g_fact(2) = 1.d0  ! nu_bar_e

      ! only done for electron (anit-)neutrios ?
	  {^SP&
       !call subroutine get_neutrino_temp(wprim,ixI^L,ixO^L,idim,eta,T_vu)
	   !Black_nr(ixO^S,^SP) = g_fact(^SP) * Units * T_vu(ixO^S,^SP)**3 * fermi_dirac(eta,2)
       !Black_er(ixO^S,^SP) = g_fact(^SP) * Units * T_vu(ixO^S,^SP)**4 * fermi_dirac(eta,3)	  
      
	  Black_nr(ixO^S,^SP) = g_fact(^SP) * Units * Temp(ixO^S,^SP)**3 * fermi_dirac(eta,2)
      Black_er(ixO^S,^SP) = g_fact(^SP) * Units * Temp(ixO^S,^SP)**4 * fermi_dirac(eta,3)
	  \}
	  
  end subroutine blackbody

  subroutine Q_emission_corr(wprim,ixI^L,ixO^L,idim,eta,Temp,Q_er_new,Q_nr_new)	 
      {#IFDEF UNIT_TESTS
      use mod_m1_tests
      }
      {#IFNDEF UNIT_TESTS
      include "amrvacdef.f"
      }
      !use mod_m1_constants
      integer, intent(in) :: ixI^L, ixO^L,idim
	   double precision, intent(in) :: eta
       double precision, intent(in) :: wprim(ixI^S,1:nw)
	   double precision, intent(in) :: Temp(ixI^L,1:^SP)
	   double precision, intent(out) :: Q_er_new(ixI^L,1:^SP), Q_nr_new(ixI^L,1:^SP)
	   	   ! internal: 
	   double precision  :: Black_er(ixI^L,1:^SP), Black_nr(ixI^L,1:^SP)
	   double precision  :: kappa_er(ixI^L,1:^SP), kappa_nr(ixI^L,1:^SP) !M1_ToCheck do I really need grid indicees?
	   double precision  :: kappa_er_EQ(ixI^L,1:^SP), kappa_nr_EQ(ixI^L,1:^SP)
	   double precision  :: T_vu(ixI^L,1:^SP)
	   double precision  :: maxterm
       
	  {^KSP&
		
       call blackbody(wprim,ixI^L,ixO^L,idim,eta,Temp,Black_er,Black_nr)
	   
	   if(^KSP.eq.1 .or. ^KSP.eq.2) then
	   	    !M1_ToCheck TODO : calculate: kappa_er_EQ, kappa_nr_EQ at equilibrium:
			
		call subroutine get_neutrino_temp(wprim,ixI^L,ixO^L,idim,eta,T_vu)
		maxterm = 1.d0
		if(T_vu(ixO^S,^SP)/Temp(ixO^S,^SP)) then
		   maxterm = (T_vu(ixO^S,^SP)/Temp(ixO^S,^SP))**2
		end if 
		kappa_er(ixO^S,^SP) = kappa_er_EQ(ixO^S,^SP) * maxterm
		kappa_nr(ixO^S,^SP) = kappa_nr_EQ(ixO^S,^SP) * maxterm
	   Q_er_new(ixO^S,^SP) = Black_er(ixO^S,^SP) * kappa_er(ixO^S,^SP)
	   Q_nr_new(ixO^S,^SP) = Black_nr(ixO^S,^SP) * kappa_nr(ixO^S,^SP)
	   else
	   	    !M1_ToCheck TODO : calculate: Q_er, Q_nr.. for kappa_x:
	   kappa_er(ixO^S,^SP) =  Q_er_new(ixO^S,^SP) / Black_er(ixO^S,^SP)
	   kappa_nr(ixO^S,^SP) =  Q_nr_new(ixO^S,^SP) / Black_nr(ixO^S,^SP)
	   end if 
	   	   
	  \}
	  
  end subroutine Q_emission_corr	
  
  }

end module mod_fermi
