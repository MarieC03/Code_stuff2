module mod_m1

  use mod_m1_closure 
  !use mod_m1_closure, disabled => m1_update_closure 
  implicit none

  public 

  contains
  {#IFNDEF UNIT_TESTS   
  !< Startup m1 modules, make sure routines are associated.
  subroutine m1_startup()
    use mod_m1_closure, only: m1_closure_activate
    {#IFNDEF UNIT_TESTS_EXPLICIT
    {#IFNDEF M1_EXPLICIT
    use mod_m1_eas, only: m1_eas_activate,m1_eas_gray_activate
    use mod_m1_eas_test, only: m1_eas_test_activate

    use mod_Weakhub_reader
    }
    }
    {#IFDEF M1_EXPLICIT  
    include "amrvacdef.f"
    }
    call m1_closure_activate()
    {#IFDEF M1_RATES 
    ! default associate empty rates at first
     {#IFNDEF M1_RATES
     call m1_eas_test_activate(1)
     }
    ! standard gray microphysics
     {#IFDEF M1_RATES
      call m1_eas_gray_activate()
      call m1_read_Weakhub()  
      if(mype==0) write(*,*)"M1: Actived tabulated rates from Weakhub"
     }
     ! m1_eas_tests:
      {#IFDEF M1_RATES_TEST1
       call m1_eas_tests_activate(1)
      }
      {#IFDEF M1_RATES_TEST2
       call m1_eas_tests_activate(2)
      }
    ! activates the eas if m1_get_eas is associated to gray or tests
      call m1_eas_activate()
    }
    if(mype==0) write(*,*) "Rise and shine, M1 is ready!"

  end subroutine m1_startup
  }

  !-----------------------------------------------
  subroutine m1_read_Weakhub()
   use mod_Weakhub_reader
   {#IFNDEF UNIT_TESTS
   !include "amrvacdef.f"
   }
   {#IFDEF UNIT_TESTS
   character(len=512) :: fileWeakhub
   fileWeakhub = "/mnt/rafast/mcassing/Weakhub/dURCA_convenforbetaproc_DD2_grey_npe_20240723_Weakhub_grey.h5"
   }

   call read_h5_Weakhub_table_gray(fileWeakhub)

  end subroutine m1_read_Weakhub
   !-------------------------------------------

  !< Compute lambda_p/m for M1 variables. This routine is a level 1 
  !< M1 routine and is called inside of get_lambda (amrvacphys.t)
  subroutine m1_get_wavespeeds(ixI^L, ixO^L, idim, wprim, x, lambda, qtC)
    ! This routine assumes the  !
    ! primitives are up to date !
    ! ====================================================================== 
    use mod_m1_closure
	  use mod_m1_metric_interface
    ! ====================================================================== 
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    integer, intent(in) :: ixI^L, ixO^L, idim
    double precision, intent(inout) :: wprim(ixI^S,1:nw) 
	  double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: qtC
	  double precision, intent(out)   :: lambda(ixI^S, 1:ncons, 1:2)
    ! ====================================================================== 
    !internal variables
    integer :: ix^D,i,idir
    double precision, dimension(ixI^S) :: dthick, dthin
    double precision, dimension(ixI^S) :: lmthin, lmthick, lpthin, lpthick
    double precision, dimension(ixI^S) :: p_tmp, r_tmp, chi, lfact
    double precision, dimension(ixI^S) :: F_sq
    double precision, dimension(ixI^S,^NC) :: F_hi, vel
    integer :: metric_comp 
	  type(m1_metric_helper) :: metricM1  
    double precision, dimension(ixI^S,^NC) :: testa, testb
    double precision :: M1_mini_pl, M1_mini_min
    {#IFDEF UNITS_TESTS
    double precision :: TESTqtC = 1.0d+10
    }

    ! ====================================================================== 

    if(qtC .ge. TESTqtC) then 
      {#IFDEF DEBUGMEEM1WAVE 
      write(335,*)"wave start: before metric"
       }
    end if 

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wprim,x,ixI^L,ixO^L) 
    }
    if(qtC .ge. TESTqtC) then 
      {#IFDEF DEBUGMEEM1WAVE 
      write(335,*)"wave start: after metric"
       }
    end if 

    select case(idim)
    case (1)
      metric_comp = 1 
    case (2)
      metric_comp = 4 
    case (3)
      metric_comp = 6 
    end select 

    
    !###############################################################################################
    !###############################################################################################
	  {^KSP& ! Species loop 
    ! Remove normalisation from radiation moments (temporarily)
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) / metricM1%sqrtg(ixO^S)
    wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) / metricM1%sqrtg(ixO^S)
    {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) / metricM1%sqrtg(ixO^S) \}  

    if(m1_actual_speeds) then
        ! in case singularity and sqrtg = 0:
        {^D& do ix^D=ixOmin^D,ixOmax^D \}
        if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
          wprim(ix^D,nrad^KSP_) = 0.0d0
          wprim(ix^D,erad^KSP_) = 0.0d0
          {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
        end if 
        {^D& end do \}
    
        if(qtC .ge. TESTqtC) then 
        {#IFDEF DEBUGMEEM1WAVE 
         write(501,*)"wave"
        {^D& do ix^D=ixOmin^D,ixOmax^D \}
          write(335,*)"wave start: ixD",ix^D
          write(335,*)"wave start: sqrtg",metricM1%sqrtg(ix^D)
          write(335,*)"wave start: alp",metricM1%alp(ix^D)     
          write(335,*)"wave start: beta",metricM1%beta(ix^D,idim)
          {^C& write(335,*)"wave start: frad",wprim(ix^D,frad^KSP^C_)\}
        {^D& end do \}
        }
        end if
    
        if(qtC .ge. TESTqtC) then 
        {#IFDEF DEBUGMEEM1WAVE 
          write(502,*)"---------------------"
          write(502,*)"wave  wprim nrad",wprim(ixO^S,nrad^KSP_)
          write(502,*)"wave  wprim erad",wprim(ixO^S,erad^KSP_)   
          {^C& write(502,*)"wave  wprim: frad",wprim(ixO^S,frad^KSP^C_)\}
        }
        end if
    
        ! update closure 
        ! This routine in general will have reconstructed 
        ! cell face values as input, so the closure always needs to 
        ! be updated.
        call m1_update_closure(metricM1, wprim, x, ixI^L, ixO^L, ^KSP, .true., chi=chi, W=lfact, vel=vel)
    
        if(qtC .ge. TESTqtC) then 
        {#IFDEF DEBUGMEEM1WAVE 
          write(502,*)"wave 2wprim nrad",wprim(ixO^S,nrad^KSP_)
          write(502,*)"wave 2wprim erad",wprim(ixO^S,erad^KSP_)   
          {^C& write(502,*)"wave 2wprim: frad",wprim(ixO^S,frad^KSP^C_)\}
        }
        end if
    
        ! Compute F_sq 
        call metricM1%raise(ixI^L,ixO^L,wprim(ixI^S,frad^KSP1_:frad^KSP^NC_), F_hi)
    
        if(qtC .ge. TESTqtC) then 
        {#IFDEF DEBUGMEEM1WAVE  
        {^C& write(502,*)"wave 2wprim: F_hi",F_hi(ixO^S,^C)\}
        }
        end if
    
        F_sq(ixO^S)=0.0d0
        {^C& 
        F_sq(ixO^S) =F_sq(ixO^S) + wprim(ixO^S,frad^KSP^C_)*F_hi(ixO^S,^C)
        \}
    
        if(qtC .ge. TESTqtC) then 
        {#IFDEF DEBUGMEEM1WAVE 
        {^D& do ix^D=ixOmin^D,ixOmax^D \}
           write(404,*)"wave start: ixD",ix^D
          {^C& write(404,*)"wave after closure: frad",wprim(ix^D,frad^KSP^C_)\}
          {^C& write(404,*)"wave after closure: F_hi",F_hi(ix^D,^C)\}
            write(404,*)"wave after closure: F_sq",F_sq(ix^D)
        {^D& end do \}
        }
        end if
        ! ====================================================================== 
        ! lambda_m^(thin)
        lmthin(ixO^S) = -metricM1%beta(ixO^S,idim) &
        - metricM1%alp(ixO^S)*dabs(F_hi(ixO^S,idim))/dsqrt(F_sq(ixO^S)+1.d-30)
        ! lambda_p^(thin)
        lpthin(ixO^S) = -metricM1%beta(ixO^S,idim) &
        + metricM1%alp(ixO^S)*dabs(F_hi(ixO^S,idim))/dsqrt(F_sq(ixO^S)+1.d-30)
        ! Auxiliaries for thick lambdas 
        p_tmp(ixO^S) = metricM1%alp(ixO^S) * vel(ixO^S,idim)/lfact(ixO^S)
        r_tmp(ixO^S) = dsqrt( metricM1%alp(ixO^S)*metricM1%alp(ixO^S) * metricM1%gammaUPij(ixO^S,metric_comp) * &
             (2.0d0*lfact(ixO^S)*lfact(ixO^S) + 1.0d0) - 2*lfact(ixO^S)*lfact(ixO^S)*p_tmp(ixO^S)*p_tmp(ixO^S))  
        ! lambda_m^(thick)
        lmthick(ixO^S) = MIN(-metricM1%beta(ixO^S,idim)+p_tmp(ixO^S), &
             -metricM1%beta(ixO^S,idim)+ (2.d0*p_tmp(ixO^S)*lfact(ixO^S)*lfact(ixO^S) - r_tmp(ixO^S))/(2.0d0*lfact(ixO^S)*lfact(ixO^S)+1.0d0) )
        ! lambda_p^(thick) 
        lpthick(ixO^S) = MAX(-metricM1%beta(ixO^S,idim)+p_tmp(ixO^S), &
             -metricM1%beta(ixO^S,idim)+ (2.0d0*p_tmp(ixO^S)*lfact(ixO^S)*lfact(ixO^S) + r_tmp(ixO^S))/(2.0d0*lfact(ixO^S)*lfact(ixO^S)+1.0d0) )
        ! ====================================================================== 
        dthin(ixO^S) = 1.5d0*chi(ixO^S) - 0.5d0
        dthick(ixO^S) = 1.0d0-dthin(ixO^S)

        M1_mini_pl = 1.0d-30
        M1_mini_min = 1.0d-30
    
        ! Lambda for n, F are identical! 
        ! _____________ Store lambda + ____________________
        lambda(ixO^S,erad^KSP_,1) = dthick(ixO^S)*lpthick(ixO^S) + dthin(ixO^S)*lpthin(ixO^S) + M1_mini_pl
        lambda(ixO^S,nrad^KSP_,1) = dthick(ixO^S)*lpthick(ixO^S) + dthin(ixO^S)*lpthin(ixO^S)+ M1_mini_pl
        {^C& lambda(ixO^S,frad^KSP^C_,1) = dthick(ixO^S)*lpthick(ixO^S) + dthin(ixO^S)*lpthin(ixO^S) + M1_mini_pl \} 
        ! ______________ Store lambda - ______________________
        lambda(ixO^S,erad^KSP_,2) = dthick(ixO^S)*lmthick(ixO^S) + dthin(ixO^S)*lmthin(ixO^S) + M1_mini_min
        lambda(ixO^S,nrad^KSP_,2) = dthick(ixO^S)*lmthick(ixO^S) + dthin(ixO^S)*lmthin(ixO^S) + M1_mini_min
        {^C& lambda(ixO^S,frad^KSP^C_,2) = dthick(ixO^S)*lmthick(ixO^S) + dthin(ixO^S)*lmthin(ixO^S) + M1_mini_min \} 
        ! ======================================================================
    !++++++++++++++ Speeds as in Radice et al +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if(m1_radice_speeds) then
        ! _____________ Store lambda + ____________________
        lambda(ixO^S,erad^KSP_,1) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) + metricM1%beta(ixO^S,1)) + M1_mini_pl
        lambda(ixO^S,nrad^KSP_,1) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) + metricM1%beta(ixO^S,1)) + M1_mini_pl
        {^C& lambda(ixO^S,frad^KSP^C_,1) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) + metricM1%beta(ixO^S,1)) + M1_mini_pl \} 
        ! ______________ Store lambda - ______________________
        lambda(ixO^S,erad^KSP_,2) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) - metricM1%beta(ixO^S,1)) + M1_mini_min
        lambda(ixO^S,nrad^KSP_,2) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) - metricM1%beta(ixO^S,1)) + M1_mini_min
        {^C& lambda(ixO^S,frad^KSP^C_,2) = dabs( metricM1%alp(ixO^S) * sqrt(metricM1%gammaUPij(ixO^S,1)) - metricM1%beta(ixO^S,1)) + M1_mini_min\} 
        ! ====================================================================== 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else 
      lambda(ix^D,nrad^KSP_,1) = 1.0d0
      lambda(ix^D,erad^KSP_,1) = 1.0d0
      {^C& lambda(ix^D,frad^KSP^C_,1) = 1.0d0 \}
      lambda(ix^D,nrad^KSP_,2) = - 1.0d0
      lambda(ix^D,erad^KSP_,2) = - 1.0d0
      {^C& lambda(ix^D,frad^KSP^C_,2) = - 1.0d0 \}
    end if 
    !++++++++++++++ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    ! rescale back rad. variables to evolved ones
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) * metricM1%sqrtg(ixO^S) 
    wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) * metricM1%sqrtg(ixO^S) 
    {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) *metricM1%sqrtg(ixO^S) \}

    ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
     if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
       wprim(ix^D,nrad^KSP_) = 0.0d0
       wprim(ix^D,erad^KSP_) = 0.0d0
       {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
       lambda(ix^D,nrad^KSP_,1) = 0.577350269189625d0
       lambda(ix^D,erad^KSP_,1) = 0.577350269189625d0
       {^C& lambda(ix^D,frad^KSP^C_,1) = 0.577350269189625d0 \}
       lambda(ix^D,nrad^KSP_,2) = - 0.577350269189625d0
       lambda(ix^D,erad^KSP_,2) = - 0.577350269189625d0
       {^C& lambda(ix^D,frad^KSP^C_,2) = - 0.577350269189625d0 \}
     end if 
    {^D& end do \}

    ! in case wavespeeds too close to zero:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
     if(lambda(ix^D,nrad^KSP_,1) .le. 1.d-4 .and. lambda(ix^D,nrad^KSP_,2) .ge. -1.d-4) then
       lambda(ix^D,nrad^KSP_,1) = 1.0d0
       lambda(ix^D,erad^KSP_,1) = 1.0d0
       {^C& lambda(ix^D,frad^KSP^C_,1) = 1.0d0 \}
       lambda(ix^D,nrad^KSP_,2) = - 1.0d0
       lambda(ix^D,erad^KSP_,2) = - 1.0d0
       {^C& lambda(ix^D,frad^KSP^C_,2) = - 1.0d0 \}
     end if 
    {^D& end do \}

    \} ! end ^ KSP
    ! ====================================================================== 
    !###############################################################################################
    !###############################################################################################


    
    if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1WAVE 
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    write(335,*)"**********************"
    write(335,*)"wave: ixD",ix^D
    write(335,*)"wave: x",x(ix^D,1)
    write(335,*)"wave: y",x(ix^D,2)
    write(335,*)"wave: idim",idim

    write(335,*)"wave: sqrtg",metricM1%sqrtg(ix^D)
    write(335,*)"wave: alp",metricM1%alp(ix^D)     
    write(335,*)"wave: beta",metricM1%beta(ix^D,idim) 
    write(335,*)"---------------"
    write(335,*)"wave: chi dthin dthick",chi(ix^D),dthin(ix^D),dthick(ix^D)
    write(335,*)"wave: Fsq",F_sq(ix^D)
    {^C& write(335,*)"wave: F_hi",F_hi(ix^D,^C) \}
    write(335,*)"wave: dabs(F_hi)/dsqrt(F**2 + 1d-30)",dabs(F_hi(ix^D,idim))/dsqrt(F_sq(ix^D)+1.d-30)
    write(335,*)"wave: dabs(F_hi)",dabs(F_hi(ix^D,idim))
    write(335,*)"wave: 1/dsqrt(F**2 + 1d-30)",1.d0/dsqrt(F_sq(ix^D)+1.d-30)
    write(335,*)"wave: 1/dsqrt(F**2 + 1d-15)",1.d0/dsqrt(F_sq(ix^D)+1.d-15)
    write(335,*)"wave: 1d-30/dsqrt(F**2 + 1d-10)",1.d-30/dsqrt(F_sq(ix^D)+1.d-30)
    write(335,*)"wave: lfact vel",lfact(ix^D),vel(ix^D,idim)
    write(335,*)"wave: ptmp rtmp",p_tmp(ix^D),r_tmp(ix^D)
    write(335,*)"wave: Lpthin Lpthick",lpthin(ix^D),lpthick(ix^D)
    write(335,*)"wave: Lmthin Lmthick",lmthin(ix^D),lmthick(ix^D)
    write(335,*)"---------------"
    write(335,*)"wave: nrad",wprim(ix^D,nrad^KSP_)
    write(335,*)"wave: erad",wprim(ix^D,erad^KSP_)
    write(335,*)"wave: frad1",wprim(ix^D,frad^KSP1_)
    write(335,*)"wave: frad2",wprim(ix^D,frad^KSP2_)
    write(335,*)"wave: frad3",wprim(ix^D,frad^KSP3_)
    write(335,*)"wave: lambda+",lambda(ix^D,erad^KSP_,1)
    write(335,*)"wave: lambda-",lambda(ix^D,erad^KSP_,2)
     {^D& end do \}
     }
    end if
      
    call metricM1%destroy()

  end subroutine m1_get_wavespeeds
  
  subroutine m1_get_fluxes(ixI^L, ixO^L, idim, x, wprim, flux, qtC)
    use mod_m1_closure
	  use mod_m1_metric_interface
    ! ====================================================================== 
    {#IFNDEF UNIT_TESTS
    include "amrvacdef.f"
    }
    ! ====================================================================== 
    integer, intent(in) :: ixI^L, ixO^L, idim
	  double precision, intent(in)  :: qtC
	  double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout)  :: wprim(ixI^S,1:nw) 
    double precision, intent(out) :: flux(ixI^S,1:nwflux)
    ! ====================================================================== 
    ! internal variables
    ! ====================================================================== 
    integer iw, ix^D, i, idir
    double precision, dimension(ixI^S,^NC) :: F, Hup, vel
    double precision, dimension(ixI^S,m1_npress) :: P
    double precision, dimension(ixI^S)    :: Gamma, Jradi, E, Wlor
    double precision, dimension(ixI^S,^NC) :: Pimu
    double precision, dimension(ixI^S,^NC) ::  Pud
	  double precision :: zetaL
    {#IFDEF UNITS_TESTS
    double precision :: TESTqtC = 1.0d+10
    }
	  logical :: get_zeta
    logical :: check_atmo_in_getflux = .false.
    !integer, dimension(3,3) :: idimtoimu  
    ! ====================================================================== 
	  type(m1_metric_helper)   :: metricM1 
    ! ====================================================================== 

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wprim,x,ixI^L,ixO^L)
    }

    if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1FLUX 
      write(501,*)"flux"
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
      write(335,*)"flux start: ixD",ix^D
      write(335,*)"flux start: sqrtg",metricM1%sqrtg(ix^D)
      write(335,*)"flux start: alp",metricM1%alp(ix^D)     
      write(335,*)"flux start: beta",metricM1%beta(ix^D,idim)
    {^D& end do \}
    }
  end if
	  
    {^KSP& ! Species loop 
    ! Remove normalisation 
    ! IMPORTANT: always needs to be added back at the end of 
    ! the routine. 
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) / metricM1%sqrtg(ixO^S)
    wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) / metricM1%sqrtg(ixO^S)
    {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) / metricM1%sqrtg(ixO^S) \}   
    
     ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}  
    ! update closure to get E,Fi,zeta
    call m1_update_closure(metricM1, wprim, x, ixI^L, ixO^L,^KSP,.true.,Gamma,Hup,Jradi,P,W=Wlor,vel=vel)

    ! Remove Gamma from n 
    wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_)/Gamma(ixO^S) 
    ! get energy
    E(ixO^S)  = wprim(ixO^S,erad^KSP_)
    ! Pimu = ( PUPxx, PUPxy, PUPxz ) if idim = 1 
    ! Pimu = ( PUPxy, PUPyy, PUPyz ) if idim = 2
    ! Pimu = ( PUPxz, PUPyz, PUPzz ) if idim = 3 
    {^C& Pimu(ixO^S,^C) = P(ixO^S,metricM1%idimtoimu(idim,^C) ) \} 

    if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1FLUX
   write(401,*)"Pimu 1",Pimu(ixO^S,1) 
   write(401,*)"Pimu 2",Pimu(ixO^S,2) 
   write(401,*)"Pimu 3",Pimu(ixO^S,3) 
   write(401,*)"----------------"
   }
    end if

    ! Shuffle some indices up and down cause why not 
    call metricM1%lower(ixI^L,ixO^L,Pimu, Pud) !ixO
		call metricM1%raise(ixI^L,ixO^L,wprim(ixI^S,frad^KSP1_:frad^KSP^NC_), F) 

    if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1FLUX
    write(401,*)"Pud 1",Pud(ixO^S,1) 
    write(401,*)"Pud 2",Pud(ixO^S,2) 
    write(401,*)"Pud 3",Pud(ixO^S,3) 
    write(401,*)"----------------"

    write(402,*)"wprim f 1",wprim(ixO^S,frad^KSP1_) 
    write(402,*)"wprim f 2",wprim(ixO^S,frad^KSP2_) 
    write(402,*)"wprim f 3",wprim(ixO^S,frad^KSP3_) 
    write(402,*)"F 1",F(ixO^S,1) 
    write(402,*)"F 2",F(ixO^S,2) 
    write(402,*)"F 3",F(ixO^S,3) 
    write(402,*)"----------------"
    }
  end if
    
    ! Get fluxes 
    ! ====================================================================== 
    ! Nflux^idim = sqrtg * alp * n * (( W*v^idim-beta^idim/alp) + H^idim/Jrad )
    ! See Radice, Bernuzzi et al. Eqs. (23) and (28)
    flux(ixO^S,nrad^KSP_) = metricM1%sqrtg(ixO^S) * metricM1%alp(ixO^S) * wprim(ixO^S,nrad^KSP_) &
         * (Wlor(ixO^S)*( vel(ixO^S,idim) - metricM1%beta(ixO^S,idim)/metricM1%alp(ixO^S)) &
         + Hup(ixO^S,idim)/(Jradi(ixO^S) + M1_TINY))  
    ! ====================================================================== 
    ! Eflux^idim  = sqrtg * ( alp * F^idim - beta^idim* E )
    flux(ixO^S,erad^KSP_) = metricM1%sqrtg(ixO^S) * ( metricM1%alp(ixO^S) * F(ixO^S,idim) &
    - metricM1%beta(ixO^S,idim) * E(ixO^S) )
    ! ====================================================================== 
    ! Fflux^idim_i  = sqrtg * ( alp * P^idim_i - beta^idim * F_i )
    {^C&flux(ixO^S,frad^KSP^C_) = metricM1%sqrtg(ixO^S) * ( metricM1%alp(ixO^S) * Pud(ixO^S,^C) &
    - metricM1%beta(ixO^S,idim) * wprim(ixO^S,frad^KSP^C_) ) \}
    ! ====================================================================== 
    ! ====================================================================== 
    ! ====================================================================== 
    ! NB here is where we re-densitize the radiation vars 
 		 wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_)*metricM1%sqrtg(ixO^S)*Gamma(ixO^S)
		 wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_)*metricM1%sqrtg(ixO^S)
		 {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_)*metricM1%sqrtg(ixO^S) \}     
     
     ! in case singularity and sqrtg = 0:
     {^D& do ix^D=ixOmin^D,ixOmax^D \}
     if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
     flux(ix^D,nrad^KSP_) = 0.0d0
     flux(ix^D,erad^KSP_) = 0.0d0
     {^C& flux(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}


     if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1FLUX 
     {^D& do ix^D=ixOmin^D,ixOmax^D \}
     write(335,*)"flux: ixD",ix^D
     write(335,*)"flux: x",x(ix^D,1)
     write(335,*)"flux: y",x(ix^D,2)
     write(335,*)"flux: sqrtg",metricM1%sqrtg(ix^D)
     write(335,*)"flux: alp",metricM1%alp(ix^D)     
     write(335,*)"flux: beta",metricM1%beta(ix^D,idim)
     write(335,*)"flux: Hup, Jrad",Hup(ix^D,idim),Jradi(ix^D) 
     write(335,*)"flux: u0+idim",wprim(ix^D,u0_+idim)
     write(335,*)"flux: Pud",Pud(ix^D,idim)
     write(335,*)"flux: nrad",wprim(ix^D,nrad^KSP_)
     write(335,*)"flux: erad",wprim(ix^D,erad^KSP_)
     write(335,*)"flux: frad1",wprim(ix^D,frad^KSP1_)
     write(335,*)"flux: frad2",wprim(ix^D,frad^KSP2_)
     write(335,*)"flux: frad3",wprim(ix^D,frad^KSP3_)
     write(335,*)"flux: -- fc nrad",flux(ix^D,nrad^KSP_)
     write(335,*)"flux: -- fc erad",flux(ix^D,erad^KSP_)
     write(335,*)"flux: -- fc frad1",flux(ix^D,frad^KSP1_)
     write(335,*)"flux: -- fc frad2",flux(ix^D,frad^KSP2_)
     write(335,*)"flux: -- fc frad3",flux(ix^D,frad^KSP3_)
    {^D& end do \}
       }
     end if

         !---------------------------------------
    ! is set to false:
    ! check atmosphere
    if(check_atmo_in_getflux) then
      {^D& do ix^D=ixImin^D,ixImax^D \}
      ! check wprim
      if(wprim(ix^D,erad^KSP_) .lt. m1_E_atmo) then
        wprim(ix^D,nrad^KSP_) = m1_E_atmo
        wprim(ix^D,erad^KSP_) = m1_E_atmo
      end if 
      {^C&
      if((wprim(ix^D,frad^KSP^C_) .lt. m1_E_atmo) .and. (wprim(ix^D,frad^KSP^C_) .gt. -1.0d0*m1_E_atmo)) then
         wprim(ix^D,frad^KSP^C_) = m1_E_atmo 
         if(wprim(ix^D,frad^KSP^C_) .gt. -1.0d0*m1_E_atmo) wprim(ix^D,frad^KSP^C_) = -1.0d0*m1_E_atmo
      end if 
      \}
      ! check flux
      if(flux(ix^D,erad^KSP_) .lt. m1_E_atmo) then
        flux(ix^D,nrad^KSP_) = m1_E_atmo
        flux(ix^D,erad^KSP_) = m1_E_atmo
      end if 
      {^C&
      if((flux(ix^D,frad^KSP^C_) .lt. m1_E_atmo) .and. (flux(ix^D,frad^KSP^C_) .gt. -1.0d0*m1_E_atmo)) then
         flux(ix^D,frad^KSP^C_) = m1_E_atmo 
      end if 
      \}
      {^D& end do \}
   endif
   !---------------------------------------

    \} ! End of species loop 
    call metricM1%destroy()
    ! ====================================================================== 
    ! ====================================================================== 
    !                               All done! 
    ! ====================================================================== 
    ! ====================================================================== 
  end subroutine m1_get_fluxes

  !< Apply correction to hll flux 
  !< to preserve asymptotic (diffusion)
  !< limit of RTE
  subroutine m1_correct_asymptotic_fluxes(ix^D, ixI^L, tvdlfeps1, dxdim, idims, wLC, wRC, fC, fLC, fRC, cminC, cmaxC, wradimpl, qtC)
  {#IFDEF UNIT_TESTS
  use mod_m1_tests
  }
  {#IFNDEF UNIT_TESTS
  include "amrvacdef.f"
  }
  integer, intent(in) :: ix^D, ixI^L, idims  
  double precision, intent(in) :: tvdlfeps1, qtC
  double precision, intent(in), dimension(1:ndim)                   :: dxdim        !< grid spacing 
  double precision, intent(in), dimension(ixI^S,1:nwflux)           :: fLC, fRC     !< left/right fluxes 
  double precision, intent(inout), dimension(ixI^S,1:nwflux,1:ndim) :: fC           !< flux 
  double precision, intent(in), dimension(ixI^S,1:nw)               :: wLC, wRC     !< left/right state 
  double precision, intent(in), dimension(ixI^S,1:ncons)            :: cminC, cmaxC !< min / max wavespeeds 
  double precision, intent(in), dimension(ixI^S,1:nm1rad_eas)           :: wradimpl     !old state 

  ! internals 
  double precision :: A 
  {#IFDEF UNITS_TESTS
  double precision :: TESTqtC = 1.0d+10
  }
  double precision :: fc_loc
  integer          :: iw, ixR^D, i,j

  {#IFDEF UNIT_TESTS
  !Kronecker delta and Levi-Civita tensors
  INTEGER:: kr(0:3,0:3),lvc(1:3,1:3,1:3)
  }

  {#IFDEF UNIT_TESTS    
    ! Kronecker delta defined in mod_m1_tests    
   do i= 0,3
    do j=0,3
       kr(i,j)=0            
       if(i.eq.j) then
        kr(i,j)=1
       end if
    end do 
   end do 
  }
  {^D& ixR^D = ix^D + kr(idims,^D)\}
  
  {^KSP& 
    

  if(qtC .ge. TESTqtC) then 
  {#IFDEF DEBUGMEEM1CORRFLUX 
    write(335,*) "kappa",wradimpl(ix^D,kappa_a^KSP_)
    write(335,*) "kappa",wradimpl(ixR^D,kappa_a^KSP_)
    write(335,*) "ixD, ixR",ix^D,ixR^D
    write(335,*) "dxdim",dxdim(idims)
    write(335,*) "---------"
    write(335,*) "cmax erad",cmaxC(ix^D,erad^KSP_)
    write(335,*) "cmin erad",cminC(ix^D,erad^KSP_)
    
    write(335,*) "cmax nrad",cmaxC(ix^D,erad^KSP_)
    write(335,*) "cmin nrad",cminC(ix^D,erad^KSP_)
    
    write(335,*) "cmax frad1",cmaxC(ix^D,frad^KSP1_)
    write(335,*) "cmin frad1",cminC(ix^D,frad^KSP1_)
    write(335,*) "cmax frad2",cmaxC(ix^D,frad^KSP2_)
    write(335,*) "cmin frad2",cminC(ix^D,frad^KSP2_)
    write(335,*) "cmax frad3",cmaxC(ix^D,frad^KSP3_)
    write(335,*) "cmin frad3",cminC(ix^D,frad^KSP3_)
    write(335,*) "---------"
    write(335,*) "fC hll nrad",fC(ix^D,nrad^KSP_,idims)
    write(335,*) "fC hll erad",fC(ix^D,erad^KSP_,idims)
    write(335,*) "fC hll frad1",fC(ix^D,frad^KSP1_,idims)
    write(335,*) "fC hll frad2",fC(ix^D,frad^KSP2_,idims)
    write(335,*) "fC hll frad3",fC(ix^D,frad^KSP3_,idims)
     }
  end if 

   !write(92,*) wradimpl(ix^D,kappa_a^KSP_)

    ! Compute A factor 
    A = MIN(1.0d0, dabs(4.0d0/ dabs(dxdim(idims)) /(wradimpl(ixR^D,kappa_a^KSP_) + wradimpl(ix^D,kappa_a^KSP_) &
                        + wradimpl(ixR^D,kappa_s^KSP_) + wradimpl(ix^D,kappa_s^KSP_) + 1.0d-40)))

     if(A .ne. A) then
       A = 1.0d0
     end if 
     if(A .gt. 1.0d0) then 
       A = 1.0d0
     end if
    
    if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1CORRFLUX
    write(335,*) "corr flux: A",A
    write(325,*) "corr flux: -----idim",idims
    write(325,*) "corr flux: dx",dxdim(idims)
    write(325,*) "corr flux: kappaA_i",(wradimpl(ix^D,kappa_a^KSP_))
    write(325,*) "corr flux: kappaA_iR",(wradimpl(ixR^D,kappa_a^KSP_))
    write(325,*) "corr flux: kappaA_Av",(wradimpl(ixR^D,kappa_a^KSP_) + wradimpl(ix^D,kappa_a^KSP_))/2.0
    write(325,*) "corr flux: kappaS_i",(wradimpl(ix^D,kappa_a^KSP_))
    write(325,*) "corr flux: kappaS_iR",(wradimpl(ixR^D,kappa_a^KSP_))
    write(325,*) "corr flux: kappaS_Av",(wradimpl(ixR^D,kappa_a^KSP_) + wradimpl(ix^D,kappa_a^KSP_))/2.0
    write(325,*) "corr flux: argument",dabs(4.0d0/ dabs(dxdim(idims)) /(wradimpl(ixR^D,kappa_a^KSP_) + wradimpl(ix^D,kappa_a^KSP_) &
    + wradimpl(ixR^D,kappa_s^KSP_) + wradimpl(ix^D,kappa_s^KSP_) + 1.0d-40))
    write(325,*) "corr flux: A",A
    }
    end if 

    !-------------------- erad
    ! overwrite HLL like this else in fc already is cmax*cmin/(cmax-xmin)*(WR-wL) without A
    !------------------------                    
    iw = erad^KSP_    
    if(m1_erad_LLF) then
      call mpistop("m1 LLF not implemented yet")
    else if(m1_parabolic) then
    ! to test diffusion limit and equations get parabolic
    fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw))&
                  /(cmaxC(ix^D,iw)-cminC(ix^D,iw))    
    else
      fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw)&
      +A*tvdlfeps1*cminC(ix^D,iw)&
      *cmaxC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
      /(cmaxC(ix^D,iw)-cminC(ix^D,iw)) 
    end if 
    !-------------------- nrad
    ! overwrite HLL like this else in fc already is cmax*cmin/(cmax-xmin)*(WR-wL) without A
    !------------------------ 
    iw = nrad^KSP_ 
    if(m1_parabolic) then
      ! test diffusion limit:
      fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw))&
                    /(cmaxC(ix^D,iw)-cminC(ix^D,iw))
    else 
      fC(ix^D, iw,idims) = (cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw)&
                  +A*tvdlfeps1*cminC(ix^D,iw)&
                  *cmaxC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
                  /(cmaxC(ix^D,iw)-cminC(ix^D,iw))
    end if 

     !-------------------- frad-C
    ! overwrite HLL like this else in fc already is cmax*cmin/(cmax-xmin)*(WR-wL) without A
    !------------------------   
    {^C&
    iw = frad^KSP^C_ 
    if(m1_frad_A2_A3) then
       fC(ix^D, iw,idims) = A**2*(cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw)&
       +A*tvdlfeps1*cminC(ix^D,iw)&
       *cmaxC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))&
       /(cmaxC(ix^D,iw)-cminC(ix^D,iw))
   
       fC(ix^D,iw,idims) = fC(ix^D, iw,idims) +0.5d0 *(1.0d0-A**2)* ( fLC(ix^D,iw) + fRC(ix^D,iw) )
    else if(m1_frad_A2_A) then
      fC(ix^D, iw,idims) = (A**2*(cmaxC(ix^D,iw)*fLC(ix^D,iw)-cminC(ix^D,iw)*fRC(ix^D,iw))&
                          +A*tvdlfeps1*cmaxC(ix^D,iw)*cminC(ix^D,iw)*(wRC(ix^D,iw)-wLC(ix^D,iw)))/&
                          (cmaxC(ix^D,iw)-cminC(ix^D,iw)) + 0.5d0*(1.0d0-A**2)*(fLC(ix^D,iw)+fRC(ix^D,iw))
    else if(m1_frad_LLF) then
      call mpistop("LLF for frad not implemented yet")
    else if(m1_parabolic) then
       !diffusive lim
       fC(ix^D, iw,idims) =  0.5d0 * ( fLC(ix^D,iw) + fRC(ix^D,iw) )
    else 
      fC(ix^D, iw,idims) = A**2 * fC(ix^D, iw,idims) + 0.5d0 * (1.0d0 - A**2) * ( fLC(ix^D,iw) + fRC(ix^D,iw) )
    end if 
    \}

    !---------------------------------------
    ! check atmosphere
    ! check flux
    !!if(fC(ix^D,erad^KSP_,idims) .lt. m1_E_atmo) then
    !!  fC(ix^D,nrad^KSP_,idims) = 0.0d0 !m1_E_atmo
    !!  fC(ix^D,erad^KSP_,idims) = 0.0d0 !m1_E_atmo
    !!end if 
    !!{C
    !!if((fC(ix^D,frad^KSP^C_,idims) .lt. m1_E_atmo) .and. (fC(ix^D,frad^KSP^C_,idims) .gt. -1.0d0*m1_E_atmo)) then
    !!  fC(ix^D,frad^KSP^C_,idims) = 0.0d0 !m1_E_atmo 
    !!end if 
    !!}
    !------------------
    ! check for NAN :
    {#IFDEF M1_FLUX_NAN_CHECK
    if(fC(ix^D,nrad^KSP_,idims) .ne. fC(ix^D,nrad^KSP_,idims)) then
      write(777,*) "---hll: NAN detected: FC(nrad)", fC(ix^D,nrad^KSP_,idims)
      fC(ix^D,nrad^KSP_,idims) = 0.0d0
    end if 
    if(fC(ix^D,erad^KSP_,idims) .ne. fC(ix^D,erad^KSP_,idims)) then
      write(777,*) "---hll: NAN detected: FC(erad)", fC(ix^D,erad^KSP_,idims)
      fC(ix^D,erad^KSP_,idims) = 0.0d0
    end if 
    {^C&
    if(fC(ix^D,frad^KSP^C_,idims) .ne. fC(ix^D,frad^KSP^C_,idims))then
      write(777,*) "---hll: NAN detected: FC(frad)", fC(ix^D,frad^KSP^C_,idims)
      fC(ix^D,frad^KSP^C_,idims) = 0.0d0
    end if
    \}
    }
    !---------------------------------------
    if(qtC .ge. TESTqtC) then  
    {#IFDEF DEBUGMEEM1CORRFLUX        
    write(335,*)"corr flux: ixD",ix^D
    write(335,*)"corr flux: fc nrad",fC(ix^D,nrad^KSP_,idims)
    write(335,*)"corr flux: fc erad",fC(ix^D,erad^KSP_,idims)
    write(335,*)"corr flux: fc frad1",fC(ix^D,frad^KSP1_,idims)
    write(335,*)"corr flux: fc frad2",fC(ix^D,frad^KSP2_,idims)
    write(335,*)"corr flux: fc frad3",fC(ix^D,frad^KSP3_,idims)
      }
    end if 

  \}   

  end subroutine m1_correct_asymptotic_fluxes

  ! ====================================================================== 
  ! ====================================================================== 
  !< compute geometrical sources for M1
  subroutine m1_add_geometrical_sources(ixI^L,ixO^L,x,wprim,wcons, qdt,d_k_gamma_ij,Kij_local,qtC)
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    use mod_m1_closure 
	  use mod_m1_metric_interface
    {#IFNDEF UNIT_TESTS
    use mod_metric
    include "amrvacdef.f"
    }
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, qtC
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
	  double precision, intent(inout)    :: Kij_local(ixI^S,^NC,^NC)
	  double precision, intent(inout)    :: d_k_gamma_ij(ixI^S,^NC,^NC,^NC)
    double precision, intent(inout) :: wcons(ixI^S, 1:nw) 
	  double precision, intent(inout) :: wprim(ixI^S, 1:nw)
    ! internal
    integer :: inonzero, i, j, pressind, ix^D
    double precision :: symfact
    double precision, dimension(ixI^S) :: rhs
    double precision, dimension(ixI^S,m1_npress) :: P
    double precision, dimension(ixI^S,^NC) :: Pimu
    double precision, dimension(ixI^S,^NC) :: FU
	  type(m1_metric_helper) :: metricM1 
    {#IFDEF UNITS_TESTS
    double precision :: TESTqtC = 1.0d+10
    }
    
    !===============================

    {#IFDEF UNIT_TESTS2
    call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wprim,x,ixI^L,ixO^L)
    }
    
     if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1GEO 
      write(501,*)"geo"
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    write(335,*)"addgeo: ixD",ix^D
    write(335,*)"before addgeo: wcons nrad",wcons(ix^D,nrad^KSP_)
    write(335,*)"before addgeo: wcons erad",wcons(ix^D,erad^KSP_)
    write(335,*)"before addgeo: wcons frad1",wcons(ix^D,frad^KSP1_)
    write(335,*)"before addgeo: wcons frad2",wcons(ix^D,frad^KSP2_)
    write(335,*)"before addgeo: wcons frad3",wcons(ix^D,frad^KSP3_)
    {^D& end do \}
    }
    {#IFDEF DEBUGMEEM1GEO 
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
      write(335,*)"addgeo start: sqrtg",metricM1%sqrtg(ix^D)
      write(335,*)"addgeo start: alp",metricM1%alp(ix^D)     
      {^C& write(335,*)"addgeo start: beta",metricM1%beta(ix^D,^C) \}
    {^D& end do \}
    }
     end if

	{^KSP&
     rhs(ixO^S) = 0.0d0
     
     ! multiply with sqrtg at end of add_geometrical_sources
 		 wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_) &
     / metricM1%sqrtg(ixO^S)
		 wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_) / metricM1%sqrtg(ixO^S)
		 {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) &
     / metricM1%sqrtg(ixO^S) \}   

    ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}


    ! update closure to get E,Fi,zeta and pressure
    call m1_update_closure(metricM1, wprim, x, ixI^L, ixO^L,^KSP,.true.,Press=P) 

    if(qtC .ge. TESTqtC) then 
      {#IFDEF DEBUGMEEM1GEO 
      {^D& do ix^D=ixOmin^D,ixOmax^D \}
      do i = 1, m1_npress        
	     write(1003,*)"addgeo: -- ixD",ix^D 
	     write(1003,*)"addgeo: ipress",i
	     write(1003,*)"addgeo: press",P(ix^D,i) 
      end do
      {^D& end do \}
      }
      end if

    !------------  N-density -------------
	  ! geom-source is zero:
	  ! wcons(ixO^S,nrad^KSP_) =  wcons(ixO^S,nrad^KSP_) +qdt*0.0d0*metricM1%sqrtg(ixO^S)

    ! ------------ Energy ------------------------------
    ! -- 
    ! Source term is:
    ! sqrtg [ alp PK - F^i partial_i alp ]
    ! --
      {#IFDEF M1_FLAT_KIJ
      Kij_local(ixO^S,:,:) = 0.0d0
      }
      {#IFDEF M1_FLAT_DKDGIJ
       d_k_gamma_ij(ixO^S,:,:,:)=0.0d0
      }

    rhs(ixO^S) = Kij_local(ixO^S,1,1)*P(ixO^S,1) &
         + Kij_local(ixO^S,1,2)*P(ixO^S,2) &
         + Kij_local(ixO^S,1,3)*P(ixO^S,3) &
         + Kij_local(ixO^S,2,2)*P(ixO^S,4) &
         + Kij_local(ixO^S,2,3)*P(ixO^S,5) &
         + Kij_local(ixO^S,3,3)*P(ixO^S,6) & 
         + Kij_local(ixO^S,2,1)*P(ixO^S,2) & 
         + Kij_local(ixO^S,3,1)*P(ixO^S,3) & 
         + Kij_local(ixO^S,3,2)*P(ixO^S,5) 
   
    rhs(ixO^S) = rhs(ixO^S) * metricM1%alp(ixO^S)

    call metricM1%raise(ixI^L,ixO^L,wprim(ixI^S,frad^KSP1_:frad^KSP^NC_),FU)   


    if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1GEO 
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
	  write(403,*)"addgeo: -- ixD",ix^D 
      write(403,*)"addgeo: Flow frad",wprim(ix^D,frad^KSP1_),wprim(ix^D,frad^KSP2_),wprim(ix^D,frad^KSP3_)
      write(403,*)"addgeo: Fup frad",FU(ix^D,1),FU(ix^D,2),FU(ix^D,3)    
    {^D& end do \}
    }
    end if

	  do i = 1,^NC !M1_ToCheck +-FU?
       rhs(ixO^S) = rhs(ixO^S) - FU(ixO^S,i)* metricM1%dalp(ixO^S,i)  
    end do 
	
    wcons(ixO^S,erad^KSP_) = wcons(ixO^S,erad^KSP_) + qdt*rhs(ixO^S)*metricM1%sqrtg(ixO^S)

    ! ------------ Fluxes ----------------------------
    ! --
    ! Source term is:
    ! sqrtgamma [ F_i partial_k beta^i -
    ! E partial_k alp + alp/2*P^{ij} partial_k gamma_{ij} ]
    ! --
    
    {^C&
    ! first set all to zero if all metric derivs are zero 
    rhs(ixO^S) = 0.0d0
    ! - E partial alp
       rhs(ixO^S) = -1.0d0 * wprim(ixO^S,erad^KSP_) * metricM1%dalp(ixO^S,^C)
       if(qtC .ge. TESTqtC) then 
       {#IFDEF DEBUGMEEM1GEO2
       write(335,*)"-----C------",^C
       write(335,*)"addgeo: metric dalp",metricM1%dalp(ixO^S,^C)
       write(335,*)"addgeo: wp erad",wprim(ixO^S,erad^KSP_)
       }
       end if

    ! F_i partial_k beta^i
    do i=1, ^NC
       rhs(ixO^S) = rhs(ixO^S) + wprim(ixO^S,frad^KSP1_+i-1) * metricM1%dbeta(ixO^S,i,^C) !^C,i)
       if(qtC .ge. TESTqtC) then 
       {#IFDEF DEBUGMEEM1GEO2
       write(335,*)"addgeo: metric dbeta i C",metricM1%dbeta(ixO^S,i,^C)
       write(335,*)"addgeo: metric w:frad i C",wprim(ixO^S,frad^KSP1_+i-1)
       write(335,*)"addgeo: rhs",rhs(ixO^S)
       }
       end if
    end do

    pressind = 1
    symfact  = 1.0d0
    do i= 1, ^NC
      if(qtC .ge. TESTqtC) then
      {#IFDEF DEBUGMEEM1GEO2
      write(335,*)"-----i------",i
      }
      end if
       do j=i, ^NC
          if(qtC .ge. TESTqtC) then
          {#IFDEF DEBUGMEEM1GEO2
          write(335,*)"-----j------",j
          write(335,*)"addgeo: rhs",rhs(ixO^S)
          }
          end if
          if ( i .ne. j) then
             symfact = 2.0d0
          else
             symfact = 1.0d0
          end if
          rhs(ixO^S) = rhs(ixO^S) + 0.5d0 * metricM1%alp(ixO^S) * symfact &
          * P(ixO^S,pressind) * d_k_gamma_ij(ixO^S,i,j,^C)
          pressind = pressind + 1 
          if(qtC .ge. TESTqtC) then
          {#IFDEF DEBUGMEEM1GEO2
          write(335,*)"addgeo: alpha",metricM1%alp(ixO^S)
          write(335,*)"addgeo: d_k_gamma_ij",d_k_gamma_ij(ixO^S,i,j,^C)
          write(335,*)"addgeo: Press",P(ixO^S,pressind-1)
          write(335,*)"addgeo: rhs 2 ",rhs(ixO^S)
          }
          end if 
       end do
    end do
    wcons(ixO^S,frad^KSP^C_) = wcons(ixO^S,frad^KSP^C_) + qdt*rhs(ixO^S)*metricM1%sqrtg(ixO^S)
    \}   
    
  !-------------------------
     ! rescale back wprim with sqrtg 
 		 wprim(ixO^S,nrad^KSP_) = wprim(ixO^S,nrad^KSP_)*metricM1%sqrtg(ixO^S)
		 wprim(ixO^S,erad^KSP_) = wprim(ixO^S,erad^KSP_)*metricM1%sqrtg(ixO^S)
		 {^C& wprim(ixO^S,frad^KSP^C_) = wprim(ixO^S,frad^KSP^C_) *metricM1%sqrtg(ixO^S) \}   
     
    ! in case singularity and sqrtg = 0:
    {^D& do ix^D=ixOmin^D,ixOmax^D \}
    if(metricM1%sqrtg(ix^D).eq. 0.0d0) then
     wprim(ix^D,nrad^KSP_) = 0.0d0
     wprim(ix^D,erad^KSP_) = 0.0d0
     {^C& wprim(ix^D,frad^KSP^C_) = 0.0d0 \}
    end if 
    {^D& end do \}

     if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1GEO 

    {^D& do ix^D=ixOmin^D,ixOmax^D \}

    write(335,*)"addgeo: ixD",ix^D
    write(335,*)"addgeo: sqrtg",metricM1%sqrtg(ix^D)
    write(335,*)"addgeo: alp",metricM1%alp(ix^D)     
    {^C& write(335,*)"addgeo: beta",metricM1%beta(ix^D,^C) \}
    write(335,*)"addgeo: wcons nrad",wcons(ix^D,nrad^KSP_)
    write(335,*)"addgeo: wcons erad",wcons(ix^D,erad^KSP_)
    write(335,*)"addgeo: wcons frad1",wcons(ix^D,frad^KSP1_)
    write(335,*)"addgeo: wcons frad2",wcons(ix^D,frad^KSP2_)
    write(335,*)"addgeo: wcons frad3",wcons(ix^D,frad^KSP3_)

    write(335,*)"addgeo: ixD",ix^D
    write(335,*)"addgeo: wprim nrad",wprim(ix^D,nrad^KSP_)
    write(335,*)"addgeo: wprim erad",wprim(ix^D,erad^KSP_)
    write(335,*)"addgeo: wprim frad1",wprim(ix^D,frad^KSP1_)
    write(335,*)"addgeo: wprim frad2",wprim(ix^D,frad^KSP2_)
    write(335,*)"addgeo: wprim frad3",wprim(ix^D,frad^KSP3_)
    {^D& end do \}
    }
     end if

    \} ! end ^KSP
    ! All done
    call metricM1%destroy()

    !===============================

  end subroutine m1_add_geometrical_sources

  !****************************************************************************

  {#IFNDEF UNIT_TESTS_EXPLICIT
  {#IFNDEF M1_EXPLICIT
  subroutine m1_add_collisional_sources(dtfactor,qdt,qtC,x,wcons,wold,wradimpl,ixI^L)
    ! wcons = psa%w%w out, at return it will be =U+F+S+G
    ! wold = psb%w%w old state without implicit part, at return = U+F+S
    ! wradimpl = psm1%%pwrad%wradimpl

    use mod_m1_internal
    use mod_m1_collisional, only: m1_get_implicit_collisional_sources
    use mod_m1_eas
    use mod_eos, only: small_rho, small_temp
    {#IFNDEF UNIT_TESTS
    use mod_m1_metric_interface
    include "amrvacdef.f"
    }
    integer, intent(in) :: ixI^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(in)    :: dtfactor   !< Timestep factor 
    double precision, intent(inout) :: wcons(ixI^S,1:nw)
    double precision, intent(inout) :: wold(ixI^S,1:nw)
    double precision, intent(in) :: wradimpl(ixI^S,1:nm1rad_eas)
    double precision, intent(in) :: qdt,qtC
    !internal
    integer :: ix^D, ixO^L
    double precision :: N
    double precision :: eas_ixD(1:m1_num_eas)
    double precision, dimension(1:m1_numvars_internal) :: wrad
    double precision, dimension(1:m1_numvars_internal+1,^NS) :: sources
    double precision, dimension(1:^NC) :: vel_
    integer, dimension(1:6) :: signum
    type(m1_metric_helper) :: metricM1 	

    ixO^L=ixI^L^LSUBdixB;

	  !need wprim here for metric!
    {#IFDEF UNITT_TESTS2
     call fill_metric(metricM1)
    }
    {#IFNDEF UNIT_TESTS2
	  call metricM1%fill_metric(wold,x,ixI^L,ixO^L)  
    }
    
    if( ^NS .eq. 4) then
      signum(1) = 1
      signum(2) = -1
      signum(3) = 0
      signum(4) = 0
    else if( ^NS .eq. 3) then
       signum(1) = 1
       signum(2) = -1
       signum(3) = 0
    else if(^NS .eq. 2) then     
      signum(1) = 1
      signum(2) = -1
    else 
       signum(1) = 1 
    end if     

    if(qtC .ge. TESTqtC) then 
    {#IFDEF DEBUGMEEM1COLL      
    write(335,*)"m1_add_collisional_sources"
    }
    end if 

    {^KSP&
    {do ix^D=ixOmin^D,ixOmax^D \}

    if(qtC .ge. TESTqtC) then   
      {#IFDEF DEBUGMEEM1COLL      
      write(335,*)"m1coll: ixD",ix^D
      }
    end if 

    !! ps2 = wcons, ps1 = wold
    !! Ut_explicit has to be wold
    N = wold(ix^D,nrad^KSP_)  
    wrad(m1_energy_) = wold(ix^D,erad^KSP_)
    {^C& wrad(m1_flux^C_) = wold(ix^D,frad^KSP^C_) \}
    {^C& vel_(^C) = wold(ix^D,u0_+^C) \}

    if(qtC .ge. TESTqtC) then   
      {#IFDEF DEBUGMEEM1COLL      
      write(335,*)"m1coll: wradimpl(Q_er)",wradimpl(ix^D,Q_er^KSP_)
      write(335,*)"m1coll: wradimpl(Q_nr)",wradimpl(ix^D,Q_nr^KSP_)
      write(335,*)"m1coll: wradimpl(k_a)",wradimpl(ix^D,kappa_a^KSP_)
      write(335,*)"m1coll: wradimpl(k_s)",wradimpl(ix^D,kappa_s^KSP_)
      write(335,*)"m1coll: wradimpl(k_n)",wradimpl(ix^D,kappa_nr^KSP_)
      }
    end if 

    eas_ixD(Q_ems) = wradimpl(ix^D,Q_er^KSP_)
    eas_ixD(k_a) = wradimpl(ix^D,kappa_a^KSP_)
    eas_ixD(k_s) = wradimpl(ix^D,kappa_s^KSP_)
    eas_ixD(Q_ems_n) = wradimpl(ix^D,Q_nr^KSP_)
    eas_ixD(k_n) = wradimpl(ix^D,kappa_nr^KSP_)

    if(qtC .ge. TESTqtC) then   
      {#IFDEF DEBUGMEEM1COLL      
      write(335,*)"m1coll: before m1_get_implicit_collisional"
      }
    end if 

    call m1_get_implicit_collisional_sources(wrad,metricM1,vel_,eas_ixD,&
         ixI^L,ix^D,N,sources(:,^KSP),qdt,^KSP,dtfactor)

         if(qtC .ge. TESTqtC) then   
          {#IFDEF DEBUGMEEM1COLL      
          write(335,*)"m1coll: after m1_get_implicit_collisional"
          }
        end if 

    if(M1_FLUID_BACKREACT) then
      if(wcons(ix^D,rho_) > small_rho .and. wcons(ix^D,T_eps_) > small_temp ) then
       ! NOTE: qdt*dtfactor
       ! update fluid, m1 backreaction on fluid
       wcons(ix^D,tau_) = wcons(ix^D,tau_) - qdt*dtfactor*sources(1,^KSP)
       {^C& wcons(ix^D,s^C_) = wcons(ix^D,s^C_) - qdt*dtfactor*sources(1+^C,^KSP) \}
       {#IFDEF TABEOS
       wcons(ix^D,ye_) = wcons(ix^D,ye_) - qdt*dtfactor * signum(^KSP) * sources(1+^NC+1,^KSP) 
       }
      end if 
    end if  ! end M1_BAKCREACT

    if(qtC .ge. TESTqtC) then   
      {#IFDEF DEBUGMEEM1COLL      
      write(335,*)"m1coll: after m1_fluid_backreact"
      }
    end if 

    ! update rad
    wcons(ix^D, nrad^KSP_) = N
    wcons(ix^D, erad^KSP_) = wrad(m1_energy_)
    {^C& wcons(ix^D, frad^KSP^C_) = wrad(m1_flux^C_) \}

    if(qtC .ge. TESTqtC) then   
      {#IFDEF DEBUGMEEM1COLL      
      write(335,*)"m1coll: nrad",wcons(ix^D, nrad^KSP_)
      write(335,*)"m1coll: erad",wcons(ix^D, erad^KSP_)
      write(335,*)"m1coll: frad1",wcons(ix^D, frad^KSP1_)
      write(335,*)"m1coll: frad2",wcons(ix^D, frad^KSP2_)
      write(335,*)"m1coll: frad3",wcons(ix^D, frad^KSP3_)
      }
    end if 

    {enddo ^D&\}
    \} ! end KSP
    
    call metricM1%destroy()


    
  end subroutine m1_add_collisional_sources
   } ! end IFNDEF M1_EXPLICIT
   } ! end IFNDEF UNIT_TESTS 
  
end module mod_m1
