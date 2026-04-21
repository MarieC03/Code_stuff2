!program TNTextract
module fileIO

interface operator( .f. )
  module procedure file_exists
end interface

contains

function file_exists(filename) result(res)
  implicit none
  character(len=*),intent(in) :: filename
  logical                     :: res

  ! Check if the file exists
  inquire( file=trim(filename), exist=res )
end function

end module

program PhiAverages
use fileIO
       INTEGER, PARAMETER :: NDIM =400
       INTEGER, PARAMETER :: NDIM2 =150
       REAL :: OMEGA(NDIM,NDIM),Jmom(NDIM,NDIM),Rho(NDIM,NDIM)
       REAL :: Temp(NDIM,NDIM),Press(NDIM,NDIM)
	   REAL :: Gphiphi(NDIM,NDIM),Gphir(NDIM,NDIM),Gphiz(NDIM,NDIM)
	   REAL :: Vphi(NDIM,NDIM),Vr(NDIM,NDIM)
	   REAL :: BetaPhi(NDIM,NDIM),Alpha(NDIM,NDIM)
	   REAL :: Entropy(NDIM,NDIM),Bvec0(NDIM,NDIM)
	   REAL :: Bvec1(NDIM,NDIM),Bvec2(NDIM,NDIM) 
	   REAL :: Enua(NDIM,NDIM),Enue(NDIM,NDIM)
	   REAL :: Enux(NDIM,NDIM),Epsnua(NDIM,NDIM)
	   REAL :: Epsnue(NDIM,NDIM),Epsnux(NDIM,NDIM)
	   REAL :: Ye(NDIM,NDIM)
       REAL :: XGRID(NDIM),YGRID(NDIM)
       REAL :: xx,yy,r,phi,cosphi
!   new files
       REAL :: RGRID(0:NDIM2),PHIGRID(-NDIM2:NDIM2)
       REAL :: NUMP(0:NDIM2,-NDIM2:NDIM2),CO(0:NDIM2)
       REAL :: OME(0:NDIM2,-NDIM2:NDIM2),JM(0:NDIM2,-NDIM2:NDIM2)
       REAL :: TE(0:NDIM2,-NDIM2:NDIM2),PR(0:NDIM2,-NDIM2:NDIM2),PI
       REAL :: RH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: GPHIPHI2(0:NDIM2,-NDIM2:NDIM2),GPHIR2(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: GPHIZ2(0:NDIM2,-NDIM2:NDIM2),VPHI2(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: VR2(0:NDIM2,-NDIM2:NDIM2),BETAPHI2(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: ALPHA2(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: ENTROPY_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: BVEC0_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: BVEC1_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: BVEC2_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: ENUA_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: ENUE_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: ENUX_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: EPSNUA_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: EPSNUE_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: EPSNUX_PH(0:NDIM2,-NDIM2:NDIM2)	   
	   REAL :: YE_PH(0:NDIM2,-NDIM2:NDIM2)
	   REAL :: RCIRC(0:NDIM2,-NDIM2:NDIM2)
       REAL :: AVJM,AVOME,AVTE,AVPR,AVRH
	   REAL :: AVGPHIPHI,AVGPHIR,AVGPHIZ,AVVPHI,AVVZ
	   REAL :: AVBETAPHI,AVALPHA
	   Integer :: tdigit,tdezi1,tdezi2,tdezi3,tdigimax
	   Integer :: idig,jdezi,kdezi,ldezi
	   Integer :: nana,dezi
	   ! in:
	   character(150) filename,filename1,filename2,filename3
	   character(150) filename4,filename5,filename6,filename7
	   character(150) filename8,filename9,filename10,filename11
	   character(150) filename12,filename13,filename14,filename15,filename16
	   character(150) filename17,filename18,filename19,filename20,filename21
	   character(150) filename22,filename23,filename24,filename25,filename26
	   
        !out:	   
	   character(150) filename120,filename130,filename140,filename150,filename160
	   ! rphi output commented
	   ! rphi and average: 160 ---> unit=400
	   
	   tdigimax = 34
	   
	   
       do idig=0,tdigimax
	   tdigit = idig
	   do jdezi=0,9
	   tdezi1 = jdezi
	   do kdezi=0,9
	   tdezi2 = kdezi
	   do ldezi = 0,9
	   tdezi3 = ldezi
	   
		 !write(*,*) "*************************************************************"
	     !write(*,*) "--------------------",tdigit,".",tdezi1,tdezi2,tdezi3
		 !write(*,*) "-------------------------------------------------------------"
		 
	  ! if name only has 1 dezimal digit
	  if((tdezi2.eq.0).and.(tdezi3.eq.0))then
	     write(filename,'(a,i0,a,i0,a)') 'Omega_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
	  endif
	  
	  if(file_exists(filename) .and. (tdezi2.eq.0) .and. (tdezi3.eq.0)) then
	     write(*,*)"file found:",tdigit,".",tdezi1
		 write(filename1,'(a,i0,a,i0,a)') 'Jmom_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename2,'(a,i0,a,i0,a)') 'Rho_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename3,'(a,i0,a,i0,a)') 'Temp_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename4,'(a,i0,a,i0,a)') 'Press_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename5,'(a,i0,a,i0,a)') 'Gphi2_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename6,'(a,i0,a,i0,a)') 'Gphir_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename7,'(a,i0,a,i0,a)') 'Gphiz_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename8,'(a,i0,a,i0,a)') 'Vphi_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename9,'(a,i0,a,i0,a)') 'Vr_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename10,'(a,i0,a,i0,a)') 'Betaphi_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename11,'(a,i0,a,i0,a)') 'Alpha_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename12,'(a,i0,a,i0,a)') 'Entropy_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename13,'(a,i0,a,i0,a)') 'Bvec0_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename14,'(a,i0,a,i0,a)') 'Bvec1_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename15,'(a,i0,a,i0,a)') 'Bvec2_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename16,'(a,i0,a,i0,a)') 'Enua_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename17,'(a,i0,a,i0,a)') 'Enue_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename18,'(a,i0,a,i0,a)') 'Enux_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'		 
		 write(filename19,'(a,i0,a,i0,a)') 'Eps_nua_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt' 
		 write(filename20,'(a,i0,a,i0,a)') 'Eps_nue_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt' 
		 write(filename21,'(a,i0,a,i0,a)') 'Eps_nux_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
		 write(filename22,'(a,i0,a,i0,a)') 'Ye_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,'.txt'
			 		 
		 write(filename120,'(a,i0,a,i0,a)') '.\TNT_phiAv\Zout_2D_xy_TNT_',tdigit,'.',tdezi1,'00ms_fort.dat'
		 write(filename130,'(a,i0,a,i0,a)') '.\TNT_phiAv\Zout_Y0_TNT_',tdigit,'.',tdezi1,'00ms_fort.dat'
		 write(filename140,'(a,i0,a,i0,a)') '.\TNT_phiAv\Zout_X0_TNT_',tdigit,'.',tdezi1,'00ms_fort.dat'
		 write(filename150,'(a,i0,a,i0,a)') '.\TNT_phiAv\Zout_2D_rphi_TNT_',tdigit,'.',tdezi1,'00ms_fort.dat'
		 write(filename160,'(a,i0,a,i0,a)') '.\TNT_phiAv\Zout_phiAV_TNT_',tdigit,'.',tdezi1,'00ms_fort.dat'
		 
		 open(unit=3,file=filename,status="old")
         open(unit=4,file=filename1,status="old")
         open(unit=5,file=filename2,status="old")
         open(unit=7,file=filename3,status="old")
		 open(unit=8,file=filename4,status="old")
		 open(unit=9,file=filename5,status="old")
		 open(unit=10,file=filename6,status="old")
		 open(unit=11,file=filename7,status="old")
		 open(unit=12,file=filename8,status="old")
		 open(unit=13,file=filename9,status="old")
		 open(unit=14,file=filename10,status="old")
		 open(unit=15,file=filename11,status="old")
		 open(unit=16,file=filename12,status="old")
		 open(unit=17,file=filename13,status="old")
		 open(unit=18,file=filename14,status="old")
		 open(unit=19,file=filename15,status="old")
		 open(unit=20,file=filename16,status="old")
		 open(unit=21,file=filename17,status="old")
		 open(unit=22,file=filename18,status="old")
		 open(unit=23,file=filename19,status="old")
		 open(unit=24,file=filename20,status="old")
		 open(unit=25,file=filename21,status="old")
		 open(unit=26,file=filename22,status="old")
		 		 		 		 		 		 		 		 
         open(unit=200,file=filename120,status="unknown")
         open(unit=201,file=filename130,status="unknown")
         open(unit=202,file=filename140,status="unknown")
		 !open(unit=400,file=filename150,status="unknown")
         open(unit=400,file=filename160,status="unknown")

		  GOTO 757
		end if

        ! if name only has 2 dezimal digits
		if(tdezi3.eq.0) then
	     write(filename,'(a,i0,a,i0,i0,a)') 'Omega_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
	    endif
		
	  if(file_exists(filename) .and. (tdezi3.eq.0)) then
	     write(*,*)"file found:",tdigit,".",tdezi1,tdezi2		 
	     write(filename,'(a,i0,a,i0,i0,a)') 'Omega_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename1,'(a,i0,a,i0,i0,a)') 'Jmom_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename2,'(a,i0,a,i0,i0,a)') 'Rho_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename3,'(a,i0,a,i0,i0,a)') 'Temp_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename4,'(a,i0,a,i0,i0,a)') 'Press_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename5,'(a,i0,a,i0,i0,a)') 'Gphi2_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename6,'(a,i0,a,i0,i0,a)') 'Gphir_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename7,'(a,i0,a,i0,i0,a)') 'Gphiz_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename8,'(a,i0,a,i0,i0,a)') 'Vphi_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename9,'(a,i0,a,i0,i0,a)') 'Vr_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename10,'(a,i0,a,i0,i0,a)') 'Betaphi_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename11,'(a,i0,a,i0,i0,a)') 'Alpha_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename12,'(a,i0,a,i0,i0,a)') 'Entropy_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename13,'(a,i0,a,i0,i0,a)') 'Bvec0_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename14,'(a,i0,a,i0,i0,a)') 'Bvec1_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename15,'(a,i0,a,i0,i0,a)') 'Bvec2_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename16,'(a,i0,a,i0,i0,a)') 'Enua_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename17,'(a,i0,a,i0,i0,a)') 'Enue_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename18,'(a,i0,a,i0,i0,a)') 'Enux_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'		 
		 write(filename19,'(a,i0,a,i0,i0,a)') 'Eps_nua_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt' 
		 write(filename20,'(a,i0,a,i0,i0,a)') 'Eps_nue_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt' 
		 write(filename21,'(a,i0,a,i0,i0,a)') 'Eps_nux_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
		 write(filename22,'(a,i0,a,i0,i0,a)') 'Ye_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,'.txt'
			 		 
		 write(filename120,'(a,i0,a,i0,i0,a)') '.\TNT_phiAv\Zout_2D_xy_TNT_',tdigit,'.',tdezi1,tdezi2,'0ms_fort.dat'
		 write(filename130,'(a,i0,a,i0,i0,a)') '.\TNT_phiAv\Zout_Y0_TNT_',tdigit,'.',tdezi1,tdezi2,'0ms_fort.dat'
		 write(filename140,'(a,i0,a,i0,i0,a)') '.\TNT_phiAv\Zout_X0_TNT_',tdigit,'.',tdezi1,tdezi2,'0ms_fort.dat'
		 write(filename150,'(a,i0,a,i0,i0,a)') '.\TNT_phiAv\Zout_2D_rphi_TNT_',tdigit,'.',tdezi1,tdezi2,'0ms_fort.dat'
		 write(filename160,'(a,i0,a,i0,i0,a)') '.\TNT_phiAv\Zout_phiAV_TNT_',tdigit,'.',tdezi1,tdezi2,'0ms_fort.dat'
		 
		 open(unit=3,file=filename,status="old")
         open(unit=4,file=filename1,status="old")
         open(unit=5,file=filename2,status="old")
         open(unit=7,file=filename3,status="old")
		 open(unit=8,file=filename4,status="old")
		 open(unit=9,file=filename5,status="old")
		 open(unit=10,file=filename6,status="old")
		 open(unit=11,file=filename7,status="old")
		 open(unit=12,file=filename8,status="old")
		 open(unit=13,file=filename9,status="old")
		 open(unit=14,file=filename10,status="old")
		 open(unit=15,file=filename11,status="old")
		 open(unit=16,file=filename12,status="old")
		 open(unit=17,file=filename13,status="old")
		 open(unit=18,file=filename14,status="old")
		 open(unit=19,file=filename15,status="old")
		 open(unit=20,file=filename16,status="old")
		 open(unit=21,file=filename17,status="old")
		 open(unit=22,file=filename18,status="old")
		 open(unit=23,file=filename19,status="old")
		 open(unit=24,file=filename20,status="old")
		 open(unit=25,file=filename21,status="old")
		 open(unit=26,file=filename22,status="old")
		 		 		 		 		 		 		 		 
         open(unit=200,file=filename120,status="unknown")
         open(unit=201,file=filename130,status="unknown")
         open(unit=202,file=filename140,status="unknown")
		 !open(unit=400,file=filename150,status="unknown")
         open(unit=400,file=filename160,status="unknown")

		  GOTO 757
		end if
		
		! if name has 3 dezimal digits
	     write(filename,'(a,i0,a,i0,i0,i0,a)') 'Omega_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
	  
	  if(file_exists(filename)) then
	     write(*,*)"file found:",tdigit,".",tdezi1,tdezi2,tdezi3
		 write(filename,'(a,i0,a,i0,i0,i0,a)') 'Omega_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename1,'(a,i0,a,i0,i0,i0,a)') 'Jmom_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename2,'(a,i0,a,i0,i0,i0,a)') 'Rho_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename3,'(a,i0,a,i0,i0,i0,a)') 'Temp_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename4,'(a,i0,a,i0,i0,i0,a)') 'Press_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename5,'(a,i0,a,i0,i0,i0,a)') 'Gphi2_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename6,'(a,i0,a,i0,i0,i0,a)') 'Gphir_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename7,'(a,i0,a,i0,i0,i0,a)') 'Gphiz_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename8,'(a,i0,a,i0,i0,i0,a)') 'Vphi_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename9,'(a,i0,a,i0,i0,i0,a)') 'Vr_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename10,'(a,i0,a,i0,i0,i0,a)') 'Betaphi_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename11,'(a,i0,a,i0,i0,i0,a)') 'Alpha_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename12,'(a,i0,a,i0,i0,i0,a)') 'Entropy_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename13,'(a,i0,a,i0,i0,i0,a)') 'Bvec0_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename14,'(a,i0,a,i0,i0,i0,a)') 'Bvec1_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename15,'(a,i0,a,i0,i0,i0,a)') 'Bvec2_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename16,'(a,i0,a,i0,i0,i0,a)') 'Enua_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename17,'(a,i0,a,i0,i0,i0,a)') 'Enue_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename18,'(a,i0,a,i0,i0,i0,a)') 'Enux_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'		 
		 write(filename19,'(a,i0,a,i0,i0,i0,a)') 'Eps_nua_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt' 
		 write(filename20,'(a,i0,a,i0,i0,i0,a)') 'Eps_nue_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt' 
		 write(filename21,'(a,i0,a,i0,i0,i0,a)') 'Eps_nux_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
		 write(filename22,'(a,i0,a,i0,i0,i0,a)') 'Ye_2D_np_saved_TNTYST_',tdigit,'.',tdezi1,tdezi2,tdezi3,'.txt'
			 		 
		 write(filename120,'(a,i0,a,i0,i0,i0,a)') '.\TNT_phiAv\Zout_2D_xy_TNT_',tdigit,'.',tdezi1,tdezi2,tdezi3,'ms_fort.dat'
		 write(filename130,'(a,i0,a,i0,i0,i0,a)') '.\TNT_phiAv\Zout_Y0_TNT_',tdigit,'.',tdezi1,tdezi2,tdezi3,'ms_fort.dat'
		 write(filename140,'(a,i0,a,i0,i0,i0,a)') '.\TNT_phiAv\Zout_X0_TNT_',tdigit,'.',tdezi1,tdezi2,tdezi3,'ms_fort.dat'
		 write(filename150,'(a,i0,a,i0,i0,i0,a)') '.\TNT_phiAv\Zout_2D_rphi_TNT_',tdigit,'.',tdezi1,tdezi2,tdezi3,'ms_fort.dat'
		 write(filename160,'(a,i0,a,i0,i0,i0,a)') '.\TNT_phiAv\Zout_phiAV_TNT_',tdigit,'.',tdezi1,tdezi2,tdezi3,'ms_fort.dat'
		 
		 open(unit=3,file=filename,status="old")
         open(unit=4,file=filename1,status="old")
         open(unit=5,file=filename2,status="old")
         open(unit=7,file=filename3,status="old")
		 open(unit=8,file=filename4,status="old")
		 open(unit=9,file=filename5,status="old")
		 open(unit=10,file=filename6,status="old")
		 open(unit=11,file=filename7,status="old")
		 open(unit=12,file=filename8,status="old")
		 open(unit=13,file=filename9,status="old")
		 open(unit=14,file=filename10,status="old")
		 open(unit=15,file=filename11,status="old")
		 open(unit=16,file=filename12,status="old")
		 open(unit=17,file=filename13,status="old")
		 open(unit=18,file=filename14,status="old")
		 open(unit=19,file=filename15,status="old")
		 open(unit=20,file=filename16,status="old")
		 open(unit=21,file=filename17,status="old")
		 open(unit=22,file=filename18,status="old")
		 open(unit=23,file=filename19,status="old")
		 open(unit=24,file=filename20,status="old")
		 open(unit=25,file=filename21,status="old")
		 open(unit=26,file=filename22,status="old")
		 		 		 		 		 		 		 		 
         open(unit=200,file=filename120,status="unknown")
         open(unit=201,file=filename130,status="unknown")
         open(unit=202,file=filename140,status="unknown")
		 !open(unit=400,file=filename150,status="unknown")
         open(unit=400,file=filename160,status="unknown")

		  GOTO 757
		end if
		
		if(file_exists(filename).eq.(.false.)) cycle

 757    CONTINUE	 
       PI=4.*ATAN(1.)
       xmin = -20.
       xmax = 20.
       xstep = (xmax-xmin)/float(NDIM)
	   
	    DNUMP=0.0
        DOME=0.0
        DJM=0.0
        DTE=0.0
        DPR=0.0
        DRH=0.0
		DGPHIPHI=0.0
		DGPHIR=0.0
		DGPHIZ=0.0
		DVPHI=0.0
		DVR=0.0
		DBETA=0.0
		DALPHA=0.0
		DRCIRC=0.0
		
		DENTROP=0.0
		DBVEC0 = 0.0
		DBVEC1 = 0.0
		DBVEC2 = 0.0
		DENUA = 0.0
		DENUE = 0.0
		DENUX = 0.0
		DEPSNUA = 0.0
		DEPSNUE = 0.0
		DEPSNUEX = 0.0
		DYE = 0.0
		
        NUMP(:,:)=0.0
		
		
        OME(:,:)=0.0
        JM(:,:)=0.0
        TE(:,:)=0.0
        PR(:,:)=0.0
        RH(:,:)=0.0 
		GPHIPHI2(:,:)=0.0
		GPHIR2(:,:)=0.0
		GPHIZ2(:,:)=0.0
		VPHI2(:,:)=0.0
		VR2(:,:)=0.0		
		BETAPHI2(:,:)=0.0		
		ALPHA2(:,:)=0.0		
		RCIRC(:,:)=0.0
		
		ENTROPY_PH(:,:) =0.0
		BVEC0_PH(:,:) = 0.0
		BVEC1_PH(:,:) = 0.0
		BVEC2_PH(:,:) = 0.0
		ENUA_PH(:,:) = 0.0
		ENUE_PH(:,:) = 0.0
		ENUX_PH(:,:) = 0.0
		EPSNUA_PH(:,:) = 0.0
		EPSNUE_PH(:,:) = 0.0
		EPSNUX_PH(:,:) =0.0
		YE_PH(:,:) =0.0

		DO i=1,NDIM
       xx= xmin + xstep*float(i)
       XGRID(i) = xx
!       write(*,*)"------------- x(i) = ",XGRID(i)

       DO j=1,NDIM
       yy= xmin + xstep*float(j)
       YGRID(j) = yy
	   READ(3,*) Omega(i,j)
	   READ(4,*) Jmom(i,j)
	   READ(5,*) Rho(i,j)
	   READ(7,*) Temp(i,j)
	   READ(8,*) Press(i,j)
	   READ(9,*) Gphiphi(i,j)
	   READ(10,*) Gphir(i,j)
	   READ(11,*) Gphiz(i,j)
	   READ(12,*) Vphi(i,j)
	   READ(13,*) Vr(i,j)
	   READ(14,*) BetaPhi(i,j)	   
	   READ(15,*) Alpha(i,j)	   
	   READ(16,*) Entropy(i,j)   
	   READ(17,*) Bvec0(i,j)  
	   READ(18,*) Bvec1(i,j)  
	   READ(19,*) Bvec2(i,j)  
	   READ(20,*) Enua(i,j)
	   READ(21,*) Enue(i,j)
	   READ(22,*) Enux(i,j)
	   READ(23,*) Epsnua(i,j)
	   READ(24,*) Epsnue(i,j)
	   READ(25,*) Epsnux(i,j)
	   READ(26,*) Ye(i,j)
	   
	   
	   
 !       DNUMP=0.0
 !       DOME=0.0
 !       DJM=0.0
 !       DTE=0.0
 !       DPR=0.0
 !       DRH=0.0
 !       NUMP(:,:)=0.0
 !       OME(:,:)=0.0
 !       JM(:,:)=0.0
 !       TE(:,:)=0.0
 !       PR(:,:)=0.0
 !       RH(:,:)=0.0 
		
       I2=I/4
       J2=J/4
       IF((4*I2).EQ.I.AND.(4*J2).EQ.J) THEN           ! reduced output
 !      write(200,557) XGRID(i),YGRID(j),Jmom(i,j),Omega(i,j),Rho(i,j),Temp(i,j),Press(i,j),sqrt(Gphiphi(i,j))/2.0/PI,Gphiphi(i,j),&
 !	         Vphi(i,j),Vr(i,j),BetaPhi(i,j),Alpha(i,j)
	   write(200,567) XGRID(i),YGRID(j),Jmom(i,j),Omega(i,j),Rho(i,j),Temp(i,j),Press(i,j),sqrt(Gphiphi(i,j))/2.0/PI,Gphiphi(i,j),&
	         Vphi(i,j),Vr(i,j),BetaPhi(i,j),Alpha(i,j),Entropy(i,j),Bvec0(i,j),Bvec1(i,j),Bvec2(i,j),Enua(i,j),Enue(i,j),Enux(i,j),&
			 Epsnua(i,j),Epsnue(i,j),Epsnux(i,j),Ye(i,j)
       ENDIF
       !if(YGRID(j).eq.0) then
       if(j.eq.NDIM/2) then    ! slice of y=0
        !write(201,557) XGRID(i),YGRID(j),Jmom(i,j),Omega(i,j),Rho(i,j),Temp(i,j),Press(i,j),sqrt(Gphiphi(i,j))/2.0/PI,Gphiphi(i,j),Vphi(i,j),Vr(i,j),BetaPhi(i,j),Alpha(i,j)
        write(201,567) XGRID(i),YGRID(j),Jmom(i,j),Omega(i,j),Rho(i,j),Temp(i,j),Press(i,j),sqrt(Gphiphi(i,j))/2.0/PI,Gphiphi(i,j),&
	         Vphi(i,j),Vr(i,j),BetaPhi(i,j),Alpha(i,j),Entropy(i,j),Bvec0(i,j),Bvec1(i,j),Bvec2(i,j),Enua(i,j),Enue(i,j),Enux(i,j),&
			 Epsnua(i,j),Epsnue(i,j),Epsnux(i,j),Ye(i,j)
	   end if
       if(i.eq.NDIM/2) then    ! slice of x=0
       ! write(202,557) XGRID(i),YGRID(j),Jmom(i,j),Omega(i,j),Rho(i,j),Temp(i,j),Press(i,j),sqrt(Gphiphi(i,j))/2.0/PI,Gphiphi(i,j),Vphi(i,j),Vr(i,j),BetaPhi(i,j),Alpha(i,j)
       	 write(202,567) XGRID(i),YGRID(j),Jmom(i,j),Omega(i,j),Rho(i,j),Temp(i,j),Press(i,j),sqrt(Gphiphi(i,j))/2.0/PI,Gphiphi(i,j),&
	         Vphi(i,j),Vr(i,j),BetaPhi(i,j),Alpha(i,j),Entropy(i,j),Bvec0(i,j),Bvec1(i,j),Bvec2(i,j),Enua(i,j),Enue(i,j),Enux(i,j),&
			 Epsnua(i,j),Epsnue(i,j),Epsnux(i,j),Ye(i,j)
       end if
       END DO
       WRITE(200,77)
 77    FORMAT(/)
       END DO
 555   FORMAT(1X,11E15.5)
!-------------------------------------- input finished ------------------
       RCEN=0.5
!  calculate for r=0.
       DO i=1,NDIM
        DO j=1,NDIM
        xx = XGRID(i)
        yy = YGRID(j)
        r = sqrt(xx**2+yy**2)
        IF(R.GT.RCEN)                   CYCLE
        DNUMP=DNUMP+1.                    ! count entries in cell
        DOME=DOME+OMEGA(I,J)
        DJM=DJM+Jmom(I,J)
        DTE=DTE+Temp(I,J)
        DPR=DPR+Press(I,J)
        DRH=DRH+Rho(I,J)
		DGPHIPHI=DGPHIPHI+Gphiphi(I,J)
		DGPHIR=DGPHIR+Gphir(I,J)
		DGPHIZ=DGPHIZ+Gphiz(I,J)
		DVPHI=DVPHI+Vphi(I,J)
		DVR=DVR+Vr(I,J)
		DBETA=DBETA+BetaPhi(I,J)
		DALPHA=DALPHA+Alpha(I,J)
		DENTROP =DENTROP+Entropy(I,J)
		DBVEC0=DBVEC0+Bvec0(I,J)
		DBVEC1=DBVEC1+Bvec1(I,J)
		DBVEC2=DBVEC2+Bvec2(I,J)
		DENUA=DENUA+Enua(I,J)
		DENUE=DENUE+Enue(I,J)
		DENUX=DENUX+Enux(I,J)
		DEPSNUA=DEPSNUA+Epsnua(I,J)
		DEPSNUE=DEPSNUE+Epsnue(I,J)
		DEPSNUX=DEPSNUX+Epsnux(I,J)
		DYE=DYE+Ye(I,J)
		DRCIRC=DRCIRC+sqrt(Gphiphi(I,J))/2.0/PI
        END DO
       END DO
!      normalize
        DOME=DOME/DNUMP
        DJM=DJM/DNUMP
        DTE=DTE/DNUMP
        DPR=DPR/DNUMP
        DRH=DRH/DNUMP
		DGPHIPHI=DGPHIPHI/DNUMP
		DGPHIR=DGPHIR/DNUMP
		DGPHIZ=DGPHIZ/DNUMP
		DVPHI=DVPHI/DNUMP
		DVR=DVR/DNUMP
		DBETA=DBETA/DNUMP
		DALPHA=DALPHA/DNUMP
		DENTROP=DENTROP/DNUMP
		DBVEC0=DBVEC0/DNUMP
		DBVEC1=DBVEC1/DNUMP
		DBVEC2=DBVEC2/DNUMP
		DENUA=DENUA/DNUMP
		DENUE=DENUE/DNUMP
		DENUX=DENUX/DNUMP
		DEPSNUA=DEPSNUA/DNUMP
		DEPSNUE=DEPSNUE/DNUMP
		DEPSNUX=DEPSNUX/DNUMP
		DYE=DYE/DNUMP
		DRCIRC=DRCIRC/DNUMP
!---------------------------------------------------------------------

! calc phi-r distributions from cartesian grid
       NDIM3=100

       RSTEP=20./FLOAT(NDIM3)
       PHISTEP=1./FLOAT(NDIM3)*PI
       NUMP(:,:)=0.

       DO i=1,NDIM
        DO j=1,NDIM
        xx = XGRID(i)
        yy = YGRID(j)
        r = sqrt(xx**2+yy**2)
        cosphi = xx/r
        phi = acos(cosphi)
        IR=NINT(r/RSTEP)
!        IPHI=NINT(PHI/PHISTEP)
        IPHI=0                                   ! average over phi
        IF(IR.GT.NDIM3)      CYCLE
        NUMP(IR,IPHI)=NUMP(IR,IPHI)+1.                    ! count entries in cell
        OME(IR,IPHI)=OME(IR,IPHI)+OMEGA(I,J)
        JM(IR,IPHI)=JM(IR,IPHI)+Jmom(I,J)
        TE(IR,IPHI)=TE(IR,IPHI)+Temp(I,J)
        PR(IR,IPHI)=PR(IR,IPHI)+Press(I,J)
        RH(IR,IPHI)=RH(IR,IPHI)+Rho(I,J)
		GPHIPHI2(IR,IPHI)=GPHIPHI2(IR,IPHI)+Gphiphi(I,J)
		GPHIR2(IR,IPHI)=GPHIR2(IR,IPHI)+Gphir(I,J)
		GPHIZ2(IR,IPHI)=GPHIZ2(IR,IPHI)+Gphiz(I,J)
		VPHI2(IR,IPHI)=VPHI2(IR,IPHI)+Vphi(I,J)
		VR2(IR,IPHI)=VR2(IR,IPHI)+Vr(I,J)
		BETAPHI2(IR,IPHI)=BETAPHI2(IR,IPHI)+BetaPhi(I,J)
		ALPHA2(IR,IPHI)=ALPHA2(IR,IPHI)+Alpha(I,J)
		ENTROPY_PH(IR,IPHI)=ENTROPY_PH(IR,IPHI)+Entropy(I,J)
		BVEC0_PH(IR,IPHI)=BVEC0_PH(IR,IPHI)+Bvec0(I,J)
		BVEC1_PH(IR,IPHI)=BVEC1_PH(IR,IPHI)+Bvec1(I,J)
		BVEC2_PH(IR,IPHI)=BVEC2_PH(IR,IPHI)+Bvec2(I,J)
		ENUA_PH(IR,IPHI)=ENUA_PH(IR,IPHI)+Enua(I,J)
		ENUE_PH(IR,IPHI)=ENUE_PH(IR,IPHI)+Enue(I,J)		
		ENUX_PH(IR,IPHI)=ENUX_PH(IR,IPHI)+Enux(I,J)	
		EPSNUA_PH(IR,IPHI)=EPSNUA_PH(IR,IPHI)+Epsnua(I,J)
		EPSNUE_PH(IR,IPHI)=EPSNUE_PH(IR,IPHI)+Epsnue(I,J)
		EPSNUX_PH(IR,IPHI)=EPSNUX_PH(IR,IPHI)+Epsnux(I,J)
		Ye_PH(IR,IPHI)=Ye_PH(IR,IPHI)+Ye(I,J)
		! circumferential radius
		RCIRC(IR,IPHI)=RCIRC(IR,IPHI)+sqrt(Gphiphi(I,J)) !/2.0/PI
        END DO
       END DO
!
        RR0=0.                                        ! center
 !       WRITE(400,559) RR0,RR0,DJM,DOME,DRH,DTE,DPR,DNUMP
 !		WRITE(400,559) RR0,RR0,DJM,DOME,DRH,DTE,DPR,DRCIRC,DGPHIPHI,DVPHI,DVR,DBETA,DALPHA,DNUMP
		WRITE(400,569) RR0,RR0,DJM,DOME,DRH,DTE,DPR,DRCIRC,DGPHIPHI,DVPHI,DVR,DBETA,DALPHA,&
		          DENTROP,DBVEC0,DBVEC1,DBVEC2,DENUA,DENUE,DENUX,DEPSNUA,DEPSNUE,DEPSNUX,DYE,DNUMP
		
!    normalize
        CO(:)=0.
        DO I=1,NDIM3
        RR=FLOAT(I)*RSTEP                   ! r loop
        J=0
        PHIS=FLOAT(J)*PHISTEP

        IF(NUMP(I,J).NE.0.) THEN
        OME(I,J)=OME(I,J)/NUMP(I,J)
        JM(I,J)=JM(I,J)/NUMP(I,J)
        TE(I,J)=TE(I,J)/NUMP(I,J)
        PR(I,J)=PR(I,J)/NUMP(I,J)
        RH(I,J)=RH(I,J)/NUMP(I,J)
		GPHIPHI2(I,J)=GPHIPHI2(I,J)/NUMP(I,J)
		GPHIR2(I,J)=GPHIR2(I,J)/NUMP(I,J)
		GPHIZ2(I,J)=GPHIZ2(I,J)/NUMP(I,J)
		VPHI2(I,J)=VPHI2(I,J)/NUMP(I,J)
		VR2(I,J)=VR2(I,J)/NUMP(I,J)
		BETAPHI2(I,J)=BETAPHI2(I,J)/NUMP(I,J)
		ALPHA2(I,J)=ALPHA2(I,J)/NUMP(I,J)
		ENTROPY_PH(I,J)=ENTROPY_PH(I,J)/NUMP(I,J)
		BVEC0_PH(I,J)=BVEC0_PH(I,J)/NUMP(I,J)
		BVEC1_PH(I,J)=BVEC1_PH(I,J)/NUMP(I,J)
		BVEC2_PH(I,J)=BVEC2_PH(I,J)/NUMP(I,J)
		ENUA_PH(I,J)=ENUA_PH(I,J)/NUMP(I,J)
		ENUE_PH(I,J)=ENUE_PH(I,J)/NUMP(I,J)
		ENUX_PH(I,J)=ENUX_PH(I,J)/NUMP(I,J)
		EPSNUA_PH(I,J)=EPSNUA_PH(I,J)/NUMP(I,J)
		EPSNUE_PH(I,J)=EPSNUE_PH(I,J)/NUMP(I,J)
		EPSNUX_PH(I,J)=EPSNUX_PH(I,J)/NUMP(I,J)
		YE_PH(I,J)=YE_PH(I,J)/NUMP(I,J)
	    ! circumferential radius	
        RCIRC(I,J)=RCIRC(I,J)/NUMP(I,J)	 !sqrt(GPHIPHI2(I,J))/2.0/PI !RCIRC(I,J)/NUMP(I,J)		
!        WRITE(400,559) RR,PHIS,JM(I,J),OME(I,J),RH(I,J),TE(I,J),PR(I,J),NUMP(I,J)
!        WRITE(400,559) RR,PHIS,JM(I,J),OME(I,J),RH(I,J),TE(I,J),PR(I,J),RCIRC(I,J),&
!		   GPHIPHI2(I,J),VPHI2(I,J),VR2(I,J),BETAPHI2(I,J),ALPHA2(I,J),NUMP(I,J)
        WRITE(400,569) RR,PHIS,JM(I,J),OME(I,J),RH(I,J),TE(I,J),PR(I,J),RCIRC(I,J),&
		   GPHIPHI2(I,J),VPHI2(I,J),VR2(I,J),BETAPHI2(I,J),ALPHA2(I,J),ENTROPY_PH(I,J),&
		   BVEC0_PH(I,J),BVEC1_PH(I,J),BVEC2_PH(I,J),ENUA_PH(I,J),ENUE_PH(I,J),ENUX_PH(I,J),&
		   EPSNUA_PH(I,J),EPSNUE_PH(I,J),EPSNUX_PH(I,J),YE_PH(I,J),NUMP(I,J)

        CO(I)=CO(I)+1.                       ! number of nonzero cells
        ENDIF

        ENDDO
        WRITE(400,77)
        IF(NDIM3.EQ.20) WRITE(*,*) 'Files written out on unit 400'
		
	end do !dezi3
	end do !dezi2
	end do !dezi1
	end do !digit
!-------------------------------------------------------------------
 99     CONTINUE                      ! end of NDIM3 loop

 557 FORMAT(13(3XE12.5)) 
 567 FORMAT(24(3XE12.5)) 
 
 558 FORMAT(7(3XE12.5)) 
! 559 FORMAT(8(3XE12.5))
 559 FORMAT(14(3XE12.5))
 569 FORMAT(25(3XE12.5))
 
 556   FORMAT(1X,13E14.5)
       STOP
       end
