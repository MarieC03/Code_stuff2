module mod_m1_eas_microphysical_gray

  implicit none
 
contains

  ! at a given (rho,T,Ye) this routine finds the equilibrium opacities
  ! A trilinear interpolation is used to find exact kappa between
  ! the saved values in weakhub tables
  subroutine m1_get_eas_microphysical_gray_Weakhub(wrad,speciesKSP,ix^D,eas_eq,fluid_Prim) 
    use mod_m1_internal
    use mod_m1_eas_param
    use mod_Weakhub_reader
    use mod_interpolate, only: TrilinearInterpolation, QuadrilinearInterpolation
    {#IFDEF UNIT_TESTS
    use mod_m1_tests
    }
    {#IFNDEF UNIT_TESTS
    ! include "amrvacdef.f"
     }
     integer, intent(in) :: ix^D, speciesKSP
     double precision, intent(in)    :: wrad(1:m1_numvars_internal) !not used
     double precision, intent(inout) :: fluid_Prim(1:fluid_vars)
     double precision, intent(inout) :: eas_eq(1:m1_num_eas)
     ! internal
     double precision :: Temp_fluid, rho_fluid, ye_fluid, ymu_fluid
     double precision :: Temp_fluid_units, Rho_fluid_units, Ye_fluid_units
     !double precision :: eta_nu
     double precision :: LogRho, LogTemp, Yeloc, LogYmuloc
     double precision :: LogRhoTab1, LogRhoTab2
     double precision :: YeTab1, YeTab2
     double precision :: LogTempTab1, LogTempTab2, LogYmuTab1, LogYmuTab2
     double precision :: f111, f211, f121, f221, f112, f212, f122, f222, result
     double precision :: f1111, f2111, f1211, f2211, f1121, f2121, f1221, f2221
     double precision :: f1112, f2112, f1212, f2212, f1122, f2122, f1222, f2222
     double precision :: UNITS_rho, UNITS_temp, UNITS_ye
     integer :: irho1, iye1, itemp1, iimu
     integer :: i,J

     rho_fluid  = fluid_Prim(idx_rho)
     ye_fluid   = fluid_Prim(idx_ye)
     Temp_fluid = fluid_Prim(idx_T)
     ymu_fluid  = 0.0d0
     if (m1_use_muons) ymu_fluid = fluid_Prim(idx_Ymu)

     !-------------
     UNITS_rho  = 1 
     UNITS_ye  = 1
     UNITS_temp = 1 
     !-------------

     !-------------- treating rho
     rho_fluid_units  = rho_fluid * UNITS_rho ! now in units of table
     LogRho = DLOG(rho_fluid_units)
     if(LogRho>=logrho_max_IV) then
        !write(99,*)"LogRho_fluid is larger than logrho_max in table"
        LogRho = logrho_max_IV * 0.9999999d0
     else if(LogRho<logrho_min_IV) then
        !write(99,*)"LogRho_fluid is smaller than logrho_min in table"
        LogRho = logrho_min_IV    
     end if 
     ! 
     do i = 1,IVrho-1
      if((LogRho >= logrho_IVtable(i)) .and. (LogRho < logrho_IVtable(i+1))) then
        irho1 = i
        LogRhoTab1 = logrho_IVtable(i)
        LogRhoTab2 = logrho_IVtable(i+1)
        !frac_rhoTabNext = (DLOG(rho_fluid_units)-logrho_IVtable(i))/(logrho_IVtable(i+1)-logrho_IVtable(i))
      end if
     end do 
     !--------------
     !-------------- treating temp
     temp_fluid_units  = temp_fluid * UNITS_temp ! now in units of table
     LogTemp = DLOG(temp_fluid_units)
     if(LogTemp>=logtemp_max_IV) then
        !write(99,*)"LogTemp_fluid is larger than logtemp_max in table"
        LogTemp = logtemp_max_IV * 0.9999999d0
     else if(LogTemp<logtemp_min_IV) then
        !write(99,*)"LogTemp_fluid is smaller than logtemp_min in table"
        LogTemp = logtemp_min_IV    
     end if 
     ! 
     do i = 1,IVtemp-1
      if((LogTemp >= logTemp_IVtable(i)) .and. (LogTemp < logtemp_IVtable(i+1))) then
        itemp1 = i
        LogTempTab1 = logtemp_IVtable(i)
        LogTempTab2 = logtemp_IVtable(i+1)
        !frac_rhoTabNext = (DLOG(rho_fluid_units)-logrho_IVtable(i))/(logrho_IVtable(i+1)-logrho_IVtable(i))
      end if
     end do 
     !-------------- treating ye
     ye_fluid_units  = ye_fluid * UNITS_ye ! now in units of table
     Yeloc =ye_fluid_units
     if(Yeloc>=ye_max_IV) then
        !write(99,*)"Ye_fluid is larger than logrho_max in table"
        Yeloc = ye_max_IV * 0.9999999d0
     else if(Yeloc<ye_min_IV) then
        !write(99,*)"Ye_fluid is smaller than logrho_min in table"
        Yeloc = ye_min_IV    
     end if 
     ! 
     do i = 1,IVye-1
      if((Yeloc >= Ye_IVtable(i)) .and. (Yeloc < ye_IVtable(i+1))) then
        iye1 = i
        YeTab1 = ye_IVtable(i)
        YeTab2 = ye_IVtable(i+1)
        !frac_rhoTabNext = (DLOG(rho_fluid_units)-logrho_IVtable(i))/(logrho_IVtable(i+1)-logrho_IVtable(i))
      end if
     end do 
     !--------------
     ! ----  do Trilinear interpolation to find rates on grid of rho,ye,temp ------
     if(m1_use_muons) then
       if (IVymu < 2) call mpistop("Muonic weakhub table needs at least two ymu points")
       LogYmuloc = DLOG(max(ymu_fluid, exp(logymu_min_IV)))
       if(LogYmuloc >= logymu_max_IV) then
          LogYmuloc = logymu_max_IV * 0.9999999d0
       else if(LogYmuloc < logymu_min_IV) then
          LogYmuloc = logymu_min_IV
       end if
       iimu = 1
       do i = 1, IVymu-1
        if((LogYmuloc >= logymu_IVtable(i)) .and. (LogYmuloc < logymu_IVtable(i+1))) then
          iimu = i
          LogYmuTab1 = logymu_IVtable(i)
          LogYmuTab2 = logymu_IVtable(i+1)
        end if
       end do
     else
       iimu = 1
     end if 

     ! get kappa_a energy rate
     if (m1_use_muons) then
       f1111 = kappa_a_en_grey_table(irho1, itemp1, iye1, iimu, speciesKSP)
       f2111 = kappa_a_en_grey_table(irho1+1, itemp1, iye1, iimu, speciesKSP)
       f1211 = kappa_a_en_grey_table(irho1, itemp1+1, iye1, iimu, speciesKSP)
       f2211 = kappa_a_en_grey_table(irho1+1, itemp1+1, iye1, iimu, speciesKSP)
       f1121 = kappa_a_en_grey_table(irho1, itemp1, iye1+1, iimu, speciesKSP)
       f2121 = kappa_a_en_grey_table(irho1+1, itemp1, iye1+1, iimu, speciesKSP)
       f1221 = kappa_a_en_grey_table(irho1, itemp1+1, iye1+1, iimu, speciesKSP)
       f2221 = kappa_a_en_grey_table(irho1+1, itemp1+1, iye1+1, iimu, speciesKSP)
       f1112 = kappa_a_en_grey_table(irho1, itemp1, iye1, iimu+1, speciesKSP)
       f2112 = kappa_a_en_grey_table(irho1+1, itemp1, iye1, iimu+1, speciesKSP)
       f1212 = kappa_a_en_grey_table(irho1, itemp1+1, iye1, iimu+1, speciesKSP)
       f2212 = kappa_a_en_grey_table(irho1+1, itemp1+1, iye1, iimu+1, speciesKSP)
       f1122 = kappa_a_en_grey_table(irho1, itemp1, iye1+1, iimu+1, speciesKSP)
       f2122 = kappa_a_en_grey_table(irho1+1, itemp1, iye1+1, iimu+1, speciesKSP)
       f1222 = kappa_a_en_grey_table(irho1, itemp1+1, iye1+1, iimu+1, speciesKSP)
       f2222 = kappa_a_en_grey_table(irho1+1, itemp1+1, iye1+1, iimu+1, speciesKSP)
       call QuadrilinearInterpolation(LogRho, LogTemp, Yeloc, LogYmuloc, LogRhoTab1, LogRhoTab2, &
              LogTempTab1, LogTempTab2, YeTab1, YeTab2, LogYmuTab1, LogYmuTab2, &
              f1111, f2111, f1211, f2211, f1121, f2121, f1221, f2221, &
              f1112, f2112, f1212, f2212, f1122, f2122, f1222, f2222, result)
     else
       f111 = kappa_a_en_grey_table( irho1, itemp1, iye1,iimu, speciesKSP)
       f211 = kappa_a_en_grey_table( irho1+1, itemp1, iye1,iimu, speciesKSP)
       f121 = kappa_a_en_grey_table( irho1, itemp1+1, iye1,iimu, speciesKSP)
       f221 = kappa_a_en_grey_table( irho1+1, itemp1+1, iye1,iimu, speciesKSP)
       f112 = kappa_a_en_grey_table( irho1, itemp1, iye1+1,iimu, speciesKSP)
       f212 = kappa_a_en_grey_table( irho1+1, itemp1, iye1+1,iimu, speciesKSP)
       f122 = kappa_a_en_grey_table( irho1, itemp1+1, iye1+1, iimu, speciesKSP)
       f222 = kappa_a_en_grey_table( irho1+1, itemp1+1, iye1+1,iimu, speciesKSP)
       call TrilinearInterpolation(LogRho, LogTemp, Yeloc, LogRhoTab1, LogRhoTab2, LogTempTab1, LogTempTab2, &
              YeTab1, YeTab2, f111, f211, f121, f221, f112, f212, f122, f222, result)
     endif
     eas_eq(k_a) = result

     ! get kappa_s energy rate
     if (m1_use_muons) then
       f1111 = kappa_s_grey_table(irho1, itemp1, iye1, iimu, speciesKSP)
       f2111 = kappa_s_grey_table(irho1+1, itemp1, iye1, iimu, speciesKSP)
       f1211 = kappa_s_grey_table(irho1, itemp1+1, iye1, iimu, speciesKSP)
       f2211 = kappa_s_grey_table(irho1+1, itemp1+1, iye1, iimu, speciesKSP)
       f1121 = kappa_s_grey_table(irho1, itemp1, iye1+1, iimu, speciesKSP)
       f2121 = kappa_s_grey_table(irho1+1, itemp1, iye1+1, iimu, speciesKSP)
       f1221 = kappa_s_grey_table(irho1, itemp1+1, iye1+1, iimu, speciesKSP)
       f2221 = kappa_s_grey_table(irho1+1, itemp1+1, iye1+1, iimu, speciesKSP)
       f1112 = kappa_s_grey_table(irho1, itemp1, iye1, iimu+1, speciesKSP)
       f2112 = kappa_s_grey_table(irho1+1, itemp1, iye1, iimu+1, speciesKSP)
       f1212 = kappa_s_grey_table(irho1, itemp1+1, iye1, iimu+1, speciesKSP)
       f2212 = kappa_s_grey_table(irho1+1, itemp1+1, iye1, iimu+1, speciesKSP)
       f1122 = kappa_s_grey_table(irho1, itemp1, iye1+1, iimu+1, speciesKSP)
       f2122 = kappa_s_grey_table(irho1+1, itemp1, iye1+1, iimu+1, speciesKSP)
       f1222 = kappa_s_grey_table(irho1, itemp1+1, iye1+1, iimu+1, speciesKSP)
       f2222 = kappa_s_grey_table(irho1+1, itemp1+1, iye1+1, iimu+1, speciesKSP)
       call QuadrilinearInterpolation(LogRho, LogTemp, Yeloc, LogYmuloc, LogRhoTab1, LogRhoTab2, &
              LogTempTab1, LogTempTab2, YeTab1, YeTab2, LogYmuTab1, LogYmuTab2, &
              f1111, f2111, f1211, f2211, f1121, f2121, f1221, f2221, &
              f1112, f2112, f1212, f2212, f1122, f2122, f1222, f2222, result)
     else
       f111 = kappa_s_grey_table( irho1, itemp1, iye1,iimu, speciesKSP)
       f211 = kappa_s_grey_table( irho1+1, itemp1, iye1,iimu, speciesKSP)
       f121 = kappa_s_grey_table( irho1, itemp1+1, iye1,iimu, speciesKSP)
       f221 = kappa_s_grey_table( irho1+1, itemp1+1, iye1,iimu, speciesKSP)
       f112 = kappa_s_grey_table( irho1, itemp1, iye1+1,iimu, speciesKSP)
       f212 = kappa_s_grey_table( irho1+1, itemp1, iye1+1,iimu, speciesKSP)
       f122 = kappa_s_grey_table( irho1, itemp1+1, iye1+1, iimu, speciesKSP)
       f222 = kappa_s_grey_table( irho1+1, itemp1+1, iye1+1,iimu, speciesKSP)
       call TrilinearInterpolation(LogRho, LogTemp, Yeloc, LogRhoTab1, LogRhoTab2, LogTempTab1, LogTempTab2, &
              YeTab1, YeTab2, f111, f211, f121, f221, f112, f212, f122, f222, result)
     endif
     eas_eq(k_s) = result 

     ! get kappa_a number rate
     if (m1_use_muons) then
       f1111 = kappa_a_num_grey_table(irho1, itemp1, iye1, iimu, speciesKSP)
       f2111 = kappa_a_num_grey_table(irho1+1, itemp1, iye1, iimu, speciesKSP)
       f1211 = kappa_a_num_grey_table(irho1, itemp1+1, iye1, iimu, speciesKSP)
       f2211 = kappa_a_num_grey_table(irho1+1, itemp1+1, iye1, iimu, speciesKSP)
       f1121 = kappa_a_num_grey_table(irho1, itemp1, iye1+1, iimu, speciesKSP)
       f2121 = kappa_a_num_grey_table(irho1+1, itemp1, iye1+1, iimu, speciesKSP)
       f1221 = kappa_a_num_grey_table(irho1, itemp1+1, iye1+1, iimu, speciesKSP)
       f2221 = kappa_a_num_grey_table(irho1+1, itemp1+1, iye1+1, iimu, speciesKSP)
       f1112 = kappa_a_num_grey_table(irho1, itemp1, iye1, iimu+1, speciesKSP)
       f2112 = kappa_a_num_grey_table(irho1+1, itemp1, iye1, iimu+1, speciesKSP)
       f1212 = kappa_a_num_grey_table(irho1, itemp1+1, iye1, iimu+1, speciesKSP)
       f2212 = kappa_a_num_grey_table(irho1+1, itemp1+1, iye1, iimu+1, speciesKSP)
       f1122 = kappa_a_num_grey_table(irho1, itemp1, iye1+1, iimu+1, speciesKSP)
       f2122 = kappa_a_num_grey_table(irho1+1, itemp1, iye1+1, iimu+1, speciesKSP)
       f1222 = kappa_a_num_grey_table(irho1, itemp1+1, iye1+1, iimu+1, speciesKSP)
       f2222 = kappa_a_num_grey_table(irho1+1, itemp1+1, iye1+1, iimu+1, speciesKSP)
       call QuadrilinearInterpolation(LogRho, LogTemp, Yeloc, LogYmuloc, LogRhoTab1, LogRhoTab2, &
              LogTempTab1, LogTempTab2, YeTab1, YeTab2, LogYmuTab1, LogYmuTab2, &
              f1111, f2111, f1211, f2211, f1121, f2121, f1221, f2221, &
              f1112, f2112, f1212, f2212, f1122, f2122, f1222, f2222, result)
     else
       f111 = kappa_a_num_grey_table( irho1, itemp1, iye1,iimu, speciesKSP)
       f211 = kappa_a_num_grey_table( irho1+1, itemp1, iye1,iimu, speciesKSP)
       f121 = kappa_a_num_grey_table( irho1, itemp1+1, iye1,iimu, speciesKSP)
       f221 = kappa_a_num_grey_table( irho1+1, itemp1+1, iye1,iimu, speciesKSP)
       f112 = kappa_a_num_grey_table( irho1, itemp1, iye1+1,iimu, speciesKSP)
       f212 = kappa_a_num_grey_table( irho1+1, itemp1, iye1+1,iimu, speciesKSP)
       f122 = kappa_a_num_grey_table( irho1, itemp1+1, iye1+1, iimu, speciesKSP)
       f222 = kappa_a_num_grey_table( irho1+1, itemp1+1, iye1+1,iimu, speciesKSP)
       call TrilinearInterpolation(LogRho, LogTemp, Yeloc, LogRhoTab1, LogRhoTab2, LogTempTab1, LogTempTab2, &
              YeTab1, YeTab2, f111, f211, f121, f221, f112, f212, f122, f222, result)
     endif
     eas_eq(k_n) = result

     !if(m1_use_muons .ne. .true.) then
     !if((speciesKSP .eq. 4) .or. (speciesKSP .eq. 5)) then
     ! eas_eq(k_a) = 0.0d0
     ! eas_eq(k_s) = 0.0d0
     ! eas_eq(k_n) = 0.0d0
     !end if 
     !end if 
     !--------------
     ! Weakhub reader should be called once at beginning of bhac
     !- ---> logtemp,logrho, logye, kappa_a,kappa_s,kappa_e arrays need to be saved/stored 

  end subroutine m1_get_eas_microphysical_gray_Weakhub

end module mod_m1_eas_microphysical_gray
