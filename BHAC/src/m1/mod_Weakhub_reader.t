! n_ebin = n_groups and n_spec = n_species when creating tables
! Grey table is following Federico 2023
! weakhub_table%energies(0:weakhub_table%n_groups+1)), 0 index: -ve buffer zone, weakhub_table%n_groups+1 index: outer buffer zone
! logrho: log_e(rho);  logT: log_e(T)
! use_emissivity_table --> output emissivities eta which computed by numerical integration 
! kappa_a_st_numu  is nu_x when not using muonic
module mod_Weakhub_reader
  use mod_weakhub_parameters
  use HDF5
  {#IFNDEF UNIT_TESTS
  include 'amrvacdef.f' 
  }
  !implicit none
  public
 
  ! #########################################################
  !  Table only for stellarcollapse, not for Compose yet
  !  UNITS:
  !  all Kernels and Opacities should be in cgs units
  !  Kernels table store the values directly without log or 1 over sth
  !  all absorption opacity is corrected by stimulated absorption: kappa_a_st
  ! #########################################################

  integer, save :: n_total_nunuPP, n_total_PP, n_total_IS
  character(len=512) :: weakhub_table_filename, Table_info_filename, base, date, time
  ! ---- local ----
  ! ---- Grey tables ----, assumed eta_nu = eta_nu_eq for neutrinos, following Federico 2023
  double precision, allocatable :: kappa_a_en_grey_table(:,:,:,:,:)     ! rho, T, ye, ymu, n_spec
  double precision, allocatable :: kappa_a_num_grey_table(:,:,:,:,:)    ! rho, T, ye, ymu, n_spec
  double precision, allocatable :: kappa_s_grey_table(:,:,:,:,:)        ! rho, T, ye, ymu, n_spec

  ! ---- Spectral tables ----
  ! absorption opacities
  double precision, allocatable :: kappa_a_st_table_nue(:,:,:,:,:)       ! rho, T, ye, ymu, n_ebin
  double precision, allocatable :: kappa_a_st_table_nue_bar(:,:,:,:,:)   ! rho, T, ye, ymu, n_ebin
  ! numu kappa_a table could be nux for three species scheme !
  double precision, allocatable :: kappa_a_st_table_numu(:,:,:,:,:)      ! rho, T, ye, ymu, n_ebin  
  double precision, allocatable :: kappa_a_st_table_numu_bar(:,:,:,:,:)  ! rho, T, ye, ymu, n_ebin
  double precision, allocatable :: kappa_a_st_table_nutautaubar(:,:,:,:,:)  ! rho, T, ye, ymu, n_ebin
  ! Scattering opacities
  double precision, allocatable :: kappa_s_table_nu(:,:,:,:,:,:)        ! rho, T, ye, ymu, n_ebin, n_spec
  ! Emissivities
  double precision, allocatable :: eta_table_nue(:,:,:,:,:)       ! rho, T, ye, ymu, n_ebin
  double precision, allocatable :: eta_table_nue_bar(:,:,:,:,:)   ! rho, T, ye, ymu, n_ebin
  double precision, allocatable :: eta_table_numu(:,:,:,:,:)      ! rho, T, ye, ymu, n_ebin  
  double precision, allocatable :: eta_table_numu_bar(:,:,:,:,:)  ! rho, T, ye, ymu, n_ebin
  double precision, allocatable :: eta_table_nutautaubar(:,:,:,:,:)  ! rho, T, ye, ymu, n_ebin

  double precision, allocatable :: IS_Phi0_out_table(:,:,:) ! T, eta_e, n_total_IS
  double precision, allocatable :: IS_Phi1_out_table(:,:,:) ! T, eta_e, n_total_IS
  double precision, allocatable :: PP_Phi0_ann_table(:,:,:) ! T, eta_e, n_total_PP 
  double precision, allocatable :: PP_Phi1_ann_table(:,:,:) ! T, eta_e, n_total_PP 
  double precision, allocatable :: PP_Phi0_pro_table(:,:,:) ! T, eta_e, n_total_PP 
  double precision, allocatable :: PP_Phi1_pro_table(:,:,:) ! T, eta_e, n_total_PP 

  ! For NNbrem and de-excitation, etc 
  double precision, allocatable :: PP_Phi0_ann_table_ThreeD(:,:,:,:) ! rho, T, yp, n_total_PP 
  double precision, allocatable :: PP_Phi1_ann_table_ThreeD(:,:,:,:) ! rho, T, yp, n_total_PP 
  ! For nunubar pair process
  double precision, allocatable :: nunuPP_Phi0_ann_table(:,:,:,:) ! rho, T, yp, n_total_nunuPP 
  double precision, allocatable :: nunuPP_Phi1_ann_table(:,:,:,:) ! rho, T, yp, n_total_nunuPP 

  ! 4D table for opacity
  double precision, allocatable :: logtemp_IVtable(:) ! in ln(MeV)
  double precision, allocatable :: logrho_IVtable(:)  ! in ln(code unit)
  double precision, allocatable :: ye_IVtable(:)      ! if no muonic table using, ye = yp => ye_Ktable = yp_table
  double precision, allocatable :: logymu_IVtable(:)  ! in ln(ymu)
  ! 3D table for opacity
  double precision, allocatable :: logtemp_IIItable(:) ! in ln(MeV)
  double precision, allocatable :: logrho_IIItable(:)  ! in ln(code unit)
  double precision, allocatable :: yp_IIItable(:)      ! Proton fraction
  ! 2D table for Kernel
  double precision, allocatable :: logtemp_IItable(:) ! in ln(MeV)
  double precision, allocatable :: logeta_e_IItable(:)! in ln(eta_e)

  !.....................................
  logical :: use_muonic_table 
  !integer :: mype
  integer :: n_ebin
  !.....................................
  
contains


  !###########################################################################
  !###########################################################################
  !###########################################################################
  subroutine read_h5_Weakhub_table_gray(filename)
    !use mod_weakhub_parameters
    !use HDF5 !hdf5
    implicit none
    character(len=512), intent(in) :: filename
    character(len=512) :: energy_bin_type
    integer :: error, rank, cerror
    integer(HID_T) :: file_id, dset_id, dspace_id,i,j,k,index
    integer(HSIZE_T) :: dims1(1), dims3(3), dims4(4), dims5(5), dims6(6), dims1_char(1)
    character(8) :: date
    integer :: values(8), itemp, irho, iye, iymu, i_spec, i_ebin
    ! Lookup table
    double precision :: x2_min(2), x2_max(2)
    double precision :: x3_min(3), x3_max(3)
    double precision :: x4_min(4), x4_max(4)
   
    !......................
    use_muonic_table =.false.
    n_ebin = 1
    !......................

     cerror = 0
     call h5open_f(error)
     if (error.ne.0) then
        {#IFNDEF UNIT_TESTS
        call mpistop("Error reading in weakhub file")
        }
        {#IFDEF UNIT_TESTS
          write(99,*)"Error reading in weakhub file"
        }
     end if
     call h5fopen_f(trim(adjustl(filename)), &
          H5F_ACC_RDONLY_F,file_id,error)
     if (error.ne.0) then
        if (mype == 0) write(*,*) trim(adjustl(filename))
        {#IFNDEF UNIT_TESTS
        call mpistop("Error reading in weakhub table")
        }
        {#IFDEF UNIT_TESTS
        write(99,*)"Error reading in weakhub table"
        }
     end if

     rank = 1
     dims1(1) = 1
 
     ! 1. Read the basic info of the table
     call h5dopen_f(file_id, "weakhub_table_type",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, weakhub_table_type, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     select case(weakhub_table_type)
     case(1) ! all D
       use_IVtable  = .true.
       use_IItable  = .true.
       use_IIItable = .true.
       if (mype == 0) write(*,*) 'weakhub_table_type : table_allD'
     case(2) ! 2D+4D
       use_IVtable  = .true.
       use_IItable  = .true.
       use_IIItable = .false.
       if (mype == 0) write(*,*) 'weakhub_table_type : table_IVandII'
     case(3) ! 4D 
       use_IVtable  = .true.
       use_IItable  = .false.
       use_IIItable = .false.
       if (mype == 0) write(*,*) 'weakhub_table_type : table_IV'
     case(4) ! 2D 
       use_IVtable  = .false.
       use_IItable  = .true.
       use_IIItable = .false.
       if (mype == 0) write(*,*) 'weakhub_table_type : table_II'
     case(5) ! 3D
       use_IVtable  = .false.
       use_IItable  = .false.
       use_IIItable = .true.
       if (mype == 0) write(*,*) 'weakhub_table_type : table_III'
     case(100)
       use_IVtable  = .true.
       use_IItable  = .false.
       use_IIItable = .false.
       if (mype == 0) write(*,*) 'weakhub_table_type : table_grey'
       if (mype == 0) write(*,*) "file",filename
       if (weakhub_grey_method .ne. 1) then
        {#IFNDEF UNIT_TESTS
         call mpistop('Please set weakhub_grey_method = 1 &
                                          when reading grey table')
        }
        {#IFDEF UNIT_TESTS
        write(99,*)'Please set weakhub_grey_method = 1 &
                                          when reading grey table'
        }
       end if 
     case default
       {#IFNDEF UNIT_TESTS
       call mpistop('Reading Weakhub table which has no weakhub_table_type')
       }
       {#IFDEF UNIT_TESTS
        write(99,*) 'Reading Weakhub table which has no weakhub_table_type'
       }
     end select

     if (mype == 0) then
       if (weakhub_table_type == table_grey) then
         write(*,*) 'Grey table only read the 4D table'
       endif
     endif

     call h5dopen_f(file_id, "n_species",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, weakhub_table%n_species, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     if (weakhub_table%n_species == 6 .and. .not. use_muonic_table) then
      {#IFNDEF UNIT_TESTS
            call mpistop('If you dont use muon, why you use 6 species neutrino scheme')
      }
      {#IFDEF UNIT_TESTS
          write(99,*) 'If you dont use muon, why you use 6 species neutrino scheme'
      }
    end if 
     if (weakhub_table%n_species == 3 .and.  use_muonic_table) then
      {#IFNDEF UNIT_TESTS
            call mpistop('If you use muon, why you use 3 species neutrino scheme')
      }
      {#IFDEF UNIT_TESTS
         write(99,*)'If you use muon, why you use 3 species neutrino scheme'
      }
     end if

     call h5dopen_f(file_id, "n_groups",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, weakhub_table%n_groups, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     call h5dopen_f(file_id, "eps_min",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, weakhub_table%eps_min, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     call h5dopen_f(file_id, "eps_max",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, weakhub_table%eps_max, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     dims1_char(1) = 512
     call h5dopen_f(file_id, "energy_bin_type",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_CHARACTER, energy_bin_type, dims1_char, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     if (energy_bin_type == 'centroid') then 
       weakhub_table%energy_bin_centroid = .true.
     else
       weakhub_table%energy_bin_centroid = .false.
     endif

     call h5dopen_f(file_id, "energy_bin_spacing",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_CHARACTER, weakhub_table%energy_bin_spacing, dims1_char, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     allocate(weakhub_table%energies(0:weakhub_table%n_groups+1))
     rank = 1
     dims1(1) = n_ebin+2
     call h5dopen_f(file_id, "energies_table",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, weakhub_table%energies, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     allocate(weakhub_table%energies_i(0:weakhub_table%n_groups+1))
     rank = 1
     dims1(1) = n_ebin+2
     call h5dopen_f(file_id, "energies_i",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, weakhub_table%energies_i, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     allocate(weakhub_table%dV_eps(1:weakhub_table%n_groups))
     rank = 1
     dims1(1) = n_ebin+2
     call h5dopen_f(file_id, "energy_volume",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, weakhub_table%dV_eps, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     ! Read 4D table: Opacity
     if (use_IVtable) then
       rank = 1
       dims1(1) = 1
       call h5dopen_f(file_id, "IVrho",dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_INTEGER, IVrho, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error
   
       rank = 1
       dims1(1) = 1
       call h5dopen_f(file_id, "IVtemp",dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_INTEGER, IVtemp, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error
    
       dims1(1) = 1
       call h5dopen_f(file_id, "IVye",dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_INTEGER, IVye, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error
  
       dims1(1) = 1
       call h5dopen_f(file_id, "IVymu",dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_INTEGER, IVymu, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error
    
       if (mype == 0) write(*,*)  'IVrho, IVtemp, IVye, IVymu in the 4D table =',  &
                                   IVrho, IVtemp, IVye, IVymu
       allocate(logrho_IVtable(IVrho))
       allocate(logtemp_IVtable(IVtemp))
       allocate(ye_IVtable(IVye))
       allocate(logymu_IVtable(IVymu))
  
       rank = 1
       dims1(1) = IVrho
       call h5dopen_f(file_id, "logrho_IVtable", dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,logrho_IVtable, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error
  
       rank = 1
       dims1(1) = IVtemp
       call h5dopen_f(file_id, "logtemp_IVtable", dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,logtemp_IVtable, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error
    
       rank = 1
       dims1(1) = IVye
       call h5dopen_f(file_id, "ye_IVtable", dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,ye_IVtable, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error
  
       rank = 1
       dims1(1) = IVymu
       call h5dopen_f(file_id, "logymu_IVtable", dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,logymu_IVtable, dims1, error)
       call h5dclose_f(dset_id, error)
       cerror = cerror + error

       ! Upper and Lower values for the lookup table routines
       if (use_muonic_table) then
          x4_min(1) = logrho_IVtable(1);  x4_max(1) = logrho_IVtable(IVrho)
          x4_min(2) = logtemp_IVtable(1); x4_max(2) = logtemp_IVtable(IVtemp)
          x4_min(3) = ye_IVtable(1);      x4_max(3) = ye_IVtable(IVye)
          x4_min(4) = logymu_IVtable(1);  x4_max(4) = logymu_IVtable(IVymu)
       else
          x3_min(1) = logrho_IVtable(1);  x3_max(1) = logrho_IVtable(IVrho)
          x3_min(2) = logtemp_IVtable(1); x3_max(2) = logtemp_IVtable(IVtemp)
          x3_min(3) = ye_IVtable(1);      x3_max(3) = ye_IVtable(IVye)
       endif
  
     endif ! endif of use_IVtable for both spectral or grey tables
 
     if (weakhub_table_type == table_grey) then
       if (use_IVtable) then
         allocate(kappa_a_en_grey_table(IVrho,IVtemp,IVye,IVymu,weakhub_table%n_species))
         allocate(kappa_a_num_grey_table(IVrho,IVtemp,IVye,IVymu,weakhub_table%n_species))
         allocate(kappa_s_grey_table(IVrho,IVtemp,IVye,IVymu,weakhub_table%n_species))

         rank = 5
         dims5(1) = IVrho
         dims5(2) = IVtemp
         dims5(3) = IVye
         dims5(4) = IVymu
         dims5(5) = weakhub_table%n_species

         call h5dopen_f(file_id, "kappa_a_en_grey_table", dset_id, error)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kappa_a_en_grey_table, dims5, error)
         call h5dclose_f(dset_id, error)
         cerror = cerror + error

         call h5dopen_f(file_id, "kappa_a_num_grey_table", dset_id, error)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kappa_a_num_grey_table, dims5, error)
         call h5dclose_f(dset_id, error)
         cerror = cerror + error

         call h5dopen_f(file_id, "kappa_s_grey_table", dset_id, error)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kappa_s_grey_table, dims5, error)
         call h5dclose_f(dset_id, error)
         cerror = cerror + error

         !!if (use_muonic_table) then
         !!  weakhub_table%FourD_kappa_a_en_grey_table  = LT4_create_from_data(x4_min,x4_max,kappa_a_en_grey_table)
         !!  weakhub_table%FourD_kappa_a_num_grey_table  = LT4_create_from_data(x4_min,x4_max,kappa_a_num_grey_table)
         !!  weakhub_table%FourD_kappa_s_grey_table  = LT4_create_from_data(x4_min,x4_max,kappa_s_grey_table)
         !!else
         !!  weakhub_table%ThreeD_kappa_a_en_grey_table = &
         !!             LT3_create_from_data(x3_min,x3_max,&
         !!                             kappa_a_en_grey_table(1:IVrho,1:IVtemp,1:IVye,1,1:weakhub_table%n_species))
         !!  weakhub_table%ThreeD_kappa_a_num_grey_table = &
         !!             LT3_create_from_data(x3_min,x3_max,&
         !!                             kappa_a_num_grey_table(1:IVrho,1:IVtemp,1:IVye,1,1:weakhub_table%n_species))
         !!  weakhub_table%ThreeD_kappa_s_grey_table = &
         !!             LT3_create_from_data(x3_min,x3_max,&
         !!                             kappa_s_grey_table(1:IVrho,1:IVtemp,1:IVye,1,1:weakhub_table%n_species))
         !!endif

         !if (allocated(kappa_a_en_grey_table))      deallocate(kappa_a_en_grey_table)
         !if (allocated(kappa_a_num_grey_table))     deallocate(kappa_a_num_grey_table)
         !if (allocated(kappa_s_grey_table))         deallocate(kappa_s_grey_table)
       endif ! endif of use_IVtable
    endif !endif of grey table or not

    !must close h5 files, check for error
    if (cerror.ne.0) then
       if (mype == 0) write (*,*) "We have errors on reading HDF5 file", cerror
       {#IFNDEF UNIT_TESTS
       call mpistop('problem in weakhub reader')
       }
       {#IFDEF UNIT_TESTS
        write(99,*) 'problem in weakhub reader'
       }
    end if

    call h5fclose_f(file_id,error)
    call h5close_f(error)

    if (use_IVtable) then
      logrho_max_IV  = logrho_IVtable(IVrho) 
      logrho_min_IV  = logrho_IVtable(1)
      logtemp_max_IV = logtemp_IVtable(IVtemp) 
      logtemp_min_IV = logtemp_IVtable(1)
      logymu_max_IV  = logymu_IVtable(IVymu) 
      logymu_min_IV  = logymu_IVtable(1)
      ye_max_IV      = ye_IVtable(IVye) 
      ye_min_IV      = ye_IVtable(1)
      !if (allocated(logtemp_IVtable))            deallocate(logtemp_IVtable)
      !if (allocated(logrho_IVtable))             deallocate(logrho_IVtable)
      !if (allocated(ye_IVtable))                 deallocate(ye_IVtable)
      !if (allocated(logymu_IVtable))             deallocate(logymu_IVtable)
      if (mype == 0 ) then
        write(*,*) 'logrho_max_IV   ', logrho_max_IV
        write(*,*) 'logrho_min_IV   ', logrho_min_IV
        write(*,*) 'logtemp_max_IV  ', logtemp_max_IV
        write(*,*) 'logtemp_min_IV  ', logtemp_min_IV
        write(*,*) 'ye_max_IV       ', ye_max_IV
        write(*,*) 'ye_min_IV       ', ye_min_IV
        write(*,*) 'logymu_max_IV   ', logymu_max_IV
        write(*,*) 'logymu_min_IV   ', logymu_min_IV
      endif
    endif

    if (mype == 0) then
      write(*,*) '-----------Reading Weakhub tables successfully--------------!'
    endif

  end subroutine read_h5_Weakhub_table_gray

  !###########################################################################
  !###########################################################################
  !###########################################################################

end module mod_Weakhub_reader

