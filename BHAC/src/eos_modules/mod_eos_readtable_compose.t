module mod_eos_readtable_compose
  use mod_eos
  use mod_eos_tabulated_parameters
  implicit none
  public

contains

   subroutine activate_tablereader_compose
      table_type = compose
      call readtable_compose(eos_table_name)
   end subroutine activate_tablereader_compose

   subroutine readtable_compose(eos_filename)
   ! This routine reads the table and initializes
   ! all variables in the module. 
     use hdf5 
     {#IFDEF UNIT_TESTS
     use amrvacdef
     }
     {#IFNDEF UNIT_TESTS
     include 'amrvacdef.f'
     }

   
     character(*), intent(in)      :: eos_filename
   
     character(len=100) :: message
   
   ! HDF5 vars
     integer(HID_T) :: file_id, dset_id, dspace_id, parameters_id, thermo_id, comp_id, av_id
     integer(HSIZE_T) :: dims1(1), dims3(3), dims4(4)
     logical :: dataset_there
     integer :: error,rank,accerr
     integer :: i,j,k,iv
     integer :: status_e, status_comp, nullptr
 
     double precision :: amu_cgs_andi, eps_min_local, eos_prsmin, dummy
     double precision :: buffer1,buffer2,buffer3,buffer4
     double precision, allocatable      :: enthalpy(:,:,:)
     double precision, allocatable      :: thermo_table(:,:,:,:)
     double precision, allocatable      :: yi_table(:,:,:,:)
     double precision, allocatable      :: zav_table(:,:,:)
     double precision, allocatable      :: yav_table(:,:,:)
     double precision, allocatable      :: aav_table(:,:,:)
     integer, allocatable               :: thermo_index(:)
     integer, allocatable               :: index_yi(:)
    
     accerr=0
   
     if (mype==0) write(*,*) "Reading Compose EOS Table"
   
     call h5open_f(error)
   
     call h5fopen_f (trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)
   
     if (mype==0) write(6,*) trim(adjustl(eos_filename))

     ! use H5Gopen_f to get sub group of the file_id
     call h5gopen_f (file_id, "/Parameters", parameters_id, error)
     ! use H5Aopen_f instead of H5Dopen_f
     dims1(1)=1
     call h5Aopen_f(parameters_id, "pointsnb", dset_id, error)
     call h5Aread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
     call h5Aclose_f(dset_id,error)
  
     if(error.ne.0) then
        stop "Could not read EOS table file"
     endif
   
     dims1(1)=1
     call h5Aopen_f(parameters_id, "pointst", dset_id, error)
     call h5Aread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
     call h5Aclose_f(dset_id,error)
   
     if(error.ne.0) then
        stop "Could not read EOS table file"
     endif
   
     dims1(1)=1
     call h5Aopen_f(parameters_id, "pointsyq", dset_id, error)
     call h5Aread_f(dset_id, H5T_NATIVE_INTEGER, nye, dims1, error)
     call h5Aclose_f(dset_id,error)
   
     if(error.ne.0) then
        stop "Could not read EOS table file"
     endif

     if (mype==0) write(message,"(a25,i5,i5,i5)") "We have nrho ntemp nye: ", nrho,ntemp,nye
     if (mype==0) write(*,*) message
   
     allocate(logrho_table(nrho))
     dims1(1)=nrho
     call h5dopen_f(parameters_id, "nb", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho_table, dims1, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error

     allocate(logtemp_table(ntemp))
     dims1(1)=ntemp
     call h5dopen_f(parameters_id, "t", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp_table, dims1, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error
   
     allocate(ye_table(nye))
     dims1(1)=nye
     call h5dopen_f(parameters_id, "yq", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye_table, dims1, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error

     ! thermo table

     ! use H5Gopen_f to get sub group of the file_id
     call h5gopen_f (file_id, "/Thermo_qty", thermo_id, error)
     ! use H5Aopen_f instead of H5Dopen_f
     dims1(1)=1
     call h5Aopen_f(thermo_id, "pointsqty", dset_id, error)
     call h5Aread_f(dset_id, H5T_NATIVE_INTEGER, nthermo, dims1, error)
     call h5Aclose_f(dset_id,error)

     if (mype==0) write(message,"(a25,i5,i5,i5)") "We have nthermo", nthermo
     if (mype==0) write(*,*) message

     allocate(thermo_index(nthermo))
     dims1(1)=nthermo
     call h5dopen_f(thermo_id, "index_thermo", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, thermo_index, dims1, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error

     allocate(thermo_table(nrho, ntemp, nye, nthermo))
     dims4(1)=nrho
     dims4(2)=ntemp
     dims4(3)=nye
     dims4(4)=nthermo
     call h5dopen_f(thermo_id, "thermo", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, thermo_table(:,:,:,:), dims4, error)
     call h5dclose_f(dset_id,error)
     accerr=accerr+error

     ! Compositions
     call h5eset_auto_f(0, status_e)

     ! use H5Gopen_f to get sub group of the file_id
     call h5gopen_f (file_id, "/Composition_pairs", comp_id, error)
     ! use H5Aopen_f instead of H5Dopen_f
     dims1(1)=1
     call h5Aopen_f(comp_id, "pointspairs", dset_id, error)
     call h5Aread_f(dset_id, H5T_NATIVE_INTEGER, ncomp, dims1, error)
     call h5Aclose_f(dset_id,error)

     if (mype==0) write(message,"(a25,i5,i5,i5)") "We have ncomp: ", ncomp
     if (mype==0) write(*,*) message

     if (ncomp > 0 ) then
        allocate(index_yi(ncomp)) 
        dims1(1)=ncomp
        call h5dopen_f(comp_id, "index_yi", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, index_yi, dims1, error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        
        allocate(yi_table(nrho,ntemp,nye,ncomp)) 
        dims4(1)=nrho
        dims4(2)=ntemp
        dims4(3)=nye
        dims4(4)=ncomp
        call h5dopen_f(comp_id, "yi", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, yi_table(:,:,:,:), dims4, error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
     endif  

     nav = 0
     ! we always assume nav = 1 as this is NOT standardized!

     ! use H5Gopen_f to get sub group of the file_id
     call h5gopen_f (file_id, "/Composition_quadrupels", av_id, error)
     !call h5gopen_f (file_id, "/Composition_quadruples", av_id, error)
     ! use H5Aopen_f instead of H5Dopen_f
     dims1(1)=nav
     call h5Aopen_f(av_id, "pointsav", dset_id, error)
     call h5Aread_f(dset_id, H5T_NATIVE_INTEGER, nav, dims1, error)
     call h5Aclose_f(dset_id,error)
     if (mype==0) write(message,"(a25,i5,i5,i5)") "We have nav: ", nav
     if (mype==0) write(*,*) message

     if (nav > 0 ) then
        if (nav .ne. 1) then
          if (mype==0) write(message,*) "nav != 1 in this table &
                               so there is none or more than &
                               one definition of an average nucleus."
          if (mype==0) write(*,*) message
     
          if (mype==0) write(message,*) "  Please check and &
                                generalize accordingly."
     
          if (mype==0) write(*,*) message
        endif

        allocate(zav_table(nrho,ntemp,nye))
        dims3(1)=nrho
        dims3(2)=ntemp
        dims3(3)=nye
        call h5dopen_f(av_id, "zav", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, zav_table(:,:,:), dims3, error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        allocate(yav_table(nrho,ntemp,nye))
        dims3(1)=nrho
        dims3(2)=ntemp
        dims3(3)=nye
        call h5dopen_f(av_id, "yav", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, yav_table(:,:,:), dims3, error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        allocate(aav_table(nrho,ntemp,nye))
        dims3(1)=nrho
        dims3(2)=ntemp
        dims3(3)=nye
        call h5dopen_f(av_id, "aav", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, aav_table(:,:,:), dims3, error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
     else if (nav == 0) then
          if (mype==0) write(message,*) "nav = 0 in this table &
                               No nuclei composition provided "
          if (mype==0) write(*,*) message
 
     endif

     if(accerr.ne.0) then
       stop "Problem reading Compose EOS table file"
     endif

     call h5fclose_f (file_id,error)
     call h5close_f (error)

     ! allocate eos_tables and fill the values
     allocate(eos_tables(nrho,ntemp,nye,nvars))

     do i = 1, nrho
     do j = 1, ntemp
     do k = 1, nye
           ! for the thermo table
           ! thermo quantities and its index
           ! press: 1;    s: 2;    mu_b-m_n: 3;   mu_q: 4;   mu_l: 5
           ! (e/m_n - 1): 7;    dp/dn_b |E : 10;    dp/dE|n_b: 11;
           ! cs2: 12;      
     
        do iv = 1, nthermo 
           if (thermo_index(iv) == 1) then      
             eos_tables(i,j,k,i_logpress)  = thermo_table(i,j,k,iv)
           endif
           if (thermo_index(iv) == 2) then      
             if (use_realistic_mp_table) eos_tables(i,j,k,i_entropy)   = thermo_table(i,j,k,iv) 
           endif
           if (thermo_index(iv) == 3) then      
             eos_tables(i,j,k,i_mu_n)      = thermo_table(i,j,k,iv) 
           endif
           if (thermo_index(iv) == 4) then      
             eos_tables(i,j,k,i_mu_p)      = thermo_table(i,j,k,iv) 
           endif
           if (thermo_index(iv) == 5) then      
             eos_tables(i,j,k,i_mu_e)      = thermo_table(i,j,k,iv) 
           endif
           if (thermo_index(iv) == 7) then      
             eos_tables(i,j,k,i_logenergy) = thermo_table(i,j,k,iv) 
           endif
           if (thermo_index(iv) == 10) then      
             if (use_realistic_mp_table) eos_tables(i,j,k,i_dpdrhoe)   = thermo_table(i,j,k,iv) 
           endif
           if (thermo_index(iv) == 11) then      
             if (use_realistic_mp_table) eos_tables(i,j,k,i_dpderho)   = thermo_table(i,j,k,iv) 
           endif
           if (thermo_index(iv) == 12) then      
             eos_tables(i,j,k,i_cs2)       = thermo_table(i,j,k,iv) 
           endif
        enddo

        if (use_realistic_mp_table) then    
          if ( nav == 1 ) then
                eos_tables(i,j,k,i_abar) = aav_table(i,j,k)  
                eos_tables(i,j,k,i_zbar) = zav_table(i,j,k)  
                eos_tables(i,j,k,i_xh)   = aav_table(i,j,k) * yav_table(i,j,k)
          endif
              
          if ( ncomp > 0 ) then
                ! for the yi_table
                ! particle fraction and its corresponding index
                ! Xn: 10;   Xp: 11;  Xa: 4002;  Xe: 0;  X_quark_u: 500;  X_quark_d: 501
                ! X_quark_s: 502;
             do iv = 1, ncomp
                if (index_yi(iv) == 10) then
                  eos_tables(i,j,k,i_xn) = yi_table(i,j,k,iv) 
                endif
                if (index_yi(iv) == 11) then
                  eos_tables(i,j,k,i_xp) = yi_table(i,j,k,iv) 
                endif
                if (index_yi(iv) == 4002) then
                  eos_tables(i,j,k,i_xa) = 4.0d0 * yi_table(i,j,k,iv) 
                endif
                if (index_yi(iv) == 500) then
                  eos_tables(i,j,k,i_quark_u) = yi_table(i,j,k,iv) 
                endif
                if (index_yi(iv) == 501) then
                  eos_tables(i,j,k,i_quark_d) = yi_table(i,j,k,iv) 
                endif
                if (index_yi(iv) == 502) then
                  eos_tables(i,j,k,i_quark_s) = yi_table(i,j,k,iv) 
                endif
             enddo
          endif !endif of ncomp > 0
        endif !endif of use_realistic_mp_table
     enddo
     enddo
     enddo

     deallocate(thermo_index)
     deallocate(thermo_table)

     if (nav > 0 ) then
       deallocate(zav_table)
       deallocate(yav_table)
       deallocate(aav_table)
     endif
     if (ncomp > 0) then
       deallocate(yi_table)
       deallocate(index_yi)
     endif

     !------------------------------------------------------------------------------!
     !------------------------------------------------------------------------------!
     !--------------------------  Conversion of units   ----------------------------!
     !------------------------------------------------------------------------------!
     !------------------------------------------------------------------------------!

     ! assume baryon_mass is neutron mass (not optimal, but we only have this in compose) : massn_cgs
     ! convert n_b to rho_b
     if (mype==0) write(*, *) "The baryon mass is set to be massn_cgs = ", massn_cgs, "g, when reading compose EoS table."
     logrho_table(1:nrho)   = log10(logrho_table(1:nrho) * massn_cgs * cm3_to_fm3 * rho_gf)

     logtemp_table(1:ntemp) =  log10(logtemp_table(1:ntemp)) 

     ! Pressure
     eos_tables(1:nrho,1:ntemp,1:nye,i_logpress)  = log10(eos_tables(1:nrho,1:ntemp,1:nye,i_logpress) * &
                                                     Mev_to_erg * cm3_to_fm3 * press_gf)

     ! Specific internal energy (without rest mass)
     energy_shift = 0.0d0

     ! get the epsmin first
     eps_min_local = minval(eos_tables(:,:,:,i_logenergy))

     ! if eps in table has -ve values, we need to assign a energy shift for it
     if (eps_min_local < 0.0d0) then
        energy_shift = -2.0d0 * eps_min_local
        eos_tables(1:nrho,1:ntemp,1:nye,i_logenergy) = eos_tables(1:nrho,1:ntemp,1:nye,i_logenergy) + energy_shift
        do i = 1, nrho
        do j = 1, nye
        do k = 1, ntemp
           if (eos_tables(i,k,j,i_logenergy) .lt. 0.0d0) then
               stop 'eos_tables (logeps) still have -ve values after + energyshift'
           endif
        enddo
        enddo
        enddo
     endif

     ! It is aleady dimensionless, no need to change unit
     ! Fixme
     ! Question : page 62. of compose docu:  this is not specific internal energy and not store rest mass of neutron
     ! So, the internal energy definition should include rest mass in nuclear physics, but astro physics assumes 
     ! already no rest mass contribution inside internal energy? 
     eos_tables(1:nrho,1:ntemp,1:nye,i_logenergy) = log10(eos_tables(1:nrho,1:ntemp,1:nye,i_logenergy)) 

     ! Squared speed of sound, it is already dimensionless, but it needs to limit
     eos_tables(1:nrho,1:ntemp,1:nye,i_cs2) = max(min(0.999999d0,eos_tables(1:nrho,1:ntemp,1:nye,i_cs2)), 0.0d0)

     if (use_realistic_mp_table) then
       ! Derivatives
       ! these two are wrong, later will handle them
       eos_tables(1:nrho,1:ntemp,1:nye,i_dpdrhoe) = eos_tables(1:nrho,1:ntemp,1:nye,i_dpdrhoe)/&
                                                    (massn_cgs * cm3_to_fm3 * rho_gf) * (press_gf)
       eos_tables(1:nrho,1:ntemp,1:nye,i_dpderho) = eos_tables(1:nrho,1:ntemp,1:nye,i_dpderho)/(eps_gf) * press_gf
     endif

     ! Chemical potentials -------------------------------------------------------!
     !  all the chemical potentials in compose table are including rest mass contributions, see Eq(3.22)
     !  mu_e = mu_le - mu_q :  Eq(3.23),  defintion of mu_l : Eq(4.7)
     eos_tables(1:nrho,1:ntemp,1:nye,i_mu_e) = eos_tables(1:nrho,1:ntemp,1:nye,i_mu_e) - eos_tables(1:nrho,1:ntemp,1:nye,i_mu_p)
     !  mu_p = mu_b + mu_q = mu_n + mu_q
     eos_tables(1:nrho,1:ntemp,1:nye,i_mu_p) = eos_tables(1:nrho,1:ntemp,1:nye,i_mu_p) + eos_tables(1:nrho,1:ntemp,1:nye,i_mu_n)

     eos_tables(1:nrho,1:ntemp,1:nye,i_mu_n) = eos_tables(1:nrho,1:ntemp,1:nye,i_mu_n)

     if (use_realistic_mp_table) then
       ! Fix munu in stellarcollapse tables, now properly included rest mass
       eos_tables(1:nrho,1:ntemp,1:nye,i_munu) = eos_tables(1:nrho,1:ntemp,1:nye,i_mu_e) + &
                                                  eos_tables(1:nrho,1:ntemp,1:nye,i_mu_p) - &
                                                  eos_tables(1:nrho,1:ntemp,1:nye,i_mu_n) + mp_mev - mn_mev
     endif

     ! remarks: We make the explicit assumption that we have no muons....
     !         I know we have to fix this later, but for most EOS this is ok
     !         And if we had muons, I'm pretty sure the neutrino microphysics would become
     !         inconsistent... 
     ! *** it is because mu_l = (mu_le * n_le + mu_lnu * n_lnu)/ n_l,  if there is muon contributions,
     !      mu_e is not equal to mu_l 
     ! Page 8: Assumptions on the relation between
     !         the electron and muon chemical potentials are discussed
     !         in the description of each model separately.
     ! Page 10: In this case, the balance between the
     !          electron and muon densities depends on the assumed relation of
     !          the electron and muon chemical potentials.
     ! end chemical potentials ----------------------------------------------------!

     ! Enthalpy

     allocate(enthalpy(nrho,ntemp,nye))
     do k=1,nye
        do j=1,ntemp
           do i=1,nrho
              enthalpy(i,j,k) = (1.0d0 + 10.0d0**eos_tables(i,j,k,i_logenergy) - energy_shift + &
                   10.0d0**eos_tables(i,j,k,i_logpress)/10.0d0**logrho_table(i))
           enddo
        enddo
     enddo
   
     ! set min-max values:
     eos_rhomin = 10.0d0**logrho_table(1)
     eos_rhomax = 10.0d0**logrho_table(nrho)
   
     eos_yemin = ye_table(1)
     eos_yemax = ye_table(nye)
   
     eos_tempmin = 10.0d0**logtemp_table(1)
     eos_tempmax = 10.0d0**logtemp_table(ntemp)

     ! some of the tables, it has -1.0 eps for extremely high temperature (corrupted table)
     ! so we need to use another way to find eos_epsmin, max in activate_eos_tabulated 
      eos_epsmin = 10.0d0**minval(eos_tables(:,:,:,i_logenergy)) - energy_shift
      eos_epsmax = 10.0d0**maxval(eos_tables(:,:,:,i_logenergy)) - energy_shift

      ! corrupted compose table, need another way to find eos epsmin, hmin
      if (eos_epsmin .lt. -0.999) then 
        corrupted_compose = .true.
        {#IFNDEF UNIT_TESTS
        call mpistop('Corrupted compose')
        }
      endif

      if (.not. corrupted_compose) then
        eos_hmin = minval(enthalpy(:,:,:))
        eos_hmin = min(1.0d0, eos_hmin)
      endif

     deallocate(enthalpy)
     if (mype==0) write(6,*) "Done reading compose eos table"

   end subroutine readtable_compose
end module 
