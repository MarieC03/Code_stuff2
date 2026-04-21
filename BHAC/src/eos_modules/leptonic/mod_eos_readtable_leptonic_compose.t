!================================================================================
!
!  mod_eos_readtable_leptonic_compose.t
!
!  Reads a COMPOSE-format baryon EOS table for use in the 4D leptonic EOS.
!  The third axis is interpreted as the total proton fraction yp = ye + ymu,
!  consistent with the Margherita/BHAC leptonic EOS split.
!
!  The routine maps the COMPOSE thermodynamic and composition fields onto
!    lep_tables_baryon(nrho, ntemp, nye, lep_nvars_baryon)
!  and converts them to the same conventions used by the active leptonic EOS:
!    - ln(rho) in code units
!    - ln(T / MeV)
!    - ln(P) and ln(eps + energy_shift)
!    - relativistic cs2
!
!================================================================================

module mod_eos_readtable_leptonic_compose
  use mod_eos
  use mod_eos_leptonic_parameters
  implicit none
  public

contains

  subroutine activate_baryon_tablereader_compose
    call readtable_baryon_compose(baryon_table_name)
  end subroutine activate_baryon_tablereader_compose

  subroutine readtable_baryon_compose(eos_filename)
    use hdf5
    use mod_eos
    use mod_eos_leptonic_parameters
    {#IFDEF UNIT_TESTS
    use amrvacdef
    }
    {#IFNDEF UNIT_TESTS
    include 'amrvacdef.f'
    }

    character(*), intent(in) :: eos_filename

    integer(HID_T)   :: file_id, dset_id, parameters_id, thermo_id, comp_id, av_id
    integer(HSIZE_T) :: dims1(1), dims3(3), dims4(4)
    integer          :: error, accerr, status_e
    integer          :: i, j, k, iv

    integer          :: nthermo_loc, ncomp_loc, nav_loc
    integer, allocatable :: thermo_index(:), index_yi(:)

    double precision, allocatable :: nb_axis(:), temp_axis(:), yp_axis(:)
    double precision, allocatable :: thermo_table(:,:,:,:)
    double precision, allocatable :: yi_table(:,:,:,:)
    double precision, allocatable :: zav_table(:,:,:), yav_table(:,:,:), aav_table(:,:,:)
    double precision, allocatable :: enthalpy(:,:,:)
    double precision :: eps_min_local, h_local, compose_baryon_mass

    accerr = 0
    if (mype == 0) write(*,*) "Reading baryon (compose) EOS table: ", trim(adjustl(eos_filename))

    call h5open_f(error)
    call h5fopen_f(trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)
    if (error /= 0) call mpistop("readtable_baryon_compose: cannot open file")

    call h5gopen_f(file_id, "/Parameters", parameters_id, error)
    if (error /= 0) call mpistop("readtable_baryon_compose: cannot open /Parameters")

    dims1(1) = 1
    call h5aopen_f(parameters_id, "pointsnb", dset_id, error)
    call h5aread_f(dset_id, H5T_NATIVE_INTEGER, lep_nrho, dims1, error)
    call h5aclose_f(dset_id, error)

    call h5aopen_f(parameters_id, "pointst", dset_id, error)
    call h5aread_f(dset_id, H5T_NATIVE_INTEGER, lep_ntemp, dims1, error)
    call h5aclose_f(dset_id, error)

    call h5aopen_f(parameters_id, "pointsyq", dset_id, error)
    call h5aread_f(dset_id, H5T_NATIVE_INTEGER, lep_nye, dims1, error)
    call h5aclose_f(dset_id, error)

    if (mype == 0) write(*,'(a,3i6)') "  nrho ntemp nye: ", lep_nrho, lep_ntemp, lep_nye

    if (allocated(lep_logrho_table)) deallocate(lep_logrho_table)
    if (allocated(lep_logtemp_table)) deallocate(lep_logtemp_table)
    if (allocated(lep_ye_table)) deallocate(lep_ye_table)
    allocate(lep_logrho_table(lep_nrho))
    allocate(lep_logtemp_table(lep_ntemp))
    allocate(lep_ye_table(lep_nye))
    allocate(nb_axis(lep_nrho))
    allocate(temp_axis(lep_ntemp))
    allocate(yp_axis(lep_nye))

    dims1(1) = lep_nrho
    call h5dopen_f(parameters_id, "nb", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nb_axis, dims1, error)
    call h5dclose_f(dset_id, error)
    accerr = accerr + error

    dims1(1) = lep_ntemp
    call h5dopen_f(parameters_id, "t", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_axis, dims1, error)
    call h5dclose_f(dset_id, error)
    accerr = accerr + error

    dims1(1) = lep_nye
    call h5dopen_f(parameters_id, "yq", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, yp_axis, dims1, error)
    call h5dclose_f(dset_id, error)
    accerr = accerr + error

    call h5gopen_f(file_id, "/Thermo_qty", thermo_id, error)
    if (error /= 0) call mpistop("readtable_baryon_compose: cannot open /Thermo_qty")

    dims1(1) = 1
    call h5aopen_f(thermo_id, "pointsqty", dset_id, error)
    call h5aread_f(dset_id, H5T_NATIVE_INTEGER, nthermo_loc, dims1, error)
    call h5aclose_f(dset_id, error)
    accerr = accerr + error

    allocate(thermo_index(nthermo_loc))
    dims1(1) = nthermo_loc
    call h5dopen_f(thermo_id, "index_thermo", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, thermo_index, dims1, error)
    call h5dclose_f(dset_id, error)
    accerr = accerr + error

    allocate(thermo_table(lep_nrho, lep_ntemp, lep_nye, nthermo_loc))
    dims4(1) = lep_nrho
    dims4(2) = lep_ntemp
    dims4(3) = lep_nye
    dims4(4) = nthermo_loc
    call h5dopen_f(thermo_id, "thermo", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, thermo_table, dims4, error)
    call h5dclose_f(dset_id, error)
    accerr = accerr + error

    call h5eset_auto_f(0, status_e)

    ncomp_loc = 0
    call h5gopen_f(file_id, "/Composition_pairs", comp_id, error)
    if (error == 0) then
      dims1(1) = 1
      call h5aopen_f(comp_id, "pointspairs", dset_id, error)
      call h5aread_f(dset_id, H5T_NATIVE_INTEGER, ncomp_loc, dims1, error)
      call h5aclose_f(dset_id, error)
      accerr = accerr + error

      if (ncomp_loc > 0) then
        allocate(index_yi(ncomp_loc))
        dims1(1) = ncomp_loc
        call h5dopen_f(comp_id, "index_yi", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, index_yi, dims1, error)
        call h5dclose_f(dset_id, error)
        accerr = accerr + error

        allocate(yi_table(lep_nrho, lep_ntemp, lep_nye, ncomp_loc))
        dims4(1) = lep_nrho
        dims4(2) = lep_ntemp
        dims4(3) = lep_nye
        dims4(4) = ncomp_loc
        call h5dopen_f(comp_id, "yi", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, yi_table, dims4, error)
        call h5dclose_f(dset_id, error)
        accerr = accerr + error
      end if
      call h5gclose_f(comp_id, error)
    end if

    nav_loc = 0
    call h5gopen_f(file_id, "/Composition_quadrupels", av_id, error)
    if (error /= 0) then
      call h5gopen_f(file_id, "/Composition_quadruples", av_id, error)
    end if
    if (error == 0) then
      dims1(1) = 1
      call h5aopen_f(av_id, "pointsav", dset_id, error)
      call h5aread_f(dset_id, H5T_NATIVE_INTEGER, nav_loc, dims1, error)
      call h5aclose_f(dset_id, error)
      accerr = accerr + error

      if (nav_loc > 0) then
        allocate(zav_table(lep_nrho, lep_ntemp, lep_nye))
        allocate(yav_table(lep_nrho, lep_ntemp, lep_nye))
        allocate(aav_table(lep_nrho, lep_ntemp, lep_nye))

        dims3(1) = lep_nrho
        dims3(2) = lep_ntemp
        dims3(3) = lep_nye

        call h5dopen_f(av_id, "zav", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, zav_table, dims3, error)
        call h5dclose_f(dset_id, error)
        accerr = accerr + error

        call h5dopen_f(av_id, "yav", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, yav_table, dims3, error)
        call h5dclose_f(dset_id, error)
        accerr = accerr + error

        call h5dopen_f(av_id, "aav", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, aav_table, dims3, error)
        call h5dclose_f(dset_id, error)
        accerr = accerr + error
      end if
      call h5gclose_f(av_id, error)
    end if

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    if (accerr /= 0) call mpistop("readtable_baryon_compose: HDF5 read error")

    if (allocated(lep_tables_baryon)) deallocate(lep_tables_baryon)
    allocate(lep_tables_baryon(lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon))
    lep_tables_baryon = 0.0d0

    do i = 1, lep_nrho
    do j = 1, lep_ntemp
    do k = 1, lep_nye
      do iv = 1, nthermo_loc
        if (thermo_index(iv) == 1) then
          lep_tables_baryon(i,j,k,i_lep_logpress) = thermo_table(i,j,k,iv)
        end if
        if (thermo_index(iv) == 2) then
          lep_tables_baryon(i,j,k,i_lep_entropy) = thermo_table(i,j,k,iv)
        end if
        if (thermo_index(iv) == 3) then
          lep_tables_baryon(i,j,k,i_lep_mu_n) = thermo_table(i,j,k,iv)
        end if
        if (thermo_index(iv) == 4) then
          lep_tables_baryon(i,j,k,i_lep_mu_p) = thermo_table(i,j,k,iv)
        end if
        if (thermo_index(iv) == 5) then
          lep_tables_baryon(i,j,k,i_lep_mu_e) = thermo_table(i,j,k,iv)
        end if
        if (thermo_index(iv) == 7) then
          lep_tables_baryon(i,j,k,i_lep_logenergy) = thermo_table(i,j,k,iv)
        end if
        if (thermo_index(iv) == 12) then
          lep_tables_baryon(i,j,k,i_lep_cs2) = thermo_table(i,j,k,iv)
        end if
      end do

      if (nav_loc > 0) then
        lep_tables_baryon(i,j,k,i_lep_abar) = aav_table(i,j,k)
        lep_tables_baryon(i,j,k,i_lep_zbar) = zav_table(i,j,k)
        lep_tables_baryon(i,j,k,i_lep_xh)   = aav_table(i,j,k) * yav_table(i,j,k)
      end if

      if (ncomp_loc > 0) then
        do iv = 1, ncomp_loc
          if (index_yi(iv) == 10) then
            lep_tables_baryon(i,j,k,i_lep_xn) = yi_table(i,j,k,iv)
          end if
          if (index_yi(iv) == 11) then
            lep_tables_baryon(i,j,k,i_lep_xp) = yi_table(i,j,k,iv)
          end if
          if (index_yi(iv) == 4002) then
            lep_tables_baryon(i,j,k,i_lep_xa) = 4.0d0 * yi_table(i,j,k,iv)
          end if
        end do
      end if
    end do
    end do
    end do

    deallocate(thermo_index, thermo_table)
    if (allocated(yi_table)) deallocate(yi_table)
    if (allocated(index_yi)) deallocate(index_yi)
    if (allocated(zav_table)) deallocate(zav_table)
    if (allocated(yav_table)) deallocate(yav_table)
    if (allocated(aav_table)) deallocate(aav_table)

    compose_baryon_mass = 1.674927211d-24
    lep_baryon_mass = compose_baryon_mass
    lep_have_rel_cs2 = 1

    if (mype == 0) then
      write(*,'(a,es14.6)') "  compose baryon mass used [g]: ", compose_baryon_mass
      if (abs(massn_cgs - compose_baryon_mass) > 1.0d-12*compose_baryon_mass) then
        write(*,'(a,es14.6,a,es14.6)') "  overriding massn_cgs = ", massn_cgs, &
             " with FIL-compatible neutron mass = ", compose_baryon_mass
      end if
    end if

    do i = 1, lep_nrho
      lep_logrho_table(i) = log(nb_axis(i) * compose_baryon_mass * cm3_to_fm3 * rho_gf)
    end do
    do j = 1, lep_ntemp
      lep_logtemp_table(j) = log(temp_axis(j))
    end do
    lep_ye_table(1:lep_nye) = yp_axis(1:lep_nye)

    lep_tables_baryon(:,:,:,i_lep_logpress) = log(lep_tables_baryon(:,:,:,i_lep_logpress) * &
         Mev_to_erg * cm3_to_fm3 * press_gf)

    lep_energy_shift = 0.0d0
    eps_min_local = minval(lep_tables_baryon(:,:,:,i_lep_logenergy))
    if (eps_min_local < 0.0d0) then
      lep_energy_shift = -2.0d0 * eps_min_local
      lep_tables_baryon(:,:,:,i_lep_logenergy) = lep_tables_baryon(:,:,:,i_lep_logenergy) + lep_energy_shift
    end if
    if (minval(lep_tables_baryon(:,:,:,i_lep_logenergy)) <= 0.0d0) then
      call mpistop("readtable_baryon_compose: non-positive eps after energy shift")
    end if
    lep_tables_baryon(:,:,:,i_lep_logenergy) = log(lep_tables_baryon(:,:,:,i_lep_logenergy))

    lep_tables_baryon(:,:,:,i_lep_cs2) = max(min(0.9999999d0, lep_tables_baryon(:,:,:,i_lep_cs2)), 0.0d0)

    lep_tables_baryon(:,:,:,i_lep_mu_e) = lep_tables_baryon(:,:,:,i_lep_mu_e) - lep_tables_baryon(:,:,:,i_lep_mu_p)
    lep_tables_baryon(:,:,:,i_lep_mu_p) = lep_tables_baryon(:,:,:,i_lep_mu_p) + lep_tables_baryon(:,:,:,i_lep_mu_n)

    allocate(enthalpy(lep_nrho, lep_ntemp, lep_nye))
    do k = 1, lep_nye
    do j = 1, lep_ntemp
    do i = 1, lep_nrho
      h_local = 1.0d0 + exp(lep_tables_baryon(i,j,k,i_lep_logenergy)) - lep_energy_shift + &
           exp(lep_tables_baryon(i,j,k,i_lep_logpress)) / exp(lep_logrho_table(i))
      enthalpy(i,j,k) = h_local
    end do
    end do
    end do

    lep_eos_rhomin  = exp(lep_logrho_table(1))
    lep_eos_rhomax  = exp(lep_logrho_table(lep_nrho))
    lep_eos_tempmin = exp(lep_logtemp_table(1))
    lep_eos_tempmax = exp(lep_logtemp_table(lep_ntemp))
    lep_eos_yemin   = lep_ye_table(1)
    lep_eos_yemax   = lep_ye_table(lep_nye)
    lep_eos_epsmin  = exp(minval(lep_tables_baryon(:,:,:,i_lep_logenergy))) - lep_energy_shift
    lep_eos_epsmax  = exp(maxval(lep_tables_baryon(:,:,:,i_lep_logenergy))) - lep_energy_shift
    lep_eos_hmin    = minval(enthalpy)

    deallocate(enthalpy, nb_axis, temp_axis, yp_axis)

    if (mype == 0) then
      write(*,'(a,2es14.6)') "  rho  range  [code]: ", lep_eos_rhomin,  lep_eos_rhomax
      write(*,'(a,2es14.6)') "  temp range  [MeV]:  ", lep_eos_tempmin, lep_eos_tempmax
      write(*,'(a,2es14.6)') "  yp   range:         ", lep_eos_yemin,   lep_eos_yemax
      write(*,'(a,es14.6)')  "  h_min:              ", lep_eos_hmin
      write(*,*) "  Baryon (compose) table loaded successfully."
    end if

  end subroutine readtable_baryon_compose

end module mod_eos_readtable_leptonic_compose
