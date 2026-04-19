!================================================================================
!
!  mod_eos_readtable_leptonic_scollapse.t
!
!  Reads the stellarcollapse.org baryon EOS table for the 4D leptonic EOS.
!  The table is a standard 3D stellarcollapse table parametrised by
!  (logrho, logtemp, yp) where yp = yle + ymu.
!
!  After reading the baryon sub-table is stored in
!    lep_tables_baryon(nrho, ntemp, nye, lep_nvars_baryon)
!  and global bounds lep_eos_rhomin/max, lep_eos_tempmin/max,
!  lep_eos_yemin/max  (= bounds on yp) are set.
!
!  The routine MUST be called BEFORE readtable_leptonic so that the
!  shared axis arrays are initialised first.
!
!  Copyright (C) 2024  Harry Ho-Yin Ng
!  Based on BHAC mod_eos_readtable_scollapse and the Margherita
!  leptonic_eos_readtable_scollapse.hh by Harry Ho-Yin Ng.
!
!================================================================================

module mod_eos_readtable_leptonic_scollapse
  use mod_eos_leptonic_parameters
  implicit none
  public

contains

  subroutine activate_baryon_tablereader_scollapse
    call readtable_baryon_scollapse(baryon_table_name)
  end subroutine activate_baryon_tablereader_scollapse

  !============================================================================
  subroutine readtable_baryon_scollapse(eos_filename)
    use hdf5
    use mod_eos_leptonic_parameters
    include 'amrvacdef.f'

    character(*), intent(in) :: eos_filename

    integer(HID_T)   :: file_id, dset_id
    integer(HSIZE_T) :: dims1(1), dims1d(1)
    logical          :: dataset_there
    integer          :: error, accerr, i, j, k, iv, nrho_loc, ntemp_loc, nye_loc

    ! Unit-conversion factors (same as in standard BHAC reader)
    double precision, parameter :: rho_gf   = 1.61887093132742d-18  ! g/cm^3 -> CU
    double precision, parameter :: press_gf = 1.80123683248503d-39  ! dyn/cm^2 -> CU
    double precision, parameter :: eps_gf   = 1.11265005605362d-21  ! -> CU
    double precision, parameter :: len_gf   = 6.77269222552442d-06  ! cm -> CU
    double precision, parameter :: time_gf  = 2.03040204956746d+05  ! s -> CU
    double precision, parameter :: c_cgs    = 2.99792458d+10        ! cm/s

    double precision, allocatable :: logrho_raw(:), logtemp_raw(:), ye_raw(:)
    double precision, allocatable :: flat_tmp(:)  ! HDF5 read buffer (NTABLES * nrho*ntemp*nye)
    double precision, allocatable :: enthalpy(:,:,:)

    integer :: indold, indnew
    integer, parameter :: NTABLES = 13  ! = lep_nvars_baryon
    double precision   :: h_local

    ! ---
    accerr = 0
    if (mype == 0) write(*,*) "Reading baryon (scollapse) EOS table: ", trim(adjustl(eos_filename))

    call h5open_f(error)
    call h5fopen_f(trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)
    if (error /= 0) call mpistop("readtable_baryon_scollapse: cannot open file")

    !--------------------------------------------------------------------------
    ! Table dimensions
    !--------------------------------------------------------------------------
    dims1(1) = 1
    call h5dopen_f(file_id, "pointsrho",  dset_id, error); call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho_loc,  dims1, error); call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "pointstemp", dset_id, error); call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp_loc, dims1, error); call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, "pointsye",   dset_id, error); call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nye_loc,   dims1, error); call h5dclose_f(dset_id, error)

    if (mype == 0) write(*,'(a,3i6)') "  nrho ntemp nye: ", nrho_loc, ntemp_loc, nye_loc

    lep_nrho  = nrho_loc
    lep_ntemp = ntemp_loc
    lep_nye   = nye_loc

    !--------------------------------------------------------------------------
    ! Check for relativistic cs2 flag
    !--------------------------------------------------------------------------
    call h5lexists_f(file_id, "have_rel_cs2", dataset_there, error)
    if (dataset_there) then
      call h5dopen_f(file_id, "have_rel_cs2", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, lep_have_rel_cs2, dims1, error)
      call h5dclose_f(dset_id, error)
    else
      lep_have_rel_cs2 = 0
    end if
    if (mype == 0 .and. lep_have_rel_cs2 /= 0) write(*,*) "  cs2 is already relativistic in table"

    !--------------------------------------------------------------------------
    ! Allocate axis arrays (also needed by leptonic reader)
    !--------------------------------------------------------------------------
    if (.not. allocated(lep_logrho_table))  allocate(lep_logrho_table(lep_nrho))
    if (.not. allocated(lep_logtemp_table)) allocate(lep_logtemp_table(lep_ntemp))
    allocate(lep_ye_table(lep_nye))
    allocate(logrho_raw(lep_nrho))
    allocate(logtemp_raw(lep_ntemp))
    allocate(ye_raw(lep_nye))

    dims1d(1) = lep_nrho
    call h5dopen_f(file_id, "logrho",  dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho_raw,  dims1d, error)
    call h5dclose_f(dset_id, error)

    dims1d(1) = lep_ntemp
    call h5dopen_f(file_id, "logtemp", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp_raw, dims1d, error)
    call h5dclose_f(dset_id, error)

    dims1d(1) = lep_nye
    call h5dopen_f(file_id, "ye",      dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye_raw,      dims1d, error)
    call h5dclose_f(dset_id, error)

    dims1(1) = 1
    call h5dopen_f(file_id, "energy_shift", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_energy_shift, dims1, error)
    call h5dclose_f(dset_id, error)

    !> Read baryon mass if present
    call h5lexists_f(file_id, "mass_factor", dataset_there, error)
    if (dataset_there) then
      call h5dopen_f(file_id, "mass_factor", dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_baryon_mass, dims1, error)
      call h5dclose_f(dset_id, error)
      if (mype == 0) write(*,'(a,es14.6)') "  baryon mass from table: ", lep_baryon_mass
    else
      if (mype == 0) write(*,*) "  Using default baryon mass"
    end if

    !--------------------------------------------------------------------------
    ! Read flat table data  layout: (iv, irho*itemp*iye)  (C order from HDF5)
    !--------------------------------------------------------------------------
    allocate(flat_tmp(NTABLES * lep_nrho * lep_ntemp * lep_nye))

    ! Helper macro expanded inline:
    ! READ_EOSTABLE_HDF5 reads a hyperslab at offset iv into the flat buffer
    ! Here we read each dataset separately into consecutive slabs.

    call read_slab(file_id, "logpress",  flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_logpress,  error); accerr=accerr+error
    call read_slab(file_id, "logenergy", flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_logenergy, error); accerr=accerr+error
    call read_slab(file_id, "entropy",   flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_entropy,   error); accerr=accerr+error
    call read_slab(file_id, "cs2",       flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_cs2,       error); accerr=accerr+error
    call read_slab(file_id, "mu_e",      flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_mu_e,      error); accerr=accerr+error
    call read_slab(file_id, "mu_p",      flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_mu_p,      error); accerr=accerr+error
    call read_slab(file_id, "mu_n",      flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_mu_n,      error); accerr=accerr+error
    call read_slab(file_id, "Xa",        flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_xa,        error); accerr=accerr+error
    call read_slab(file_id, "Xh",        flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_xh,        error); accerr=accerr+error
    call read_slab(file_id, "Xn",        flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_xn,        error); accerr=accerr+error
    call read_slab(file_id, "Xp",        flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_xp,        error); accerr=accerr+error
    call read_slab(file_id, "Abar",      flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_abar,      error); accerr=accerr+error
    call read_slab(file_id, "Zbar",      flat_tmp, lep_nrho, lep_ntemp, lep_nye, NTABLES, i_lep_zbar,      error); accerr=accerr+error

    if (accerr /= 0) call mpistop("readtable_baryon_scollapse: HDF5 read error")

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    !--------------------------------------------------------------------------
    ! Convert axes to code units (log -> natural log, apply GF)
    !--------------------------------------------------------------------------
    lep_energy_shift = lep_energy_shift * eps_gf   ! apply unit conversion

    do i = 1, lep_nrho
      logrho_raw(i) = logrho_raw(i) * log(10.0d0) + log(rho_gf)
    end do
    do j = 1, lep_ntemp
      logtemp_raw(j) = logtemp_raw(j) * log(10.0d0)  ! stays in ln(MeV)
    end do
    lep_logrho_table(1:lep_nrho)   = logrho_raw(1:lep_nrho)
    lep_logtemp_table(1:lep_ntemp) = logtemp_raw(1:lep_ntemp)
    lep_ye_table(1:lep_nye)        = ye_raw(1:lep_nye)

    !--------------------------------------------------------------------------
    ! Reindex from C (iv, irho*itemp*iye) to Fortran (irho, itemp, iye, iv)
    ! and apply unit conversions to table entries
    !--------------------------------------------------------------------------
    allocate(lep_tables_baryon(lep_nrho, lep_ntemp, lep_nye, lep_nvars_baryon))
    allocate(enthalpy(lep_nrho, lep_ntemp, lep_nye))

    do iv = 1, NTABLES
    do k  = 1, lep_nye
    do j  = 1, lep_ntemp
    do i  = 1, lep_nrho
      indold = (iv-1) * (lep_nrho*lep_ntemp*lep_nye) + (i-1) + lep_nrho*((j-1) + lep_ntemp*(k-1)) + 1
      lep_tables_baryon(i, j, k, iv) = flat_tmp(indold)
    end do; end do; end do; end do
    deallocate(flat_tmp)

    ! --- unit conversions ---
    do k = 1, lep_nye
    do j = 1, lep_ntemp
    do i = 1, lep_nrho
      ! pressure: log10 -> ln  +  CGS -> code
      lep_tables_baryon(i,j,k,i_lep_logpress) = &
          lep_tables_baryon(i,j,k,i_lep_logpress) * log(10.0d0) + log(press_gf)

      ! energy: log10 -> ln  +  CGS -> code
      lep_tables_baryon(i,j,k,i_lep_logenergy) = &
          lep_tables_baryon(i,j,k,i_lep_logenergy) * log(10.0d0) + log(eps_gf)

      ! cs2: CGS -> code   ( len_gf^2 / time_gf^2 )
      lep_tables_baryon(i,j,k,i_lep_cs2) = &
          lep_tables_baryon(i,j,k,i_lep_cs2) * (len_gf*len_gf) / (time_gf*time_gf)
      if (lep_tables_baryon(i,j,k,i_lep_cs2) < 0.0d0) &
          lep_tables_baryon(i,j,k,i_lep_cs2) = 0.0d0

      ! Compute specific enthalpy h = 1 + eps + p/rho  for cs2 correction
      ! and for determining h_min
      h_local = 1.0d0 &
          + (exp(lep_tables_baryon(i,j,k,i_lep_logenergy)) - lep_energy_shift) &
          + exp(lep_tables_baryon(i,j,k,i_lep_logpress)) / exp(logrho_raw(i))
      enthalpy(i,j,k) = h_local

      ! Make cs2 relativistic (divide by h) if not already done
      if (lep_have_rel_cs2 == 0) then
        lep_tables_baryon(i,j,k,i_lep_cs2) = lep_tables_baryon(i,j,k,i_lep_cs2) / h_local
      end if
      lep_tables_baryon(i,j,k,i_lep_cs2) = min(0.9999999d0, lep_tables_baryon(i,j,k,i_lep_cs2))
    end do; end do; end do

    !--------------------------------------------------------------------------
    ! Set global bounds
    !--------------------------------------------------------------------------
    lep_eos_rhomin  = exp(logrho_raw(1))
    lep_eos_rhomax  = exp(logrho_raw(lep_nrho))
    lep_eos_tempmin = exp(logtemp_raw(1))
    lep_eos_tempmax = exp(logtemp_raw(lep_ntemp))
    lep_eos_yemin   = ye_raw(1)
    lep_eos_yemax   = ye_raw(lep_nye)
    lep_eos_epsmin  = exp(minval(lep_tables_baryon(:,:,:,i_lep_logenergy))) - lep_energy_shift
    lep_eos_epsmax  = exp(maxval(lep_tables_baryon(:,:,:,i_lep_logenergy))) - lep_energy_shift
    lep_eos_hmin    = minval(enthalpy)

    deallocate(enthalpy, logrho_raw, logtemp_raw, ye_raw)

    if (mype == 0) then
      write(*,'(a,2es14.6)') "  rho  range  [code]: ", lep_eos_rhomin,  lep_eos_rhomax
      write(*,'(a,2es14.6)') "  temp range  [MeV]:  ", lep_eos_tempmin, lep_eos_tempmax
      write(*,'(a,2es14.6)') "  yp   range:         ", lep_eos_yemin,   lep_eos_yemax
      write(*,'(a,es14.6)')  "  h_min:              ", lep_eos_hmin
      write(*,*) "  Baryon (scollapse) table loaded successfully."
    end if

  end subroutine readtable_baryon_scollapse

  !============================================================================
  !> Helper: read a single 1D dataset slice from the flat HDF5 table array
  !> at offset (iv0-1)*npoints into flat_out.
  !> The HDF5 dataset NAME has shape [nrho*ntemp*nye].
  !============================================================================
  subroutine read_slab(file_id, dset_name, flat_out, nrho, ntemp, nye, ntables, iv0, error)
    use hdf5
    implicit none
    integer(HID_T),   intent(in)    :: file_id
    character(*),     intent(in)    :: dset_name
    double precision, intent(inout) :: flat_out(ntables * nrho * ntemp * nye)
    integer,          intent(in)    :: nrho, ntemp, nye, ntables, iv0
    integer,          intent(out)   :: error

    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims1d(1)
    integer          :: npts, base

    npts = nrho * ntemp * nye
    base = (iv0-1) * npts + 1

    dims1d(1) = npts
    call h5dopen_f(file_id, trim(dset_name), dset_id, error)
    if (error /= 0) return
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, flat_out(base:base+npts-1), dims1d, error)
    call h5dclose_f(dset_id, error)
  end subroutine read_slab

end module mod_eos_readtable_leptonic_scollapse
