!================================================================================
!
!  mod_eos_readtable_leptonic.t
!
!  Reads the dedicated leptonic HDF5 table produced by the Margherita
!  leptonic EOS library.  The file contains:
!
!    logrho_table        [nrho]                  ln(rho) in code units
!    logtemp_table       [ntemp]                 ln(T/MeV)
!    yle_table           [nyle]                  yle  (linear)
!    ymu_table           [nymu]                  ln(ymu)  (log scale!)
!    eos_ymumin/max      scalar bounds
!    eos_ylemin/max      scalar bounds
!    electronic_eos_tables [nrho*ntemp*nyle * NUM_VARS_ELE]
!    muonic_eos_tables     [nrho*ntemp*nymu * NUM_VARS_MUON]
!
!  After reading, data are stored in
!    lep_tables_ele (nrho, ntemp, nyle, lep_nvars_ele)
!    lep_tables_muon(nrho, ntemp, nymu, lep_nvars_muon)
!  with units already in code units.
!
!  Copyright (C) 2024  Harry Ho-Yin Ng
!
!================================================================================

module mod_eos_readtable_leptonic
  use mod_eos_leptonic_parameters
  implicit none
  public

contains

  subroutine activate_tablereader_leptonic
    call readtable_leptonic(leptonic_table_name)
  end subroutine activate_tablereader_leptonic

  !============================================================================
  subroutine readtable_leptonic(eos_filename)
    use hdf5
    use mod_eos_leptonic_parameters
    include 'amrvacdef.f'

    character(*), intent(in) :: eos_filename

    ! HDF5 handles
    integer(HID_T) :: file_id, dset_id
    integer(HSIZE_T) :: dims1(1), dims1d(1)
    integer :: error, accerr

    ! Temporary flat arrays (C ordering from HDF5)
    double precision, allocatable :: flat_ele(:)   ! nrho*ntemp*nyle*nvars_ele
    double precision, allocatable :: flat_muon(:)  ! nrho*ntemp*nymu*nvars_muon
    double precision, allocatable :: logrho_tmp(:), logtemp_tmp(:)
    double precision :: rho_axis_diff, temp_axis_diff

    integer :: nrho_loc, ntemp_loc, nyle_loc, nymu_loc
    integer :: i, j, k, iv, indold, indnew

    !--------------------------------------------------------------------------
    accerr = 0
    if (mype == 0) write(*,*) "Reading leptonic EOS table: ", trim(adjustl(eos_filename))

    call h5open_f(error)
    call h5fopen_f(trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)
    if (error /= 0) call mpistop("mod_eos_readtable_leptonic: cannot open HDF5 file")

    !--------------------------------------------------------------------------
    ! Read table dimensions
    !--------------------------------------------------------------------------
    dims1(1) = 1
    call h5dopen_f(file_id, "nrho",  dset_id, error); accerr=accerr+error
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho_loc,  dims1, error); accerr=accerr+error
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "ntemp", dset_id, error); accerr=accerr+error
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp_loc, dims1, error); accerr=accerr+error
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "nyle",  dset_id, error); accerr=accerr+error
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nyle_loc,  dims1, error); accerr=accerr+error
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "nymu",  dset_id, error); accerr=accerr+error
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nymu_loc,  dims1, error); accerr=accerr+error
    call h5dclose_f(dset_id, error)

    if (accerr /= 0) call mpistop("mod_eos_readtable_leptonic: error reading dimensions")

    if (mype == 0) write(*,'(a,4i6)') "  nrho ntemp nyle nymu: ", nrho_loc, ntemp_loc, nyle_loc, nymu_loc

    !--------------------------------------------------------------------------
    ! Store in parameter module globals (shared with baryon reader result)
    ! Only overwrite if not already set (baryon reader sets nrho/ntemp first)
    !--------------------------------------------------------------------------
    if (lep_nrho == 0) lep_nrho = nrho_loc
    if (lep_ntemp == 0) lep_ntemp = ntemp_loc
    if (lep_nrho /= nrho_loc .or. lep_ntemp /= ntemp_loc) then
      call mpistop("mod_eos_readtable_leptonic: inconsistent nrho/ntemp with baryon table")
    end if
    lep_nyle  = nyle_loc
    lep_nymu  = nymu_loc

    !--------------------------------------------------------------------------
    ! Allocate axis arrays
    !--------------------------------------------------------------------------
    allocate(logrho_tmp(lep_nrho))
    allocate(logtemp_tmp(lep_ntemp))
    if (.not. allocated(lep_yle_table))     allocate(lep_yle_table(lep_nyle))
    if (.not. allocated(lep_logymu_table))  allocate(lep_logymu_table(lep_nymu))

    dims1d(1) = lep_nrho
    call h5dopen_f(file_id, "logrho_table",  dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho_tmp,  dims1d, error)
    call h5dclose_f(dset_id, error)

    dims1d(1) = lep_ntemp
    call h5dopen_f(file_id, "logtemp_table", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp_tmp, dims1d, error)
    call h5dclose_f(dset_id, error)

    dims1d(1) = lep_nyle
    call h5dopen_f(file_id, "yle_table",     dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_yle_table,     dims1d, error)
    call h5dclose_f(dset_id, error)

    dims1d(1) = lep_nymu
    call h5dopen_f(file_id, "ymu_table",     dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_logymu_table,  dims1d, error)
    call h5dclose_f(dset_id, error)

    if (.not. allocated(lep_logrho_table) .or. .not. allocated(lep_logtemp_table)) then
      call mpistop("mod_eos_readtable_leptonic: baryon axes must be loaded first")
    end if
    rho_axis_diff = maxval(abs(logrho_tmp - lep_logrho_table))
    temp_axis_diff = maxval(abs(logtemp_tmp - lep_logtemp_table))
    if (rho_axis_diff > 1.0d-10) then
      if (mype == 0) write(*,'(a,es14.6)') "  max |d log rho| = ", rho_axis_diff
      call mpistop("mod_eos_readtable_leptonic: rho axis mismatch")
    end if
    if (temp_axis_diff > 1.0d-10) then
      if (mype == 0) write(*,'(a,es14.6)') "  max |d log temp| = ", temp_axis_diff
      call mpistop("mod_eos_readtable_leptonic: temp axis mismatch")
    end if

    !--------------------------------------------------------------------------
    ! Read scalar bounds
    !--------------------------------------------------------------------------
    dims1(1) = 1
    call h5dopen_f(file_id, "eos_ylemin", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_eos_ylemin, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "eos_ylemax", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_eos_ylemax, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "eos_ymumin", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_eos_ymumin, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "eos_ymumax", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lep_eos_ymumax, dims1, error)
    call h5dclose_f(dset_id, error)

    !--------------------------------------------------------------------------
    ! Read flat electronic table  [nrho * ntemp * nyle * nvars_ele]
    ! HDF5 layout (C order):  iv + nvars * (i + nrho * (j + ntemp * k))
    ! Target Fortran layout:  ft(i, j, k, iv)  i.e. Fortran column-major
    !--------------------------------------------------------------------------
    allocate(flat_ele(lep_nrho * lep_ntemp * lep_nyle * lep_nvars_ele))
    dims1d(1) = size(flat_ele)
    call h5dopen_f(file_id, "electronic_eos_tables", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, flat_ele, dims1d, error)
    call h5dclose_f(dset_id, error)

    if (allocated(lep_tables_ele)) deallocate(lep_tables_ele)
    allocate(lep_tables_ele(lep_nrho, lep_ntemp, lep_nyle, lep_nvars_ele))
    do iv = 1, lep_nvars_ele
    do k  = 1, lep_nyle
    do j  = 1, lep_ntemp
    do i  = 1, lep_nrho
      ! C index (0-based):  iv-1 + nvars * (i-1 + nrho * (j-1 + ntemp * (k-1)))
      indold = (iv-1) + lep_nvars_ele * ((i-1) + lep_nrho * ((j-1) + lep_ntemp * (k-1))) + 1
      lep_tables_ele(i, j, k, iv) = flat_ele(indold)
    end do; end do; end do; end do
    deallocate(flat_ele)

    !--------------------------------------------------------------------------
    ! Read flat muonic table  [nrho * ntemp * nymu * nvars_muon]
    !--------------------------------------------------------------------------
    allocate(flat_muon(lep_nrho * lep_ntemp * lep_nymu * lep_nvars_muon))
    dims1d(1) = size(flat_muon)
    call h5dopen_f(file_id, "muonic_eos_tables", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, flat_muon, dims1d, error)
    call h5dclose_f(dset_id, error)

    if (allocated(lep_tables_muon)) deallocate(lep_tables_muon)
    allocate(lep_tables_muon(lep_nrho, lep_ntemp, lep_nymu, lep_nvars_muon))
    do iv = 1, lep_nvars_muon
    do k  = 1, lep_nymu
    do j  = 1, lep_ntemp
    do i  = 1, lep_nrho
      indold = (iv-1) + lep_nvars_muon * ((i-1) + lep_nrho * ((j-1) + lep_ntemp * (k-1))) + 1
      lep_tables_muon(i, j, k, iv) = flat_muon(indold)
    end do; end do; end do; end do
    deallocate(flat_muon)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    deallocate(logrho_tmp, logtemp_tmp)

    if (mype == 0) then
      write(*,'(a,2es14.6)') "  yle  range: ", lep_eos_ylemin, lep_eos_ylemax
      write(*,'(a,2es14.6)') "  ymu  range: ", lep_eos_ymumin, lep_eos_ymumax
      write(*,*) "  Leptonic table loaded successfully."
    end if

  end subroutine readtable_leptonic

end module mod_eos_readtable_leptonic
