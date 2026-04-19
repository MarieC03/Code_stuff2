module mod_c2b_cart3d
    use HDF5
    implicit none

    ! note that this module need API of hdf5, so specific version of hdf5 may be needed, version 1.10.5 and 1.12.2 has been proved to work
    public
    logical                                           :: debug_in_mode_c2b_cart3d = .false. !  open it for verbose output for debuging
    type carpet_single_datablock
        double precision, allocatable                 :: data(:, :, :)
        double precision                              :: coord_origin(3), delta_x(3)
        double precision                              :: coord_min(3), coord_max(3)
        double precision                              :: interp_range_min(3), interp_range_max(3)
        character*256                                 :: dset_path
        integer                                       :: shape(3), nghostzone(3), boundary_box(6)
    end type carpet_single_datablock

    type carpet_single_reflevel
        type(carpet_single_datablock), allocatable    :: all_datablocks(:) ! there are many datablocks separately stored in many different files
        integer                                       :: num_datablocks, ref_level
        double precision                              :: coord_min(3), coord_max(3), delta_x(3)
        double precision                              :: interp_range_min(3), interp_range_max(3)
    end type carpet_single_reflevel

    type carpet_variable
        type(carpet_single_reflevel), allocatable     :: all_refine_levels(:) ! real data is stored here or in volatile found_dabk is dependent on save_cpu_or_mem
        integer                                       :: save_cpu_or_mem      ! 1 for save cpu, 2 for save memory
        character*256                                 :: fpath, vfname, vname ! name of the variable in hdf5 dataset and filename
        character*3                                   :: staggered            ! list of 'T' or 'F' to specify whether staggered or not in each dimension, eg, 'TFT'
        integer                                       :: num_refine_levels    ! total number of refinement level
        integer                                       :: iter_num             ! iteration number
        integer                                       :: found_ref_level      ! found point in which refinement level
        integer                                       :: found_dabk_num       ! found point in which datablock
        type(carpet_single_datablock)                 :: found_dabk           ! datablock which must contains the point to interpolate, exist and store data only in the case that save_cpu_or_mem=2
    end type carpet_variable

    contains


    subroutine init_single_variable(cp3dvb, fpath, variable_fname, variable_name, num_refine_levels, num_datablocks, iter_num, staggered, save_cpu_or_mem)
        type(carpet_variable), intent(inout)          :: cp3dvb
        character(len=256), intent(in)                :: fpath ! path to directory where hdf5 files are stored
        character(len=256), intent(in)                :: variable_fname ! name of the file that store the variable, only first few characters are needed
        character(len=256), intent(in)                :: variable_name ! name of the variable in the dataset, can also be found in the ETK par file
        integer, intent(in)                           :: num_refine_levels, num_datablocks, iter_num
        integer, intent(in)                           :: save_cpu_or_mem
        character(len=3), intent(in)                  :: staggered ! whether staggered or not in each direction
        ! local variables
        character(len=256)                            :: local_fp
        integer                                       :: i, len_fp

        ! check allocation state, if allocated, deallocate and allocate again
        if (allocated(cp3dvb%all_refine_levels)) then
            write(*, *) "Variable ", cp3dvb%vname, " already initialized before, deallocate and ini again!"
            call destroy_single_variable(cp3dvb)
        endif
        allocate(cp3dvb%all_refine_levels(0: num_refine_levels-1))
        ! add a trailing '/' in fpath if there does not have one
        len_fp = len_trim(fpath)
        if (fpath(len_fp:len_fp)/='/') then
            local_fp = trim(fpath)//'/'
        else
            local_fp = fpath
        end if
        ! initialize
        do i=0, num_refine_levels-1
            call init_single_reflevel(cp3dvb%all_refine_levels(i), num_datablocks, i, local_fp, variable_fname, variable_name, iter_num, staggered, save_cpu_or_mem)
        enddo
        cp3dvb%save_cpu_or_mem = save_cpu_or_mem
        cp3dvb%fpath = local_fp
        cp3dvb%vfname = variable_fname
        cp3dvb%vname = variable_name
        cp3dvb%staggered = staggered
        cp3dvb%num_refine_levels = num_refine_levels
        cp3dvb%iter_num = iter_num
        cp3dvb%found_ref_level = -1
        cp3dvb%found_dabk_num = -1
    end subroutine init_single_variable

    subroutine init_single_reflevel(sgrefl, num_datablocks, ref_level, fpath, variable_fname, variable_name, iter_num, staggered, save_cpu_or_mem)
        type(carpet_single_reflevel), intent(inout)   :: sgrefl
        integer, intent(in)                           :: num_datablocks, ref_level, iter_num
        integer, intent(in)                           :: save_cpu_or_mem
        character(len=256), intent(in)                :: fpath, variable_fname, variable_name
        character(len=3), intent(in)                  :: staggered ! whether staggered or not in each direction
        ! local variables below
        integer                                       :: i, j

        if (allocated(sgrefl%all_datablocks)) then
            write(*, "(A17, I3, A54)") "Refinement level ", sgrefl%ref_level, " already initialized before, deallocate and ini again!"
            call destroy_single_reflevel(sgrefl)
        endif
        allocate(sgrefl%all_datablocks(0: num_datablocks-1))
        sgrefl%interp_range_max = -1.0d+99
        sgrefl%interp_range_min = 1.0d+99
        sgrefl%coord_max = -1.0d+99
        sgrefl%coord_min = 1.0d+99
        do i=0, num_datablocks-1
            if (save_cpu_or_mem==1) then
                call init_single_datablock(sgrefl%all_datablocks(i), fpath, variable_fname, variable_name, iter_num, ref_level, i, staggered, .true.)
            else
                call init_single_datablock(sgrefl%all_datablocks(i), fpath, variable_fname, variable_name, iter_num, ref_level, i, staggered, .false.)
            endif
            loop_dim: do j=1, 3
                sgrefl%interp_range_max(j) = max(sgrefl%interp_range_max(j), sgrefl%all_datablocks(i)%interp_range_max(j))
                sgrefl%interp_range_min(j) = min(sgrefl%interp_range_min(j), sgrefl%all_datablocks(i)%interp_range_min(j))
                sgrefl%coord_max(j) = max(sgrefl%coord_max(j), sgrefl%all_datablocks(i)%coord_max(j))
                sgrefl%coord_min(j) = min(sgrefl%coord_min(j), sgrefl%all_datablocks(i)%coord_min(j))
                sgrefl%delta_x(j) = sgrefl%all_datablocks(i)%delta_x(j)   ! assign the same value multiple times but still ok
            enddo loop_dim
        enddo
        sgrefl%num_datablocks = num_datablocks
        sgrefl%ref_level = ref_level
    end subroutine init_single_reflevel

    ! Note: if any of the default behaviors of the naming system of ETK is changed, this part should be updated accordingly 
    subroutine init_single_datablock(sgdabk, fpath, variable_fname, variable_name, iter_num, ref_level, block_id, staggered, read_data)
        type(carpet_single_datablock), intent(inout)  :: sgdabk
        character(len=256), intent(in)                :: fpath, variable_fname, variable_name
        integer, intent(in)                           :: ref_level, block_id, iter_num
        character(len=3), intent(in)                  :: staggered ! whether staggered or not in each direction
        logical, intent(in)                           :: read_data ! whether to read data or just get the structure information
        ! local variables
        character(len=256)                            :: buff1, buff2, buff_f, fname
        integer                                       :: i
        logical                                       :: file_exist

        ! create filename
        write(buff_f, *) block_id
        fname = trim(adjustl(fpath))//trim(adjustl(variable_fname))//".xyz.file_"//trim(adjustl(buff_f))//".h5" ! Note: this could be '.hdf5' in the future
        if (debug_in_mode_c2b_cart3d) write(*, *) "Reading file: ", fname
        inquire(file=fname, exist=file_exist)
        if (.not. file_exist) call mpistop("ETK file "//fname//" not exist, please check the path, variable spelling, and the data block numbers.")
        ! create dataset name based on many things
        write(buff1, *) iter_num
        ! Note: I do not know what means 'tl', so the next line may have to be modified in the future to adapt to the real situation
        buff2 = trim(adjustl(variable_name))//" it="//trim(adjustl(buff1))//" tl=0"
        write(buff1, *) ref_level
        buff2 = trim(buff2)//" rl="//trim(adjustl(buff1))
        write(buff1, *) block_id
        buff2 = trim(buff2)//" c="//trim(adjustl(buff1))
        sgdabk%dset_path = buff2
        ! read in data and important attributes
        call hdf5_read_attr_single_datablock(fname, sgdabk)
        if (read_data) call hdf5_read_single_dataset(fname, sgdabk%dset_path, sgdabk%data)
        ! calculate boundaries and many other properties
        do i=1, 3
            if(staggered(i:i)=='T') sgdabk%coord_origin(i) = sgdabk%coord_origin(i)+0.5*sgdabk%delta_x(i)
            sgdabk%coord_min(i) = sgdabk%coord_origin(i)
            sgdabk%coord_max(i) = sgdabk%coord_origin(i)+sgdabk%delta_x(i)*sgdabk%shape(i)
            if (sgdabk%boundary_box(i*2-1)==1) then
                sgdabk%interp_range_min(i) = sgdabk%coord_min(i)
            else
                sgdabk%interp_range_min(i) = sgdabk%coord_min(i)+sgdabk%nghostzone(i)*sgdabk%delta_x(i)
            endif
            if (sgdabk%boundary_box(i*2)==1) then
                sgdabk%interp_range_max(i) = sgdabk%coord_max(i)-sgdabk%delta_x(i)
            else
                sgdabk%interp_range_max(i) = sgdabk%coord_max(i)-sgdabk%nghostzone(i)*sgdabk%delta_x(i)
            endif
        enddo
    end subroutine init_single_datablock

    subroutine load_selected_ref_and_datablock(cp3dvb, nref, ndabk, mask_array)
        type(carpet_variable), intent(inout)          :: cp3dvb
        integer, intent(in)                           :: nref, ndabk
        integer, intent(in)                           :: mask_array(0:nref-1, 0:ndabk-1)
        ! local
        integer                                       :: i, j

        do i=0, nref-1
            if (allocated(cp3dvb%all_refine_levels)) then
                associate(sgrefl=>cp3dvb%all_refine_levels(i))
                    do j=0, ndabk-1
                        associate(sgdabk=>sgrefl%all_datablocks(j))
                            if (mask_array(i, j)==1 .and. (.not. allocated(sgdabk%data))) then
                                call init_single_datablock(sgdabk, cp3dvb%fpath, cp3dvb%vfname, cp3dvb%vname, &
                                    cp3dvb%iter_num, i, j, cp3dvb%staggered, .true.)
                            endif
                            ! the following will release un-necessary data blocks
                            if (mask_array(i, j)==0 .and. allocated(sgdabk%data)) deallocate(sgdabk%data)
                        end associate
                    enddo
                end associate
            endif
        enddo
    end subroutine load_selected_ref_and_datablock

    subroutine destroy_single_variable(cp3dvb)
        type(carpet_variable), intent(inout)          :: cp3dvb
        integer                                       :: i, status
        
        if (allocated(cp3dvb%all_refine_levels)) then
            do i=0, cp3dvb%num_refine_levels-1
                if (debug_in_mode_c2b_cart3d) write(*, *) "Destroying variable "//trim(adjustl(cp3dvb%vname))//" in refinement level: ", i
                call destroy_single_reflevel(cp3dvb%all_refine_levels(i))
            enddo
            deallocate(cp3dvb%all_refine_levels, stat=status)
        endif
        if (allocated(cp3dvb%found_dabk%data)) deallocate(cp3dvb%found_dabk%data)
        cp3dvb%num_refine_levels = 0
        cp3dvb%iter_num = -1
        cp3dvb%found_ref_level = -1
        cp3dvb%found_dabk_num = -1
    end subroutine destroy_single_variable

    subroutine destroy_single_reflevel(sgrefl)
        type(carpet_single_reflevel), intent(inout)   :: sgrefl
        integer                                       :: i, status
        if (allocated(sgrefl%all_datablocks)) then
            do i=0, sgrefl%num_datablocks-1
                if (debug_in_mode_c2b_cart3d) write(*, "(A25, I4, A10, I4)") "    Destroying ref_level ", sgrefl%ref_level, ", datablock", i
                if (allocated(sgrefl%all_datablocks(i)%data)) deallocate(sgrefl%all_datablocks(i)%data, stat=status)
            enddo
            deallocate(sgrefl%all_datablocks, stat=status)
        endif
        sgrefl%num_datablocks = 0
        sgrefl%ref_level = -1
    end subroutine destroy_single_reflevel

    ! This subroutine read the useful attribute of a single datablock
    subroutine hdf5_read_attr_single_datablock(fname, sgdabk)
        use hdf5
        type(carpet_single_datablock), intent(inout)  :: sgdabk         ! single data block in one refinement level
        character(len=256), intent(in)                :: fname
        ! local variables below
        integer(hsize_t)                              :: dims(1) = [1], maxdims(3), ds_dim(3)
        integer(hid_t)                                :: file_id, attr1_id, attr2_id, attr3_id, attr4_id
        integer(hid_t)                                :: dataset_id, dataspace_id
        integer(hid_t)                                :: attr1_type, attr2_type, attr3_type, attr4_type
        integer                                       :: hdferr, i

        ! open
        call H5Open_f(hdferr)
        call H5FOpen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferr)
        call H5DOpen_f(file_id, sgdabk%dset_path, dataset_id, hdferr)
        ! read attribute delta_x
        call H5Aopen_f(dataset_id, "delta", attr1_id, hdferr)
        call H5AGet_type_f(attr1_id, attr1_type, hdferr)
        call H5ARead_f(attr1_id, attr1_type, sgdabk%delta_x, dims, hdferr) ! type: H5T_NATIVE_DOUBLE
        call H5AClose_f(attr1_id, hdferr)
        ! read attribute origin
        call H5Aopen_f(dataset_id, "origin", attr2_id, hdferr)
        call H5AGet_type_f(attr2_id, attr2_type, hdferr)
        call H5ARead_f(attr2_id, attr2_type, sgdabk%coord_origin, dims, hdferr)
        call H5AClose_f(attr2_id, hdferr)
        ! read attribute of ghostzones
        call H5Aopen_f(dataset_id, "cctk_nghostzones", attr3_id, hdferr)
        call H5AGet_type_f(attr3_id, attr3_type, hdferr)
        call H5ARead_f(attr3_id, attr3_type, sgdabk%nghostzone, maxdims, hdferr)
        call H5AClose_f(attr3_id, hdferr)
        ! read attribute of boundary of box
        call H5Aopen_f(dataset_id, "cctk_bbox", attr4_id, hdferr)
        call H5AGet_type_f(attr4_id, attr4_type, hdferr)
        call H5ARead_f(attr4_id, attr4_type, sgdabk%boundary_box, maxdims, hdferr)
        call H5AClose_f(attr4_id, hdferr)
        ! read dataset shape
        call H5DGet_space_f(dataset_id, dataspace_id, hdferr)
        call H5SGet_simple_extent_dims_f(dataspace_id, ds_dim, maxdims, hdferr)
        do i=1, 3
            sgdabk%shape(i) = ds_dim(i)
        enddo
        call H5SClose_f(dataspace_id, hdferr)
        ! close files
        call H5FClose_f(file_id, hdferr)
        call H5Close_f(hdferr) 
    end subroutine hdf5_read_attr_single_datablock

    subroutine hdf5_read_single_dataset(fname, dset_path, dset)
        double precision, allocatable, intent(inout)  :: dset(:, :, :)
        character(len=256), intent(in)                :: fname, dset_path 
        ! local variables below
        integer(hid_t)                                :: file_id, dataset_id, dataspace_id, dataset_type
        integer(hsize_t)                              :: dims(3), maxdims(3)
        integer                                       :: hdferr, alloc_status
        
        call H5Open_f(hdferr)     ! initialize the hdf5
        call H5FOpen_f(fname, H5F_ACC_RDONLY_F, file_id, hdferr)
        call H5DOpen_f(file_id, dset_path, dataset_id, hdferr)
        ! get dimension
        call H5DGet_space_f(dataset_id, dataspace_id, hdferr)
        call H5SGet_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
        if (debug_in_mode_c2b_cart3d) then
            write(*, 1000) trim(dset_path), dims
            1000 format( "    reading dataset named '", A, "' which has dimension: ", 3I5)
        end if
        ! allocate and read dataset
        if (allocated(dset)) then 
            write(*, *) "dataset already exist, deallocate it: ", dset_path
            deallocate(dset)
        endif
        allocate(dset(dims(1), dims(2), dims(3)), STAT=alloc_status)
        !call H5DGet_type_f(dataset_id, dataset_type, hdferr) ! can't use original, should be transformed to aim array type
        call H5DRead_f(dataset_id, H5T_NATIVE_DOUBLE, dset, dims, hdferr)
        ! close files
        call H5DClose_f(dataset_id, hdferr)
        call H5FClose_f(file_id, hdferr)
        call H5Close_f(hdferr)    ! release and clean
    end subroutine hdf5_read_single_dataset

    subroutine output_basic_info_single_reflevel(sgrefl)
        type(carpet_single_reflevel), intent(in)        :: sgrefl
        integer                                         :: i, j
        do i=0, sgrefl%num_datablocks-1
            associate(sgdabk=>sgrefl%all_datablocks(i))
                do j=1, 3
                    write(*, 2000) j, sgdabk%shape(j), sgdabk%interp_range_min(j), sgdabk%interp_range_max(j), sgdabk%delta_x(j)
                    2000 format ( "        direction ", I2, ": shape ", I5, ", interp_range_min ", ES12.4, ", interp_range_max", ES12.4, &
                                 ", delta_x ", ES12.4)
                enddo
            end associate
        enddo
        write(*, *) NEW_LINE('a')
    end subroutine output_basic_info_single_reflevel

    integer function find_point_in_which_level(cp3dvb, coord_x, coord_y, coord_z)
        type(carpet_variable)                         :: cp3dvb
        double precision                              :: coord_x, coord_y, coord_z
        integer                                       :: i
        logical                                       :: point_is_here

        ! start from the finest refinement level to find the point
        loop_reflevel: do i=cp3dvb%num_refine_levels-1, 0, -1
            associate(curent_lv=>cp3dvb%all_refine_levels(i))
                point_is_here = .true.
                point_is_here = point_is_here .and. ((curent_lv%interp_range_min(1)<=coord_x) .and. (curent_lv%interp_range_max(1)>=coord_x))
                point_is_here = point_is_here .and. ((curent_lv%interp_range_min(2)<=coord_y) .and. (curent_lv%interp_range_max(2)>=coord_y))
                point_is_here = point_is_here .and. ((curent_lv%interp_range_min(3)<=coord_z) .and. (curent_lv%interp_range_max(3)>=coord_z))
                if (debug_in_mode_c2b_cart3d) write(*, *) 'level: ', i, curent_lv%interp_range_min(1), curent_lv%interp_range_max(1), curent_lv%interp_range_min(2), curent_lv%interp_range_max(2), curent_lv%interp_range_min(3), curent_lv%interp_range_max(3)
                if (point_is_here) then
                    find_point_in_which_level = i
                    exit loop_reflevel
                endif
            end associate
        enddo loop_reflevel
        ! point is not in the whole grid of this variable
        if (.not. point_is_here) then
            find_point_in_which_level = -1
        endif
    end function find_point_in_which_level

    integer function find_point_in_which_datablock(sgrefl, coord_x, coord_y, coord_z)
        type(carpet_single_reflevel)                  :: sgrefl
        double precision                              :: coord_x, coord_y, coord_z
        integer                                       :: i
        logical                                       :: point_is_here

        ! start from the finest refinement level to find the point
        loop_datablock: do i=0, sgrefl%num_datablocks-1
            associate(curent_db=>sgrefl%all_datablocks(i))
                point_is_here = .true.
                point_is_here = point_is_here .and. ((curent_db%interp_range_min(1)<=coord_x) .and. (curent_db%interp_range_max(1)>=coord_x))
                point_is_here = point_is_here .and. ((curent_db%interp_range_min(2)<=coord_y) .and. (curent_db%interp_range_max(2)>=coord_y))
                point_is_here = point_is_here .and. ((curent_db%interp_range_min(3)<=coord_z) .and. (curent_db%interp_range_max(3)>=coord_z))
                if (debug_in_mode_c2b_cart3d) write(*, *) 'block: ', i,  curent_db%interp_range_min(1), curent_db%interp_range_max(1), curent_db%interp_range_min(2), curent_db%interp_range_max(2), curent_db%interp_range_min(3), curent_db%interp_range_max(3)
                if (point_is_here) then
                    find_point_in_which_datablock = i
                    exit loop_datablock
                endif
            end associate
        enddo loop_datablock
        ! point is not in the whole grid of in this refinement level
        if (.not. point_is_here) then 
            find_point_in_which_datablock = -1
        endif
    end function find_point_in_which_datablock

    subroutine find_point(coord_x, coord_y, coord_z, cp3dvb, which_level, which_datablock, my_error_type)
        integer, intent(out)                          :: which_level, which_datablock
        integer, intent(out)                          :: my_error_type ! 0 for normal, 1 for out of bound
        type(carpet_variable), intent(inout)          :: cp3dvb
        double precision, intent(in)                  :: coord_x, coord_y, coord_z
        ! local variables

        ! initialize
        which_level=-1
        which_datablock=-1
        my_error_type = 0
        ! find
        which_level = find_point_in_which_level(cp3dvb, coord_x, coord_y, coord_z)
        if (which_level==-1) then
            my_error_type = 1  ! not in any case
        else
            loop_find_point: do while (which_datablock==-1)
                which_datablock = find_point_in_which_datablock(cp3dvb%all_refine_levels(which_level), coord_x, coord_y, coord_z)
                if (which_datablock==-1) then
                    if (debug_in_mode_c2b_cart3d) then
                        write(*, "(A, ES12.4, A2, ES12.4, A2, ES12.4, A, I4)") "Warning: When interpolating variable "//trim(adjustl(cp3dvb%vname))//" at point: (", &
                                coord_x, ', ', coord_y, ', ', coord_z , ", going down to one lower refinement level because it is not in level: ", which_level
                    endif
                    which_level = which_level-1
                endif
                if (which_level<0) then
                    exit loop_find_point
                endif
            enddo loop_find_point
            if (which_level==-1) then
                my_error_type = 2 ! in some boxes, but loop do not find it
            else
                if (which_level/=cp3dvb%found_ref_level .or. which_datablock/=cp3dvb%found_dabk_num) then
                    if (cp3dvb%save_cpu_or_mem==2) call init_single_datablock(cp3dvb%found_dabk, cp3dvb%fpath, cp3dvb%vfname, cp3dvb%vname, cp3dvb%iter_num, &
                        which_level, which_datablock, cp3dvb%staggered, .true.)
                    cp3dvb%found_ref_level = which_level
                    cp3dvb%found_dabk_num = which_datablock
                endif
                !write(*, "(A10, I3, A12, I3)") "In level: ", which_level, ", in block: ", which_datablock
                my_error_type = 0
            endif
        endif
    end subroutine find_point

    subroutine interp_cart3d_carpet_variable_one_point(point_value, coord_x, coord_y, coord_z, cp3dvb, type_interp, ob_value)
        double precision, intent(out)                 :: point_value
        type(carpet_variable), intent(inout)          :: cp3dvb
        double precision, intent(in)                  :: coord_x, coord_y, coord_z
        double precision, intent(in)                  :: ob_value ! out of bound value, if not find in any grid, return it
        character*256, intent(in)                     :: type_interp
        ! local variables
        integer                                       :: found_error_type
        integer                                       :: found_level, found_datablock

        call find_point(coord_x, coord_y, coord_z, cp3dvb, found_level, found_datablock, found_error_type)
        if (found_error_type==0) then
            select case(cp3dvb%save_cpu_or_mem)
            case (1)
                call interp_one_point_from_dataset(point_value, coord_x, coord_y, coord_z, &
                    cp3dvb%all_refine_levels(found_level)%all_datablocks(found_datablock), type_interp)
            case (2)
                call interp_one_point_from_dataset(point_value, coord_x, coord_y, coord_z, cp3dvb%found_dabk, type_interp)
            case default ! default the best choice
                if (.not. allocated(cp3dvb%all_refine_levels(found_level)%all_datablocks(found_datablock)%data)) then
                    call init_single_datablock(cp3dvb%all_refine_levels(found_level)%all_datablocks(found_datablock), &
                        cp3dvb%fpath, cp3dvb%vfname, cp3dvb%vname, cp3dvb%iter_num, found_level, found_datablock, cp3dvb%staggered, .true.)
                    if (debug_in_mode_c2b_cart3d) write(*, *) "Unefficient warning: dynamically loading datablock is very unefficient!", cp3dvb%vname
                endif
                call interp_one_point_from_dataset(point_value, coord_x, coord_y, coord_z, &
                    cp3dvb%all_refine_levels(found_level)%all_datablocks(found_datablock), type_interp)
            end select
            ! NaN check
            if (point_value .ne. point_value) then 
                write(*, "(A, ES12.4, A2, ES12.4, A2, ES12.4, A13, I4, A13, I4)") &
                    "NaN occured when interpolating variable "//trim(adjustl(cp3dvb%vname))//" at point: (", coord_x, ', ', &
                    coord_y, ', ', coord_z , ", ref level: ", found_level, ", datablock: ", found_datablock
                point_value = ob_value
            endif
        else
             write(*, "(A, ES12.4, A2, ES12.4, A2, ES12.4, A, ES12.4, A14, I4)") "Warning: for variable "//trim(adjustl(cp3dvb%vname))//", can't find point(", &
                coord_x, ', ', coord_y, ', ', coord_z , " in any refinement level, return out of bound value: ", ob_value, ", Error type: ", found_error_type
            point_value = ob_value
        endif
    end subroutine interp_cart3d_carpet_variable_one_point

    subroutine interp_one_point_from_dataset(point_value, coord_x, coord_y, coord_z, sgdabk, type_interp)
        double precision, intent(out)                 :: point_value
        double precision, intent(in)                  :: coord_x, coord_y, coord_z
        type(carpet_single_datablock), intent(in)     :: sgdabk
        character*256, intent(in)                     :: type_interp
        ! local variables
        integer                                       :: idx_x, idx_y, idx_z
        double precision                              :: xdbuff, ydbuff, zdbuff, xrd, yrd, zrd   ! distance to the floor point normalized by the grid size
        double precision                              :: average, rel_diff
        ! tunable parameters
        double precision                              :: threshold_4th=0.01, threshold_lin=0.1

        xdbuff = (coord_x-sgdabk%coord_min(1))/sgdabk%delta_x(1)+1
        ydbuff = (coord_y-sgdabk%coord_min(2))/sgdabk%delta_x(2)+1
        zdbuff = (coord_z-sgdabk%coord_min(3))/sgdabk%delta_x(3)+1
        idx_x = floor(xdbuff)
        idx_y = floor(ydbuff)
        idx_z = floor(zdbuff)
        xrd = xdbuff-idx_x
        yrd = ydbuff-idx_y
        zrd = zdbuff-idx_z
        select case (type_interp)
        case('Linear')
            call linear_interp_3d(point_value, sgdabk%data(idx_x:idx_x+1, idx_y:idx_y+1, idx_z:idx_z+1), xrd, yrd, zrd)
        case('Lagrangian2')
            if ((idx_x-1<=0) .or. (idx_y-1<=0) .or. (idx_z-1<=0)) then
                call mpistop("Please try to lower down the interpolation order because you are at the boundary!")
            else
                call lagrangian_interp_high_order_3d(point_value, sgdabk%data(idx_x-1:idx_x+1, idx_y-1:idx_y+1, idx_z-1:idx_z+1), xrd, yrd, zrd, -1, 1)
            endif
        case('Lagrangian3')
            if ((idx_x-1<=0) .or. (idx_x+2)>sgdabk%shape(1) .or. (idx_y-1<=0) .or. (idx_y+2)>sgdabk%shape(2) .or. (idx_z-1<=0) .or. (idx_z+2)>sgdabk%shape(3)) then
                !call mpistop("Please try to lower down the interpolation order because you are at the boundary!")
                call linear_interp_3d(point_value, sgdabk%data(idx_x:idx_x+1, idx_y:idx_y+1, idx_z:idx_z+1), xrd, yrd, zrd)
            else
                call lagrangian_interp_high_order_3d(point_value, sgdabk%data(idx_x-1:idx_x+2, idx_y-1:idx_y+2, idx_z-1:idx_z+2), xrd, yrd, zrd, -1, 2)
            endif
        case('Lagrangian4') ! Warning: higher order method needs more points, ghost zone should be larger enough to store these points (at least 3)
            if ((idx_x-2<=0) .or. (idx_x+2)>sgdabk%shape(1) .or. (idx_y-2<=0) .or. (idx_y+2)>sgdabk%shape(2) .or. (idx_z-2<=0) .or. (idx_z+2)>sgdabk%shape(3)) then
                call mpistop("Please try to lower down the interpolation order because you are at the boundary!")
            else
                call lagrangian_interp_high_order_3d(point_value, sgdabk%data(idx_x-2:idx_x+2, idx_y-2:idx_y+2, idx_z-2:idx_z+2), xrd, yrd, zrd, -2, 2)
            endif
        case('Lagrangian5') ! Warning: higher order method needs more points, ghost zone should be larger enough to store these points (at least 4)
            if ((idx_x-2<=0) .or. (idx_x+3)>sgdabk%shape(1) .or. (idx_y-2<=0) .or. (idx_y+3)>sgdabk%shape(2) .or. (idx_z-2<=0) .or. (idx_z+3)>sgdabk%shape(3)) then
                !call mpistop("Please try to lower down the interpolation order because you are at the boundary!")
                call linear_interp_3d(point_value, sgdabk%data(idx_x:idx_x+1, idx_y:idx_y+1, idx_z:idx_z+1), xrd, yrd, zrd)
            else
                call lagrangian_interp_high_order_3d(point_value, sgdabk%data(idx_x-2:idx_x+3, idx_y-2:idx_y+3, idx_z-2:idx_z+3), xrd, yrd, zrd, -2, 3)
            endif
        case('Hermitian3')
            if ((idx_x-1<=0) .or. (idx_x+2)>sgdabk%shape(1) .or. (idx_y-1<=0) .or. (idx_y+2)>sgdabk%shape(2) .or. (idx_z-1<=0) .or. (idx_z+2)>sgdabk%shape(3)) then
                !call mpistop("Please try to lower down the interpolation order because you are at the boundary!")
                call linear_interp_3d(point_value, sgdabk%data(idx_x:idx_x+1, idx_y:idx_y+1, idx_z:idx_z+1), xrd, yrd, zrd)
            else
                call Hermitian_interp_3order_3d(point_value, sgdabk%data(idx_x-1:idx_x+2, idx_y-1:idx_y+2, idx_z-1:idx_z+2), & 
                        xrd, yrd, zrd, sgdabk%delta_x(1), sgdabk%delta_x(2), sgdabk%delta_x(3))
            endif
        case('Automatic') ! Strategy of Armengol et, al., 22 [2112.09817v2]
            ! calculate the difference of the point with the surrounding average
            call linear_interp_3d(point_value, sgdabk%data(idx_x:idx_x+1, idx_y:idx_y+1, idx_z:idx_z+1), xrd, yrd, zrd)
            average = sum(sgdabk%data(idx_x-1:idx_x+2, idx_y-1:idx_y+2, idx_z-1:idx_z+2))/64.
            rel_diff = abs(point_value-average)/(average+tiny(1.0d0))
            ! select interpolation regime according to the relative difference
            if ((idx_x-2<=0) .or. (idx_x+2)>sgdabk%shape(1) .or. (idx_y-2<=0) .or. (idx_y+2)>sgdabk%shape(2) .or. (idx_z-2<=0) .or. (idx_z+2)>sgdabk%shape(3)) then
                call mpistop("Please try to lower down the interpolation order because you are at the boundary!")
            else
                if(rel_diff<threshold_4th) then
                    call lagrangian_interp_high_order_3d(point_value, sgdabk%data(idx_x-2:idx_x+2, idx_y-2:idx_y+2, idx_z-2:idx_z+2), xrd, yrd, zrd, -2, 2)
                else if (rel_diff<threshold_lin) then
                    call lagrangian_interp_high_order_3d(point_value, sgdabk%data(idx_x-1:idx_x+1, idx_y-1:idx_y+1, idx_z-1:idx_z+1), xrd, yrd, zrd, -1, 1)
                endif
            endif
        case default
            call mpistop("Unknown interpolation type in subroutine interp_cart3d_carpet_variable_one_point")
        end select
        ! NaN check
        if (point_value .ne. point_value) then
            write(*, "(A24, ES12.4, A2, ES12.4, A2, ES12.4, A15, 3ES12.4, A15, 3ES12.4, A16, 3ES12.4, A2)") &
                "NaN occured: at point: (", coord_x, ', ', coord_y, ', ', coord_z &
                , "), coord_min: (", sgdabk%coord_min(1), sgdabk%coord_min(2), sgdabk%coord_min(3) &
                , "), coord_max: (", sgdabk%coord_max(1), sgdabk%coord_max(2), sgdabk%coord_max(3) &
                , "), interp_min: (", sgdabk%interp_range_min(1), sgdabk%interp_range_min(2), sgdabk%interp_range_min(3), ")."
            write(*, "(A17, 3ES12.4, A16, 3I4, A13, 3ES12.4, A9, 3I4, A2)") &
                , "    interp_max: (", sgdabk%interp_range_max(1), sgdabk%interp_range_max(2), sgdabk%interp_range_max(3) &
                , "), data Shape: (", sgdabk%shape(1), sgdabk%shape(2), sgdabk%shape(3) &
                , "), delta_x: (", sgdabk%delta_x(1), sgdabk%delta_x(2), sgdabk%delta_x(3) &
                , "), idx: (", idx_x, idx_y, idx_z, ")."
        endif
    end subroutine interp_one_point_from_dataset

    double precision function derivatives_1st(u, h) ! error: ~h^2
        double precision, intent(in)             :: u(-1:1)
        double precision, intent(in)             :: h

        derivatives_1st = (u(1)-u(-1))/(2*h)
    end function derivatives_1st

    double precision function derivatives_2nd(u, h) ! error: ~h^2
        double precision, intent(in)             :: u(-1:1)
        double precision, intent(in)             :: h

        derivatives_2nd = (u(1)+u(-1)-2*u(0))/h**2
    end function derivatives_2nd

    double precision function derivatives_mixed_2nd(u, hx, hy)
        double precision, intent(in)             :: u(-1:1, -1:1)
        double precision, intent(in)             :: hx, hy

        derivatives_mixed_2nd = (u(1, 1)-u(-1, 1)-u(1, -1)+u(-1, -1))/(4*hx*hy)
    end function derivatives_mixed_2nd

    double precision function derivatives_mixed_3rd(u, hx, hy, hz)
        double precision, intent(in)             :: u(-1:1, -1:1, -1:1)
        double precision, intent(in)             :: hx, hy, hz

        derivatives_mixed_3rd = (u(1, 1, 1)+u(1, -1, -1)+u(-1, -1, 1)+u(-1, 1, -1)-u(-1, 1, 1)-u(1, -1, 1)-u(1, 1, -1)-u(-1, -1, -1))/(8*hx*hy*hz)
    end function derivatives_mixed_3rd

    subroutine linear_interp_3d(out_value, u, nmd_x, nmd_y, nmd_z)
        double precision, intent(out)            :: out_value
        double precision, intent(in)             :: u(0:1, 0:1, 0:1)
        double precision, intent(in)             :: nmd_x, nmd_y, nmd_z ! distance to the center grid normalized by grid size, must in (0, 1)

        out_value = u(0, 0, 0)*(1-nmd_x)*(1-nmd_y)*(1-nmd_z)+u(1, 0, 0)*nmd_x*(1-nmd_y)*(1-nmd_z)+&
                    u(0, 1, 0)*(1-nmd_x)*nmd_y*(1-nmd_z)+u(0, 0, 1)*(1-nmd_x)*(1-nmd_y)*nmd_z+u(1, 1, 0)*nmd_x*nmd_y*(1-nmd_z)+&
                    u(1, 0, 1)*nmd_x*(1-nmd_y)*nmd_z+u(0, 1, 1)*(1-nmd_x)*nmd_y*nmd_z+u(1, 1, 1)*nmd_x*nmd_y*nmd_z

    end subroutine linear_interp_3d

    ! see Table in https://en.wikipedia.org/wiki/Cubic_Hermite_spline, the input x (0-1) should be normalized by grid size
    subroutine Hermitian_base_3order(out_array, x)
        double precision, intent(out)            :: out_array(0:1, 0:1)
        double precision, intent(in)             :: x

        out_array(0, 0) = (1+2*x)*(1-x)**2  ! h00
        out_array(0, 1) = x**2*(3-2*x)      ! h01
        out_array(1, 0) = x*(1-x)**2        ! h10
        out_array(1, 1) = x**2*(x-1)        ! h11
    end subroutine

    subroutine Hermitian_interp_3order_3d(out_value, u, nmd_x, nmd_y, nmd_z, hx, hy, hz)
        double precision, intent(out)               :: out_value
        double precision, intent(in)                :: u(-1:2, -1:2, -1:2) ! the point must lie in the cell (0:1, 0:1, 0:1)
        double precision, intent(in)                :: nmd_x, nmd_y, nmd_z ! distance to the center grid normalized by grid size, must in (0, 1)
        double precision, intent(in)                :: hx, hy, hz
        ! local variables   
        double precision, dimension(0:1, 0:1)       :: weight_x, weight_y, weight_z ! weight in each dimension
        double precision, dimension(0:1, 0:1, 0:1)  :: f_x, f_y, f_z, f_xy, f_xz, f_yz, f_xyz ! derivatives of different order
        integer                                     :: i, j, k

        ! calculate the Hermitian weight in different direction
        call Hermitian_base_3order(weight_x, nmd_x)
        call Hermitian_base_3order(weight_y, nmd_y)
        call Hermitian_base_3order(weight_z, nmd_z)
        do i=0, 1
            do j=0, 1
                do k=0, 1
                    f_x(i, j, k) = derivatives_1st(u(i-1:i+1, j, k), hx)
                    f_y(i, j, k) = derivatives_1st(u(i, j-1:j+1, k), hy)
                    f_z(i, j, k) = derivatives_1st(u(i, j, k-1:k+1), hz)
                    f_xy(i, j, k) = derivatives_mixed_2nd(u(i-1:i+1, j-1:j+1, k), hx, hy)
                    f_xz(i, j, k) = derivatives_mixed_2nd(u(i-1:i+1, j, k-1:k+1), hx, hz)
                    f_yz(i, j, k) = derivatives_mixed_2nd(u(i, j-1:j+1, k-1:k+1), hy, hz)
                    f_xyz(i, j, k) = derivatives_mixed_3rd(u(i-1:i+1, j-1:j+1, k-1:k+1), hx, hy, hz)
                end do
            end do
        end do

        ! calculate the final interpolated value
        out_value = 0.
        do i=0, 1
            do j=0, 1
                do k=0, 1
                    out_value = out_value+weight_x(0, i)*weight_y(0, j)*weight_z(0, k)*u(i, j, k)+&
                                hx*weight_x(1, i)*weight_y(0, j)*weight_z(0, k)*f_x(i, j, k)+&
                                hy*weight_x(0, i)*weight_y(1, j)*weight_z(0, k)*f_y(i, j, k)+&
                                hz*weight_x(0, i)*weight_y(0, j)*weight_z(1, k)*f_z(i, j, k)+&
                                hx*hy*weight_x(1, i)*weight_y(1, j)*weight_z(0, k)*f_xy(i, j, k)+&
                                hx*hz*weight_x(1, i)*weight_y(0, j)*weight_z(1, k)*f_xz(i, j, k)+&
                                hy*hz*weight_x(0, i)*weight_y(1, j)*weight_z(1, k)*f_yz(i, j, k)+&
                                hx*hy*hz*weight_x(1, i)*weight_y(1, j)*weight_z(1, k)*f_xyz(i, j, k)
                end do
            end do
        end do
    end subroutine Hermitian_interp_3order_3d

    subroutine lagrangian_interp_high_order_3d(out_value, u, nmd_x, nmd_y, nmd_z, idx_low, idx_high) ! already assumed uniform grid
        double precision, intent(out)            :: out_value
        double precision, intent(in)             :: u(idx_low:idx_high, idx_low:idx_high, idx_low:idx_high) ! point must in the grid cell (0:1, 0:1, 0:1)
        double precision, intent(in)             :: nmd_x, nmd_y, nmd_z
        integer, intent(in)                      :: idx_low, idx_high
        ! local variables
        integer                                  :: i, j, k, l
        double precision                         :: prod = 1.

        out_value = 0.
        do i=idx_low, idx_high
            do j=idx_low, idx_high
                do k=idx_low, idx_high
                    prod = 1.
                    do l=idx_low, idx_high
                        if (i/=l) prod = prod*((nmd_x-l)/(i-l))
                        if (j/=l) prod = prod*((nmd_y-l)/(j-l))
                        if (k/=l) prod = prod*((nmd_z-l)/(k-l))
                    enddo
                    out_value = out_value+prod*u(i, j, k)
                enddo
            enddo
        enddo

    end subroutine lagrangian_interp_high_order_3d

end module mod_c2b_cart3d

