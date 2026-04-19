module mod_healpix_det
    ! aim: calculate the outflow on the healpix spheres
    ! note that the default radius here are isotropic radius, change this module accordingly if FullGR is implemented in the future
    ! note on implement new outflows: go to cal_flux_single_detector, implement it similiar with Poynting flux and mass outflow
    use HDF5
    implicit none

    public

    {#IFDEF M1
    integer                                       :: where_is_flux(1:8) ! position of the fluxes in all the variables
    character(len=20), dimension(8)               :: allowed_flux_names = ["Poynting", "Mass", "Ub_mass",&
                                                                           "M1_nue","M1_nuebar","M1_nux", &
                                                                           "M1_mu","M1_mubar"] 
                                                                          !["Poynting            ", "Mass                ", "Ub_mass             ",&
                                                                          ! "M1_Flux_nue         ","M1_Flux_nuebar      ","M1_Flux_nux         ", &
                                                                          ! "M1_Flux_mu         ","M1_Flux_mubar       "] ! currently implemented outflows
    }
    {#IFNDEF M1
    integer                                       :: where_is_flux(1:3) ! position of the fluxes in all the variables
    character(len=20), dimension(3)               :: allowed_flux_names = ["Poynting  ", "Mass      ", "Ub_mass   "] ! currently implemented outflows
    }
    integer                                       :: healpix_flux_num = 0
    character(len=1024)                           :: flux_string_for_logfile_header

    ! healpix spherical detectors, use the RING number system, see: https://healpix.sourceforge.io/pdf/intro.pdf
    type healpix_detector
        integer                                       :: current_it ! avoid calculated twice
        integer                                       :: nside ! healpix parameter: nside
        logical                                       :: zsymm ! if true, only use the northern hemisphere
        integer                                       :: pix_end ! largest pix number allowed by the nside
        integer                                       :: north_end ! end index for the northern hemisphere
        integer                                       :: ncount ! number of pixels in this mpirank
        double precision                              :: radius_det   ! radius of the detector
        double precision, allocatable                 :: coord_det^D(:) ! coordinates of detector
        integer, allocatable                          :: indexes(:)  ! healpix index that in this mpi subprocess
        integer, allocatable                          :: which_igrid(:) ! store position of each index
        double precision, allocatable                 :: outflow(:, :) ! store radial outflows


    end type healpix_detector

    type(healpix_detector), allocatable     :: detector_list(:)

	contains

{^IFTHREED

    ! initialize all the detectors
    subroutine initialize_healpix_detectors()
        include 'amrvacdef.f'
        integer                          :: i, var_pos

        if (allocated(detector_list))  call destroy_healpix_detectors()
        if (mype==0) write(*, *) "Initialize outflow detectors!"
        allocate(detector_list(1:healpix_det_num))
        ! find where are the flux parameters and init log file output
        call init_for_log_file_output
        ! initialize each detectors
        do i=1, healpix_det_num
            call initialize_single_detector(healpix_det_radii(i), detector_list(i), i)
        enddo
        if (healpix_det_open) then ! avoid bad parameter setting
            if (healpix_var_num==0) then ! implement more if you want
                call mpistop("Don't set number of healpix variables to be 0!")
            endif
        endif
    end subroutine initialize_healpix_detectors


    ! get healpix_flux_num and the string for header of log file output
    subroutine init_for_log_file_output()
        include 'amrvacdef.f'
        character(len=20)   :: string_j
        integer             :: i, j

        healpix_flux_num = 0
        flux_string_for_logfile_header = ""
        do i=1, size(where_is_flux)
            call find_var_position_in_string(healpix_vname_str, trim(adjustl(allowed_flux_names(i))), healpix_var_num, where_is_flux(i))
            if (where_is_flux(i)/=-1) then
                do j=1, healpix_det_num
                    write(string_j, *) j
                    flux_string_for_logfile_header = trim(adjustl(flux_string_for_logfile_header))&
                        //", flux_"//trim(adjustl(allowed_flux_names(i)))//"_det"//trim(adjustl(string_j))
                end do
                healpix_flux_num = healpix_flux_num+1
                !if (mype==0) write(*, *) where_is_flux(i), allowed_flux_names(i), flux_string_for_logfile_header
            endif
        enddo
    end subroutine init_for_log_file_output


    ! destroy the detectors
    subroutine destroy_healpix_detectors()
        include 'amrvacdef.f'
        ! this happended only in the case of AMR structure changed and the program finish
        integer                         :: i
        do i=1, size(detector_list)
            if (detector_list(i)%ncount/=0) then
                {^D& deallocate(detector_list(i)%coord_det^D) \}
                deallocate(detector_list(i)%indexes)
                deallocate(detector_list(i)%which_igrid)
                deallocate(detector_list(i)%outflow)
            endif
        enddo
        deallocate(detector_list)
        if (mype==0) write(*, *) "Finish destroy outflow detectors!"
    end subroutine destroy_healpix_detectors


    subroutine initialize_single_detector(r_det, detector, det_i)
        include 'amrvacdef.f'
        double precision, intent(in)             :: r_det
        integer, intent(in)                      :: det_i
        type(healpix_detector), intent(inout)    :: detector
        ! local
        double precision                                   :: x^D
        double precision, dimension(1:igridstail, 1:ndim)  :: x_min, x_max
        integer                                            :: count ! total number contained in this mpi rank
        integer                                            :: idx, i, iigrid, igrid
        character(len=10)                                  :: names(1:healpix_var_num)

        count = 0
        detector%nside = healpix_nside(det_i)
        if (healpix_det_zsymm) then
            detector%zsymm = .true.
            detector%pix_end = 6*detector%nside**2+2*detector%nside
        else
            detector%zsymm = .false.
            detector%pix_end = 12*detector%nside**2
        endif
        detector%north_end = 6*detector%nside**2-1 ! end index for the northern hemisphere
        detector%radius_det = r_det
        detector%current_it = -1
        ! find the extremes in coordinate, do not use cc min/max, use block min/max
        do iigrid=1, igridstail
            igrid = igrids(iigrid)
            {^D& x_min(iigrid, ^D) = rnode(rpxmin^D_, igrid)
                x_max(iigrid, ^D) = rnode(rpxmax^D_, igrid) \}
        enddo
        ! count how many num of pixles in this mpi rank
        do idx=0, detector%pix_end-1
            call get_coord_from_healpix_index(detector%nside, r_det, idx, x^D)
            do iigrid=1, igridstail
                !if ({^C& x^D>=x_min(^D, iigrid) .and.} .and. {^C& x^D<x_max(^D, iigrid) .and.})
                   if ((x1>=x_min(iigrid, 1)) .and. (x2>=x_min(iigrid, 2)) .and. (x3>=x_min(iigrid, 3))&
                    .and. (x1<x_max(iigrid, 1)) .and. (x2<x_max(iigrid, 2)) .and. (x3<x_max(iigrid, 3))) then
                    count = count+1
                endif
            enddo
        enddo
        ! allocate the cooresponding array
        detector%ncount = count
        if (count/=0) then
            {^D& allocate(detector%coord_det^D(1:count)) \}
            allocate(detector%indexes(1:count))
            allocate(detector%which_igrid(1:count))
            ! calculate info of these pixels
            i = 0
            do idx=0, detector%pix_end-1
                call get_coord_from_healpix_index(detector%nside, r_det, idx, x^D)
                loopgrid: do iigrid=1, igridstail
                    ! how to use index language bellow ?
                    if ((x1>=x_min(iigrid, 1)) .and. (x2>=x_min(iigrid, 2)) .and. (x3>=x_min(iigrid, 3))&
                        .and. (x1<x_max(iigrid, 1)) .and. (x2<x_max(iigrid, 2)) .and. (x3<x_max(iigrid, 3))) then
                        i = i+1
                        {^D& detector%coord_det^D(i) = x^D \}
                        detector%indexes(i) = idx
                        detector%which_igrid(i) = igrids(iigrid)
                        exit loopgrid
                    endif
                enddo loopgrid
            enddo
            allocate(detector%outflow(1:detector%ncount, 1:healpix_var_num))
        endif
    end subroutine initialize_single_detector


    ! transform index to coordinate, fixme: add surport for 2 dimension case later
    subroutine get_coord_from_healpix_index(nside, radius, idx, x^D)
    ! get the coordinate from the healpix ring index form, see arXiv:astro-ph/0409513
    ! idx should run from 0 to 6N**2+2N-1 (include) for the northern hemisphere+equitor if zsymm
    ! idx should run from 0 to 12N**2-1 (include) if not zsymm
        use mod_cfc_parameters
        use mod_metric
        include 'amrvacdef.f'
        integer, intent(in)                   :: idx, nside
        double precision, intent(in)          :: radius
        double precision, intent(out)         :: x^D
        ! local variables
        integer                :: i, j, ss
        double precision       :: ph, z, phi, theta, aux

        ! calculate the north polar cap
        if (idx<2*(nside-1)*nside) then
            ph = (idx+1.)/2
            aux = int(ph)
            aux = dsqrt(ph-dsqrt(aux))
            i = int(aux)+1
            j = idx+1-2*i*(i-1)
            z = 1-i**2/(3.0d0*nside**2)
            ss = 1
            phi = dpi*(j-ss/2.0)/(2*i)
        ! calculate the north and south belt
        else if (idx<10*nside**2+2*nside) then
            ph = idx-2*nside*(nside-1)
            i = int(ph/(4.*nside))+nside
            j = mod(ph, 4.*nside)+1
            z = 4.0d0/3.0d0-2*i/(3.0d0*nside)
            ss = mod(i-nside+1, 2)
            phi = dpi*(j+ss/2.0d0-1.0d0)/(2*nside)
        ! calculate the south polar cap: just a mirror of the north
        else if (idx<12*nside**2) then
            ph = ((12*nside**2-idx-1)+1.)/2
            aux = int(ph)
            aux = dsqrt(ph-dsqrt(aux))
            i = int(aux)+1
            j = idx+1-2*i*(i-1)
            z = i**2/(3.0d0*nside**2)-1
            ss = 1
            phi = 2*dpi-dpi*(j-ss/2.0)/(2*i)
        else
            call mpistop("The index larger than 12*nside**2 is not allowed in healpix outflow!")
        end if
        ! output coordinate according to coordinate type
        select case(coordinate)
        case(cartesian)
            select case(ndir)
            case(3)
                theta = dacos(z)
                x1 = dsin(theta) * dcos(phi) * radius
                x2 = dsin(theta) * dsin(phi) * radius
                x3 = radius*z
            case default
                call mpistop("2/1 dimensional outflow detector in cartesian coordinate not implemented yet!")
            end select
        case(spherical)
            select case(ndir)
            case(3)
                theta = dacos(z)
                x1 = radius
                x2 = theta
                x3 = phi
            case default
                call mpistop("2/1 dimensional outflow detector in spherical coordinate not implemented yet!")
            end select
        case(cylindrical)
            select case(ndir)
            case(3)
                x1 = radius
                x2 = z
                x3 = phi
            case default
                call mpistop("2/1 dimensional outflow detector in cylindrical coordinate not implemented yet!")
            end select
        case default
            call mpistop("Unknown type of coordinate for the flux detector!")
        end select
    end subroutine get_coord_from_healpix_index


    subroutine cal_flux_single_detector(detector, noutflow)
        ! note that the geometric factors is not included here and must be multiplied when calculating the total flux
        use mod_c2b_cart3d, only: linear_interp_3d
        include 'amrvacdef.f'
        integer, intent(in)                   :: noutflow
        type(healpix_detector), intent(inout) :: detector
        ! local
        integer            :: i_var, i_pixel, igrid, idx_var
        integer            :: nbmin^D, nbmax^D   ! index of the neighbor
        double precision   :: local_dx^D, nmd^D  ! distance to the neighbor
        double precision   :: sqrt_gamma({^D& 1:2 ,}), sqrt_gamma_intp
        character(len=20)  :: names(1:healpix_var_num)
        integer            :: speciesKSP ! speciesKSP of neutrino

        if (detector%ncount==0) return
        call separate_name_string(healpix_vname_str, healpix_var_num, names)
        do i_pixel=1, detector%ncount
            ! Find in the neighbors
            igrid = detector%which_igrid(i_pixel)
            {^D& local_dx^D = (px(igrid)%x(ixGhi^DD, ^D)-px(igrid)%x(ixGlo^DD, ^D))/(ixGhi^D-ixGlo^D) !rnode(rpdx^D_, igrid)
                 nmd^D = MOD(detector%coord_det^D(i_pixel)-px(igrid)%x(ixGlo^DD, ^D), local_dx^D)/local_dx^D
                 nbmin^D = int((detector%coord_det^D(i_pixel)-px(igrid)%x(ixGlo^DD, ^D))/local_dx^D)+ixGlo^D
                 nbmax^D = nbmin^D+1
            \}
            call get_sqrt_gamma_hat(px(igrid)%x(nb^S, :), {^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, sqrt_gamma({^D& 1:2 ,}))
            {#IFDEF DY_SP
                sqrt_gamma({^D& 1:2 ,}) = sqrt_gamma({^D& 1:2 ,})*(pw(igrid)%w(nb^S, psi_metric_))**6
            }
            ! note: interpolation and projection can commute, result will not change much
            call linear_interp_3d(sqrt_gamma_intp, sqrt_gamma({^D& 1:2 ,}), {^D& nmd^D ,})
            do i_var=1, healpix_var_num
                if (names(i_var)=="Poynting") then
                    call calculate_flux_poynting(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="Mass") then
                    call calculate_flux_mass(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="Ub_mass") then
                    call calculate_flux_unbound_mass(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="Spec_j") then
                    call calculate_shell_specific_j(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="Velo_r") then
                    call calculate_shell_velo_r(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="Velo_infty") then
                    call calculate_shell_velo_asymptotic(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var), .true.)
                else if (names(i_var)=="Parker_c") then
                    call calculate_shell_velo_parker_criterion(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="Omega") then
                    call calculate_shell_omega_xy(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="Inv_beta") then
                    call calculate_shell_inverse_beta(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var))
                else if (names(i_var)=="M1_nue") then
                    speciesKSP = 1
                    call calculate_flux_m1(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var), speciesKSP)
                else if (names(i_var)=="M1_nuebar") then
                    if(^NS .gt. 1) then
                      speciesKSP = 2
                      call calculate_flux_m1(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var), speciesKSP)
                    else
                      call mpistop("You don't have enough neutrino species for this!")
                    end if 
                else if (names(i_var)=="M1_nux") then
                    if(^NS .gt. 2) then
                      speciesKSP = 3
                      call calculate_flux_m1(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var), speciesKSP)
                    else
                      call mpistop("You don't have enough neutrino species for this!")
                    end if 
                else if (names(i_var)=="M1_mu") then
                    if(^NS .gt. 3) then
                      speciesKSP = 4
                      call calculate_flux_m1(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var), speciesKSP)
                    else
                      call mpistop("You don't have enough neutrino species for this!")
                    end if
                else if (names(i_var)=="M1_mubar") then
                    if(^NS .gt. 4) then
                      speciesKSP = 5
                      call calculate_flux_m1(nb^L, nmd^D, detector, i_pixel, igrid, detector%outflow(i_pixel, i_var), speciesKSP)
                    else
                      call mpistop("You don't have enough neutrino species for this!")
                    end if                    
                else if ( (INDEX(wnames, TRIM(names(i_var)))/=0) ) then
                    call find_var_position_in_string(wnames, names(i_var), nw+nwauxio, idx_var)
                    call linear_interp_3d(detector%outflow(i_pixel, i_var), pw(igrid)%w(nb^S, idx_var), {^D& nmd^D ,})
                else
                    call mpistop("variable name for calculating the flux not valid: "//names(i_var))
                endif
            enddo
            ! times the sqrt gamma
            detector%outflow(i_pixel, :) = detector%outflow(i_pixel, :)*sqrt_gamma_intp
        end do
    end subroutine cal_flux_single_detector

    subroutine calculate_flux_m1(nb^L, nmd^D, detector, pixel_index, igrid, flux, speciesKSP)  ! d
        use mod_post_process, only: m1_outflow_flux
        {#IFDEF M1
        use mod_m1_outflow
        }
        use mod_c2b_cart3d, only: linear_interp_3d
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: flux
        integer, intent(in)                   :: speciesKSP
        ! local
        double precision                      :: intp_outf(1:ndim), outflux({^D& 1:2 ,}, 1:ndim)

        flux = 0.0d0
        call m1_outflow_flux({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, &
            pw(igrid)%w(nb^S, 1:nw), px(igrid)%x(nb^S, :), &
            outflux({^D& 1:2 ,}, 1:ndim), speciesKSP)
        {^D& call linear_interp_3d(intp_outf(^D), outflux({^D& 1:2 ,}, ^D), {^DD& nmd^DD ,}) \}
        ! projection onto radial direction
        call projection_flux_on_sphere({^D& detector%coord_det^D(pixel_index) ,}, intp_outf(1:ndim), flux)
    end subroutine calculate_flux_m1

    
    subroutine calculate_flux_poynting(nb^L, nmd^D, detector, pixel_index, igrid, flux)  ! d
        use mod_post_process, only: poynting_flux
        use mod_c2b_cart3d, only: linear_interp_3d
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: flux
        ! local
        double precision                      :: intp_outf(1:ndim), poynt({^D& 1:2 ,}, 1:ndim)

        flux = 0.0d0
        call poynting_flux({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, pw(igrid)%w(nb^S, 1:nw), &
            px(igrid)%x(nb^S, :), poynt({^D& 1:2 ,}, 1:ndim))
        {^D& call linear_interp_3d(intp_outf(^D), poynt({^D& 1:2 ,}, ^D), {^DD& nmd^DD ,}) \}
        ! projection onto radial direction
        call projection_flux_on_sphere({^D& detector%coord_det^D(pixel_index) ,}, intp_outf(1:ndim), flux)
    end subroutine calculate_flux_poynting


    subroutine calculate_shell_inverse_beta(nb^L, nmd^D, detector, pixel_index, igrid, shell_value)
        ! Parker instability criterion
        use mod_c2b_cart3d, only: linear_interp_3d
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: shell_value
        ! local
        double precision                      :: inv_beta({^D& 1:2 ,}), b2({^D& 1:2 ,}), P_th({^D& 1:2 ,})

        call imhd_get_intermediate_variables({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, &
            pw(igrid)%w(nb^S, 1:nw), px(igrid)%x(nb^S, :), b2=b2({^D& 1:2 ,}), P_th=P_th({^D& 1:2 ,}))

        shell_value = 0.0d0
        inv_beta = b2/P_th
        call linear_interp_3d(shell_value, inv_beta, {^D& nmd^D ,})
    end subroutine calculate_shell_inverse_beta


    subroutine calculate_shell_omega_xy(nb^L, nmd^D, detector, pixel_index, igrid, shell_value)
        ! Parker instability criterion
        use mod_c2b_cart3d, only: linear_interp_3d
        use mod_post_process, only: cartesian_3d_get_angular_omega
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: shell_value
        ! local
        double precision                      :: my_omega({^D& 1:2 ,})

        shell_value = 0.0d0
        call cartesian_3d_get_angular_omega({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, pw(igrid)%w(nb^S, 1:nw), px(igrid)%x(nb^S, :), my_omega)
        call linear_interp_3d(shell_value, my_omega, {^D& nmd^D ,})
    end subroutine calculate_shell_omega_xy


    subroutine calculate_shell_velo_parker_criterion(nb^L, nmd^D, detector, pixel_index, igrid, shell_value)
        ! Parker instability criterion
        use mod_c2b_cart3d, only: linear_interp_3d
        use mod_post_process, only: cal_parker_criterion
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: shell_value
        ! local
        double precision                      :: parker_c({^D& 1:2 ,})

        shell_value = 0.0d0
        call cal_parker_criterion({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, pw(igrid)%w(nb^S, 1:nw), px(igrid)%x(nb^S, :), parker_c)
        call linear_interp_3d(shell_value, parker_c, {^D& nmd^D ,})
    end subroutine calculate_shell_velo_parker_criterion


    subroutine calculate_shell_velo_asymptotic(nb^L, nmd^D, detector, pixel_index, igrid, shell_value, include_mag)
        ! calculate the asymptotic velocity in infinity by assuming W_inf = -hu_t
        ! h depend on whether to include mag contribution or not (include_mag)
        use mod_c2b_cart3d, only: linear_interp_3d
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: shell_value
        logical, intent(in)                   :: include_mag
        ! local
        double precision                      :: htot_tmp({^D& 1:2 ,}), u_t({^D& 1:2 ,}), velo_tmp({^D& 1:2 ,}), lfac_tmp({^D& 1:2 ,})

        shell_value = 0.0d0
        if (include_mag) then
            call imhd_get_intermediate_variables({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, pw(igrid)%w(nb^S, 1:nw), &
                px(igrid)%x(nb^S, :), htot=htot_tmp({^D& 1:2 ,}), lfac=lfac_tmp({^D& 1:2 ,}))
        else 
            call imhd_get_intermediate_variables({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, pw(igrid)%w(nb^S, 1:nw), &
                px(igrid)%x(nb^S, :), h_th=htot_tmp({^D& 1:2 ,}), lfac=lfac_tmp({^D& 1:2 ,}))
        endif
        u_t({^D& 1:2 ,}) = -lfac_tmp({^D& 1:2 ,})*pw(igrid)%w(nb^S, alp_metric_)
        {^C& u_t({^D& 1:2 ,}) = u_t({^D& 1:2 ,})+pw(igrid)%w(nb^S, psi_metric_)**4*pw(igrid)%w(nb^S, beta_metric^C_)*pw(igrid)%w(nb^S, u^C_) \}
        velo_tmp = dsqrt(1.-1./(htot_tmp*u_t)**2)
        call linear_interp_3d(shell_value, velo_tmp({^D& 1:2 ,}), {^D& nmd^D ,})
    end subroutine calculate_shell_velo_asymptotic


    subroutine calculate_shell_velo_r(nb^L, nmd^D, detector, pixel_index, igrid, shell_value)
        ! calculate the Eulerian velocity in the R direction
        use mod_c2b_cart3d, only: linear_interp_3d
        use mod_imhd_intermediate, only: dysp_get_lfac
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: shell_value
        ! local
        double precision                      :: intp_outf(1:ndim), aux_velo({^D& 1:2 ,}, 1:ndim), lfac_tmp({^D& 1:2 ,})

        shell_value = 0.0d0
        call dysp_get_lfac({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, pw(igrid)%w(nb^S, 1:nw), &
            px(igrid)%x(nb^S, :), lfac_tmp({^D& 1:2 ,}))
        {^D&
            aux_velo({^D& 1:2 ,}, ^D) = pw(igrid)%w(nb^S, u^D_)/lfac_tmp({^D& 1:2 ,})
            call linear_interp_3d(intp_outf(^D), aux_velo({^D& 1:2 ,}, ^D), {^DD& nmd^DD ,}) 
        \}
        ! projection onto radial direction
        call projection_flux_on_sphere({^D& detector%coord_det^D(pixel_index) ,}, intp_outf(1:ndim), shell_value)
    end subroutine calculate_shell_velo_r


    subroutine calculate_shell_specific_j(nb^L, nmd^D, detector, pixel_index, igrid, shell_value)  ! detector%outflow(j, outf_count)
        ! calculate the specific angular momentum on the sphere
        use mod_c2b_cart3d, only: linear_interp_3d
        use mod_post_process, only: cartesian_3d_specific_angular_momentum
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: shell_value
        ! local
        double precision                      :: spec_j({^D& 1:2 ,})

        shell_value = 0.0d0
        call cartesian_3d_specific_angular_momentum({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, pw(igrid)%w(nb^S, 1:nw), &
            px(igrid)%x(nb^S, :), spec_j({^D& 1:2 ,}), 0.0d0, 0.0d0)
            call linear_interp_3d(shell_value, spec_j, {^D& nmd^D ,})
    end subroutine calculate_shell_specific_j


    subroutine calculate_flux_mass(nb^L, nmd^D, detector, pixel_index, igrid, flux)  ! detector%outflow(j, outf_count)
        ! see equation (13) of Kelly17, 10.1103/PhysRevD.96.123003
        use mod_c2b_cart3d, only: linear_interp_3d
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: flux
        ! local
        double precision                      :: intp_outf(1:ndim)

        flux = 0.0d0
        {#IFDEF DY_SP
            !note: rho*alpha*(WU)-D*beta = D*(alpha*U-beta)
            {^D& call linear_interp_3d(intp_outf(^D), &
                pw(igrid)%w(nb^S, rho_)*pw(igrid)%w(nb^S, alp_metric_)*pw(igrid)%w(nb^S, u^D_)-&
                pw(igrid)%w(nb^S, D_)*pw(igrid)%w(nb^S, beta_metric^D_), {^DD& nmd^DD ,}) \}
        }
        {#IFNDEF DY_SP
            {^D& call linear_interp_3d(intp_outf(^D), &
                pw(igrid)%w(nb^S, rho_)*pw(igrid)%w(nb^S, u^D_), {^DD& nmd^DD ,}) \}
        }
        ! projection onto radial direction
        call projection_flux_on_sphere({^D& detector%coord_det^D(pixel_index) ,}, intp_outf(1:ndim), flux)
    end subroutine calculate_flux_mass


    subroutine calculate_flux_unbound_mass(nb^L, nmd^D, detector, pixel_index, igrid, flux)
        ! using Bernnouli Criteria for unbound condition: -h_tot*u_t-eoshmin
        use mod_c2b_cart3d, only: linear_interp_3d
        use mod_imhd_intermediate, only: imhd_get_intermediate_variables
        include 'amrvacdef.f'
        integer, intent(in)                   :: nbmin^D, nbmax^D   ! index for the neighbors of this pixel
        double precision, intent(in)          :: nmd^D
        type(healpix_detector), intent(inout) :: detector
        integer, intent(in)                   :: pixel_index, igrid
        double precision, intent(out)         :: flux
        ! local
        double precision                      :: intp_outf(1:ndim), ub_condition, ub_cond_srd({^D& 1:2 ,})

        flux = 0.0d0
        call imhd_get_intermediate_variables({^D& 1 ,}, {^D& 2 ,}, {^D& 1 ,}, {^D& 2 ,}, &
            pw(igrid)%w(nb^S, 1:nw), px(igrid)%x(nb^S, :), ub_bernoulli=ub_cond_srd({^D& 1:2 ,}))
        call linear_interp_3d(ub_condition, ub_cond_srd({^D& 1:2 ,}), {^D& nmd^D ,})
        if (ub_condition>0) then
            {#IFDEF DY_SP
                !note: rho*alpha*(WU)-D*beta = D*(alpha*U-beta)
                {^D& call linear_interp_3d(intp_outf(^D), &
                    pw(igrid)%w(nb^S, rho_)*pw(igrid)%w(nb^S, alp_metric_)*pw(igrid)%w(nb^S, u^D_)-&
                    pw(igrid)%w(nb^S, D_)*pw(igrid)%w(nb^S, beta_metric^D_), {^DD& nmd^DD ,}) \}
            }
            {#IFNDEF DY_SP
                {^D& call linear_interp_3d(intp_outf(^D), &
                    pw(igrid)%w(nb^S, rho_)*pw(igrid)%w(nb^S, u^D_), {^DD& nmd^DD ,}) \}
            }
            ! projection onto radial direction
            call projection_flux_on_sphere({^D& detector%coord_det^D(pixel_index) ,}, intp_outf(1:ndim), flux)
        endif
    end subroutine calculate_flux_unbound_mass


    ! get the radial direction of flux
    ! fixme: for more generic metric other than CFC, the radial direction dS_r may has contributions
    ! from other directions in flat spacetime, please extend it accordingly
    subroutine projection_flux_on_sphere(x^D, flux_j, j_out)
        use mod_cfc_parameters
        include 'amrvacdef.f'
        double precision, intent(in)                    :: x^D
        double precision, intent(in)                    :: flux_j(1:ndim)
        double precision, intent(out)                   :: j_out
        ! local
        integer                                         :: i
        double precision                                :: radius

        select case(coordinate)
        case(cartesian)
            radius = dsqrt({^D& x^D**2 +})
            if (radius>10*smalldouble) then
                j_out = {^D& flux_j(^D)*x^D/radius +}
            else
                j_out = 0.0d0
            endif
        case(spherical)
            j_out = flux_j(1)
        case(cylindrical) ! in the order of: r, phi, z
            radius = x1**2+x3**2
            if (radius>10*smalldouble) then
                j_out = flux_j(1)*x1/radius+flux_j(3)*x3/radius
            else
                j_out = 0.0d0
            endif
        case default
            call mpistop("Unknown type of coordinate for the flux detector!")
        end select
    end subroutine projection_flux_on_sphere


    ! surface integral, this must act on a global scale, called by printlog_special and sum by mpi_reduce
    ! please note all the results are integrated in the whole sphere, even if you use z-symmetry
    ! because the flux of non-zsymm is not just simply 2 times of the zsymm one
    subroutine surface_integral_outflow(it_local, ndetector, nflow_var, outf)
        include 'amrvacdef.f'
        integer, intent(in)               :: ndetector, nflow_var, it_local
        double precision, intent(out)     :: outf(1:nflow_var, 1:ndetector)
        ! local
        integer                           :: i, j, k
        double precision                  :: geofac

        call calculate_all_outflow(it_local)
        do i=1, ndetector
            if (detector_list(i)%ncount/=0) then
                geofac = 4*dpi*(detector_list(i)%radius_det)**2/(12*(detector_list(i)%nside)**2)
                k = 0
                do j=1, size(where_is_flux)
                    if (where_is_flux(j)/=-1) then
                        k = k+1
                        if (k<=nflow_var) then
                            if (detector_list(i)%zsymm) then
                                ! north*2+equatorial
                                outf(k, i) = (sum(detector_list(i)%outflow(:, where_is_flux(j)), mask= & 
                                    (detector_list(i)%indexes<=detector_list(i)%north_end))*2+ &
                                    sum(detector_list(i)%outflow(:, where_is_flux(j)), mask= &
                                    (detector_list(i)%indexes>detector_list(i)%north_end)))*geofac
                            else
                                outf(k, i) = sum(detector_list(i)%outflow(:, where_is_flux(j)))*geofac
                            endif
                        else
                            call mpistop("Internal number of fluxes not consistent, please check settings!")
                        endif
                    endif
                enddo
            else
                outf(:, i) = 0.0d0
            endif
        enddo
    end subroutine surface_integral_outflow


    ! calculate all kinds of outflow allowed
    subroutine calculate_all_outflow(it_local)
        include 'amrvacdef.f'
        integer, intent(in)               :: it_local
        integer                           :: i
        ! loop for detectors, count the flux need and allocate the flux array
        if (healpix_var_num/=0) then
            do i=1, healpix_det_num
                ! avoid calculate multiple times
                if (detector_list(i)%current_it/=it_local) then
                    call cal_flux_single_detector(detector_list(i), healpix_var_num)
                    detector_list(i)%current_it = it_local
                endif
            enddo
        else
            call mpistop("Detector specified but the name of outflow variables not allowed or not implemented!")
        endif
    end subroutine calculate_all_outflow


    ! save a snapshot for all the detectors to hdf5 file
    subroutine save_outflow_snap_to_hdf5(it_local)
        ! each group: detector, each dataset: it, each sub dataset: idx, flux1, flux2, ....
        include 'amrvacdef.f'
        integer, intent(in)               :: it_local
        ! local
        integer                           :: i, j, nvarp1, ncount, ncount_recv, rank, status, ipe, pos
        character(len=256)                :: group_name, dset_name, filename, var_names
        double precision, allocatable     :: da_store(:, :), idx_send(:)
        integer                           :: intstatus(MPI_STATUS_SIZE)
        integer                           :: store_dims(1:2), start_dset(1:2) ! shape of array
        integer(HID_T)                    :: file_id, group_id, attr_space_id, attr_id, hdf_dim=1
        logical                           :: lk_exist, at_exist

        ! get the output file name
        write(filename, "(a,a)") TRIM(filenameout), "_outflow.hdf5"
        nvarp1 = healpix_var_num+1
        if (npe>1) call MPI_BARRIER(icomm, ierrmpi)
        if (mype==0) then
            ! Open or create the HDF5 file
            call H5Open_f(status)
            call H5FOpen_f(filename, H5F_ACC_RDWR_F, file_id, status)
            if (status /= 0) then
                call H5FClose_f(file_id, status)
                call h5FCreate_f(filename, H5F_ACC_TRUNC_F, file_id, status)
            end if        
        endif
        loopdet: do i=1, healpix_det_num
            pos = 1
            write(group_name, "(I)") i
            group_name = "det:"//trim(adjustl(group_name))
            write(dset_name, "(I)") it_local
            dset_name = "it:"//trim(adjustl(dset_name))
            ncount = detector_list(i)%ncount
            if (mype/=0) then ! send data if needed
                call MPI_SEND(ncount, 1, MPI_INTEGER, 0, 3*i, icomm, ierrmpi) ! send shape first
                if (ncount/=0) then
                    if (allocated(idx_send)) deallocate(idx_send)
                    allocate(idx_send(1:ncount))
                    do j=1, ncount
                        idx_send(j) = detector_list(i)%indexes(j)
                    enddo
                    call MPI_SEND(idx_send, ncount, &
                        MPI_DOUBLE_PRECISION, 0, 3*i+1, icomm, ierrmpi)
                    call MPI_SEND(detector_list(i)%outflow, healpix_var_num*ncount, &
                        MPI_DOUBLE_PRECISION, 0, 3*i+2, icomm, ierrmpi)
                endif
            else
                if (allocated(da_store)) deallocate(da_store)
                allocate(da_store(1:detector_list(i)%pix_end, 1:nvarp1))
                ! open group
                call H5Lexists_f(file_id, group_name, lk_exist, status)
                if (lk_exist) then
                    call H5GOpen_f(file_id, group_name, group_id, status)
                else
                    call H5GCreate_f(file_id, group_name, group_id, status)
                end if
                ! add attribute radius
                call H5AExists_by_name_f(group_id, group_name, "radius", at_exist, status)
                call h5AOpen_name_f(group_id, "radius", attr_id, status)
                if (status==0) then
                    call H5AWrite_f(attr_id, H5T_NATIVE_DOUBLE, [detector_list(i)%radius_det], [hdf_dim], status)
                    call h5AClose_f(attr_id, status)
                else
                    call H5SCreate_simple_f(H5S_SCALAR_F, [hdf_dim], attr_space_id, status)
                    call H5ACreate_f(group_id, "radius", H5T_NATIVE_DOUBLE, attr_space_id, attr_id, status)
                    call H5AWrite_f(attr_id, H5T_NATIVE_DOUBLE, [detector_list(i)%radius_det], [hdf_dim], status)
                    call H5AClose_f(attr_id, status)
                    call h5SClose_f(attr_space_id, status)
                endif
                ! add attribute nside
                call H5AExists_by_name_f(group_id, group_name, "nside", at_exist, status)
                call h5AOpen_name_f(group_id, "nside", attr_id, status)
                if (status==0) then
                    call H5AWrite_f(attr_id, H5T_NATIVE_INTEGER, [detector_list(i)%nside], [hdf_dim], status)
                    call h5AClose_f(attr_id, status)
                else
                    call H5SCreate_simple_f(H5S_SCALAR_F, [hdf_dim], attr_space_id, status)
                    call H5ACreate_f(group_id, "nside", H5T_NATIVE_INTEGER, attr_space_id, attr_id, status)
                    call H5AWrite_f(attr_id, H5T_NATIVE_INTEGER, [detector_list(i)%nside], [hdf_dim], status)
                    call H5AClose_f(attr_id, status)
                    call h5SClose_f(attr_space_id, status)
                endif
                ! create dataset and save flux for mype=0
                if (ncount/=0) then
                    ! allocate array to store
                    da_store(1:ncount, 1) = detector_list(i)%indexes(1:ncount)
                    da_store(1:ncount, 2:nvarp1) = detector_list(i)%outflow(1:ncount, 1:healpix_var_num)
                    pos = pos+ncount
                endif
            endif
            if (mype==0) then ! collect outflow from other npe 
                if(npe>1) then
                    do ipe=1, npe-1
                        call MPI_RECV(ncount_recv, 1, MPI_INTEGER, ipe, 3*i, icomm, intstatus, ierrmpi)
                        if (ncount_recv/=0) then
                            call MPI_RECV(da_store(pos:pos+ncount_recv-1, 1), ncount_recv, MPI_DOUBLE_PRECISION, &
                                ipe, 3*i+1, icomm, intstatus, ierrmpi)
                            call MPI_RECV(da_store(pos:pos+ncount_recv-1, 2:nvarp1), healpix_var_num*ncount_recv, &
                                MPI_DOUBLE_PRECISION, ipe, 3*i+2, icomm, intstatus, ierrmpi)
                            pos = pos+ncount_recv
                        endif
                    enddo
                endif
                !output to hdf5 file
                call write_single_dataset(group_id, dset_name, "HP_idx "//TRIM(adjustl(healpix_vname_str)), &
                     nvarp1, detector_list(i)%pix_end, da_store)
                call H5GClose_f(group_id, status)
            endif
            if (npe>1) then
                call MPI_BARRIER(icomm, ierrmpi)
            endif
        enddo loopdet
        if (mype==0) then
            call h5FClose_f(file_id, status)
            call H5Close_f(status)
            if (allocated(da_store)) deallocate(da_store)
        else
            if (allocated(idx_send)) deallocate(idx_send)
        endif
    end subroutine save_outflow_snap_to_hdf5


    subroutine write_single_dataset(group_id, dset_name, var_names, total_shape1, total_shape2, da_out)
        include 'amrvacdef.f'
        integer(HID_T), intent(in)      :: group_id
        character(*), intent(in)        :: dset_name, var_names
        integer, intent(in)             :: total_shape1, total_shape2
        double precision, intent(in)    :: da_out(1:total_shape2, 1:total_shape1)
        ! local
        integer(HID_T)                  :: dataset_id, dataspace_id, type_id, attr_space_id, attr_id, hdf_dim=1
        integer(HSIZE_T)                :: init_shape(1:2), len_char
        integer                         :: status
        logical                         :: lk_exist, at_exist

        init_shape(1) = total_shape1
        init_shape(2) = total_shape2
        call H5LExists_f(group_id, dset_name, lk_exist, status)
        if (lk_exist) then
            call H5LDelete_f(group_id, dset_name, status)
        endif
        call H5SCreate_simple_f(2, init_shape, dataspace_id, status)
        call H5DCreate_f(group_id, dset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, status)
        call H5DWrite_f(dataset_id, H5T_NATIVE_DOUBLE, TRANSPOSE(da_out), init_shape, status)
        ! add attribute time
        call h5AOpen_name_f(dataset_id, "time", attr_id, status)
        if (status==0) then
            call H5AWrite_f(attr_id, H5T_NATIVE_DOUBLE, [t], [hdf_dim], status)
            call h5AClose_f(attr_id, status)
        else
            call H5SCreate_simple_f(H5S_SCALAR_F, [hdf_dim], attr_space_id, status)
            call H5ACreate_f(dataset_id, "time", H5T_NATIVE_DOUBLE, attr_space_id, attr_id, status)
            call H5AWrite_f(attr_id, H5T_NATIVE_DOUBLE, [t], [hdf_dim], status)
            call H5AClose_f(attr_id, status)
            call h5SClose_f(attr_space_id, status)
        endif
        ! add atrribute name of all variables
        call h5AOpen_name_f(dataset_id, "var_names", attr_id, status)
        if (status==0) then
            call H5AWrite_f(attr_id, H5T_STRING, [var_names], [hdf_dim], status)
            call h5AClose_f(attr_id, status)
        else
            ! Define a variable-length string datatype
            call H5Tcopy_f(H5T_C_S1, type_id, status)       ! Start with the default string datatype
            len_char = len(var_names)
            call H5Tset_size_f(type_id, len_char, status)
            call H5SCreate_simple_f(H5S_SCALAR_F, [hdf_dim], attr_space_id, status)
            call H5ACreate_f(dataset_id, "var_names", type_id, attr_space_id, attr_id, status)
            call H5AWrite_f(attr_id, type_id, [var_names], [hdf_dim], status)
            call H5AClose_f(attr_id, status)
            call h5SClose_f(attr_space_id, status)
        endif
        call H5DClose_f(dataset_id, status)
        call H5SClose_f(dataspace_id, status)
    end subroutine write_single_dataset
}

! auxilary subroutines
    subroutine separate_name_string(input_string, num_str, var_names)
        character(*), intent(in)         :: input_string
        integer, intent(in)              :: num_str
        character(*), intent(out)        :: var_names(1:num_str)
        ! Local variables
        integer                          :: iw, pos1, pos2
        character(len=len(input_string)) :: scanstring  
    
        scanstring = TRIM(adjustl(input_string))  ! Normalize input string
        pos1 = 1
        iw = 0
        ! Initialize output var_names array to empty strings
        var_names = ''
        do
            ! Skip leading spaces
            do while (pos1 <= len(scanstring) .and. scanstring(pos1:pos1) == " ")
                pos1 = pos1 + 1
            enddo
            if (pos1 > len(scanstring)) exit  ! End of string
            ! Find the position of the next space or end of the string
            pos2 = INDEX(scanstring(pos1:), " ")
            if (pos2 == 0) then  ! Last word in the string
                iw = iw + 1
                if (iw <= num_str) var_names(iw) = scanstring(pos1:)
                exit
            endif
            ! Extract the word and store it
            iw = iw + 1
            if (iw <= num_str) var_names(iw) = scanstring(pos1:pos1+pos2-2)
            ! Move to the next word
            pos1 = pos1 + pos2
        enddo
        ! Fill remaining slots with empty strings if input had fewer var_names
        if (iw < num_str) var_names(iw+1:num_str) = ''
    end subroutine separate_name_string
    
    subroutine find_var_position_in_string(input_string, name_var, num_str, var_pos)
        character(*), intent(in)         :: input_string
        character(*), intent(in)         :: name_var
        integer, intent(in)              :: num_str
        integer, intent(out)             :: var_pos
        ! Local variables
        integer                          :: pos1, pos2
        integer                          :: word_count
        character(len=len(input_string)) :: scanstring
    
        scanstring = TRIM(adjustl(input_string))  ! Normalize input string
        pos1 = 1
        word_count = 0
        var_pos = -1  ! Initialize as not found
        do
            ! Skip leading spaces before the current word
            do while (pos1 <= len(scanstring) .and. scanstring(pos1:pos1) == " ")
                pos1 = pos1 + 1
            enddo
            if (pos1 > len(scanstring)) exit  ! End of string
            ! Find the position of the next space or the end of the string
            pos2 = INDEX(scanstring(pos1:), " ")
            if (pos2 == 0) then  ! Last word in the string
                word_count = word_count + 1
                if (TRIM(scanstring(pos1:)) == TRIM(name_var)) then
                    var_pos = word_count
                    return
                endif
                exit
            endif
            word_count = word_count + 1
            if (TRIM(scanstring(pos1:pos1+pos2-2)) == TRIM(name_var)) then
                var_pos = word_count
                return
            endif
            pos1 = pos1 + pos2
        enddo
        !if (var_pos==-1) write(*, *) "Warning: can't find string in subroutine find_var_position_in_string."
    end subroutine find_var_position_in_string

end module mod_healpix_det
