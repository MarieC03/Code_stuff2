!=============================================================================
!INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
!INCLUDE:amrvacnul/usrflags.t
!INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================

!========================================================================================
!==========                   module to be included in NSimport                ==========
!========================================================================================


module Init_c2b_cart_3d_interface
    use mod_c2b_cart3d
    use mod_indices
    !include 'amrvacdef.f'
    public
    
    type c2b_3d_profile
        type(carpet_variable)                         :: rho
        type(carpet_variable)                         :: velox, veloy, veloz
        type(carpet_variable)                         :: alpha, psi, betax, betay, betaz
        type(carpet_variable)                         :: temp, ye, ymu, eps, pressure
        type(carpet_variable)                         :: Ax, Ay, Az, Bx, By, Bz
        type(carpet_variable)                         :: Enue, Enue_bar, Enux, Enumu, Enumu_bar
        type(carpet_variable)                         :: Nnue, Nnue_bar, Nnux, Nnumu, Nnumu_bar
        type(carpet_variable)                         :: Fnue_x, Fnue_bar_x, Fnux_x, Fnumu_x, Fnumu_bar_x
        type(carpet_variable)                         :: Fnue_y, Fnue_bar_y, Fnux_y, Fnumu_y, Fnumu_bar_y
        type(carpet_variable)                         :: Fnue_z, Fnue_bar_z, Fnux_z, Fnumu_z, Fnumu_bar_z
    end type c2b_3d_profile

    type(c2b_3d_profile)                              :: my_c2b_3d_profile


    contains

    subroutine init_c2b_3d_profile(fpath, fnbase_list, dsname_list, num_refine_levels, num_datablocks, iter_num, &
            save_cpu_or_mem, init_magnetic_field, init_B_from_vecA, use_eps_init, init_m1_vars)
        use mod_eos
        character(len=256), intent(in)                :: fpath ! path to directory where hdf5 files are stored
        character(len=256), dimension(45), intent(in) :: fnbase_list ! name of the file that store the variable
        character(len=256), dimension(45), intent(in) :: dsname_list ! name of the variable in the attribute

        !character(len=256), dimension(29), intent(in) :: fnbase_list ! name of the file that store the variable
        !character(len=256), dimension(29), intent(in) :: dsname_list ! name of the variable in the attribute

        !character(len=256), dimension(24), intent(in) :: fnbase_list ! name of the file that store the variable
        !character(len=256), dimension(24), intent(in) :: dsname_list ! name of the variable in the attribute

        integer, intent(in)                           :: num_refine_levels, num_datablocks, iter_num, save_cpu_or_mem
        logical, intent(in)                           :: init_magnetic_field, init_B_from_vecA, use_eps_init, init_m1_vars

        integer, parameter                            :: FIL_rho       = 1 
        integer, parameter                            :: FIL_velocx    = 2 
        integer, parameter                            :: FIL_velocy    = 3 
        integer, parameter                            :: FIL_velocz    = 4 
        integer, parameter                            :: FIL_psi       = 5
        integer, parameter                            :: FIL_alpha     = 6
        integer, parameter                            :: FIL_betax     = 7
        integer, parameter                            :: FIL_betay     = 8
        integer, parameter                            :: FIL_betaz     = 9
        integer, parameter                            :: FIL_temp      = 10
        integer, parameter                            :: FIL_ye        = 11
        integer, parameter                            :: FIL_ymu       = 12
        integer, parameter                            :: FIL_eps       = 13
        integer, parameter                            :: FIL_pressure  = 14
        integer, parameter                            :: FIL_Ax        = 15
        integer, parameter                            :: FIL_Ay        = 16
        integer, parameter                            :: FIL_Az        = 17
        integer, parameter                            :: FIL_Bx        = 18
        integer, parameter                            :: FIL_By        = 19
        integer, parameter                            :: FIL_Bz        = 20

        integer, parameter                            :: FIL_Nnue      = 21
        integer, parameter                            :: FIL_Enue      = 22
        integer, parameter                            :: FIL_Fnue_x    = 23
        integer, parameter                            :: FIL_Fnue_y    = 24
        integer, parameter                            :: FIL_Fnue_z    = 25

        integer, parameter                            :: FIL_Nnue_bar  = 26
        integer, parameter                            :: FIL_Enue_bar  = 27
        integer, parameter                            :: FIL_Fnue_bar_x = 28
        integer, parameter                            :: FIL_Fnue_bar_y = 29
        integer, parameter                            :: FIL_Fnue_bar_z = 30
        
        integer, parameter                            :: FIL_Nnux      = 31
        integer, parameter                            :: FIL_Enux      = 32
        integer, parameter                            :: FIL_Fnux_x    = 33
        integer, parameter                            :: FIL_Fnux_y    = 34
        integer, parameter                            :: FIL_Fnux_z    = 35
        
        integer, parameter                            :: FIL_Nnumu     = 36
        integer, parameter                            :: FIL_Enumu     = 37
        integer, parameter                            :: FIL_Fnumu_x   = 38
        integer, parameter                            :: FIL_Fnumu_y   = 39
        integer, parameter                            :: FIL_Fnumu_z   = 40
        
        integer, parameter                            :: FIL_Nnumu_bar  = 41
        integer, parameter                            :: FIL_Enumu_bar  = 42
        integer, parameter                            :: FIL_Fnumu_bar_x = 43
        integer, parameter                            :: FIL_Fnumu_bar_y = 44
        integer, parameter                            :: FIL_Fnumu_bar_z = 45

        integer :: Nspecies

        if (mype==0) write(*, *) "Initializing c2b_3d_profile"
 
        call init_single_variable(my_c2b_3d_profile%rho, fpath, fnbase_list(FIL_rho), dsname_list(FIL_rho), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%velox, fpath, fnbase_list(FIL_velocx), dsname_list(FIL_velocx), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%veloy, fpath, fnbase_list(FIL_velocy), dsname_list(FIL_velocy), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%veloz, fpath, fnbase_list(FIL_velocz), dsname_list(FIL_velocz), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%psi, fpath, fnbase_list(FIL_psi), dsname_list(FIL_psi), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%alpha, fpath, fnbase_list(FIL_alpha), dsname_list(FIL_alpha), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%betax, fpath, fnbase_list(FIL_betax), dsname_list(FIL_betax), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%betay, fpath, fnbase_list(FIL_betay), dsname_list(FIL_betay), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
        call init_single_variable(my_c2b_3d_profile%betaz, fpath, fnbase_list(FIL_betaz), dsname_list(FIL_betaz), &
            num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)

        if (init_magnetic_field) then
            if (init_B_from_vecA) then
                call init_single_variable(my_c2b_3d_profile%Ax, fpath, fnbase_list(FIL_Ax), dsname_list(FIL_Ax), &
                    num_refine_levels, num_datablocks, iter_num, 'FTT', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Ay, fpath, fnbase_list(FIL_Ay), dsname_list(FIL_Ay), &
                    num_refine_levels, num_datablocks, iter_num, 'TFT', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Az, fpath, fnbase_list(FIL_Az), dsname_list(FIL_Az), &
                    num_refine_levels, num_datablocks, iter_num, 'TTF', save_cpu_or_mem)
            else
                call init_single_variable(my_c2b_3d_profile%Bx, fpath, fnbase_list(FIL_Bx), dsname_list(FIL_Bx), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%By, fpath, fnbase_list(FIL_By), dsname_list(FIL_By), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Bz, fpath, fnbase_list(FIL_Bz), dsname_list(FIL_Bz), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            endif
        end if

        if (eos_uses_ye()) then
            call init_single_variable(my_c2b_3d_profile%ye, fpath, fnbase_list(FIL_ye), dsname_list(FIL_ye), &
                num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            if (eos_has_ymu()) then
                call init_single_variable(my_c2b_3d_profile%ymu, fpath, fnbase_list(FIL_ymu), dsname_list(FIL_ymu), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            endif
            if (use_eps_init) then
                call init_single_variable(my_c2b_3d_profile%eps, fpath, fnbase_list(FIL_eps), dsname_list(FIL_eps), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            else
                call init_single_variable(my_c2b_3d_profile%temp, fpath, fnbase_list(FIL_temp), dsname_list(FIL_temp), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            endif
        else
            if (use_eps_init) then
                call init_single_variable(my_c2b_3d_profile%eps, fpath, fnbase_list(FIL_eps), dsname_list(FIL_eps), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            else
                call init_single_variable(my_c2b_3d_profile%pressure, fpath, fnbase_list(FIL_pressure), dsname_list(FIL_pressure), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            endif
        endif

        if(init_m1_vars)then
            Nspecies = ^NS
          ! write(*,*)"init_c2b_3D: Number of species : NS",^NS
            ! species 1
            if(Nspecies .ge. 1) then
             !write(*,*)"Nspecies >= 1",Nspecies
                call init_single_variable(my_c2b_3d_profile%Nnue, fpath, fnbase_list(FIL_Nnue), dsname_list(FIL_Nnue), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Enue, fpath, fnbase_list(FIL_Enue), dsname_list(FIL_Enue), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnue_x, fpath, fnbase_list(FIL_Fnue_x), dsname_list(FIL_Fnue_x), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnue_y, fpath, fnbase_list(FIL_Fnue_y), dsname_list(FIL_Fnue_y), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnue_z, fpath, fnbase_list(FIL_Fnue_z), dsname_list(FIL_Fnue_z), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            end if 
            ! species 2
            if(Nspecies .ge. 2) then
             !write(*,*)"Nspecies >= 2",Nspecies
                call init_single_variable(my_c2b_3d_profile%Nnue_bar, fpath, fnbase_list(FIL_Nnue_bar), dsname_list(FIL_Nnue_bar), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Enue_bar, fpath, fnbase_list(FIL_Enue_bar), dsname_list(FIL_Enue_bar), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnue_bar_x, fpath, fnbase_list(FIL_Fnue_bar_x), dsname_list(FIL_Fnue_bar_x), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnue_bar_y, fpath, fnbase_list(FIL_Fnue_bar_y), dsname_list(FIL_Fnue_bar_y), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnue_bar_z, fpath, fnbase_list(FIL_Fnue_bar_z), dsname_list(FIL_Fnue_bar_z), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            end if
            ! species 3
            if(Nspecies .ge. 3) then
             !write(*,*)"Nspecies >= 3",Nspecies
                call init_single_variable(my_c2b_3d_profile%Nnux, fpath, fnbase_list(FIL_Nnux), dsname_list(FIL_Nnux), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Enux, fpath, fnbase_list(FIL_Enux), dsname_list(FIL_Enux), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnux_x, fpath, fnbase_list(FIL_Fnux_x), dsname_list(FIL_Fnux_x), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnux_y, fpath, fnbase_list(FIL_Fnux_y), dsname_list(FIL_Fnux_y), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnux_z, fpath, fnbase_list(FIL_Fnux_z), dsname_list(FIL_Fnux_z), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            end if 
            if(Nspecies .ge. 4) then
                call init_single_variable(my_c2b_3d_profile%Nnumu, fpath, fnbase_list(FIL_Nnumu), dsname_list(FIL_Nnumu), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Enumu, fpath, fnbase_list(FIL_Enumu), dsname_list(FIL_Enumu), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnumu_x, fpath, fnbase_list(FIL_Fnumu_x), dsname_list(FIL_Fnumu_x), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnumu_y, fpath, fnbase_list(FIL_Fnumu_y), dsname_list(FIL_Fnumu_y), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnumu_z, fpath, fnbase_list(FIL_Fnumu_z), dsname_list(FIL_Fnumu_z), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            end if
            if(Nspecies .ge. 5) then
                call init_single_variable(my_c2b_3d_profile%Nnumu_bar, fpath, fnbase_list(FIL_Nnumu_bar), dsname_list(FIL_Nnumu_bar), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Enumu_bar, fpath, fnbase_list(FIL_Enumu_bar), dsname_list(FIL_Enumu_bar), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnumu_bar_x, fpath, fnbase_list(FIL_Fnumu_bar_x), dsname_list(FIL_Fnumu_bar_x), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnumu_bar_y, fpath, fnbase_list(FIL_Fnumu_bar_y), dsname_list(FIL_Fnumu_bar_y), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
                call init_single_variable(my_c2b_3d_profile%Fnumu_bar_z, fpath, fnbase_list(FIL_Fnumu_bar_z), dsname_list(FIL_Fnumu_bar_z), &
                    num_refine_levels, num_datablocks, iter_num, 'FFF', save_cpu_or_mem)
            end if
        end if 


    end subroutine init_c2b_3d_profile


    subroutine detect_level_and_blocks_to_be_loaded(ixI^L, x, cp3dvb, nref, ndabk, mask_array, mirrorz)
        include 'amrvacdef.f'
        integer, intent(in)                        :: ixI^L
        double precision, intent(in)               :: x(ixI^S, 1:ndim)
        type(carpet_variable), intent(inout)       :: cp3dvb
        integer, intent(in)                        :: nref, ndabk
        integer, intent(out)                       :: mask_array(0:nref-1, 0:ndabk-1)
        logical, intent(in)                        :: mirrorz
        ! local 
        integer                                    :: which_level, which_datablock, my_error_type
        integer                                    :: ix^D
        double precision                           :: x3

        mask_array = 0
        {do ix^D=ixImin^D, ixImax^D \}
            if (mirrorz) then
                x3 = abs(x(ix^D, 3)) ! z-symmetry
            else
                x3 = x(ix^D, 3) ! z-symmetry
            endif
            call find_point(x(ix^D, 1), x(ix^D, 2), x3, cp3dvb, which_level, which_datablock, my_error_type)
            if (my_error_type==0) mask_array(which_level, which_datablock) = 1
        {end do^D&\}
    end subroutine detect_level_and_blocks_to_be_loaded

    subroutine destroy_c2b_3d_profile()
        include 'amrvacdef.f'
        integer :: Nspecies
        if (mype==0) write(*, *) "Destroying c2b_3d_profile"
        Nspecies = ^NS
        call destroy_single_variable(my_c2b_3d_profile%rho)
        call destroy_single_variable(my_c2b_3d_profile%velox)
        call destroy_single_variable(my_c2b_3d_profile%veloy)
        call destroy_single_variable(my_c2b_3d_profile%veloz)
        call destroy_single_variable(my_c2b_3d_profile%psi)
        call destroy_single_variable(my_c2b_3d_profile%alpha)
        call destroy_single_variable(my_c2b_3d_profile%betax)
        call destroy_single_variable(my_c2b_3d_profile%betay)
        call destroy_single_variable(my_c2b_3d_profile%betaz)
        call destroy_single_variable(my_c2b_3d_profile%temp)
        call destroy_single_variable(my_c2b_3d_profile%ye)
        call destroy_single_variable(my_c2b_3d_profile%ymu)
        call destroy_single_variable(my_c2b_3d_profile%eps)
        call destroy_single_variable(my_c2b_3d_profile%pressure)
        call destroy_single_variable(my_c2b_3d_profile%Ax)
        call destroy_single_variable(my_c2b_3d_profile%Ay)
        call destroy_single_variable(my_c2b_3d_profile%Az)
        call destroy_single_variable(my_c2b_3d_profile%Bx)
        call destroy_single_variable(my_c2b_3d_profile%By)
        call destroy_single_variable(my_c2b_3d_profile%Bz)
        if(Nspecies .ge. 1) then
        call destroy_single_variable(my_c2b_3d_profile%Nnue)
        call destroy_single_variable(my_c2b_3d_profile%Enue)
        call destroy_single_variable(my_c2b_3d_profile%Fnue_x)
        call destroy_single_variable(my_c2b_3d_profile%Fnue_y)
        call destroy_single_variable(my_c2b_3d_profile%Fnue_z)
        end if 
        if(Nspecies .ge. 2) then
        call destroy_single_variable(my_c2b_3d_profile%Nnue_bar)
        call destroy_single_variable(my_c2b_3d_profile%Enue_bar)
        call destroy_single_variable(my_c2b_3d_profile%Fnue_bar_x)
        call destroy_single_variable(my_c2b_3d_profile%Fnue_bar_y)
        call destroy_single_variable(my_c2b_3d_profile%Fnue_bar_z)
        end if 
        if(Nspecies .ge. 3) then
        call destroy_single_variable(my_c2b_3d_profile%Nnux)
        call destroy_single_variable(my_c2b_3d_profile%Enux)
        call destroy_single_variable(my_c2b_3d_profile%Fnux_x)
        call destroy_single_variable(my_c2b_3d_profile%Fnux_y)
        call destroy_single_variable(my_c2b_3d_profile%Fnux_z)
        end if 
        if(Nspecies .ge. 4) then
        call destroy_single_variable(my_c2b_3d_profile%Nnumu)
        call destroy_single_variable(my_c2b_3d_profile%Enumu)
        call destroy_single_variable(my_c2b_3d_profile%Fnumu_x)
        call destroy_single_variable(my_c2b_3d_profile%Fnumu_y)
        call destroy_single_variable(my_c2b_3d_profile%Fnumu_z)
        end if 
        if(Nspecies .ge. 5) then
        call destroy_single_variable(my_c2b_3d_profile%Nnumu_bar)
        call destroy_single_variable(my_c2b_3d_profile%Enumu_bar)
        call destroy_single_variable(my_c2b_3d_profile%Fnumu_bar_x)
        call destroy_single_variable(my_c2b_3d_profile%Fnumu_bar_y)
        call destroy_single_variable(my_c2b_3d_profile%Fnumu_bar_z)
        end if 

    end subroutine destroy_c2b_3d_profile

end module Init_c2b_cart_3d_interface


module NSimport
    implicit none
    ! Global quantities 
    integer                           :: profile_type         = 3       ! 1: ctwob, 2: xns, 3: ctwob3d_cart
    double precision                  :: shift_origin         = 0.0d0       ! shift origin when initializing
      !> quantities related to refine methods
      integer                           :: refine_type          = 3       
      !double precision, parameter       :: r_entire             = 6.0d0  !15.0d0 !30.0d0  ! radius that must cover the entire star
      !double precision, parameter       :: r_outer              = 15.0d0  ! Outside the radius that refinement lv would be the lowest
      double precision, parameter       :: r_entire             = 20.0d0 !15.0d0 !30.0d0  ! radius that must cover the entire star
      double precision, parameter       :: r_outer              = 30.0d0 !20.0d0  ! Outside the radius that refinement lv would be the lowest
      !double precision, parameter       :: r_entire             = 3.0d0  !15.0d0 !30.0d0  ! radius that must cover the entire star
      !double precision, parameter       :: r_outer              = 15.0d0  ! Outside the radius that refinement lv would be the lowest
      double precision, parameter       :: r_core               = 1.0d0   ! radius of the core that you want to coarsen
      logical                           :: cut_r_core           = .false. ! want to cut r_core to increase dt?
      !> quantities related to B-field
      logical                           :: init_B_from_vecA     = .false. !.true.
      !> m1
      logical                           :: init_m1_vars = .true. !KEN CHANGE dont forget

    ! parameters used only in cases that init with XNS results
    integer                           :: xns_profile_type     = 3       ! 1: pure hydro, 2: torodial stored in hydroeq.dat, 
    character*128                     :: profile_path_xns     = '/mnt/rafast/hng/BHAC-FIL_handover/xns/purep_smallm'
    !character*128                     :: profile_path_xns     = '/mnt/rafast/hng/BHAC-FIL_handover/xns/mag_diff_NS_2048x128'
    ! parameters used only in cases that init with Carpet2BHAC data in 2d case
    character*128                     :: profile_path_met     = './' 
    character*128                     :: profile_path_hydro   = './' 
    character*128                     :: hydro_type           = 'prim'
    logical                           :: read_C2B_metric      = .false.

    ! ------------------ parameters used only in cases that init with 3D FIL output data -------------------------!
    ! 4th order, 73 or 58ms
    !character(len=256)                :: fpath_c2b_3dcart="/lustre/hpe/ws11/ws11.1/ws/xfpchaba-handoff_new/handoff_copy/handoff-runs/0.2M-res-q1-2.55M/mhd_full_seed_late_cy/output-0005/data_hdf5_3D"
    !integer                           :: iteration_num = 255744

    ! 4th order, -0-43ms
    !character(len=256)                :: fpath_c2b_3dcart = "/lustre/hpe/ws10/ws10.1/ws/xfpchaba-handoff/handoff-runs/0.2M-res-q1-2.55M/mhd_full_seed_late/output-0002/data_hdf5_3D/"
    !integer                           :: iteration_num = 101376   

    !character(len=256) :: fpath_c2b_3dcart = "/lustre/hpe/ws10/ws10.3/ws/xfpjiang-jiangbns/init_data/FIL/0.2-mhd-hybrid/3D/"
    !character(len=256) :: fpath_c2b_3dcart = "/lustre/hpe/ws10/ws10.3/ws/xfpchaba-handoff-project/0.2M-res-q1-2.55M-2ndorder/mhd/output-0003/data_hdf5_3D/"
       
    ! character(len=256) :: fpath_c2b_3dcart = "/lustre/hpe/ws11/ws11.1/ws/xfpjiang-jiang-ho/init_data/FIL/mhd_full_seed_late_cy/3D"
   
  ! character(len=256) :: fpath_c2b_3dcart =  "/mnt/raarchive/cassing/FIL_runs/NS_HOT_M1_plasma_NEW/output-0000/data_hdf5_3D" !"/mnt/rafast/miler/bhac_import_data/KEN_TOV/"

 ! character(len=256) :: fpath_c2b_3dcart = "/mnt/raarchive/cassing/FIL_runs/M130_TNTYST_LOWRES_BNS/Supermuc/output0009/Reduced3D_it196608" !data_hdf5_3D

 ! "/mnt/raarchive/cassing/FIL_runs/NS_HOT_M1_plasma_NEW/output-0001/data_hdf5_3D" !"/mnt/rafast/miler/bhac_import_data/KEN_TOV/"
 
    !character(len=256) :: fpath_c2b_3dcart = "/mnt/rafast/miler/FIL_runs/NS_TOV_HOT/output-0000/data_hdf5_3D/"
    !"/mnt/raarchive/cassing/FIL_runs/NS_3D_testM1/NS_simTEST/output-0000/data_hdf5_3D"
    !"/mnt/raarchive/cassing/FIL_runs/HARRY_DATA/TNTYST_postmerger/3D"
 
    !  character(len=256) :: fpath_c2b_3dcart = "/scratch/astro/mcassing/FIL_data/TNTYST_postmerger/3D"
 
    ! character(len=256) :: fpath_c2b_3dcart = "/gpfs/scratch/ehpc590/Marie/FIL_runs/NS_HOT_M1_plasma_NEW/output-0001/data_hdf5_3D"

character(len=256) :: fpath_c2b_3dcart = "/mnt/raarchive/hng/FIL_data/new_03npe_m1_SFHo_midres/3D_data_output-0014/data_hdf5_3D"
!"/mnt/raarchive/hng/FIL_data/new_03npe_m1_DD2_midres/output-0001/data_hdf5_3D"

      integer                           :: iteration_num = 358912 !531456 ! 140032 !512 !1024 !200448 ! the iteration number you want to read, before reading, please use h5ls to check what iteration you want exists in that file
   
    !************************
   !! setup for M1 1 species:
   !************************
   !!! list of strings that contain the name of the saved variables, like 'HYDROBASE::rho', "HYDROBASE::vel[1]", and so on
   !! character(len=256), dimension(24) :: variable_name_list =  [ character(len=256):: "HYDROBASE::rho", "HYDROBASE::vel[0]", &
   !!     "HYDROBASE::vel[1]", "HYDROBASE::vel[2]", "ANTELOPE::W_z4c", "ADMBASE::alp", "ADMBASE::betax", "ADMBASE::betay", "ADMBASE::betaz", &
   !!     "HYDROBASE::temperature", "HYDROBASE::Y_e", "HYDROBASE::eps", "ILLINOISGRMHD::P", &
   !!     "ILLINOISGRMHD::Ax", "ILLINOISGRMHD::Ay", "ILLINOISGRMHD::Az", "HYDROBASE::Bvec[0]", "HYDROBASE::Bvec[1]", "HYDROBASE::Bvec[2]", &
   !!     "FRANKFURT_M1::Nnue","FRANKFURT_M1::Enue","FRANKFURT_M1::Fnue_x", "FRANKFURT_M1::Fnue_y", "FRANKFURT_M1::Fnue_z"]
   !!     ! list of strings that contain the name of the file where the cooresponding variable are saved, like 'rho', 'vel[1]', and so on
   !! character(len=256), dimension(24) :: variable_fname_list = [ character(len=256):: "rho", "vel[0]", &
   !!     "vel[1]", "vel[2]", "W_z4c", "alp", "betax", "betay", "betaz", "temperature", 'Y_e', 'eps', 'P', &
   !!     "Ax", "Ay", "Az", "Bvec[0]", "Bvec[1]", "Bvec[2]","Nnue","Enue","Fnue_x","Fnue_y","Fnue_z"]

    !************************
   !! setup for M1 2 species:
   !************************
   !!! list of strings that contain the name of the saved variables, like 'HYDROBASE::rho', "HYDROBASE::vel[1]", and so on
   ! character(len=256), dimension(29) :: variable_name_list =  [ character(len=256):: "HYDROBASE::rho", "HYDROBASE::vel[0]", &
   !     "HYDROBASE::vel[1]", "HYDROBASE::vel[2]", "ANTELOPE::W_z4c", "ADMBASE::alp", "ADMBASE::betax", "ADMBASE::betay", "ADMBASE::betaz", &
   !     "HYDROBASE::temperature", "HYDROBASE::Y_e", "HYDROBASE::eps", "ILLINOISGRMHD::P", &
   !     "ILLINOISGRMHD::Ax", "ILLINOISGRMHD::Ay", "ILLINOISGRMHD::Az", "HYDROBASE::Bvec[0]", "HYDROBASE::Bvec[1]", "HYDROBASE::Bvec[2]", &
   !     "FRANKFURT_M1::Nnue","FRANKFURT_M1::Enue","FRANKFURT_M1::Fnue_x", "FRANKFURT_M1::Fnue_y", "FRANKFURT_M1::Fnue_z", &
   !     "FRANKFURT_M1::Nnue_bar","FRANKFURT_M1::Enue_bar","FRANKFURT_M1::Fnue_bar_x","FRANKFURT_M1::Fnue_bar_y","FRANKFURT_M1::Fnue_bar_z"]
   !     ! list of strings that contain the name of the file where the cooresponding variable are saved, like 'rho', 'vel[1]', and so on
   ! character(len=256), dimension(29) :: variable_fname_list = [ character(len=256):: "rho", "vel[0]", &
   !     "vel[1]", "vel[2]", "W_z4c", "alp", "betax", "betay", "betaz", "temperature", 'Y_e', 'eps', 'P', &
   !     "Ax", "Ay", "Az", "Bvec[0]", "Bvec[1]", "Bvec[2]","Nnue","Enue","Fnue_x","Fnue_y","Fnue_z", &
   !     "Nnue_bar","Enue_bar","Fnue_bar_x","Fnue_bar_y","Fnue_bar_z"]
   !************************
    !! setup for M1 5 species:
   !************************
   ! list of strings that contain the name of the saved variables, like 'HYDROBASE::rho', "HYDROBASE::vel[1]", and so on
    character(len=256), dimension(45) :: variable_name_list =  [ character(len=256):: "ILLINOISGRMHD::rho_b", "HYDROBASE::vel[0]", &
        "HYDROBASE::vel[1]", "HYDROBASE::vel[2]", "ANTELOPE::psi_z4c", "ADMBASE::alp", "ADMBASE::betax", "ADMBASE::betay", "ADMBASE::betaz", &
        "ILLINOISGRMHD::temp", "ILLINOISGRMHD::ye", "ILLINOISGRMHD::ymu", "HYDROBASE::eps", "HYDROBASE::press", &
        "ILLINOISGRMHD::Ax", "ILLINOISGRMHD::Ay", "ILLINOISGRMHD::Az", "ILLINOISGRMHD::Bx", "ILLINOISGRMHD::By", "ILLINOISGRMHD::Bz", &
        "FRANKFURT_M1::Nnue","FRANKFURT_M1::Enue","FRANKFURT_M1::Fnue_x", "FRANKFURT_M1::Fnue_y", "FRANKFURT_M1::Fnue_z", &
        "FRANKFURT_M1::Nnue_bar","FRANKFURT_M1::Enue_bar","FRANKFURT_M1::Fnue_bar_x","FRANKFURT_M1::Fnue_bar_y","FRANKFURT_M1::Fnue_bar_z", &
        "FRANKFURT_M1::Nnux","FRANKFURT_M1::Enux", "FRANKFURT_M1::Fnux_x", "FRANKFURT_M1::Fnux_y", "FRANKFURT_M1::Fnux_z", &
        "FRANKFURT_M1::Nnumu","FRANKFURT_M1::Enumu", "FRANKFURT_M1::Fnumu_x", "FRANKFURT_M1::Fnumu_y", "FRANKFURT_M1::Fnumu_z", &
        "FRANKFURT_M1::Nnumu_bar","FRANKFURT_M1::Enumu_bar", "FRANKFURT_M1::Fnumu_bar_x", "FRANKFURT_M1::Fnumu_bar_y", "FRANKFURT_M1::Fnumu_bar_z"]
        ! list of strings that contain the name of the file where the cooresponding variable are saved, like 'rho', 'vel[1]', and so on
    character(len=256), dimension(45) :: variable_fname_list = [ character(len=256):: "rho_b", "vel[0]", &
        "vel[1]", "vel[2]", "psi_z4c", "alp", "betax", "betay", "betaz", "temp", 'ye', 'ymu', 'eps', 'press', &
        "Ax", "Ay", "Az", "Bx", "By", "Bz","Nnue","Enue","Fnue_x","Fnue_y","Fnue_z", &
        "Nnue_bar","Enue_bar","Fnue_bar_x","Fnue_bar_y","Fnue_bar_z",&
        "Nnux", "Enux", "Fnux_x","Fnux_y","Fnux_z", &
        "Nnumu", "Enumu", "Fnumu_x","Fnumu_y","Fnumu_z", &
        "Nnumu_bar", "Enumu_bar", "Fnumu_bar_x","Fnumu_bar_y","Fnumu_bar_z"] !KEN changed temp

    integer                           :: num_ref_levels = 8           ! number of refinement levels, to check this, you can use h5ls to check how many refinement levels are there
    integer                           :: num_data_blocks = 40 !80!16!96!80        ! number of blocks that the data have been separated into
    integer                           :: save_cpu_or_mem = 3          ! 1 to save more cpu time for faster, 2 to save more memory, 3 is the best
    logical                           :: mirror_zplane = .false.      ! used in the case that: originally in ETK data there is a z-symmetry, but you don't want it in BHAC
    logical                           :: init_magnetic_field = .true. ! whether to initialize the magnetic field 
    logical                           :: use_eps_init = .false. !KEN.true.        ! if there is no eps, temperature should be provide in the case of realistic eos, and should provide pressure in other eos cases
    character*256                     :: interp_fluid = 'Linear', interp_metric='Linear', interp_vectA = 'Linear', interp_m1 = 'Linear' ! allowed: 'Linear', 'Hermitian3', "Lagrangian2", "Lagrangian3", "Lagrangian4", "Lagrangian5", "Automatic"
    !character*256                     :: interp_fluid = 'Lagrangian3', interp_metric='Lagrangian3', interp_vectA = 'Hermitian3' ! allowed: 'Linear', 'Hermitian3', "Lagrangian2", "Lagrangian3", "Lagrangian4", "Lagrangian5", "Automatic"
    ! ------------------ End of parameters used only in cases that init with 3D FIL output data ------------------!

    
    contains


    subroutine set_atmo(ixI^L,ixO^L,w,x)
        ! The atmosphere treatment
        use mod_eos

        include 'amrvacdef.f'

        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: x(ixI^S,1:ndim)
        double precision, intent(inout)  :: w(ixI^S,1:nw)
        ! .. local ..
        double precision                 :: xBL(ixI^S,1:ndim)
        logical, dimension(ixI^S)        :: patchw
        double precision, parameter      :: tol_atmo = 0.001d0
        double precision, parameter      :: radial_shift = 1.d0
        {#IFDEF ELECTRONS
              double precision                 :: Pe(ixI^S), sefl(ixI^S)
        }
        !-----------------------------------------------------------------------------
        patchw(ixO^S) = .false.

        if (eos_uses_ye()) then
          where(w(ixO^S,rho_) .lt. small_rho_thr)
            w(ixO^S,rho_)    = small_rho
            w(ixO^S,u1_)     = 0.0d0
            w(ixO^S,u2_)     = 0.0d0
            w(ixO^S,u3_)     = 0.0d0
            w(ixO^S,ye_)     = big_ye
            w(ixO^S,T_eps_)  = small_temp
            patchw(ixO^S)    = .true.
          end where
          if (eos_has_ymu()) then
            where(patchw(ixO^S))
              w(ixO^S,ymu_) = eos_ymumin
            end where
          else
            where(patchw(ixO^S))
              w(ixO^S,ymu_) = 0.0d0
            end where
          endif
        else
          where(w(ixO^S,rho_) .lt. small_rho_thr)
            w(ixO^S,rho_)    = small_rho
            w(ixO^S,u1_)     = 0.0d0
            w(ixO^S,u2_)     = 0.0d0
            w(ixO^S,u3_)     = 0.0d0
            w(ixO^S,T_eps_)  = small_eps
            patchw(ixO^S)    = .true.
          end where
        endif


    end subroutine set_atmo
    !=============================================================================

    subroutine read_NS_CTWOB(ixI^L,ixO^L,w,x)
        ! import of NS data from phi-averaged c2b file
        !use mod_carpetgrid
        use mod_oneblock, only: interpolate_oneblock
        use mod_xns
        use mod_c2b
        use mod_eos
  
        include 'amrvacdef.f'
        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: x(ixI^S,1:ndim)
        double precision, intent(inout)  :: w(ixI^S,1:nw)
        !--------------------------------------
        double precision :: rBL, phi
        integer          :: ix^D, i
        double precision :: xCart(ixI^S, 1:ndim)
        double precision :: uCart(ixI^S, 1:ndir)
        double precision :: uBL(ixI^S, 1:ndir)
        double precision :: vD(ixI^S, 1:ndir)
        double precision :: vD_perturb(ixI^S, 1:ndir)
        double precision :: vU_perturb(ixI^S, 1:ndir)
        double precision, dimension(ixI^S) :: sqrV
        double precision :: eps, rho, press, vel1, vel2, vel3, ye, temp, lfac
        double precision :: rad, theta
  
        !---------------------------------------------------------
        ! Not change yet
        {do ix^D=ixOmin^D, ixOmax^D \}
            rad = x(ix^D,r_)
            theta = x(ix^D, ^Z)
            call mod_c2b_map_2D(w(ix^D,rho_),rad,theta,c2b_prof%rho(:,:), &
                                c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
            w(ix^D, rho_) = max(eos_rhomin, min(w(ix^D, rho_), eos_rhomax))
            if (eos_uses_ye()) then
              call mod_c2b_map_2D(w(ix^D,ye_),rad,theta,c2b_prof%ye(:,:), &
                                  c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
              w(ix^D, ye_) = max(eos_yemin, min(w(ix^D, ye_), eos_yemax))
              if (eos_has_ymu()) then
                w(ix^D, ymu_) = eos_ymumin
              else
                w(ix^D, ymu_) = 0.0d0
              endif
              call mod_c2b_map_2D(w(ix^D,T_eps_),rad,theta,c2b_prof%temp(:,:), &
                                  c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
            endif
            w(ix^D, T_eps_) = max(eos_tempmin, min(w(ix^D, T_eps_), eos_tempmax))
    
            call mod_c2b_map_2D(w(ix^D,u1_),rad,theta,c2b_prof%v1(:,:), &
                                c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
            call mod_c2b_map_2D(w(ix^D,u2_),rad,theta,c2b_prof%v2(:,:), &
                                c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
            call mod_c2b_map_2D(w(ix^D,u3_),rad,theta,c2b_prof%v3(:,:), &
                                c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)

        {end do \}

    !call check_data_correctness(ixI^L, ixO^L, w, x)

    end subroutine read_NS_CTWOB

    subroutine read_NS_XNS_cart(ixI^L,ixO^L,w,x)
        use mod_xns
        use mod_eos
        include 'amrvacdef.f'
        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: x(ixI^S,1:ndim)
        double precision, intent(inout)  :: w(ixI^S,1:nw)
        !--------------------------------------
        integer          :: ix^D, idir
        double precision :: rad, theta, rad_xy, sin_phi, cos_phi, sin_theta, shifted_x, shifted_y
        double precision :: vphi(ixI^S), betaphi(ixI^S), v(ixI^S,1:ndir), gamma(ixI^S,1:3,1:3), lfac(ixI^S), prs_tmp(ixI^S)
        if (shift_origin .ne. 0.0d0) call mpistop('shift origin must be zero')
        if (eos_uses_ye()) call mpistop('read_NS_XNS_cart requires an EOS without Ye/Ymu dependence')
  
        {do ix^D=ixOmin^D, ixOmax^D \}
            
            !shifted_x = x(ix^D, 1:ndim)-shift_origin
            !rad = dsqrt( sum(shifted_x(1:ndim)**2) )
            !theta = acos( shifted_x(3) / rad )
            !rad_xy = dsqrt( sum(shifted_x(1:2)**2) )

            shifted_x = x(ix^D, 1) - shift_origin
            shifted_y = x(ix^D, 2) - shift_origin
            rad = dsqrt(x(ix^D,3)**2 + shifted_x**2 + shifted_y**2 )
            rad_xy = dsqrt( shifted_x**2 + shifted_y**2 )
            theta = acos( x(ix^D,3) / rad )
            !theta  = dmod(2.0d0*dpi + datan2(rad_xy,x(ix^D,3)),    2.0d0*dpi)  this leads to dots in staggered grid
            sin_phi = shifted_y / rad_xy
            cos_phi = shifted_x / rad_xy
            sin_theta = rad_xy / rad
   
            !assume z_symm
            if ( theta > 0.5d0 * dpi ) theta = dpi - theta
   
            call mod_XNS_map_2D(w(ix^D,rho_),rad, theta,prof%rho(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            call mod_XNS_map_2D(prs_tmp(ix^D),rad, theta,prof%press(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            call mod_XNS_map_2D(w(ix^D,alp_metric_),rad, theta,prof%alp(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            call mod_XNS_map_2D(w(ix^D,psi_metric_),rad, theta,prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            call mod_XNS_map_2D(betaphi(ix^D),rad, theta,prof%beta3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            call mod_XNS_map_2D(vphi(ix^D),rad, theta,prof%v3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)

            !call mod_XNS_map_2D(w(ix^D,b1_),rad, theta,prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            !call mod_XNS_map_2D(w(ix^D,b2_),rad, theta,prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            !call mod_XNS_map_2D(w(ix^D,b3_),rad, theta,prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
   
            w(ix^D, alp_metric_) = min(w(ix^D, alp_metric_), 1.0d0) 
            w(ix^D, psi_metric_) = max(w(ix^D, psi_metric_), 1.0d0) 
   
            vphi(ix^D) = vphi(ix^D) * rad * sin_theta
            betaphi(ix^D) = betaphi(ix^D) * rad * sin_theta

            ! fixme: v not used here
            !v(ix^D, 1) = - sin_phi * vphi(ix^D)
            !v(ix^D, 2) =   cos_phi * vphi(ix^D)
            ! this u^C is without W
            w(ix^D, u1_) = - sin_phi * vphi(ix^D)
            w(ix^D, u2_) =   cos_phi * vphi(ix^D)
            w(ix^D, beta_metric1_) = - sin_phi * betaphi(ix^D)
            w(ix^D, beta_metric2_) =   cos_phi * betaphi(ix^D)
   
            ! B-field
            !call mod_XNS_map_2D(Br(ix^D),rad, theta,prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            !call mod_XNS_map_2D(Btheta(ix^D),rad, theta,prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
            !call mod_XNS_map_2D(Bphi(ix^D),rad, theta,prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)

            call eos_get_eps_one_grid(prs_tmp(ix^D),w(ix^D, rho_),w(ix^D, T_eps_))
            w(ix^D, T_eps_) = max( w(ix^D, T_eps_), small_eps )
        {end do^D&\}

    end subroutine read_NS_XNS_cart

    subroutine read_NS_3DCarpet_cart(ixI^L,ixO^L,w,x)
        use mod_eos
        use mod_c2b_cart3d
        use Init_c2b_cart_3d_interface
        use mod_imhd_intermediate
        use mod_m1_metric_interface
        use mod_m1_closure
        include 'amrvacdef.f'

        integer, intent(in)              :: ixI^L, ixO^L
        double precision, intent(in)     :: x(ixI^S,1:ndim)
        double precision, intent(inout)  :: w(ixI^S,1:nw)
        !--------------------------------------
        integer          :: ix^D, idir
        double precision :: gamma(ixI^S, 1:3, 1:3), x1, x2, x3, lfac(ixI^S), prs_tmp(ixI^S), eps_tmp(ixI^S)
        integer          :: mask_array(0:num_ref_levels-1, 0:num_data_blocks-1)
        integer          :: Nspecies
        double precision :: Gamma_M1(ixI^S) !KEN
        integer          :: i !KEN
        
	    type(m1_metric_helper) :: metricM1 
        type(m1_closure_helpers) :: stateM1 !KEN used this
    
        Nspecies = ^NS
        !write(*,*)"read_NS_3D: Number of species : NS",^NS

        ! use rho as an example to detect how many data blocks should be loaded
        call detect_level_and_blocks_to_be_loaded(ixI^L, x, my_c2b_3d_profile%rho, num_ref_levels, num_data_blocks, mask_array, mirror_zplane)
        !write(*, *) "mask array: ", mask_array
        ! load all the variables that is needed
        call load_selected_ref_and_datablock(my_c2b_3d_profile%rho, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%velox, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%veloy, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%veloz, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%psi, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%alpha, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%betax, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%betay, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%betaz, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%ye, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%eps, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%temp, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%pressure, num_ref_levels, num_data_blocks, mask_array)
        if(init_m1_vars)then
            if(Nspecies .ge. 1) then
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Nnue, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Enue, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnue_x, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnue_y, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnue_z, num_ref_levels, num_data_blocks, mask_array)
            end if 
            if(Nspecies .ge. 2) then
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Nnue_bar, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Enue_bar, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnue_bar_x, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnue_bar_y, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnue_bar_z, num_ref_levels, num_data_blocks, mask_array)
            end if 
            if(Nspecies .ge. 3) then
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Nnux, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Enux, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnux_x, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnux_y, num_ref_levels, num_data_blocks, mask_array)
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Fnux_z, num_ref_levels, num_data_blocks, mask_array)
            end if 
        end if
        if (.not. (allocated(my_c2b_3d_profile%rho%all_refine_levels) .or. (allocated(my_c2b_3d_profile%psi%all_refine_levels)))) then
            call mpistop("FIL reader data not initialize successfully, consider wrong path file or variable name!")
        endif

        ! interpolation
        w(ixI^S, :) = 0.0d0
        !w(ixI^S, lfac_) = 1.0d0
        w(ixI^S, psi_metric_) = 1.0d0

        {do ix^D=ixOmin^D, ixOmax^D \}
            x1 = x(ix^D, 1)
            x2 = x(ix^D, 2)
            if (mirror_zplane) then
                x3 = abs(x(ix^D, 3)) ! z-symmetry
            else
                x3 = x(ix^D, 3) ! z-symmetry
            endif
            call interp_cart3d_carpet_variable_one_point(w(ix^D, rho_), x1, x2, x3, my_c2b_3d_profile%rho, interp_fluid, eos_rhomin+10*tiny(1.0d0))
            call interp_cart3d_carpet_variable_one_point(w(ix^D, u1_), x1, x2, x3, my_c2b_3d_profile%velox, interp_fluid, 0.0d0)
            call interp_cart3d_carpet_variable_one_point(w(ix^D, u2_), x1, x2, x3, my_c2b_3d_profile%veloy, interp_fluid, 0.0d0)
            call interp_cart3d_carpet_variable_one_point(w(ix^D, u3_), x1, x2, x3, my_c2b_3d_profile%veloz, interp_fluid, 0.0d0)
            call interp_cart3d_carpet_variable_one_point(w(ix^D, psi_metric_), x1, x2, x3, my_c2b_3d_profile%psi, interp_metric, 1.0d0)
            call interp_cart3d_carpet_variable_one_point(w(ix^D, alp_metric_), x1, x2, x3, my_c2b_3d_profile%alpha, interp_metric, 1.0d0)
            call interp_cart3d_carpet_variable_one_point(w(ix^D, beta_metric1_), x1, x2, x3, my_c2b_3d_profile%betax, interp_metric, 0.0d0)
            call interp_cart3d_carpet_variable_one_point(w(ix^D, beta_metric2_), x1, x2, x3, my_c2b_3d_profile%betay, interp_metric, 0.0d0)
            call interp_cart3d_carpet_variable_one_point(w(ix^D, beta_metric3_), x1, x2, x3, my_c2b_3d_profile%betaz, interp_metric, 0.0d0)
            if (eos_uses_ye()) then
                call interp_cart3d_carpet_variable_one_point(w(ix^D, ye_), x1, x2, x3, my_c2b_3d_profile%ye, interp_fluid, eos_yemin+10*tiny(1.0d0))
                if (eos_has_ymu()) then
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, ymu_), x1, x2, x3, my_c2b_3d_profile%ymu, interp_fluid, eos_ymumin)
                else
                    w(ix^D, ymu_) = 0.0d0
                endif
                if (use_eps_init) then
                    call interp_cart3d_carpet_variable_one_point(eps_tmp(ix^D), x1, x2, x3, my_c2b_3d_profile%eps, interp_fluid, eos_epsmin+10*tiny(1.0d0))
                    if (eos_has_ymu()) then
                        call eos_get_temp_one_grid(w(ix^D, rho_), eps_tmp(ix^D), w(ix^D, T_eps_), w(ix^D, ye_), ymu=w(ix^D, ymu_))
                    else
                        call eos_get_temp_one_grid(w(ix^D, rho_), eps_tmp(ix^D), w(ix^D, T_eps_), w(ix^D, ye_))
                    endif
                else
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, T_eps_), x1, x2, x3, my_c2b_3d_profile%temp, interp_fluid, eos_tempmin+10*tiny(1.0d0))
                endif
            else  ! single star with idealgas?
                if (use_eps_init) then
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, T_eps_), x1, x2, x3, my_c2b_3d_profile%eps, interp_fluid, small_eps)
                else
                    call interp_cart3d_carpet_variable_one_point(prs_tmp(ix^D), x1, x2, x3, my_c2b_3d_profile%pressure, interp_fluid, small_press)
                    call eos_get_eps_one_grid(prs_tmp(ix^D), w(ix^D, rho_), w(ix^D, T_eps_))
                endif
            endif
            if(init_m1_vars) then
                {^KSP&
                if(^KSP.eq. 1) then
                    !write(*,*)"KSP =1 ",^KSP
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, erad^KSP_), x1, x2, x3, my_c2b_3d_profile%Enue, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, nrad^KSP_), x1, x2, x3, my_c2b_3d_profile%Nnue, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP1_), x1, x2, x3, my_c2b_3d_profile%Fnue_x, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP2_), x1, x2, x3, my_c2b_3d_profile%Fnue_y, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP3_), x1, x2, x3, my_c2b_3d_profile%Fnue_z, interp_m1, 0.0d0)
                end if 
                if(^KSP .eq. 2) then
                    !write(*,*)"KSP =2 ",^KSP
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, erad^KSP_), x1, x2, x3, my_c2b_3d_profile%Enue_bar, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, nrad^KSP_), x1, x2, x3, my_c2b_3d_profile%Nnue_bar, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP1_), x1, x2, x3, my_c2b_3d_profile%Fnue_bar_x, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP2_), x1, x2, x3, my_c2b_3d_profile%Fnue_bar_y, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP3_), x1, x2, x3, my_c2b_3d_profile%Fnue_bar_z, interp_m1, 0.0d0)
                end if 
                if(^KSP .eq. 3) then
                    !write(*,*)"KSP =3 ",^KSP
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, erad^KSP_), x1, x2, x3, my_c2b_3d_profile%Enux, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, nrad^KSP_), x1, x2, x3, my_c2b_3d_profile%Nnux, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP1_), x1, x2, x3, my_c2b_3d_profile%Fnux_x, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP2_), x1, x2, x3, my_c2b_3d_profile%Fnux_y, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP3_), x1, x2, x3, my_c2b_3d_profile%Fnux_z, interp_m1, 0.0d0)
                end if 
                if(^KSP .eq. 4) then
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, erad^KSP_), x1, x2, x3, my_c2b_3d_profile%Enumu, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, nrad^KSP_), x1, x2, x3, my_c2b_3d_profile%Nnumu, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP1_), x1, x2, x3, my_c2b_3d_profile%Fnumu_x, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP2_), x1, x2, x3, my_c2b_3d_profile%Fnumu_y, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP3_), x1, x2, x3, my_c2b_3d_profile%Fnumu_z, interp_m1, 0.0d0)
                end if 
                if(^KSP .eq. 5) then
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, erad^KSP_), x1, x2, x3, my_c2b_3d_profile%Enumu_bar, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, nrad^KSP_), x1, x2, x3, my_c2b_3d_profile%Nnumu_bar, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP1_), x1, x2, x3, my_c2b_3d_profile%Fnumu_bar_x, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP2_), x1, x2, x3, my_c2b_3d_profile%Fnumu_bar_y, interp_m1, 0.0d0)
                    call interp_cart3d_carpet_variable_one_point(w(ix^D, frad^KSP3_), x1, x2, x3, my_c2b_3d_profile%Fnumu_bar_z, interp_m1, 0.0d0)
                end if 
                \}
            end if 

            if (mirror_zplane .and. x(ix^D, 3)<0.d0) then
                w(ix^D, u3_) = -w(ix^D, u3_)
                w(ix^D, beta_metric3_) = -w(ix^D, beta_metric3_)
                !if(init_m1_vars)then
                !    {^KSP&
                !    if(Nspecies .ge. 1) then
                !    w(ix^D, frad^KSP3_) = - w(ix^D, frad^KSP3_)
                !    end if 
                !    if(Nspecies .ge. 2) then
                !    w(ix^D, frad^KSP3_) = - w(ix^D, frad^KSP3_)
                !    end if 
                !    if(Nspecies .ge. 3) then
                !    w(ix^D, frad^KSP3_) = - w(ix^D, frad^KSP3_)
                !    end if 
                !    \}
                !end if
            endif
        {end do^D& \}

        ! Warning: convert conformal factor from FIL to BHAC, please confirm this formula with the person who give you the FIL data
        !w(ixO^S, psi_metric_) = dsqrt(1.0d0/w(ixO^S, psi_metric_))
        ! psi_z4c already psi, W_z4c need convert

        ! calculate lorentz factor and velocity, note here we can not use dysp_get_lfac which assume w(.., ui)=WUi already
        if (eos_type/=hybrid) then
            call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S, 1:3, 1:3))
            lfac(ixO^S) = 1.d0
            do idir = 1, ndir
                gamma(ixO^S, idir, idir) = gamma(ixO^S, idir, idir)*w(ixO^S, psi_metric_)**4
                lfac(ixO^S) = lfac(ixO^S)-gamma(ixO^S, idir, idir)*w(ixO^S, u1_+idir-1)**2
            enddo
            lfac(ixO^S) = dsqrt(1.0d0/lfac(ixO^S))
            {^C&w(ixO^S, u^C_) = lfac(ixO^S)*w(ixO^S, u^C_)\}
        endif

        call metricM1%fill_metric(w,x,ixI^L,ixO^L)  
           
        {^KSP&

           {do ix^D=ixOmin^D, ixOmax^D \}     
           	  stateM1%E = w(ix^D, erad^KSP_ )
              {^C& stateM1%F_low(^C) = w(ix^D, frad^KSP^C_ ) \} 
              {^C& stateM1%vel(^C)   = w(ix^D, u0_ + ^C ) \}         
	          call m1_update_closure_ixD(stateM1,metricM1,ix^D,.true.,get_vel_impl=.true.)
              Gamma_M1(ix^D) = stateM1%Gamma
           {end do^D& \}

 

           ! These are the gauge independent quantities and not the conserved 
              w(ixO^S, erad^KSP_ ) = w(ixO^S, erad^KSP_ ) * metricM1%sqrtg(ixO^S) * metricM1%sqrtg(ixO^S)
              w(ixO^S, nrad^KSP_ ) = w(ixO^S, nrad^KSP_ ) * Gamma_M1(ixO^S) * metricM1%sqrtg(ixO^S) * metricM1%sqrtg(ixO^S)
              w(ixO^S, frad^KSP1_ ) = w(ixO^S, frad^KSP1_ ) * metricM1%sqrtg(ixO^S) * metricM1%sqrtg(ixO^S)
              w(ixO^S, frad^KSP2_ ) = w(ixO^S, frad^KSP2_ ) * metricM1%sqrtg(ixO^S) * metricM1%sqrtg(ixO^S)
              w(ixO^S, frad^KSP3_ ) = w(ixO^S, frad^KSP3_ ) * metricM1%sqrtg(ixO^S) * metricM1%sqrtg(ixO^S)

           ! where ((w(ixO^S, erad^KSP_) < m1_E_atmo * (1.0d0 +1.0d-2)) .or. (w(ixO^S, nrad^KSP_) < m1_E_atmo * (1.0d0 +1.0d-2)))
           !     w(ixO^S, erad^KSP_ ) = m1_E_atmo
           !     w(ixO^S, nrad^KSP_ ) = m1_E_atmo
           !     w(ixO^S, frad^KSP1_ ) = 0.0d0
           !     w(ixO^S, frad^KSP2_ ) = 0.0d0
           !     w(ixO^S, frad^KSP3_ ) = 0.0d0
           ! end where

        \} !end KSP

        
    end subroutine read_NS_3DCarpet_cart

end module NSimport



!========================================================================================
!==========                            initialization                          ==========
!========================================================================================



subroutine initglobaldata_usr
    use mod_eos
    use NSimport
    use mod_xns
    use mod_c2b
    use Init_c2b_cart_3d_interface
    use mod_imhd_con2prim
    use mod_cfc, only : cfc_solver_activate
    use mod_cfc_parameters
    include 'amrvacdef.f'


    if (snapshotini==-1) then
        select case(profile_type)
        case(1) 
            call mod_c2b_read_profile(read_C2B_metric, hydro_type, profile_path_hydro, profile_path_met)
        case(2)
            call mod_XNS_read_profile(profile_path_xns, xns_profile_type) 
        case(3)
            if (eos_type==hybrid) then ! the variable would be different for hybrid case
                variable_name_list(2) = "POST::Wvel[0]"
                variable_name_list(3) = "POST::Wvel[1]"
                variable_name_list(4) = "POST::Wvel[2]"
                variable_fname_list(2) = "Wvel[0]"
                variable_fname_list(3) = "Wvel[1]"
                variable_fname_list(4) = "Wvel[2]"
            endif
            call init_c2b_3d_profile(fpath_c2b_3dcart, variable_fname_list, variable_name_list, num_ref_levels, &
                num_data_blocks, iteration_num, save_cpu_or_mem, init_magnetic_field, &
                init_B_from_vecA, use_eps_init, init_m1_vars)
        case(4)
            write(33,*)"initglobaldata"
        case default
            call mpistop('This profile type does not support yet')
        end select
    endif
    call eos_atmo_activate()
    call eos_initialize_atmo()
    call imhd_con2prim_setup()

end subroutine initglobaldata_usr


subroutine metric_initialize()
    use mod_cfc, only: cfc_metric_init
    use mod_cfc_parameters
    include 'amrvacdef.f'

    if (initialize_metric) then 
        if (.not. useprimitiveRel) call mpistop('CFC initialize must use useprimitiveRel')
        call cfc_metric_init()
    endif

end subroutine metric_initialize


subroutine initonegrid_usr(ixI^L, ixO^L, s)
    use NSimport
    use mod_eos
    use mod_xns
    use mod_c2b
    use Init_c2b_cart_3d_interface
    use mod_metric
    use mod_imhd_intermediate
    include 'amrvacdef.f'

    integer, intent(in)               :: ixI^L, ixO^L
    type(state), intent(inout)        :: s
    ! .. local variables ..
    double precision                  :: rho, press, vel1, vel2, vel3
    integer                           :: ix^D,i,ix_next^D
    double precision                  :: rad, theta
    integer                           :: idir 
    double precision                  :: vD(ixI^S, 1:ndir), v(ixI^S,1:ndir), lfac(ixI^S)
    double precision                  :: gamma(ixI^S,1:3,1:3)
    double precision                  :: pressure(ixI^S),var_eps(ixI^S),var_temp(ixI^S)  ! <- ADD THIS LINE
    integer :: iivar
    logical  :: check_info
    integer :: iw
    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})
        w(ixI^S, :)           = zero
        !w(ixI^S, lfac_)       = 1.0d0
        w(ixI^S, alp_metric_) = 1.0d0
        w(ixI^S, psi_metric_) = 1.0d0
        v(ixI^S, :)           = zero
    
        ! Fill the hydro and metric data at this grid
        select case(profile_type)
        case(1)
            call read_NS_CTWOB(ixI^L,ixO^L,w,x)

            call Eos_update_one_grid(ixI^L, ixO^L, w)
    
            ! first clean the ID 
            call set_atmo(ixI^L,ixO^L,w,x)
    
            if (read_C2B_metric) then
                {do ix^D=ixImin^D, ixImax^D \}
                    theta = x(ix^D, ^Z)
                    call mod_c2b_map_2D(w(ix^D,alp_metric_),x(ix^D,r_),theta,c2b_prof%alp(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    w(ix^D, alp_metric_) = min(w(ix^D, alp_metric_), 1.0d0) 
                    call mod_c2b_map_2D(w(ix^D,psi_metric_),x(ix^D,r_),theta,c2b_prof%psi(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    w(ix^D, psi_metric_) = max(w(ix^D, psi_metric_), 1.0d0) 
                    call mod_c2b_map_2D(w(ix^D,beta_metric1_),x(ix^D,r_),theta,c2b_prof%beta1(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    call mod_c2b_map_2D(w(ix^D,beta_metric2_),x(ix^D,r_),theta,c2b_prof%beta2(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    call mod_c2b_map_2D(w(ix^D,beta_metric3_),x(ix^D,r_),theta,c2b_prof%beta3(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                {end do^D&\}
        
                    ! get the metric
                    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
                    do idir = 1, ndir
                        gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
                    end do
                    ! get W
                    lfac(ixO^S) = 1.0d0
                    lfac(ixO^S) = 1.0d0 - ({^C& gamma(ixO^S,^C,^C) * w(ixO^S, u^C_)**2 +})
                    lfac(ixO^S) = dsqrt( 1.0d0 / lfac(ixO^S) )
                    !w(ixO^S, lfac_) = lfac(ixO^S)
        
                    {^C&w(ixO^S, u^C_)  = w(ixO^S, u^C_) * lfac(ixO^S) \}
        
            endif
    
        case(2)
            call read_NS_XNS_cart(ixI^L, ixO^L, w, x)
    
            ! get the metric
            call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
            do idir = 1, ndir
                gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
            end do
    
            ! get W
            lfac(ixO^S) = 1.0d0
            ! the u^C is only v^C here
            lfac(ixO^S) = 1.0d0 - ({^C& gamma(ixO^S,^C,^C) * w(ixO^S, u^C_)**2 +})
            lfac(ixO^S) = dsqrt( 1.0d0 / lfac(ixO^S) )
            !w(ixO^S, lfac_) = lfac(ixO^S)
    
            {^C&w(ixO^S, u^C_)  = w(ixO^S, u^C_) * lfac(ixO^S) \}
            ! first clean the ID 
            call set_atmo(ixI^L,ixO^L,w,x)
    
        case(3)
            call read_NS_3DCarpet_cart(ixI^L, ixO^L, w, x)
        case(4)
            !!
            w(ixI^S, alp_metric_) = 1.0d0
            w(ixI^S, psi_metric_) = 1.0d0
            w(ixI^S, beta_metric1_) = 0.0d0
            w(ixI^S, beta_metric2_) = 0.0d0
            w(ixI^S, beta_metric3_) = 0.0d0
            w(ixI^S, u1_) = 0.0d0
            w(ixI^S, u2_) = 0.0d0
            w(ixI^S, u3_) = 0.0d0
            w(ixI^S, rho_) = small_rho !1.000000000000000E-010 !1.000000000000000E-004 !0.0d0
            w(ixI^S, T_eps_) = small_temp !1.0d0 !10.0d0
            w(ixI^S, ye_) = big_ye !0.4d0 !1.000000000000000E-001
            if (eos_has_ymu()) w(ixI^S, ymu_) = eos_ymumin
            !w(ixI^S, press_) = 0.0d0
        end select

        call Eos_update_one_grid(ixI^L, ixO^L, w)
        call set_atmo(ixI^L,ixO^L,w,x)

            !TEST
       !     w(ixI^S, alp_metric_) = 1.0d0
       !     w(ixI^S, psi_metric_) = 1.0d0
       !     w(ixI^S, beta_metric1_) = 0.0d0
       !     w(ixI^S, beta_metric2_) = 0.0d0
       !     w(ixI^S, beta_metric3_) = 0.0d0

     
    !   do iivar=1,nw
    !     write(91,*)"check var i",iivar,w(ixO^S, iivar)
    !     if(iivar .eq. 7 .or. iivar .eq. 8 .or. iivar .eq. 9 ) then
    !     write(92,*)"check var i",iivar,w(ixI^S, iivar)
    !     end if  
    !   end do

       ! call check_data_correctness(ixI^L, ixO^L, w, x)
       ! write(*, *) "Congratulation! Passed check_data_correctness in initonegrid_usr!"

       ! 
     !   write(*,*) "Going to debug data correctness"
     !   call debug_check_data_correctness(.true.,check_info)
     !   write(*,*)"check_info",check_info   

! do  iw = 1,nw
!   write(66,*)"iw=",iw,w(ixO^S,iw)
! end do 

    end associate
end subroutine initonegrid_usr


subroutine initonegrid_usr2(ixI^L, ixO^L, s)
    use NSimport
    use mod_eos
    use mod_xns
    use mod_c2b
    use Init_c2b_cart_3d_interface
    use mod_metric
    use mod_imhd_intermediate
    include 'amrvacdef.f'

    integer, intent(in)               :: ixI^L, ixO^L
    type(state), intent(inout)        :: s
    ! .. local variables ..
    double precision                  :: rho, press, vel1, vel2, vel3
    integer                           :: ix^D,i,ix_next^D
    double precision                  :: rad, theta
    integer                           :: idir
    double precision                  :: vD(ixI^S, 1:ndir), v(ixI^S,1:ndir), lfac(ixI^S)
    double precision                  :: gamma(ixI^S,1:3,1:3)
    double precision                  :: r_star, delta(0:2) ! for perturbation
    double precision                  :: pressure(ixI^S),var_eps(ixI^S),var_temp(ixI^S)  ! <- ADD THIS LINE

    logical  :: check_info
    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})

        w(ixI^S, :)           = zero
        !w(ixI^S, lfac_)       = 1.0d0
        w(ixI^S, alp_metric_) = 1.0d0
        w(ixI^S, psi_metric_) = 1.0d0
        v(ixI^S, :)           = zero
    
        ! Fill the hydro and metric data at this grid
        select case(profile_type)
        case(1)
            call read_NS_CTWOB(ixI^L,ixO^L,w,x)

            call Eos_update_one_grid(ixI^L, ixO^L, w)
    
            ! first clean the ID 
            call set_atmo(ixI^L,ixO^L,w,x)
    
            if (read_C2B_metric) then
                {do ix^D=ixImin^D, ixImax^D \}
                    theta = x(ix^D, ^Z)
                    call mod_c2b_map_2D(w(ix^D,alp_metric_),x(ix^D,r_),theta,c2b_prof%alp(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    w(ix^D, alp_metric_) = min(w(ix^D, alp_metric_), 1.0d0) 
                    call mod_c2b_map_2D(w(ix^D,psi_metric_),x(ix^D,r_),theta,c2b_prof%psi(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    w(ix^D, psi_metric_) = max(w(ix^D, psi_metric_), 1.0d0) 
                    call mod_c2b_map_2D(w(ix^D,beta_metric1_),x(ix^D,r_),theta,c2b_prof%beta1(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    call mod_c2b_map_2D(w(ix^D,beta_metric2_),x(ix^D,r_),theta,c2b_prof%beta2(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                    call mod_c2b_map_2D(w(ix^D,beta_metric3_),x(ix^D,r_),theta,c2b_prof%beta3(:,:), &
                                        c2b_prof%radius, c2b_prof%theta, c2b_prof%Nr, c2b_prof%Nth)
                {end do^D&\}

                    ! get the metric
                    call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
                    do idir = 1, ndir
                        gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
                    end do
                    ! get W
                    lfac(ixO^S) = 1.0d0
                    lfac(ixO^S) = 1.0d0 - ({^C& gamma(ixO^S,^C,^C) * w(ixO^S, u^C_)**2 +})
                    lfac(ixO^S) = dsqrt( 1.0d0 / lfac(ixO^S) )
                    !w(ixO^S, lfac_) = lfac(ixO^S)
        
                    {^C&w(ixO^S, u^C_)  = w(ixO^S, u^C_) * lfac(ixO^S) \}
        
            endif
    
        case(2)
            call read_NS_XNS_cart(ixI^L, ixO^L, w, x)
    
            ! get the metric
            call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma(ixI^S,1:3,1:3))
            do idir = 1, ndir
                gamma(ixO^S,idir,idir) = gamma(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
            end do
    
            ! get W
            lfac(ixO^S) = 1.0d0
            ! the u^C is only v^C here
            lfac(ixO^S) = 1.0d0 - ({^C& gamma(ixO^S,^C,^C) * w(ixO^S, u^C_)**2 +})
            lfac(ixO^S) = dsqrt( 1.0d0 / lfac(ixO^S) )
            !w(ixO^S, lfac_) = lfac(ixO^S)
    
            {^C&w(ixO^S, u^C_)  = w(ixO^S, u^C_) * lfac(ixO^S) \}
            ! first clean the ID 
            call set_atmo(ixI^L,ixO^L,w,x)
    
        case(3)
            call read_NS_3DCarpet_cart(ixI^L, ixO^L, w, x)
        case(4)
            !!
            w(ixI^S, alp_metric_) = 1.0d0
            w(ixI^S, psi_metric_) = 1.0d0
            w(ixI^S, beta_metric1_) = 0.0d0
            w(ixI^S, beta_metric2_) = 0.0d0
            w(ixI^S, beta_metric3_) = 0.0d0
            w(ixI^S, u1_) = 0.0d0
            w(ixI^S, u2_) = 0.0d0
            w(ixI^S, u3_) = 0.0d0
            w(ixI^S, rho_) = small_rho !1.000000000000000E-010 !1.000000000000000E-004 !0.0d0
            w(ixI^S, T_eps_) = small_temp !1.0d0 !10.0d0
            w(ixI^S, ye_) = big_ye !0.4d0 !1.000000000000000E-001
            if (eos_has_ymu()) w(ixI^S, ymu_) = eos_ymumin
            !w(ixI^S, press_) = 0.0d0
        end select

        call Eos_update_one_grid(ixI^L, ixO^L, w)
        call set_atmo(ixI^L,ixO^L,w,x)

            !TEST
        !    w(ixI^S, alp_metric_) = 1.0d0
        !    w(ixI^S, psi_metric_) = 1.0d0
        !    w(ixI^S, beta_metric1_) = 0.0d0
        !    w(ixI^S, beta_metric2_) = 0.0d0
        !    w(ixI^S, beta_metric3_) = 0.0d0

       !call check_data_correctness(ixI^L, ixO^L, w, x)
       !write(*, *) "Congratulation! Passed check_data_correctness in initonegrid_usr2!"
        
       ! write(*,*) "Going to debug data correctness"
       ! call debug_check_data_correctness(.true.,check_info)
      !
 ! write(*,*)"check_info",check_info   
 
    end associate
end subroutine initonegrid_usr2

!=============================================================================
subroutine Set_Radiation_Circle(x0,y0,Rcirc,Qems,kappaA,kappaS,kappaN,Qnr)
double precision, intent(out) :: x0,y0,Rcirc,Qems,kappaA,kappaS,kappaN,Qnr

     
     Rcirc = 20.0d0 !0.5d0
     x0 =0.0d0 !1.5d0
     y0 =0.0d0

     Qems = 10.d0 !1.0d+5
     kappaA = 10.0d0 !1.0d+5
     kappaS = 0.0d0 !1.0d2
     kappaN = 0.0d0
     Qnr = 0.0d0
end subroutine Set_Radiation_Circle

 !=============================================================================
! it is inputing ixG, ixM, ps(igrid)
subroutine initonegrid_M1_usr(ixI^L,ixO^L,s,sm1)

  ! initialize one grid within ixO^L
  use mod_oneblock, only: interpolate_oneblock
  use mod_metric, only: lower3
  use mod_eos
  use amrvacpar  
  include 'amrvacdef.f'

  integer, intent(in)              :: ixI^L, ixO^L
  type(m1_impl_state), intent(inout)       :: sm1
  type(state)                      :: s
  ! .. local variables ..
  integer          :: ix^D, i
  double precision :: x0,y0
  double precision :: Rcirc, dist 
  double precision :: Qems, kappaA,kappaS,kappaN,Qnr
  logical :: stopflag
  logical :: set_box = .false.
  logical :: set_circle = .true.
  !-----------------------------------------------------------------------------
   !associate(wradimpl=>sm1%pwrad%wradimpl) ! don't use this 
  associate(x=>s%x%x)
   !sm1%pwrad%wradimpl(ixI^S, :)=0.0d0
     
    ! Rcirc = 0.5d0
    ! x0 = 1.5d0
    ! y0 = 2.0d0

     ! set values:
     !wradimpl(ixI^S,kappa_a^KSP_)=0.0d0 ! don't use this 
  
   {^KSP&
     sm1%pwrad%wradimpl(ixI^S,kappa_a^KSP_)=1.0d-30
     sm1%pwrad%wradimpl(ixI^S,kappa_nr^KSP_)=1.0d-30
     sm1%pwrad%wradimpl(ixI^S,kappa_s^KSP_)=1.0d-30
     sm1%pwrad%wradimpl(ixI^S,Q_er^KSP_)=1.0d-30
     sm1%pwrad%wradimpl(ixI^S,Q_nr^KSP_)=1.0d-30
   \}

  call Set_Radiation_Circle(x0,y0,Rcirc,Qems,kappaA,kappaS,kappaN,Qnr)

  {^KSP&
  {do ix^D = ixImin^D, ixImax^D \}
   ! circle:
   if(set_circle) then
    dist = (x(ix^D,1) - x0)**2 + (x(ix^D,2) - y0)**2
    if(sqrt(dist) .le. Rcirc) then     
       sm1%pwrad%wradimpl(ix^D,Q_nr^KSP_)= Qnr
       sm1%pwrad%wradimpl( ix^D, Q_er^KSP_) = Qems !1.0d0 !1.0d+15
       sm1%pwrad%wradimpl( ix^D, kappa_a^KSP_) = kappaA !1.0d0
       sm1%pwrad%wradimpl( ix^D, kappa_s^KSP_) = kappaS !1.0d0
       sm1%pwrad%wradimpl( ix^D, kappa_nr^KSP_) = kappaN !1.0d0
    end if
   !!if(dist < Rcirc**2) then
   !!  sm1%pwrad%wradimpl( ix^D, kappa_a^KSP_) = 1.0d+10
   !!end if 
   else if(set_box) then
   if((x(ix^D,1).ge. 1.0d0) .and. (x(ix^D,1).le.2.0d0) ) then
    if((x(ix^D,2).ge. 1.5d0) .and. (x(ix^D,2).le. 2.5d0)) then
      sm1%pwrad%wradimpl( ix^D, kappa_a^KSP_) = 1.0d+4
    end if
   end if 
   end if
 
  {end do^D&\} 
  \}

   sm1%pwrad%wradimpl(ixI^S, :)=0.0d0

     stopflag = .false.
  end associate
end subroutine initonegrid_M1_usr
!=========================================

subroutine debug_check_data_correctness(include_ghost, check_info)
    use mod_small_values, only: Simple_check_data_correctness_pt
    include 'amrvacdef.f'
    logical, intent(in)              :: include_ghost
    character(len=*), intent(in)     :: check_info
    ! local
    character(len=256)                 :: new_info
    integer                          :: ix^D, myix_min^D, myix_max^D, iigrid, igrid

if (include_ghost) then 
    {^D& 
        myix_min^D = ixGlo^D
        myix_max^D = ixGhi^D
    }    
else 
    {^D& 
        myix_min^D = ixMlo^D
        myix_max^D = ixMhi^D
    }    
endif

do iigrid = 1, igridstail
    igrid = igrids(iigrid)
    write(new_info, "(a,I,a,a)") "For igrid: ", igrid, ", ", trim(adjustl(check_info))
    associate(x=>px(igrid)%x,w=>pw(igrid)%w)
        {do ix^D=myix_min^D, myix_max^D \}
            call Simple_check_data_correctness_pt(w(ix^D, 1:nw), x(ix^D, 1:ndim), w(ix^D, 1:nw), check_info)
        {end do \}
    end associate
enddo
write(*, *) mype, " pass: ", check_info
end subroutine debug_check_data_correctness
!==================================================

!========================================
subroutine fixer_after_usr2()
  use Init_c2b_cart_3d_interface
  use NSimport

  if (profile_type == 3) call destroy_c2b_3d_profile()

end subroutine fixer_after_usr2

subroutine fix_M1_after_usr2(ixI^L, ixO^L, s)
  use Init_c2b_cart_3d_interface
  use NSimport
  use mod_m1_metric_interface
  use mod_m1_closure

  include 'amrvacdef.f'
  integer, intent(in)               :: ixI^L, ixO^L
  type(state), intent(inout)        :: s

  double precision :: Gamma_M1(ixI^S)
  integer          :: ix^D,i
  
  type(m1_metric_helper) :: metricM1 
  type(m1_closure_helpers) :: stateM1 

  associate(x=>s%x%x,w=>s%w%w)



  if(init_m1_vars)then
             
	       call metricM1%fill_metric(w,x,ixI^L,ixO^L)  
           
           {^KSP&

           !{do ix^D=ixOmin^D, ixOmax^D \}     
           	  !stateM1%E = w(ix^D, erad^KSP_ )
              !{C& stateM1%F_low(^C) = w(ix^D, frad^KSP^C_ ) } 
              !{C& stateM1%vel(^C)   = w(ix^D, u0_ + ^C ) }         
	          !call m1_update_closure_ixD(stateM1,metricM1,ix^D,.true.,get_vel_impl=.true.)
              !Gamma_M1(ix^D) = stateM1%Gamma
           !{end do^D& \}

              !w(ixO^S, erad^KSP_ ) = w(ixO^S, erad^KSP_ ) * metricM1%sqrtg(ixO^S)
              !w(ixO^S, nrad^KSP_ ) = w(ixO^S, nrad^KSP_ ) * Gamma_M1(ixO^S) * metricM1%sqrtg(ixO^S)
              !w(ixO^S, frad^KSP1_ ) = w(ixO^S, frad^KSP1_ ) * metricM1%sqrtg(ixO^S)
              !w(ixO^S, frad^KSP2_ ) = w(ixO^S, frad^KSP2_ ) * metricM1%sqrtg(ixO^S)
              !w(ixO^S, frad^KSP3_ ) = w(ixO^S, frad^KSP3_ ) * metricM1%sqrtg(ixO^S)

              !These are now the conserved qunatities
              w(ixO^S, erad^KSP_ ) = (w(ixO^S, erad^KSP_ ) / metricM1%sqrtg(ixO^S))
              w(ixO^S, nrad^KSP_ ) = (w(ixO^S, nrad^KSP_ ) / metricM1%sqrtg(ixO^S))
              w(ixO^S, frad^KSP1_ ) = w(ixO^S, frad^KSP1_ ) / metricM1%sqrtg(ixO^S)
              w(ixO^S, frad^KSP2_ ) = w(ixO^S, frad^KSP2_ ) / metricM1%sqrtg(ixO^S)
              w(ixO^S, frad^KSP3_ ) = w(ixO^S, frad^KSP3_ ) / metricM1%sqrtg(ixO^S)
            \}
        end if 
        end associate

end subroutine fix_M1_after_usr2

subroutine get_B_field(ixI^L, ixO^L, s)
    use mod_interpolate
    use Init_c2b_cart_3d_interface
    use NSimport
    include 'amrvacdef.f'

    integer, intent(in)               :: ixI^L, ixO^L
    type(state), intent(inout)        :: s
    ! .. local variables ..
    double precision                  :: rho, press, vel1, vel2, vel3
    integer                           :: ix^D,i,ix_next^D, idims, ixCmin^D, ixCmax^D, idir
    double precision                  :: rad, theta, x1, x2, x3
    double precision                  :: vD(ixI^S, 1:ndir), v(ixI^S,1:ndir), lfac(ixI^S)
    double precision                  :: wi(ixI^S,1:nw), xi(ixI^S,1:ndim), Br(ixI^S), Btheta(ixI^S), Bphi(ixI^S), psi_tmp
    double precision                  :: psi6(ixI^S)
    double precision                  :: psi6_surfC^D(ixI^S)

    !-----------------------------------------------------------------------------

    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w},mygeo=>s%geo)
        ! Before this, you need psi
        do idims= 1, ^ND
            xi(ixI^S,1:ndim) = x(ixI^S,1:ndim)
            xi(ixI^S,idims)  = x(ixI^S,idims) + 0.5d0* mygeo%dx(ixI^S,idims)

            call metric_interpolation(ixI^L,ixO^L^LADD3,idims,w,x,wi,xi)
            select case (idims)
            {case (^D)
                psi6_surfC^D(ixI^S) = wi(ixI^S, psi_metric_)**6
                mygeo%surfaceC^D(ixI^S)  = mygeo%surfaceC^D(ixI^S) * psi6_surfC^D(ixI^S)
            \}
            end select
        enddo
        psi6(ixI^S) = w(ixI^S,psi_metric_)**6
        mygeo%dvolume(ixI^S) = mygeo%dvolume(ixI^S) * psi6(ixI^S)
    {#IFDEF STAGGERED 
        ws = 0.0d0
    }
    select case(profile_type)
    case(1)
        call mpistop('No B-field with profile_type = 1')
    case(2)  
        {#IFDEF STAGGERED                                             
            ws = 0.0d0                                                
        }   
        call get_B_field_XNS(ixI^L, ixO^L, s)
    case(3)
        call get_B_field_carpet3D_cart(ixI^L, ixO^L, s)
    case(4)
        do idims=1,ndim
           w(ixI^S, b0_+idims) = 0.0d0
         !  ws(ixI^S, bs0_+idims) = 0.0d0
        end do
    end select
  
        mygeo%dvolume(ixI^S)           = mygeo%dvolume(ixI^S) / psi6(ixI^S)
        {^D&  mygeo%surfaceC^D(ixI^S)  = mygeo%surfaceC^D(ixI^S) / psi6_surfC^D(ixI^S) \}
        {#IFDEF STAGGERED
            call faces2centers(ixO^L, s)

        ! code test
        do idims=1,ndim
              if (w(ix^D, b0_+idims) .ne. w(ix^D, b0_+idims) ) then
                select case(idims)
                {case(^D)
                   !m1_TODO test
                  !write(93,*) 'w(ix^DD, b0_+idims), mygeo%dvolume(ix^DD), w(ix^DD,psi_metric_)'
               !   write(93,*)  "b0+idims, geo_vol, psi"
               !   write(93,*) w(ix^DD, b0_+idims), mygeo%dvolume(ix^DD), w(ix^DD,psi_metric_)
               !   !write(93,*) 'surface'
               !   write(93,*)"surface"
               !   write(93,*) mygeo%surfaceC^D(ixO^S)
               !   write(93,*) 'ws'
               !   write(93,*) "ws"
               !   write(93,*) ws(ixO^S, bs0_+idims)
               !   !call mpistop('nan')
               !    !m1

                  !write(*,*) 'w(ix^DD, b0_+idims), mygeo%dvolume(ix^DD), w(ix^DD,psi_metric_)'
                  !write(*,*) w(ix^DD, b0_+idims), mygeo%dvolume(ix^DD), w(ix^DD,psi_metric_)
                  !write(*,*) 'surface'
                  !write(*,*) mygeo%surfaceC^D(ixO^S)
                  !write(*,*) 'ws'
                  !write(*,*) ws(ixO^S, bs0_+idims)
                  !call mpistop('nan')
                \}
                end select
              endif
        enddo
        }

        !me test
       ! w(ixI^S,b1_) = 0.0d0
       ! w(ixI^S,b2_) = 0.0d0
       ! w(ixI^S,b3_) = 0.0d0

{#IFDEF STAGGERED        
      !  ws(ixI^S,bs1_) = 0.0d0
      ! ws(ixI^S,bs2_) = 0.0d0
      ! ws(ixI^S,bs3_) = 0.0d0
}

    call check_data_correctness(ixI^L, ixO^L, w, x)

    end associate
end subroutine get_B_field

subroutine get_B_field_XNS(ixI^L, ixO^L, s)
    use mod_interpolate
    use mod_XNS
    use NSimport
    include 'amrvacdef.f'

    integer, intent(in)               :: ixI^L, ixO^L
    type(state), intent(inout)        :: s
    ! .. local variables ..
    integer                           :: ix^D, idims, ixCmin^D, ixCmax^D, idir
    double precision                  :: wi(ixI^S,1:nw), xi(ixI^S,1:ndim), Br(ixI^S), Btheta(ixI^S), Bphi(ixI^S), psi_tmp
    double precision                  :: rad, theta, rad_xy, sin_phi, cos_phi, sin_theta, cos_theta, shifted_x, shifted_y

    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w},mygeo=>s%geo)

    if (init_B_from_vecA) then

      {#IFNDEF STAGGERED
        call b_from_vectorpotential(ixI^L,ixO^L,w,x)
      }{#IFDEF STAGGERED
        call b_from_vectorpotential(s%ws%ixG^L,ixI^L,ixO^L,ws,x)
      }
        {#IFDEF STAGGERED
          do idir=1,ndim
            ixCmin^D=ixImin^D;
            ixCmax^D=ixImax^D-kr(idir,^D);
            xi = x
            xi(ixI^S,idir)=xi(ixI^S, idir)+0.5d0*mygeo%dx(ixI^S, idir)
            ! cell-face B_idir
            {do ix^D = ixC^LIM^D \}
               ! get the interface coordinates and use XNS to map
               rad    = dsqrt( sum(xi(ix^D,1:ndim)**2) )
               rad_xy = dsqrt( sum(xi(ix^D,1:2)**2) )
               ! this is bad
               !theta  = dmod(2.0d0*dpi + datan2(x(ix^D,1),x(ix^D,2)), 2.0d0*dpi)  ! this is not good
               theta  = dmod(2.0d0*dpi + datan2(rad_xy,x(ix^D,3)),    2.0d0*dpi) ! this is okay

               !theta = acos( xi(ix^D,3) / rad )
               call mod_XNS_map_2D(psi_tmp,rad, theta,prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               psi_tmp = max(psi_tmp,1.0d0)
               ws(ix^D, idir) = ws(ix^D, idir) * psi_tmp**6
            {end do^D&\}
          enddo
        } 
    else
        {#IFNDEF STAGGERED
        {do ix^D=ixOmin^D, ixOmax^D \}
          shifted_x = x(ix^D, 1) - shift_origin
          shifted_y = x(ix^D, 2) - shift_origin
          rad = dsqrt(x(ix^D,3)**2 + shifted_x**2 + shifted_y**2 )
          rad_xy = dsqrt( shifted_x**2 + shifted_y**2 )
          theta  = dmod(2.0d0*dpi + datan2(rad_xy,x(ix^D,3)),    2.0d0*dpi)
          !theta = acos( x(ix^D,3) / rad )
          sin_phi = shifted_y / rad_xy
          cos_phi = shifted_x / rad_xy
          sin_theta = rad_xy / rad
          !assume z_symm
          if ( theta > 0.5d0 * dpi ) theta = dpi - theta
          call mod_XNS_map_2D(Br(ix^D),rad, theta,prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
          call mod_XNS_map_2D(Btheta(ix^D),rad, theta,prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
          call mod_XNS_map_2D(Bphi(ix^D),rad, theta,prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)

          Btheta(ix^D) = Btheta(ix^D) * rad
          Bphi(ix^D)   = Bphi(ix^D)   * rad_xy

          sin_theta = rad_xy     / rad
          cos_theta = x(ix^D, 3) / rad
          sin_phi   = x(ix^D, 2) / rad_xy
          cos_phi   = x(ix^D, 1) / rad_xy

          w(ix^D, b1_) = Br(ix^D) * sin_theta * cos_phi + Btheta(ix^D) * cos_theta * cos_phi - Bphi(ix^D) * sin_phi
          w(ix^D, b2_) = Br(ix^D) * sin_theta * sin_phi + Btheta(ix^D) * cos_theta * sin_phi + Bphi(ix^D) * cos_phi
          w(ix^D, b3_) = Br(ix^D) * cos_theta           - Btheta(ix^D) * sin_theta
        {end do \}

        }
 
        {#IFDEF STAGGERED
          ws = 0.0d0
          do idir=1,ndim
            ixCmin^D=ixImin^D;
            ixCmax^D=ixImax^D-kr(idir,^D);
            xi = x
            xi(ixI^S,idir)=xi(ixI^S, idir)+0.5d0*mygeo%dx(ixI^S, idir)

            ! cell-face B_idir
            {do ix^D = ixC^LIM^D \}
               ! get the interface coordinates and use XNS to map
               rad    = dsqrt( sum(xi(ix^D,1:ndim)**2) )
               rad_xy = dsqrt( sum(xi(ix^D,1:2)**2) )
               ! this is bad
               !theta  = dmod(2.0d0*dpi + datan2(x(ix^D,1),x(ix^D,2)), 2.0d0*dpi)  ! this is not good
               theta  = dmod(2.0d0*dpi + datan2(rad_xy,xi(ix^D,3)),    2.0d0*dpi) ! this is okay

               !theta  = dmod(2.0d0*dpi + datan2(rad_xy,xi(ix^D,3)),    2.0d0*dpi)
               !theta = acos( xi(ix^D,3) / rad )

               call mod_XNS_map_2D(psi_tmp,rad, theta,prof%psi(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               call mod_XNS_map_2D(Br(ix^D),rad, theta,prof%B1(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               call mod_XNS_map_2D(Btheta(ix^D),rad, theta,prof%B2(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               call mod_XNS_map_2D(Bphi(ix^D),rad, theta,prof%B3(:,:),prof%radius,prof%theta,prof%Nr,prof%Nth)
               psi_tmp = max(psi_tmp,1.0d0)

               !Br(ix^D)     = Br(ix^D)
               Btheta(ix^D) = Btheta(ix^D) * rad
               Bphi(ix^D)   = Bphi(ix^D)   * rad_xy

               sin_theta = rad_xy      / rad
               cos_theta = xi(ix^D, 3) / rad
               sin_phi   = xi(ix^D, 2) / rad_xy
               cos_phi   = xi(ix^D, 1) / rad_xy

               select case (idir)
               case (1)
                  ws(ix^D,idir) = Br(ix^D) * sin_theta * cos_phi + Btheta(ix^D) * cos_theta * cos_phi - Bphi(ix^D) * sin_phi
               case (2)
                  ws(ix^D,idir) = Br(ix^D) * sin_theta * sin_phi + Btheta(ix^D) * cos_theta * sin_phi + Bphi(ix^D) * cos_phi
               case (3)
                  ws(ix^D,idir) = Br(ix^D) * cos_theta           - Btheta(ix^D) * sin_theta
               end select
               ws(ix^D, idir) = ws(ix^D, idir) * psi_tmp**6
            {end do^D&\}
          end do
        }
    endif
    end associate
end subroutine get_B_field_XNS

subroutine get_B_field_carpet3D_cart(ixI^L, ixO^L, s)
    use mod_interpolate
    use mod_XNS
    use mod_c2b_cart3d, only: load_selected_ref_and_datablock
    use Init_c2b_cart_3d_interface
    use NSimport
    include 'amrvacdef.f'

    integer, intent(in)               :: ixI^L, ixO^L
    type(state), intent(inout)        :: s
    ! .. local variables ..
    integer                           :: ix^D, idims, ixCmin^D, ixCmax^D, idir
    double precision                  :: wi(ixI^S,1:nw), xi(ixI^S,1:ndim)
    double precision                  :: shifted_x, shifted_y, x1, x2, x3, psi_tmp
    integer                           :: mask_array(0:num_ref_levels-1, 0:num_data_blocks-1)

    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w},mygeo=>s%geo)
    call detect_level_and_blocks_to_be_loaded(ixI^L, x, my_c2b_3d_profile%rho, num_ref_levels, num_data_blocks, mask_array, mirror_zplane)
    call load_selected_ref_and_datablock(my_c2b_3d_profile%psi, num_ref_levels, num_data_blocks, mask_array)
    ! sanity check
    if (.not. allocated(my_c2b_3d_profile%psi%all_refine_levels)) then
        call mpistop("Error, psi should be initialized, but it is not initialized or already &
            been destroyed in subroutine get_B_field_carpet3D_cart!")
    endif

    if (init_B_from_vecA) then
        {#IFDEF STAGGERED
          call b_from_vectorpotential(s%ws%ixG^L,ixI^L,ixO^L,ws,x)
          do idir=1,ndim
            ixCmin^D=ixImin^D;
            ixCmax^D=ixImax^D-kr(idir,^D);
            xi = x
            xi(ixI^S,idir)=xi(ixI^S, idir)+0.5d0*mygeo%dx(ixI^S, idir)
            ! cell-face B_idir
            {do ix^D = ixC^LIM^D \}
                x1 = xi(ix^D, 1) - shift_origin
                x2 = xi(ix^D, 2) - shift_origin
                if (mirror_zplane) then
                    x3 = abs(xi(ix^D, 3)) ! z-symmetry
                else
                    x3 = xi(ix^D, 3) ! z-symmetry
                endif
                call interp_cart3d_carpet_variable_one_point(psi_tmp, x1, x2, x3, my_c2b_3d_profile%psi, interp_metric, 1.0d0)
                psi_tmp = dsqrt(1.0d0/psi_tmp)
                psi_tmp = max(psi_tmp,1.0d0)
                ws(ix^D, idir) = ws(ix^D, idir) * psi_tmp**6
            {end do^D&\}
          enddo
        }
        {#IFNDEF STAGGERED
            call b_from_vectorpotential(ixI^L,ixO^L,w,x)
        }
    else
        call load_selected_ref_and_datablock(my_c2b_3d_profile%Bx, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%By, num_ref_levels, num_data_blocks, mask_array)
        call load_selected_ref_and_datablock(my_c2b_3d_profile%Bz, num_ref_levels, num_data_blocks, mask_array)
        ! sanity check
        if (.not. (allocated(my_c2b_3d_profile%Bx%all_refine_levels) &
            .and. allocated(my_c2b_3d_profile%By%all_refine_levels) .and. allocated(my_c2b_3d_profile%Bz%all_refine_levels))) then
            call mpistop("Error, Bx, By, and Bz should all be initialized, but one of them is not initialized or already &
                been destroyed in subroutine get_B_field_carpet3D_cart!")
        endif
        {#IFNDEF STAGGERED
            {do ix^D=ixOmin^D, ixOmax^D \}
                x1 = x(ix^D, 1) - shift_origin
                x2 = x(ix^D, 2) - shift_origin
                if (mirror_zplane) then
                    x3 = abs(x(ix^D, 3)) ! z-symmetry
                else
                    x3 = x(ix^D, 3) ! z-symmetry
                endif

                call interp_cart3d_carpet_variable_one_point(w(ix^D, b1_), x1, x2, x3, my_c2b_3d_profile%Bx, interp_fluid, 0.0d0)
                call interp_cart3d_carpet_variable_one_point(w(ix^D, b2_), x1, x2, x3, my_c2b_3d_profile%By, interp_fluid, 0.0d0)
                call interp_cart3d_carpet_variable_one_point(w(ix^D, b3_), x1, x2, x3, my_c2b_3d_profile%Bz, interp_fluid, 0.0d0)

                if (mirror_zplane .and. x(ix^D, 3)<0.d0) then
                    w(ix^D, b1_) = -w(ix^D, b1_)
                    w(ix^D, b2_) = -w(ix^D, b2_)
                endif
                ! Note that the EM unit in FIL and BHAC are different, this lead to the factor dsqrt(4*dpi)
                {^C& w(ix^D, b^C_) = w(ix^D, b^C_)/dsqrt(4.0d0*dpi)\}
            {end do \}
        }
 
        {#IFDEF STAGGERED
            ws = 0.0d0
            do idir=1,ndim
                ixCmin^D=ixImin^D;
                ixCmax^D=ixImax^D-kr(idir,^D);
                xi = x
                xi(ixI^S,idir)=xi(ixI^S, idir)+0.5d0*mygeo%dx(ixI^S, idir)
    
                ! cell-face B_idir
                {do ix^D = ixC^LIM^D \}
                    x1 = xi(ix^D, 1) - shift_origin
                    x2 = xi(ix^D, 2) - shift_origin
                    if (mirror_zplane) then
                        x3 = abs(xi(ix^D, 3)) ! z-symmetry
                    else
                        x3 = xi(ix^D, 3) ! z-symmetry
                    endif
                    call interp_cart3d_carpet_variable_one_point(psi_tmp, x1, x2, x3, my_c2b_3d_profile%psi, interp_metric, 1.0d0)
                    ! contradictory to intuition here, Bx, and By should be inversed in the mirror plane of z, but Bz don't need
                    select case (idir)
                    case (1)
                        call interp_cart3d_carpet_variable_one_point(ws(ix^D,idir), x1, x2, x3, my_c2b_3d_profile%Bx, interp_fluid, 0.0d0)
                        if (mirror_zplane .and. xi(ix^D, 3)<0.d0) ws(ix^D, idir) = -ws(ix^D, idir)
                    case (2)
                        call interp_cart3d_carpet_variable_one_point(ws(ix^D,idir), x1, x2, x3, my_c2b_3d_profile%By, interp_fluid, 0.0d0)
                        if (mirror_zplane .and. xi(ix^D, 3)<0.d0) ws(ix^D, idir) = -ws(ix^D, idir)
                    case (3)
                        call interp_cart3d_carpet_variable_one_point(ws(ix^D,idir), x1, x2, x3, my_c2b_3d_profile%Bz, interp_fluid, 0.0d0)
                    end select
                    psi_tmp = dsqrt(1.0d0/psi_tmp)
                    psi_tmp = max(psi_tmp,1.0d0)
                    ! Note that the EM unit in FIL and BHAC are different, this lead to the factor dsqrt(4*dpi)
                    ws(ix^D, idir) = (ws(ix^D, idir)*psi_tmp**6)/dsqrt(4.0d0*dpi)
                {end do^D&\}
            end do
        }
    endif
    end associate
end subroutine get_B_field_carpet3D_cart

subroutine initvecpot_usr(ixI^L, ixO^L, x, A, idir)
    ! initialize the vectorpotential on the corners
    ! used by b_from_vectorpotential()

    ! use mod_metric, only: CoordToBL
    use mod_interpolate
    use mod_XNS
    use mod_eos
    use NSimport
    use Init_c2b_cart_3d_interface
    use mod_c2b_cart3d, only: load_selected_ref_and_datablock
    include 'amrvacdef.f'

    integer, intent(in)                :: ixI^L, ixO^L, idir
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    ! local
    double precision                   :: rad, theta, rad_xy_sqr, rad_xy
    integer                            :: ix^D, i
    double precision                   :: x1, x2, x3
    integer                            :: mask_array(0:num_ref_levels-1, 0:num_data_blocks-1)
    !-----------------------------------------------------------------------------
        A = 0.d0
        select case (profile_type)
        case (1)
            if (idir==3) return ! only A_phi exist
            {do ix^D = ixO^LIM^D \}
               rad_xy_sqr = sum(x(ix^D,1:2)**2)
               if ( rad_xy_sqr < tiny(0.0d0) ) cycle
               rad   = dsqrt( sum(x(ix^D,1:ndim)**2) )
               rad_xy = dsqrt( x(ix^D,1)**2 + x(ix^D,2)**2 )
               !theta  = dmod(2.0d0*dpi + datan2(x(ix^D,1),x(ix^D,2)), 2.0d0*dpi)   
               !theta  = dmod(2.0d0*dpi + datan2(rad_xy,x(ix^D,3)),    2.0d0*dpi)  
               theta = acos( x(ix^D,3) / rad )
        
               call mod_XNS_map_2D(A(ix^D),rad, theta,prof%A3,prof%radius,prof%theta,prof%Nr,prof%Nth)
        
               select case ( idir )
               case (1)
                  A(ix^D) = - A(ix^D) * x(ix^D, 2) / rad_xy_sqr
               case (2)
                  A(ix^D) =   A(ix^D) * x(ix^D, 1) / rad_xy_sqr
               end select
            {end do^D&\}
        case (3)
            call detect_level_and_blocks_to_be_loaded(ixI^L, x, my_c2b_3d_profile%rho, num_ref_levels, num_data_blocks, mask_array, mirror_zplane)
            if (idir==1) then
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Ax, num_ref_levels, num_data_blocks, mask_array)
            else if (idir==2) then
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Ay, num_ref_levels, num_data_blocks, mask_array)
            else
                call load_selected_ref_and_datablock(my_c2b_3d_profile%Az, num_ref_levels, num_data_blocks, mask_array)
            endif
            ! sanity check
            if (.not. (allocated(my_c2b_3d_profile%Ax%all_refine_levels) &
                .and. allocated(my_c2b_3d_profile%Ay%all_refine_levels) .and. allocated(my_c2b_3d_profile%Az%all_refine_levels))) then
                call mpistop("Error, Ax, Ay, and Az should all be initialized, but one of them is not initialized or already &
                    been destroyed in subroutine get_B_field_carpet3D_cart!")
            endif
            {do ix^D = ixO^LIM^D  \}
                x1 = x(ix^D, 1)
                x2 = x(ix^D, 2)
                if (mirror_zplane) then
                    x3 = abs(x(ix^D, 3)) ! z-symmetry
                else
                    x3 = x(ix^D, 3) ! z-symmetry
                endif
                if (idir==1) then
                    call interp_cart3d_carpet_variable_one_point(A(ix^D), x1, x2, x3, my_c2b_3d_profile%Ax, interp_vectA, 0.d0)
                else if (idir==2) then
                    call interp_cart3d_carpet_variable_one_point(A(ix^D), x1, x2, x3, my_c2b_3d_profile%Ay, interp_vectA, 0.d0)
                else if (idir==3) then
                    call interp_cart3d_carpet_variable_one_point(A(ix^D), x1, x2, x3, my_c2b_3d_profile%Az, interp_vectA, 0.d0)
                    if (mirror_zplane .and. x(ix^D, 3)<0.d0) A(ix^D) = -A(ix^D)
                else
                    write(*, *) "Error of dimension when calling initvecpot_usr, please reset parameter idir"
                    A(ix^D) = 0.d0
                endif
                ! Note that the EM unit in FIL and BHAC are different, this lead to the following relation
                A(ix^D) = A(ix^D)/dsqrt(4.0d0*dpi)
            {end do^D& \}
        case(4) 
            A(ixI^S) = 0.0d0
        case default
            write(*, *) "Initialize B field from vector potential not implemented yet in this profile_type", profile_type
        end select

end subroutine initvecpot_usr
!=============================================================================

subroutine set_conserve()
  include 'amrvacdef.f'
  ! .. local ..
  integer                   :: iigrid, igrid, ix^D
  !-----------------------------------------------------------------------------
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     associate(x=>ps(igrid)%x%x,w=>ps(igrid)%w%w{#IFDEF STAGGERED ,ws=>ps(igrid)%ws%w})
     call set_tmpGlobals(igrid)
     ! Fixme: after metric_initialization, psi6 changed, but psi6 inside ws was old!, how to deal with?
     ! {#IFDEF STAGGERED
     !call faces2centers(ixM^LL,ps(igrid))
     ! }
     call conserve(ixG^LL,ixM^LL,w,x,patchfalse)
     end associate
  end do
end subroutine set_conserve
  !=============================================================================

subroutine check_data_correctness(ixI^L, ixO^L, w, x)
    use mod_eos
    include 'amrvacdef.f'
    integer, intent(in)              :: ixI^L, ixO^L
    double precision, intent(in)     :: x(ixI^S, 1:ndim)
    double precision, intent(in)     :: w(ixI^S, 1:nw)
    ! local variables
    integer                          :: ix^D, i, idim
    logical                          :: stopflag
  
    stopflag = .false. 
  
    ! EOS check
    {do ix^D=ixOmin^D, ixOmax^D \}
        if (eos_uses_ye()) then
            if (w(ix^D, ye_)  .lt. eos_yemin) then
                write(*,*) w(ix^D, ye_), eos_yemin, 'ye'
                call mpistop("ye < eosmin in check correctness")
            endif
            if (w(ix^D, rho_) .lt. eos_rhomin)  then
                write(*,*)  w(ix^D, rho_), eos_rhomin, 'rho'
                call mpistop('rho < eosmin in data check correctness')
            endif
            if (w(ix^D, T_eps_) .lt. eos_tempmin) then
                write(*,*) w(ix^D, T_eps_), 'temp'
                call mpistop('temp < eosmin in data check correctness')
            endif
            if (w(ix^D, ye_)   .gt. eos_yemax)   call mpistop ("ye > eosmax in check correctness"   ) 
            if (eos_has_ymu()) then
                if (w(ix^D, ymu_) .lt. eos_ymumin) call mpistop ("ymu < eosmin in check correctness")
                if (w(ix^D, ymu_) .gt. eos_ymumax) call mpistop ("ymu > eosmax in check correctness")
            endif
            if (w(ix^D, rho_)  .gt. eos_rhomax)  call mpistop ("rho > eosmax in check correctness"  )
            if (w(ix^D, T_eps_) .gt. eos_tempmax) call mpistop ("temp > eosmax in check correctness" )
        else
            !nothing to do
        endif
    {end do \}

    {do ix^D=ixOmin^D, ixOmax^D \}
        ! NaN check
        do i = 1, nw
            if (w(ix^D, i) .ne. w(ix^D, i)) then
                write(*,*) w(ix^D, i), " appeared in check correctness at position: ", x(ix^D,:), " in variable ", i

                stopflag = .true.
            endif
        enddo



        ! deal with stop
        if (stopflag) then
            write(*,*) ' x(ix^D,:)', x(ix^D,:)
            call mpistop('Sth wrong in data check correctness')
        endif

        ! Lorentz factor check
       ! if (w(ix^D,lfac_) .lt. 1.0d0 ) then
       !     write(*,*) w(ix^D,lfac_), 'lfac < 1 in check correctness, at position: ', x(ix^D,:)
       !     call mpistop('Sth wrong in data check correctness')
       ! endif

        ! metric check
        {#IFDEF DY_SP
            if (w(ix^D, alp_metric_) > 1.0d0 ) then
                write(*,*) w(ix^D, alp_metric_), 'alp_metric > 1 in check correctness, at position: ', x(ix^D,:)
                call mpistop('Sth wrong in data check correctness')
            endif
            if (w(ix^D, alp_metric_) <= 0.0d0 ) then
                write(*,*) w(ix^D, alp_metric_), 'alp_metric <= 0 in check correctness, at position: ', x(ix^D,:)
                call mpistop('Sth wrong in data check correctness')
            endif
            if ((w(ix^D, psi_metric_) < 1.0d0) .or. (w(ix^D, psi_metric_) > 3.0d0)) then
                write(*,*) w(ix^D, psi_metric_), 'psi_metric < 1 or > 3 in check correctness, at position: ', x(ix^D,:)
                call mpistop('Sth wrong in data check correctness')
            endif
        }
  
    {end do \}
end subroutine check_data_correctness



!========================================================================================
!==========                         personalized output                        ==========
!========================================================================================

subroutine specialvar_output(ixI^L,ixO^L,nwmax,w,s,normconv,dxgrid,level,sold,sm1)

   ! this subroutine can be used in convert, to add auxiliary variables to the
   ! converted output file, for further analysis using tecplot, paraview, ....
   ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
   !
   ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
   ! corresponding normalization values (default value 1)

   use mod_physaux, only: get_b2, get_u4, get_b4
   use mod_metric, only: lower3, lower4   !, u3CoordToCart
   use mod_transform, only: matvec, boost
   use mod_imhd_intermediate
    use mod_m1_metric_interface

   use mod_m1_eas
   use mod_eos

   include 'amrvacdef.f'

   integer, intent(in)                :: ixI^L,ixO^L,nwmax
   double precision, intent(inout)    :: w(ixI^S,nwmax)
   type(state)                        :: s
   type(m1_impl_state), intent(in)       :: sm1
   type(state), intent(in)            :: sold
   double precision                   :: normconv(0:nwmax)
   double precision, intent(in)       :: dxgrid(ixI^S,1:ndim)
   integer, intent(in)                :: level
   ! .. local ..
   integer                            :: ix^D, idir
   double precision, dimension(ixI^S) :: b2,lfact
   double precision, dimension(ixI^S,0:ndir) :: u4, u4d, b4, Lmom, FluxSum
   double precision, dimension(ixG^T)        :: rhoh
   double precision, dimension(ixI^S, 1:ndim) :: xCart
   double precision, dimension(ixI^S, 1:ndir) :: uCoord
   double precision, dimension(ixI^S, 1:ndir) :: betaCoord
   double precision, dimension(ixI^S, 1:ndir) :: betaCart
   double precision, dimension(ixI^S)         :: vR, Omega, dwdt
   double precision, dimension(ixImin1:ixImax1,1:3,&
      1:3)         :: Jac_dummy
   double precision, dimension(ixI^S, 1:ndir) :: uCart
   double precision                   ::x_carpet1, x_carpet2, x_carpet3
   double precision :: wradimpl1(ixI^S,1:nm1rad_eas)   
   double precision :: wradimpl2(ixI^S,1:nm1rad_eas)
   double precision :: tmp1 = 1.0d0
   double precision, dimension(ixI^S,1:3,1:3) :: gamma_ij
   type(m1_metric_helper) :: metricM1 	
  type(m1_closure_helpers) :: stateM1 

   logical :: calc_eqRates = .true.
   integer :: igrid1 = 1
   double precision, dimension(ixI^S) :: J_out, Gamma_out
   double precision, dimension(ixI^S) :: Gamma1, Gamma2, Gamma3
   integer :: ispec
   
   double precision, dimension(ixI^S) :: eos_ent, eos_prs
   double precision, dimension(ixI^S) :: eos_rho, eos_eps, eos_temp, eos_ye, eos_ymu
   !-----------------------------------------------------------------------------
   associate(x=>s%x%x{#IFDEF STAGGERED ,ws=>s%ws%w},wold=>sold%w%w)   

	  call metricM1%fill_metric(w,x,ixI^L,ixO^L)  

   call imhd_get_intermediate_variables(ixI^L, ixO^L, w(ixI^S,1:nw), x, &
                                        b2=b2(ixI^S),lfac=lfact(ixI^S))

   if (nwmax-nw .gt. 0) then
      w(ixO^S,nw+1) = b2(ixO^S)
   end if

   if (nwmax-nw .gt. 1) then
      {#IFDEF STAGGERED
      call div_staggered(ixI^L,ixO^L,s,w(ixO^S,nw+2))
      }{#IFNDEF STAGGERED
      ! Reduce output array size, +1 was added for eventual pointdata output
      call get_divb(ixI^L,ixO^L^LSUB1,w(ixI^S,1:nw),w(ixI^S,nw+2))
      }
   end if

   if (nwmax-nw .gt. 2) &
       w(ixI^S,nw+3) = dble(level)

   if (nwmax-nw .gt. 3) &
       w(ixI^S,nw+4) = 0.0d0

    call m1_add_eas(tmp1,t,igrid1,x,w,wradimpl1,ixI^L,calc_eqRates)

    if (nwmax-nw .gt. 4) &
       w(ixI^S,nw+5) = wradimpl1(ixI^S,Q_er1_)

     if (nwmax-nw .gt. 5) &
       w(ixI^S,nw+6) = wradimpl1(ixI^S,kappa_a1_)

     if (nwmax-nw .gt. 6) &
       w(ixI^S,nw+7) = wradimpl1(ixI^S,kappa_s1_)

     if (nwmax-nw .gt. 7) &
       w(ixI^S,nw+8) = wradimpl1(ixI^S,kappa_nr1_)

        if (nwmax-nw .gt. 8) &
       w(ixI^S,nw+9) = wradimpl1(ixI^S,Q_nr1_)

     if (nwmax-nw .gt. 9) &
       w(ixI^S,nw+10) = wradimpl1(ixI^S,Q_er2_)

     if (nwmax-nw .gt. 10) &
       w(ixI^S,nw+11) = wradimpl1(ixI^S,kappa_a2_)

       if (nwmax-nw .gt. 11) &
       w(ixI^S,nw+12) = wradimpl1(ixI^S,kappa_s2_)

       if (nwmax-nw .gt. 12) &
       w(ixI^S,nw+13) = wradimpl1(ixI^S,kappa_nr2_)

       if (nwmax-nw .gt. 13) &
       w(ixI^S,nw+14) = wradimpl1(ixI^S,Q_nr2_)

       if (nwmax-nw .gt. 14) &
       w(ixI^S,nw+15) = wradimpl1(ixI^S,Q_er3_)

     if (nwmax-nw .gt. 15) &
       w(ixI^S,nw+16) = wradimpl1(ixI^S,kappa_a3_)

       if (nwmax-nw .gt. 16) &
       w(ixI^S,nw+17) = wradimpl1(ixI^S,kappa_s3_)

       if (nwmax-nw .gt. 17) &
       w(ixI^S,nw+18) = wradimpl1(ixI^S,kappa_nr3_)

       if (nwmax-nw .gt. 18) &
       w(ixI^S,nw+19) = wradimpl1(ixI^S,Q_nr3_)

       if (nwmax-nw .gt. 19) &
       w(ixI^S,nw+20) = lfact(ixI^S) 
       
      {^C& FluxSum(ixI^S,^C) = 0.0d0 \}
      {^KSP& {^C& FluxSum(ixI^S,^C) = FluxSum(ixI^S,^C) + w(ixI^S,frad^KSP^C_) \} \}

     Lmom(ixI^S,1) = x(ixI^S,2) * FluxSum(ixI^S,3) - FluxSum(ixI^S,2) * x(ixI^S,3)
     Lmom(ixI^S,2) = - (x(ixI^S,1) * FluxSum(ixI^S,3) - FluxSum(ixI^S,1) * x(ixI^S,1)  )
     Lmom(ixI^S,3) = x(ixI^S,1) * FluxSum(ixI^S,2) - FluxSum(ixI^S,1) *  x(ixI^S,2)

       if(nwmax-nw .gt. 20) &
       w(ixI^S,nw+21) = Lmom(ixI^S,1)

       if(nwmax-nw .gt. 21) &
       w(ixI^S,nw+22) = Lmom(ixI^S,2)

       if(nwmax-nw .gt. 22) &
       w(ixI^S,nw+23) = Lmom(ixI^S,3)

       ! Get the 3-metric components
   call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, gamma_ij(ixI^S,1:3,1:3))
   
   ! Include conformal factor
   do idir = 1, ndir
       gamma_ij(ixO^S,idir,idir) = gamma_ij(ixO^S,idir,idir) * w(ixO^S, psi_metric_)**4
   end do
   
   if (nwmax-nw .gt. 23) &
      w(ixI^S,nw+24) = gamma_ij(ixI^S,1,1)  ! gxx
   
   if (nwmax-nw .gt. 24) &
      w(ixI^S,nw+25) = gamma_ij(ixI^S,2,2)  ! gyy
   
   if (nwmax-nw .gt. 25) &
      w(ixI^S,nw+26) = gamma_ij(ixI^S,3,3)  ! gzz
      
   if (nwmax-nw .gt. 26) &
      w(ixI^S,nw+27) = gamma_ij(ixI^S,1,2)  ! gxy
    
    {do ix^D=ixOmin^D, ixOmax^D \}     
           	  stateM1%E = w(ix^D, erad1_ ) / metricM1%sqrtg(ix^D)
              {^C& stateM1%F_low(^C) = w(ix^D, frad1^C_ ) / metricM1%sqrtg(ix^D) \} 
              {^C& stateM1%vel(^C)   = w(ix^D, u0_ + ^C ) \}         
	          call m1_update_closure_ixD(stateM1,metricM1,ix^D,.true.,get_vel_impl=.true.)
              Gamma1(ix^D) = stateM1%Gamma
              J_out(ix^D) = stateM1%J
    {end do^D& \}

   ! Compute J and Gamma for each species
   !call m1_update_closure(metricM1, w, x, ixI^L, ixO^L, &
   !                       species=1, get_zeta=.true., &
   !                       Gamma=Gamma1, Jrad=J_out)
   
   if (nwmax-nw .gt. 27) &
      w(ixI^S,nw+28) = J_out(ixI^S)  ! J1
   
   if (nwmax-nw .gt. 28) &
      w(ixI^S,nw+29) = Gamma1(ixI^S)  ! Gamma1


   {do ix^D=ixOmin^D, ixOmax^D \}     
           	  stateM1%E = w(ix^D, erad2_ ) / metricM1%sqrtg(ix^D)
              {^C& stateM1%F_low(^C) = w(ix^D, frad2^C_ ) / metricM1%sqrtg(ix^D) \} 
              {^C& stateM1%vel(^C)   = w(ix^D, u0_ + ^C ) \}         
	          call m1_update_closure_ixD(stateM1,metricM1,ix^D,.true.,get_vel_impl=.true.)
              Gamma2(ix^D) = stateM1%Gamma
              J_out(ix^D) = stateM1%J
    {end do^D& \}
   !call m1_update_closure(metricM1, w, x, ixI^L, ixO^L, &
   !                       species=2, get_zeta=.true., &
   !                       Gamma=Gamma2, Jrad=J_out)
   
   if (nwmax-nw .gt. 29) &
      w(ixI^S,nw+30) = J_out(ixI^S)  ! J2
   
   if (nwmax-nw .gt. 30) &
      w(ixI^S,nw+31) = Gamma2(ixI^S)  ! Gamma2

   {do ix^D=ixOmin^D, ixOmax^D \}     
           	  stateM1%E = w(ix^D, erad3_ ) / metricM1%sqrtg(ix^D)
              {^C& stateM1%F_low(^C) = w(ix^D, frad3^C_ ) / metricM1%sqrtg(ix^D) \} 
              {^C& stateM1%vel(^C)   = w(ix^D, u0_ + ^C ) \}         
	          call m1_update_closure_ixD(stateM1,metricM1,ix^D,.true.,get_vel_impl=.true.)
              Gamma3(ix^D) = stateM1%Gamma
              J_out(ix^D) = stateM1%J
    {end do^D& \}
   
   !call m1_update_closure(metricM1, w, x, ixI^L, ixO^L, &
   !                       species=3, get_zeta=.true., &
   !                       Gamma=Gamma3, Jrad=J_out)
   
   if (nwmax-nw .gt. 31) &
      w(ixI^S,nw+32) = J_out(ixI^S)  ! J3
   
   if (nwmax-nw .gt. 32) &
      w(ixI^S,nw+33) = Gamma3(ixI^S)  ! Gamma3

   ! Output sqrtg
   if (nwmax-nw .gt. 33) &
      w(ixI^S,nw+34) = metricM1%sqrtg(ixI^S)  ! sqrtg

   ! Species 1 primitives
   if (nwmax-nw .gt. 34) &
      w(ixI^S,nw+35) = w(ixI^S, nrad1_) / (metricM1%sqrtg(ixI^S) * Gamma1(ixI^S))  ! nrad1_prim
   
   if (nwmax-nw .gt. 35) &
      w(ixI^S,nw+36) = w(ixI^S, erad1_) / metricM1%sqrtg(ixI^S)  ! erad1_prim
   
   ! Species 2 primitives
   if (nwmax-nw .gt. 36) &
      w(ixI^S,nw+37) = w(ixI^S, nrad2_) / (metricM1%sqrtg(ixI^S) * Gamma2(ixI^S))  ! nrad2_prim
   
   if (nwmax-nw .gt. 37) &
      w(ixI^S,nw+38) = w(ixI^S, erad2_) / metricM1%sqrtg(ixI^S)  ! erad2_prim
   
   ! Species 3 primitives
   if (nwmax-nw .gt. 38) &
      w(ixI^S,nw+39) = w(ixI^S, nrad3_) / (metricM1%sqrtg(ixI^S) * Gamma3(ixI^S))  ! nrad3_prim
   
   if (nwmax-nw .gt. 39) &
      w(ixI^S,nw+40) = w(ixI^S, erad3_) / metricM1%sqrtg(ixI^S)  ! erad3_prim

    {do ix^D=ixOmin^D, ixOmax^D \}
      eos_rho(ix^D) = w(ix^D, rho_)
      eos_temp(ix^D) = w(ix^D, T_eps_)
      eos_ye(ix^D)  = w(ix^D, ye_)
      if (eos_has_ymu()) then
           eos_ymu(ix^D) = w(ix^D, ymu_)
           call eos_temp_get_all_one_grid(eos_rho(ix^D), eos_temp(ix^D), &
                eos_ye(ix^D), eos_eps(ix^D),&
                 prs = eos_prs(ix^D), &
                 ent = eos_ent(ix^D), ymu=eos_ymu(ix^D))
      else
           call eos_temp_get_all_one_grid(eos_rho(ix^D), eos_temp(ix^D), &
                eos_ye(ix^D), eos_eps(ix^D),&
                 prs = eos_prs(ix^D), &
                 ent = eos_ent(ix^D))
      endif
    {end do^D& \}

   if (nwmax-nw .gt. 40) &
      w(ixI^S,nw+41) = eos_ent(ixI^S)  ! entropy

   if (nwmax-nw .gt. 41) &
      w(ixI^S,nw+42) = eos_prs(ixI^S)  ! pressure

 end associate
 end subroutine specialvar_output
 !=============================================================================
 subroutine specialvarnames_output

   ! newly added variables need to be concatenated with the varnames/primnames string

   include 'amrvacdef.f'
   integer                            :: iw
   !-----------------------------------------------------------------------------
   write(primnames,"(a,a)") trim(primnames),' B2'
   write(wnames,"(a,a)") trim(wnames),' B2'
   write(primnames,"(a,a)") trim(primnames),' divB'
   write(wnames,"(a,a)") trim(wnames),' divB'
   write(primnames,"(a,a)") trim(primnames),' lv'
   write(wnames,"(a,a)") trim(wnames),' lv'
   write(primnames,"(a,a)") trim(primnames),' Omega'
   write(wnames,"(a,a)") trim(wnames),' Omega'
   write(primnames,"(a,a)") trim(primnames),' Q_er1'
   write(wnames,"(a,a)") trim(wnames),' Q_er1'
   write(primnames,"(a,a)") trim(primnames),' kappa_a1'
   write(wnames,"(a,a)") trim(wnames),' kappa_a1'
   write(primnames,"(a,a)") trim(primnames),' kappa_s1'
   write(wnames,"(a,a)") trim(wnames),' kappa_s1'
   write(primnames,"(a,a)") trim(primnames),' kappa_n1'
   write(wnames,"(a,a)") trim(wnames),' kappa_n1'
   write(primnames,"(a,a)") trim(primnames),' Q_nr1'
   write(wnames,"(a,a)") trim(wnames),' Q_nr1'
   write(primnames,"(a,a)") trim(primnames),' Q_er2'
   write(wnames,"(a,a)") trim(wnames),' Q_er2'
   write(primnames,"(a,a)") trim(primnames),' kappa_a2'
   write(wnames,"(a,a)") trim(wnames),' kappa_a2'
   write(primnames,"(a,a)") trim(primnames),' kappa_s2'
   write(wnames,"(a,a)") trim(wnames),' kappa_s2'
   write(primnames,"(a,a)") trim(primnames),' kappa_n2'
   write(wnames,"(a,a)") trim(wnames),' kappa_n2'
   write(primnames,"(a,a)") trim(primnames),' Q_nr2'
   write(wnames,"(a,a)") trim(wnames),' Q_nr2'
   write(primnames,"(a,a)") trim(primnames),' Q_er3'
   write(wnames,"(a,a)") trim(wnames),' Q_er3'
   write(primnames,"(a,a)") trim(primnames),' kappa_a3'
   write(wnames,"(a,a)") trim(wnames),' kappa_a3'
   write(primnames,"(a,a)") trim(primnames),' kappa_s3'
   write(wnames,"(a,a)") trim(wnames),' kappa_s3'
   write(primnames,"(a,a)") trim(primnames),' kappa_n3'
   write(wnames,"(a,a)") trim(wnames),' kappa_n3'
   write(primnames,"(a,a)") trim(primnames),' Q_nr3'
   write(wnames,"(a,a)") trim(wnames),' Q_nr3'
   write(primnames,"(a,a)") trim(primnames),' lfac'
   write(wnames,"(a,a)") trim(wnames),' lfac'
   write(primnames,"(a,a)") trim(primnames),' Lmom1'
   write(wnames,"(a,a)") trim(wnames),' Lmom1'
   write(primnames,"(a,a)") trim(primnames),' Lmom2'
   write(wnames,"(a,a)") trim(wnames),' Lmom2'
   write(primnames,"(a,a)") trim(primnames),' Lmom3'
   write(wnames,"(a,a)") trim(wnames),' Lmom3'
   write(primnames,"(a,a)") trim(primnames),' gxx'
   write(wnames,"(a,a)") trim(wnames),' gxx'
   write(primnames,"(a,a)") trim(primnames),' gyy'
   write(wnames,"(a,a)") trim(wnames),' gyy'
   write(primnames,"(a,a)") trim(primnames),' gzz'
   write(wnames,"(a,a)") trim(wnames),' gzz'
   write(primnames,"(a,a)") trim(primnames),' gxy'
   write(wnames,"(a,a)") trim(wnames),' gxy'
   write(primnames,"(a,a)") trim(primnames),' J1'
   write(wnames,"(a,a)") trim(wnames),' J1'
   write(primnames,"(a,a)") trim(primnames),' Gamma1'
   write(wnames,"(a,a)") trim(wnames),' Gamma1'
   write(primnames,"(a,a)") trim(primnames),' J2'
   write(wnames,"(a,a)") trim(wnames),' J2'
   write(primnames,"(a,a)") trim(primnames),' Gamma2'
   write(wnames,"(a,a)") trim(wnames),' Gamma2'
   write(primnames,"(a,a)") trim(primnames),' J3'
   write(wnames,"(a,a)") trim(wnames),' J3'
   write(primnames,"(a,a)") trim(primnames),' Gamma3'
   write(wnames,"(a,a)") trim(wnames),' Gamma3'
   write(primnames,"(a,a)") trim(primnames),' sqrtg'
   write(wnames,"(a,a)") trim(wnames),' sqrtg'
   write(primnames,"(a,a)") trim(primnames),' nrad1_prim'
   write(wnames,"(a,a)") trim(wnames),' nrad1_prim'
   write(primnames,"(a,a)") trim(primnames),' erad1_prim'
   write(wnames,"(a,a)") trim(wnames),' erad1_prim'
   write(primnames,"(a,a)") trim(primnames),' nrad2_prim'
   write(wnames,"(a,a)") trim(wnames),' nrad2_prim'
   write(primnames,"(a,a)") trim(primnames),' erad2_prim'
   write(wnames,"(a,a)") trim(wnames),' erad2_prim'
   write(primnames,"(a,a)") trim(primnames),' nrad3_prim'
   write(wnames,"(a,a)") trim(wnames),' nrad3_prim'
   write(primnames,"(a,a)") trim(primnames),' erad3_prim'
   write(wnames,"(a,a)") trim(wnames),' erad3_prim'
   write(primnames,"(a,a)") trim(primnames),' ent'
   write(wnames,"(a,a)") trim(wnames),' ent'
   write(primnames,"(a,a)") trim(primnames),' prs'
   write(wnames,"(a,a)") trim(wnames),' prs'

 end subroutine specialvarnames_output


subroutine printlog_special
    ! printlog: calculates volume averaged mean values
    use mod_physaux, only: get_u4
    use mod_eos
    use mod_imhd_intermediate
    use mod_post_process
    use mod_gw_br_Iij_Qij
    use mod_timing
    {#IFDEF D3
        use mod_healpix_det
    }
    use mod_forest,only:nleafs,nleafs_active,nleafs_level
    include 'amrvacdef.f'

    double precision     :: rho_max_out, T_eps_max_out ,ye_min_out, lfac_max_out
    double precision     :: alp_min_out, psi_max_out, temp_local, rho_local
    double precision     :: lfac_max_local, lfac(ixG^T), divB(ixG^T)
    double precision     :: Dpsi6_max_local, Dpsi6_max
    double precision     :: rest_mass_local, rest_mass
    double precision     :: adm_mass_local, adm_mass, adm_vol(ixG^T)
    double precision     :: center_mass_x_local, center_mass_x, center_mass_y_local, center_mass_y
    double precision     :: int_energy_local, int_energy
    double precision     :: em_energy_local, em_energy
    double precision     :: kinetic_energy_local, kinetic_energy
    double precision     :: em_energy_comv_local, em_energy_comv
    double precision     :: em_energy_tor_local, em_energy_tor, em_energy_zz_local
    double precision     :: em_energy_pol_local, em_energy_pol, em_energy_zz
    double precision     :: Q2_ij_local(1:6), Q2_ij(1:6)
    double precision     :: divB_max_local, divB_max_out
    double precision     :: h_00_max_local, h_00_max_out
    double precision     :: delta_max_local, delta_max_out
    double precision     :: sphoutflow(1:healpix_var_num, 1:healpix_det_num), sphhealpix_local(1:healpix_var_num, 1:healpix_det_num)
    logical              :: fileopen
    integer              :: iigrid, igrid, level, iw, i, idir, jdir, kdir, ii, jj
    double precision     :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
    double precision     :: volumeflat(1:nlevelshi)
    integer              :: numlevels, nx^D, nc, ncells, dit
    double precision     :: dtTimeLast, now, cellupdatesPerSecond, activeBlocksPerCore, wctPerCodeTime, timeToFinish
    integer, dimension(1:nlevelshi) :: isum_send, isum_recv
    double precision, dimension(1:nw+1+nlevelshi) :: dsum_send, dsum_recv
    character(len=80)    :: filename, string_i
    character(len=2048)  :: line, line_h, line_hbase
    logical, save        :: opened=.false., already_written
    integer              :: amode, status(MPI_STATUS_SIZE), ix^D
    double precision, dimension(ixG^T, 1:igridstail) :: bsq_comoving, E_em_eul, dvolume, E_em_eul_tor, E_kinetic, E_internal, E_em_eul_pol, E_em_eul_zz
    double precision, dimension(ixG^T, 1:6, 1:igridstail) :: Q2_ij_grid
    double precision     :: Dpsi6(ixG^T, 1:igridstail)
    double precision     :: Wv(ixG^T,1:3)


    ! find extreme values
    rho_max_out = find_global_extreme(rho_, "MAX")
    T_eps_max_out = find_global_extreme(T_eps_, "MAX")
    ye_min_out = 0.0d0
    if (eos_uses_ye()) ye_min_out = find_global_extreme(ye_, "MIN")
    alp_min_out = find_global_extreme(alp_metric_, "MIN")
    psi_max_out = find_global_extreme(psi_metric_, "MAX")
    {#IFDEF GW_BR
      h_00_max_out = find_global_extreme(h_00_, "MAX")
      {#IFDEF DELTA
      delta_max_out = find_global_extreme(delta_br_, "MAX")
      }
    }

    {#IFDEF STAGGERED
    ! fixme: not a good way to print divB max, it contracts to the vtu output, why?
    ! find max divB
    divB_max_local = 0.0d0
    do iigrid = 1, igridstail
       igrid = igrids(iigrid)
       call div_staggered(ixG^LL,ixM^LL,ps(igrid),divB(ixG^T))
       divB_max_local = max(divB_max_local, maxval(divB(ixM^T)) )
    end do
    call MPI_ALLREDUCE(divB_max_local, divB_max_out, 1, mpi_double_precision, MPI_MAX, icomm, ierrmpi)
    }

    ! find max of D*phi6
    ! and calculate the rest mass, mass center at the same time
    Dpsi6_max_local = 0.0d0
    center_mass_x_local = 0.0d0
    center_mass_y_local = 0.0d0
    rest_mass_local = 0.0d0

    do iigrid = 1, igridstail
        igrid = igrids(iigrid)
        Dpsi6(ixM^T, iigrid) = pw(igrid)%w(ixM^T, D_)
        {#IFDEF DY_SP
            Dpsi6(ixM^T, iigrid) = Dpsi6(ixM^T, iigrid)*(pw(igrid)%w(ixM^T, psi_metric_))**6
        }
        Dpsi6_max_local = max(Dpsi6_max_local, maxval(Dpsi6(ixM^T, iigrid)))
        rest_mass_local = rest_mass_local+sum(Dpsi6(ixM^T, iigrid)*pgeo(igrid)%dvolume(ixM^T))
        center_mass_x_local = center_mass_x_local+sum(px(igrid)%x(ixM^T, 1)*Dpsi6(ixM^T, iigrid)*pgeo(igrid)%dvolume(ixM^T))
        center_mass_y_local = center_mass_y_local+sum(px(igrid)%x(ixM^T, 2)*Dpsi6(ixM^T, iigrid)*pgeo(igrid)%dvolume(ixM^T))
    enddo
    call MPI_ALLREDUCE(Dpsi6_max_local, Dpsi6_max, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)
    call MPI_ALLREDUCE(rest_mass_local, rest_mass, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(center_mass_x_local, center_mass_x, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(center_mass_y_local, center_mass_y, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    center_mass_x = center_mass_x/rest_mass
    center_mass_y = center_mass_y/rest_mass

    ! calculate volume integral and average
    rest_mass_local = 0.0d0
    adm_mass_local = 0.0d0
    int_energy_local = 0.0d0
    em_energy_local = 0.0d0
    em_energy_comv_local = 0.0d0
    em_energy_tor_local = 0.0d0
    em_energy_pol_local = 0.0d0
    kinetic_energy_local = 0.0d0
    em_energy_zz_local = 0.0d0
    Q2_ij_local = 0.0d0
    do iigrid = 1, igridstail
        igrid = igrids(iigrid)
        call imhd_get_intermediate_variables(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T, 1:nw), px(igrid)%x(ixG^T, 1:ndim), &
                                             lfac=lfac(ixG^T), b2=bsq_comoving(ixG^T, iigrid))
        lfac_max_local = max(lfac_max_local, maxval(lfac(ixM^T)) )
        dvolume(ixM^T, iigrid) = pgeo(igrid)%dvolume(ixM^T)*pw(igrid)%w(ixM^T,psi_metric_)**6
        call cal_adm_mass_volume_contribution(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, adm_vol)
        adm_mass_local = adm_mass_local+sum(adm_vol(ixM^T)*dvolume(ixM^T, iigrid))
        if (eos_uses_ye()) then
            {do ix^D = ixMlo^D, ixMhi^D \}
                temp_local = pw(igrid)%w(ix^D, T_eps_)
                rho_local = pw(igrid)%w(ix^D, rho_)
                if (eos_has_ymu()) then
                    call eos_temp_get_all_one_grid(rho_local, temp_local, pw(igrid)%w(ix^D, ye_), E_internal(ix^D, iigrid), &
                                                   ymu=pw(igrid)%w(ix^D, ymu_))
                else
                    call eos_temp_get_all_one_grid(rho_local, temp_local, pw(igrid)%w(ix^D, ye_), E_internal(ix^D, iigrid))
                endif
            {enddo ^D&\}
        else
            E_internal(ixM^T, iigrid) = pw(igrid)%w(ixM^T, T_eps_)
        endif
        int_energy_local = int_energy_local+sum(E_internal(ixM^T, iigrid)*&
                pw(igrid)%w(ixM^T,rho_)*lfac(ixM^T)*dvolume(ixM^T, iigrid))
        ! calculate manetic energy
        call imhd_get_total_em_energy(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T, 1:nw), px(igrid)%x(ixG^T, 1:ndim), E_em_eul(ixG^T, iigrid))
        call imhd_get_em_energy_cylind_decomp(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T, 1:nw), px(igrid)%x(ixG^T, 1:ndim), &
            E_em_eul_pol(ixG^T, iigrid), E_em_eul_tor(ixG^T, iigrid), E_em_eul_zz(ixG^T, iigrid), center_mass_x, center_mass_y)
        call total_kinetic_energy(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T, 1:nw), px(igrid)%x(ixG^T, 1:ndim), E_kinetic(ixG^T, iigrid))
        em_energy_comv_local = em_energy_comv_local+sum(bsq_comoving(ixM^T, iigrid)/2*dvolume(ixM^T, iigrid))
        em_energy_local = em_energy_local+sum(E_em_eul(ixM^T, iigrid)*dvolume(ixM^T, iigrid))
        em_energy_tor_local = em_energy_tor_local+sum(E_em_eul_tor(ixM^T, iigrid)*dvolume(ixM^T, iigrid))
        em_energy_pol_local = em_energy_pol_local+sum(E_em_eul_pol(ixM^T, iigrid)*dvolume(ixM^T, iigrid))
        em_energy_zz_local = em_energy_zz_local+sum(E_em_eul_zz(ixM^T, iigrid)*dvolume(ixM^T, iigrid))
        kinetic_energy_local = kinetic_energy_local+sum(E_kinetic(ixM^T, iigrid)*dvolume(ixM^T, iigrid))
        {#IFDEF GW_BR
        ! GW analysis: quadruople
        call compute_Q2_ij(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, &
                           Q2_ij_grid(ixG^T,1:6,iigrid), lfac(ixG^T), dvolume(ixG^T,iigrid))
        do idir = 1, 6
          Q2_ij_local(idir) = Q2_ij_local(idir) + sum(Q2_ij_grid(ixM^T,idir,iigrid))
        enddo
        }
    enddo

    ! calculate healpix flux output
    {#IFDEF D3
    if (healpix_det_open) then
        sphoutflow = 0.0d0
        call surface_integral_outflow(it, healpix_det_num, healpix_flux_num, sphhealpix_local(1:healpix_flux_num, 1:healpix_det_num))
        call MPI_ALLREDUCE(sphhealpix_local(1:healpix_flux_num, 1:healpix_det_num), sphoutflow(1:healpix_flux_num, 1:healpix_det_num), &
            healpix_flux_num*healpix_det_num, mpi_double_precision, MPI_SUM, icomm, ierrmpi)
    endif
    }

    ! gather all the results!
    call MPI_ALLREDUCE(lfac_max_local, lfac_max_out, 1, mpi_double_precision, &
          MPI_MAX, icomm, ierrmpi)
    call MPI_ALLREDUCE(adm_mass_local, adm_mass, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(int_energy_local, int_energy, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(em_energy_local, em_energy, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(em_energy_comv_local, em_energy_comv, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(em_energy_tor_local, em_energy_tor, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(em_energy_pol_local, em_energy_pol, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(em_energy_zz_local, em_energy_zz, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(kinetic_energy_local, kinetic_energy, 1, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)
    call MPI_ALLREDUCE(Q2_ij_local, Q2_ij, 6, mpi_double_precision, &
          MPI_SUM, icomm, ierrmpi)

    if (mype==0) then
        ! To compute cell updates per second, we do the following:
        nx^D=ixMhi^D-ixMlo^D+1;
        nc={nx^D*}
        ncells = nc * nleafs_active
        ! assumes the number of active leafs have not changed since last compute.
        now        = MPI_WTIME()
        dit        = it - itTimeLast
        dtTimeLast = now - timeLast
        itTimeLast = it
        timeLast   = now
        cellupdatesPerSecond = dble(ncells) * dble(nstep) * dble(dit) / (dtTimeLast * dble(npe))

        ! write file header
        if (.not.opened) then
            ! generate filename
            write(filename,"(a,a)") TRIM(filenamelog),".log"

            amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
            amode=ior(amode,MPI_MODE_APPEND)
            call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, &
                               MPI_INFO_NULL,log_fh,ierrmpi)
            opened=.true.

           ! if (snapshotini==-1) then
                call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout), MPI_CHARACTER,status,ierrmpi)
                call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)
                line_hbase = "rho_max, T_eps_max, ye_min, lfac_max, alp_min, psi_max, h_00_max, delta_max, Q2_11, Q2_12, Q2_13, Q2_22, Q2_23, Q2_33, divB_max, &
                    Dpsi6_max, rest_mass, adm_mass, int_energy, em_energy_tot, em_energy_pol, em_energy_tor, em_energy_zz, em_energy_comv, kinetic_energy, &
                    cmx, cmy"
                {#IFDEF D3
                    if (healpix_det_open) line_hbase = trim(adjustl(line_hbase))//trim(adjustl(flux_string_for_logfile_header))
                }
                do i=1, mxnest
                    write(string_i, *) i
                    line_hbase = trim(adjustl(line_hbase))//", nleafs_level" //trim(adjustl(string_i))
                end do
                write(line_h,*) "it, t, dt, "//trim(line_hbase)//", cellupdatesPerSecond"
                call MPI_FILE_WRITE(log_fh, line_h, len_trim(line_h), MPI_CHARACTER, status,ierrmpi)
            !endif
        end if
        call MPI_FILE_WRITE(log_fh, achar(10), 1, MPI_CHARACTER, status, ierrmpi)

        ! write one line
        already_written = .false.
        {#IFDEF D3
            if (healpix_det_open) then
                write(line,*)it, t, dt, rho_max_out, T_eps_max_out, ye_min_out, lfac_max_out, alp_min_out, psi_max_out, &
                    h_00_max_out, delta_max_out, Q2_ij(1:6), divB_max_out, Dpsi6_max, &
                    rest_mass, adm_mass, int_energy, em_energy, em_energy_pol, em_energy_tor, em_energy_zz, em_energy_comv, kinetic_energy, &
                    center_mass_x, center_mass_y, transpose(sphoutflow(1:healpix_flux_num, 1:healpix_det_num)), nleafs_level(1:mxnest), cellupdatesPerSecond
                already_written = .true.
            endif
        }
        if (.not. already_written) then
            write(line,*)it, t, dt, rho_max_out, T_eps_max_out, ye_min_out, lfac_max_out, alp_min_out, psi_max_out, &
                h_00_max_out, delta_max_out, Q2_ij(1:6), divB_max_out, Dpsi6_max, &
                rest_mass, adm_mass, int_energy, em_energy, em_energy_pol, em_energy_tor, em_energy_zz, em_energy_comv, kinetic_energy, &
                center_mass_x, center_mass_y, nleafs_level(1:mxnest), cellupdatesPerSecond
        endif
        call MPI_FILE_WRITE(log_fh, line, len_trim(line), MPI_CHARACTER, status, ierrmpi)
    end if


end subroutine printlog_special



!========================================================================================
!==========                  user specified behavior                           ==========
!========================================================================================



!> before the main loop, we improve the initial data here --> initialize metric and calculate the lfac --> conserve
subroutine usr_improve_initial_condition
    use mod_metric, only: lower3
    use mod_eos
    use mod_xns
    use mod_c2b
    !use mod_cfc
    include 'amrvacdef.f'

    integer              :: ix^D, igrid, iigrid
    double precision     :: vD(ixG^T, 1:ndir)
    double precision     :: sqrV(ixG^T)
    
end subroutine usr_improve_initial_condition

  !=============================================================================
subroutine perturb(ixI^L,ixO^L, w, x)

    !USE IFPORT
    include 'amrvacdef.f'

    integer, intent(in)              :: ixI^L, ixO^L
    double precision, intent(in)     :: x(ixI^S,1:ndim)
    double precision, intent(inout)  :: w(ixI^S,1:nw)
    ! .. local ..
    integer                          :: ix^D
    double precision                 :: r(ixI^S)
    double precision                 :: Xr, rdble
    double precision, parameter      :: amp=0.04d0
    !-----------------------------------------------------------------------------

end subroutine perturb
!=============================================================================
subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,s)

    use mod_metric,only: get_g_component, get_alpha, get_beta, lower3

    ! special boundary types, user defined
    include 'amrvacdef.f'

    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt
    type(state), intent(inout)   :: s
    !-----------------------------------------------------------------------------

    ! .. local variables ..
    integer                           :: ix^D, iiw, ixD^L, ixDB^L
    double precision,dimension(ixI^S) :: a,b,c,d, ut
    !-----------------------------------------------------------------------------
    call mpistop("Special boundary conditions are not defined.")
end subroutine specialbound_usr

!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

    ! this subroutine is ONLY to be used for computing auxiliary variables
    ! which happen to be non-local (like div v), and are in no way used for
    ! flux computations. As auxiliaries, they are also not advanced

    include 'amrvacdef.f'

    integer, intent(in):: igrid,level,ixI^L,ixO^L
    double precision, intent(in):: qt,x(ixI^S,1:ndim)
    double precision, intent(inout):: w(ixI^S,1:nw)
    !-----------------------------------------------------------------------------


end subroutine process_grid_usr
!=============================================================================
subroutine correctaux_usr(ixI^L,ixO^L,w,x,patchierror,subname)

    use mod_metric
    include 'amrvacdef.f'

    integer, intent(in)            :: ixI^L, ixO^L
    integer, intent(inout)         :: patchierror(ixG^T)
    character(len=*), intent(in)   :: subname
    double precision, intent(inout):: w(ixI^S,1:nw)
    double precision, intent(in)   :: x(ixI^S,1:ndim)
    ! .. local ..
    double precision               :: wtmp(ixI^S,1:nw)
    double precision               :: xBL(ixI^S,1:ndim)
    logical                        :: patchw(ixG^T)
    !-----------------------------------------------------------------------------
    call mpistop("no need correctaux")

end subroutine correctaux_usr
!=============================================================================
subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)

    ! internal boundary, user defined
    !
    ! This subroutine can be used to artificially overwrite ALL conservative
    ! variables in a user-selected region of the mesh, and thereby act as
    ! an internal boundary region. It is called just before external (ghost cell)
    ! boundary regions will be set by the BC selection. Here, you could e.g.
    ! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
    ! which can be used to identify the internal boundary region location.
    ! Its effect should always be local as it acts on the mesh.
    !

    include 'amrvacdef.f'

    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)

    ! .. local ..
    !----------------------------------------------------------------------------
    
end subroutine bc_int
  !=============================================================================
subroutine fixp_usr(ixI^L,ixO^L,w,x)
    include 'amrvacdef.f'

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(inout)    :: w(ixI^S,1:nw)
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    ! .. local ..
    !----------------------------------------------------------------------------
end subroutine fixp_usr

subroutine fixp_usr_pt(ix^D,w_pt,x_pt)
    include 'amrvacdef.f'
    integer, intent(in)                :: ix^D
    double precision, intent(inout)    :: w_pt(1:nw)
    double precision, intent(in)       :: x_pt(1:ndim)
                   
end subroutine fixp_usr_pt

!=============================================================================
subroutine normalize_initial()

    include 'amrvacdef.f'

    ! .. local ..
!--------------------------------


end subroutine normalize_initial
!=============================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

    include 'amrvacdef.f'

    integer, intent(in)             :: ixG^L, ixO^L
    integer, intent(inout)          :: flag
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision, intent(in)    :: x(ixG^S,1:ndim)

    ! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
    ! flag=0  : Treat as normal domain
    ! flag=1  : Treat as passive, but reduce by safety belt
    ! flag=2  : Always treat as passive

    !-----------------------------------------------------------------------------

end subroutine flag_grid_usr
  !=============================================================================
subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = 1001 !getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
end subroutine init_random_seed
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

! do iw= iw^LIM
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixI^L,ixO^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ixO^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================


subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen
use mod_eos
use NSimport
use mod_c2b_cart3d
use Init_c2b_cart_3d_interface
include 'amrvacdef.f'

integer, intent(in) :: igrid, level, ixI^L, ixO^L
double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
integer, intent(inout) :: refine, coarsen
double precision             :: ratio(ixI^S)
integer                      :: i_level, to_level, found_level
double precision             :: phi_grav ! relativistic gravitational potential :=  1 - alp
double precision             :: phi_grav_cut
double precision, parameter  :: phi_grav_max = 0.2d0
double precision             :: rho_cut = 1.d-5
double precision, parameter  :: r_in = 12.0d5 * length_gf ! 100 km
double precision             :: radius(ixI^S), radius2(ixI^S)
double precision             :: r_cut, r_cut2, r_max, r_min, rho_max, rela_diff(3)
double precision             :: x_middle(1:3), shifted_x(ixI^S), shifted_y(ixI^S)
double precision             :: aux_denstiy, density_lev, rho_threshold, auxv1, auxv2


    refine  = 0
    coarsen = 0
    shifted_x(ixI^S) = x(ixI^S,1)-shift_origin
    shifted_y(ixI^S) = x(ixI^S,2)-shift_origin

    select case(refine_type)
    case(1)
    !noting
        if (cut_r_core) then
          if ( (minval(shifted_x(ixO^S)) .lt. r_core) ) then
             refine  = -1
             coarsen = 1
          endif
        endif
    case(2)
        phi_grav = minval( w(ixO^S, alp_metric_) )
        phi_grav = 1.0d0 - phi_grav

        to_level = mxnest
        phi_grav_cut = phi_grav_max
        do i_level = mxnest-1, 1, -1
           phi_grav_cut = phi_grav_cut / 2.0d0
           if ( phi_grav < phi_grav_cut ) then
              to_level = i_level
           end if
        end do

        if ( level > to_level ) then
           refine = -1
           coarsen = 1
        else if ( level < to_level ) then
           refine = 1
           coarsen = -1
        end if

    case (3)
         !r_cut = 20.0d0
        radius(ixO^S) = dsqrt( shifted_x(ixO^S)**2 + shifted_y(ixO^S)**2 + x(ixO^S,3)**2  )
        if ( minval(radius(ixO^S)) > r_outer) then
            refine = -1
            coarsen = 1
            return
        endif

        r_cut = r_entire
        if ( minval(radius(ixO^S)) < r_cut) then
            refine = 1
            coarsen = -1
            return
        else
            refine = -1
            coarsen = 1
        endif

    case (4) ! new alp refine
        ! check if the density is almost/already
        ! above rho_cut
        rho_max = maxval( w(ixO^S, rho_) )
        if (rho_max >= rho_cut) then
           ! always at the finest level
           refine = 1
           coarsen = -1
           return
        end if

        ! check if this is the most outer block
        if ( (maxval(shifted_x(ixO^S))+dx(1,level)) >= xprobmax1 .or. &
             (minval(shifted_x(ixO^S))-dx(1,level)) <= xprobmin1 .or. &
             (maxval(shifted_y(ixO^S))+dx(2,level)) >= xprobmax2 .or. &
             (minval(shifted_y(ixO^S))-dx(2,level)) <= xprobmin2 .or. &
             (maxval(x(ixO^S,3))+dx(3,level)) >= xprobmax3 .or. &
            ((minval(x(ixO^S,3))-dx(3,level)) <= xprobmin3 .and. abs(xprobmin3) <= smalldouble) &
             ) then
           ! always at the lowest level
           refine = -1
           coarsen = 1
           return
        end if

        phi_grav = minval( w(ixO^S, alp_metric_) )
        phi_grav = 1.0d0 - phi_grav

        to_level = mxnest
        phi_grav_cut = phi_grav_max
        do i_level = mxnest-1, 1, -1
           phi_grav_cut = phi_grav_cut / 2.0d0
           if ( phi_grav < phi_grav_cut ) then
              to_level = i_level
           end if
        end do

        if ( level > to_level ) then
           refine = -1
           coarsen = 1

        ! if the current level is higher than the maximum allowed level,
        ! we strictly limit the level.
        ! However, we do not restrict if the level is too low.
        ! The goal here is just put the upper bound of the level.

        ! if refine_criterion == 0, purely depents on user define level,
        ! then do the following
        else if ( level < to_level ) then
           refine = 1
           coarsen = -1
        end if

    case (5)
        x_middle(1) = (shifted_x(ixOmax^D) + shifted_x(ixOmin^D))/2.d0
        x_middle(2) = (shifted_y(ixOmax^D) + shifted_y(ixOmin^D))/2.d0
        x_middle(3) = (x(ixOmax^D,3) + x(ixOmin^D,3))/2.d0
        if (mirror_zplane) then
            found_level = find_point_in_which_level(my_c2b_3d_profile%rho, x_middle(1), x_middle(2), abs(x_middle(3))) ! consider z-symm in original carpet data
        else
            found_level = find_point_in_which_level(my_c2b_3d_profile%rho, x_middle(1), x_middle(2), x_middle(3)) ! consider z-symm in original carpet data
        endif
        if (found_level/=-1) then
            rela_diff(1) = (my_c2b_3d_profile%rho%all_refine_levels(found_level)%delta_x(1)-dx(1, level))/dx(1, level)
            rela_diff(2) = (my_c2b_3d_profile%rho%all_refine_levels(found_level)%delta_x(2)-dx(2, level))/dx(2, level)
            rela_diff(3) = (my_c2b_3d_profile%rho%all_refine_levels(found_level)%delta_x(3)-dx(3, level))/dx(3, level)
            write(*, *) "found in level: ", found_level, ", bhac level: ", level, "  ", dx(1, level), "  ", dx(2, level), "  ", dx(3, level), "  ", my_c2b_3d_profile%rho%all_refine_levels(found_level)%delta_x(1)
            if (rela_diff(1)<-0.1 .and. rela_diff(2)<-0.1 .and. rela_diff(3)<-0.1) then
                write(*, *) "   force refinement"
                refine = 1
                coarsen = -1
            else if  (rela_diff(1)>0.1 .and. rela_diff(2)>0.1 .and. rela_diff(3)>0.1) then
                write(*, *) "   force coarsen"
                refine = -1
                coarsen = 1
            else
                write(*, *) "   not need to be refine or coarsen"
                refine = -1
                coarsen = -1
            endif
        endif
    case (6) ! AMR using density profile for BNS merger
        aux_denstiy = maxval(w(ixO^S, rho_))
        density_lev = log10(aux_denstiy/rho_gf)
        write(*, *) "debug amr: ", aux_denstiy, density_lev
        if (density_lev<10) then
            write(*, *) "  AMR force coarsen"
            refine = -1
            coarsen = 1
            return
        endif
        ! aim: (>14, mxnest), (>13, mxnest-1), ....
        to_level = mxnest
        rho_threshold = 14
        findlevel: do i_level = mxnest, 1, -1
            if ( density_lev > rho_threshold ) then
                to_level = i_level
                exit findlevel
            end if
            rho_threshold = rho_threshold-0.80
        enddo findlevel
        if ( level > to_level ) then
            refine = -1
            coarsen = 1
        else if ( level < to_level ) then
            refine = 1
            coarsen = -1
        end if
    case (7) ! AMR using density profile for BNS merger
        aux_denstiy = maxval(w(ixO^S, rho_))
        density_lev = log10(aux_denstiy/rho_gf)
        radius(ixI^S) = dsqrt( shifted_x(ixI^S)**2 + shifted_y(ixI^S)**2)
        auxv1 = minval(radius)
        auxv2 = minval(x(ixO^S, 3)**2)
        if (auxv1>80 .or. auxv2>50.) then
            refine = -1
            coarsen = 1
            return
        endif
        if (auxv1<25. .and. auxv2<15.) then
            if (level<mxnest) then
                write(*, *) "trigger 1"
                refine = 1
                coarsen = -1
                return
            endif
        else if (auxv1<40. .and. auxv2<25.) then 
            if (level<mxnest-1) then 
                write(*, *) "trigger 2"
                refine = 1
                coarsen = -1 
            else if (level>mxnest-1) then
                refine = -1
                coarsen = 1
            else
                ! do nothing
            endif
        !else if (auxv1<80. .and. auxv2<40.) then 
        !    if (level<mxnest-2) then 
        !        write(*, *) "trigger 3"
        !        refine = 1
        !        coarsen = -1 
        !    else if (level>mxnest-2) then
        !        refine = -1
        !        coarsen = 1
        !    else
        !        ! do nothing
        !    endif
        endif
     case default
        call mpistop('pls specify a special refine case, otherwise, set it to one')
     end select

end subroutine specialrefine_grid




!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

include 'amrvacdef.f'

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0

subroutine usr_atmo_pt(ix^D, w_pt, x_pt)
   ! The atmosphere treatment
   include 'amrvacdef.f'
   integer, intent(In)              :: ix^D
   double precision, intent(in)     :: x_pt(1:ndim)
   double precision, intent(inout)  :: w_pt(1:nw)

end subroutine usr_atmo_pt
