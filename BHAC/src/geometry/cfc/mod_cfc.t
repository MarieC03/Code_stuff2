!> Module for cfc
!< Note that this metric solver support with mod_grhd only!
module mod_cfc
  use mod_cfc_parameters
  use mod_cfc_alp
  use mod_cfc_beta
  use mod_cfc_psi
  use mod_gw_br
  use mod_eos
  implicit none
  public

  ! public methods
  public :: cfc_metric_init
  public :: cfc_solve_metric
  public :: cfc_check_metric
  public :: cfc_solver_activate
  public :: cfc_update_flag
!  public :: cfc_metric_interpolation

  contains

  logical function cfc_update_flag()
    include 'amrvacdef.f'
    cfc_update_flag = .False.

    if (use_cfc .and. cfc_evolve) then

       ! check now in which stage to adjust cfc_dit_update/cfc_dt_update
       if (t < cfc_t_end_of_stage1) then
       ! stage 1
          !keep default one
       else if (t < cfc_t_end_of_stage2 .and. t .ge. cfc_t_end_of_stage1) then
       ! stage 2
          !if (mype == 0) print*, 'CFC arrive stage2'
          cfc_dit_update = cfc_dit_update_stage2
          cfc_dt_update  = cfc_dt_update_stage2
       else if (t .ge. cfc_t_end_of_stage2) then     
       ! stage 3
          !if (mype == 0) print*, 'CFC arrive stage3'
          cfc_dit_update = cfc_dit_update_stage3
          cfc_dt_update  = cfc_dt_update_stage3
       else
          call mpistop('You did not set cfc_time_of_stage1,2,3 in correct way')
       endif

       if ( it >= cfc_it_last_update - 1 + cfc_dit_update ) then
          ! note that metric will not be updated if the time difference is smaller than cfc_smallest_dt.
          if ( t >= cfc_t_last_update + cfc_smallest_dt .or. t == 0.0d0) then
             cfc_update_flag=.True.
          end if
       end if

       if ( t >= cfc_t_last_update + cfc_dt_update ) then
          cfc_update_flag=.True.
       end if

       if (.not. cfc_update_flag .and. use_check_psi) call cfc_check_metric(cfc_update_flag)

       if (cfc_update_flag) then
          cfc_t_last_update = t
          cfc_it_last_update = it
       end if

    end if
  end function cfc_update_flag

  subroutine cfc_solver_activate()
    use mod_multigrid_coupling
    include 'amrvacdef.f'
    integer                      :: n

    if (mype==0) then
       print*,'-----------------------------------------------------------------------------'
       print*,'-------------------------------CFC activated---------------------------------'
       print*,'-----------------------------------------------------------------------------'
    end if

    ! set the parameters here
    cfc_print       = 1000
    !cfc_it_max      = 1000000
    ! have put to parameters
    !cfc_n_cycle(1)  = 2!15
    !cfc_n_cycle(2)  = 2!15
!    cfc_psi_tol_init = 1.0d-3
!    cfc_dt_update = 2.030017d0 ! time_gf = 2.03d5
!    reconstruct_cfc = .false.
!    cfc_redblack = .true.
!    cfc_n_interpolation = 4
!    cfc_beta_robin_bc = .true.
!    cfc_test_alp = .false.

    if (mype == 0) then
       write(*,*) 'cfc_it_max', cfc_it_max
       write(*,*) 'cfc_evolve', cfc_evolve
       write(*,*) 'cfc_tol :', cfc_tol(:)
       write(*,*) 'cfc_tol_evolve :', cfc_tol_evolve(:)
       write(*,*) 'cfc_smallest_dt', cfc_smallest_dt 
       write(*,*) 'cfc_t_end_of_stage1', cfc_t_end_of_stage1
       write(*,*) 'cfc_t_end_of_stage2', cfc_t_end_of_stage2

       write(*,*) 'cfc_dit_update', cfc_dit_update
       write(*,*) 'cfc_dit_update_stage2', cfc_dit_update_stage2 
       write(*,*) 'cfc_dit_update_stage3', cfc_dit_update_stage3

       write(*,*) 'cfc_dt_update', cfc_dt_update
       write(*,*) 'cfc_dt_update_stage2', cfc_dt_update_stage2 
       write(*,*) 'cfc_dt_update_stage3', cfc_dt_update_stage3
  
       write(*,*) 'cfc_n_cycle(1:2)', cfc_n_cycle(:)
       write(*,*) 'cfc_beta_robin_bc', cfc_beta_robin_bc
    endif

    if ( cfc_dt_update < 0.0d0 ) then
       call mpistop(" cfc_dt_update can not be negative")
    end if
    if ( cfc_n_interpolation < 2 ) then
       call mpistop(" cfc_n_interpolation < 2")
    end if

    if ( mod(cfc_n_interpolation,2) == 0 ) then
       n = cfc_n_interpolation / 2
    else
       n = (cfc_n_interpolation-1) / 2
    endif
    ! fixme: this seems no needed, as the code not yet read dixB
!    dixB = max(dixB, n)
!    if ( n > dixB) then
!       write(*,*) "cfc_n_interpolation = ", cfc_n_interpolation
!       write(*,*) "number of ghostcells needed = ", n
!       write(*,*) "dixB = ", dixB
!       call mpistop(" cfc_n_interpolation is not match with dixB")
!    end if

    cfc_t_last_update = t
    cfc_it_last_update = it
    
    ! activate the solver
    use_multigrid = .True.
    use_cfc = .True.
  end subroutine cfc_solver_activate

  subroutine cfc_metric_init()
    use mod_multigrid_coupling
    use mod_forest
    use mod_imhd_intermediate
    use mod_small_values
    include 'amrvacdef.f'

    integer                        :: id, iigrid, igrid
    double precision               :: psi_err = huge(0.0d0)
    double precision               :: psi_err_old = huge(0.0d0)
    double precision               :: local_err = 0.0d0
    integer, parameter             :: mg_it_max = 10000
    integer                        :: nc, lvl, mg_it
    type(tree_node), pointer       :: pnode

    double precision, allocatable  :: Aij(:^D&,:,:,:)
    double precision, allocatable  :: A2(:^D&,:)
    double precision, allocatable  :: lfac(:^D&,:)
    double precision, allocatable  :: lfac2(:^D&,:)
    integer                        :: idims, ixC^L, ix^D

    double precision               :: psi_xns_err, psi_xns_err_local
    double precision               :: alp_xns_err, alp_xns_err_local
    double precision               :: array_max_local(1:3), array_max(1:3), lfac_max_local, lfacmax

    integer :: i

    ! limitation:
    ! before here, must fill lfac and xi to get correct conserve variables

    allocate(Aij(ixG^T,1:3,1:3,1:igridstail))
    allocate(A2(ixG^T,1:igridstail))
    allocate(lfac(ixG^T,1:igridstail))
    allocate(lfac2(ixG^T,1:igridstail))
    Aij  = 0.0d0
    A2   = 0.0d0
    lfac = 1.0d0
    lfac2 = 1.0d0
 
    ! initialize the solver
    if (cfc_redblack) then
       mg%smoother_type = mg_smoother_gsrb
    else
       mg%smoother_type = mg_smoother_gs
    end if

    mg%n_cycle_up   = cfc_n_cycle(1)
    mg%n_cycle_down = cfc_n_cycle(2)

    ! Before initialization, make sure you have psi6 cons and correct v^i, lfac already

    call getbc(t,ps,psCoarse)

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
 
       !if (.not. use_lfac_profile_to_init) then
       !! store vi instead of Wvi inside pwold       
       !{^C& pwold(igrid)%w(ixM^T, u^C_) = pw(igrid)%w(ixM^T, u^C_) &
       !      / pwold(igrid)%w(ixM^T,lfac_)   \}
       !endif
       ! make sure all psi and alp are physical
       pw(igrid)%w(ixG^T, psi_metric_) = max( pw(igrid)%w(ixG^T, psi_metric_), 1.0d0 )
       pw(igrid)%w(ixG^T, alp_metric_) = min( max( pw(igrid)%w(ixG^T, alp_metric_), 0.0d0 ), 1.0d0)

       !TEST
      ! pw(igrid)%w(ixG^T, psi_metric_) = 1.0d0
      ! pw(igrid)%w(ixG^T, alp_metric_) = 1.0d0
      ! pw(igrid)%w(ixG^T, beta_metric1_) = 0.0d0
      ! pw(igrid)%w(ixG^T, beta_metric2_) = 0.0d0
      ! pw(igrid)%w(ixG^T, beta_metric3_) = 0.0d0

       {#IFDEF M1
       {^KSP&         
       {#IFNDEF VAR1
       pwold(igrid)%w(ixG^T,nrad^KSP_) = pw(igrid)%w(ixG^T,nrad^KSP_)
       pwold(igrid)%w(ixG^T,erad^KSP_) = pw(igrid)%w(ixG^T,erad^KSP_)
       {^C& pwold(igrid)%w(ixG^T,frad^KSP^C_) = pw(igrid)%w(ixG^T,frad^KSP^C_) \}
       }
       !pw(igrid)%w(ixG^T,frad^KSP1_) = pwold(igrid)%w(ixG^T,frad^KSP1_) 
       !pw(igrid)%w(ixG^T,frad^KSP2_) = pwold(igrid)%w(ixG^T,frad^KSP2_) 
       !pw(igrid)%w(ixG^T,frad^KSP3_) = pwold(igrid)%w(ixG^T,frad^KSP3_) 
       \}
       }

       if (initialize_with_cons) then
              ! you cant, since you dont have correct gammaij yet
              ! call primitive(ixG^LL,ixG^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
              ! fixed initial cons
         !     pwold(igrid)%w(ixG^T,D_)          = pw(igrid)%w(ixG^T,D_)
         !     pwold(igrid)%w(ixG^T,tau_)        = pw(igrid)%w(ixG^T,tau_)
         !{^C& pwold(igrid)%w(ixG^T,S^C_)        = pw(igrid)%w(ixG^T,S^C_) \}
         !     pwold(igrid)%w(ixG^T,psi_metric_) = pw(igrid)%w(ixG^T,psi_metric_)
              ! every conserved quantity needed to times psi_metric_ **6
              call conformal_transformation(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
              pwold(igrid)%w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
              call mpistop('stop using initialize with cons')
       else
              ! warning: fix cons, does not work, you changed the actually prim profiles after update psi
              ! --> change T, ye, eps
              ! Therefore, we fix prim
       
             ! call conserve(ixG^LL,ixG^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim),patchfalse)
             ! ! every conserved quantity needed to times psi_metric_ **6
             ! call conformal_transformation(ixG^LL, ixG^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
          if (.not. fix_prim_init) then
              pw(igrid)%w(ixG^T,D_)   = 0.0d0
              if (eos_type == tabulated) pw(igrid)%w(ixG^T,Dye_) = 0.0d0
               pw(igrid)%w(ixG^T,tau_) = 0.0d0
          {^C& pw(igrid)%w(ixG^T,S^C_) = 0.0d0 \}
              ! me B-field test:
           !{C pw(igrid)%w(ixG^T,b^C_) = 0.0d0 }

             !dont set pw-M1 to zero else amrvacusr is overwritten and conserve does not calculate them!
             !this is only fine in case we do pwold = pw before and pw = pwold here
              {#IFDEF M1 
               {^KSP&
               {#IFNDEF VAR1
               pw(igrid)%w(ixG^T,nrad^KSP_) = 0.0d0
               pw(igrid)%w(ixG^T,erad^KSP_) = 0.0d0
               {^C& pw(igrid)%w(ixG^T,frad^KSP^C_) = 0.0d0 \}
               }
               \}               
              }
              !M1_TODO: check here for NAN in any var!
              {#IFDEF NAN_CHECK
              write(75,*) "igrid",igrid
              do i=1,nw
                write(75,*) "i",i,pw(igrid)%w(ixG^T,i)
              end do
              }
              call conserve(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim),patchfalse)
              ! every conserved quantity needed to times psi_metric_ **6
              call conformal_transformation(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
              
              !M1_TODO: check here for NAN in any var!
              {#IFDEF NAN_CHECK
              write(76,*) "igrid",igrid
              do i=1,nw
                write(76,*) "i",i,pw(igrid)%w(ixG^T,i)
              end do
              }
              pwold(igrid)%w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
          else
              pwold(igrid)%w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
              ! assume you have a guess of Wv first
              call imhd_get_intermediate_variables(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, lfac=lfac(ixG^T,iigrid))
              ! store the current v in pw
          {^C& pwold(igrid)%w(ixM^T, u^C_) = pw(igrid)%w(ixM^T, u^C_) &
                                             / lfac(ixM^T,iigrid)    \}
          endif


       endif

    end do
    !$OMP END PARALLEL DO

    !M1_TODO: check here for NAN in any var!
    {#IFDEF NAN_CHECK   
    do iigrid=1,igridstail; igrid=igrids(iigrid);
    write(77,*) "igrid",igrid
    do i=1,nw
      write(77,*) "i",i,pw(igrid)%w(ixG^T,i)
    end do
    end do 
    }

    call getbc(t,ps,psCoarse)
    
    !M1_TODO: check here for NAN in any var!
    {#IFDEF NAN_CHECK 
    write(78,*) "igrid",igrid
    do i=1,nw
      write(78,*) "i",i,pw(igrid)%w(ixG^T,i)
    end do
    }

    ! this was commented:
    !do iigrid=1,igridstail; igrid=igrids(iigrid)
    !  call Simple_check_data_correctness_pt(pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim), pw(igrid)%w(ixG^T,1:nw), 'b4 metric init')
    !enddo
    if (mype==0) then
      write(*, '(A,ES9.2,A)') ' Start initializating metric'
      write(*, '(A4,A10,A12,A12,A12)') '  #', 'mg_it', 'psi_err'
    end if

    if (init_psi) then
      do mg_it = 1, mg_it_max
  
         psi_err_old = psi_err
         
  
         !$OMP PARALLEL DO PRIVATE(igrid)
         do iigrid=1,igridstail; igrid=igrids(iigrid)
            ! copy psi
            pw1(igrid)%w(ixM^T, psi_metric_) = pw(igrid)%w(ixM^T, psi_metric_)
 
 
            !if (.not. use_lfac_profile_to_init) then
            !   ! at each iteration, keep vi instead of Wvi, so we need to map pwold of vi here
            !   {^C&pw(igrid)%w(ixM^T, u^C_) = pwold(igrid)%w(ixM^T, u^C_)  \}
            !   ! get new lfac with current psi, note that the input now is v instead of Wv,
            !   ! the following subroutine return 1+v^2
            !   call imhd_get_intermediate_variables(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, lfac2=lfac2(ixG^T,iigrid))
            !   lfac(ixM^T,iigrid) = 1.0d0 / dsqrt(2.0d0-lfac2(ixM^T,iigrid))
            !   ! store the current Wv in pw
            !   {^C&pw(igrid)%w(ixM^T, u^C_) = pw(igrid)%w(ixM^T, u^C_) &
            !                                  * lfac(ixM^T,iigrid)    \}
            !endif

            if (fix_prim_init) then
            ! it has very large initial pert, since u fix the prim, not cons
            ! everytime recover the v^i not Wv^i
            {^C& pw(igrid)%w(ixM^T,u^C_) = pwold(igrid)%w(ixM^T,u^C_) \}
                 call imhd_get_intermediate_variables(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, lfac2=lfac2(ixG^T,iigrid))
                 lfac(ixM^T,iigrid) = 1.0d0 / dsqrt(2.0d0-lfac2(ixM^T,iigrid))
            {^C& pw(igrid)%w(ixM^T, u^C_) = pw(igrid)%w(ixM^T, u^C_) &
                                            * lfac(ixM^T,iigrid)    \}
   
               call conserve(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim),patchfalse)
               ! every conserved quantity needed to times psi_metric_ **6
               call conformal_transformation(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))

               
               {#IFDEF NAN_CHECK
               write(79,*) "igrid",igrid
               do i=1,nw
                 write(79,*) "i",i,pw(igrid)%w(ixG^T,i)
               end do
               }
            !!me:
            !!   {do ix^D=ixMlo^D, ixMhi^D \}
            !!  ! if (pw(igrid)%w(ix^D, psi_metric_) .ne. pw(igrid)%w(ix^D, psi_metric_)) then
            !!     write(*,*) "after conformal trafo"
            !!     write(*,*) pw(igrid)%w(ixM^T, psi_metric_)
            !!     write(24,*) "after conformal trafo, ixD",ix^D
            !!     write(24,*) pw(igrid)%w(ixM^T, psi_metric_)
            !!  ! end if 
            !!   {enddo^D&\}
            !!!

            endif
  
         end do
         !$OMP END PARALLEL DO
  
         ! update the boundaries ensure the ghostzones also have conformal transed values
         call getbc(t,ps,psCoarse)

         {#IFDEF NAN_CHECK        
         do iigrid=1,igridstail; igrid=igrids(iigrid)
         write(80,*) "igrid",igrid
         do i=1,nw
           write(80,*) "i",i,pw(igrid)%w(ixG^T,i)
         end do
         end do
         }

         call cfc_solve_beta(.False.,Aij(ixG^T,1:3,1:3,1:igridstail))

         {#IFDEF NAN_CHECK
         write(81,*) "igrid"
         do i=1,nw
           write(81,*) "i",i,Aij(ixG^T,1:3,1:3,1:igridstail)
         end do
         }
  
         !$OMP PARALLEL DO PRIVATE(igrid)
         do iigrid = 1, igridstail;  igrid = igrids(iigrid);
            ! Geometry subroutines expect this to be set
            !block => ps(igrid)
            ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
            !call cfc_get_Aij_grid(ps(igrid), ixG^LL, ixM^LL, &
            !       A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
            call cfc_get_Aij_grid(pw(igrid)%w, px(igrid)%x, ixG^LL, ixM^LL, &
                   A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))


               {#IFDEF NAN_CHECK
                 write(82,*) "iigrid",iigrid
                 write(83,*) "iigrid",iigrid
                 do i=1,nw
                  write(82,*) "i",i,Aij(ixG^T,1:3,1:3,iigrid)
                  write(83,*) "i",i,A2(ixG^T,iigrid)
                 end do
               }
            

             !!me:
             !!      {do ix^D=ixMlo^D, ixMhi^D \}
             !!     ! if (pw(igrid)%w(ix^D, psi_metric_) .ne. pw(igrid)%w(ix^D, psi_metric_)) then
             !!        write(*,*) "after cfc get Aij"
             !!        write(*,*) pw(igrid)%w(ixM^T, psi_metric_)
             !!        write(23,*) "after cfc get Aij, ixD",ix^D
             !!        write(23,*) pw(igrid)%w(ixM^T, psi_metric_)
             !!     ! end if 
             !!      {enddo^D&\}
             !!!
 
         end do
         !$OMP END PARALLEL DO
         
         !M1_TODO
         !write(*,*)"before solve cfc psi"
         {#IFDEF NAN_CHECK
         do iigrid = 1, igridstail;  igrid = igrids(iigrid);
         write(85,*) "igrid",igrid
         do i=1,nw
          write(85,*) "i",i,pw(igrid)%w(ixG^T, i)
         end do
         end do 
         }
         
         call cfc_solve_psi(A2(ixG^T,1:igridstail))
         local_err = 0.0d0

         {#IFDEF NAN_CHECK
         do iigrid = 1, igridstail;  igrid = igrids(iigrid);
         write(84,*) "igrid",igrid
         do i=1,nw
          write(84,*) "i",i,pw(igrid)%w(ixG^T, i)
         end do
         end do 
          }

         !!!!!me      
        !!$OMP PARALLEL DO PRIVATE(igrid)
        !! do iigrid=1,igridstail; igrid=igrids(iigrid)      
        !!    write(22,*)"-----------------in cfc psi: igrid",igrid
        !!    {do ix^D=ixMlo^D, ixMhi^D \}
        !!    !if (pw(igrid)%w(ix^D, psi_metric_) .ne. pw(igrid)%w(ix^D, psi_metric_)) then
        !!        !write(*,*) pw(igrid)%w(ixM^T, psi_metric_)
        !!        write(22,*) "ixD---",ix^D
        !!        write(22,*) "psi",pw(igrid)%w(ix^D, psi_metric_)
        !!    ! end if
        !!    {enddo^D&\}
        !!    !!
        !! end do
        !!$OMP END PARALLEL DO
         !!!!!


         !$OMP PARALLEL DO PRIVATE(igrid)
         do iigrid=1,igridstail; igrid=igrids(iigrid)
            !!mme
            !write(*,*)"in cfc psi: igrid",igrid
            ! me: Simple_check_data_correctness was commented
           !call Simple_check_data_correctness_pt(pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim), pw(igrid)%w(ixG^T,1:nw), 'iterating psi')
           {do ix^D=ixMlo^D, ixMhi^D \}
              if (pw(igrid)%w(ix^D, psi_metric_) .ne. pw(igrid)%w(ix^D, psi_metric_)) then
                  write(*,*) pw(igrid)%w(ixM^T, psi_metric_)
                  write(*,*) pw(igrid)%w(ix^D, psi_metric_), px(igrid)%x(ix^D,:), 'psi ix^D x ix^D'
                  call mpistop('pw(igrid)%w(ix^D, psi_metric_) NaN when initialization')
              endif
              pw(igrid)%w(ix^D, psi_metric_) = max(pw(igrid)%w(ix^D, psi_metric_), 1.0d0)
           {enddo^D&\}

           local_err = max( local_err, maxval( &
                     dabs( (pw(igrid)%w(ixM^T, psi_metric_) - pw1(igrid)%w(ixM^T, psi_metric_) ) &
                           /pw1(igrid)%w(ixM^T, psi_metric_) ) ) )
         end do
         !$OMP END PARALLEL DO
         call mpi_allreduce(local_err, psi_err, 1, mpi_double, mpi_max, icomm, ierrmpi)
  
  
         if ( (psi_err <= cfc_psi_tol_init)  ) then
            if (mype==0) then
               write(*,*) "! psi initialization completed! "
               write(*,*) "! mg_it = ", mg_it, "psi_err = ", psi_err
               write(*,*) " "
               print*,'-------------------------------------------------------------------------------'
               write(*,*) " "
            end if
            exit
         else if ( psi_err >= psi_err_old ) then
            if (mype==0) then
               write(*,*) "! Warning: metric is not converging."
               write(*,*) "! mg_it = ", mg_it, "psi_err = ", psi_err
               write(*,*) "! Stop initializating here."
               !$OMP PARALLEL DO PRIVATE(igrid)
               do iigrid=1,igridstail; igrid=igrids(iigrid)
                  ! use last working psi
                  pw(igrid)%w(ixM^T, psi_metric_) = pw1(igrid)%w(ixM^T, psi_metric_)
               end do
               !$OMP END PARALLEL DO
            end if
            exit
         else if (mg_it >= mg_it_max) then
            if (mype==0) then
               write(*,*) "! Warning: metric failed to converge. mg_it = ", mg_it, "psi_err = ", psi_err
               write(*,*) "! Stop initializating here."
            end if
            exit
         end if
  
         if (mype==0) then
            write(*, '(A4,I10,ES12.3,ES12.3,ES12.3)') " #", &
                     mg_it, psi_err
         end if
  
      end do ! after psi loop
    endif

    ! psi trick
    if (psi_trick) then
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
          call inverse_conformal_trans(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
          call primitive(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))

          ! store the current v in pw
          {^C& pwold(igrid)%w(ixM^T, u^C_) = pw(igrid)%w(ixM^T, u^C_) &
                                             / pw(igrid)%w(ixM^T,lfac_)    \}
          ! Overshoot the psi
          pw(igrid)%w(ixM^T,psi_metric_) = pw(igrid)%w(ixM^T,psi_metric_) * psi_trick_value
          pw(igrid)%w(ixM^T,psi_metric_) = max(pw(igrid)%w(ixM^T,psi_metric_), 1.0d0)

          ! take v^i back
          {^C& pw(igrid)%w(ixM^T,u^C_) = pwold(igrid)%w(ixM^T,u^C_) \}
          call imhd_get_intermediate_variables(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, lfac2=lfac2(ixG^T,iigrid))
          lfac(ixM^T,iigrid) = 1.0d0 / dsqrt(2.0d0-lfac2(ixM^T,iigrid))
          {^C& pw(igrid)%w(ixM^T, u^C_) = pw(igrid)%w(ixM^T, u^C_) &
                                          * lfac(ixM^T,iigrid)    \}

          ! do p2c
          pw(igrid)%w(ixG^T,D_)       = 0.0d0
          if (eos_type == tabulated) pw(igrid)%w(ixG^T,Dye_)     = 0.0d0
          pw(igrid)%w(ixG^T,tau_)     = 0.0d0
          {^C& pw(igrid)%w(ixG^T,S^C_) = 0.0d0 \}
          !me test Bfield:
          !{C pw(igrid)%w(ixG^T,b^C_) = 0.0d0 }
          call conserve(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim),patchfalse)
          ! every conserved quantity needed to times psi_metric_ **6
          call conformal_transformation(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
          pwold(igrid)%w(ixM^T,1:nw) = pw(igrid)%w(ixM^T,1:nw)
      end do
      !$OMP END PARALLEL DO
      call getbc(t,ps,psCoarse)

          do mg_it = 1, mg_it_max
  
             psi_err_old = psi_err
             
  
             !$OMP PARALLEL DO PRIVATE(igrid)
             do iigrid=1,igridstail; igrid=igrids(iigrid)
                ! copy psi
                pw1(igrid)%w(ixM^T, psi_metric_) = pw(igrid)%w(ixM^T, psi_metric_)
  
             end do
             !$OMP END PARALLEL DO
  
             ! update the boundaries ensure the ghostzones also have conformal transed values
             call getbc(t,ps,psCoarse)
  
             call cfc_solve_beta(.False.,Aij(ixG^T,1:3,1:3,1:igridstail))
  
             !$OMP PARALLEL DO PRIVATE(igrid)
             do iigrid = 1, igridstail;  igrid = igrids(iigrid);
                ! Geometry subroutines expect this to be set
                !block => ps(igrid)
                ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
                !call cfc_get_Aij_grid(ps(igrid), ixG^LL, ixM^LL, &
                !       A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
                call cfc_get_Aij_grid(pw(igrid)%w, px(igrid)%x, ixG^LL, ixM^LL, &
                       A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
  
             end do
             !$OMP END PARALLEL DO
             call cfc_solve_psi(A2(ixG^T,1:igridstail))
  
             local_err = 0.0d0
             !$OMP PARALLEL DO PRIVATE(igrid)
             do iigrid=1,igridstail; igrid=igrids(iigrid)
                local_err = max( local_err, maxval( &
                          dabs( (pw(igrid)%w(ixM^T, psi_metric_) - pw1(igrid)%w(ixM^T, psi_metric_) ) &
                                /pw1(igrid)%w(ixM^T, psi_metric_) ) ) )
             end do
             !$OMP END PARALLEL DO
             call mpi_allreduce(local_err, psi_err, 1, mpi_double, mpi_max, icomm, ierrmpi)
  
  
             if ( (psi_err <= cfc_psi_tol_init)  ) then
                if (mype==0) then
                   write(*,*) "! psi initialization completed! "
                   write(*,*) "! mg_it = ", mg_it, "psi_err = ", psi_err
                   write(*,*) " "
                   print*,'-------------------------------------------------------------------------------'
                   write(*,*) " "
                end if
                exit
             else if ( psi_err >= psi_err_old ) then
                if (mype==0) then
                   write(*,*) "! Warning: metric is not converging."
                   write(*,*) "! mg_it = ", mg_it, "psi_err = ", psi_err
                   write(*,*) "! Stop initializating here."
                   !$OMP PARALLEL DO PRIVATE(igrid)
                   do iigrid=1,igridstail; igrid=igrids(iigrid)
                      ! use last working psi
                      pw(igrid)%w(ixM^T, psi_metric_) = pw1(igrid)%w(ixM^T, psi_metric_)
                   end do
                   !$OMP END PARALLEL DO
                end if
                exit
             else if (mg_it >= mg_it_max) then
                if (mype==0) then
                   write(*,*) "! Warning: metric failed to converge. mg_it = ", mg_it, "psi_err = ", psi_err
                   write(*,*) "! Stop initializating here."
                end if
                exit
             end if
  
             if (mype==0) then
                write(*, '(A4,I10,ES12.3,ES12.3,ES12.3)') " #", &
                         mg_it, psi_err
             end if
  
          end do ! after psi loop
          call getbc(t,ps,psCoarse)
          ! this trick may not converge, bound it
          do iigrid=1,igridstail; igrid=igrids(iigrid)
             pw(igrid)%w(ixG^T, psi_metric_) = max(pw(igrid)%w(ixG^T, psi_metric_), 1.0d0)
          end do
    endif  ! end of psi trick

    ! solve rest of the metric variables
    if (cfc_test_alp) then
       mg%timers(:)%t = 0.0d0 ! reset all timers
       do iigrid=1,igridstail; igrid=igrids(iigrid)
          ! flat space as an inital guess
          pw(igrid)%w(ixG^T, alp_metric_) = 1.0d0 
       end do
       if (mype==0) write(*,*) 'reset all timers, and set alpha = 1'
    end if


! this is to abandon "keeping imported prim to be absolutely unchanged in initialization" --> extremely large initial perturbation
! So, go for c2p with newest updated metric --> decrease initial pert
    if (.not. fix_prim_init) then
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ! save psi6 cons
    !      pwold(igrid)%w(ixG^T,D_)   = pw(igrid)%w(ixG^T,D_)   
    !      pwold(igrid)%w(ixG^T,Dye_) = pw(igrid)%w(ixG^T,Dye_)   
    !      pwold(igrid)%w(ixG^T,tau_) = pw(igrid)%w(ixG^T,tau_) 
    ! {^C& pwold(igrid)%w(ixG^T,S^C_) = pw(igrid)%w(ixG^T,S^C_)  \}

          call inverse_conformal_trans(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
          call primitive(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
          ! restore psi6 cons to find new alp beta
          pw(igrid)%w(ixG^T,D_)   = pwold(igrid)%w(ixG^T,D_)  
          pw(igrid)%w(ixG^T,tau_) = pwold(igrid)%w(ixG^T,tau_)
     {^C& pw(igrid)%w(ixG^T,S^C_) = pwold(igrid)%w(ixG^T,S^C_) \}
     {^C& pw(igrid)%w(ixG^T,b^C_) = pwold(igrid)%w(ixG^T,b^C_) \}
          if (eos_type == tabulated) pw(igrid)%w(ixG^T,Dye_) = pwold(igrid)%w(ixG^T,Dye_)
          {#IFDEF M1
          {^KSP&         
          {#IFNDEF VAR1
          pw(igrid)%w(ixG^T,nrad^KSP_) = pwold(igrid)%w(ixG^T,nrad^KSP_)
          pw(igrid)%w(ixG^T,erad^KSP_) = pwold(igrid)%w(ixG^T,erad^KSP_)
          {^C& pw(igrid)%w(ixG^T,frad^KSP^C_) = pwold(igrid)%w(ixG^T,frad^KSP^C_) \}
          }
          !pw(igrid)%w(ixG^T,frad^KSP1_) = pwold(igrid)%w(ixG^T,frad^KSP1_) 
          !pw(igrid)%w(ixG^T,frad^KSP2_) = pwold(igrid)%w(ixG^T,frad^KSP2_) 
          !pw(igrid)%w(ixG^T,frad^KSP3_) = pwold(igrid)%w(ixG^T,frad^KSP3_) 
          \}
          }
       end do
       !$OMP END PARALLEL DO
    endif

    call cfc_solve_alp(A2(ixG^T,1:igridstail))
    if (.not. cfc_test_alp) &
       call cfc_solve_beta(.True., Aij(ixG^T,1:3,1:3,1:igridstail))

    ! once all the metric variables are solved, Aij is not needed anymore
    deallocate(Aij)
    deallocate(A2)

    psi_xns_err = 0.0d0
    psi_xns_err_local = 0.0d0
    alp_xns_err = 0.0d0
    alp_xns_err_local = 0.0d0
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid)
       psi_xns_err_local = maxval(dabs((pw(igrid)%w(ixM^T, psi_metric_) - pwold(igrid)%w(ixM^T, psi_metric_) ) &
                       /pwold(igrid)%w(ixM^T, psi_metric_) ))
       alp_xns_err_local = maxval(dabs((pw(igrid)%w(ixM^T, alp_metric_) - pwold(igrid)%w(ixM^T, alp_metric_) ) &
                       /pwold(igrid)%w(ixM^T, alp_metric_) ))
    end do
    !$OMP END PARALLEL DO

    call mpi_allreduce(psi_xns_err_local, psi_xns_err, 1, mpi_double, mpi_max, icomm, ierrmpi)
    call mpi_allreduce(alp_xns_err_local, alp_xns_err, 1, mpi_double, mpi_max, icomm, ierrmpi)

    if (mype==0) then
       write(*,*) "! metric initialization completed! "
       write(*,*) " "
       write(*,*) " psi_xns_err = ", psi_xns_err
       write(*,*) " alp_xns_err = ", alp_xns_err
       write(*,*) " "
       print*,'-------------------------------------------------------------------------------'
       write(*,*) " "
    end if


   ! ! clear up cons 
   ! !$OMP PARALLEL DO PRIVATE(igrid)
   ! do iigrid=1,igridstail; igrid=igrids(iigrid);
   !    pw(igrid)%w(ixM^T, 1:ncons) = 0.0d0
   ! end do
   ! !$OMP END PARALLEL DO

    ! recover all original hydro vars (no metric vars)

!then temp changed largely.
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);

!       if (.not. initialize_with_cons) then
!         ! restore the hydro variables only prim!! 
        if (fix_prim_init) then
          pw(igrid)%w(ixM^T,nprim_lo:nprim_hi) = pwold(igrid)%w(ixM^T,nprim_lo:nprim_hi)
     {^C& pw(igrid)%w(ixM^T,b^C_)              = pwold(igrid)%w(ixM^T,b^C_)   \}
           {#IFDEF M1
          {^KSP&         
          {#IFNDEF VAR1
          pw(igrid)%w(ixM^T,nrad^KSP_) = pwold(igrid)%w(ixM^T,nrad^KSP_)
          pw(igrid)%w(ixM^T,erad^KSP_) = pwold(igrid)%w(ixM^T,erad^KSP_)
          {^C& pw(igrid)%w(ixM^T,frad^KSP^C_) = pwold(igrid)%w(ixM^T,frad^KSP^C_) \}
          }
          !pw(igrid)%w(ixG^T,frad^KSP1_) = pwold(igrid)%w(ixG^T,frad^KSP1_) 
          !pw(igrid)%w(ixG^T,frad^KSP2_) = pwold(igrid)%w(ixG^T,frad^KSP2_) 
          !pw(igrid)%w(ixG^T,frad^KSP3_) = pwold(igrid)%w(ixG^T,frad^KSP3_) 
          \}
          }
          ! the pwold(u) is v^i only, we need to find Wv
          call imhd_get_intermediate_variables(ixG^LL, ixM^LL, pw(igrid)%w, px(igrid)%x, lfac2=lfac2(ixG^T,iigrid))
          lfac(ixM^T,iigrid) = 1.0d0 / dsqrt(2.0d0-lfac2(ixM^T,iigrid))
          ! store the current Wv in ps
     {^C& pw(igrid)%w(ixM^T, u^C_) = pw(igrid)%w(ixM^T, u^C_) * lfac(ixM^T,iigrid)   \}
        endif

        ! with the current prim to recover corresponding cons without using psi6 cons
        call conserve(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim),patchfalse)
    end do
    !$OMP END PARALLEL DO

   deallocate(lfac)
   deallocate(lfac2)

  end subroutine cfc_metric_init

  subroutine cfc_check_metric(update)
    use mod_multigrid_coupling
    use mod_forest
    use mod_imhd_intermediate, only: conformal_transformation, inverse_conformal_trans
    include 'amrvacdef.f'

    logical, intent(out)           :: update
    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl
    type(tree_node), pointer       :: pnode
    double precision               :: res = 0.0d0

    double precision, allocatable  :: Aij(:^D&,:,:,:)
    double precision, allocatable  :: A2(:^D&,:)

    update = .False.

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid =  igrids(iigrid);
       pw1(igrid)%w(ixG^T, 1:nw) = pw(igrid)%w(ixG^T,1:nw)
       !pwold(igrid)%w(ixG^T, 1:nw) = pw(igrid)%w(ixG^T,1:nw)
       call conformal_transformation(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
    end do
    !$OMP END PARALLEL DO

    allocate(Aij(ixG^T,1:3,1:3,1:igridstail))
    allocate(A2(ixG^T,1:igridstail))
    Aij = 0.0d0
    A2 = 0.0d0
    ! fixme: should I solve vecX for metric check?
    !call cfc_solve_beta(.False.,Aij(ixG^T,1:3,1:3,1:igridstail))
    ! obtain the conformal extrinsic curvature
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid =  igrids(iigrid);
       ! Geometry subroutines expect this to be set
       call set_tmpGlobals(igrid)

       !call cfc_get_Aij_grid(ps(igrid), ixG^LL, ixM^LL, &
       !          A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
       call cfc_get_Aij_grid(pw(igrid)%w, px(igrid)%x, ixG^LL, ixM^LL, &
                 A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
    end do
    !$OMP END PARALLEL DO

    call cfc_check_psi(res, A2(ixG^T,1:igridstail))
    !if (mype==0) write(*,*) 'res of psi:', res
    if ( res >= cfc_psi_tol_max ) then
       update = .true.
    end if
    deallocate(Aij)
    deallocate(A2)

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid =  igrids(iigrid);
       !call inverse_conformal_trans(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
       pw(igrid)%w(ixG^T, 1:nw) = pw1(igrid)%w(ixG^T,1:nw)
       !pw(igrid)%w(ixG^T, 1:nw) = pwold(igrid)%w(ixG^T,1:nw)

    end do
    !$OMP END PARALLEL DO
  end subroutine cfc_check_metric

  subroutine cfc_solve_metric()
    use mod_multigrid_coupling
    use mod_forest
    use mod_imhd_intermediate, only: conformal_transformation, inverse_conformal_trans
    include 'amrvacdef.f'

    integer                        :: id, iigrid, igrid
    integer                        :: nc, lvl
    type(tree_node), pointer       :: pnode

    double precision, allocatable  :: Aij(:^D&,:,:,:)
    double precision, allocatable  :: A2(:^D&,:)
    double precision               :: array_max_local(1:4), array_max(1:4), alp_min_local, alp_min
    double precision :: max_local(10), max_global(10)
    double precision :: min_local(1), min_global(1)

    ! b4 here, shared ghostzones for both nwflux + aux
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid =  igrids(iigrid);

       pw1(igrid)%w(ixG^T,psi_metric_) = pw(igrid)%w(ixG^T,psi_metric_)
       pw1(igrid)%w(ixG^T,D_)          = pw(igrid)%w(ixG^T,D_)
       if (eos_type == tabulated) pw1(igrid)%w(ixG^T,Dye_)        = pw(igrid)%w(ixG^T,Dye_)
       pw1(igrid)%w(ixG^T,tau_)        = pw(igrid)%w(ixG^T,tau_)
  {^C& pw1(igrid)%w(ixG^T,S^C_)        = pw(igrid)%w(ixG^T,S^C_) \}
  {^C& pw1(igrid)%w(ixG^T,b^C_)        = pw(igrid)%w(ixG^T,b^C_) \}
     {#IFDEF M1 
      {^KSP&
      {#IFNDEF VAR1
      pw1(igrid)%w(ixG^T,nrad^KSP_)        = pw(igrid)%w(ixG^T,nrad^KSP_)
      pw1(igrid)%w(ixG^T,erad^KSP_)        = pw(igrid)%w(ixG^T,erad^KSP_)
      {^C& pw1(igrid)%w(ixG^T,frad^KSP^C_)        = pw(igrid)%w(ixG^T,frad^KSP^C_) \}
      }
      !!M1_TODO theoretically M1 coming in here in cons already (nothing would happen in conformal-trafo)
      \}
     }
       call conformal_transformation(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
    end do
    !$OMP END PARALLEL DO

    allocate(Aij(ixG^T,1:3,1:3,1:igridstail))
    allocate(A2(ixG^T,1:igridstail))
    Aij = 0.0d0
    A2 = 0.0d0

    ! initialize the solver
    if (cfc_redblack) then
       mg%smoother_type = mg_smoother_gsrb
    else
       mg%smoother_type = mg_smoother_gs
    end if

    mg%n_cycle_up   = cfc_n_cycle(1)
    mg%n_cycle_down = cfc_n_cycle(2)

    !-----------------------------------------------------------------------
    ! Step 1: Solve the PDE for X and thus the conformal extrinsic curvature
    !-----------------------------------------------------------------------
    call cfc_solve_beta(.False.,Aij(ixG^T,1:3,1:3,1:igridstail))

    ! obtain the conformal extrinsic curvature
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid = 1, igridstail;  igrid =  igrids(iigrid);
       ! Geometry subroutines expect this to be set
       !block => ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       !call cfc_get_Aij_grid(ps(igrid), ixG^LL, ixM^LL, &
       !          A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
       call cfc_get_Aij_grid(pw(igrid)%w, px(igrid)%x, ixG^LL, ixM^LL, &
                 A2(ixG^T,iigrid), Aij(ixG^T,1:3,1:3,iigrid))
    end do
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    ! Step 2: After having the extrinsic curvature, we can solve psi
    !-----------------------------------------------------------------------
    call cfc_solve_psi(A2(ixG^T,1:igridstail))

    !-----------------------------------------------------------------------
    ! Step 3: Before solving the alp, we need to know S_star first 
    !             (That means update the prim variables)
    !-----------------------------------------------------------------------
    !  before c2p, you need to recover prim for solving alp 
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! for safety:
       pw(igrid)%w(ixG^T,psi_metric_) = max(1.0d0,pw(igrid)%w(ixG^T,psi_metric_))
       ! psi has ghostzones here already

       !using new psi to do c2p
       call inverse_conformal_trans(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
       call primitive(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))

!fixme, semi-cowling way not working
!          pw(igrid)%w(ixG^T,D_)    = pwold(igrid)%w(ixG^T,D_)
!          pw(igrid)%w(ixG^T,Dye_)  = pwold(igrid)%w(ixG^T,Dye_)
!          pw(igrid)%w(ixG^T,tau_)  = pwold(igrid)%w(ixG^T,tau_)
!     {^C& pw(igrid)%w(ixG^T,S^C_)  = pwold(igrid)%w(ixG^T,S^C_) \}
!       call primitive(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim))
       ! restore the original cons with psi6 to solve alp beta
       
    end do
    !$OMP END PARALLEL DO
   
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pw(igrid)%w(ixG^T,D_)    = pw1(igrid)%w(ixG^T,D_)   * pw1(igrid)%w(ixG^T,psi_metric_)**6
       if (eos_type == tabulated) pw(igrid)%w(ixG^T,Dye_)  = pw1(igrid)%w(ixG^T,Dye_) * pw1(igrid)%w(ixG^T,psi_metric_)**6
       pw(igrid)%w(ixG^T,tau_)  = pw1(igrid)%w(ixG^T,tau_) * pw1(igrid)%w(ixG^T,psi_metric_)**6
  {^C& pw(igrid)%w(ixG^T,S^C_)  = pw1(igrid)%w(ixG^T,S^C_) * pw1(igrid)%w(ixG^T,psi_metric_)**6 \}
  {^C& pw(igrid)%w(ixG^T,b^C_)  = pw1(igrid)%w(ixG^T,b^C_) * pw1(igrid)%w(ixG^T,psi_metric_)**6 \}
    {#IFDEF M1    
     {^KSP&
     ! TODO actually without psi6 because is cons all time?
       {#IFNDEF VAR1
       pw(igrid)%w(ixG^T,nrad^KSP_)  = pw1(igrid)%w(ixG^T,nrad^KSP_) !* pw1(igrid)%w(ixG^T,psi_metric_)**6
       pw(igrid)%w(ixG^T,erad^KSP_)  = pw1(igrid)%w(ixG^T,erad^KSP_) !* pw1(igrid)%w(ixG^T,psi_metric_)**6
       {^C& pw(igrid)%w(ixG^T,frad^KSP^C_)  = pw1(igrid)%w(ixG^T,frad^KSP^C_)  \} !* pw1(igrid)%w(ixG^T,psi_metric_)**6
       }
      \}
    }
    end do
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    ! Step 4: Solve alp
    !-----------------------------------------------------------------------
    call cfc_solve_alp(A2(ixG^T,1:igridstail))
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! for safety:
       pw(igrid)%w(ixG^T,alp_metric_) = min(1.0d0,pw(igrid)%w(ixG^T,alp_metric_))
    end do
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    ! Step 5: Solve shift vector beta
    !-----------------------------------------------------------------------
    call cfc_solve_beta(.True., Aij(ixG^T,1:3,1:3,1:igridstail))

    ! U is currently psi6 U, --> change it back to U, but we should not use inverse conformal trans, 
    ! as we need to consistent with prim as
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call inverse_conformal_trans(ixG^LL, ixM^LL, pw(igrid)%w(ixG^T,1:nw), px(igrid)%x(ixG^T,1:ndim))
       call set_tmpGlobals(igrid)
       call conserve(ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,1:nw),px(igrid)%x(ixG^T,1:ndim),patchfalse)
    end do
    !$OMP END PARALLEL DO

    
    !do iigrid=1,igridstail; igrid=igrids(iigrid);
    !  write(333,*) "igrid,alp",igrid,pw(igrid)%w(ixG^T,alp_metric_)
    !  write(333,*) "igrid,psi",igrid,pw(igrid)%w(ixG^T,psi_metric_)
    !end do 
    ! once all the metric variables are solved, Aij is not needed anymore
    deallocate(Aij)
    deallocate(A2)

  end subroutine cfc_solve_metric

  !> calculate Aij on a single grid
  subroutine cfc_get_Aij_grid( w, x, ixI^L, ixO^L, A2, Aij)
    use mod_metric
    include 'amrvacdef.f'
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(out)   :: A2(ixI^S)
    double precision, intent(out)   :: Aij(ixI^S,1:3,1:3)

    integer                         :: jxO^L, hxO^L
    integer                         :: idir, jdir, kdir, ldir
    double precision                :: fij(ixI^S, 1:3, 1:3)
    double precision                :: dvecX(ixI^S,1:3,1:3)
    {^NOONED
    double precision                :: cot_theta(ixI^S)
    }


    call get_gamma_ij_hat(x(ixI^S,1:ndim),ixI^L,ixO^L,fij(ixI^S, 1:3, 1:3))

    dvecX = 0.0d0
    do idir = 1, ndim
       {^C&
          call partial_d( w(ixI^S,Xvec_metric^C_), ixI^L, ixO^L, x(ixI^S,1:ndim), idir, dvecX(ixI^S, ^C, idir))
       \}
    end do

    Aij = 0.0d0
    select case (coordinate)
    case (cartesian)
       Aij(ixO^S, 1,1) = 2.0d0/3.0d0 * ( 2.0d0 * dvecX(ixO^S,1,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } & 
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 2,2) = 2.0d0/3.0d0 * ( - dvecX(ixO^S,1,1) &
                         {^NOONED + 2.0d0 * dvecX(ixO^S,2,2) } &
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 3,3) = 2.0d0/3.0d0 * ( - dvecX(ixO^S,1,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } &
                         {^IFTHREED + 2.0d0 * dvecX(ixO^S,3,3) } )
    case (cylindrical)
       Aij(ixO^S, 1,1) = 2.0d0/3.0d0 * ( 2.0d0 * dvecX(ixO^S,1,1) - w(ixO^S,Xvec_metric1_)/x(ixO^S,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } & 
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 2,2) = 2.0d0/3.0d0 * ( -dvecX(ixO^S,1,1) - w(ixO^S,Xvec_metric1_)/x(ixO^S,1) &
                         {^NOONED + 2.0d0 * dvecX(ixO^S,2,2) } &
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 3,3) = 2.0d0/3.0d0/fij(ixO^S,3,3) * ( -dvecX(ixO^S,1,1) + 2.0d0 * w(ixO^S,Xvec_metric1_)/x(ixO^S,1) &
                         {^NOONED - dvecX(ixO^S,2,2) } &
                         {^IFTHREED + 2.0d0 * dvecX(ixO^S,3,3) } )
    case (spherical)
       {^NOONED   cot_theta(ixO^S) = dcos(x(ixO^S,2))/dsin(x(ixO^S,2))    }
       Aij(ixO^S, 1,1) = 2.0d0/3.0d0 * ( 2.0d0 * ( dvecX(ixO^S,1,1) - w(ixO^S,Xvec_metric1_)/x(ixO^S,1) )&
                         {^NOONED - dvecX(ixO^S,2,2) - cot_theta(ixO^S) * w(ixO^S,Xvec_metric2_) } & 
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 2,2) = 2.0d0/3.0d0/fij(ixO^S,2,2) * ( -dvecX(ixO^S,1,1) + w(ixO^S,Xvec_metric1_)/x(ixO^S,1) &
                         {^NOONED + 2.0d0 * dvecX(ixO^S,2,2) - cot_theta(ixO^S) * w(ixO^S,Xvec_metric2_) } &
                         {^IFTHREED - dvecX(ixO^S,3,3) } )
       Aij(ixO^S, 3,3) = 2.0d0/3.0d0/fij(ixO^S,3,3) * ( -dvecX(ixO^S,1,1) + w(ixO^S,Xvec_metric1_)/x(ixO^S,1) &
                         {^NOONED - dvecX(ixO^S,2,2) + 2.0d0 * cot_theta(ixO^S) * w(ixO^S,Xvec_metric2_) } &
                         {^IFTHREED + 2.0d0 * dvecX(ixO^S,3,3) } )
   
    end select
    {^NOONED
    Aij(ixO^S, 1,2) =  dvecX(ixO^S,1,2)/fij(ixO^S,2,2) + dvecX(ixO^S,2,1) 
    Aij(ixO^S, 1,3) =  dvecX(ixO^S,3,1)                + dvecX(ixO^S,1,3)/fij(ixO^S,3,3)
    Aij(ixO^S, 2,3) =  dvecX(ixO^S,3,2)/fij(ixO^S,2,2) + dvecX(ixO^S,2,3)/fij(ixO^S,3,3)
    do idir = 1, 2
       do jdir = idir+1, 3
          Aij(ixO^S, jdir, idir) = Aij(ixO^S, idir, jdir) 
       end do
    end do
    }

    ! get A2
    !fixme:
    A2=0.0d0
    do jdir = 1, 3; do idir = 1, 3
       do kdir = 1, 3
          do ldir = 1, 3
             A2(ixO^S) = A2(ixO^S) &
              + fij(ixO^S,idir,kdir) * fij(ixO^S,jdir,ldir) &
                  * Aij(ixO^S,idir,jdir) * Aij(ixO^S,kdir,ldir) 
          end do
       end do
    end do; end do

  end subroutine cfc_get_Aij_grid

end module mod_cfc


