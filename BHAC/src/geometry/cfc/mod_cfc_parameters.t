module mod_cfc_parameters

  public 

  logical                                   :: use_cfc = .False.
  logical                                   :: use_multigrid = .false.
  logical                                   :: initialize_metric = .false.
  logical                                   :: restart_init_metric = .false.

  integer                                   :: coordinate         =-1
  integer, parameter                        :: Cartesian          = 0
  integer, parameter                        :: Cartesian_stretched= 1
  integer, parameter                        :: cylindrical        = 2
  integer, parameter                        :: spherical          = 3
  
  !-------------------------------------------------------------------!
  !          solver related parameters
  !-------------------------------------------------------------------!
  !> tolerence for initializing psi
  double precision                          ::  cfc_psi_tol_init = 1.0d-8
  
  !> provide conservative and conformal factor for initialization
  logical                                   ::  initialize_with_cons = .false.
  !> fix prim to initialize 
  logical                                   ::  fix_prim_init    = .false.
  !> Initialize psi
  logical                                   ::  init_psi         = .true.
  !> evolve cfc every rk step
  logical                                   ::  evolve_every_rk_step = .false.
  !> Dirty trick to consistent with Fil psi
  logical                                   ::  psi_trick        = .false.
  double precision                          ::  psi_trick_value  = 1.0005d0
  !> Activate checking psi every time step
  logical                                   ::  use_check_psi    = .true.

  !> tolerence for 1: alp, 2: psi, 3: beta/X
  double precision, dimension(1:3)          ::  cfc_tol = 1.0d-6
  !> tolerence for 1: alp, 2: psi, 3: beta/X when evolution
  double precision, dimension(1:3)          ::  cfc_tol_evolve = 1.0d-6
  !> maximun iteration for each metric solver
  integer                                   ::  cfc_it_max = 1000
  !> Print out the status of the solver at every cfc_print cfc iteration
  integer                                   ::  cfc_print = 10000000
  !> N smooths will be applied for 1: cycling up, 2: cycling down
  integer, dimension(1:2)                   ::  cfc_n_cycle = (/2,2/)
  logical                                   ::  cfc_redblack = .True.

  !> outer boundary of beta, either robin (default) or dirichlet
  ! If cartesian, it should be .false.
  logical                                   ::  cfc_beta_robin_bc = .True.
  !> When initialize metric, you use the profile lorentz factor
  logical                                   ::  use_lfac_profile_to_init = .true.

  !-------------------------------------------------------------------!
  !          update related parameters
  !-------------------------------------------------------------------!

  !> evolve metric or not
  logical                                   ::  cfc_evolve = .false.
  !> last metric update time
  double precision                          ::  cfc_t_last_update = 0.0d0
  integer                                   ::  cfc_it_last_update = 0
  !> time of different stages that will change the freq. of update of cfc solver
  double precision                          ::  cfc_t_end_of_stage1 = 1.0d+99
  double precision                          ::  cfc_t_end_of_stage2 = 1.0d+99

  !> solve the metric at every N steps
  integer                                   ::  cfc_dit_update = 1000000000
  integer                                   ::  cfc_dit_update_stage2 = 1000000000
  integer                                   ::  cfc_dit_update_stage3 = 1000000000
  !> allowed smallest time step for metric solver (for dit_update only)
  double precision                          ::  cfc_smallest_dt = 0.0d0
  !> solve the metric at every dt
  double precision                          ::  cfc_dt_update = 1.0d+99
  double precision                          ::  cfc_dt_update_stage2 = 1.0d+99
  double precision                          ::  cfc_dt_update_stage3 = 1.0d+99 
  !> Maximum tolerence for psi
  double precision                          ::  cfc_psi_tol_max = 1.0d-3

  !-------------------------------------------------------------------!
  !          interpolation related parameters
  !-------------------------------------------------------------------!
  !> number of points used to interpolate metric at the cell interface
  integer                                   ::  cfc_n_interpolation = 4
  !> Instead of using interpolation, metric variables can also be reconstructed by limiters, currently only wenozp5 is used.
  logical                                   ::  reconstruct_cfc = .False.

  !-------------------------------------------------------------------!
  !          debug related parameters
  !-------------------------------------------------------------------!
  !> Debug flag, test the performance of the metric solver 
  logical                                   ::  cfc_test_alp = .False.

 
  !-------------------------------------------------------------------!
  !          GW backreaction related parameters
  !-------------------------------------------------------------------!
  logical                                   ::  use_gw_br = .False.

  !> tolerence for 1: U_br, 2: U_i, 3: R
  double precision, dimension(1:3)          ::  gw_br_tol = 1.0d-6
  !> tolerence for 1: U_br, 2: U_i, 3: R when evolution
  double precision, dimension(1:3)          ::  gw_br_tol_evolve = 1.0d-6

  logical                                   ::  gw_br_include_dU_i     = .false.
  logical                                   ::  gw_br_use_I3_ij        = .false.
  logical                                   ::  gw_br_use_wi           = .true.
  logical                                   ::  use_h00_coupling       = .true.
  logical                                   ::  use_index_contract     = .false.
  logical                                   ::  use_h00_cfc_betai      = .false.
  logical                                   ::  use_h00_cfc_alp        = .false.
  logical                                   ::  gw_br_couple_weakly    = .false.

  !> Directly modify the metric alp and psi values
  logical                                   ::  use_3rd_d_Iij              = .false.
  !> maximum tolerance of relative changes between dt and dt_old to avoid unphysical Q3_ij
  double precision, parameter               ::  tol_rel_dt          = 10.0d0 ! 10 = 1000%
  double precision                          ::  gw_br_dt_update     = 1.0d+99
  double precision                          ::  gw_br_t_last_update = 0.0d0
  double precision                          ::  gw_br_t_last_time   = 0.0d0
  double precision                          ::  dt1                 = 0.0d0
  double precision                          ::  dt2                 = 0.0d0
  double precision                          ::  dt3                 = 0.0d0

end module mod_cfc_parameters
