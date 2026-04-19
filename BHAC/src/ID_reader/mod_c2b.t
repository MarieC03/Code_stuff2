!only has 2d
module mod_c2b

  implicit none
  private

  type c2b_profile
     integer                       :: Nr, Nth, Ntot, Nphi
     double precision, allocatable :: radius(:), theta(:), theta_r(:,:)
     double precision, allocatable :: beta1(:,:), beta2(:,:), beta3(:,:)
     double precision, allocatable :: alp(:,:), psi(:,:)
     !prim
     double precision, allocatable :: rho(:,:), press(:,:), v1(:,:), v2(:,:), v3(:,:)
     double precision, allocatable :: lfac(:,:), ye(:,:), temp(:,:), eps(:,:)
     !cons
     double precision, allocatable :: D(:,:), Dye(:,:), S_1(:,:), S_2(:,:), S_3(:,:), tau(:,:)
     double precision, allocatable :: sqrtgamma(:,:)
  end type c2b_profile

  ! save the profle
  type(c2b_profile) :: c2b_prof

  public :: c2b_profile
  public :: c2b_prof
  public :: mod_c2b_allocate_var
  public :: mod_c2b_deallocate_var
  public :: mod_c2b_read_profile
  public :: mod_c2b_map_2D
  public :: mod_c2b_map_1D

  integer, parameter :: C2B_ghostcells = 11

contains

subroutine mod_c2b_allocate_var(Nr, Nth)
  integer, intent(in) :: Nth, Nr
  allocate(c2b_prof%radius(Nr))
  allocate(c2b_prof%theta(Nth))
  allocate(c2b_prof%theta_r(Nr,Nth))
  !metric
  allocate(c2b_prof%alp(Nr,Nth),c2b_prof%psi(Nr,Nth))
  allocate(c2b_prof%beta1(Nr,Nth),c2b_prof%beta2(Nr,Nth),c2b_prof%beta3(Nr,Nth))
  allocate(c2b_prof%sqrtgamma(Nr,Nth))
  !hydro prim
  allocate(c2b_prof%rho(Nr,Nth),c2b_prof%press(Nr,Nth))
  allocate(c2b_prof%temp(Nr,Nth),c2b_prof%ye(Nr,Nth))
  allocate(c2b_prof%lfac(Nr,Nth),c2b_prof%eps(Nr,Nth))
  allocate(c2b_prof%v1(Nr,Nth),c2b_prof%v2(Nr,Nth),c2b_prof%v3(Nr,Nth))
  !hydro cons
  allocate(c2b_prof%D(Nr,Nth))
  allocate(c2b_prof%Dye(Nr,Nth),c2b_prof%tau(Nr,Nth))
  allocate(c2b_prof%S_1(Nr,Nth),c2b_prof%S_2(Nr,Nth),c2b_prof%S_3(Nr,Nth))
  
end subroutine mod_c2b_allocate_var

subroutine mod_c2b_deallocate_var
  deallocate(c2b_prof%radius)
  deallocate(c2b_prof%theta)
  deallocate(c2b_prof%theta_r)
  deallocate(c2b_prof%alp,c2b_prof%psi)
  deallocate(c2b_prof%beta1,c2b_prof%beta2,c2b_prof%beta3)
  deallocate(c2b_prof%sqrtgamma)

  deallocate(c2b_prof%rho,c2b_prof%press)
  deallocate(c2b_prof%temp,c2b_prof%ye)
  deallocate(c2b_prof%lfac,c2b_prof%eps)
  deallocate(c2b_prof%v1,c2b_prof%v2,c2b_prof%v3)

  deallocate(c2b_prof%D)
  deallocate(c2b_prof%Dye,c2b_prof%tau)
  deallocate(c2b_prof%S_1,c2b_prof%S_2,c2b_prof%S_3)
end subroutine mod_c2b_deallocate_var

subroutine mod_c2b_read_profile(read_C2B_metric, hydro_type, profile_path_hydro, profile_path_met)

  character*(*)           :: profile_path_met
  character*(*)           :: profile_path_hydro
  character*(*)           :: hydro_type
  character*128 :: filename
  logical, intent(in)     :: read_C2B_metric
  double precision        :: minrho
  double precision        :: x_star, x, x_minus, y, y_minus, maxrho

  integer :: Nth, Nr, Nphi, Ntot
  integer :: i,j,k
  call mod_c2b_get_dim(profile_path_hydro, Nr, Nth, Nphi, Ntot)

  if (Nphi > 1) call mpistop('This numerical.met is not 2D, is 3D')

  c2b_prof%Nr   = Nr
  c2b_prof%Nth  = Nth
  c2b_prof%Ntot = Ntot
  call mod_c2b_allocate_var(Nr, Nth)

  ! read profile
  call mod_c2b_get_grid(profile_path_hydro, Nr, Nth, c2b_prof%radius, c2b_prof%theta_r)

  ! read metric first, since psi will be replaced later  
  if (read_C2B_metric) then
       call mod_c2b_get_met_profile(profile_path_met, Nr, Nth, &
                                  c2b_prof%alp,c2b_prof%psi,&
                                  c2b_prof%beta1,c2b_prof%beta2,c2b_prof%beta3)
  endif

  select case(hydro_type)
  case('prim')
       call mod_c2b_get_prim_profile(profile_path_hydro, Nr, Nth, &
                                     c2b_prof%rho,c2b_prof%press,&
                                     c2b_prof%temp,c2b_prof%ye,&
                                     c2b_prof%lfac,c2b_prof%eps,&
                                     c2b_prof%v1,c2b_prof%v2,c2b_prof%v3, &
                                     c2b_prof%sqrtgamma, c2b_prof%psi)
  case('cons')
  ! read psi again from the conformal factor, will replace the psi from get_met_profile
       call mod_c2b_get_cons_profile(profile_path_hydro, Nr, Nth, &
                                     c2b_prof%rho,c2b_prof%press,&
                                     c2b_prof%temp,c2b_prof%ye,&
                                     c2b_prof%lfac,c2b_prof%eps,&
                                     c2b_prof%v1,c2b_prof%v2,c2b_prof%v3, &
                                     c2b_prof%D,c2b_prof%Dye,c2b_prof%tau, &
                                     c2b_prof%S_1,c2b_prof%S_2,c2b_prof%S_3, &
                                     c2b_prof%sqrtgamma, c2b_prof%psi)
  case default
     call mpistop('C2B reader: specify a hydro type of the fluidID.csv')
  end select

  do j = 1, Nth
    do i = 1, Nr
      c2b_prof%theta(j) = c2b_prof%theta_r(i,j)
!      maxrho = max(maxval(c2b_prof%rho),maxrho)
    enddo
  enddo
!write(*,*) maxrho
!stop 'maxrho'


end subroutine mod_c2b_read_profile

subroutine mod_c2b_get_dim(lprofile_name, Nr, Nth, Nphi, Ntot)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(out) :: Nth, Nr, Nphi, Ntot
  
  filename = trim(adjustl(lprofile_name))//"/fluidID.csv"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) 
  read(666,*) Ntot, Nr, Nth, Nphi

  close(666)

end subroutine mod_c2b_get_dim

subroutine mod_c2b_get_grid(lprofile_name, Nr, Nth, radius, theta_r)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(in) :: Nth, Nr
  double precision, intent(out) :: radius(Nr), theta_r(Nr,Nth)
  
  integer i,j
  
  filename = trim(adjustl(lprofile_name))//"/fluidID.csv"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) 
  read(666,*) 
  read(666,*) 
  do j = 1, Nth
    do i = 1, Nr
      read(666,*) radius(i), theta_r(i,j)
      !write(*,*) radius(i), theta_r(i,j), i, j , 'radius, theta_r, i, j '
    enddo
  enddo
 
  close(666)

end subroutine mod_c2b_get_grid

subroutine mod_c2b_get_met_profile(lprofile_name, Nr, Nth, alp, psi, beta1, beta2, beta3)

  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(in) :: Nth, Nr
  double precision, intent(out) :: psi(Nr,Nth),alp(Nr,Nth),beta1(Nr,Nth),beta2(Nr,Nth),beta3(Nr,Nth)

  integer i,j
  
  double precision :: discrim, rho_min
  double precision :: dx,dtheta,dphi, phi_max
  double precision :: dummy
  
  !--------------------------------------------------------!       
  ! read Hydroeq.dat profile      
  !--------------------------------------------------------!       
  filename = trim(adjustl(lprofile_name))//"/numerical.met"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) 
  read(666,*) 
  read(666,*) 
  read(666,*) 
  read(666,*) 

  do j=1,Nth
     do i=1,Nr
        read(666,*) dummy, dummy, alp(i,j), beta1(i,j), beta2(i,j), &
                    beta3(i,j), psi(i,j)
   !write(*,*) alp(i,j), beta1(i,j), beta2(i,j), &
   !                 beta3(i,j), psi(i,j)

     enddo
  enddo

  close(666)

  do j=1,Nth
     do i=1,Nr
       psi(i,j) = (psi(i,j))**(1.0d0/4.0d0)
       !write(*,*) psi(i,j), i, j, 'psi  i j '
     enddo
  enddo
end subroutine mod_c2b_get_met_profile

subroutine mod_c2b_get_prim_profile(lprofile_name, Nr, Nth, rho, press,&
                                     temp, ye,&
                                     lfac, eps,&
                                     v1, v2, v3, sqrtgamma, psi)


  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(in) :: Nth, Nr
  double precision, dimension(Nr,Nth), intent(out) :: rho, press, temp, ye, lfac, eps, v1, v2, v3
  double precision, dimension(Nr,Nth), intent(out) :: sqrtgamma, psi

  integer i,j
  
  double precision :: discrim, rho_min
  double precision :: dx,dtheta,dphi, phi_max
  double precision :: dummy
  
  !--------------------------------------------------------!       
  ! read Hydroeq.dat profile      
  !--------------------------------------------------------!       
  filename = trim(adjustl(lprofile_name))//"/fluidID.csv"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) 
  read(666,*) 
  read(666,*) 

  do j=1,Nth
     do i=1,Nr
        read(666,*) dummy, dummy, rho(i,j), ye(i,j), temp(i,j), eps(i,j), &
                    v1(i,j), v2(i,j), v3(i,j), press(i,j), lfac(i,j), &
                    sqrtgamma(i,j), psi(i,j)
                    
!   write(*,*) 'rho(i,j), ye(i,j), temp(i,j), eps(i,j), v1(i,j), &
!                v2(i,j), v3(i,j), press(i,j), lfac(i,j)'
!   write(*,*) rho(i,j), ye(i,j), temp(i,j), eps(i,j), v1(i,j), v2(i,j), v3(i,j), press(i,j), lfac(i,j)
     enddo
  enddo
  close(666)

end subroutine mod_c2b_get_prim_profile

subroutine mod_c2b_get_cons_profile(lprofile_name, Nr, Nth, rho, press,&
                                     temp, ye,&
                                     lfac, eps,&
                                     v1, v2, v3, &
                                     D, Dye, tau, &
                                     s_1, s_2, s_3, &
                                     sqrtgamma, psi)


  character*(*) lprofile_name
  character*128 :: filename

  integer, intent(in) :: Nth, Nr
  double precision, dimension(Nr,Nth), intent(out) :: rho, press, temp, ye, lfac, eps, v1, v2, v3
  double precision, dimension(Nr,Nth), intent(out) :: D, Dye, tau, S_1, S_2, S_3, sqrtgamma, psi

  integer i,j
  
  double precision :: discrim, rho_min
  double precision :: dx,dtheta,dphi, phi_max
  double precision :: dummy
  
  !--------------------------------------------------------!       
  ! read Hydroeq.dat profile      
  !--------------------------------------------------------!       
  filename = trim(adjustl(lprofile_name))//"/fluidID.csv"
  open(666,file=trim(filename),status='unknown', & 
       form='formatted',action='read')
  read(666,*) 
  read(666,*) 
  read(666,*) 

  do j=1,Nth
     do i=1,Nr
        read(666,*) dummy, dummy, rho(i,j), ye(i,j), temp(i,j), eps(i,j), &
                    v1(i,j), v2(i,j), v3(i,j), press(i,j), lfac(i,j), &
!code test
!                    sqrtgamma(i,j), dummy, D(i,j), tau(i,j), Dye(i,j), &
                    sqrtgamma(i,j), psi(i,j), D(i,j), tau(i,j), Dye(i,j), &
                    S_1(i,j), S_2(i,j), S_3(i,j)
     enddo
  enddo
  close(666)

end subroutine mod_c2b_get_cons_profile

subroutine mod_c2b_map_1D(point_value, point_radius, parray, pradius ,nr)

  implicit none
  
  double precision, intent(out)  :: point_value
  double precision, intent(in)   :: point_radius
  double precision, dimension(1:nr), intent(in) :: parray
  double precision, dimension(1:nr), intent(in) :: pradius
  integer, intent(in)   :: nr 

  integer :: ir
  integer :: pm1r
  double precision   :: dr
  double precision   :: del_r

  ! find the closest points
  ir = minloc( dabs(pradius(2:nr-1)-point_radius), dim = 1 )
  ! make sure the index falls into correct range
  ir = min( max(ir, 2), nr-1 )
  del_r = point_radius - pradius(ir)
  if (del_r > 0.0d0) then
     pm1r = 1
  else
     pm1r = -1
  end if
 
  dr = ( pradius(ir + pm1r) - pradius(ir) )
  point_value = parray(ir) &
                + del_r/dr * (parray(ir + pm1r) - parray(ir))

end subroutine mod_c2b_map_1D

subroutine mod_c2b_map_2D(point_value, point_radius, point_theta,&
                             parray, &
                             pradius ,ptheta ,nr, ntheta)
  implicit none
  
  double precision, intent(out)  :: point_value
  double precision, intent(in)   :: point_radius, point_theta
  double precision, dimension(1:nr,1:ntheta), intent(in) :: parray
  double precision, dimension(1:nr), intent(in) :: pradius
  double precision, dimension(1:ntheta), intent(in) :: ptheta
  integer, intent(in)   :: nr, ntheta

  integer :: ir, itheta
  integer :: pm1r, pm1theta
  double precision   :: dr, dtheta
  double precision   :: del_r, del_theta

  double precision   :: fh(4), a(4)

  ! find the closest points
  ir = minloc( dabs(pradius(2:nr-1)-point_radius), dim = 1 )
  itheta = minloc( dabs(ptheta(2:ntheta-1)-point_theta), dim = 1 )
  ! make sure the index falls into correct range
  ir = min( max(ir, 2), nr-1 )
  itheta = min( max(itheta, 2), ntheta-1 )

  del_r = point_radius - pradius(ir)
  del_theta = point_theta - ptheta(itheta)

  if (del_r >= 0.0d0) then
     pm1r = 1
  else
     pm1r = -1
  end if
  if (del_theta >= 0.0d0) then
     pm1theta = 1
  else
     pm1theta = -1
  end if
 
  dr = ( pradius(ir + pm1r) - pradius(ir) )
  dtheta = ( ptheta(itheta + pm1theta) - ptheta(itheta) )

  fh(1) = parray(ir, itheta)
  fh(2) = parray(ir + pm1r, itheta)
  fh(3) = parray(ir, itheta + pm1theta)
  fh(4) = parray(ir + pm1r, itheta + pm1theta)

  a(1) = fh(1)
  a(2) = (fh(2)-fh(1))/dr
  a(3) = (fh(3)-fh(1))/dtheta
  a(4) = (fh(4)-fh(2)-fh(3)+fh(1))/dtheta/dr

  point_value  = a(1) &
                 + a(2) * del_r &
                 + a(3) * del_theta &
                 + a(4) * del_theta * del_r
end subroutine mod_c2b_map_2D


end module mod_c2b

