
! ====================================================================== 
module mod_m1_metric_interface
! ====================================================================== 

    implicit none 
    ! ====================================================================== 
    integer, parameter :: XX_=1, XY_=2, XZ_=3, YY_=4, YZ_=5, ZZ_=6 
    ! ====================================================================== 
    ! In the no-CFC case this is unfortunately 
    ! strictly a locally allocated copy of bits and pieces 
    ! of the myM array. So long as the grid patches are not 
    ! too big this should not be an issue since it gets deallocated 
    ! right away. 
    type m1_metric_helper 
        logical :: valid 
        double precision, allocatable :: sqrtg({^D&:,}), alp({^D&:,})
        double precision, allocatable :: gammaij({^D&:,},:), gammaUPij({^D&:,},:), beta({^D&:,}, :)
        double precision, allocatable :: dalp({^D&:,}, :), dbeta({^D&:,}, :, :)
        integer, allocatable :: idimtoimu(:,:)!1:3,1:3)
        contains 
        procedure :: init 
        procedure :: destroy 
		procedure :: fill_metric
        procedure :: raise 
        procedure :: lower 
		procedure :: raise_ixD
		procedure :: lower_ixD
    end type m1_metric_helper 
    ! ====================================================================== 
    contains 
    ! ====================================================================== 
    ! Allocates metric object 
    subroutine init(self, ix^L)
        class(m1_metric_helper), intent(inout) :: self 
        integer, intent(in) :: ix^L 
        self%valid = .false. 
        if( .not. allocated(self%sqrtg) ) then
            allocate(self%sqrtg(ix^S))
            allocate(self%alp(ix^S))
            allocate(self%gammaij(ix^S, 6))
            allocate(self%gammaUPij(ix^S, 6))
            allocate(self%beta(ix^S, ^NC))
            allocate(self%dalp(ix^S, ^NC))
            allocate(self%dbeta(ix^S, ^NC, ^NC))
            allocate(self%idimtoimu(^NC, ^NC))
            self%valid = .true.
        end if 
    end subroutine init 
    ! ====================================================================== 
    ! Deallocates metric object 
    subroutine destroy(self)
        class(m1_metric_helper), intent(inout) :: self 
        if( allocated(self%sqrtg) ) then 
            deallocate(self%sqrtg)
            deallocate(self%alp)
            deallocate(self%gammaij)
            deallocate(self%gammaUPij)
            deallocate(self%beta)
            deallocate(self%dalp)
            deallocate(self%dbeta)
            deallocate(self%idimtoimu)
        end if 
        self%valid = .false.
    end subroutine destroy 
    ! ====================================================================== 
    {#IFNDEF UNIT_TESTS
    ! Fills the metricH structure used inside m1 
    !KEN added Flag get_deriv
    subroutine fill_metric(metricH,wprim,x,ixI^L,ixO^L, get_deriv)
        
        use mod_metric
        {#IFNDEF UNIT_TESTS
        include "amrvacdef.f"
        }
        ! ====================================================================== 
        ! ====================================================================== 
        class(m1_metric_helper), intent(inout) :: metricH 
        integer, intent(in) :: ixI^L,ixO^L
        double precision, intent(in)       :: x(ixI^S,1:ndim)
        double precision, intent(inout)    :: wprim(ixI^S,1:nw)
        logical, intent(in), optional :: get_deriv
        ! ====================================================================== 
        logical :: get_deriv_impl                            ! <-- local resolved variable
        double precision  :: dummy(ixI^S, 1:3,1:3) 
        double precision  :: dummy2 
        double precision  :: betaD(ixI^S,1:ndir) 
        double precision, dimension(ixI^S,1:^NC,1:^NC,1:^NC) :: dgammainvdk
        double precision, dimension(ixI^S,0:^NC,1:^NC)       :: dbetaDownjDk
        double precision  :: sqrtgHat(ixI^S)
        double precision  :: epsilon = 1.0d-30 
        integer :: idir, jdir, ncomp 
        integer :: i,j,k,inonzero,ix^D
        double precision :: psi4(ixI^S)
        double precision :: psi6(ixI^S)
        double precision :: psi4inv(ixI^S)

        ! ====================================================================== 
        if (present(get_deriv)) then
          get_deriv_impl = get_deriv
        else
          get_deriv_impl = .false.
        end if

        ! ====================================================================== 
        call metricH%init(ixI^L)

        metricH%idimtoimu(1,1) = 1 ! XX 
        metricH%idimtoimu(1,2) = 2 ! XY 
        metricH%idimtoimu(1,3) = 3 ! XZ 
    
        metricH%idimtoimu(2,1) = 2 ! XY
        metricH%idimtoimu(2,2) = 4 ! YY 
        metricH%idimtoimu(2,3) = 5 ! YZ 
    
        metricH%idimtoimu(3,1) = 3 ! XZ 
        metricH%idimtoimu(3,2) = 5 ! YZ 
        metricH%idimtoimu(3,3) = 6 ! ZZ 

        ! ======================================================================   
        ! initialize dalp with zero for different ndim
        ! KEN changed idir to :
        ! KEN this pre init is not needed for the dysp case. I copy it over to non dysp
        ! KEN moved dalp init into deriv section
        !metricH%dalp(ixI^S,:) = 0.0d0
        ! initialize dbeta with zero for different ndim
        !metricH%dbeta(ixI^S,:,:) = 0.0d0 
        ! initialize metric with flat metric
        !metricH%gammaij(ixO^S,:)   = 0.0d0 
        !metricH%gammaUPij(ixO^S,:) = 0.0d0 
        ! ====================================
        !metricH%gammaij(ixO^S,XX_)     = 1.0d0
        !metricH%gammaij(ixO^S,YY_)     = 1.0d0
        !metricH%gammaij(ixO^S,ZZ_)     = 1.0d0 
        !metricH%gammaij(ixO^S,XY_)     = 0.0d0 
        !metricH%gammaij(ixO^S,XZ_)     = 0.0d0 
        !metricH%gammaij(ixO^S,YZ_)     = 0.0d0 
        ! ====================================
        !metricH%gammaUPij(ixO^S,XX_)   = 1.0d0
        !metricH%gammaUPij(ixO^S,YY_)   = 1.0d0
        !metricH%gammaUPij(ixO^S,ZZ_)   = 1.0d0 
        !metricH%gammaUPij(ixO^S,XY_)   = 0.0d0 
        !metricH%gammaUPij(ixO^S,XZ_)   = 0.0d0 
        !metricH%gammaUPij(ixO^S,YZ_)   = 0.0d0 
        ! ======================================================================   
        {#IFDEF DY_SP 
        !KEN this gets completely overwritten in get_sqrt_gamma_hat
        !sqrtgHat(ixI^S) = 1.0d0
        ! ======================================================================
        ! CFC case 
        ! ======================================================================
        ! 1) sqrtgamma_hat 
        call get_sqrt_gamma_hat(x(ixI^S, 1:^ND), ixI^L, ixI^L, sqrtgHat(ixI^S))
        ! ======================================================================
        ! Recover physical sqrtg = sqrtg_hat * pow(psi,-4) 
        ! ======================================================================

        !KEN did this to avoid expensive **
        psi4(ixI^S)    = wprim(ixI^S, psi_metric_)**4
        psi6(ixI^S)   = psi4(ixI^S) * wprim(ixI^S, psi_metric_)**2
        psi4inv(ixI^S) = 1.0d0 / psi4

        metricH%sqrtg(ixI^S) = sqrtgHat(ixI^S) * psi6(ixI^S)

        {#IFDEF METRIC_CHECK
        {^D& do ix^DB=ixImin^DB, ixImax^DB\}
        write(5001,*)"-----ix",ix^D
        write(5001,*)"x",x(ix^D, 1) 
        write(5001,*)"abs(x)",dabs( x(ix^D, 1) )
        write(5001,*)"sqrtg",sqrtgHat(ix^D)
        write(5001,*)"psi",wprim(ix^D, psi_metric_)
        write(5001,*)"sqrtg w psi",metricH%sqrtg(ix^D)
        {^D& end do \}
        }
        ! ======================================================================
        ! 2) gamma_ij_hat 
        call get_gamma_ij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, dummy(ixI^S,1:3,1:3))

        metricH%gammaij(ixO^S,XX_) = dummy(ixO^S,1,1) * psi4(ixO^S) 
        metricH%gammaij(ixO^S,XY_) = dummy(ixO^S,1,2) * psi4(ixO^S)
        metricH%gammaij(ixO^S,XZ_) = dummy(ixO^S,1,3) * psi4(ixO^S) 
        metricH%gammaij(ixO^S,YY_) = dummy(ixO^S,2,2) * psi4(ixO^S) 
        metricH%gammaij(ixO^S,YZ_) = dummy(ixO^S,2,3) * psi4(ixO^S)
        metricH%gammaij(ixO^S,ZZ_) = dummy(ixO^S,3,3) * psi4(ixO^S)
        ! Recover physical metricH gamma_ij = gamma_ij_hat * pow(psi,-4) 
        {#IFDEF METRIC_ZERO_COMP
        if(^ND .lt. 3 ) then
            metricH%gammaij(ixO^S,XZ_) = 0.0d0
            metricH%gammaij(ixO^S,YZ_) = 0.0d0
            metricH%gammaij(ixO^S,ZZ_) = 0.0d0
            if(^ND .lt. 2) then
                metricH%gammaij(ixO^S,XY_) = 0.0d0
                metricH%gammaij(ixO^S,YY_) = 0.0d0
                metricH%gammaij(ixO^S,YZ_) = 0.0d0
                metricH%gammaij(ixO^S,ZZ_) = 0.0d0
            end if 
        end if 
        }
        ! ======================================================================
        ! 3) gamma_inv_ij_hat 
        call get_gammainvij_hat(x(ixI^S, 1:ndim), ixI^L, ixO^L, dummy(ixI^S,1:3,1:3))

        select case (coordinate)
        case (cartesian)
          
        case (cylindrical)
          ! note: r_ = 1,z_ = 2, phi_ = 3.
            {do ix^D=ixImin^D, ixImax^D \}    
            if (dabs(x(ix^D, 1)) < epsilon) then
                dummy(ix^D,3,3) = 1.0d0 / (epsilon**2)
            end if
            {end do^D&\}
        case (spherical)
          ! note: r_ = 1,theta_=2, phi_ = 3.
            {do ix^D=ixImin^D, ixImax^D \}    
            if (dabs(x(ix^D, 1)) < epsilon) then
                dummy(ix^D,2,2) = 1.0d0 / (epsilon**2)
                if((dabs(x(ix^D, 2)) < epsilon) .or. (dabs(x(ix^D, 2) - 3.14159265358979323846d0) < epsilon)) then
                    dummy(ix^D,3,3) = 1.0d0 / (epsilon**2 {^NOONED * epsilon**2 })
                end if
            end if
            {end do^D&\}
        end select


        metricH%gammaUPij(ixO^S,XX_) = dummy(ixO^S,1,1) * psi4inv(ixO^S)
        metricH%gammaUPij(ixO^S,XY_) = dummy(ixO^S,1,2) * psi4inv(ixO^S)
        metricH%gammaUPij(ixO^S,XZ_) = dummy(ixO^S,1,3) * psi4inv(ixO^S)
        metricH%gammaUPij(ixO^S,YY_) = dummy(ixO^S,2,2) * psi4inv(ixO^S)
        metricH%gammaUPij(ixO^S,YZ_) = dummy(ixO^S,2,3) * psi4inv(ixO^S)
        metricH%gammaUPij(ixO^S,ZZ_) = dummy(ixO^S,3,3) * psi4inv(ixO^S)
        ! Recover inverse physical metricH gamma^ij = gamma^ij_hat * pow(psi,4) 
        {#IFDEF METRIC_ZERO_COMP
        if(^ND .lt. 3 ) then
            metricH%gammaUPij(ixO^S,XZ_) = 0.0d0
            metricH%gammaUPij(ixO^S,YZ_) = 0.0d0
            metricH%gammaUPij(ixO^S,ZZ_) = 0.0d0
            if(^ND .lt. 2) then
                metricH%gammaUPij(ixO^S,XY_) = 0.0d0
                metricH%gammaUPij(ixO^S,YY_) = 0.0d0
                metricH%gammaUPij(ixO^S,YZ_) = 0.0d0
                metricH%gammaUPij(ixO^S,ZZ_) = 0.0d0
            end if 
        end if 
        }
        ! ======================================================================
        ! Get lapse 
        metricH%alp(ixO^S) = wprim(ixO^S, alp_metric_)
        ! ======================================================================
        ! Get shift 
        {^C& metricH%beta(ixO^S,^C) = wprim(ixO^S, beta_metric^C_) \}
		! ======================================================================
		! Get lapse partial derivative
        if (get_deriv_impl) then
        !KEN put the deriv in only when needed
        metricH%dalp(ixI^S,:) = 0.0d0
        ! initialize dbeta with zero for different ndim
        metricH%dbeta(ixI^S,:,:) = 0.0d0
		do idir = 1,ndim
          call partial_d( wprim(ixI^S, alp_metric_) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,metricH%dalp(ixI^S,idir) )
		end do 
		! ======================================================================
		! Get shift partial derivative: d_i beta^j
		do idir = 1,ndim
          {^C&
            jdir = ^C
            call partial_d( wprim(ixI^S, beta_metric^C_) ,ixI^L,ixO^L,x(ixI^S,1:ndim),idir,metricH%dbeta(ixI^S,jdir,idir) )
           \}
		end do
        endif !get_deriv_impl
        ! ======================================================================
        ! ======================================================================
        }
        {#IFNDEF DY_SP 

        !KEN moved this here, it was written even befor ifdef dysp
        metricH%dalp(ixI^S,:) = 0.0d0
        ! initialize dbeta with zero for different ndim
        metricH%dbeta(ixI^S,:,:) = 0.0d0 
        ! initialize metric with flat metric
        metricH%gammaij(ixO^S,:)   = 0.0d0 
        metricH%gammaUPij(ixO^S,:) = 0.0d0 
        ! ====================================
        metricH%gammaij(ixO^S,XX_)     = 1.0d0
        metricH%gammaij(ixO^S,YY_)     = 1.0d0
        metricH%gammaij(ixO^S,ZZ_)     = 1.0d0 
        ! ====================================
        metricH%gammaUPij(ixO^S,XX_)   = 1.0d0
        metricH%gammaUPij(ixO^S,YY_)   = 1.0d0
        metricH%gammaUPij(ixO^S,ZZ_)   = 1.0d0 
        ! ======================================================================   
        ! ======================================================================
        ! no CFC case, we are just copying stuff around 
        ! ======================================================================
        ! 1) sqrtgamma
        metricH%sqrtg(ixO^S) = myM%sqrtgamma(ixO^S)
        !call get_sqrtgamma({x(ix^DD,^D)},m%sqrtgamma(ix^D)) ! is actually same
        ! ======================================================================
        ! 2) gamma_ij and gamma^ij 

        do inonzero = 1, myM%nnonzeroDgDk
            i = myM%nonzeroDgDk(inonzero)%i
            j = myM%nonzeroDgDk(inonzero)%j
             ncomp = metricH%idimtoimu(i,j) 
            metricH%gammaij(ixO^S, ncomp) = myM%g(i,j)%elem(ixO^S) 
            metricH%gammaUPij(ixO^S, ncomp) = myM%gammainv(i,j)%elem(ixO^S)
        enddo
        ! ======================================================================
        ! 3) lapse
        metricH%alp(ixO^S) = myM%alpha(ixO^S)
        ! ======================================================================
        ! 4) shift 
        do idir=1,3
           metricH%beta(ixO^S,idir) = myM%beta(idir)%elem(ixO^S) 
        end do
	    ! ======================================================================
	    ! 5) Get lapse partial derivative
        ! PREVIOUS VERSION: does not work for non-cartesian:
        !do idir = 1,^NC
        !   !metricH%dalp(ixO^S,idir) = myM%dalphaDj(idir)%elem(ixO^S)
        !end do

        ! NOTE: unfortunately this does not work and myM%nonzerodalpha is non-existend 
        ! although dalpjha_is_zero = F :
        !do idir = 1,^NC       
        ! if (.not. dalphadj_is_zero(idir)) then     
        !    metricH%dalp(ixO^S,idir) = myM%dalphaDj(idir)%elem(ixO^S)
        !    metricH%dalp(ixO^S,idir) = myM%nonzeroDalphaDj(idir)%elem(ixO^S)
        ! end if
        !end do 

	    ! ======================================================================
	    ! 6) Get shift partial derivative
        ! PREVIOUS VERSION: does not work for non-cartesian coordinates
	    ! do idir = 1,3
        !     do jdir = 1,3
        !      !if (dbetaidj_is_zero(idir,jdir)) cycle
	    !      !metricH%dbeta(ixO^S,idir,jdir) = myM%dbetaidj(idir,jdir)%elem(ixO^S) 
        !     end do
	    ! end do 
        !==========================================================

        ! Thus we calcuated the derivatives on our own again:
        !---------------------------------------------
        if (get_deriv_impl) then
        
        if (.not. init_from_g4) then
            !5) lapse alpha: dalphadj
            do j=1, ^NC
                {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
                call get_alpha({x(ix^DD,^D)},dummy2,dalphadj=metricH%dalp(ix^D,j),jdir=j)
                !call get_alpha({x(ix^DD,^D)},metricH%alp(ix^DB),dalphadj=metricH%dalp(ix^D,j),jdir=j)
                {^D& end do \}
             end do
             !6) shift beta: dbetaidj
             do i=1,^NC
                do j=1,^NC
                {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
                call get_beta(i,{x(ix^DD,^D)},dummy2,dbetaidj= metricH%dbeta(ix^D,i,j),jdir=j)
                !call get_beta(i,{x(ix^DD,^D)},metricH%beta(ix^DB,i),dbetaidj= metricH%dbeta(ix^D,i,j),jdir=j)
                {^D& end do \}
                end do 
             end do
        else !---------------------------------------------
            ! for init_from_g4 spacetimes like cks this does not work,
            ! some arrays don't exist even though nonzeroDbetaidj = T

            call lower3(ixI^L,ixO^L,myM,metricH%beta,betaD)

            do k=1,^NC
               {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
               call get_dgammainvdk({x(ix^DD,^D)},dgammainvdk(ix^D,1:^NC,1:^NC,k),k)
               {^D& end do \}
            end do

            do j=0,^NC
               do k=1,^NC
                  {^D& do ix^DB=ixOmin^DB, ixOmax^DB\}
                  call get_g_component(0,j,{x(ix^DD,^D)},dummy2,dgdk=dbetaDownjDk(ix^D,j,k),kdir=k)
                  !call get_g_component(0,j,{x(ix^DD,^D)},metricH%beta(ix^DB,j),dgdk=dbetaDownjDk(ix^D,j,k),kdir=k)
                  {^D& end do \}
               end do
            end do     

            ! Get derivative of contravariant shift:
            ! d_j beta^i = d_j (gamma^{ik}*beta_k)
            !----------------------------------------------
            do inonzero=1,myM%nnonzeroDbetaiDj
                i = myM%nonzeroDbetaiDj(inonzero)%i
                j = myM%nonzeroDbetaiDj(inonzero)%j
                metricH%dbeta(ixO^S,i,j) = &
                     + {^C& betaD(ixO^S,^C)*dgammainvdk(ixO^S,i,^C,j)|+} &
                     + {^C& myM%gammainv(i,^C)%elem(ixO^S)*dbetaDownjDk(ixO^S,^C,j)|+}
             end do
            !! !-----------------------------------------
            ! Derivative of lapse:
            ! d_j alpha=1/2*(beta^k*d_j beta_k + beta_k*d_j beta^k - d_j g_{00})
            !---------------------------------------
            do inonzero=1, myM%nnonzeroDalphaDj
                j = myM%nonzeroDalphaDj(inonzero)%j
                metricH%dalp(ixO^S,j) = &
                     0.5d0*({^C& metricH%beta(ixO^S,^C)*dbetaDownjDk(ixO^S,^C,j) |+} &
                     + {^C& betaD(ixO^S,^C)*metricH%dbeta(ixO^S,^C,j)|+} &
                     - dbetaDownjDk(ixO^S,0,j) &
                     )
             end do
            !***********************************************
        end if 
        end if ! get_deriv_impl
        !--------------------------------------------- end if init_from_g4

        } ! #IFNDEF DY_SP
        ! ======================================================================
        ! ====================================================================== 

        !KEN changed only DB to D
        ! Ken removed his part, because he is mean
        !***********************************************
        {^D& do ix^D=ixImin^D,ixImax^D \}
           if(metricH%sqrtg(ix^D) .le. 1.0d-7 .and. metricH%sqrtg(ix^D) .ge. 0.0d0 ) then
              metricH%sqrtg(ix^D) = 1.0d-7
           else if(metricH%sqrtg(ix^D) .ge. -1.0d-7 .and. metricH%sqrtg(ix^D) .lt. 0.0d0 ) then
              metricH%sqrtg(ix^D) = -1.0d-7
           end if 
        {^D& end do \}

        
        ! ====================================================================== 
        ! ---- For testing and switching on flat metric parts in m1 -------
        {#IFDEF M1_DBETA_ZERO
          do j = 1,^NC
          {^C& metricH%dbeta(ixO^S,j,^C) = 0.0d0 \}
          end do
        }
        {#IFDEF M1_DALPHA_ZERO
          do j=1,^NC
            metricH%dalp(ixO^S,j) = 0.0d0
          end do
        }

        {#IFDEF M1_FLAT_METRIC
            ! sqrt(det(gamma))
            metricH%sqrtg(ixO^S) = 1.0d0  
            ! ==================================== 
            ! Lapse 
            metricH%alp(ixO^S)   = 1.0d0
            {^C& metricH%dalp(ixO^S,^C) = 0.0d0 \} 
            ! ====================================
            ! Beta 
            {^C& metricH%beta(ixO^S,^C) = 0.0d0 \} 
            do idir = 1, ndim   
                {^C& metricH%dbeta(ixO^S,idir,^C) = 0.0d0 \}
            end do 
            ! Gamma and Gamma up 
            ! ==================================== 
            metricH%gammaij(ixO^S,:)   = 0.0d0 
            metricH%gammaUPij(ixO^S,:) = 0.0d0 
            ! ====================================
            metricH%gammaij(ixO^S,XX_)     = 1.0d0
            metricH%gammaij(ixO^S,YY_)     = 1.0d0
            metricH%gammaij(ixO^S,ZZ_)     = 1.0d0 
            ! ====================================
            metricH%gammaUPij(ixO^S,XX_)   = 1.0d0
            metricH%gammaUPij(ixO^S,YY_)   = 1.0d0
            metricH%gammaUPij(ixO^S,ZZ_)   = 1.0d0 
        }
        ! ====================================================================== 

    end subroutine fill_metric 
    ! ====================================================================== 
    }
    {#IFDEF UNIT_TESTS2 
      ! flat metric used for unit-tests done with fruit-library
     ! ====================================================================== 
    subroutine fill_metric(metricH)
        integer, parameter :: ndim = (^ND)
        integer :: idir, jdir, ncomp 
        ! ====================================================================== 
        ! ====================================================================== 
        class(m1_metric_helper), intent(inout) :: metricH 
        ! ====================================================================== 
        ! ====================================================================== 
        integer :: ixI^L, ixO^L
        ! For testing we use one point 
        {^D& ixImin^D = 1 \}
        {^D& ixImax^D = 2 \} 
        {^D& ixOmin^D = 1 \}
        {^D& ixOmax^D = 2 \}
        ! Initialize 
        call metricH%init(ixI^L)
        ! treat indices
        metricH%idimtoimu(1,1) = 1 ! XX 
        metricH%idimtoimu(1,2) = 2 ! XY 
        metricH%idimtoimu(1,3) = 3 ! XZ 
    
        metricH%idimtoimu(2,1) = 2 ! XY
        metricH%idimtoimu(2,2) = 4 ! YY 
        metricH%idimtoimu(2,3) = 5 ! YZ 
    
        metricH%idimtoimu(3,1) = 3 ! XZ 
        metricH%idimtoimu(3,2) = 5 ! YZ 
        metricH%idimtoimu(3,3) = 6 ! ZZ 
        ! Fill 
        ! ====================================================================== 
        ! For now unit tests are flat 
        ! sqrt(det(gamma))
        metricH%sqrtg(ixO^S) = 1.0d0  
        ! ====================================================================== 
        ! Lapse 
        metricH%alp(ixO^S)   = 1.0d0
        {^C& metricH%dalp(ixO^S,^C) = 4.0d0 \} ! false for flat, but to test sources
        ! ====================================================================== 
        ! Beta 
        {^C& metricH%beta(ixO^S,^C) = 0.0d0 \} 
        do idir = 1, ndim   ! false for flat, but to test sources react to some num
            {^C& metricH%dbeta(ixO^S,idir,^C) = 7.0d0 \}
        end do 
        ! Gamma and Gamma up 
        ! ====================================================================== 
        metricH%gammaij(ixO^S,:)   = 0.0d0 
        metricH%gammaUPij(ixO^S,:) = 0.0d0 
        ! ====================================================================== 
        metricH%gammaij(ixO^S,XX_)     = 1.0d0
        metricH%gammaij(ixO^S,YY_)     = 1.0d0
        metricH%gammaij(ixO^S,ZZ_)     = 1.0d0 
        ! ====================================================================== 
        metricH%gammaUPij(ixO^S,XX_)   = 1.0d0
        metricH%gammaUPij(ixO^S,YY_)   = 1.0d0
        metricH%gammaUPij(ixO^S,ZZ_)   = 1.0d0 
        ! ====================================================================== 
    end subroutine fill_metric
    }
    ! ====================================================================== 
    subroutine raise(self,ixI^L,ixO^L,covec, vec)
        class(m1_metric_helper), intent(in) :: self
        integer, intent(in) :: ixI^L, ixO^L
        double precision, dimension(ixI^S, ^NC), intent(in) :: covec 
        double precision, dimension(ixI^S, ^NC), intent(out) :: vec
        integer :: idims,jdims, icomp 

        do idims=1,^NC
        vec(ixO^S,idims) = 0.0d0
        do jdims=1,^NC
            icomp = self%idimtoimu(idims,jdims) 
            vec(ixO^S,idims) = vec(ixO^S,idims) + self%gammaUPij(ixO^S, icomp) * covec(ixO^S,jdims)
        end do
        end do

    end subroutine raise 
    ! ====================================================================== 
    subroutine lower(self,ixI^L,ixO^L,vec, covec)
        class(m1_metric_helper), intent(in) :: self
        integer, intent(in) :: ixI^L, ixO^L
        double precision, dimension(ixI^S, ^NC), intent(in) :: vec 
        double precision, dimension(ixI^S, ^NC), intent(out) :: covec
        integer :: idims, jdims,icomp 

        do idims=1,^NC
        covec(ixO^S,idims) = 0.0d0
        do jdims=1,^NC
            icomp = self%idimtoimu(idims,jdims)
            covec(ixO^S,idims) = covec(ixO^S,idims) + self%gammaij(ixO^S, icomp) * vec(ixO^S,jdims)
        end do
        end do

    end subroutine lower 
    ! ====================================================================== 
    subroutine raise_ixD(self,ix^D,covec, vec)
        class(m1_metric_helper), intent(in) :: self
        integer, intent(in) :: ix^D
        double precision, dimension(^NC), intent(in) :: covec 
        double precision, dimension(^NC), intent(out) :: vec
        integer :: idims,jdims, icomp 
        do idims=1,^NC
        vec(idims) = 0.0d0
        do jdims=1,^NC
            icomp = self%idimtoimu(idims,jdims)
            vec(idims) = vec(idims) + self%gammaUPij(ix^D, icomp) * covec(jdims)
        end do
        end do

    end subroutine raise_ixD
    ! ====================================================================== 
    subroutine lower_ixD(self,ix^D,vec, covec)
        class(m1_metric_helper), intent(in) :: self
        integer, intent(in) :: ix^D
        double precision, dimension(^NC), intent(in) :: vec 
        double precision, dimension(^NC), intent(out) :: covec
        integer :: idims,jdims, icomp 

        do idims=1,^NC
        covec(idims) = 0.0d0
        do jdims=1,^NC
            icomp = self%idimtoimu(idims,jdims)
            covec(idims) = covec(idims) + self%gammaij(ix^D, icomp) * vec(jdims)
        end do
        end do

    end subroutine lower_ixD
 	! ====================================================================== \
! ====================================================================== 
end module mod_m1_metric_interface
! ====================================================================== 