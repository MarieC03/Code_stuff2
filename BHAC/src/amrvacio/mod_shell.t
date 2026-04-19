!================================================================================
!
!    BHAC (The Black Hole Accretion Code) solves the equations of
!    general relativistic magnetohydrodynamics and other hyperbolic systems
!    in curved spacetimes.
!
!    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
!    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
!
!    This file is part of BHAC.
!
!    BHAC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BHAC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

!=============================================================================
{#IFNDEF D1
module mod_shell
  !
  ! Contains routines for building an arbitrary shell and integrate on it
  ! Only for 2D and 3D.
  ! 03.Oct. 2017 Oliver Porth
  !
  use mod_physicaldata
  implicit none


  integer, parameter                               :: nshellshi = 512

  type tshell
     double precision                             :: r=0.0d0, {dxShell^DE=0.0d0}
     integer                                      :: {^DE& ixGmin^DE=-1,ixGmax^DE=-1}
     integer                                      :: nwmin=-1,nwmax=-1 ! allocated index-range
     double precision, dimension(:^D&,:), pointer :: w=>Null()
     double precision, dimension(:^D&,:), pointer :: x=>Null(), xShell=>Null()
     logical                                       :: allocated=.false.
  end type tshell

  type(walloc), save,dimension(ngridshi)           :: gridvars
  type(tshell), save,dimension(nshellshi)          :: shell
  
  !=============================================================================
contains
  !=============================================================================
  subroutine alloc_shell(is,ixGmax^DE,nwmax)

    include 'amrvacdef.f'

    integer, intent(in)                            :: is,ixGmax^DE,nwmax
    ! .. local ..
    integer                                        :: ixGmin^DE
    !-----------------------------------------------------------------------------

    ixGmin^DE=1;
    
    if (.not. shell(is)%allocated) then
       allocate(shell(is)%w(1,{ixGmin^DE:ixGmax^DE},1:nwmax))
       allocate(shell(is)%x(1,{ixGmin^DE:ixGmax^DE},1:ndim))
       allocate(shell(is)%xShell(1,{ixGmin^DE:ixGmax^DE},1:ndim))
       shell(is)%ixGmin^DE=ixGmin^DE;
       shell(is)%ixGmax^DE=ixGmax^DE;
       shell(is)%nwmin=1;       shell(is)%nwmax=nwmax
       shell(is)%allocated = .true.
    else
       call mpistop('Trying to allocate already allocated shell')
    end if
    
  end subroutine alloc_shell
  !=============================================================================
  subroutine dealloc_shell(is)

    include 'amrvacdef.f'
    
    integer, intent(in)                            :: is
    !-----------------------------------------------------------------------------
    
    if (shell(is)%allocated) then
       deallocate(shell(is)%w)
       deallocate(shell(is)%x)
       deallocate(shell(is)%xShell)
       shell(is)%ixGmin^DE=-1;
       shell(is)%ixGmax^DE=-1;
       shell(is)%nwmin=-1;       shell(is)%nwmax=-1
       shell(is)%allocated = .false.
       shell(is)%r=0.0d0
       shell(is)%dxShell^DE=0.0d0;
    else
       call mpistop('Trying to deallocate already deallocated shell')
    end if
    
  end subroutine dealloc_shell
  !=============================================================================
  subroutine init_gridvars(ixG^L,nwmax)
    
    include 'amrvacdef.f'

    integer, intent(in)                             :: ixG^L, nwmax
    ! .. local ..
    integer                                         :: igrid, iigrid
    double precision,dimension(0:nwmax)             :: normconv 
    !-----------------------------------------------------------------------------

    if (nwmax .lt. nw) call mpistop('init_gridvars: nwmax needs to be .ge. than nw')
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);         
       call set_tmpGlobals(igrid)
       
       call alloc_pw(gridvars(igrid),ixG^L,1,nwmax)
       gridvars(igrid)%w(ixG^S,1:nw) = pw(igrid)%w(ixG^T,1:nw)

       if(nwmax>nw) then 
              ! M1_test
             !call specialvar_output(ixG^LL,ixM^LL^LADD1,nwmax,&
             !gridvars(igrid)%w,ps(igrid),normconv,pgeo(igrid)%dx, node(plevel_,igrid), psold(igrid))
             call specialvar_output(ixG^LL,ixM^LL^LADD1,nwmax,&
             gridvars(igrid)%w,ps(igrid),normconv,pgeo(igrid)%dx, node(plevel_,igrid), psold(igrid), psm1(igrid))
       endif
       !if(nwmax>nw) call specialvar_output(ixG^LL,ixM^LL^LADD1,nwmax,&
       !     gridvars(igrid)%w,ps(igrid),normconv)

       !if(saveprim)  call primitive(ixG^LL,ixG^LL,&
       !     gridvars(igrid)%w(ixG^T,1:nw),px(igrid)%x)

    end do

  end subroutine init_gridvars
  !=============================================================================
  subroutine finish_gridvars()

    include 'amrvacdef.f'

    ! .. local ..
    integer             :: iigrid, igrid
    !-----------------------------------------------------------------------------
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call dealloc_pw(gridvars(igrid))
    end do

  end subroutine finish_gridvars
  !=============================================================================
  subroutine fill_shell(is,rshell,SphToCoord)

    use mod_interpolate, only: find_point_ipe, interpolate_var
    include 'amrvacdef.f'
    
    integer, intent(in)                                        :: is
    double precision, intent(in)                               :: rshell

    interface
       subroutine SphToCoord(ixI^L,ixO^L,xKS,xCKS)
         integer,intent(in)                                     :: ixI^L, ixO^L
         double precision, dimension(ixI^S,1:^ND), intent(in)   :: xKS
         double precision, dimension(ixI^S,1:^ND), intent(out)  :: xCKS
       end subroutine SphToCoord
    end interface
    
    ! .. local ..
    integer                                       :: ix^DE, N^Z, N^PHI, igrid, ipe, iw
    !-----------------------------------------------------------------------------
    associate (s=>shell(is))

      if (.not. s%allocated) call mpistop('fill_shell: trying to fill unallocated shell')

      s%r = rshell
      s%w = 0.0d0
      
      ! First get the coordinates:
      {^IFZIN
      N^Z = s%ixGmax^Z-s%ixGmin^Z+1
      s%dxshell^Z   = dpi/dble(N^Z-1)
      }{^IFPHIIN
      N^PHI = s%ixGmax^PHI-s%ixGmin^PHI+1
      s%dxshell^PHI = 2.0d0*dpi/dble(N^PHI-1)
      }      

      {do ix^DE=  s%ixG^LIM^DE \}
         {^IFZIN
         s%xShell(1,ix^DE,^Z)=(ix^Z-1)*s%dxShell^Z
         }{^IFPHIIN
         s%xShell(1,ix^DE,^PHI)=(ix^PHI-1)*s%dxShell^PHI
         }
         s%xShell(1,ix^DE,1)=rshell
      {end do^DE& \}

      call SphToCoord(1,s%ixGmin^DE,1,s%ixGmax^DE,1,s%ixGmin^DE,1,s%ixGmax^DE,s%xShell,s%x)

      
      {do ix^DE=  s%ixG^LIM^DE \}

         !Finding processor and blocks for each point
         call find_point_ipe(s%x(1,ix^DE,1:ndim),igrid,ipe)
         
         if(ipe.eq.mype)then
            
            do iw=s%nwmin,s%nwmax
               call interpolate_var(igrid,ixG^LL,gridvars(igrid)%w(ixG^T,iw),&
                    px(igrid)%x,s%x(1,ix^DE,1:ndim),s%w(1,ix^DE,iw))
            end do
                   
         end if
      {end do^DE&\}
      
      
      ! Reduce to head-node:
      if (mype==0) then
         call MPI_REDUCE(MPI_IN_PLACE,s%w,size(s%w),MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
      else
         call MPI_REDUCE(s%w,s%w,size(s%w),MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
      end if
      call MPI_BARRIER(icomm, ierrmpi)
      
    end associate
  end subroutine fill_shell
  !=============================================================================
  subroutine write_shell

    include 'amrvacdef.f'

    ! Writes a topological sphere
    ! 03.Oct 2017
    integer                                   :: is
    logical, save                             :: firstshell=.true.
    !-----------------------------------------------------------------------------

    if (firstshell) then
       ishell=shellNext
       firstshell=.false.
    end if

    call init_gridvars(ixG^LL,nw+nwauxio)
    do is=1,nshells
       call put_shell(shellcoord(is))
    end do
    call finish_gridvars

    ishell=ishell+1

  end subroutine write_shell
  !=============================================================================
  subroutine put_shell(rshell)

    use mod_metric, only: SPToCoord
    include 'amrvacdef.f'

    double precision, intent(in)               :: rshell
    !-----------------------------------------------------------------------------

    call alloc_shell(1, nxShell^DE, nw+nwauxio)

    call fill_shell(1, rshell, SPToCoord)

    select case(shell_type)
    case('csv')
       call put_shell_csv(shell(1),ishell)
    case('vtu')
       call put_shell_vtu(shell(1),ishell)       
    case default
       call mpistop('put_shell: unknown shell_type')
    end select
    
    call dealloc_shell(1)

  end subroutine put_shell
  !=============================================================================
  subroutine put_shell_csv(shell,isout)

    use mod_tocart
    include 'amrvacdef.f'

    type(tshell), intent(in)      :: shell  ! the shell structure
    integer, intent(in)           :: isout  ! the output index
    ! .. local ..
    character(len=1024)           :: filename, xlabel
    logical                       :: fileopen
    character(len=10)             :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
    character(len=1024)           :: outfilehead
    integer                       :: iw, idir, ix^DE
    character(len=1024)           :: line, data
    double precision, dimension(1,{shell%ixGmin^DE:shell%ixGmax^DE},shell%nwmin:shell%nwmax) :: wCart
    double precision, dimension(1,{shell%ixGmin^DE:shell%ixGmax^DE},1:ndim)          :: xCart
    double precision, parameter :: minvalue = 1.0d-99, maxvalue = 1.0d+99
    double precision            :: roundoff_minmax
    !-----------------------------------------------------------------------------

    if (mype==0) then

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Open the file:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       inquire(unitshell,opened=fileopen)
       if(.not.fileopen)then
          ! generate filename: 
          write(xlabel,"(D10.3)") shell%r
          if(shell%r>=zero)then
             write(xlabel(1:1),"(a)") "+"
          end if
          write(filename,"(a,i4.4,a)") TRIM(filenameout)//'_r'//trim(xlabel)//'_n',isout,'.csv'
          open(unitshell,file=filename,status='unknown',form='formatted')
       end if

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Write the header:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call getheadernames(wnamei,xandwnamei,outfilehead)
       line=''
       do iw=1,ndim+nw+nwauxio-1
          line = trim(line)//trim(xandwnamei(iw))//', '
       end do
       line = trim(line)//trim(xandwnamei(ndim+nw+nwauxio))

       write(unitshell,'(a)') trim(line)

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate Cartesian components:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call cartesian_covariant(1,shell%ixGmin^DE,1,shell%ixGmax^DE,1,shell%ixGmin^DE,1,shell%ixGmax^DE,shell%w,shell%x,&
            wCart,xCart,wCoord_is_primitive=saveprim)

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Write data to file:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       {do ix^DE=  shell%ixG^LIM^DE \}

          line=''
          do idir=1,ndim
             write(data,"(es14.6,a)") roundoff_minmax(xCart(1,ix^DE,idir),minvalue,maxvalue),', '
             line=trim(line)//trim(data)
          end do
          do iw=1,nw+nwauxio-1
             write(data,"(es14.6,a)") roundoff_minmax(wCart(1,ix^DE,iw),minvalue,maxvalue),', '
             line=trim(line)//trim(data)
          end do
          write(data,"(es14.6)") roundoff_minmax(wCart(1,ix^DE,nw+nwauxio),minvalue,maxvalue)
          line=trim(line)//trim(data)

          write(unitshell,"(a)") trim(line)
       
       {end do^DE&\}

       close(unitshell)
       
    end if! mype.eq.0
       
  end subroutine put_shell_csv
  !=============================================================================
  subroutine put_shell_vtu(shell,isout)

    use mod_tocart
    include 'amrvacdef.f'

    type(tshell), intent(in)      :: shell  ! the shell structure
    integer, intent(in)           :: isout  ! the output index
    ! .. local ..
    character(len=1024)           :: filename, xlabel
    logical                       :: fileopen
    character(len=10)             :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
    character(len=1024)           :: outfilehead
    integer                       :: iw, ix^DE
    character(len=1024)           :: line, data
    double precision, dimension(1,{shell%ixGmin^DE:shell%ixGmax^DE},shell%nwmin:shell%nwmax) :: wCart
    double precision, dimension(1,{shell%ixGmin^DE:shell%ixGmax^DE},1:ndim)          :: xCart
    double precision, parameter :: minvalue = 1.0d-99, maxvalue = 1.0d+99
    double precision            :: roundoff_minmax
    integer                     :: nx^DE, nxC^DE, nc, np, icell
    double precision            :: x_VTK(1:3)
    integer                     :: VTK_type
    !-----------------------------------------------------------------------------

    if (mype==0) then

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Open the file:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       inquire(unitshell,opened=fileopen)
       if(.not.fileopen)then
          ! generate filename: 
          write(xlabel,"(D10.3)") shell%r
          if(shell%r>=zero)then
             write(xlabel(1:1),"(a)") "+"
          end if
          write(filename,"(a,i4.4,a)") TRIM(filenameout)//'_r'//trim(xlabel)//'_n',isout,'.vtu'
          open(unitshell,file=filename,status='unknown',form='formatted')
       end if

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Get the header:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call getheadernames(wnamei,xandwnamei,outfilehead)


       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate Cartesian components:
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call cartesian_covariant(1,shell%ixGmin^DE,1,shell%ixGmax^DE,1,shell%ixGmin^DE,1,shell%ixGmax^DE,shell%w,shell%x,&
            wCart,xCart,wCoord_is_primitive=saveprim)


       
       nx^DE=shell%ixGmax^DE-shell%ixGmin^DE; ! cells
       nxC^DE=nx^DE+1;                        ! corners
       nc={nx^DE*}
       np={nxC^DE*}

       
       ! generate xml header
       write(unitshell,'(a)')'<?xml version="1.0"?>'
       write(unitshell,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
       {#IFDEF BIGENDIAN write(unitshell,'(a)')' version="0.1" byte_order="BigEndian">'}
       {#IFNDEF BIGENDIAN write(unitshell,'(a)')' version="0.1" byte_order="LittleEndian">'}
       write(unitshell,'(a)')'  <UnstructuredGrid>'
       write(unitshell,'(a)')'<FieldData>'
       write(unitshell,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
            'NumberOfTuples="1" format="ascii">'
       write(unitshell,*) real(t)
       write(unitshell,'(a)')'</DataArray>'
       write(unitshell,'(a)')'</FieldData>'
       
       
       ! write out one VTK PIECE
       write(unitshell,'(a,i7,a,i7,a)') &
            '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'

       
       !==============================
       ! output Pointdata
       !==============================
       write(unitshell,'(a)')'<PointData>'
       do iw=1,nw+nwauxio
          if(iw<=nw) then 
             if(.not.writew(iw)) cycle
          endif
          write(unitshell,'(a,a,a)')&
               '<DataArray type="Float64" Name="',TRIM(wnamei(iw)),'" format="ascii">'
          write(unitshell,'(200(1pe14.6))') {^DE&(|}roundoff_minmax(wCart(1,ix^DE,iw),minvalue,maxvalue),{ix^DE=shell%ixGmin^DE,shell%ixGmax^DE)}
          write(unitshell,'(a)')'</DataArray>'
       end do
       write(unitshell,'(a)')'</PointData>'
       !==============================
       ! Done: Output Pointdata
       !==============================

       !==============================
       ! output Cornerpoints
       !==============================
       write(unitshell,'(a)')'<Points>'
       write(unitshell,'(a)')'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
       ! write cell corner coordinates in a backward dimensional loop, always 3D output
      {^DE& do ix^DEB=shell%ixGmin^DEB,shell%ixGmax^DEB \}
            x_VTK(1:3)=zero;
            x_VTK(1:ndim)=xCart(1,ix^DE,1:ndim)
            write(unitshell,'(3(1pe14.6))') x_VTK
      {^DE&end do \}
      write(unitshell,'(a)')'</DataArray>'
      write(unitshell,'(a)')'</Points>'
      !==============================
      ! Done: output Cornerpoints
      !==============================


      !==============================
      ! cell Metainformation
      !==============================
      write(unitshell,'(a)')'<Cells>'

      ! connectivity part
      write(unitshell,'(a)')'<DataArray type="Int32" Name="connectivity" format="ascii">'

      {^DE& do ix^DEB=1,nx^DEB\}
      {^IFTWOD
      write(unitshell,'(2(i7,1x))')ix2-1,ix2
      }{^IFTHREED
      write(unitshell,'(4(i7,1x))')(ix3-1)*nxC2+ix2-1, &
           (ix3-1)*nxC2+ix2,ix3*nxC2+ix2-1,ix3*nxC2+ix2
      }{^DE& end do\}
              
      write(unitshell,'(a)')'</DataArray>'

      ! offsets data array
      write(unitshell,'(a)')'<DataArray type="Int32" Name="offsets" format="ascii">'
      do icell=1,nc
         write(unitshell,'(i7)') icell*(2**(^ND-1))
      end do
      write(unitshell,'(a)')'</DataArray>'

      ! VTK cell type data array
      write(unitshell,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
      ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
      {^IFTWOD VTK_type=3 \}
      {^IFTHREED VTK_type=8 \}
      do icell=1,nc
         write(unitshell,'(i2)') VTK_type
      enddo
      write(unitshell,'(a)')'</DataArray>'
      
      write(unitshell,'(a)')'</Cells>'
      !==============================
      ! Done: cell Metainformation
      !==============================
      write(unitshell,'(a)')'</Piece>'

      
       write(unitshell,'(a)')'</UnstructuredGrid>'
       write(unitshell,'(a)')'</VTKFile>'
       close(unitshell)

    end if

  end subroutine put_shell_vtu
  !=============================================================================
  subroutine shell_integrate(s,iw,gfac,int,imask)

    include 'amrvacdef.f'
    
    type(tshell), intent(in)      :: s      ! the shell structure
    integer, intent(in)           :: iw     ! the variable number to integrate
    integer, optional, intent(in) :: imask  ! the optional mask
    double precision, intent(out) :: int    ! the integral over the shell

    interface
       subroutine gfac(ixI^L,ixO^L,xShell,gg) ! the geometric factor
         integer,intent(in)                                     :: ixI^L, ixO^L
         double precision, dimension(ixI^S,1:^ND), intent(in)   :: xShell
         double precision, dimension(ixI^S), intent(out)        :: gg
       end subroutine gfac
    end interface

    ! .. local ..
    double precision, dimension(1,{s%ixGmin^DE:s%ixGmax^DE})   :: gg, mask
    integer                                                    :: ix^D, ip^DE, itmp^DE, icnt
    double precision                                           :: tmp
    !-----------------------------------------------------------------------------

    if (mype.ne.0) then

       call mpistop("shell_integrate: called for mype.ne.0, currently only headnode has the full shell!")

    else

       if (present(imask)) then
          mask(:^D&) = s%w(:^D&,imask)
       else
          mask(:^D&) = 1.0d0
       end if
       
       call gfac(1,s%ixGmin^DE,1,s%ixGmax^DE,1,s%ixGmin^DE,1,s%ixGmax^DE,s%xshell,gg)


       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       int = 0.0d0
       {do ix^DE= s%ixGmin^DE,s%ixGmax^DE-1\}

       ! Average to center (not very fast, but does it matter?):
       icnt=0
       tmp=0.0d0
       {do ip^DE= 0,1\}
       itmp^DE=ix^DE+ip^DE;
       icnt=icnt+1
       tmp = tmp + s%w(1,itmp^DE,iw) * mask(1,itmp^DE) * gg(1,itmp^DE)
       {end do^DE&\}
       tmp = tmp/dble(icnt)

       int = int + tmp * {s%dxShell^DE|*}

       {end do^DE&\}
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    end if

  end subroutine shell_integrate
  !=============================================================================
end module mod_shell
}
{^IFONED
! Shell only makes sense in 2D or 3D.
}
!=============================================================================
