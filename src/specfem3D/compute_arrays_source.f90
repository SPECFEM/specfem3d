!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine compute_arrays_source(ispec_selected_source,sourcearray, &
                                   xi_source,eta_source,gamma_source, &
                                   Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                   xigll,yigll,zigll,nspec)

  use constants

  implicit none

  integer :: ispec_selected_source,nspec

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision :: xi_source,eta_source,gamma_source
  double precision :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xix,xiy,xiz,etax,etay,etaz, &
        gammax,gammay,gammaz

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  ! local parameters
  double precision :: xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  double precision :: hlagrange
  double precision :: dsrc_dx, dsrc_dy, dsrc_dz
  double precision :: dxis_dx, detas_dx, dgammas_dx
  double precision :: dxis_dy, detas_dy, dgammas_dy
  double precision :: dxis_dz, detas_dz, dgammas_dz

  integer :: k,l,m

! compute Lagrange polynomials at the source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

  dxis_dx = ZERO
  dxis_dy = ZERO
  dxis_dz = ZERO
  detas_dx = ZERO
  detas_dy = ZERO
  detas_dz = ZERO
  dgammas_dx = ZERO
  dgammas_dy = ZERO
  dgammas_dz = ZERO

  do m = 1,NGLLZ
     do l = 1,NGLLY
        do k = 1,NGLLX

           xixd    = dble(xix(k,l,m,ispec_selected_source))
           xiyd    = dble(xiy(k,l,m,ispec_selected_source))
           xizd    = dble(xiz(k,l,m,ispec_selected_source))
           etaxd   = dble(etax(k,l,m,ispec_selected_source))
           etayd   = dble(etay(k,l,m,ispec_selected_source))
           etazd   = dble(etaz(k,l,m,ispec_selected_source))
           gammaxd = dble(gammax(k,l,m,ispec_selected_source))
           gammayd = dble(gammay(k,l,m,ispec_selected_source))
           gammazd = dble(gammaz(k,l,m,ispec_selected_source))

           hlagrange = hxis(k) * hetas(l) * hgammas(m)

           dxis_dx = dxis_dx + hlagrange * xixd
           dxis_dy = dxis_dy + hlagrange * xiyd
           dxis_dz = dxis_dz + hlagrange * xizd

           detas_dx = detas_dx + hlagrange * etaxd
           detas_dy = detas_dy + hlagrange * etayd
           detas_dz = detas_dz + hlagrange * etazd

           dgammas_dx = dgammas_dx + hlagrange * gammaxd
           dgammas_dy = dgammas_dy + hlagrange * gammayd
           dgammas_dz = dgammas_dz + hlagrange * gammazd

       enddo
     enddo
  enddo

! calculate source array
  sourcearrayd(:,:,:,:) = ZERO
  do m = 1,NGLLZ
     do l = 1,NGLLY
        do k = 1,NGLLX

           dsrc_dx = (hpxis(k)*dxis_dx)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dx)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dx)
           dsrc_dy = (hpxis(k)*dxis_dy)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dy)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dy)
           dsrc_dz = (hpxis(k)*dxis_dz)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dz)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dz)

           sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + (Mxx*dsrc_dx + Mxy*dsrc_dy + Mxz*dsrc_dz)
           sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + (Mxy*dsrc_dx + Myy*dsrc_dy + Myz*dsrc_dz)
           sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + (Mxz*dsrc_dx + Myz*dsrc_dy + Mzz*dsrc_dz)

       enddo
     enddo
  enddo

  ! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:), kind=CUSTOM_REAL)

  end subroutine compute_arrays_source

! =======================================================================

! compute array for acoustic source
  subroutine compute_arrays_source_acoustic(sourcearray,hxis,hetas,hgammas,factor_source)

  use constants

  implicit none

  real(kind=CUSTOM_REAL) :: factor_source
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

! local parameters
! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX) :: hxis
  double precision, dimension(NGLLY) :: hetas
  double precision, dimension(NGLLZ) :: hgammas
  integer :: i,j,k

! initializes
  sourcearray(:,:,:,:) = 0._CUSTOM_REAL
  sourcearrayd(:,:,:,:) = 0.d0

! calculates source array for interpolated location
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ! identical source array components in x,y,z-direction
        sourcearrayd(:,i,j,k) = hxis(i)*hetas(j)*hgammas(k)*dble(factor_source)
      enddo
    enddo
  enddo

! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:),kind=CUSTOM_REAL)

  end subroutine compute_arrays_source_acoustic

!================================================================

  subroutine compute_arrays_adjoint_source(adj_source_file,irec)

  use specfem_par, only: myrank,xi_receiver,eta_receiver,gamma_receiver,adj_sourcearray, &
                          xigll,yigll,zigll,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC,it

  use constants

  implicit none

! input
  integer irec
  character(len=*) adj_source_file

! Gauss-Lobatto-Legendre points of integration and weights
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)

  real(kind=CUSTOM_REAL), dimension(NTSTEP_BETWEEN_READ_ADJSRC,NDIM) :: adj_src

  integer icomp, itime, i, j, k, ier, it_start, it_end, it_sub_adj
  double precision :: junk
  ! note: should have same order as orientation in write_seismograms_to_file()
  character(len=3),dimension(NDIM) :: comp
  character(len=MAX_STRING_LEN) :: filename

  ! gets channel names
  do icomp=1,NDIM
    call write_channel_name(icomp,comp(icomp))
  enddo

  ! range of the block we need to read
  it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )
  it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC + 1
  it_end   = it_start + NTSTEP_BETWEEN_READ_ADJSRC - 1

  adj_src(:,:) = 0._CUSTOM_REAL

  ! loops over components
  do icomp = 1, NDIM

    filename = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'/../SEM/'//trim(adj_source_file)//'.'//comp(icomp)//'.adj'
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat = ier)
    ! cycles to next file (this might be more error prone)
    !if (ier /= 0) cycle
    ! requires adjoint files to exist (users will have to be more careful in setting up adjoint runs)
    if (ier /= 0) call exit_MPI(myrank, ' file '//trim(filename)//' does not exist - required for adjoint runs')

    ! reads in adjoint source trace
    !! skip unused blocks
    do itime = 1, it_start-1
      read(IIN,*,iostat=ier) junk, junk
      if (ier /= 0) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration (1111)')
    enddo
    !! read the block we need
    do itime = it_start, it_end
      read(IIN,*,iostat=ier) junk, adj_src(itime-it_start+1,icomp)
      !!! used to check whether we read the correct block
      ! if (icomp==1)      print *, junk, adj_src(itime-it_start+1,icomp)
      if (ier /= 0) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration (2222)')
    enddo
    close(IIN)

  enddo

  ! lagrange interpolators for receiver location
  call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

  ! interpolates adjoint source onto GLL points within this element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
      enddo
    enddo
  enddo

end subroutine compute_arrays_adjoint_source

!================================================================

subroutine compute_arrays_adjoint_source_SU()

  use specfem_par, only: myrank,xi_receiver,eta_receiver,gamma_receiver,adj_sourcearray,adj_sourcearrays, &
                          xigll,yigll,zigll,it,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC,nrec_local,number_receiver_global
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION
  use specfem_par_elastic, only: ELASTIC_SIMULATION
  use constants

  implicit none

! Gauss-Lobatto-Legendre points of integration and weights
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)
  real(kind=CUSTOM_REAL), dimension(NTSTEP_BETWEEN_READ_ADJSRC,NDIM) :: adj_src
  real(kind=CUSTOM_REAL), dimension(NTSTEP_BETWEEN_READ_ADJSRC) :: adj_temp
  integer i, j, k, ier, irec_local, irec, it_start, it_sub_adj

  ! note: should have same order as orientation in write_seismograms_to_file()
  character(len=MAX_STRING_LEN) :: procname, filename

  ! range of the block we need to read
  it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )
  it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC

  write(procname,"(i4)") myrank

  if (ACOUSTIC_SIMULATION) then

    filename = trim(OUTPUT_FILES)//'../SEM/'//trim(adjustl(procname))//'_dp_SU.adj'
    open(unit=IIN_SU1,file=filename,status='old',access='stream',iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
    do irec_local = 1,nrec_local
       irec = number_receiver_global(irec_local)
       read(IIN_SU1,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       adj_src(:,1)=adj_temp(:)
       adj_src(:,2)=0.0  !TRIVIAL
       adj_src(:,3)=0.0  !TRIVIAL

       ! lagrange interpolators for receiver location
       call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
       call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
       call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

       ! interpolates adjoint source onto GLL points within this element
       do k = 1, NGLLZ
         do j = 1, NGLLY
           do i = 1, NGLLX
             adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
           enddo
         enddo
       enddo

       adj_sourcearrays(irec_local,:,:,:,:,:) = adj_sourcearray(:,:,:,:,:)

    enddo
    close(IIN_SU1)

  else if (ELASTIC_SIMULATION) then

    filename = trim(OUTPUT_FILES)//'../SEM/'//trim(adjustl(procname))//'_dx_SU.adj'
    open(unit=IIN_SU1,file=filename,status='old',access='stream',iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    filename = trim(OUTPUT_FILES)//'../SEM/'//trim(adjustl(procname))//'_dy_SU.adj'
    open(unit=IIN_SU2,file=filename,status='old',access='stream',iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    filename = trim(OUTPUT_FILES)//'../SEM/'//trim(adjustl(procname))//'_dz_SU.adj'
    open(unit=IIN_SU3,file=filename,status='old',access='stream',iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    do irec_local = 1,nrec_local
       irec = number_receiver_global(irec_local)

       read(IIN_SU1,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       adj_src(:,1)=adj_temp(:)
       read(IIN_SU2,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       adj_src(:,2)=adj_temp(:)
       read(IIN_SU3,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       adj_src(:,3)=adj_temp(:)

       ! lagrange interpolators for receiver location
       call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
       call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
       call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

       ! interpolates adjoint source onto GLL points within this element
       do k = 1, NGLLZ
         do j = 1, NGLLY
           do i = 1, NGLLX
             adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
           enddo
         enddo
       enddo

       adj_sourcearrays(irec_local,:,:,:,:,:) = adj_sourcearray(:,:,:,:,:)

    enddo

    close(IIN_SU1)
    close(IIN_SU2)
    close(IIN_SU3)

  else

    call exit_MPI(myrank,'SU_FORMAT not implemented for adjoint poroelastic simulations yet')

  endif

end subroutine compute_arrays_adjoint_source_SU

