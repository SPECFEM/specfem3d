!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


  subroutine compute_arrays_source_cmt(ispec_selected_source,sourcearray, &
                                       hxis,hetas,hgammas,hpxis,hpetas,hpgammas, &
                                       Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                                       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,nspec)

  use constants

  implicit none

  integer :: ispec_selected_source,nspec

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas
  double precision :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! local parameters
  double precision :: xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd

  double precision :: hlagrange
  double precision :: dsrc_dx, dsrc_dy, dsrc_dz
  double precision :: dxis_dx, detas_dx, dgammas_dx
  double precision :: dxis_dy, detas_dy, dgammas_dy
  double precision :: dxis_dz, detas_dz, dgammas_dz

  integer :: k,l,m

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

           dsrc_dx = (hpxis(k)*dxis_dx)*hetas(l)*hgammas(m) &
                    + hxis(k)*(hpetas(l)*detas_dx)*hgammas(m) &
                    + hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dx)

           dsrc_dy = (hpxis(k)*dxis_dy)*hetas(l)*hgammas(m) &
                    + hxis(k)*(hpetas(l)*detas_dy)*hgammas(m) &
                    + hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dy)

           dsrc_dz = (hpxis(k)*dxis_dz)*hetas(l)*hgammas(m) &
                    + hxis(k)*(hpetas(l)*detas_dz)*hgammas(m) &
                    + hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dz)

           sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + (Mxx*dsrc_dx + Mxy*dsrc_dy + Mxz*dsrc_dz)
           sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + (Mxy*dsrc_dx + Myy*dsrc_dy + Myz*dsrc_dz)
           sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + (Mxz*dsrc_dx + Myz*dsrc_dy + Mzz*dsrc_dz)

       enddo
     enddo
  enddo

  ! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:), kind=CUSTOM_REAL)

  end subroutine compute_arrays_source_cmt

!
!-------------------------------------------------------------------------------------------------
!

! compute array for acoustic source

  subroutine compute_arrays_source_forcesolution(sourcearray,hxis,hetas,hgammas,factor_source,comp_x,comp_y,comp_z,nu_source)

  use constants

  implicit none

  real(kind=CUSTOM_REAL) :: factor_source
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray
  double precision, dimension(NGLLX) :: hxis
  double precision, dimension(NGLLY) :: hetas
  double precision, dimension(NGLLZ) :: hgammas
  double precision :: comp_x,comp_y,comp_z
  double precision, dimension(NDIM,NDIM) :: nu_source

! local parameters
  integer :: i,j,k
  double precision :: hlagrange

! initializes
  sourcearray(:,:,:,:) = 0._CUSTOM_REAL

! calculates source array for interpolated location
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        hlagrange = hxis(i) * hetas(j) * hgammas(k) * dble(factor_source)
        ! identical source array components in x,y,z-direction
        sourcearray(:,i,j,k) =  hlagrange * ( nu_source(1,:) * comp_x + &
                                              nu_source(2,:) * comp_y + &
                                              nu_source(3,:) * comp_z )
      enddo
    enddo
  enddo

  end subroutine compute_arrays_source_forcesolution

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_arrays_adjoint_source(adj_source_file,irec_local)

  use specfem_par, only: myrank,source_adjoint,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC,it

  use constants

  implicit none

! input
  integer irec_local
  character(len=*) adj_source_file

! local
  integer icomp, itime, ier, it_start, it_end, it_sub_adj
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
      read(IIN,*,iostat=ier) junk, source_adjoint(icomp,irec_local,itime-it_start+1)
      !!! used to check whether we read the correct block
      ! if (icomp==1)      print *, junk, adj_src(itime-it_start+1,icomp)
      if (ier /= 0) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration (2222)')
    enddo
    close(IIN)

  enddo

  end subroutine compute_arrays_adjoint_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_arrays_adjoint_source_SU()

  use specfem_par, only: myrank,source_adjoint,it,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC,nrec_local
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION
  use specfem_par_elastic, only: ELASTIC_SIMULATION
  use constants

  implicit none
  real(kind=CUSTOM_REAL), dimension(NTSTEP_BETWEEN_READ_ADJSRC) :: adj_temp
  integer :: ier, irec_local, it_start, it_sub_adj

  ! note: should have same order as orientation in write_seismograms_to_file()
  character(len=MAX_STRING_LEN) :: procname, filename

  ! range of the block we need to read
  it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )
  it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC

  write(procname,"(i4)") myrank

  if (ACOUSTIC_SIMULATION) then

    filename = trim(OUTPUT_FILES)//'../SEM/'//trim(adjustl(procname))//'_p_SU.adj'
    open(unit=IIN_SU1,file=filename,status='old',access='stream',iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
    do irec_local = 1,nrec_local
       read(IIN_SU1,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       source_adjoint(1,irec_local,:) = adj_temp(:)
       source_adjoint(2,irec_local,:) = 0.0  !TRIVIAL
       source_adjoint(3,irec_local,:) = 0.0  !TRIVIAL
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

       read(IIN_SU1,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       source_adjoint(1,irec_local,:) = adj_temp(:)
       read(IIN_SU2,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       source_adjoint(2,irec_local,:) = adj_temp(:)
       read(IIN_SU3,pos=4*((irec_local-1)*(60+NSTEP) + 60 + it_start)+1 ) adj_temp
       source_adjoint(3,irec_local,:) = adj_temp(:)

    enddo

    close(IIN_SU1)
    close(IIN_SU2)
    close(IIN_SU3)

  else

    call exit_MPI(myrank,'SU_FORMAT not implemented for adjoint poroelastic simulations yet')

  endif

  end subroutine compute_arrays_adjoint_source_SU

!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_arrays_source_forcesolution_fluid(ispec_selected_source,sourcearray, &
                                                       hxis,hetas,hgammas,hpxis,hpetas,hpgammas, &
                                                       factor_source, comp_x,comp_y,comp_z, nu_source, &
                                                       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,nspec)
  use constants

  implicit none

  integer :: ispec_selected_source,nspec

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas
  double precision, dimension(NDIM,NDIM) :: nu_source
  double precision :: comp_x,comp_y,comp_z
  real(kind=CUSTOM_REAL) :: factor_source

  double precision :: FX, FY, FZ

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! local parameters
  double precision :: xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd

  double precision :: hlagrange
  double precision :: dsrc_dx, dsrc_dy, dsrc_dz
  double precision :: dxis_dx, detas_dx, dgammas_dx
  double precision :: dxis_dy, detas_dy, dgammas_dy
  double precision :: dxis_dz, detas_dz, dgammas_dz

  integer :: k,l,m

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

  FX = factor_source *(nu_source(1,1)*comp_x + nu_source(1,2)*comp_y +  nu_source(1,3)*comp_z)
  FY = factor_source *(nu_source(2,1)*comp_x + nu_source(2,2)*comp_y +  nu_source(2,3)*comp_z)
  FZ = factor_source *(nu_source(3,1)*comp_x + nu_source(3,2)*comp_y +  nu_source(3,3)*comp_z)

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
           !! for now fixed force direction and stf is defined after
           sourcearrayd(:,k,l,m) = sourcearrayd(:,k,l,m) + (FX*dsrc_dx + FY*dsrc_dy + FZ*dsrc_dz)

           !! to do :
           !!  this is for time changing force direction
           !! sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) +  dsrc_dx  * (stf_comp_x(t))
           !! sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) +  dsrc_dy  * (stf_comp_y(t))
           !! sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) +  dsrc_dz  * (stf_comp_z(t))
           !!
           !! after we need to add : sum(sourcearrayd(:,k,l,m)) to acoustic potential
           !! or sourcearrayd(1,k,l,m) * stf_comp_x(t) +
           !!    sourcearrayd(2,k,l,m) * stf_comp_y(t) +
           !!    sourcearrayd(3,k,l,m) * stf_comp_z(t)
           !!

       enddo
     enddo
  enddo

  ! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:), kind=CUSTOM_REAL)

  end subroutine compute_arrays_source_forcesolution_fluid
