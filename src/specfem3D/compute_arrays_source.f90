!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
       Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
       hxis,hpxis,hetas,hpetas,hgammas,hpgammas,nspec)

  implicit none

  include "constants.h"

  integer ispec_selected_source
  integer nspec

  double precision Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xix,xiy,xiz,etax,etay,etaz, &
        gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,l,m

! calculate G_ij for general source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  do m=1,NGLLZ
    do l=1,NGLLY
      do k=1,NGLLX

        xixd    = dble(xix(k,l,m,ispec_selected_source))
        xiyd    = dble(xiy(k,l,m,ispec_selected_source))
        xizd    = dble(xiz(k,l,m,ispec_selected_source))
        etaxd   = dble(etax(k,l,m,ispec_selected_source))
        etayd   = dble(etay(k,l,m,ispec_selected_source))
        etazd   = dble(etaz(k,l,m,ispec_selected_source))
        gammaxd = dble(gammax(k,l,m,ispec_selected_source))
        gammayd = dble(gammay(k,l,m,ispec_selected_source))
        gammazd = dble(gammaz(k,l,m,ispec_selected_source))

        G11(k,l,m) = Mxx*xixd+Mxy*xiyd+Mxz*xizd
        G12(k,l,m) = Mxx*etaxd+Mxy*etayd+Mxz*etazd
        G13(k,l,m) = Mxx*gammaxd+Mxy*gammayd+Mxz*gammazd
        G21(k,l,m) = Mxy*xixd+Myy*xiyd+Myz*xizd
        G22(k,l,m) = Mxy*etaxd+Myy*etayd+Myz*etazd
        G23(k,l,m) = Mxy*gammaxd+Myy*gammayd+Myz*gammazd
        G31(k,l,m) = Mxz*xixd+Myz*xiyd+Mzz*xizd
        G32(k,l,m) = Mxz*etaxd+Myz*etayd+Mzz*etazd
        G33(k,l,m) = Mxz*gammaxd+Myz*gammayd+Mzz*gammazd

      enddo
    enddo
  enddo

! calculate source array for interpolated location
  do m=1,NGLLZ
    do l=1,NGLLY
      do k=1,NGLLX
        call multiply_arrays_source(sourcearrayd,G11,G12,G13,G21,G22,G23, &
                  G31,G32,G33,hxis,hpxis,hetas,hpetas,hgammas,hpgammas,k,l,m)
      enddo
    enddo
  enddo

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    sourcearray(:,:,:,:) = sngl(sourcearrayd(:,:,:,:))
  else
    sourcearray(:,:,:,:) = sourcearrayd(:,:,:,:)
  endif

  end subroutine compute_arrays_source

!=============================================================================

  subroutine compute_arrays_adjoint_source(myrank, adj_source_file, &
                                           xi_receiver,eta_receiver,gamma_receiver, adj_sourcearray, &
                                           xigll,yigll,zigll, &
                                           it_sub_adj,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC)

  implicit none

  include 'constants.h'

! input
  integer myrank, NSTEP, it_sub_adj, NTSTEP_BETWEEN_READ_ADJSRC

  double precision xi_receiver, eta_receiver, gamma_receiver

  character(len=*) adj_source_file

! output
  real(kind=CUSTOM_REAL),dimension(NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ) :: adj_sourcearray

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)

  real(kind=CUSTOM_REAL), dimension(NTSTEP_BETWEEN_READ_ADJSRC,NDIM) :: adj_src

  integer icomp, itime, i, j, k, ios, it_start, it_end
  double precision :: junk
  ! note: should have same order as orientation in write_seismograms_to_file()
  character(len=3),dimension(NDIM) :: comp != (/ "BHE", "BHN", "BHZ" /)
  character(len=256) :: filename

  ! gets channel names
  do icomp=1,NDIM
    call write_channel_name(icomp,comp(icomp))
  enddo

  ! range of the block we need to read
  it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC + 1
  it_end   = it_start + NTSTEP_BETWEEN_READ_ADJSRC - 1

  adj_src(:,:) = 0._CUSTOM_REAL

  ! loops over components
  do icomp = 1, NDIM

    filename = OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))//'/../SEM/'//trim(adj_source_file)//'.'//comp(icomp)//'.adj'
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat = ios)
    ! cycles to next file (this might be more error prone)
    !if (ios /= 0) cycle
    ! requires adjoint files to exist (users will have to be more careful in setting up adjoint runs)
    if (ios /= 0) call exit_MPI(myrank, ' file '//trim(filename)//' does not exist - required for adjoint runs')

    ! reads in adjoint source trace
    !! skip unused blocks
    do itime = 1, it_start-1
      read(IIN,*,iostat=ios) junk, junk
      if( ios /= 0 ) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration (1111)')
    enddo
    !! read the block we need
    do itime = it_start, it_end
      read(IIN,*,iostat=ios) junk, adj_src(itime-it_start+1,icomp)
      !!! used to check whether we read the correct block
      ! if (icomp==1)      print *, junk, adj_src(itime-it_start+1,icomp)
      if( ios /= 0 ) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration (2222)')
    enddo
    close(IIN)

  enddo

  ! lagrange interpolators for receiver location
  call lagrange_any(xi_receiver,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_receiver,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammar,hpgammar)

  ! interpolates adjoint source onto GLL points within this element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
      enddo
    enddo
  enddo

end subroutine compute_arrays_adjoint_source


! =======================================================================

! compute array for acoustic source
  subroutine compute_arrays_source_acoustic(sourcearray,hxis,hetas,hgammas,factor_source)

  implicit none

  include "constants.h"

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
  if(CUSTOM_REAL == SIZE_REAL) then
    sourcearray(:,:,:,:) = sngl(sourcearrayd(:,:,:,:))
  else
    sourcearray(:,:,:,:) = sourcearrayd(:,:,:,:)
  endif

  end subroutine compute_arrays_source_acoustic


! testing read in adjoint sources block by block

!!!the original version
!!!
!!!subroutine compute_arrays_adjoint_source(myrank, adj_source_file, &
!!!                    xi_receiver,eta_receiver,gamma_receiver, adj_sourcearray, &
!!!                    xigll,yigll,zigll,NSTEP)
!!!
!!!
!!!  implicit none
!!!
!!!  include 'constants.h'
!!!
!!!! input
!!!  integer myrank, NSTEP
!!!
!!!  double precision xi_receiver, eta_receiver, gamma_receiver
!!!
!!!  character(len=*) adj_source_file
!!!
!!!! output
!!!  real(kind=CUSTOM_REAL),dimension(NSTEP,NDIM,NGLLX,NGLLY,NGLLZ) :: adj_sourcearray
!!!
!!!! Gauss-Lobatto-Legendre points of integration and weights
!!!  double precision, dimension(NGLLX) :: xigll
!!!  double precision, dimension(NGLLY) :: yigll
!!!  double precision, dimension(NGLLZ) :: zigll
!!!
!!!  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
!!!        hgammar(NGLLZ), hpgammar(NGLLZ)
!!!
!!!  real(kind=CUSTOM_REAL) :: adj_src(NSTEP,NDIM)
!!!
!!!  integer icomp, itime, i, j, k, ios
!!!  double precision :: junk
!!!  ! note: should have same order as orientation in write_seismograms_to_file()
!!!  character(len=3),dimension(NDIM) :: comp = (/ "BHE", "BHN", "BHZ" /)
!!!  character(len=256) :: filename
!!!
!!!  !adj_sourcearray(:,:,:,:,:) = 0.
!!!  adj_src = 0._CUSTOM_REAL
!!!
!!!  ! loops over components
!!!  do icomp = 1, NDIM
!!!
!!!    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
!!!    open(unit=IIN,file=trim(filename),status='old',action='read',iostat = ios)
!!!    if (ios /= 0) cycle ! cycles to next file
!!!    !if (ios /= 0) call exit_MPI(myrank, ' file '//trim(filename)//'does not exist')
!!!
!!!    ! reads in adjoint source trace
!!!    do itime = 1, NSTEP
!!!
!!!      ! things become a bit tricky because of the Newmark time scheme at
!!!      ! the very beginning of the time loop. however, when we read in the backward/reconstructed
!!!      ! wavefields at the end of the first time loop, we can use the adjoint source index from 1 to NSTEP
!!!      ! (and then access it in reverse NSTEP-it+1  down to 1, for it=1,..NSTEP; see compute_add_sources*.f90).
!!!      read(IIN,*,iostat=ios) junk, adj_src(itime,icomp)
!!!      if( ios /= 0 ) &
!!!        call exit_MPI(myrank, &
!!!          'file '//trim(filename)//' has wrong length, please check with your simulation duration')
!!!    enddo
!!!    close(IIN)
!!!
!!!  enddo
!!!
!!!  ! lagrange interpolators for receiver location
!!!  call lagrange_any(xi_receiver,NGLLX,xigll,hxir,hpxir)
!!!  call lagrange_any(eta_receiver,NGLLY,yigll,hetar,hpetar)
!!!  call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammar,hpgammar)
!!!
!!!  ! interpolates adjoint source onto GLL points within this element
!!!  do k = 1, NGLLZ
!!!    do j = 1, NGLLY
!!!      do i = 1, NGLLX
!!!        adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
!!!      enddo
!!!    enddo
!!!  enddo
!!!
!!!end subroutine compute_arrays_adjoint_source

