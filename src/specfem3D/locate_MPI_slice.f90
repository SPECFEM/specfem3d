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


!--------------------------------------------------------------------------------------------------------------------
!  locate MPI slice which contains the point and bcast to all
!--------------------------------------------------------------------------------------------------------------------

  subroutine locate_MPI_slice_and_bcast_to_all(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
       xi, eta, gamma, ispec_selected, islice_selected, distance_from_target, domain, nu)

  use constants, only: HUGEVAL
  use specfem_par, only: NPROC,myrank

  integer,                                        intent(inout)  :: ispec_selected, islice_selected, domain
  double precision,                               intent(in)     :: x_to_locate, y_to_locate, z_to_locate
  double precision,                               intent(inout)  :: x_found,  y_found,  z_found
  double precision,                               intent(inout)  :: xi, eta, gamma, distance_from_target
  double precision, dimension(3,3),               intent(inout)  :: nu

  double precision,   dimension(:,:), allocatable                :: distance_from_target_all
  double precision,   dimension(:,:), allocatable                :: xi_all, eta_all, gamma_all
  double precision,   dimension(:,:,:), allocatable              :: nu_all
  double precision,   dimension(:,:), allocatable                :: x_found_all, y_found_all, z_found_all
  integer,            dimension(:,:), allocatable                :: ispec_selected_all, domain_all
  integer                                                        :: iproc

  !! to avoid compler error when calling gather_all*
  double precision,  dimension(1)                                :: distance_from_target_dummy
  double precision,  dimension(1)                                :: xi_dummy, eta_dummy, gamma_dummy
  double precision,  dimension(1)                                :: x_found_dummy, y_found_dummy, z_found_dummy
  integer,           dimension(1)                                :: ispec_selected_dummy, islice_selected_dummy, domain_dummy

  allocate(distance_from_target_all(1,0:NPROC-1), &
           xi_all(1,0:NPROC-1), &
           eta_all(1,0:NPROC-1), &
           gamma_all(1,0:NPROC-1), &
           x_found_all(1,0:NPROC-1), &
           y_found_all(1,0:NPROC-1), &
           z_found_all(1,0:NPROC-1), &
           nu_all(3,3,0:NPROC-1))

  allocate(ispec_selected_all(1,0:NPROC-1),domain_all(1,0:NPROC-1))

  distance_from_target = dsqrt( (x_to_locate - x_found)**2&
                               +(y_to_locate - y_found)**2&
                               +(z_to_locate - z_found)**2)

  !! it's just to avoid compiler error
  distance_from_target_dummy(1)=distance_from_target
  xi_dummy(1)=xi
  eta_dummy(1)=eta
  gamma_dummy(1)=gamma
  ispec_selected_dummy(1)=ispec_selected
  x_found_dummy(1)=x_found
  y_found_dummy(1)=y_found
  z_found_dummy(1)=z_found
  domain_dummy(1)=domain

  ! gather all on myrank=0
  call gather_all_dp(distance_from_target_dummy, 1, distance_from_target_all, 1, NPROC)
  call gather_all_dp(xi_dummy,    1,  xi_all,    1,  NPROC)
  call gather_all_dp(eta_dummy,   1,  eta_all,   1,  NPROC)
  call gather_all_dp(gamma_dummy, 1,  gamma_all, 1,  NPROC)
  call gather_all_dp(x_found_dummy, 1,  x_found_all, 1,  NPROC)
  call gather_all_dp(y_found_dummy, 1,  y_found_all, 1,  NPROC)
  call gather_all_dp(z_found_dummy, 1,  z_found_all, 1,  NPROC)
  call gather_all_dp(nu, 3*3,  nu_all, 3*3,  NPROC)
  call gather_all_i(ispec_selected_dummy, 1, ispec_selected_all, 1, NPROC)
  call gather_all_i(domain_dummy, 1, domain_all, 1, NPROC)

  ! find the slice and element to put the source
  if (myrank == 0) then

     distance_from_target = HUGEVAL

     do iproc=0, NPROC-1
       if (distance_from_target > distance_from_target_all(1,iproc)) then
         distance_from_target =  distance_from_target_all(1,iproc)
         islice_selected_dummy(1) = iproc
         ispec_selected_dummy(1) = ispec_selected_all(1,iproc)
         domain_dummy(1) = domain_all(1,iproc)
         xi_dummy(1)    = xi_all(1,iproc)
         eta_dummy(1)   = eta_all(1,iproc)
         gamma_dummy(1) = gamma_all(1,iproc)
         distance_from_target_dummy(1)=distance_from_target
         x_found_dummy(1)=x_found_all(1,iproc)
         y_found_dummy(1)=y_found_all(1,iproc)
         z_found_dummy(1)=z_found_all(1,iproc)
         nu(:,:)=nu_all(:,:,iproc)
       endif
     enddo

  endif

  ! bcast from myrank=0
  call bcast_all_i(islice_selected_dummy,1)
  call bcast_all_i(domain_dummy,1)
  call bcast_all_i(ispec_selected_dummy,1)
  call bcast_all_dp(xi_dummy,1)
  call bcast_all_dp(eta_dummy,1)
  call bcast_all_dp(gamma_dummy,1)
  call bcast_all_dp(distance_from_target_dummy,1)
  call bcast_all_dp(nu,3*3)
  call bcast_all_dp(x_found_dummy,1)
  call bcast_all_dp(y_found_dummy,1)
  call bcast_all_dp(z_found_dummy,1)

  !! it was just to avoid compler error
  islice_selected=islice_selected_dummy(1)
  domain=domain_dummy(1)
  ispec_selected=ispec_selected_dummy(1)
  xi=xi_dummy(1)
  eta=eta_dummy(1)
  gamma=gamma_dummy(1)
  x_found=x_found_dummy(1)
  y_found=y_found_dummy(1)
  z_found=z_found_dummy(1)
  distance_from_target=distance_from_target_dummy(1)

  deallocate(distance_from_target_all, xi_all, eta_all, gamma_all, x_found_all, y_found_all, z_found_all, nu_all)
  deallocate(ispec_selected_all,domain_all)

  end subroutine locate_MPI_slice_and_bcast_to_all
