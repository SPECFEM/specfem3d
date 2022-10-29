!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

  subroutine locate_MPI_slice(npoints_subset,ipoin_already_done, &
                              ispec_selected_subset, &
                              x_found_subset, y_found_subset, z_found_subset, &
                              xi_subset,eta_subset,gamma_subset, &
                              idomain_subset,nu_subset,final_distance_subset, &
                              npoints_total, ispec_selected, islice_selected, &
                              x_found,y_found,z_found, &
                              xi_point, eta_point, gamma_point, &
                              idomain,nu_point,final_distance)

! locates subset of points in all slices

  use constants, only: HUGEVAL,NDIM
  use specfem_par, only: NPROC,myrank

  integer, intent(in) :: npoints_subset,ipoin_already_done
  integer, dimension(npoints_subset), intent(in)  :: ispec_selected_subset,idomain_subset
  double precision, dimension(npoints_subset), intent(in) :: x_found_subset, y_found_subset, z_found_subset
  double precision, dimension(npoints_subset), intent(in) :: xi_subset, eta_subset, gamma_subset
  double precision, dimension(NDIM,NDIM,npoints_subset), intent(in) :: nu_subset
  double precision, dimension(npoints_subset), intent(in) :: final_distance_subset

  integer, intent(in) :: npoints_total
  integer, dimension(npoints_total), intent(inout)  :: ispec_selected, islice_selected, idomain
  double precision, dimension(npoints_total), intent(inout)  :: x_found, y_found, z_found
  double precision, dimension(npoints_total), intent(inout)  :: xi_point, eta_point, gamma_point
  double precision, dimension(NDIM,NDIM,npoints_total), intent(inout)  :: nu_point
  double precision, dimension(npoints_total), intent(inout)  :: final_distance

  ! local parameters
  integer :: ipoin,ipoin_in_this_subset,iproc
  double precision :: distmin

  ! gather arrays
  integer, dimension(npoints_subset,0:NPROC-1) :: ispec_selected_all,idomain_all
  double precision, dimension(npoints_subset,0:NPROC-1) :: xi_all,eta_all,gamma_all
  double precision, dimension(npoints_subset,0:NPROC-1) :: x_found_all,y_found_all,z_found_all
  double precision, dimension(npoints_subset,0:NPROC-1) :: final_distance_all
  double precision, dimension(NDIM,NDIM,npoints_subset,0:NPROC-1) :: nu_all

  ! initializes with dummy values
  ispec_selected_all(:,:) = -1
  idomain_all(:,:) = -1000
  xi_all(:,:) = 0.d0
  eta_all(:,:) = 0.d0
  gamma_all(:,:) = 0.d0
  x_found_all(:,:) = 0.d0
  y_found_all(:,:) = 0.d0
  z_found_all(:,:) = 0.d0
  final_distance_all(:,:) = HUGEVAL

  ! gather all (on main process)
  call gather_all_i(ispec_selected_subset,npoints_subset,ispec_selected_all,npoints_subset,NPROC)
  call gather_all_i(idomain_subset,npoints_subset,idomain_all,npoints_subset,NPROC)

  call gather_all_dp(x_found_subset,npoints_subset,x_found_all,npoints_subset,NPROC)
  call gather_all_dp(y_found_subset,npoints_subset,y_found_all,npoints_subset,NPROC)
  call gather_all_dp(z_found_subset,npoints_subset,z_found_all,npoints_subset,NPROC)

  call gather_all_dp(xi_subset,npoints_subset,xi_all,npoints_subset,NPROC)
  call gather_all_dp(eta_subset,npoints_subset,eta_all,npoints_subset,NPROC)
  call gather_all_dp(gamma_subset,npoints_subset,gamma_all,npoints_subset,NPROC)

  call gather_all_dp(nu_subset,NDIM*NDIM*npoints_subset,nu_all,NDIM*NDIM*npoints_subset,NPROC)
  call gather_all_dp(final_distance_subset,npoints_subset,final_distance_all,npoints_subset,NPROC)

  ! find the slice and element to put the source
  if (myrank == 0) then

    ! loops over subset
    do ipoin_in_this_subset = 1,npoints_subset

      ! mapping from station/source number in current subset to real station/source number in all the subsets
      ipoin = ipoin_in_this_subset + ipoin_already_done

      distmin = HUGEVAL
      do iproc = 0,NPROC-1
        if (final_distance_all(ipoin_in_this_subset,iproc) < distmin) then
          distmin =  final_distance_all(ipoin_in_this_subset,iproc)

          islice_selected(ipoin) = iproc
          ispec_selected(ipoin) = ispec_selected_all(ipoin_in_this_subset,iproc)
          idomain(ipoin) = idomain_all(ipoin_in_this_subset,iproc)

          xi_point(ipoin)    = xi_all(ipoin_in_this_subset,iproc)
          eta_point(ipoin)   = eta_all(ipoin_in_this_subset,iproc)
          gamma_point(ipoin) = gamma_all(ipoin_in_this_subset,iproc)
          nu_point(:,:,ipoin) = nu_all(:,:,ipoin_in_this_subset,iproc)

          x_found(ipoin) = x_found_all(ipoin_in_this_subset,iproc)
          y_found(ipoin) = y_found_all(ipoin_in_this_subset,iproc)
          z_found(ipoin) = z_found_all(ipoin_in_this_subset,iproc)
        endif
      enddo
      final_distance(ipoin) = distmin
    enddo

  endif ! end of section executed by main process only

  end subroutine locate_MPI_slice


!
!------------------------------------------------------------------------------------------
!

  subroutine locate_MPI_slice_and_bcast_to_all_single(x_to_locate, y_to_locate, z_to_locate, &
                                                      x_found, y_found, z_found, &
                                                      xi, eta, gamma, ispec_selected, islice_selected, &
                                                      distance_from_target, domain, nu)

! locates MPI slice of single point and broadcasts result

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
  integer                                                        :: iproc, ier

  !! to avoid compler error when calling gather_all*
  double precision,  dimension(1)                                :: distance_from_target_dummy
  double precision,  dimension(1)                                :: xi_dummy, eta_dummy, gamma_dummy
  double precision,  dimension(1)                                :: x_found_dummy, y_found_dummy, z_found_dummy
  integer,           dimension(1)                                :: ispec_selected_dummy, islice_selected_dummy, domain_dummy

  allocate(distance_from_target_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1782')
  allocate(xi_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1783')
  allocate(eta_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1784')
  allocate(gamma_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1785')
  allocate(x_found_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1786')
  allocate(y_found_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1787')
  allocate(z_found_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1788')
  allocate(nu_all(3,3,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1789')

  allocate(ispec_selected_all(1,0:NPROC-1),domain_all(1,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1790')

  distance_from_target = dsqrt( (x_to_locate - x_found)**2 &
                               +(y_to_locate - y_found)**2 &
                               +(z_to_locate - z_found)**2)

  !! it's just to avoid compiler error
  distance_from_target_dummy(1) = distance_from_target
  xi_dummy(1) = xi
  eta_dummy(1) = eta
  gamma_dummy(1) = gamma
  ispec_selected_dummy(1) = ispec_selected
  x_found_dummy(1) = x_found
  y_found_dummy(1) = y_found
  z_found_dummy(1) = z_found
  domain_dummy(1) = domain

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
  islice_selected = islice_selected_dummy(1)
  domain = domain_dummy(1)
  ispec_selected = ispec_selected_dummy(1)
  xi = xi_dummy(1)
  eta = eta_dummy(1)
  gamma = gamma_dummy(1)
  x_found = x_found_dummy(1)
  y_found = y_found_dummy(1)
  z_found = z_found_dummy(1)
  distance_from_target = distance_from_target_dummy(1)

  deallocate(distance_from_target_all, xi_all, eta_all, gamma_all, x_found_all, y_found_all, z_found_all, nu_all)
  deallocate(ispec_selected_all,domain_all)

  end subroutine locate_MPI_slice_and_bcast_to_all_single









