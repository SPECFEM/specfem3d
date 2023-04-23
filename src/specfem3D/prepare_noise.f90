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


  subroutine prepare_noise()

! prepares noise simulations

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_noise
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: ier,num

  ! for noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "preparing noise arrays"
      call flush_IMAIN()
    endif

    ! checks if free surface is defined
    call sum_all_i(num_free_surface_faces,num)
    if (myrank == 0) then
      if (num == 0) then
        print *,"Error: noise simulation doesn't have a free_surface_face "
        call exit_mpi(myrank,'Error: noise simulations needs a free surface')
      endif
    endif

    ! allocates arrays
    allocate(noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1975')
    if (ier /= 0) call exit_mpi(myrank,'error allocating noise source array')

    allocate(normal_x_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1976')
    if (ier /= 0) stop 'error allocating array normal_x_noise'
    allocate(normal_y_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1977')
    if (ier /= 0) stop 'error allocating array normal_y_noise'
    allocate(normal_z_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1978')
    if (ier /= 0) stop 'error allocating array normal_z_noise'
    allocate(mask_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1979')
    if (ier /= 0) stop 'error allocating array mask_noise'
    allocate(noise_surface_movie(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1980')
    if (ier /= 0) stop 'error allocating array noise_surface_movie'

    ! initializes
    noise_sourcearray(:,:,:,:,:) = 0._CUSTOM_REAL
    normal_x_noise(:)            = 0._CUSTOM_REAL
    normal_y_noise(:)            = 0._CUSTOM_REAL
    normal_z_noise(:)            = 0._CUSTOM_REAL
    mask_noise(:)                = 0._CUSTOM_REAL
    noise_surface_movie(:,:,:) = 0._CUSTOM_REAL

    ! sets up noise source for main receiver station
    call read_parameters_noise(nrec,NSTEP,NGLLSQUARE*num_free_surface_faces, &
                               islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver,nu_rec, &
                               noise_sourcearray,xigll,yigll,zigll, &
                               ibool, &
                               xstore,ystore,zstore, &
                               irec_main_noise,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                               NSPEC_AB,NGLOB_AB, &
                               num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                               ispec_is_acoustic)

    ! read noise distribution and direction
    call read_noise_distribution_direction(NGLLSQUARE*num_free_surface_faces, &
                                           normal_x_noise,normal_y_noise,normal_z_noise,mask_noise)

    ! checks flags for noise simulation
    call check_parameters_noise(NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                LOCAL_PATH,num_free_surface_faces,NSTEP)
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_noise
