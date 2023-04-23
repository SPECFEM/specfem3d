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

!==============================================================================
!> \file save_forward_arrays_adios.F90
!!
!! \author MPBL
!==============================================================================

#include "config.fh"


  subroutine save_forward_arrays_adios()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  use pml_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  !--- Local parameters for ADIOS ---
  character(len=MAX_STRING_LEN) :: output_name,group_name
  integer(kind=8) :: group_size_inc
  integer(kind=8) :: local_dim

  !--- Variables to allreduce - wmax stands for world_max
  integer :: NGLOB_wmax, NSPEC_ATTENUATION_wmax, NSPEC_STRAIN_wmax, N_SLS_wmax
  integer, parameter :: num_vars = 4
  integer, dimension(num_vars) :: max_global_values

  !-----------------------------------------------------------------.
  ! Get maximum value for each variable used to define a local_dim. |
  ! ADIOS write equally sized chunks for each processor.            |
  !-----------------------------------------------------------------'
  ! Filling a temporary array to avoid doing allreduces for each var.
  max_global_values(1) = NGLOB_AB
  max_global_values(2) = NSPEC_ATTENUATION_AB
  max_global_values(3) = NSPEC_STRAIN_ONLY
  max_global_values(4) = N_SLS

  call max_allreduce_i(max_global_values,num_vars)

  NGLOB_wmax                   = max_global_values(1)
  NSPEC_ATTENUATION_wmax       = max_global_values(2)
  NSPEC_STRAIN_wmax            = max_global_values(3)
  N_SLS_wmax                   = max_global_values(4)

  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  group_size_inc = 0
  output_name = get_adios_filename(trim(LOCAL_PATH) // "/forward_arrays")

  ! initializes i/o group
  group_name = "SPECFEM3D_FORWARD_ARRAYS"
  call init_adios_group(myadios_fwd_group,group_name)

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(myadios_fwd_group, group_size_inc, '', STRINGIFY_VAR(ngllx))
  call define_adios_scalar(myadios_fwd_group, group_size_inc, '', STRINGIFY_VAR(nglly))
  call define_adios_scalar(myadios_fwd_group, group_size_inc, '', STRINGIFY_VAR(ngllz))

  call define_adios_scalar(myadios_fwd_group, group_size_inc, '', "nglob", NGLOB_AB)
  call define_adios_scalar(myadios_fwd_group, group_size_inc, '', STRINGIFY_VAR(NSPEC_ATTENUATION_AB))
  call define_adios_scalar(myadios_fwd_group, group_size_inc, '', STRINGIFY_VAR(NSPEC_STRAIN_ONLY))
  call define_adios_scalar(myadios_fwd_group, group_size_inc, "", STRINGIFY_VAR(N_SLS))

  ! acoustic wavefields
  if (ACOUSTIC_SIMULATION) then
    local_dim = NGLOB_wmax
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(potential_acoustic))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(potential_dot_acoustic))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(potential_dot_dot_acoustic))
  endif

  ! elastic wavefields
  if (ELASTIC_SIMULATION) then
    local_dim = NDIM * NGLOB_wmax
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(displ))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(veloc))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(accel))
    if (ATTENUATION) then
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_wmax * N_SLS_wmax
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(R_xx))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(R_yy))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(R_xy))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(R_xz))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(R_yz))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, "", &
                                       STRINGIFY_VAR(R_trace))

      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_wmax
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_xx))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_yy))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_xy))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_xz))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_yz))
      call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, "", &
                                       STRINGIFY_VAR(epsilondev_trace))
    endif
  endif

  ! poroelastic wavefields
  if (POROELASTIC_SIMULATION) then
    local_dim = NDIM * NGLOB_wmax
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(displs_poroelastic))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(velocs_poroelastic))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(accels_poroelastic))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(displw_poroelastic))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(velocw_poroelastic))
    call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                     STRINGIFY_VAR(accelw_poroelastic))
  endif

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  ! opens file for writing
  call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,output_name,group_name)

  ! sets group size
  call set_adios_group_size(myadios_fwd_file,group_size_inc)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call write_adios_scalar(myadios_fwd_file,myadios_fwd_group, STRINGIFY_VAR(ngllx))
  call write_adios_scalar(myadios_fwd_file,myadios_fwd_group, STRINGIFY_VAR(nglly))
  call write_adios_scalar(myadios_fwd_file,myadios_fwd_group, STRINGIFY_VAR(ngllz))

  call write_adios_scalar(myadios_fwd_file,myadios_fwd_group, "nglob",NGLOB_AB)

  call write_adios_scalar(myadios_fwd_file,myadios_fwd_group, STRINGIFY_VAR(NSPEC_ATTENUATION_AB))
  call write_adios_scalar(myadios_fwd_file,myadios_fwd_group, STRINGIFY_VAR(NSPEC_STRAIN_ONLY))
  call write_adios_scalar(myadios_fwd_file,myadios_fwd_group, STRINGIFY_VAR(N_SLS))

  if (ACOUSTIC_SIMULATION) then
    local_dim = NGLOB_wmax
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_acoustic))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_dot_acoustic))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_dot_dot_acoustic))
  endif
  if (ELASTIC_SIMULATION) then
    local_dim = NDIM * NGLOB_wmax
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displ))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(veloc))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accel))
    if (ATTENUATION) then
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_wmax * N_SLS_wmax
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xx))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_yy))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xy))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xz))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_yz))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_trace))

      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_wmax
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_xx))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_yy))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_xy))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_xz))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_yz))
      call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_trace))
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    local_dim = NDIM * NGLOB_wmax
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displs_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(velocs_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accels_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displw_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(velocw_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file,myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accelw_poroelastic))
  endif

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call write_adios_perform(myadios_fwd_file)

  ! closes file
  call close_file_adios(myadios_fwd_file)

  end subroutine save_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simulations when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios() except than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.

  subroutine save_forward_arrays_undoatt_adios()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Local parameters
  integer :: iteration_on_subset_tmp
  character(len=MAX_STRING_LEN) :: file_name
  integer(kind=8),save :: group_size_inc
  integer(kind=8) :: local_dim
  ! ADIOS variables
  character(len=MAX_STRING_LEN) :: group_name
  ! multiple/single file for storage of snapshots
  logical :: do_open_file,do_close_file,do_init_group

  !--- Variables to allreduce - wmax stands for world_max
  integer,save :: NGLOB_wmax, NSPEC_ATTENUATION_wmax, N_SLS_wmax

  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  !-----------------------------------------------------------------.
  ! Get maximum value for each variable used to define a local_dim. |
  ! ADIOS write equally sized chunks for each processor.            |
  !-----------------------------------------------------------------'
  ! only need to do this once for the first iteration and save the max values
  if (iteration_on_subset_tmp == 1) then
    ! Filling a temporary array to avoid doing allreduces for each var.
    call max_allreduce_singlei(NGLOB_AB,NGLOB_wmax)
    call max_allreduce_singlei(NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_wmax)
    call max_allreduce_singlei(N_SLS,N_SLS_wmax)        ! probably not needed, N_SLS will be the same for all partitions
  endif

  ! file handling
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) then
    ! single file for all steps
    do_open_file = .false.
    do_close_file = .false.
    do_init_group = .false.

    ! single file
    file_name = get_adios_filename(trim(LOCAL_PATH) // "/forward_arrays_undoatt",ADIOS2_ENGINE_UNDO_ATT)

    group_name = "SPECFEM3D_FORWARD_ARRAYS_UNDOATT"

    ! open file at first call of this routine
    if (is_adios_version1) then
      ! adds steps by appending to file
      do_open_file = .true.
      do_close_file = .true.
      if (iteration_on_subset_tmp == 1) do_init_group = .true. ! only needs to initialize group once
    else
      ! adds steps by commands
      if (.not. is_initialized_fwd_group) then
        do_open_file = .true.
        do_init_group = .true.
      endif
    endif
  else
    ! for each step a single file
    do_open_file = .true.
    do_close_file = .true.
    do_init_group = .true.

    ! files for each iteration step
    write(file_name,'(a, a, i6.6)') trim(LOCAL_PATH), '/save_frame_at', iteration_on_subset_tmp
    file_name = get_adios_filename(trim(file_name))

    write(group_name, '(a, i6)') "SPECFEM3D_FORWARD_ARRAYS_UNDOATT", iteration_on_subset_tmp
  endif

  ! debug
  !if (myrank == 0) print *,'debug: undoatt adios: save forward step iteration_on_subset_tmp = ',iteration_on_subset_tmp, &
  !                         'open/close file',do_open_file,do_close_file,do_init_group

  ! opens file for writing
  if (do_open_file) then
    ! prepares group & metadata
    !
    ! note, see adios manual:
    ! "These routines prepare ADIOS metadata construction,
    ! for example, setting up groups, variables, attributes and IO transport method,
    ! and hence must be called before any other ADIOS I/O operations,
    ! i.e., adios_open, adios_group_size, adios_write, adios_close."
    if (do_init_group) then
      call init_adios_group_undo_att(myadios_fwd_group,group_name)

      ! adds wavefield compression
      if (ADIOS_COMPRESSION_ALGORITHM /= 0) then
        ! sets adios flag to add compression operation for the following define_adios_** function calls
        call define_adios_compression()
      endif

      ! defines ADIOS variables
      group_size_inc = 0
      ! iteration number
      call define_adios_scalar(myadios_fwd_group, group_size_inc, '', "iteration", iteration_on_subset_tmp)

      ! wavefields (displ/veloc/accel) for all domains
      ! acoustic wavefields
      if (ACOUSTIC_SIMULATION) then
        local_dim = NGLOB_wmax
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(potential_acoustic))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(potential_dot_acoustic))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(potential_dot_dot_acoustic))
      endif

      ! elastic wavefields
      if (ELASTIC_SIMULATION) then
        local_dim = NDIM * NGLOB_wmax
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(displ))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(veloc))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(accel))
        if (ATTENUATION) then
          local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_wmax * N_SLS_wmax
          call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                           STRINGIFY_VAR(R_xx))
          call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                           STRINGIFY_VAR(R_yy))
          call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                           STRINGIFY_VAR(R_xy))
          call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                           STRINGIFY_VAR(R_xz))
          call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                           STRINGIFY_VAR(R_yz))
          call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, "", &
                                           STRINGIFY_VAR(R_trace))

          ! strain not needed anymore, to save file diskspace - will be re-constructed based on b_displ...
          !local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_wmax
          !call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
          !                                 STRINGIFY_VAR(epsilondev_xx))
          !call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
          !                                 STRINGIFY_VAR(epsilondev_yy))
          !call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
          !                                 STRINGIFY_VAR(epsilondev_xy))
          !call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
          !                                 STRINGIFY_VAR(epsilondev_xz))
          !call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
          !                                 STRINGIFY_VAR(epsilondev_yz))
          !call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, "", &
          !                                 STRINGIFY_VAR(epsilondev_trace))
        endif
      endif

      ! poroelastic wavefields
      if (POROELASTIC_SIMULATION) then
        local_dim = NDIM * NGLOB_wmax
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(displs_poroelastic))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(velocs_poroelastic))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(accels_poroelastic))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(displw_poroelastic))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(velocw_poroelastic))
        call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(accelw_poroelastic))
      endif

      ! re-sets compression flag (in case other routines will call the define_adios_** function calls)
      if (ADIOS_COMPRESSION_ALGORITHM /= 0) use_adios_compression = .false.
    endif

    ! Open an ADIOS handler to the restart file.
    if (is_adios_version1) then
      ! checks if we open for first time or append
      if (iteration_on_subset_tmp == 1) then
        ! creates new file
        call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,file_name,group_name)
      else
        ! append to existing file
        call open_file_adios_write_append(myadios_fwd_file,myadios_fwd_group,file_name,group_name)

        ! debug: note, do not call as the inquiry on the appended file handle will seg-fault
        ! call show_adios_file_variables(myadios_fwd_file,myadios_fwd_group,file_name)
      endif

      ! debug
      !if (myrank == 0) print *,'debug: undoatt adios: save forward step = ',iteration_on_subset_tmp,' handle ',myadios_fwd_file

    else
      ! version 2, only opens once at beginning
      call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,file_name,group_name)
    endif

    call set_adios_group_size(myadios_fwd_file,group_size_inc)

    ! sets flag
    is_initialized_fwd_group = .true.
  endif

  ! indicate new step section
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) call write_adios_begin_step(myadios_fwd_file)

  ! iteration number
  call write_adios_scalar(myadios_fwd_file, myadios_fwd_group, "iteration", iteration_on_subset_tmp)

  ! write the previously defined variable to the ADIOS file
  if (ACOUSTIC_SIMULATION) then
    local_dim = NGLOB_wmax
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_acoustic))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_dot_acoustic))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_dot_dot_acoustic))
  endif
  if (ELASTIC_SIMULATION) then
    local_dim = NDIM * NGLOB_wmax
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displ))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(veloc))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accel))
    if (ATTENUATION) then
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_wmax * N_SLS_wmax
      call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xx))
      call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_yy))
      call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xy))
      call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xz))
      call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_yz))
      call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_trace))

      ! strain not needed anymore, to save file diskspace - will be re-constructed based on b_displ...
      !local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_wmax
      !call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
      !                                 STRINGIFY_VAR(epsilondev_xx))
      !call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
      !                                 STRINGIFY_VAR(epsilondev_yy))
      !call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
      !                                 STRINGIFY_VAR(epsilondev_xy))
      !call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
      !                                 STRINGIFY_VAR(epsilondev_xz))
      !call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
      !                                 STRINGIFY_VAR(epsilondev_yz))
      !call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
      !                                 STRINGIFY_VAR(epsilondev_trace))
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    local_dim = NDIM * NGLOB_wmax
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displs_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(velocs_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accels_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displw_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(velocw_poroelastic))
    call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accelw_poroelastic))
  endif

  ! perform writing
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) then
    ! end step to indicate output is completed. ADIOS2 can do I/O
    call write_adios_end_step(myadios_fwd_file)
  else
    ! Reset the path to its original value to avoid bugs and write out arrays.
    call write_adios_perform(myadios_fwd_file)
  endif

  ! Close ADIOS handler to the restart file.
  if (do_close_file) then
    ! flushes all engines (makes sure i/o is all written out)
    call flush_adios_group_all(myadios_fwd_group)
    ! closes file
    call close_file_adios(myadios_fwd_file)
    ! re-sets flag
    is_initialized_fwd_group = .false.

    ! debug
    !if (myrank == 0) print *,'debug: undoatt adios: close file save forward step = ',iteration_on_subset_tmp, &
    !                         ' handle ',myadios_fwd_file
  endif

  end subroutine save_forward_arrays_undoatt_adios
