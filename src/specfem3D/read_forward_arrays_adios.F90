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


  subroutine read_forward_arrays_adios()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  use pml_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: database_name,group_name
  ! ADIOS variables
  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num,i
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count

  integer(kind=8) :: local_dim_potential_acoustic,local_dim_potential_dot_acoustic,local_dim_potential_dot_dot_acoustic, &
    local_dim_displ,local_dim_veloc,local_dim_accel, &
    local_dim_R_xx,local_dim_R_yy,local_dim_R_xy,local_dim_R_xz,local_dim_R_yz, &
    local_dim_epsilondev_xx,local_dim_epsilondev_yy,local_dim_epsilondev_xy, &
    local_dim_epsilondev_xz,local_dim_epsilondev_yz, &
    local_dim_R_trace, &
    local_dim_epsilondev_trace, &
    local_dim_displs_poroelastic,local_dim_velocs_poroelastic,local_dim_accels_poroelastic, &
    local_dim_displw_poroelastic,local_dim_velocw_poroelastic,local_dim_accelw_poroelastic

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  sel_num = 0

  database_name = get_adios_filename(trim(LOCAL_PATH) // "/forward_arrays")

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  reading forward array file: ',trim(database_name)
#if defined(USE_ADIOS)
    write(IMAIN,*) '  using ADIOS1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '  using ADIOS2 file format'
#endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initializes i/o group
  group_name = "SPECFEM3D_FORWARD_ARRAYS"
  call init_adios_group(myadios_group,group_name)

  ! setup the ADIOS library to read the file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,database_name)

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  if (ACOUSTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "potential_acoustic", &
                          local_dim_potential_acoustic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "potential_dot_acoustic", &
                          local_dim_potential_dot_acoustic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "potential_dot_dot_acoustic", &
                          local_dim_potential_dot_dot_acoustic)
  endif
  if (ELASTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "displ", local_dim_displ)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "veloc", local_dim_veloc)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "accel", local_dim_accel)
    if (ATTENUATION) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "R_xx", local_dim_R_xx)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "R_yy", local_dim_R_yy)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "R_xy", local_dim_R_xy)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "R_xz", local_dim_R_xz)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "R_yz", local_dim_R_yz)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "epsilondev_xx", local_dim_epsilondev_xx)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "epsilondev_yy", local_dim_epsilondev_yy)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "epsilondev_xy", local_dim_epsilondev_xy)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "epsilondev_xz", local_dim_epsilondev_xz)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "epsilondev_yz", local_dim_epsilondev_yz)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "R_trace", local_dim_R_trace)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "epsilondev_trace", local_dim_epsilondev_trace)
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "displs_poroelastic", local_dim_displs_poroelastic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "velocs_poroelastic", local_dim_velocs_poroelastic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "accels_poroelastic", local_dim_accels_poroelastic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "displw_poroelastic", local_dim_displw_poroelastic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "velocw_poroelastic", local_dim_velocw_poroelastic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "accelw_poroelastic", local_dim_accelw_poroelastic)
  endif

  !-----------------------------------.
  ! Read arrays from forward_arrays.bp |
  !-----------------------------------'
  if (ACOUSTIC_SIMULATION) then
    start(1) = local_dim_potential_acoustic * myrank
    count(1) = NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "potential_acoustic/array", b_potential_acoustic)

    start(1) = local_dim_potential_dot_acoustic * myrank
    count(1) = NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "potential_dot_acoustic/array", b_potential_dot_acoustic)

    start(1) = local_dim_potential_dot_dot_acoustic * myrank
    count(1) = NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "potential_dot_dot_acoustic/array", b_potential_dot_dot_acoustic)
  endif

  ! elastic wavefields
  if (ELASTIC_SIMULATION) then
    start(1) = local_dim_displ * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "displ/array", b_displ)

    start(1) = local_dim_veloc * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "veloc/array", b_veloc)

    start(1) = local_dim_accel * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "accel/array", b_accel)

    ! memory variables if attenuation
    if (ATTENUATION) then
      start(1) = local_dim_R_xx * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xx/array", b_R_xx)

      start(1) = local_dim_R_yy * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yy/array", b_R_yy)

      start(1) = local_dim_R_xy * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xy/array", b_R_xy)

      start(1) = local_dim_R_xz * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xz/array", b_R_xz)

      start(1) = local_dim_R_yz * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yz/array", b_R_yz)

      start(1) = local_dim_R_trace * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_trace/array", b_R_trace)

      start(1) = local_dim_epsilondev_xx * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     "epsilondev_xx/array", b_epsilondev_xx)

      start(1) = local_dim_epsilondev_yy * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     "epsilondev_yy/array", b_epsilondev_yy)

      start(1) = local_dim_epsilondev_xy * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     "epsilondev_xy/array", b_epsilondev_xy)

      start(1) = local_dim_epsilondev_xz * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     "epsilondev_xz/array", b_epsilondev_xz)

      start(1) = local_dim_epsilondev_yz * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     "epsilondev_yz/array", b_epsilondev_yz)

      start(1) = local_dim_epsilondev_trace * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     "epsilondev_trace/array", b_epsilondev_trace)
    endif
  endif

  ! poroelastic wavefields
  if (POROELASTIC_SIMULATION) then
    start(1) = local_dim_displs_poroelastic * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "displs_poroelastic/array", b_displs_poroelastic)

    start(1) = local_dim_velocs_poroelastic * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "velocs_poroelastic/array", b_velocs_poroelastic)

    start(1) = local_dim_accels_poroelastic * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "accels_poroelastic/array", b_accels_poroelastic)

    start(1) = local_dim_displw_poroelastic * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "displw_poroelastic/array", b_displw_poroelastic)

    start(1) = local_dim_velocw_poroelastic * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "velocw_poroelastic/array", b_velocw_poroelastic)

    start(1) = local_dim_accelw_poroelastic * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "accelw_poroelastic/array", b_accelw_poroelastic)
  endif

  !---------------------------------------------------------------.
  ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
  !---------------------------------------------------------------'
  call read_adios_perform(myadios_file)

  ! frees selection
  do i = 1,sel_num
    sel => selections(i)
    call delete_adios_selection(sel)
  enddo

  ! closes default file and finalizes read method
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,group_name)

  end subroutine read_forward_arrays_adios


!-------------------------------------------------------------------------------
!> \brief Read forward arrays for undo attenuation from an ADIOS file.

  subroutine read_forward_arrays_undoatt_adios(iteration_on_subset_tmp)

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  use adios_helpers_mod
  use manager_adios

  implicit none
  ! Arguments
  integer, intent(in) :: iteration_on_subset_tmp

  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name,group_name
  ! ADIOS variables
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: local_dim_potential_acoustic,local_dim_displ,local_dim_R_xx, &
                     local_dim_displs_poroelastic

  ! shorten the name of iteration variable and make it integer*8
  integer(kind=8) :: step
  integer :: t_tmp
  ! selection
  !integer(kind=8)         :: sel
  ! selections array
  integer(kind=8), dimension(8),target :: selections
  integer :: sel_num, i
  integer(kind=8), pointer :: sel => null()
  ! multiple/single file for storage of snapshots
  logical :: do_open_file,do_close_file

  ! selections array
  sel_num = 0

  ! file handling
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) then
    ! iterations here go down from N to 1, but ADIOS files has steps 0..N-1
    step = iteration_on_subset_tmp - 1

    ! single file for all steps
    do_open_file = .false.
    do_close_file = .false.

    ! single file
    file_name = get_adios_filename(trim(LOCAL_PATH) // "/forward_arrays_undoatt",ADIOS2_ENGINE_UNDO_ATT)

    group_name = "SPECFEM3D_FORWARD_ARRAYS_UNDOATT"

    ! open file at first call of this routine
    if (.not. is_initialized_fwd_group) do_open_file = .true.
    ! close at last step (step counting down from N-1 to 0)
    if (step == 0) do_close_file = .true.

  else
    ! single ADIOS files for each subset with step 0 entry
    step = 0

    ! for each step a single file
    do_open_file = .true.
    do_close_file = .true.

    ! files for each iteration step
    write(file_name,'(a,a,i6.6)') trim(LOCAL_PATH), '/save_frame_at',iteration_on_subset_tmp
    file_name = get_adios_filename(trim(file_name))

    write(group_name, '(a, i6)') "SPECFEM3D_FORWARD_ARRAYS_UNDOATT", iteration_on_subset_tmp
  endif

  ! debug
  !if (myrank == 0) print *,'debug: undoatt adios: read forward step = ',step,'open/close file',do_open_file,do_close_file

  ! opens file for reading
  if (do_open_file) then
    ! creates adios group
    call init_adios_group_undo_att(myadios_fwd_group,group_name)

    ! opens adios file
    call open_file_adios_read(myadios_fwd_file,myadios_fwd_group,file_name)

    ! debug
    !if (myrank == 0) call show_adios_file_variables(myadios_fwd_file,myadios_fwd_group,file_name)

    ! sets flag
    is_initialized_fwd_group = .true.
  endif

  ! daniel todo: see if calling read_adios_perform() only once at the end would increase performance.
  !              this will need to have sel variables stored in an array.
  !
  !              note also that reading scalars is always blocking/synchronized.
  !
  ! checks if correct snapshot number of wavefields
  call read_adios_scalar(myadios_fwd_file,myadios_fwd_group,myrank,"iteration",t_tmp,step)
  ! checks step and iteration number
  if (t_tmp /= iteration_on_subset_tmp) then
    print *,'Error: invalid iteration step found in reading undoatt arrays: found ',t_tmp,' instead of ',iteration_on_subset_tmp
    call exit_mpi(myrank,'Invalid iteration step read in read_forward_arrays_undoatt_adios() routine')
  endif

  ! gets the 'chunks' sizes
  !
  ! note: for SPECFEM3D_Cartesian, the different (parallel) slices/chunks can have different sizes after partitioning.
  !       this is simpler for SPECFEM3D_GLOBE where all slices have the same number of elements after meshing.
  !
  !       here, we need first to read the overall maximum size value given by local_dim_* to estimate and select the bounding box
  !       for reading in the corresponding data arrays.
  !       we can assume that the local_dim_* is the same for displ/veloc/accel wavefields, for potential/potential_dot/.. etc.
  !       and only need to read in one local_dim_* value for similar arrays.
  if (ACOUSTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_fwd_file, myadios_fwd_group, myrank, "potential_acoustic", &
                                     local_dim_potential_acoustic)
  endif
  if (ELASTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_fwd_file, myadios_fwd_group, myrank, "displ", local_dim_displ)
    if (ATTENUATION) then
      call read_adios_scalar_local_dim(myadios_fwd_file, myadios_fwd_group, myrank, &
                                       "R_xx", local_dim_R_xx)
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_fwd_file, myadios_fwd_group, myrank, &
                                     "displs_poroelastic", local_dim_displs_poroelastic)
  endif

  ! reads in arrays
  ! acoustic wavefield
  if (ACOUSTIC_SIMULATION) then
    start(1) = local_dim_potential_acoustic * myrank
    count(1) = NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "potential_acoustic/array", b_potential_acoustic, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "potential_dot_acoustic/array", b_potential_dot_acoustic, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "potential_dot_dot_acoustic/array", b_potential_dot_dot_acoustic, step)
  endif

  ! elastic wavefield
  if (ELASTIC_SIMULATION) then
    start(1) = local_dim_displ * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "displ/array", b_displ, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "veloc/array", b_veloc, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "accel/array", b_accel, step)

    ! memory variables if attenuation
    if (ATTENUATION) then
      start(1) = local_dim_R_xx * myrank
      count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "R_xx/array", b_R_xx, step)
      call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "R_yy/array", b_R_yy, step)
      call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "R_xy/array", b_R_xy, step)
      call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "R_xz/array", b_R_xz, step)
      call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "R_yz/array", b_R_yz, step)
      call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, "R_trace/array", b_R_trace, step)

      ! strain not needed anymore, to save file diskspace - will be re-constructed based on b_displ...
      !start(1) = local_dim_epsilondev_xx * myrank
      !count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      !sel_num = sel_num+1
      !sel => selections(sel_num)
      !call set_selection_boundingbox(sel, start, count)
      !
      !call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
      !                               "epsilondev_xx/array", b_epsilondev_xx, step)
      !call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
      !                               "epsilondev_yy/array", b_epsilondev_yy, step)
      !call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
      !                               "epsilondev_xy/array", b_epsilondev_xy, step)
      !call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
      !                               "epsilondev_xz/array", b_epsilondev_xz, step)
      !call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
      !                               "epsilondev_yz/array", b_epsilondev_yz, step)
      !call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
      !                               "epsilondev_trace/array", b_epsilondev_trace, step)
    endif
  endif

  ! poroelastic wavefields
  if (POROELASTIC_SIMULATION) then
    start(1) = local_dim_displs_poroelastic * myrank
    count(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "displs_poroelastic/array", b_displs_poroelastic, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "velocs_poroelastic/array", b_velocs_poroelastic, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "accels_poroelastic/array", b_accels_poroelastic, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "displw_poroelastic/array", b_displw_poroelastic, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "velocw_poroelastic/array", b_velocw_poroelastic, step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "accelw_poroelastic/array", b_accelw_poroelastic, step)
  endif

  ! perform actual reading
  call read_adios_perform(myadios_fwd_file)

  ! free selection structures
  do i = 1, sel_num
    sel => selections(i)
    call delete_adios_selection(sel)
  enddo

  ! Close ADIOS handler to the restart file.
  ! note: for single file, only close at the very end
  if (do_close_file) then
    call close_file_adios_read_and_finalize_method(myadios_fwd_file)
    call delete_adios_group(myadios_fwd_group,group_name)

    is_initialized_fwd_group = .false.

    ! debug
    !if (myrank == 0) print *,'debug: undoatt adios: close file read forward step = ',step,' handle ',myadios_fwd_file
  endif

  end subroutine read_forward_arrays_undoatt_adios
