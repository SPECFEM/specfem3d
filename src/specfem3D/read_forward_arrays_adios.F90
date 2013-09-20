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
!
! United States and French Government Sponsorship Acknowledged.

!==============================================================================
subroutine read_forward_arrays_adios()

  use mpi
  use adios_read_mod

  use pml_par

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  real(kind=CUSTOM_REAL):: minl,maxl,min_all,max_all
  integer :: ier,inum

  character(len=256) :: database_name
  integer(kind=8) :: handle

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count_ad

  integer :: local_dim_potential_acoustic,            &
             local_dim_potential_dot_acoustic,        &
             local_dim_potential_dot_dot_acoustic,    &
             local_dim_displ,                         &
             local_dim_veloc,                         &
             local_dim_accel,                         &
             local_dim_R_xx,                          &
             local_dim_R_yy,                          &
             local_dim_R_xy,                          &
             local_dim_R_xz,                          &
             local_dim_R_yz,                          &
             local_dim_epsilondev_xx,                 &
             local_dim_epsilondev_yy,                 &
             local_dim_epsilondev_xy,                 &
             local_dim_epsilondev_xz,                 &
             local_dim_epsilondev_yz,                 &
             local_dim_R_trace,                       &
             local_dim_epsilondev_trace,              &
             local_dim_displs_poroelastic,            &
             local_dim_velocs_poroelastic,            &
             local_dim_accels_poroelastic,            &
             local_dim_displw_poroelastic,            &
             local_dim_velocw_poroelastic,            &
             local_dim_accelw_poroelastic

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  sel_num = 0

  database_name = adjustl(LOCAL_PATH)
  database_name = database_name(1:len_trim(database_name)) // "/forward_arrays.bp"

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                               "verbose=1", ier)
  call adios_read_open_file (handle, database_name, 0, MPI_COMM_WORLD, ier)

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  if (ACOUSTIC_SIMULATION) then
    call adios_get_scalar(handle, "potential_acoustic/local_dim",&
                          local_dim_potential_acoustic,ier)
    call adios_get_scalar(handle, "potential_dot_acoustic/local_dim",&
                          local_dim_potential_dot_acoustic,ier)
    call adios_get_scalar(handle, "potential_dot_dot_acoustic/local_dim",&
                          local_dim_potential_dot_dot_acoustic,ier)
  endif
  if (ELASTIC_SIMULATION) then
    call adios_get_scalar(handle, "displ/local_dim",&
                          local_dim_displ,ier)
    call adios_get_scalar(handle, "veloc/local_dim",&
                          local_dim_veloc,ier)
    call adios_get_scalar(handle, "accel/local_dim",&
                          local_dim_accel,ier)
    if (ATTENUATION) then
      call adios_get_scalar(handle, "R_xx/local_dim",&
                            local_dim_R_xx,ier)
      call adios_get_scalar(handle, "R_yy/local_dim",&
                            local_dim_R_yy,ier)
      call adios_get_scalar(handle, "R_xy/local_dim",&
                            local_dim_R_xy,ier)
      call adios_get_scalar(handle, "R_xz/local_dim",&
                            local_dim_R_xz,ier)
      call adios_get_scalar(handle, "R_yz/local_dim",&
                            local_dim_R_yz,ier)
      call adios_get_scalar(handle, "epsilondev_xx/local_dim",&
                            local_dim_epsilondev_xx,ier)
      call adios_get_scalar(handle, "epsilondev_yy/local_dim",&
                            local_dim_epsilondev_yy,ier)
      call adios_get_scalar(handle, "epsilondev_xy/local_dim",&
                            local_dim_epsilondev_xy,ier)
      call adios_get_scalar(handle, "epsilondev_xz/local_dim",&
                            local_dim_epsilondev_xz,ier)
      call adios_get_scalar(handle, "epsilondev_yz/local_dim",&
                            local_dim_epsilondev_yz,ier)
      if (FULL_ATTENUATION_SOLID) then
        call adios_get_scalar(handle, "R_trace/local_dim",&
                              local_dim_R_trace,ier)
        call adios_get_scalar(handle, "epsilondev_trace/local_dim",&
                              local_dim_epsilondev_trace,ier)
      endif
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    call adios_get_scalar(handle, "displs_poroelastic/local_dim",&
                          local_dim_displs_poroelastic,ier)
    call adios_get_scalar(handle, "velocs_poroelastic/local_dim",&
                          local_dim_velocs_poroelastic,ier)
    call adios_get_scalar(handle, "accels_poroelastic/local_dim",&
                          local_dim_accels_poroelastic,ier)
    call adios_get_scalar(handle, "displw_poroelastic/local_dim",&
                          local_dim_displw_poroelastic,ier)
    call adios_get_scalar(handle, "velocw_poroelastic/local_dim",&
                          local_dim_velocw_poroelastic,ier)
    call adios_get_scalar(handle, "accelw_poroelastic/local_dim",&
                          local_dim_accelw_poroelastic,ier)
  endif

  !-----------------------------------.
  ! Read arrays from forward_arrays.bp |
  !-----------------------------------'
  if( ACOUSTIC_SIMULATION ) then
    start(1) = local_dim_potential_acoustic * myrank
    count_ad(1) = NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "potential_acoustic/array", 0, 1, &
                             b_potential_acoustic, ier)

    start(1) = local_dim_potential_dot_acoustic * myrank
    count_ad(1) = NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "potential_dot_acoustic/array", 0, 1, &
                             b_potential_dot_acoustic, ier)

    start(1) = local_dim_potential_dot_dot_acoustic * myrank
    count_ad(1) = NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "potential_dot_dot_acoustic/array", 0, 1, &
                             b_potential_dot_dot_acoustic, ier)

  endif

  ! elastic wavefields
  if( ELASTIC_SIMULATION ) then
    start(1) = local_dim_displ * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "displ/array", 0, 1, &
                             b_displ, ier)

    start(1) = local_dim_veloc * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "veloc/array", 0, 1, &
                             b_veloc, ier)

    start(1) = local_dim_accel * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "accel/array", 0, 1, &
                             b_accel, ier)

    ! memory variables if attenuation
    if( ATTENUATION ) then
      start(1) = local_dim_R_xx * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "R_xx/array", 0, 1, &
                               b_R_xx, ier)

      start(1) = local_dim_R_yy * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "R_yy/array", 0, 1, &
                               b_R_yy, ier)

      start(1) = local_dim_R_xy * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "R_xy/array", 0, 1, &
                               b_R_xy, ier)

      start(1) = local_dim_R_xz * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "R_xz/array", 0, 1, &
                               b_R_xz, ier)

      start(1) = local_dim_R_yz * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB * N_SLS
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "R_yz/array", 0, 1, &
                               b_R_yz, ier)

      start(1) = local_dim_epsilondev_xx * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "epsilondev_xx/array", 0, 1, &
                               b_epsilondev_xx, ier)

      start(1) = local_dim_epsilondev_yy * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "epsilondev_yy/array", 0, 1, &
                               b_epsilondev_yy, ier)

      start(1) = local_dim_epsilondev_xy * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "epsilondev_xy/array", 0, 1, &
                               b_epsilondev_xy, ier)

      start(1) = local_dim_epsilondev_xz * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "epsilondev_xz/array", 0, 1, &
                               b_epsilondev_xz, ier)

      start(1) = local_dim_epsilondev_yz * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_ONLY
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "epsilondev_yz/array", 0, 1, &
                               b_epsilondev_yz, ier)

      if(FULL_ATTENUATION_SOLID) then
        start(1) = local_dim_R_trace * myrank
        count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB_kappa * N_SLS
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "R_trace/array", 0, 1, &
                                 b_R_trace, ier)

        start(1) = local_dim_epsilondev_trace * myrank
        count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_AB_kappa
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "epsilondev_trace/array", 0, 1, &
                                 b_epsilondev_trace, ier)
      endif
    endif
  endif

  ! poroelastic wavefields
  if( POROELASTIC_SIMULATION ) then
    start(1) = local_dim_displs_poroelastic * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "displs_poroelastic/array", 0, 1, &
                             b_displs_poroelastic, ier)

    start(1) = local_dim_velocs_poroelastic * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "velocs_poroelastic/array", 0, 1, &
                             b_velocs_poroelastic, ier)

    start(1) = local_dim_accels_poroelastic * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "accels_poroelastic/array", 0, 1, &
                             b_accels_poroelastic, ier)

    start(1) = local_dim_displw_poroelastic * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "displw_poroelastic/array", 0, 1, &
                             b_displw_poroelastic, ier)

    start(1) = local_dim_velocw_poroelastic * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "velocw_poroelastic/array", 0, 1, &
                             b_velocw_poroelastic, ier)

    start(1) = local_dim_accelw_poroelastic * myrank
    count_ad(1) = NDIM * NGLOB_ADJOINT
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "accelw_poroelastic/array", 0, 1, &
                             b_accelw_poroelastic, ier)
  endif

  !---------------------------------------------------------------.
  ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
  !---------------------------------------------------------------'
  call adios_perform_reads(handle, ier)
  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)

end subroutine read_forward_arrays_adios
