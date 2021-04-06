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

module combine_vol_data_adios_mod

  use adios_read_mod

  implicit none

contains

!=============================================================================

!> Print help message.
  subroutine print_usage_adios()

  implicit none

  print *, 'Usage: '
  print *, '   xcombine_data start_slice end_slice varname var_file ' // &
           'mesh_file output_dir high/low-resolution'
  print *, '   or '
  print *, '   xcombine_data slice_list varname var_file mesh_file ' // &
           'output_dir high/low-resolution'
  print *
  print *, '* possible varnames are '
  print *, '   rho_vp, rho_vs, kappastore, mustore, alpha_kernel, etc'
  print *
  print *, '   that are stored in the local directory as ' // &
           'real(kind=CUSTOM_REAL) varname(NGLLX,NGLLY,NGLLZ,NSPEC_AB)  '
  print *, '   in var_file.bp'
  print *
  print *, '* mesh_files are used to link variable to the topology'
  print *, '* output_dir indicates where var_name.vtk will be written'
  print *, '* give 0 for low resolution and 1 for high resolution'
  print *

  stop ' Reenter command line options'

  end subroutine print_usage_adios

!=============================================================================

!> Interpret command line arguments
  subroutine read_args_adios(arg, MAX_NUM_NODES, node_list, num_node, &
                             var_name, value_file_name, mesh_file_name, &
                             outdir, ires)

  use constants, only: MAX_STRING_LEN

  implicit none
  ! Arguments
  character(len=*), intent(in) :: arg(:)
  integer, intent(in) :: MAX_NUM_NODES
  integer, intent(out) :: node_list(:)
  integer, intent(out) :: num_node, ires
  character(len=*), intent(out) :: var_name, value_file_name, mesh_file_name, outdir

  ! Variables
  character(len=MAX_STRING_LEN) :: sline
  integer :: it, iproc, proc1, proc2, ier, njunk

  if (command_argument_count() == 6) then
    num_node = 0
    open(unit = 20, file = trim(arg(1)), status = 'unknown',iostat = ier)
    if (ier /= 0) then
      print *,'Error opening ',trim(arg(1))
      stop
    endif
    do while (1 == 1)
      read(20,'(a)',iostat=ier) sline
      if (ier /= 0) exit
      read(sline,*,iostat=ier) njunk
      if (ier /= 0) exit
      num_node = num_node + 1
      if (num_node > MAX_NUM_NODES) &
          stop 'error number of slices exceeds MAX_NUM_NODES...'
      node_list(num_node) = njunk
    enddo
    close(20)
    var_name = arg(2)
    value_file_name = arg(3)
    mesh_file_name = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
  else if (command_argument_count() == 7) then
    read(arg(1),*) proc1
    read(arg(2),*) proc2
    do iproc = proc1, proc2
      it = iproc - proc1 + 1
      if (it > MAX_NUM_NODES) &
          stop 'error number of slices exceeds MAX_NUM_NODES...'
      node_list(it) = iproc
    enddo
    num_node = proc2 - proc1 + 1
    var_name = arg(3)
    value_file_name= arg(4)
    mesh_file_name = arg(5)
    outdir = arg(6)
    read(arg(7),*) ires
  else
    call print_usage_adios()
  endif

  end subroutine read_args_adios


!=============================================================================

!> Open ADIOS value and mesh files, read mode
  subroutine init_adios(value_file_name, mesh_file_name, value_handle, mesh_handle)

  use adios_manager_mod, only: adios_setup,comm_adios,ADIOS_VERBOSITY

  implicit none
  ! Parameters
  character(len=*), intent(in) :: value_file_name, mesh_file_name
  integer(kind=8), intent(out) :: value_handle, mesh_handle
  ! Variables
  integer :: ier
  integer :: comm

  ! initializes adios
  call adios_setup()

  ! gets MPI communicator
  comm = comm_adios

  call adios_read_init_method(ADIOS_READ_METHOD_BP, comm, ADIOS_VERBOSITY, ier)

  print *,'ADIOS opening mesh file: ',trim(mesh_file_name)
  call adios_read_open_file(mesh_handle, trim(mesh_file_name), 0, comm, ier)
  if (ier /= 0) call abort_mpi()

  print *,'ADIOS opening data file: ',trim(value_file_name)
  call adios_read_open_file(value_handle, trim(value_file_name), 0, comm, ier)
  if (ier /= 0) call abort_mpi()

  end subroutine init_adios


!=============================================================================

!> Open ADIOS value and mesh files, read mode
  subroutine clean_adios(value_handle, mesh_handle)

  use adios_manager_mod, only: adios_cleanup
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: value_handle, mesh_handle
  ! Variables
  integer :: ier

  call adios_read_close(mesh_handle,ier)
  call adios_read_close(value_handle,ier)

  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0 ) stop 'Error adios read finalize'

  ! finalize adios
  call adios_cleanup()

  end subroutine clean_adios

!=============================================================================

  subroutine read_scalars_adios_mesh(mesh_handle, iproc, NGLOB_AB, NSPEC_AB, ibool_offset, x_global_offset)

  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: mesh_handle
  integer, intent(in) :: iproc
  integer, intent(out) :: NGLOB_AB, NSPEC_AB
  integer(kind=8) :: ibool_offset, x_global_offset
  ! Variables
  integer(kind=8) :: sel
  integer :: ier

  call adios_selection_writeblock(sel, iproc)
  call adios_schedule_read(mesh_handle, sel, "nglob", 0, 1, NGLOB_AB, ier)
  call adios_schedule_read(mesh_handle, sel, "nspec", 0, 1, NSPEC_AB, ier)
  call adios_schedule_read(mesh_handle, sel, "ibool/offset", 0, 1, ibool_offset, ier)
  call adios_schedule_read(mesh_handle, sel, "x_global/offset", 0, 1, x_global_offset, ier)

  call adios_perform_reads(mesh_handle, ier)
  if (ier /= 0) call abort_mpi()

  end subroutine read_scalars_adios_mesh

!=============================================================================

  subroutine read_ibool_adios_mesh(mesh_handle, ibool_offset, NGLLX, NGLLY, NGLLZ, NSPEC_AB, ibool)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: mesh_handle,ibool_offset
  integer, intent(in) :: NGLLX, NGLLY, NGLLZ, NSPEC_AB
  integer, dimension(:,:,:,:), intent(inout) :: ibool
  ! Variables
  integer(kind=8), dimension(1) :: start, count_ad
  integer(kind=8) :: sel
  integer :: ier

  start(1) = ibool_offset
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(mesh_handle, sel, "ibool/array", 0, 1, ibool, ier)

  call adios_perform_reads(mesh_handle, ier)
  if (ier /= 0) call abort_mpi()

  end subroutine read_ibool_adios_mesh


!=============================================================================

  subroutine read_coordinates_adios_mesh(mesh_handle, x_global_offset, NGLOB_AB, xstore, ystore, zstore)

  use constants, only: CUSTOM_REAL

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: mesh_handle,x_global_offset
  integer, intent(in) :: NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(:), intent(inout) :: xstore, ystore, zstore
  ! Variables
  integer(kind=8), dimension(1) :: start, count_ad
  integer(kind=8) :: sel
  integer :: ier

  start(1) = x_global_offset
  count_ad(1) = NGLOB_AB
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(mesh_handle, sel, "x_global/array", 0, 1, xstore, ier)
  call adios_schedule_read(mesh_handle, sel, "y_global/array", 0, 1, ystore, ier)
  call adios_schedule_read(mesh_handle, sel, "z_global/array", 0, 1, zstore, ier)

  call adios_perform_reads(mesh_handle, ier)
  if (ier /= 0) call abort_mpi()

  end subroutine read_coordinates_adios_mesh

!=============================================================================

  subroutine read_double_values_adios(value_handle, var_name, ibool_offset, NSPEC_AB, dat)

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: value_handle,ibool_offset
  character(len=*), intent(in) :: var_name
  integer, intent(in) :: NSPEC_AB
  double precision, dimension(:,:,:,:), intent(inout) :: dat
  ! Variables
  integer(kind=8), dimension(1) :: start, count_ad
  integer(kind=8) :: sel
  integer :: ier

  start(1) = ibool_offset
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(value_handle, sel, trim(var_name) // "/array", 0, 1, dat, ier)

  call adios_perform_reads(value_handle, ier)
  if (ier /= 0) call abort_mpi()

  end subroutine read_double_values_adios

!=============================================================================

  subroutine read_float_values_adios(value_handle, var_name, ibool_offset, NSPEC_AB, dat)

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: value_handle,ibool_offset
  character(len=*), intent(in) :: var_name
  integer, intent(in) :: NSPEC_AB
  real, dimension(:,:,:,:), intent(inout) :: dat
  ! Variables
  integer(kind=8), dimension(1) :: start, count_ad
  integer(kind=8) :: sel
  integer :: ier

  start(1) = ibool_offset
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(value_handle, sel, trim(var_name) // "/array", 0, 1, dat, ier)

  call adios_perform_reads(value_handle, ier)
  if (ier /= 0) call abort_mpi()

  end subroutine read_float_values_adios

end module combine_vol_data_adios_mod
