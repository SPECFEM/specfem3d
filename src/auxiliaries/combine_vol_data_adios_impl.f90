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

module combine_vol_data_adios_mod

  use adios_helpers_mod
  use manager_adios

  implicit none

contains

!=============================================================================

!> Print help message.
  subroutine print_usage_adios()

  implicit none

  print *, 'Usage: '
  print *, '   xcombine_data start_slice end_slice varname var_file mesh_file output_dir high/low-resolution'
  print *, '   or '
  print *, '   xcombine_data slice_list varname var_file mesh_file output_dir high/low-resolution'
  print *
  print *, ' with'
  print *, '   start_slice/end_slice - start/end range of slice numbers to combine'
  print *, '   slice_list            - text file containing slice numbers to combine (or use name "all" for all slices)'
  print *, '   varname      - possible varnames are: '
  print *, '                    rho, vp, vs, kappastore, mustore, alpha_kl, etc.'
  print *, '   var_file     - datafile that holds array, as real(kind=CUSTOM_REAL):: varname(NGLLX,NGLLY,NGLLZ,NSPEC),'
  print *, '                  (e.g. OUTPUT_FILES/kernels.bp)'
  print *, '   mesh_file    - are used to link variable to the topology (e.g. DATABASES_MPI/external_mesh.bp)'
  print *, '   output_dir   - indicates where var_name.vtk will be written'
  print *, '   high/low res - give 0 for low resolution and 1 for high resolution'

  print *

  stop ' Reenter command line options'

  end subroutine print_usage_adios

!=============================================================================

!> Interpret command line arguments
  subroutine read_args_adios(arg, MAX_NUM_NODES, node_list, num_node, &
                             var_name, value_file_name, mesh_file_name, &
                             outdir, ires, NPROCTOT)

  use constants, only: IIN,MAX_STRING_LEN

  implicit none
  ! Arguments
  character(len=*), intent(in) :: arg(:)
  integer, intent(in) :: MAX_NUM_NODES
  integer, intent(out) :: node_list(:)
  integer, intent(out) :: num_node, ires
  character(len=*), intent(out) :: var_name, value_file_name, mesh_file_name, outdir
  integer, intent(in) :: NPROCTOT

  ! Variables
  character(len=MAX_STRING_LEN) :: sline,slice_list_name
  integer :: i, it, iproc, proc1, proc2, ier, njunk

  if (command_argument_count() == 6) then
    slice_list_name = arg(1)
    var_name = arg(2)
    value_file_name = arg(3)
    mesh_file_name = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
    ! gets slice list
    if (trim(slice_list_name) == 'all' .or. trim(slice_list_name) == '-1') then
      ! combines all slices
      num_node = 0
      do i = 0,NPROCTOT-1
        num_node = num_node + 1
        if (num_node > MAX_NUM_NODES ) stop 'Error number of slices exceeds MAX_NUM_NODES...'
        node_list(num_node) = i
      enddo
    else
      ! reads in slices files
      open(unit = IIN, file = trim(slice_list_name), status = 'unknown',iostat = ier)
      if (ier /= 0) then
        print *,'Error opening slice file ',trim(slice_list_name)
        stop
      endif
      num_node = 0
      do while (1 == 1)
        read(IIN,'(a)',iostat=ier) sline
        if (ier /= 0) exit
        read(sline,*,iostat=ier) njunk
        if (ier /= 0) exit
        num_node = num_node + 1
        if (num_node > MAX_NUM_NODES) &
            stop 'error number of slices exceeds MAX_NUM_NODES...'
        node_list(num_node) = njunk
      enddo
      close(IIN)
    endif
  else if (command_argument_count() == 7) then
    read(arg(1),*) proc1
    read(arg(2),*) proc2
    var_name = arg(3)
    value_file_name= arg(4)
    mesh_file_name = arg(5)
    outdir = arg(6)
    read(arg(7),*) ires
    ! sets slice range
    do iproc = proc1, proc2
      it = iproc - proc1 + 1
      if (it > MAX_NUM_NODES) &
          stop 'error number of slices exceeds MAX_NUM_NODES...'
      node_list(it) = iproc
    enddo
    num_node = proc2 - proc1 + 1
  else
    call print_usage_adios()
  endif

  end subroutine read_args_adios


!=============================================================================

!> Open ADIOS value and mesh files, read mode
  subroutine init_adios(value_file_name, mesh_file_name)

  implicit none
  ! Parameters
  character(len=*), intent(in) :: value_file_name, mesh_file_name

  ! initializes adios
  call initialize_adios()

  ! initializes read method and opens mesh file (using default handle)
  print *,'ADIOS opening mesh file: ',trim(mesh_file_name)
  call init_adios_group(myadios_group,"MeshReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,mesh_file_name)

  ! opens second adios file for reading data values
  print *,'ADIOS opening data file: ',trim(value_file_name)
  call init_adios_group(myadios_val_group,"ValReader")
  call open_file_adios_read(myadios_val_file,myadios_val_group,value_file_name)

  end subroutine init_adios


!=============================================================================

!> Open ADIOS value and mesh files, read mode
  subroutine clean_adios()

  implicit none

  ! closes file with data values
  call close_file_adios_read(myadios_val_file)

  ! closes default file and finalizes read method
  call close_file_adios_read_and_finalize_method(myadios_file)

  ! finalizes adios
  call finalize_adios()

  end subroutine clean_adios

!=============================================================================

  subroutine read_scalars_adios_mesh(iproc, NGLOB_AB, NSPEC_AB)

  implicit none

  ! Parameters
  integer, intent(in) :: iproc
  integer, intent(out) :: NGLOB_AB, NSPEC_AB

  ! reads nglob & nspec
  ! the following adios calls allow to have different nglob/nspec values for different processes (iproc)
  call read_adios_scalar(myadios_file,myadios_group,iproc,"nglob",NGLOB_AB)
  call read_adios_scalar(myadios_file,myadios_group,iproc,"nspec",NSPEC_AB)

  end subroutine read_scalars_adios_mesh

!=============================================================================

  subroutine read_ibool_adios_mesh(iproc, NGLLX, NGLLY, NGLLZ, NSPEC_AB, ibool)

  implicit none
  ! Parameters
  integer, intent(in) :: iproc, NGLLX, NGLLY, NGLLZ, NSPEC_AB
  integer, dimension(:,:,:,:), intent(inout) :: ibool
  ! Variables
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel
  integer(kind=8) :: offset_ibool

  call read_adios_scalar(myadios_file,myadios_group,iproc,"ibool/offset",offset_ibool)

  start(1) = offset_ibool
  count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "ibool/array", ibool)

  ! perform reading
  call read_adios_perform(myadios_file)

  call delete_adios_selection(sel)

  end subroutine read_ibool_adios_mesh


!=============================================================================

  subroutine read_coordinates_adios_mesh(iproc, NGLOB_AB, xstore, ystore, zstore)

  use constants, only: CUSTOM_REAL

  implicit none
  ! Parameters
  integer, intent(in) :: iproc, NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(:), intent(inout) :: xstore, ystore, zstore
  ! Variables
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel
  integer(kind=8) :: offset_coord

  call read_adios_scalar(myadios_file,myadios_group,iproc,"x_global/offset",offset_coord)

  start(1) = offset_coord
  count(1) = NGLOB_AB
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "x_global/array", xstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "y_global/array", ystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "z_global/array", zstore)

  ! perform reading
  call read_adios_perform(myadios_file)

  call delete_adios_selection(sel)

  end subroutine read_coordinates_adios_mesh

!=============================================================================

  subroutine read_values_adios(var_name, iproc, NSPEC_AB, data)

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  implicit none
  ! Parameters
  character(len=*), intent(in) :: var_name
  integer, intent(in) :: iproc, NSPEC_AB
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(inout) :: data
  ! Variables
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel
  integer(kind=8) :: offset

  ! gets data values
  if (.true.) then
    ! default
    ! assumes GLL type array size (NGLLX,NGLLY,NGLLZ,nspec)
    call read_adios_array(myadios_val_file,myadios_val_group,iproc,NSPEC_AB,trim(var_name),data)
  else
    ! reads in data offset
    call read_adios_scalar(myadios_val_file,myadios_val_group,iproc,trim(var_name) // "/offset",offset)

    ! reads in data array
    start(1) = offset
    count(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
    call set_selection_boundingbox(sel , start, count)

    call read_adios_schedule_array(myadios_val_file, myadios_val_group, sel, start, count, trim(var_name) // "/array", data)

    ! perform reading
    call read_adios_perform(myadios_val_file)

    call delete_adios_selection(sel)
  endif

  end subroutine read_values_adios

end module combine_vol_data_adios_mod
