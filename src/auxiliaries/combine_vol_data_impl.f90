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

module combine_vtk_par

    ! global point data
    real,dimension(:),allocatable :: total_dat

    ! maximum number of slices
    integer,parameter :: MAX_NUM_NODES = 600

end module combine_vtk_par

!=============================================================================

module combine_vol_data_mod

  implicit none

contains

!=============================================================================

!> Print help message.
  subroutine print_usage()

  implicit none
  print *, 'Usage: '
  print *, '   xcombine_data start_slice end_slice filename input_dir ' // &
           'output_dir high/low-resolution'
  print *, '    or '
  print *, '   xcombine_data slice_list filename input_dir output_dir ' // &
           'high/low-resolution'
  print *
  print *, ' possible filenames are '
  print *, '   rho_vp, rho_vs, kappastore, mustore, alpha_kernel, etc'
  print *
  print *, '   that are stored in the local directory as ' // &
           'real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,NSPEC_AB)  '
  print *, '   in filename.bin'
  print *
  print *, ' files have been collected in input_dir, ' // &
           'output mesh file goes to output_dir '
  print *, ' give 0 for low resolution and 1 for high resolution'
  print *

  stop ' Reenter command line options'

  end subroutine print_usage

!=============================================================================

!> Interpret command line arguments
  subroutine read_args(arg, MAX_NUM_NODES, node_list, num_node, filename, indir, outdir, ires)

  use constants, only: MAX_STRING_LEN

  implicit none
  ! Arguments
  character(len=*), intent(in) :: arg(:)
  integer, intent(in) :: MAX_NUM_NODES
  integer, intent(out) :: node_list(:)
  integer, intent(out) :: num_node, ires
  character(len=*), intent(out) :: filename, indir, outdir
  ! Variables
  character(len=MAX_STRING_LEN) :: sline
  integer :: it, iproc, proc1, proc2, ier, njunk

  if (command_argument_count() == 5) then
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
      if (num_node > MAX_NUM_NODES) stop 'error number of slices exceeds MAX_NUM_NODES...'
      node_list(num_node) = njunk
    enddo
    close(20)
    filename = arg(2)
    indir= arg(3)
    outdir = arg(4)
    read(arg(5),*) ires
  else if (command_argument_count() == 6) then
    read(arg(1),*) proc1
    read(arg(2),*) proc2
    do iproc = proc1, proc2
      it = iproc - proc1 + 1
      if (it > MAX_NUM_NODES) stop 'error number of slices exceeds MAX_NUM_NODES...'
      node_list(it) = iproc
    enddo
    num_node = proc2 - proc1 + 1
    filename = arg(3)
    indir = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
  else
    call print_usage()
  endif

  end subroutine read_args

end module combine_vol_data_mod
