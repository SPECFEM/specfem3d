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

program xdecompose_mesh

  use constants, only: MAX_STRING_LEN

  use decompose_mesh_par, only: nparts,localpath_name,outputpath_name,ADIOS_FOR_DATABASES

  implicit none

  integer :: i,myrank
  logical :: BROADCAST_AFTER_READ
  character(len=MAX_STRING_LEN) :: arg(3)

  ! user output
  print *
  print *,'**********************'
  print *,'Serial mesh decomposer'
  print *,'**********************'
  print *

  ! check usage
  do i = 1,3
    call get_command_argument(i,arg(i))
    if (i <= 3 .and. trim(arg(i)) == '') then
      print *, 'Usage: ./xdecompose_mesh  nparts  input_directory output_directory'
      print *
      print *, '  where'
      print *, '      nparts = number of partitions'
      print *, '      input_directory = directory containing mesh files mesh_file,nodes_coords_file,..'
      print *, '      output_directory = directory for output files proc***_Databases'
      print *
      stop ' Reenter command line options'
    endif
  enddo

  read(arg(1),*) nparts
  localpath_name = arg(2)
  outputpath_name = arg(3)

 ! needs local_path for mesh files
  myrank = 0
  BROADCAST_AFTER_READ = .false.
  call read_parameter_file(BROADCAST_AFTER_READ)

  ! checks adios parameters
  if (ADIOS_FOR_DATABASES) then
    ! note: this decomposer runs as a single, serial process.
    !       writing out ADIOS files would require a parallel section, for each process - not clear yet how to do that...
    print *, 'Error: ADIOS_FOR_DATABASES set to .true. in Par_file'
    print *, 'ADIOS format for databases stored by xdecompose_mesh not implemented yet, please check your Par_file settings...'
    stop 'Error ADIOS_FOR_DATABASES setting; Reenter command line'
  endif

  ! reads in (CUBIT) mesh files: mesh_file,nodes_coord_file, ...
  call read_mesh_files()

  ! checks valence of nodes
  call check_valence()

  ! sets up elements for local time stepping
  call lts_setup_elements()

  ! partitions mesh (using scotch, metis, or patoh partitioners via constants.h)
  call decompose_mesh()

  ! writes out database files
  call write_mesh_databases()

  ! user output
  print *
  print *,'finished successfully'
  print *

end program xdecompose_mesh

