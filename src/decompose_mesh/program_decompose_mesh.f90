!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

program xdecompose_mesh

  use constants, only: MAX_STRING_LEN

  use decompose_mesh, only: nparts,localpath_name,outputpath_name,read_mesh_files,check_valence, &
                                  scotch_partitioning,write_mesh_databases,ADIOS_FOR_DATABASES

  implicit none

  integer :: i,myrank
  logical :: BROADCAST_AFTER_READ
  character(len=MAX_STRING_LEN) :: arg(3)

! check usage
  do i=1,3
    call get_command_argument(i,arg(i))
    if (i <= 3 .and. trim(arg(i)) == "") then
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
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! checks adios parameters
  if (ADIOS_FOR_DATABASES) then
    print *, 'Error: ADIOS_FOR_DATABASES set to .true. in Par_file'
    print *, 'ADIOS format for databases stored by xdecompose_mesh not implemented yet, please check your Par_file settings...'
    stop 'Error ADIOS_FOR_DATABASES setting; Reenter command line'
  endif

! reads in (CUBIT) mesh files: mesh_file,nodes_coord_file, ...
  call read_mesh_files()

! checks valence of nodes
  call check_valence()

! partitions mesh
  call scotch_partitioning()

! writes out database files
  call write_mesh_databases()

end program xdecompose_mesh

