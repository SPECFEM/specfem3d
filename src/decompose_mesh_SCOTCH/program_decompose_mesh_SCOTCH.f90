!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

program pre_meshfem3D

  use decompose_mesh_SCOTCH,only: nparts,localpath_name, outputpath_name,&
                                  read_mesh_files, &
                                  check_valence, &
                                  scotch_partitioning, &
                                  write_mesh_databases
  implicit none
  integer :: i
  character(len=256) :: arg(3)

!  include './constants_decompose_mesh_SCOTCH.h'

! check usage
  do i=1,3
    call getarg(i,arg(i))
    if (i <= 3 .and. trim(arg(i)) == "") then
      print *, 'Usage: ./decompose_mesh_SCOTCH  nparts  input_directory output_directory'
      print *
      print *, '  where'
      print *, '      nparts = number of partitons'
      print *, '      input_directory = directory containing mesh files mesh_file,nodes_coords_file,..'
      print *, '      output_directory = directory for output files proc***_Databases'
      print *
      stop ' Reenter command line options'
    endif
  enddo

  read(arg(1),*) nparts
  localpath_name = arg(2)
  outputpath_name = arg(3)

! reads in (CUBIT) mesh files: mesh_file,nodes_coord_file, ...
  call read_mesh_files()

! checks valence of nodes
  call check_valence()

! partitions mesh
  call scotch_partitioning()

! writes out database files
  call write_mesh_databases()

end program pre_meshfem3D

