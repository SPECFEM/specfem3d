!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
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

program pre_meshfem3D

  use constants, only: MAX_STRING_LEN
  use decompose_mesh,only: nparts,localpath_name, outputpath_name, &
                                  read_mesh_files, &
                                  check_valence, &
                                  scotch_partitioning, &
                                  write_mesh_databases, &
                                  DT, HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,NPROC, &
                                  NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP, &
                                  UTM_PROJECTION_ZONE,SIMULATION_TYPE,NGNOD,NGNOD2D, &
                                  NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                                  NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,MOVIE_TYPE, &
                                  MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                                  USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION, &
                                  ATTENUATION,USE_OLSEN_ATTENUATION, &
                                  APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,USE_FORCE_POINT_SOURCE, &
                                  STACEY_ABSORBING_CONDITIONS,SAVE_FORWARD,STACEY_INSTEAD_OF_FREE_SURFACE, &
                                  ANISOTROPY,SAVE_MESH_FILES,USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION, &
                                  LOCAL_PATH,TOMOGRAPHY_PATH,PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE, &
                                  f0_FOR_PML,IMODEL,SEP_MODEL_DIRECTORY,FULL_ATTENUATION_SOLID,TRACTION_PATH, &
                                  COUPLE_WITH_EXTERNAL_CODE,EXTERNAL_CODE_TYPE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none

  integer :: i
  character(len=MAX_STRING_LEN) :: arg(3)
  ! ADIOS parameters
  logical :: ADIOS_ENABLED, ADIOS_FOR_DATABASES, ADIOS_FOR_MESH, &
             ADIOS_FOR_FORWARD_ARRAYS, ADIOS_FOR_KERNELS

! check usage
  do i=1,3
    call get_command_argument(i,arg(i))
    if (i <= 3 .and. trim(arg(i)) == "") then
      print *, 'Usage: ./xdecompose_mesh  nparts  input_directory output_directory'
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

 ! needs local_path for mesh files
  call read_parameter_file(NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT,NGNOD,NGNOD2D, &
                           UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,TOMOGRAPHY_PATH, &
                           ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                           APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ANISOTROPY,STACEY_ABSORBING_CONDITIONS,MOVIE_TYPE, &
                           MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                           NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                           SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                           NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                           NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                           USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE, &
                           USE_RICKER_TIME_FUNCTION,OLSEN_ATTENUATION_RATIO,PML_CONDITIONS, &
                           PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL, SEP_MODEL_DIRECTORY, &
                           FULL_ATTENUATION_SOLID,TRACTION_PATH,COUPLE_WITH_EXTERNAL_CODE,EXTERNAL_CODE_TYPE, &
                           MESH_A_CHUNK_OF_THE_EARTH)

  ! gets adios parameters in Par_file
  call read_adios_parameters(ADIOS_ENABLED, ADIOS_FOR_DATABASES, &
                             ADIOS_FOR_MESH, ADIOS_FOR_FORWARD_ARRAYS, &
                             ADIOS_FOR_KERNELS)

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

end program pre_meshfem3D

