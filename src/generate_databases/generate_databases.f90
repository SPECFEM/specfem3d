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
!
! United States and French Government Sponsorship Acknowledged.

!=============================================================================!
!                                                                             !
!  generate_databases produces a spectral element grid                        !
!  for a local or regional model.                                             !
!  The mesher uses the UTM projection                                         !
!                                                                             !
!=============================================================================!
!
! Please find in the header of specfem3D.F90 further code informations.
!
! ************** PROGRAM STARTS HERE **************

  program xgenerate_databases

  use manager_adios
  use generate_databases_par

  implicit none

  include 'version.fh'

  ! local parameters
  ! timing
  double precision, external :: wtime

  ! MPI initialization
  call init_mpi()

  ! sizeprocs returns number of processes started (should be equal to NPROC).
  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_generate_databases.txt',status='unknown')

  ! get MPI starting time
  time_start = wtime()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '*****************************************'
    write(IMAIN,*) '*** Specfem3D MPI database generation ***'
    write(IMAIN,*) '*****************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Running Git package version of the code: ', git_package_version
    write(IMAIN,*) 'which is Git ', git_commit_version
    write(IMAIN,*) 'dating ', git_date_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! read the parameter file
  call read_parameters()

  ! reads topography and bathymetry file
  call read_topography()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '************************************'
    write(IMAIN,*) 'reading partition files in the model'
    write(IMAIN,*) '************************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! Initialize ADIOS I/O
  if (ADIOS_ENABLED) then
    call initialize_adios()
  endif

  ! reads Databases files
  if (ADIOS_FOR_DATABASES) then
    call read_partition_files_adios()
  else
    call read_partition_files()
  endif

  ! external mesh creation
  call setup_mesh()

  ! finalize mesher
  call finalize_databases()

  if (ADIOS_ENABLED) then
    call finalize_adios()
  endif

  ! MPI finish
  call finalize_mpi()

  end program xgenerate_databases

