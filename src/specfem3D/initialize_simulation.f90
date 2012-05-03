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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine initialize_simulation()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_movie
  implicit none

  integer :: ier

  ! read the parameter file
  call read_parameter_file( NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        OCEANS,TOPOGRAPHY,ANISOTROPY,ABSORBING_CONDITIONS, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY)

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)))

  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_rank(myrank)

  ! checks flags
  call initialize_simulation_check()

  ! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_solver.txt',status='unknown')
  ! user output
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*) '**** Specfem 3-D Solver - MPI version f90 ****'
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*)
    write(IMAIN,*)
    if(FIX_UNDERFLOW_PROBLEM) write(IMAIN,*) 'Fixing slow underflow trapping problem using small initial field'
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NPROC,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',NPROC-1
    write(IMAIN,*)
    write(IMAIN,*) 'There is a total of ',NPROC,' slices'
    write(IMAIN,*)
    write(IMAIN,*) ' NDIM = ',NDIM
    write(IMAIN,*)
    write(IMAIN,*) ' NGLLX = ',NGLLX
    write(IMAIN,*) ' NGLLY = ',NGLLY
    write(IMAIN,*) ' NGLLZ = ',NGLLZ
    write(IMAIN,*)
    ! write information about precision used for floating-point operations
    if(CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',&
                   tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)
  endif

  ! reads in numbers of spectral elements and points for this process' domain
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',&
        action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error: could not open database '
    print*,'path: ',prname(1:len_trim(prname))//'external_mesh.bin'
    call exit_mpi(myrank,'error opening database')
  endif
  read(27) NSPEC_AB
  read(27) NGLOB_AB
  close(27)

  ! attenuation arrays size
  if( ATTENUATION ) then
    !pll
    NSPEC_ATTENUATION_AB = NSPEC_AB
  else
    ! if attenuation is off, set dummy size of arrays to one
    NSPEC_ATTENUATION_AB = 1
  endif
  ! needed for attenuation and/or kernel computations
  if( ATTENUATION .or. SIMULATION_TYPE == 3 ) then
    COMPUTE_AND_STORE_STRAIN = .true.
    NSPEC_STRAIN_ONLY = NSPEC_AB
  else
    COMPUTE_AND_STORE_STRAIN = .false.
    NSPEC_STRAIN_ONLY = 1
  endif

  ! anisotropy arrays size
  if( ANISOTROPY ) then
    NSPEC_ANISO = NSPEC_AB
  else
    ! if off, set dummy size
    NSPEC_ANISO = 1
  endif

  ! allocate arrays for storing the databases
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          xix(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          xiy(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          xiz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          etax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          etay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          etaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          gammax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          gammay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for databases'
  ! mesh node locations
  allocate(xstore(NGLOB_AB), &
          ystore(NGLOB_AB), &
          zstore(NGLOB_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for mesh nodes'
  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for material properties'
  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB), &
          ispec_is_elastic(NSPEC_AB), &
          ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for material flags'

  ! initializes adjoint simulations
  call initialize_simulation_adjoint()

  end subroutine initialize_simulation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_simulation_check()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_movie
  implicit none

  integer :: sizeprocs
  integer :: ier

  ! sizeprocs returns number of processes started
  ! (should be equal to NPROC)
  call world_size(sizeprocs)

  ! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) then
    if( myrank == 0 ) then
      write(IMAIN,*) 'error: number of processors supposed to run on: ',NPROC
      write(IMAIN,*) 'error: number of MPI processors actually run on: ',sizeprocs
      print*, 'error: number of processors supposed to run on: ',NPROC
      print*, 'error: number of MPI processors actually run on: ',sizeprocs      
    endif
    call exit_MPI(myrank,'wrong number of MPI processes')
  endif

  ! check that we have at least one source
  if(NSOURCES < 1) call exit_MPI(myrank,'need at least one source')

  ! check simulation type
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
        call exit_mpi(myrank,'SIMULATION_TYPE can only be 1, 2, or 3')

  ! check that optimized routines from Deville et al. (2002) can be used
  if( USE_DEVILLE_PRODUCTS) then
    if(NGLLX < 5 .or. NGLLY < 5 .or. NGLLZ < 5 .or. NGLLX > 10 .or. NGLLY > 10 .or. NGLLZ > 10) &
      stop 'Deville et al. (2002) routines can only be used if NGLLX = NGLLY = NGLLZ is in [5-10]'
  endif

  ! absorbing surfaces
  if( ABSORBING_CONDITIONS ) then
    ! for arbitrary orientation of elements, which face belongs to xmin,xmax,etc... -
    ! does it makes sense to have different NGLLX,NGLLY,NGLLZ?
    ! there is a problem with absorbing boundaries for faces with different NGLLX,NGLLY,NGLLZ values
    ! just to be sure for now..
    if( NGLLX /= NGLLY .and. NGLLY /= NGLLZ ) &
      stop 'ABSORBING_CONDITIONS must have NGLLX = NGLLY = NGLLZ'
  endif

  ! exclusive movie flags
  if( EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP ) then
    if( EXTERNAL_MESH_MOVIE_SURFACE .and. EXTERNAL_MESH_CREATE_SHAKEMAP ) &
      stop 'EXTERNAL_MESH_MOVIE_SURFACE and EXTERNAL_MESH_MOVIE_SURFACE cannot be both true'
    if( MOVIE_SURFACE ) &
      stop 'MOVIE_SURFACE cannot be used when EXTERNAL_MESH_MOVIE_SURFACE or EXTERNAL_MESH_CREATE_SHAKEMAP is true'
    if( CREATE_SHAKEMAP ) &
      stop 'CREATE_SHAKEMAP cannot be used when EXTERNAL_MESH_MOVIE_SURFACE or EXTERNAL_MESH_CREATE_SHAKEMAP is true'
  endif

  ! checks directories
  if( myrank == 0 ) then
    ! tests if OUTPUT_FILES directory exists
    call get_value_string(dummystring, 'OUTPUT_FILES', OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)))
    ! note: inquire behaves differently when using intel ifort or gfortran compilers
    !INQUIRE( FILE = dummystring(1:len_trim(dummystring))//'/.', EXIST = exists )
    open(IOUT,file=trim(dummystring)//'/dummy.txt',status='unknown',iostat=ier)
    if( ier /= 0 ) then
      print*,"OUTPUT_FILES directory does not work: ",trim(dummystring)
      call exit_MPI(myrank,'error OUTPUT_FILES directory')
    endif
    close(IOUT,status='delete')

    ! tests if LOCAL_PATH directory exists
    dummystring = adjustl(LOCAL_PATH)
    !INQUIRE( FILE = dummystring(1:len_trim(dummystring))//'/.', EXIST = exists )
    open(IOUT,file=trim(dummystring)//'/dummy.txt',status='unknown',iostat=ier)
    if( ier /= 0 ) then
      print*,"LOCAL_PATH directory does not work: ",trim(dummystring)
      call exit_MPI(myrank,'error LOCAL_PATH directory')
    endif
    close(IOUT,status='delete')
  endif

  end subroutine initialize_simulation_check
!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_simulation_adjoint()

! initialization for ADJOINT simulations

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  implicit none

  ! check simulation parameters
  if (SIMULATION_TYPE /= 1 .and. NSOURCES > 1000) &
    call exit_mpi(myrank, 'for adjoint simulations, NSOURCES <= 1000')

  ! snapshot file names: ADJOINT attenuation
  if (ATTENUATION .and. ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3)) &
    call create_name_database(prname_Q,myrank,LOCAL_PATH_Q)

  ! number of elements and points for adjoint arrays
  if( SIMULATION_TYPE == 3 ) then
    NSPEC_ADJOINT = NSPEC_AB
    NGLOB_ADJOINT = NGLOB_AB
  else
    ! dummy array size
    NSPEC_ADJOINT = 1
    NGLOB_ADJOINT =  1
  endif

  ! moho boundary
  if( SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3 ) then
    NSPEC_BOUN = NSPEC_AB
  else
    NSPEC_BOUN = 1
  endif

  end subroutine initialize_simulation_adjoint
