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


  subroutine initialize_simulation()

  use manager_adios
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_movie

  implicit none

  include 'version.fh'

  integer :: ier
  logical :: BROADCAST_AFTER_READ
  character(len=MAX_STRING_LEN) :: path_to_add
  logical :: simul_run_flag

  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_rank(myrank)

  ! open main output file, only written to by process 0
  if (myrank == 0) then
    if (IMAIN /= ISTANDARD_OUTPUT) then
      open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_solver.txt',status='unknown',action='write',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file output_solver.txt for writing output info')
    endif

    write(IMAIN,*) '**********************************************'
    write(IMAIN,*) '**** Specfem 3-D Solver - MPI version f90 ****'
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Running Git package version of the code: ', git_package_version
    write(IMAIN,*) 'which is Git ', git_commit_version
    write(IMAIN,*) 'dating ', git_date_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! read the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(BROADCAST_AFTER_READ)

  ! checks flags
  call initialize_simulation_check()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*)
    if (FIX_UNDERFLOW_PROBLEM) write(IMAIN,*) 'Fixing slow underflow trapping problem using small initial field'
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
    if (CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    if (FORCE_VECTORIZATION_VAL) write(IMAIN,*) 'using force vectorization'
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ', &
                   tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)

    write(IMAIN,'(a)',advance='no') ' velocity model: '
    select case (IMODEL)
    case (IMODEL_DEFAULT)
    write(IMAIN,'(a)',advance='yes') '  default '
    case (IMODEL_GLL)
    write(IMAIN,'(a)',advance='yes') '  gll'
    case (IMODEL_1D_PREM)
    write(IMAIN,'(a)',advance='yes') '  1d_prem'
    case (IMODEL_1D_CASCADIA)
    write(IMAIN,'(a)',advance='yes') '  1d_cascadia'
    case (IMODEL_1D_SOCAL)
    write(IMAIN,'(a)',advance='yes') '  1d_socal'
    case (IMODEL_SALTON_TROUGH)
    write(IMAIN,'(a)',advance='yes') '  salton_trough'
    case (IMODEL_TOMO)
    write(IMAIN,'(a)',advance='yes') '  tomo'
    case (IMODEL_USER_EXTERNAL)
    write(IMAIN,'(a)',advance='yes') '  external'
    case (IMODEL_IPATI)
    write(IMAIN,'(a)',advance='yes') '  ipati'
    case (IMODEL_IPATI_WATER)
    write(IMAIN,'(a)',advance='yes') '  ipati_water'
    case (IMODEL_SEP)
    write(IMAIN,'(a)',advance='yes') '  SEP'
    case (IMODEL_COUPLED)
    write(IMAIN,'(a)',advance='yes') '  model coupled with injection method'
    end select

    write(IMAIN,*)
    call flush_IMAIN()
  endif
  ! synchronizes processes
  call synchronize_all()

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    simul_run_flag = .true.
  else
    simul_run_flag = .false.
  endif

  ! initializes ADIOS
  if (ADIOS_ENABLED) then
    call initialize_adios()
  endif

  if ((SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) .and. READ_ADJSRC_ASDF) then
    call asdf_setup(current_asdf_handle, path_to_add, simul_run_flag)
  endif

  ! reads in numbers of spectral elements and points for the part of the mesh handled by this process
  call create_name_database(prname,myrank,LOCAL_PATH)

  ! read the value of NSPEC_AB, NGLOB_AB and NSPEC_IRREGULAR
  ! (we need it to define some array sizes below)
  if (ADIOS_FOR_MESH) then
    call read_mesh_for_init_ADIOS()
  else
    call read_mesh_for_init()
  endif

  ! attenuation arrays size
  if (ATTENUATION) then
    NSPEC_ATTENUATION_AB = NSPEC_AB
  else
    ! if attenuation is off, set dummy size of arrays to one
    NSPEC_ATTENUATION_AB = 1
  endif

  ! needed for attenuation and/or kernel computations
  if (ATTENUATION .or. SIMULATION_TYPE == 3) then
    COMPUTE_AND_STORE_STRAIN = .true.
    NSPEC_STRAIN_ONLY = NSPEC_AB
  else
    COMPUTE_AND_STORE_STRAIN = .false.
    NSPEC_STRAIN_ONLY = 1
  endif

  ! anisotropy arrays size
  if (ANISOTROPY) then
    NSPEC_ANISO = NSPEC_AB
  else
    ! if off, set dummy size
    NSPEC_ANISO = 1
  endif

  allocate(irregular_element_number(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2388')
  if (ier /= 0) stop 'error allocating arrays for irregular element numbering'
  irregular_element_number(:) = 0

  ! allocate arrays for storing the databases
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2389')
  if (ier /= 0) stop 'error allocating ibool'
  ibool(:,:,:,:) = 0

  if (NSPEC_IRREGULAR > 0) then
     allocate(xixstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2390')
     allocate(xiystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2391')
     allocate(xizstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2392')
     allocate(etaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2393')
     allocate(etaystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2394')
     allocate(etazstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2395')
     allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2396')
     allocate(gammaystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2397')
     allocate(gammazstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2398')
     allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2399')
  else
    allocate(xixstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2400')
    allocate(xiystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2401')
    allocate(xizstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2402')
    allocate(etaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2403')
    allocate(etaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2404')
    allocate(etazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2405')
    allocate(gammaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2406')
    allocate(gammaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2407')
    allocate(gammazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2408')
    allocate(jacobianstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2409')
  endif
  xixstore(:,:,:,:) = 0.0_CUSTOM_REAL; xiystore(:,:,:,:) = 0.0_CUSTOM_REAL; xizstore(:,:,:,:) = 0.0_CUSTOM_REAL
  etaxstore(:,:,:,:) = 0.0_CUSTOM_REAL; etaystore(:,:,:,:) = 0.0_CUSTOM_REAL; etazstore(:,:,:,:) = 0.0_CUSTOM_REAL
  gammaxstore(:,:,:,:) = 0.0_CUSTOM_REAL; gammaystore(:,:,:,:) = 0.0_CUSTOM_REAL; gammazstore(:,:,:,:) = 0.0_CUSTOM_REAL

  ! mesh node locations
  allocate(xstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2410')
  allocate(ystore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2411')
  allocate(zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2412')
  if (ier /= 0) stop 'error allocating arrays for mesh nodes'
  xstore(:) = 0.0_CUSTOM_REAL; ystore(:) = 0.0_CUSTOM_REAL; zstore(:) = 0.0_CUSTOM_REAL

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2413')
  allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2414')
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating rho array 2414')
  if (ier /= 0) stop 'error allocating arrays for material properties'
  kappastore(:,:,:,:) = 0.0_CUSTOM_REAL; mustore(:,:,:,:) = 0.0_CUSTOM_REAL; rhostore(:,:,:,:) = 0.0_CUSTOM_REAL

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2415')
  allocate(ispec_is_elastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2416')
  allocate(ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2417')
  if (ier /= 0) stop 'error allocating arrays for material flags'
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! initializes adjoint simulations
  call initialize_simulation_adjoint()

  ! initializes GPU cards
  if (GPU_MODE) call initialize_GPU()

  ! output info for possible OpenMP
  call init_openmp()

  ! synchronizes processes
  call synchronize_all()

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

  integer :: ier,ioutputs

  character(len=MAX_STRING_LEN) :: HEADER_FILE
  logical :: ABSORB_FREE_SURFACE_VAL

  NAMELIST/MESHER/ABSORB_FREE_SURFACE_VAL

  ! sizeprocs returns number of processes started
  ! (should be equal to NPROC)
  call world_size(sizeprocs)

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      write(IMAIN,*) 'error: number of processors supposed to run on: ',NPROC
      write(IMAIN,*) 'error: number of MPI processors actually run on: ',sizeprocs
      print *
      print *, 'error specfem3D: number of processors supposed to run on: ',NPROC
      print *, 'error specfem3D: number of MPI processors actually run on: ',sizeprocs
      print *
    endif
    call exit_MPI(myrank,'wrong number of MPI processes')
  endif

  ! check that we have at least one source
  if (NSOURCES < 1 .and. (.not. HAS_FINITE_FAULT_SOURCE .and. .not. INVERSE_FWI_FULL_PROBLEM)) &
    call exit_MPI(myrank,'need at least one source')

  ! check simulation type
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    call exit_mpi(myrank,'SIMULATION_TYPE can only be 1, 2, or 3')

  ! gravity only on GPU supported
  if (.not. GPU_MODE .and. GRAVITY) stop 'GRAVITY only supported in GPU mode'

  if (NGLLX /= NGLLY .or. NGLLY /= NGLLZ) stop 'Methods that can handle unstructured meshes require NGLLX = NGLLY = NGLLZ'

  ! absorbing surfaces
  if (STACEY_ABSORBING_CONDITIONS) then
    if (PML_CONDITIONS) then
      print *, 'please modify Par_file and recompile solver'
      stop 'STACEY_ABSORBING_CONDITIONS and PML_CONDITIONS are both set to .true.'
    else if (PML_INSTEAD_OF_FREE_SURFACE) then
      print *, 'please modify Par_file and recompile solver'
      stop 'PML_INSTEAD_OF_FREE_SURFACE = .true. is incompatible with STACEY_ABSORBING_CONDITIONS = .true.'
    endif
  else
    if (STACEY_INSTEAD_OF_FREE_SURFACE) then
      print *, 'please modify Par_file and recompile solver'
      stop 'STACEY_ABSORBING_CONDITIONS must be activated when STACEY_INSTEAD_OF_FREE_SURFACE is set to .true.'
    endif
  endif

  if (PML_CONDITIONS) then
    if (STACEY_INSTEAD_OF_FREE_SURFACE) then
      print *, 'please modify Par_file and recompile solver'
      stop 'STACEY_INSTEAD_OF_FREE_SURFACE = .true. is incompatible with PML_CONDITIONS = .true.'
    !else if (.not. SUPPRESS_UTM_PROJECTION) then
    !  print *, 'please modify Par_file and recompile solver'
    !  stop 'SUPPRESS_UTM_PROJECTION must be activated when PML_CONDITIONS is set to .true.'
    endif
  else
    if (PML_INSTEAD_OF_FREE_SURFACE) &
      stop 'PML_CONDITIONS must be activated when PML_INSTEAD_OF_FREE_SURFACE is set to .true.'
  endif

  if (STACEY_INSTEAD_OF_FREE_SURFACE .and. PML_INSTEAD_OF_FREE_SURFACE) then
    print *, 'please modify Par_file and recompile solver'
    stop 'error: STACEY_INSTEAD_OF_FREE_SURFACE and PML_INSTEAD_OF_FREE_SURFACE are both set to .true.'
  endif

  ! checks the MOVIE_TYPE parameter
  if (MOVIE_TYPE < 1 .or. MOVIE_TYPE > 3) then
    stop 'error: MOVIE_TYPE must be either 1, 2 or 3! Please modify Par_file and recompile solver'
  endif

  ! check that the code has been compiled with the right values
  if (myrank == 0) then
    HEADER_FILE = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'/values_from_mesher.h'

    open(unit=IOUT,file=trim(HEADER_FILE),status='old',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(HEADER_FILE)
      print *
      print *,'please check if xgenerate_databases has been run before this solver, exiting now...'
      stop 'error opening file values_from_mesher.h'
    endif
    read(IOUT,NML=MESHER)
    close(IOUT)

    if (STACEY_INSTEAD_OF_FREE_SURFACE .neqv. ABSORB_FREE_SURFACE_VAL) then
      write(IMAIN,*) 'STACEY_INSTEAD_OF_FREE_SURFACE:',STACEY_INSTEAD_OF_FREE_SURFACE,ABSORB_FREE_SURFACE_VAL
      call exit_MPI(myrank,'it seems you have changed STACEY_INSTEAD_OF_FREE_SURFACE, you need to rerun xgenerate_databases')
    endif
  endif

  call synchronize_all()

  ! checks directories
  if (myrank == 0) then
    ! tests if OUTPUT_FILES directory exists
    ! note: inquire behaves differently when using intel ifort or gfortran compilers
    !INQUIRE( FILE = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'/.', EXIST = exists)
    open(IOUT,file=trim(OUTPUT_FILES)//'/dummy.txt',status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,"OUTPUT_FILES directory does not work: ",trim(OUTPUT_FILES)
      call exit_MPI(myrank,'error OUTPUT_FILES directory')
    endif
    close(IOUT,status='delete')

    ! tests if LOCAL_PATH directory exists
    !INQUIRE( FILE = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/.', EXIST = exists)
    open(IOUT,file=trim(LOCAL_PATH)//'/dummy.txt',status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,"LOCAL_PATH directory does not work: ",trim(LOCAL_PATH)
      call exit_MPI(myrank,'error LOCAL_PATH directory')
    endif
    close(IOUT,status='delete')
  endif

  ! safety check
  if (NB_RUNS_ACOUSTIC_GPU > 1) then
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) &
      stop 'NB_RUNS_ACOUSTIC_GPU > 1 not compatible with NUMBER_OF_SIMULTANEOUS_RUNS > 1 yet'
    if (SIMULATION_TYPE /= 1) &
      stop 'NB_RUNS_ACOUSTIC_GPU > 1 not compatible with SIMULATION_TYPE /= 1 yet'
    if (SAVE_SEISMOGRAMS_DISPLACEMENT .or. SAVE_SEISMOGRAMS_VELOCITY .or. SAVE_SEISMOGRAMS_ACCELERATION) &
      stop 'Invalid seismogram output for NB_RUNS_ACOUSTIC_GPU > 1, only pressure output implemented yet'
    if (.not. SAVE_SEISMOGRAMS_PRESSURE ) &
      stop 'NB_RUNS_ACOUSTIC_GPU > 1 not compatible with elastic wavefield seismograms yet'
    if (.not. GPU_MODE ) &
      stop 'NB_RUNS_ACOUSTIC_GPU > 1 only applies with GPU_MODE'
    if (INVERSE_FWI_FULL_PROBLEM) &
      stop 'NB_RUNS_ACOUSTIC_GPU > 1 not compatible with INVERSE_FWI_FULL_PROBLEM yet'
  endif

  ! file output
  if (SU_FORMAT .and. ASDF_FORMAT) &
    stop 'Please choose either SU_FORMAT or ASDF_FORMAT, both outputs together are not implemented yet...'

  ! ASDF for 1 output component-type only
  if (ASDF_FORMAT) then
    ! counts output types
    ioutputs = 0
    if (SAVE_SEISMOGRAMS_DISPLACEMENT) ioutputs = ioutputs + 1
    if (SAVE_SEISMOGRAMS_VELOCITY) ioutputs = ioutputs + 1
    if (SAVE_SEISMOGRAMS_ACCELERATION) ioutputs = ioutputs + 1
    if (SAVE_SEISMOGRAMS_PRESSURE) ioutputs = ioutputs + 1
    ! check
    if (ioutputs > 1 .or. ioutputs == 0) &
      stop 'Please save only a single type (disp/veloc/accel/pressure) when using ASDF_FORMAT...'
  endif

  ! acoustic kernel simulations
  if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
    ! for USE_TRICK_FOR_BETTER_PRESSURE, the potential_acoustic becomes the potential_dot_dot_acoustic:
    !  "use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
    !   use the second derivative of the source for the source time function instead of the source itself,
    !   and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
    !   this is mathematically equivalent, but numerically significantly more accurate because in the explicit
    !   Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
    !   thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
    !   is accurate at second order and thus contains significantly less numerical noise."
    ! however, for kernels expressions, we need both b_potential_acoustic and b_potential_dot_dot_acoustic as defined
    ! from the original acoustic potential definition.
    if (USE_TRICK_FOR_BETTER_PRESSURE) &
      stop 'For acoustic kernels, please set USE_TRICK_FOR_BETTER_PRESSURE to .false.'
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

  ! checks number of sources for adjoint simulations
  ! The limit below is somewhat arbitrary. For pure adjoint simulations (SIMULATION_TYPE == 2),
  ! the code outputs displacement (NT.S00001.BXX.semd,..) and strains (NT.S00001.SEE.semd,..)
  ! as well as source derivative kernels (src_frechet.00001,..) all for each point source.
  ! The naming convention for these files uses (.., i6.6,..), which limits the number of sources to 999999.
  ! If that is still too low, you can increase it further (if so, change all the occurrences of (.., i6.6,..) in the code).
  if (SIMULATION_TYPE /= 1 .and. NSOURCES > 999999) &
    call exit_MPI(myrank,'for adjoint simulations, NSOURCES <= 999999, if you need more change i6.6 in write_seismograms.f90')

  if (SIMULATION_TYPE /= 1 .and. POROELASTIC_SIMULATION) stop 'poroelastic simulations for adjoint runs not supported yet'

  ! snapshot file names: ADJOINT attenuation
  if (ATTENUATION .and. ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3)) &
    call create_name_database(prname_Q,myrank,OUTPUT_FILES)

  ! number of elements and points for adjoint arrays
  if (SIMULATION_TYPE == 3) then
    NSPEC_ADJOINT = NSPEC_AB
    NGLOB_ADJOINT = NGLOB_AB
  else
    ! dummy array size
    NSPEC_ADJOINT = 1
    NGLOB_ADJOINT =  1
  endif

  ! moho boundary
  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    NSPEC_BOUN = NSPEC_AB
  else
    NSPEC_BOUN = 1
  endif

  end subroutine initialize_simulation_adjoint

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_GPU()

! initialization for GPU cards

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer :: ncuda_devices,num_device,ncuda_devices_min,ncuda_devices_max

  ! GPU_MODE now defined in Par_file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "GPU_MODE Active."
    call flush_IMAIN()
  endif

  if (CUSTOM_REAL /= 4) stop 'GPU mode runs only with CUSTOM_REAL == 4'

  if (SAVE_MOHO_MESH) stop 'GPU mode does not support SAVE_MOHO_MESH yet'

  if (ATTENUATION) then
    if (N_SLS /= 3) stop 'GPU mode does not support N_SLS /= 3 yet'
  endif

  if (POROELASTIC_SIMULATION) stop 'poroelastic simulations on GPUs not supported yet'

  if (NPROC == 1 .and. NUMBER_OF_SIMULTANEOUS_RUNS > 1 ) then
    num_device = mygroup
  else
    num_device = myrank
  endif

  ! initializes GPU and outputs info to files for all processes
  call initialize_gpu_device(num_device,ncuda_devices)

  ! collects min/max of local devices found for statistics
  call synchronize_all()
  call min_all_i(ncuda_devices,ncuda_devices_min)
  call max_all_i(ncuda_devices,ncuda_devices_max)

  if (myrank == 0) then
    write(IMAIN,*) "GPU number of devices per node: min =",ncuda_devices_min
    write(IMAIN,*) "                                max =",ncuda_devices_max
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine initialize_GPU
