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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine initialize_simulation()

  use adios_manager_mod
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_movie

  implicit none

  include 'version.fh'

  integer :: ier
  logical :: BROADCAST_AFTER_READ

  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_rank(myrank)

  ! read the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! checks flags
  call initialize_simulation_check()

  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_solver.txt',status='unknown')

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*) '**** Specfem 3-D Solver - MPI version f90 ****'
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Version: ', git_version
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
    end select

    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (ADIOS_ENABLED) then
    call adios_setup()
  endif

  ! reads in numbers of spectral elements and points for the part of the mesh handled by this process
  call create_name_database(prname,myrank,LOCAL_PATH)

#ifdef DEBUG_COUPLED
    include "../../../add_to_initialize_simulation.F90"
#endif

! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  if (ADIOS_FOR_MESH) then
    call read_mesh_for_init_ADIOS(NSPEC_AB, NGLOB_AB)
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
  if (ier /= 0) stop 'error allocating arrays for databases'

  ! mesh node locations
  allocate(xstore(NGLOB_AB), &
           ystore(NGLOB_AB), &
           zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB), &
           ispec_is_elastic(NSPEC_AB), &
           ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'error allocating arrays for material flags'
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! initializes adjoint simulations
  call initialize_simulation_adjoint()

  ! initializes GPU cards
  if (GPU_MODE) call initialize_GPU()

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

  integer :: ier

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
  if (NSOURCES < 1) call exit_MPI(myrank,'need at least one source')

  ! check simulation type
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    call exit_mpi(myrank,'SIMULATION_TYPE can only be 1, 2, or 3')

  ! gravity only on GPU supported
  if (.not. GPU_MODE .and. GRAVITY) &
    stop 'GRAVITY only supported in GPU mode'

  ! absorbing surfaces
  if (STACEY_ABSORBING_CONDITIONS) then
    ! for arbitrary orientation of elements, which face belongs to xmin,xmax,etc... -
    ! does it makes sense to have different NGLLX,NGLLY,NGLLZ?
    ! there is a problem with absorbing boundaries for faces with different NGLLX,NGLLY,NGLLZ values
    ! just to be sure for now..
    if (NGLLX /= NGLLY .and. NGLLY /= NGLLZ) &
      stop 'STACEY_ABSORBING_CONDITIONS must have NGLLX = NGLLY = NGLLZ'
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
    else if (.not. SUPPRESS_UTM_PROJECTION) then
      print *, 'please modify Par_file and recompile solver'
      stop 'SUPPRESS_UTM_PROJECTION must be activated when PML_CONDITIONS is set to .true.'
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
  if (MOVIE_TYPE /= 1 .and. MOVIE_TYPE /= 2) then
    stop 'error: MOVIE_TYPE must be either 1 or 2! Please modify Par_file and recompile solver'
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
  integer :: ncuda_devices,ncuda_devices_min,ncuda_devices_max

  ! GPU_MODE now defined in Par_file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "GPU_MODE Active."
    call flush_IMAIN()
  endif

  ! check for GPU runs
  if (NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5) stop 'GPU mode can only be used if NGLLX == NGLLY == NGLLZ == 5'

  if (CUSTOM_REAL /= 4) stop 'GPU mode runs only with CUSTOM_REAL == 4'

  if (SAVE_MOHO_MESH) stop 'GPU mode does not support SAVE_MOHO_MESH yet'

  if (ATTENUATION) then
    if (N_SLS /= 3) stop 'GPU mode does not support N_SLS /= 3 yet'
  endif

  if (POROELASTIC_SIMULATION) stop 'poroelastic simulations on GPUs not supported yet'

  ! initializes GPU and outputs info to files for all processes
  call initialize_cuda_device(myrank,ncuda_devices)

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
