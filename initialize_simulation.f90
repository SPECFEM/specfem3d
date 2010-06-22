!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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
                        OCEANS,ANISOTROPY,ABSORBING_CONDITIONS, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD)

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

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

  ! anisotropy arrays size
  if( ANISOTROPY ) then
    NSPEC_ANISO = NSPEC_AB
  else
    ! if off, set dummy size
    NSPEC_ANISO = 1
  endif

  ! allocate arrays for storing the databases
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(xix(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(xiy(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(xiz(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(etax(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(etay(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(etaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(gammax(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(gammay(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  ! mesh node locations  
  allocate(xstore(NGLOB_AB))
  allocate(ystore(NGLOB_AB))
  allocate(zstore(NGLOB_AB))
  ! material properties  
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB))
  allocate(ispec_is_elastic(NSPEC_AB))
  allocate(ispec_is_poroelastic(NSPEC_AB))
    
  ! ocean mass matrix
  allocate(rmass_ocean_load(NGLOB_AB))  
  
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
  
  ! sizeprocs returns number of processes started
  ! (should be equal to NPROC)
  call world_size(sizeprocs)

  ! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) call exit_MPI(myrank,'wrong number of MPI processes')

  ! check that we have at least one source
  if(NSOURCES < 1) call exit_MPI(myrank,'need at least one source')

  ! check simulation type
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
        call exit_mpi(myrank,'SIMULATION_TYPE can only be 1, 2, or 3')

  ! check that optimized routines from Deville et al. (2002) can be used
  if( USE_DEVILLE_PRODUCTS) then
    if(NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5) &
      stop 'Deville et al. (2002) routines can only be used if NGLLX = NGLLY = NGLLZ = 5'
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

  ! strain/attenuation
  if( ATTENUATION .and. SIMULATION_TYPE == 3 ) then
    NSPEC_ATT_AND_KERNEL = NSPEC_AB
  else
    NSPEC_ATT_AND_KERNEL = 1
  endif

  ! moho boundary
  if( SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3 ) then    
    NSPEC_BOUN = NSPEC_AB
  else
    NSPEC_BOUN = 1
  endif
  
  end subroutine initialize_simulation_adjoint
  
  
  
  
