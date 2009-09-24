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
  
! sizeprocs returns number of processes started
! (should be equal to NPROC)
! myrank is the rank of each process, between 0 and sizeprocs-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

! read the parameter file
  call read_parameter_file( &
        NPROC,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT, &
        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
        ATTENUATION,USE_OLSEN_ATTENUATION,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        OCEANS,ANISOTROPY,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD)

  if (sizeprocs == 1 .and. (NPROC_XI /= 1 .or. NPROC_ETA /= 1)) then
    stop 'must have NPROC_XI = NPROC_ETA = 1 for a serial run'
  endif

! check simulation type
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
        call exit_mpi(myrank,'SIMULATION_TYPE can only be 1, 2, or 3')

! check simulation parameters
  if (SIMULATION_TYPE /= 1 .and. NSOURCES > 1000) call exit_mpi(myrank, 'for adjoint simulations, NSOURCES <= 1000')
! LQY -- note: kernel simulations with attenuation turned on has been implemented

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! check that optimized routines from Deville et al. (2002) can be used
  if(NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5) &
    stop 'optimized routines from Deville et al. (2002) such as mxm_m1_m2_5points can only be used if NGLL = 5'

! info about external mesh simulation
! nlegoff -- should be put in compute_parameters and read_parameter_file for clarity
  NPROC = sizeprocs
! chris: DT_ext_mesh & NSTE_ext_mesh were in constants.h, I suppressed it, now it is Par_file & read in 
! read_parameters_file.f90
!  DT = DT_ext_mesh
!  NSTEP = NSTEP_ext_mesh
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',action='read',form='unformatted')

  read(27) NSPEC_AB
  read(27) NGLOB_AB
  close(27)

  if( ATTENUATION ) then
    !pll
    NSPEC_ATTENUATION_AB = NSPEC_AB
  else
    ! if attenuation is off, set dummy size of arrays to one
    NSPEC_ATTENUATION_AB = 1
  endif

! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_solver.txt',status='unknown')

  if(myrank == 0) then

  write(IMAIN,*)
  write(IMAIN,*) '**********************************************'
  write(IMAIN,*) '**** Specfem 3-D Solver - MPI version f90 ****'
  write(IMAIN,*) '**********************************************'
  write(IMAIN,*)
  write(IMAIN,*)

  if(FIX_UNDERFLOW_PROBLEM) write(IMAIN,*) 'Fixing slow underflow trapping problem using small initial field'

  write(IMAIN,*)
  write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
  write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
  write(IMAIN,*)

  write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi'
  write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta'
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
  write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
  write(IMAIN,*)

  endif

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) call exit_MPI(myrank,'wrong number of MPI processes')

! check that we have at least one source
  if(NSOURCES < 1) call exit_MPI(myrank,'need at least one source')


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
  allocate(xstore(NGLOB_AB))
  allocate(ystore(NGLOB_AB))
  allocate(zstore(NGLOB_AB))
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(not_fully_in_bedrock(NSPEC_AB))
  allocate(flag_sediments(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(idoubling(NSPEC_AB))
  allocate(rmass(NGLOB_AB))
  allocate(rmass_ocean_load(NGLOB_AB))
  allocate(updated_dof_ocean_load(NGLOB_AB))
  allocate(displ(NDIM,NGLOB_AB))
  allocate(veloc(NDIM,NGLOB_AB))
  allocate(accel(NDIM,NGLOB_AB))
  allocate(iflag_attenuation_store(NGLLX,NGLLY,NGLLZ,NSPEC_AB))




  end subroutine