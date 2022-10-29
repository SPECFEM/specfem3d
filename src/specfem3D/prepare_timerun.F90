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


  subroutine prepare_timerun()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie

  implicit none

  ! local parameters
  double precision :: tCPU,tstart
  double precision, external :: wtime

  ! synchonizes
  call synchronize_all()

  ! get MPI starting time
  tstart = wtime()

  ! user output infos
  call prepare_timerun_user_output()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()

  ! sets up time increments
  call prepare_timerun_constants()

  ! initializes arrays
  call prepare_wavefields()

  ! initializes fault rupture arrays
  call prepare_timerun_faults()

  ! prepares attenuation arrays
  call prepare_attenuation()

  ! prepares gravity arrays
  call prepare_gravity()

  ! ZN I do not use if (USE_LDDRK) call prepare_timerun_lddrk()
  ! ZN in order to avoid the error of using unallocated arrays later on in the code,
  ! ZN since R_**_lddrk are arguments in subroutine compute_forces_viscoelastic
  call prepare_timerun_lddrk()

  ! prepares C-PML arrays
  if (PML_CONDITIONS) call prepare_timerun_pml()

  ! Stacey boundaries
  call prepare_timerun_stacey()

  ! prepares ADJOINT simulations
  call prepare_timerun_adjoint()

  ! prepares noise simulations
  call prepare_noise()

  ! prepares coupling with injection boundary
  call couple_with_injection_prepare_boundary()

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_GPU()

  ! optimizes array memory layout for better performance
  call prepare_optimized_arrays()

  ! synchronize all the processes
  call synchronize_all()

  ! elapsed time since beginning of preparation
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for preparing timerun in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) '************'
    write(IMAIN,*) ' time loop'
    write(IMAIN,*) '************'
    if (USE_LDDRK) then
      write(IMAIN,*) '              scheme:         LDDRK with',NSTAGE_TIME_SCHEME,'stages'
    else
      write(IMAIN,*) '              scheme:         Newmark'
    endif
    write(IMAIN,*)
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(NSTEP*DT),' seconds'
    write(IMAIN,*) 'start time:',sngl(-t0),' seconds'
    write(IMAIN,*)

    ! flushes file buffer for main output file (IMAIN)
    call flush_IMAIN()

    !daniel debug: total time estimation
    !  average time per element per time step:
    !     elastic elements    ~ dt = 1.17e-05 s (intel xeon 2.6GHz, stand 2013)
    !                              = 3.18e-07 s (Kepler K20x, stand 2013)
    !
    !  total time per time step:
    !     T_total = dt * nspec_total
    !
    !  total time using nproc processes (slices) for NSTEP time steps:
    !     T_simulation = T_total * NSTEP / nproc

  endif

  ! synchronize all the processes
  call synchronize_all()

  end subroutine prepare_timerun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_user_output()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie

  implicit none

  ! user info
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Simulation setup:'
    write(IMAIN,*)
    if (NOISE_TOMOGRAPHY /= 0) then
      write(IMAIN,*) 'noise simulation:'
      write(IMAIN,*) '  simulation type = ',NOISE_TOMOGRAPHY
      write(IMAIN,*)
    endif

    if (ACOUSTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating acoustic simulation'
    else
      write(IMAIN,*) '  no acoustic simulation'
    endif

    if (ELASTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating elastic simulation'
    else
      write(IMAIN,*) '  no elastic simulation'
    endif

    if (POROELASTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating poroelastic simulation'
    else
      write(IMAIN,*) '  no poroelastic simulation'
    endif
    write(IMAIN,*)

    if (ATTENUATION) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
      if (UNDO_ATTENUATION_AND_OR_PML) &
        write(IMAIN,*) 'using undo_attenuation scheme'

      if (USE_OLSEN_ATTENUATION) then
        write(IMAIN,*) 'using attenuation from Olsen et al.'
      else
        write(IMAIN,*) '  not using attenuation from Olsen et al.'
      endif
    else
      write(IMAIN,*) '  no attenuation'
    endif

    if (ANISOTROPY) then
      write(IMAIN,*) 'incorporating anisotropy'
    else
      write(IMAIN,*) '  no anisotropy'
    endif

    if (APPROXIMATE_OCEAN_LOAD) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) '  no oceans'
    endif

    if (GRAVITY) then
      write(IMAIN,*) 'incorporating gravity'
    else
      write(IMAIN,*) '  no gravity'
    endif

    if (MOVIE_SIMULATION) then
      write(IMAIN,*) 'incorporating movie simulation'
    else
      write(IMAIN,*) '  no movie simulation'
    endif

    write(IMAIN,*)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! safety checks
  if (ACOUSTIC_SIMULATION) then
    if (USE_TRICK_FOR_BETTER_PRESSURE .and. SIMULATION_TYPE == 3) &
      stop 'for SIMULATION_TYPE == 3, acoustic kernels need to have USE_TRICK_FOR_BETTER_PRESSURE set to .false.'
  endif
  if (UNDO_ATTENUATION_AND_OR_PML) then
    if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION .and. SIMULATION_TYPE == 3) &
      stop 'for SIMULATION_TYPE == 3, UNDO_ATT for coupled elastic/acoustic simulations not implemented yet'
  endif

  end subroutine prepare_timerun_user_output
!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_mass_matrices()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing mass matrices"
    call flush_IMAIN()
  endif

  ! synchronize all the processes before assembling the mass matrix
  ! to make sure all the nodes have finished to read their databases
  call synchronize_all()

  ! the mass matrix needs to be assembled with MPI here once and for all
  if (ACOUSTIC_SIMULATION) then
    ! adds contributions
    if (STACEY_ABSORBING_CONDITIONS) then
      if (USE_LDDRK) then
        rmass_acoustic(:) = rmass_acoustic(:)
      else
        ! adds boundary contributions for newmark scheme
        rmass_acoustic(:) = rmass_acoustic(:) + rmassz_acoustic(:)
      endif
      ! not needed anymore
      deallocate(rmassz_acoustic)
    endif

    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_acoustic, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmass_acoustic <= 0._CUSTOM_REAL) rmass_acoustic = 1._CUSTOM_REAL
    rmass_acoustic(:) = 1._CUSTOM_REAL / rmass_acoustic(:)

  endif

  if (ELASTIC_SIMULATION) then
    ! switches to three-component mass matrix
    if (STACEY_ABSORBING_CONDITIONS) then
      if (USE_LDDRK) then
        rmassx(:) = rmass(:)
        rmassy(:) = rmass(:)
        rmassz(:) = rmass(:)
      else
        ! adds boundary contributions for Newmark scheme
        rmassx(:) = rmass(:) + rmassx(:)
        rmassy(:) = rmass(:) + rmassy(:)
        rmassz(:) = rmass(:) + rmassz(:)
      endif
    else
      rmassx(:) = rmass(:)
      rmassy(:) = rmass(:)
      rmassz(:) = rmass(:)
    endif
    ! not needed anymore
    deallocate(rmass)

    ! assemble mass matrix
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassx, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassy, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassz, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmassx <= 0._CUSTOM_REAL) rmassx = 1._CUSTOM_REAL
    where(rmassy <= 0._CUSTOM_REAL) rmassy = 1._CUSTOM_REAL
    where(rmassz <= 0._CUSTOM_REAL) rmassz = 1._CUSTOM_REAL
    rmassx(:) = 1._CUSTOM_REAL / rmassx(:)
    rmassy(:) = 1._CUSTOM_REAL / rmassy(:)
    rmassz(:) = 1._CUSTOM_REAL / rmassz(:)

    ! ocean load
    if (APPROXIMATE_OCEAN_LOAD) then
      call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_ocean_load, &
                                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                        my_neighbors_ext_mesh)
      where(rmass_ocean_load <= 0._CUSTOM_REAL) rmass_ocean_load = 1._CUSTOM_REAL
      rmass_ocean_load(:) = 1._CUSTOM_REAL / rmass_ocean_load(:)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_solid_poroelastic, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_fluid_poroelastic, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    ! fills mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmass_solid_poroelastic <= 0._CUSTOM_REAL) rmass_solid_poroelastic = 1._CUSTOM_REAL
    where(rmass_fluid_poroelastic <= 0._CUSTOM_REAL) rmass_fluid_poroelastic = 1._CUSTOM_REAL
    rmass_solid_poroelastic(:) = 1._CUSTOM_REAL / rmass_solid_poroelastic(:)
    rmass_fluid_poroelastic(:) = 1._CUSTOM_REAL / rmass_fluid_poroelastic(:)
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!


  subroutine prepare_timerun_constants()

! precomputes constants for time integration

  use specfem_par

  implicit none

  if (myrank == 0) then
    write(IMAIN,*) "preparing constants"
    call flush_IMAIN()
  endif

  ! time scheme
  if (.not. USE_LDDRK) then
    ! Newmark time scheme, only single update needed
    NSTAGE_TIME_SCHEME = 1
  else
    ! LDDRK time scheme with 6-stages
    NSTAGE_TIME_SCHEME = 6
  endif

  ! distinguish between single and double precision for reals
  deltat = real(DT,kind=CUSTOM_REAL)
  deltatover2 = deltat/2._CUSTOM_REAL
  deltatsqover2 = deltat*deltat/2._CUSTOM_REAL

  ! adjoint runs: timing
  if (SIMULATION_TYPE == 3) then
    if (UNDO_ATTENUATION_AND_OR_PML) then
      ! moves forward
      b_deltat = deltat
      b_deltatover2 = deltatover2
      b_deltatsqover2 = deltatsqover2
    else
      ! reconstructed wavefield moves backward in time from last snapshot
      ! backward/reconstructed wavefields: time stepping is in time-reversed sense
      ! (negative time increments)
      b_deltat = - real(DT,kind=CUSTOM_REAL)
      b_deltatover2 = b_deltat/2._CUSTOM_REAL
      b_deltatsqover2 = b_deltat*b_deltat/2._CUSTOM_REAL
    endif
  else
    ! will not be used, but initialized
    b_deltat = 0._CUSTOM_REAL
    b_deltatover2 = 0._CUSTOM_REAL
    b_deltatsqover2 = 0._CUSTOM_REAL
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_constants

!
!-------------------------------------------------------------------------------------------------
!


  subroutine prepare_timerun_lddrk()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  integer :: ier

  if (USE_LDDRK) then
    NGLOB_AB_LDDRK = NGLOB_AB
    NSPEC_ATTENUATION_AB_LDDRK = NSPEC_ATTENUATION_AB
  else
    NGLOB_AB_LDDRK = 1
    NSPEC_ATTENUATION_AB_LDDRK = 1
  endif

  ! user output
  if (USE_LDDRK) then
    if (myrank == 0) then
      write(IMAIN,*) "preparing LDDRK"
      call flush_IMAIN()
    endif
  endif

  if (ACOUSTIC_SIMULATION) then
    allocate(potential_acoustic_lddrk(NGLOB_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2120')
    if (ier /= 0) stop 'Error allocating array potential_acoustic_lddrk'
    allocate(potential_dot_acoustic_lddrk(NGLOB_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2121')
    if (ier /= 0) stop 'Error allocating array potential_dot_acoustic_lddrk'

    potential_acoustic_lddrk(:) = 0._CUSTOM_REAL
    potential_dot_acoustic_lddrk(:) = 0._CUSTOM_REAL
    if (FIX_UNDERFLOW_PROBLEM) then
      potential_acoustic_lddrk(:) = VERYSMALLVAL
      potential_dot_acoustic_lddrk(:) = VERYSMALLVAL
    endif
  endif

  if (ELASTIC_SIMULATION) then
    allocate(displ_lddrk(NDIM,NGLOB_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2122')
    if (ier /= 0) stop 'Error allocating array displ_lddrk'
    allocate(veloc_lddrk(NDIM,NGLOB_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2123')
    if (ier /= 0) stop 'Error allocating array veloc_lddrk'

    displ_lddrk(:,:) = 0._CUSTOM_REAL
    veloc_lddrk(:,:) = 0._CUSTOM_REAL
    if (FIX_UNDERFLOW_PROBLEM) then
      displ_lddrk(:,:) = VERYSMALLVAL
      veloc_lddrk(:,:) = VERYSMALLVAL
    endif

    ! note: currently, they need to be defined, as they are used in some subroutine arguments
    allocate(R_xx_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2124')
    allocate(R_yy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2125')
    allocate(R_xy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2126')
    allocate(R_xz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2127')
    allocate(R_yz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2128')
    if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

    allocate(R_trace_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2129')
    if (ier /= 0) stop 'Error allocating array R_trace_lddrk etc.'

    if (SIMULATION_TYPE == 3) then
      allocate(b_R_xx_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2130')
      allocate(b_R_yy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2131')
      allocate(b_R_xy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2132')
      allocate(b_R_xz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2133')
      allocate(b_R_yz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2134')
      if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

      allocate(b_R_trace_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2135')
      if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'
    endif

    ! initializes
    if (ATTENUATION) then
      R_xx_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
      R_yy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
      R_xy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
      R_xz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
      R_yz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
      R_trace_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
      if (FIX_UNDERFLOW_PROBLEM) then
        R_xx_lddrk(:,:,:,:,:) = VERYSMALLVAL
        R_yy_lddrk(:,:,:,:,:) = VERYSMALLVAL
        R_xy_lddrk(:,:,:,:,:) = VERYSMALLVAL
        R_xz_lddrk(:,:,:,:,:) = VERYSMALLVAL
        R_yz_lddrk(:,:,:,:,:) = VERYSMALLVAL
        R_trace_lddrk(:,:,:,:,:) = VERYSMALLVAL
      endif

      if (SIMULATION_TYPE == 3) then
        b_R_xx_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_yy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_xy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_xz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_yz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_trace_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        if (FIX_UNDERFLOW_PROBLEM) then
          b_R_xx_lddrk(:,:,:,:,:) = VERYSMALLVAL
          b_R_yy_lddrk(:,:,:,:,:) = VERYSMALLVAL
          b_R_xy_lddrk(:,:,:,:,:) = VERYSMALLVAL
          b_R_xz_lddrk(:,:,:,:,:) = VERYSMALLVAL
          b_R_yz_lddrk(:,:,:,:,:) = VERYSMALLVAL
          b_R_trace_lddrk(:,:,:,:,:) = VERYSMALLVAL
        endif
      endif
    endif
  endif

  ! safety stop
  if (POROELASTIC_SIMULATION) then
    if (USE_LDDRK) &
      stop 'LDDRK has not been implemented for POROELASTIC_SIMULATION'
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_pml()

  use constants, only: IMAIN,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: myrank,deltat,SIMULATION_TYPE,GPU_MODE, &
    UNDO_ATTENUATION_AND_OR_PML,PML_CONDITIONS,SAVE_MESH_FILES

  use pml_par

  implicit none

  ! local parameters
  integer :: ispec,ispec_CPML,NSPEC_CPML_GLOBAL
  integer :: i,j,k,ier
  real(kind=CUSTOM_REAL) :: deltatpow2,deltatpow3,deltatpow4,deltat_half
  real(kind=CUSTOM_REAL) :: coef0_1,coef1_1,coef2_1, &
                            coef0_2,coef1_2,coef2_2, &
                            coef0_3,coef1_3,coef2_3
  real(kind=CUSTOM_REAL) :: kappa_x,d_x,alpha_x, &
                            kappa_y,d_y,alpha_y, &
                            kappa_z,d_z,alpha_z
  real(kind=CUSTOM_REAL) :: beta_x,beta_y,beta_z

  ! checks if anything to do
  if (.not. PML_CONDITIONS) return

  if (myrank == 0) then
    write(IMAIN,*) "preparing CPML"
    call flush_IMAIN()
  endif

  ! safety stops
  if (SIMULATION_TYPE /= 1 .and. .not. UNDO_ATTENUATION_AND_OR_PML) &
    stop 'Error: PMLs for adjoint runs require the flag UNDO_ATTENUATION_AND_OR_PML to be set'
  if (GPU_MODE) &
    stop 'Error: PMLs only supported in CPU mode for now'

  ! total number of PML elements
  call sum_all_i(NSPEC_CPML,NSPEC_CPML_GLOBAL)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  incorporating C-PML  '
    write(IMAIN,*)
    write(IMAIN,*) '  number of C-PML spectral elements in the global mesh: ', NSPEC_CPML_GLOBAL
    write(IMAIN,*)
    write(IMAIN,*) '  thickness of C-PML layer in X direction: ', CPML_width_x
    write(IMAIN,*) '  thickness of C-PML layer in Y direction: ', CPML_width_y
    write(IMAIN,*) '  thickness of C-PML layer in Z direction: ', CPML_width_z
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

!! DK DK Feb 2017: I do not think there is any particular reason for that any more, thus commenting it out
! ! checks that 8-node mesh elements are used (27-node elements are not supported)
! if (NGNOD /= NGNOD_EIGHT_CORNERS) &
!   stop 'error: the C-PML code works for 8-node bricks only; should be made more general'

  ! allocates and initializes C-PML arrays
  call pml_allocate_arrays()

  ! defines C-PML spectral elements local indexing
  do ispec_CPML = 1,NSPEC_CPML
    ispec = CPML_to_spec(ispec_CPML)
    spec_to_CPML(ispec) = ispec_CPML
  enddo

  ! defines C-PML element type array: 1 = face, 2 = edge, 3 = corner
  do ispec_CPML = 1,NSPEC_CPML
    ! X_surface C-PML
    if (CPML_regions(ispec_CPML) == 1) then
      CPML_type(ispec_CPML) = 1

    ! Y_surface C-PML
    else if (CPML_regions(ispec_CPML) == 2) then
      CPML_type(ispec_CPML) = 1

    ! Z_surface C-PML
    else if (CPML_regions(ispec_CPML) == 3) then
      CPML_type(ispec_CPML) = 1

    ! XY_edge C-PML
    else if (CPML_regions(ispec_CPML) == 4) then
      CPML_type(ispec_CPML) = 2

    ! XZ_edge C-PML
    else if (CPML_regions(ispec_CPML) == 5) then
      CPML_type(ispec_CPML) = 2

    ! YZ_edge C-PML
    else if (CPML_regions(ispec_CPML) == 6) then
      CPML_type(ispec_CPML) = 2

    ! XYZ_corner C-PML
    else if (CPML_regions(ispec_CPML) == 7) then
      CPML_type(ispec_CPML) = 3
    endif
  enddo

  ! useful constants
  deltatpow2 = deltat**2
  deltatpow3 = deltat**3
  deltatpow4 = deltat**4
  deltat_half = deltat * 0.5_CUSTOM_REAL

  ! arrays for coefficients
  allocate(convolution_coef_acoustic_alpha(9,NGLLX,NGLLY,NGLLZ,NSPEC_CPML), &
           convolution_coef_acoustic_beta(9,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating coef_acoustic array')

  ! initializes
  convolution_coef_acoustic_alpha(:,:,:,:,:) = 0._CUSTOM_REAL
  convolution_coef_acoustic_beta(:,:,:,:,:) = 0._CUSTOM_REAL

  ! pre-computes convolution coefficients
  do ispec_CPML = 1,NSPEC_CPML
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          kappa_x = k_store_x(i,j,k,ispec_CPML)
          kappa_y = k_store_y(i,j,k,ispec_CPML)
          kappa_z = k_store_z(i,j,k,ispec_CPML)
          d_x = d_store_x(i,j,k,ispec_CPML)
          d_y = d_store_y(i,j,k,ispec_CPML)
          d_z = d_store_z(i,j,k,ispec_CPML)
          alpha_x = alpha_store_x(i,j,k,ispec_CPML)
          alpha_y = alpha_store_y(i,j,k,ispec_CPML)
          alpha_z = alpha_store_z(i,j,k,ispec_CPML)

          ! gets recursive convolution coefficients
          ! alpha coefficients
          call compute_convolution_coef(alpha_x, coef0_1, coef1_1, coef2_1)
          call compute_convolution_coef(alpha_y, coef0_2, coef1_2, coef2_2)
          call compute_convolution_coef(alpha_z, coef0_3, coef1_3, coef2_3)

          convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML) = coef0_1
          convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML) = coef1_1
          convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML) = coef2_1

          convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML) = coef0_2
          convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML) = coef1_2
          convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML) = coef2_2

          convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML) = coef0_3
          convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML) = coef1_3
          convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML) = coef2_3

          ! gets recursive convolution coefficients
          ! beta coefficients
          beta_x = alpha_x + d_x / kappa_x
          beta_y = alpha_y + d_y / kappa_y
          beta_z = alpha_z + d_z / kappa_z

          call compute_convolution_coef(beta_x, coef0_1, coef1_1, coef2_1)
          call compute_convolution_coef(beta_y, coef0_2, coef1_2, coef2_2)
          call compute_convolution_coef(beta_z, coef0_3, coef1_3, coef2_3)

          convolution_coef_acoustic_beta(1,i,j,k,ispec_CPML) = coef0_1
          convolution_coef_acoustic_beta(2,i,j,k,ispec_CPML) = coef1_1
          convolution_coef_acoustic_beta(3,i,j,k,ispec_CPML) = coef2_1

          convolution_coef_acoustic_beta(4,i,j,k,ispec_CPML) = coef0_2
          convolution_coef_acoustic_beta(5,i,j,k,ispec_CPML) = coef1_2
          convolution_coef_acoustic_beta(6,i,j,k,ispec_CPML) = coef2_2

          convolution_coef_acoustic_beta(7,i,j,k,ispec_CPML) = coef0_3
          convolution_coef_acoustic_beta(8,i,j,k,ispec_CPML) = coef1_3
          convolution_coef_acoustic_beta(9,i,j,k,ispec_CPML) = coef2_3
        enddo
      enddo
    enddo
  enddo

  ! outputs informations about C-PML elements in VTK-file format
  if (SAVE_MESH_FILES) call pml_output_VTKs()

  ! synchonizes
  call synchronize_all()

  contains

    subroutine compute_convolution_coef(bb,coef0,coef1,coef2)

    use constants, only: CUSTOM_REAL
    use specfem_par, only: deltat
    use pml_par, only: FIRST_ORDER_CONVOLUTION,min_distance_between_CPML_parameter

    implicit none

    real(kind=CUSTOM_REAL),intent(in) :: bb
    real(kind=CUSTOM_REAL),intent(out) :: coef0, coef1, coef2

    real(kind=CUSTOM_REAL) :: bbpow2,bbpow3 !,c0
    real(kind=CUSTOM_REAL) :: temp

    ! permanent factors (avoids divisions which are computationally expensive)
    ! note: compilers precompute these constant factors (thus division in these statemenets are still fine)
    real(kind=CUSTOM_REAL),parameter :: ONE_OVER_8 = 0.125_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: ONE_OVER_48 = 1._CUSTOM_REAL / 48._CUSTOM_REAL
    !real(kind=CUSTOM_REAL),parameter :: ONE_OVER_128 = 0.0078125_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: ONE_OVER_384 = 1._CUSTOM_REAL / 384._CUSTOM_REAL

    real(kind=CUSTOM_REAL),parameter :: FACTOR_A = 3._CUSTOM_REAL / 8._CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: FACTOR_B = 7._CUSTOM_REAL / 48._CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: FACTOR_C = 5._CUSTOM_REAL / 128._CUSTOM_REAL

    ! unused
    !real(kind=CUSTOM_REAL),parameter :: ONE_OVER_12 = 1._CUSTOM_REAL / 12._CUSTOM_REAL
    !real(kind=CUSTOM_REAL),parameter :: ONE_OVER_24 = 1._CUSTOM_REAL / 24._CUSTOM_REAL
    !real(kind=CUSTOM_REAL),parameter :: ONE_OVER_960 = 1._CUSTOM_REAL / 960._CUSTOM_REAL
    !real(kind=CUSTOM_REAL),parameter :: ONE_OVER_1920 = 1._CUSTOM_REAL / 1920._CUSTOM_REAL
    !real(kind=CUSTOM_REAL),parameter :: SEVEN_OVER_3840 = 7._CUSTOM_REAL / 3840._CUSTOM_REAL
    !real(kind=CUSTOM_REAL),parameter :: FIVE_OVER_11520 = 5._CUSTOM_REAL/11520._CUSTOM_REAL

    ! recursive convolution coefficients
    !
    ! see Xie et al. (2014), second-order recursive scheme given by eq. (60)
    !                        and also appendix D, eq. (D6a) and (D6b) for p = 0
    !
    ! coefficients needed for the recursive scheme are:
    !    coef0 = exp(-b delta_t)
    !          = exp(- 1/2 b delta_t) * exp(- 1/2 b delta_t)
    !
    !    coef1 = 1/b (1 - exp( - 1/2 b delta_t)                             -> see also factor xi_0^(n+1) in eq. D6b
    !
    !    coef2 = 1/b (1 - exp( - 1/2 b delta_t) exp(- 1/2 b delta_t)
    !          = coef1 * exp(- 1/2 b delta_t)                               -> see also factor xi_0^n in eq. D6a
    !
    ! helper variables
    temp = exp(- 0.5_CUSTOM_REAL * bb * deltat)

    ! calculates coefficients
    !
    ! note: exponentials are expensive functions
    !c0 = exp(-a)
    !
    ! cheap exponential (up to 5 terms exp(x) = 1 + x *( 1 + x/2 * (1 + x/3 * (1 + x/4 * (1 + x/5 * ..))))
    !x0 = -a
    !c0 = 1.0 + x0 * (1.0 + 0.5 * x0 * (1.0  + 0.333333333333_CUSTOM_REAL * x0 * (1.0 + 0.25 * x0 * (1.0 + 0.2 * x0))))

    !  real function my_exp(n, x) result(f)
    !  integer, intent(in) :: n, x
    !  f = 1.0
    !  do i = n-1,1,-1
    !      f = 1.0 + x * f / i
    !  enddo
    !  end function

    ! determines coefficients
    coef0 = temp*temp

    if (abs(bb) >= min_distance_between_CPML_parameter) then
      if (FIRST_ORDER_CONVOLUTION) then
        ! first-order scheme
        coef1 = (1._CUSTOM_REAL - coef0) / bb
        coef2 = 0._CUSTOM_REAL
      else
        ! second-order convolution scheme
        ! calculates coefficients
        !
        ! cheap exponential (up to 5 terms exp(x) = 1 + x *( 1 + x/2 * (1 + x/3 * (1 + x/4 * (1 + x/5 * ..))))
        !x1 = - 0.5 * bb * deltat
        !temp = 1.0 + x1 * (1.0 + 0.5 * x1 * (1.0  + 0.333333333333_CUSTOM_REAL * x1 * (1.0 + 0.25 * x1 * (1.0 + 0.2 * x1))))

        coef1 = (1._CUSTOM_REAL - temp) / bb
        coef2 = coef1 * temp
      endif
    else
      ! approximation for small beta values
      if (FIRST_ORDER_CONVOLUTION) then
        coef1 = deltat
        coef2 = 0._CUSTOM_REAL
      else
        ! Taylor expansion to third-order
        ! (writing out powers can help for speed, that is bb * bb is usually faster than bb**2)
        bbpow2 = bb * bb
        bbpow3 = bbpow2 * bb

        coef1 = deltat_half + &
                (- ONE_OVER_8 * deltatpow2 * bb + ONE_OVER_48 * deltatpow3 * bbpow2 - ONE_OVER_384 * deltatpow4 * bbpow3)
        coef2 = deltat_half + &
                (- FACTOR_A * deltatpow2 * bb + FACTOR_B * deltatpow3 * bbpow2 - FACTOR_C * deltatpow4 * bbpow3)
      endif
    endif

    end subroutine compute_convolution_coef

  end subroutine prepare_timerun_pml

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_adjoint()

! prepares adjoint simulations

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer :: ier
  integer :: ispec,ispec2D

  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    if (myrank == 0) then
      write(IMAIN,*) "preparing adjoint fields"
      call flush_IMAIN()
    endif
  endif

  ! moment tensor derivatives
  if (nrec_local > 0 .and. SIMULATION_TYPE == 2) then
    ! allocate Frechet derivatives array
    allocate(Mxx_der(nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2136')
    allocate(Myy_der(nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2137')
    allocate(Mzz_der(nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2138')
    allocate(Mxy_der(nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2139')
    allocate(Mxz_der(nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2140')
    allocate(Myz_der(nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2141')
    allocate(sloc_der(NDIM,nrec_local),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2142')
    if (ier /= 0) stop 'error allocating array Mxx_der and following arrays'
    Mxx_der = 0._CUSTOM_REAL
    Myy_der = 0._CUSTOM_REAL
    Mzz_der = 0._CUSTOM_REAL
    Mxy_der = 0._CUSTOM_REAL
    Mxz_der = 0._CUSTOM_REAL
    Myz_der = 0._CUSTOM_REAL
    sloc_der = 0._CUSTOM_REAL
  endif

  ! initializes adjoint kernels and reconstructed/backward wavefields
  if (SIMULATION_TYPE == 3) then
    ! elastic domain
    if (ELASTIC_SIMULATION) then
      rho_kl(:,:,:,:)   = 0._CUSTOM_REAL

      if (ANISOTROPIC_KL) then
        cijkl_kl(:,:,:,:,:) = 0._CUSTOM_REAL
      else
        mu_kl(:,:,:,:)    = 0._CUSTOM_REAL
        kappa_kl(:,:,:,:) = 0._CUSTOM_REAL
      endif

      if (APPROXIMATE_HESS_KL) then
        hess_kl(:,:,:,:)   = 0._CUSTOM_REAL
        hess_rho_kl(:,:,:,:)   = 0._CUSTOM_REAL
        hess_mu_kl(:,:,:,:)   = 0._CUSTOM_REAL
        hess_kappa_kl(:,:,:,:)   = 0._CUSTOM_REAL
      endif

      ! reconstructed/backward elastic wavefields
      b_displ = 0._CUSTOM_REAL
      b_veloc = 0._CUSTOM_REAL
      b_accel = 0._CUSTOM_REAL
      if (FIX_UNDERFLOW_PROBLEM) b_displ = VERYSMALLVAL

      ! moho kernels
      if (SAVE_MOHO_MESH) moho_kl(:,:) = 0._CUSTOM_REAL
    endif

    ! acoustic domain
    if (ACOUSTIC_SIMULATION) then
      rho_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
      kappa_ac_kl(:,:,:,:) = 0._CUSTOM_REAL

      if (APPROXIMATE_HESS_KL) then
         hess_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
         hess_rho_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
         hess_kappa_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
      endif

      ! reconstructed/backward acoustic potentials
      b_potential_acoustic = 0._CUSTOM_REAL
      b_potential_dot_acoustic = 0._CUSTOM_REAL
      b_potential_dot_dot_acoustic = 0._CUSTOM_REAL
      if (FIX_UNDERFLOW_PROBLEM) b_potential_acoustic = VERYSMALLVAL
    endif

    ! poroelastic domain
    if (POROELASTIC_SIMULATION) then
      rhot_kl(:,:,:,:)   = 0._CUSTOM_REAL
      rhof_kl(:,:,:,:)   = 0._CUSTOM_REAL
      sm_kl(:,:,:,:)   = 0._CUSTOM_REAL
      eta_kl(:,:,:,:)   = 0._CUSTOM_REAL
      mufr_kl(:,:,:,:)    = 0._CUSTOM_REAL
      B_kl(:,:,:,:) = 0._CUSTOM_REAL
      C_kl(:,:,:,:) = 0._CUSTOM_REAL
      M_kl(:,:,:,:) = 0._CUSTOM_REAL

      !if (APPROXIMATE_HESS_KL) &
      !  hess_kl(:,:,:,:)   = 0._CUSTOM_REAL

      ! reconstructed/backward elastic wavefields
      b_displs_poroelastic = 0._CUSTOM_REAL
      b_velocs_poroelastic = 0._CUSTOM_REAL
      b_accels_poroelastic = 0._CUSTOM_REAL
      b_displw_poroelastic = 0._CUSTOM_REAL
      b_velocw_poroelastic = 0._CUSTOM_REAL
      b_accelw_poroelastic = 0._CUSTOM_REAL
      if (FIX_UNDERFLOW_PROBLEM) b_displs_poroelastic = VERYSMALLVAL
      if (FIX_UNDERFLOW_PROBLEM) b_displw_poroelastic = VERYSMALLVAL
    endif
  endif

! initialize Moho boundary index
  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    ! note: boundary arrays are setup in mesher, here we only need to assign the boundary element index (ispec2D)
    !       to each corresponding element (ispec) for storing strain values at the right array place
    allocate(ispec2D_moho_top(NSPEC_AB),ispec2D_moho_bot(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating ispec2D_moho arrays')
    ispec2D_moho_top(:) = 0
    ispec2D_moho_bot(:) = 0

    ! sets indexing from ispec to ispec2D
    do ispec2D = 1,NSPEC2D_MOHO
      ! bottom element
      ispec = ibelm_moho_bot(ispec2D)
      if (ispec < 1 .or. ispec > NSPEC_AB) stop 'Invalid index in ibelm_moho_bot'
      ! indexing from ispec to ispec2D
      ispec2D_moho_bot(ispec) = ispec2D

      ! top element
      ispec = ibelm_moho_top(ispec2D)
      if (ispec < 1 .or. ispec > NSPEC_AB) stop 'Invalid index in ibelm_moho_bot'
      ! indexing from ispec to ispec2D
      ispec2D_moho_top(ispec) = ispec2D
    enddo
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_adjoint

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_stacey()

! prepares stacey absorbing boundary

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_coupling

  implicit none

  ! local parameters
  integer :: ier
  integer(kind=8) :: filesize

  ! initializes
  SAVE_STACEY = .false.    ! flag to indicate if we need to read/write boundary contributions from disk

  ! stacey absorbing fields will be reconstructed for adjoint simulations
  ! using snapshot files of wavefields
  if (STACEY_ABSORBING_CONDITIONS) then

    if (myrank == 0) then
      write(IMAIN,*) "preparing Stacey absorbing boundaries"
      call flush_IMAIN()
    endif

    ! sets flag to check if we need to save the stacey contributions to file
    if (UNDO_ATTENUATION_AND_OR_PML) then
      ! not needed for undo_attenuation scheme
      SAVE_STACEY = .false.
    else
      ! save in simulation type 1 with save_forward set, read back in simulation type 3
      if (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) then
        SAVE_STACEY = .true.
      else
        SAVE_STACEY = .false.
      endif
    endif

    ! opens absorbing wavefield saved/to-be-saved by forward simulations
    if (num_abs_boundary_faces > 0 .and. SAVE_STACEY) then
      b_num_abs_boundary_faces = num_abs_boundary_faces
    else
      b_num_abs_boundary_faces = 0
    endif

    ! elastic domains
    if (ELASTIC_SIMULATION) then
      ! allocates wavefield
      if (b_num_abs_boundary_faces > 0) then
        allocate(b_absorb_field(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2143')
        if (ier /= 0) stop 'error allocating array b_absorb_field'
      else
        ! dummy
        allocate(b_absorb_field(1,1,1))
      endif
      b_absorb_field(:,:,:) = 0.0_CUSTOM_REAL

      if (num_abs_boundary_faces > 0 .and. SAVE_STACEY) then
        ! size of single record
        b_reclen_field = CUSTOM_REAL * NDIM * NGLLSQUARE * num_abs_boundary_faces

        ! check integer size limit: size of b_reclen_field must fit onto an 4-byte integer
        if (num_abs_boundary_faces > int(2147483646.0 / (CUSTOM_REAL * NDIM * NGLLSQUARE))) then
          print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_field
          print *,'  ',CUSTOM_REAL, NDIM, NGLLSQUARE, num_abs_boundary_faces
          print *,'bit size Fortran: ',bit_size(b_reclen_field)
          call exit_MPI(myrank,"error b_reclen_field integer limit")
        endif

        ! total file size
        filesize = b_reclen_field
        filesize = filesize * NSTEP

        if (SIMULATION_TYPE == 3) then
          ! opens existing files
          call open_file_abs_r(IOABS,trim(prname)//'absorb_field.bin', &
                              len_trim(trim(prname)//'absorb_field.bin'), &
                              filesize)
        else
          ! opens new file
          call open_file_abs_w(IOABS,trim(prname)//'absorb_field.bin', &
                              len_trim(trim(prname)//'absorb_field.bin'), &
                              filesize)
        endif
      endif

      if (COUPLE_WITH_INJECTION_TECHNIQUE) then
        ! boundary contribution to save together with absorbing boundary array b_absorb_field
        ! for reconstructing backward wavefields in kernel simulations
        if (b_num_abs_boundary_faces > 0 .and. SIMULATION_TYPE == 1) then
          ! only needed for forward simulation to store boundary contributions
          allocate(b_boundary_injection_field(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 2198')
          if (ier /= 0) stop 'error allocating array b_boundary_injection_field'
        else
          ! dummy
          allocate(b_boundary_injection_field(1,1,1))
        endif
        b_boundary_injection_field(:,:,:) = 0.0_CUSTOM_REAL
      endif

    endif

    ! acoustic domains
    if (ACOUSTIC_SIMULATION) then
      ! allocates wavefield
      if (b_num_abs_boundary_faces > 0) then
        allocate(b_absorb_potential(NGLLSQUARE,b_num_abs_boundary_faces * NB_RUNS_ACOUSTIC_GPU),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2144')
        if (ier /= 0) stop 'error allocating array b_absorb_potential'
      else
        ! dummy
        allocate(b_absorb_potential(1,1))
      endif
      b_absorb_potential(:,:) = 0.0_CUSTOM_REAL

      if (num_abs_boundary_faces > 0 .and. SAVE_STACEY) then
        ! size of single record
        b_reclen_potential = CUSTOM_REAL * NGLLSQUARE * num_abs_boundary_faces * NB_RUNS_ACOUSTIC_GPU


        ! check integer size limit: size of b_reclen_potential must fit onto an 4-byte integer
        if (num_abs_boundary_faces > int(2147483646.0 / (CUSTOM_REAL * NGLLSQUARE))) then
          print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_potential
          print *,'  ',CUSTOM_REAL, NGLLSQUARE, num_abs_boundary_faces
          print *,'bit size Fortran: ',bit_size(b_reclen_potential)
          call exit_MPI(myrank,"error b_reclen_potential integer limit")
        endif

        ! total file size (two lines to implicitly convert to 8-byte integers)
        filesize = b_reclen_potential
        filesize = filesize*NSTEP

        ! debug check size limit
        !if (NSTEP > 2147483646 / b_reclen_potential) then
        !  print *,'file size needed exceeds integer 4-byte limit: ',b_reclen_potential,NSTEP
        !  print *,'  ',CUSTOM_REAL, NGLLSQUARE, num_abs_boundary_faces,NSTEP
        !  print *,'file size Fortran: ',filesize
        !  print *,'file bit size Fortran: ',bit_size(filesize)
        !endif

        if (SIMULATION_TYPE == 3) then
          ! opens existing files
          call open_file_abs_r(IOABS_AC,trim(prname)//'absorb_potential.bin', &
                              len_trim(trim(prname)//'absorb_potential.bin'), &
                              filesize)
        else
          ! opens new file
          call open_file_abs_w(IOABS_AC,trim(prname)//'absorb_potential.bin', &
                              len_trim(trim(prname)//'absorb_potential.bin'), &
                              filesize)
        endif
      endif
    endif

    ! poroelastic domains
    if (POROELASTIC_SIMULATION) then
      ! allocates wavefields for solid and fluid phases
      if (b_num_abs_boundary_faces > 0) then
        allocate(b_absorb_fields(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2145')
        allocate(b_absorb_fieldw(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2146')
        if (ier /= 0) stop 'error allocating array b_absorb_fields and b_absorb_fieldw'
      else
        ! dummy
        allocate(b_absorb_fields(1,1,1))
        allocate(b_absorb_fieldw(1,1,1))
      endif
      b_absorb_fields(:,:,:) = 0.0_CUSTOM_REAL
      b_absorb_fieldw(:,:,:) = 0.0_CUSTOM_REAL

      if (num_abs_boundary_faces > 0 .and. SAVE_STACEY) then
        ! size of single record
        b_reclen_field_poro = CUSTOM_REAL * NDIM * NGLLSQUARE * num_abs_boundary_faces

        ! check integer size limit: size of b_reclen_field must fit onto an
        ! 4-byte integer
        if (num_abs_boundary_faces > int(2147483646.0 / (CUSTOM_REAL * NDIM * NGLLSQUARE))) then
          print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_field_poro
          print *,'  ',CUSTOM_REAL, NDIM, NGLLSQUARE, num_abs_boundary_faces
          print *,'bit size Fortran: ',bit_size(b_reclen_field_poro)
          call exit_MPI(myrank,"error b_reclen_field_poro integer limit")
        endif

        ! total file size
        filesize = b_reclen_field_poro
        filesize = filesize*NSTEP

        if (SIMULATION_TYPE == 3) then
          ! opens existing files
          call open_file_abs_r(IOABS,trim(prname)//'absorb_fields.bin', &
                              len_trim(trim(prname)//'absorb_fields.bin'), &
                              filesize)
          call open_file_abs_r(IOABS,trim(prname)//'absorb_fieldw.bin', &
                              len_trim(trim(prname)//'absorb_fieldw.bin'), &
                              filesize)
        else
          ! opens new file
          call open_file_abs_w(IOABS,trim(prname)//'absorb_fields.bin', &
                              len_trim(trim(prname)//'absorb_fields.bin'), &
                              filesize)
          call open_file_abs_w(IOABS,trim(prname)//'absorb_fieldw.bin', &
                              len_trim(trim(prname)//'absorb_fieldw.bin'), &
                              filesize)
        endif
      endif
    endif

  else
    ! no STACEY_ABSORBING_CONDITIONS
    b_num_abs_boundary_faces = 0

    ! dummy arrays
    if (ELASTIC_SIMULATION) then
      allocate(b_absorb_field(1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2151')
      if (ier /= 0) stop 'error allocating array b_absorb_field'
    endif

    if (ACOUSTIC_SIMULATION) then
      allocate(b_absorb_potential(1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2152')
      if (ier /= 0) stop 'error allocating array b_absorb_potential'
    endif

    if (POROELASTIC_SIMULATION) then
      allocate(b_absorb_fields(1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2153')
      allocate(b_absorb_fieldw(1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2154')
      if (ier /= 0) stop 'error allocating array b_absorb_fields and b_absorb_fieldw'
    endif
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_faults()

  use specfem_par

  use fault_solver_common, only: fault_check_mesh_resolution,USE_KELVIN_VOIGT_DAMPING
  use fault_solver_dynamic, only: BC_DYNFLT_init,SIMULATION_TYPE_DYN
  use fault_solver_kinematic, only: BC_KINFLT_init,SIMULATION_TYPE_KIN

  implicit none

  ! local parameters
  logical :: FAULT_SIMULATION_all,SIMULATION_TYPE_DYN_all,SIMULATION_TYPE_KIN_all,USE_KELVIN_VOIGT_DAMPING_all

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing fault simulation"
    call flush_IMAIN()
  endif

  ! initializes kinematic and dynamic fault solvers
  ! dynamic rupture
  call BC_DYNFLT_init(prname)

  ! kinematic rupture
  call BC_KINFLT_init(prname)

  ! checks simulation flag
  if (SIMULATION_TYPE_DYN .or. SIMULATION_TYPE_KIN) then
    ! updates flag
    FAULT_SIMULATION = .true.
  else
    FAULT_SIMULATION = .false.
  endif

  ! all processes will need to have FAULT_SIMULATION flag set if any flag is .true. somewhere
  ! (needed for proper MPI assembly)
  !
  ! note: the flag USE_KELVING_VOIGT_DAMPING
  !       can change from one process to another if a process has no fault in it.
  !       thus, it is a local flag and deals only with its single MPI process.
  !
  !       the flags SIMULATION_TYPE_DYN, SIMULATION_TYPE_KIN
  !       will be the same for all processes since all processes read in the Par_file_fault
  !       and run the corresponding fault initialization routines.
  !       we will however synchronize them just to be sure they are consistent, in case the initialization changes in future.
  !
  !       FAULT_SIMULATION will be an overall flag, which must be consistent for all MPI processes.
  !       it will determine if we need to call fault routines and assembly stages.
  call any_all_l( SIMULATION_TYPE_DYN, SIMULATION_TYPE_DYN_all )
  SIMULATION_TYPE_DYN = SIMULATION_TYPE_DYN_all

  call any_all_l( SIMULATION_TYPE_KIN, SIMULATION_TYPE_KIN_all )
  SIMULATION_TYPE_KIN = SIMULATION_TYPE_KIN_all

  call any_all_l( FAULT_SIMULATION, FAULT_SIMULATION_all )
  FAULT_SIMULATION = FAULT_SIMULATION_all

  ! just to check damping and see if any damping will be used, for user output
  call any_all_l( USE_KELVIN_VOIGT_DAMPING, USE_KELVIN_VOIGT_DAMPING_all )

  ! user output
  if (FAULT_SIMULATION) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) "  Fault simulation turned on:"
      if (SIMULATION_TYPE_DYN) write(IMAIN,*) "    dynamic   rupture simulation"
      if (SIMULATION_TYPE_KIN) write(IMAIN,*) "    kinematic rupture simulation"
      if (USE_KELVIN_VOIGT_DAMPING_all) then
        write(IMAIN,*) "    using Kelvin Voigt damping"
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  else
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  no fault simulation"
      call flush_IMAIN()
    endif
  endif

  ! user output of critical time step estimation
  if (FAULT_SIMULATION) call fault_check_mesh_resolution()

  ! simulation checks
  if (SIMULATION_TYPE_DYN .and. SIMULATION_TYPE_KIN) &
    stop 'Invalid rupture flags, cannot have dynamic and kinematic ruptures combined in a simulation yet'

  if (SIMULATION_TYPE == 3) then
    ! backward simulation support needs to be checked
!    if (FAULT_SIMULATION) &
!      stop 'Fault simulation in SIMULATION_TYPE == 3 not supported yet'
    ! damping needs to be checked if doable for adjoint/kernel simulations
!    if (USE_KELVIN_VOIGT_DAMPING_all) &
!      stop 'Using Kelvin-Voigt damping in SIMULATION_TYPE == 3 not supported yet'
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_faults
