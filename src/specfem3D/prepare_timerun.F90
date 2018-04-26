!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

  subroutine prepare_timerun()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use fault_solver_dynamic, only: BC_DYNFLT_init
  use fault_solver_kinematic, only: BC_KINFLT_init

  implicit none

  ! local parameters
  double precision :: tCPU,tstart
  double precision, external :: wtime

  ! get MPI starting time
  tstart = wtime()

  ! user output infos
  call prepare_timerun_user_output()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()

  ! initializes arrays
  call prepare_wavefields()

  ! Loading kinematic and dynamic fault solvers.
  call BC_DYNFLT_init(prname,DT,myrank)

  call BC_KINFLT_init(prname,DT,myrank)

  ! sets up time increments
  call prepare_timerun_constants()

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

  ! prepares ADJOINT simulations
  call prepare_timerun_adjoint()

  ! prepares noise simulations
  call prepare_noise()

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_GPU()

#ifdef USE_OPENMP
  ! prepares arrays for OpenMP
  call prepare_timerun_OpenMP()
#endif

  ! prepars coupling with injection boundary
  call couple_with_injection_prepare_boundary()

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

  ! flag for any movie simulation
  if (MOVIE_SURFACE .or. CREATE_SHAKEMAP .or. MOVIE_VOLUME .or. PNM_IMAGE) then
    MOVIE_SIMULATION = .true.
  else
    MOVIE_SIMULATION = .false.
  endif

  ! user info
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Simulation setup:'
    write(IMAIN,*)

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

    !! CD CD !!
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
    !! CD CD

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

  end subroutine prepare_timerun_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!


  subroutine prepare_timerun_constants()

! precomputes constants for time integration

  use specfem_par

  implicit none

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
    ! backward/reconstructed wavefields: time stepping is in time-reversed sense
    ! (negative time increments)
    b_deltat = - real(DT,kind=CUSTOM_REAL)
    b_deltatover2 = b_deltat/2._CUSTOM_REAL
    b_deltatsqover2 = b_deltat*b_deltat/2._CUSTOM_REAL
  endif

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

  if (ACOUSTIC_SIMULATION) then
    allocate(potential_acoustic_lddrk(NGLOB_AB_LDDRK),stat=ier)
    if (ier /= 0) stop 'Error allocating array potential_acoustic_lddrk'
    allocate(potential_dot_acoustic_lddrk(NGLOB_AB_LDDRK),stat=ier)
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
    if (ier /= 0) stop 'Error allocating array displ_lddrk'
    allocate(veloc_lddrk(NDIM,NGLOB_AB_LDDRK),stat=ier)
    if (ier /= 0) stop 'Error allocating array veloc_lddrk'

    displ_lddrk(:,:) = 0._CUSTOM_REAL
    veloc_lddrk(:,:) = 0._CUSTOM_REAL
    if (FIX_UNDERFLOW_PROBLEM) then
      displ_lddrk(:,:) = VERYSMALLVAL
      veloc_lddrk(:,:) = VERYSMALLVAL
    endif

    ! note: currently, they need to be defined, as they are used in some subroutine arguments
    allocate(R_xx_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
             R_yy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
             R_xy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
             R_xz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
             R_yz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
    if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

    allocate(R_trace_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK))
    if (ier /= 0) stop 'Error allocating array R_trace_lddrk etc.'

    if (SIMULATION_TYPE == 3) then
      allocate(b_R_xx_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
               b_R_yy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
               b_R_xy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
               b_R_xz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK), &
               b_R_yz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
      if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

      allocate(b_R_trace_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK))
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

  end subroutine prepare_timerun_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_pml()

  use pml_par
  use specfem_par, only: myrank,SIMULATION_TYPE,GPU_MODE,UNDO_ATTENUATION_AND_OR_PML
  use constants, only: IMAIN,NGNOD_EIGHT_CORNERS

  implicit none

  ! local parameters
  integer :: ispec,ispec_CPML,NSPEC_CPML_GLOBAL

  ! safety stops
  if (SIMULATION_TYPE /= 1 .and. .not. UNDO_ATTENUATION_AND_OR_PML) &
          stop 'Error: PMLs for adjoint runs require the flag UNDO_ATTENUATION_AND_OR_PML to be set'

  if (GPU_MODE) stop 'Error: PMLs only supported in CPU mode for now'

  ! total number of PML elements
  call sum_all_i(NSPEC_CPML,NSPEC_CPML_GLOBAL)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'incorporating C-PML  '
    write(IMAIN,*)
    write(IMAIN,*) 'number of C-PML spectral elements in the global mesh: ', NSPEC_CPML_GLOBAL
    write(IMAIN,*)
    write(IMAIN,*) 'thickness of C-PML layer in X direction: ', CPML_width_x
    write(IMAIN,*) 'thickness of C-PML layer in Y direction: ', CPML_width_y
    write(IMAIN,*) 'thickness of C-PML layer in Z direction: ', CPML_width_z
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
  do ispec_CPML=1,NSPEC_CPML

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
  integer(kind=8) :: filesize

  ! moment tensor derivatives
  if (nrec_local > 0 .and. SIMULATION_TYPE == 2) then
    ! allocate Frechet derivatives array
    allocate(Mxx_der(nrec_local),Myy_der(nrec_local), &
             Mzz_der(nrec_local),Mxy_der(nrec_local), &
             Mxz_der(nrec_local),Myz_der(nrec_local), &
             sloc_der(NDIM,nrec_local),stat=ier)
    if (ier /= 0) stop 'error allocating array Mxx_der and following arrays'
    Mxx_der = 0._CUSTOM_REAL
    Myy_der = 0._CUSTOM_REAL
    Mzz_der = 0._CUSTOM_REAL
    Mxy_der = 0._CUSTOM_REAL
    Mxz_der = 0._CUSTOM_REAL
    Myz_der = 0._CUSTOM_REAL
    sloc_der = 0._CUSTOM_REAL
  endif

  ! attenuation backward memories
  if (ATTENUATION .and. SIMULATION_TYPE == 3) then
    ! precompute Runge-Kutta coefficients if attenuation
    call get_attenuation_memory_values(tau_sigma,b_deltat,b_alphaval,b_betaval,b_gammaval)
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

      ! memory variables if attenuation
      if (ATTENUATION) then
        b_R_trace = 0._CUSTOM_REAL
        b_R_xx = 0._CUSTOM_REAL
        b_R_yy = 0._CUSTOM_REAL
        b_R_xy = 0._CUSTOM_REAL
        b_R_xz = 0._CUSTOM_REAL
        b_R_yz = 0._CUSTOM_REAL
        b_epsilondev_trace = 0._CUSTOM_REAL
        b_epsilondev_xx = 0._CUSTOM_REAL
        b_epsilondev_yy = 0._CUSTOM_REAL
        b_epsilondev_xy = 0._CUSTOM_REAL
        b_epsilondev_xz = 0._CUSTOM_REAL
        b_epsilondev_yz = 0._CUSTOM_REAL
      endif

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
    ispec2D_moho_top = 0
    ispec2D_moho_bot = 0
  endif

! stacey absorbing fields will be reconstructed for adjoint simulations
! using snapshot files of wavefields
  if (STACEY_ABSORBING_CONDITIONS) then

    ! opens absorbing wavefield saved/to-be-saved by forward simulations
    if (num_abs_boundary_faces > 0 .and. (SIMULATION_TYPE == 3 .or. &
          (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

      b_num_abs_boundary_faces = num_abs_boundary_faces

      ! elastic domains
      if (ELASTIC_SIMULATION) then
        ! allocates wavefield
        allocate(b_absorb_field(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) stop 'error allocating array b_absorb_field'

        ! size of single record
        b_reclen_field = CUSTOM_REAL * NDIM * NGLLSQUARE * num_abs_boundary_faces

        ! check integer size limit: size of b_reclen_field must fit onto an 4-byte integer
        if (num_abs_boundary_faces > 2147483646 / (CUSTOM_REAL * NDIM * NGLLSQUARE)) then
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

      ! acoustic domains
      if (ACOUSTIC_SIMULATION) then
        ! allocates wavefield
        allocate(b_absorb_potential(NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) stop 'error allocating array b_absorb_potential'

        ! size of single record
        b_reclen_potential = CUSTOM_REAL * NGLLSQUARE * num_abs_boundary_faces * NB_RUNS_ACOUSTIC_GPU


        ! check integer size limit: size of b_reclen_potential must fit onto an 4-byte integer
        if (num_abs_boundary_faces > 2147483646 / (CUSTOM_REAL * NGLLSQUARE)) then
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

      ! poroelastic domains
      if (POROELASTIC_SIMULATION) then
        ! allocates wavefields for solid and fluid phases
        allocate(b_absorb_fields(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        allocate(b_absorb_fieldw(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) stop 'error allocating array b_absorb_fields and b_absorb_fieldw'

        ! size of single record
        b_reclen_field_poro = CUSTOM_REAL * NDIM * NGLLSQUARE * num_abs_boundary_faces

        ! check integer size limit: size of b_reclen_field must fit onto an
        ! 4-byte integer
        if (num_abs_boundary_faces > 2147483646 / (CUSTOM_REAL * NDIM * NGLLSQUARE)) then
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
    else
      ! num_abs_boundary_faces is zero
      ! needs dummy array
      b_num_abs_boundary_faces = 0
      if (ELASTIC_SIMULATION) then
        allocate(b_absorb_field(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) stop 'error allocating array b_absorb_field'
      endif

      if (ACOUSTIC_SIMULATION) then
        allocate(b_absorb_potential(NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) stop 'error allocating array b_absorb_potential'
      endif

      if (POROELASTIC_SIMULATION) then
        allocate(b_absorb_fields(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        allocate(b_absorb_fieldw(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
        if (ier /= 0) stop 'error allocating array b_absorb_fields and b_absorb_fieldw'
      endif
    endif
  else ! STACEY_ABSORBING_CONDITIONS
    ! needs dummy array
    b_num_abs_boundary_faces = 0
    if (ELASTIC_SIMULATION) then
      allocate(b_absorb_field(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
      if (ier /= 0) stop 'error allocating array b_absorb_field'
    endif

    if (ACOUSTIC_SIMULATION) then
      allocate(b_absorb_potential(NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
      if (ier /= 0) stop 'error allocating array b_absorb_potential'
    endif

    if (POROELASTIC_SIMULATION) then
      allocate(b_absorb_fields(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
      allocate(b_absorb_fieldw(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),stat=ier)
      if (ier /= 0) stop 'error allocating array b_absorb_fields and b_absorb_fieldw'
    endif
  endif

  end subroutine prepare_timerun_adjoint

!
!-------------------------------------------------------------------------------------------------
!


! OpenMP version uses "special" compute_forces_viscoelastic routine
! we need to set num_elem_colors_elastic arrays

#ifdef USE_OPENMP
  subroutine prepare_timerun_OpenMP()

  use specfem_par
  use specfem_par_elastic

  implicit none

  ! local parameters
  integer :: ier
  integer :: NUM_THREADS
  integer :: OMP_GET_MAX_THREADS

  ! safety stop
  print *
  print *,'Note: there is currently no OpenMP version of the code!'
  print *,'      a former implementation has been broken and moved to utils/unused_routines/ directory'
  print *,'      Please re-compile the code, disenabling OpenMP support and/or consider contributing a new version if needed.'
  print *

! the old OpenMP implementation for compute_forces_viscoelastic is in utils/unused_routines/:
! older_please_do_not_use_anymore_partial_OpenMP_port/older_not_maintained_compute_forces_viscoelastic_Dev_openmp.f90

  stop 'OpenMP version is currently not implemented.'

  ! unused below but might be helpful in future...

  NUM_THREADS = OMP_GET_MAX_THREADS()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Using:',NUM_THREADS, ' OpenMP threads'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! OpenMP for elastic simulation only supported yet
  if (ELASTIC_SIMULATION) then

!! DK DK July 2014: I do not know who wrote the OpenMP version, but it is currently broken
!! DK DK July 2014: because the arrays below are undeclared; I therefore need to comment them out
!! DK DK July 2014: for now and put a stop statement instead
    stop 'from DK DK, July 2014: the OpenMP version is currently broken here, not sure who wrote it, please fix it if possible'

!   ! allocate cfe_openmp local arrays for OpenMP version
!   allocate(dummyx_loc(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(dummyy_loc(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(dummyz_loc(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(dummyx_loc_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(dummyy_loc_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(dummyz_loc_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempx1(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempx2(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempx3(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempy1(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempy2(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempy3(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempz1(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempz2(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(newtempz3(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempx1(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempx2(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempx3(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempy1(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempy2(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempy3(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempz1(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempz2(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempz3(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempx1_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempx2_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempx3_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempy1_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempy2_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempy3_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempz1_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempz2_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))
!   allocate(tempz3_att(NGLLX,NGLLY,NGLLZ,NUM_THREADS))

    ! set num_elem_colors array in case no mesh coloring is used
    if (.not. USE_MESH_COLORING_GPU) then
      ! deallocate dummy array
      if (allocated(num_elem_colors_elastic)) deallocate(num_elem_colors_elastic)

      ! loads with corresonding values
      num_colors_outer_elastic = 1
      num_colors_inner_elastic = 1
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) stop 'error allocating num_elem_colors_elastic array'

      ! sets to all elements in inner/outer phase
      num_elem_colors_elastic(1) = nspec_outer_elastic
      num_elem_colors_elastic(2) = nspec_inner_elastic
    endif

  endif

  end subroutine prepare_timerun_OpenMP
#endif

