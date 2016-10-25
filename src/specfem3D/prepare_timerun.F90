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

  subroutine prepare_timerun()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use fault_solver_dynamic, only: BC_DYNFLT_init
  use fault_solver_kinematic, only: BC_KINFLT_init
  use gravity_perturbation, only: gravity_init

  implicit none

  ! local parameters
  double precision :: tCPU
  double precision, external :: wtime

  ! get MPI starting time
  time_start = wtime()

  ! user output infos
  call prepare_timerun_user_output()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()

  ! initializes arrays
  call prepare_timerun_init_wavefield()

  ! Loading kinematic and dynamic fault solvers.
  call BC_DYNFLT_init(prname,DT,myrank)

  call BC_KINFLT_init(prname,DT,myrank)

  ! sets up arrays for gravity field
  call gravity_init()

  ! sets up time increments
  call prepare_timerun_constants()

  ! prepares attenuation arrays
  call prepare_timerun_attenuation()

  ! prepares gravity arrays
  call prepare_timerun_gravity()

  ! ZN I do not use if (USE_LDDRK) call prepare_timerun_lddrk()
  ! ZN in order to avoid the error of using unallocated arrays later on in the code,
  ! ZN since R_**_lddrk are arguments in subroutine compute_forces_viscoelastic
  call prepare_timerun_lddrk()

  ! prepares C-PML arrays
  if (PML_CONDITIONS) call prepare_timerun_pml()

  ! prepares ADJOINT simulations
  call prepare_timerun_adjoint()

  ! prepares noise simulations
  call prepare_timerun_noise()

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_timerun_GPU()

#ifdef OPENMP_MODE
  ! prepares arrays for OpenMP
  call prepare_timerun_OpenMP()
#endif

  ! compute the gravity integrals if needed
  if (GRAVITY_INTEGRALS) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...computing gravity integrals'
      call flush_IMAIN()
    endif
    call compute_gravity_integrals()
  endif

  ! elapsed time since beginning of preparation
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for preparing timerun in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'time loop:'
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
    if (ATTENUATION) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
      if (USE_OLSEN_ATTENUATION) then
        write(IMAIN,*) 'using Olsen''s attenuation'
      else
        write(IMAIN,*) 'not using Olsen''s attenuation'
      endif
    else
      write(IMAIN,*) 'no attenuation'
    endif

    write(IMAIN,*)
    if (ANISOTROPY) then
      write(IMAIN,*) 'incorporating anisotropy'
    else
      write(IMAIN,*) 'no anisotropy'
    endif

    write(IMAIN,*)
    if (APPROXIMATE_OCEAN_LOAD) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) 'no oceans'
    endif

    write(IMAIN,*)
    if (GRAVITY) then
      write(IMAIN,*) 'incorporating gravity'
    else
      write(IMAIN,*) 'no gravity'
    endif

    write(IMAIN,*)
    if (ACOUSTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating acoustic simulation'
    else
      write(IMAIN,*) 'no acoustic simulation'
    endif

    write(IMAIN,*)
    if (ELASTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating elastic simulation'
    else
      write(IMAIN,*) 'no elastic simulation'
    endif

    write(IMAIN,*)
    if (POROELASTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating poroelastic simulation'
    else
      write(IMAIN,*) 'no poroelastic simulation'
    endif

    write(IMAIN,*)
    if (MOVIE_SIMULATION) then
      write(IMAIN,*) 'incorporating movie simulation'
    else
      write(IMAIN,*) 'no movie simulation'
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
                        my_neighbours_ext_mesh)

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
                        my_neighbours_ext_mesh)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassy, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassz, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

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
                        my_neighbours_ext_mesh)
      where(rmass_ocean_load <= 0._CUSTOM_REAL) rmass_ocean_load = 1._CUSTOM_REAL
      rmass_ocean_load(:) = 1._CUSTOM_REAL / rmass_ocean_load(:)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_solid_poroelastic, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_fluid_poroelastic, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

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

  subroutine prepare_timerun_init_wavefield()

! initializes arrays

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! initialize acoustic arrays to zero
  if (ACOUSTIC_SIMULATION) then
    potential_acoustic(:) = 0._CUSTOM_REAL
    potential_dot_acoustic(:) = 0._CUSTOM_REAL
    potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if (FIX_UNDERFLOW_PROBLEM) potential_acoustic(:) = VERYSMALLVAL
  endif

  ! initialize elastic arrays to zero/verysmallvall
  if (ELASTIC_SIMULATION) then
    displ(:,:) = 0._CUSTOM_REAL
    veloc(:,:) = 0._CUSTOM_REAL
    accel(:,:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if (FIX_UNDERFLOW_PROBLEM) displ(:,:) = VERYSMALLVAL
  endif

  ! initialize poroelastic arrays to zero/verysmallvall
  if (POROELASTIC_SIMULATION) then
    displs_poroelastic(:,:) = 0._CUSTOM_REAL
    velocs_poroelastic(:,:) = 0._CUSTOM_REAL
    accels_poroelastic(:,:) = 0._CUSTOM_REAL
    displw_poroelastic(:,:) = 0._CUSTOM_REAL
    velocw_poroelastic(:,:) = 0._CUSTOM_REAL
    accelw_poroelastic(:,:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if (FIX_UNDERFLOW_PROBLEM) displs_poroelastic(:,:) = VERYSMALLVAL
    if (FIX_UNDERFLOW_PROBLEM) displw_poroelastic(:,:) = VERYSMALLVAL
  endif

  end subroutine prepare_timerun_init_wavefield


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_constants()

! precomputes constants for time integration

  use specfem_par

  implicit none

  ! local parameters
  integer :: ier

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

  ! seismograms
  if (nrec_local > 0) then
    ! allocate seismogram array
    allocate(seismograms_d(NDIM,nrec_local,NSTEP),stat=ier)
    if (ier /= 0) stop 'error allocating array seismograms_d'
    allocate(seismograms_v(NDIM,nrec_local,NSTEP),stat=ier)
    if (ier /= 0) stop 'error allocating array seismograms_v'
    allocate(seismograms_a(NDIM,nrec_local,NSTEP),stat=ier)
    if (ier /= 0) stop 'error allocating array seismograms_a'
    allocate(seismograms_p(NDIM,nrec_local,NSTEP),stat=ier)
    if (ier /= 0) stop 'error allocating array seismograms_p'

    ! initialize seismograms
    seismograms_d(:,:,:) = 0._CUSTOM_REAL
    seismograms_v(:,:,:) = 0._CUSTOM_REAL
    seismograms_a(:,:,:) = 0._CUSTOM_REAL
    seismograms_p(:,:,:) = 0._CUSTOM_REAL
  endif

  end subroutine prepare_timerun_constants

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_attenuation()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  double precision, dimension(N_SLS) :: tau_sigma_dble
  double precision :: f_c_source
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
  real(kind=CUSTOM_REAL):: scale_factorl
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: scale_factor,scale_factor_kappa

  ! if attenuation is on, shift shear moduli to center frequency of absorption period band, i.e.
  ! rescale mu to average (central) frequency for attenuation
  if (ATTENUATION) then

    ! initializes arrays
    one_minus_sum_beta(:,:,:,:) = 1._CUSTOM_REAL
    factor_common(:,:,:,:,:) = 1._CUSTOM_REAL

    allocate( scale_factor(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error allocation scale_factor')
    scale_factor(:,:,:,:) = 1._CUSTOM_REAL

    one_minus_sum_beta_kappa(:,:,:,:) = 1._CUSTOM_REAL
    factor_common_kappa(:,:,:,:,:) = 1._CUSTOM_REAL
    allocate( scale_factor_kappa(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error allocation scale_factor_kappa')
    scale_factor_kappa(:,:,:,:) = 1._CUSTOM_REAL

    ! reads in attenuation arrays
    call create_name_database(prname,myrank,LOCAL_PATH)
    if (I_should_read_the_database) then
        open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', status='old',action='read',form='unformatted',iostat=ier)
        if (ier /= 0) then
            print *,'error: could not open ',prname(1:len_trim(prname))//'attenuation.bin'
            call exit_mpi(myrank,'error opening attenuation.bin file')
        endif
    endif

    if (I_should_read_the_database) then
        read(27) ispec
        if (ispec /= NSPEC_ATTENUATION_AB) then
            close(27)
            print *,'error: attenuation file array ',ispec,'should be ',NSPEC_ATTENUATION_AB
            call exit_mpi(myrank,'error attenuation array dimensions, please recompile and rerun generate_databases')
        endif
        read(27) one_minus_sum_beta
        read(27) factor_common
        read(27) scale_factor

        read(27) one_minus_sum_beta_kappa
        read(27) factor_common_kappa
        read(27) scale_factor_kappa

        close(27)
    endif

    call bcast_all_i_for_database(ispec, 1)
    if (size(one_minus_sum_beta) > 0) call bcast_all_cr_for_database(one_minus_sum_beta(1,1,1,1), size(one_minus_sum_beta))
    if (size(factor_common) > 0) call bcast_all_cr_for_database(factor_common(1,1,1,1,1), size(factor_common))
    if (size(scale_factor) > 0) call bcast_all_cr_for_database(scale_factor(1,1,1,1), size(scale_factor))
    call bcast_all_cr_for_database(one_minus_sum_beta_kappa(1,1,1,1), size(one_minus_sum_beta_kappa))
    call bcast_all_cr_for_database(factor_common_kappa(1,1,1,1,1), size(factor_common_kappa))
    call bcast_all_cr_for_database(scale_factor_kappa(1,1,1,1), size(scale_factor_kappa))

    ! gets stress relaxation times tau_sigma, i.e.
    ! precalculates tau_sigma depending on period band (constant for all Q_mu), and
    ! determines central frequency f_c_source of attenuation period band
    call get_attenuation_constants(min_resolved_period,tau_sigma_dble, &
              f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

    ! determines alphaval,betaval,gammaval for runge-kutta scheme
    tau_sigma(:) = real(tau_sigma_dble(:),kind=CUSTOM_REAL)

    call get_attenuation_memory_values(tau_sigma,deltat,alphaval,betaval,gammaval)

    ! shifts shear moduli
    do ispec = 1,NSPEC_AB

      ! skips non elastic elements
      if (ispec_is_elastic(ispec) .eqv. .false.) cycle

      ! determines attenuation factors for each GLL point
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

            ! scales only mu moduli
            scale_factorl = scale_factor(i,j,k,ispec)
            mustore(i,j,k,ispec) = mustore(i,j,k,ispec) * scale_factorl

            ! scales kappa moduli
            scale_factorl = scale_factor_kappa(i,j,k,ispec)
            kappastore(i,j,k,ispec) = kappastore(i,j,k,ispec) * scale_factorl

          enddo
        enddo
      enddo
    enddo

    deallocate(scale_factor)
    deallocate(scale_factor_kappa)

    ! statistics
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) "Attenuation:"
      write(IMAIN,*)
      write(IMAIN,*) "The code uses a constant Q quality factor,"
      write(IMAIN,*) "but approximated based on a series of Zener standard linear solids (SLS)."
      write(IMAIN,*) "The approximation is performed in the following frequency band:"
      write(IMAIN,*)
      write(IMAIN,*) "  Reference frequency (Hz):",sngl(ATTENUATION_f0_REFERENCE)," period (s):",sngl(1.0/ATTENUATION_f0_REFERENCE)
      write(IMAIN,*) "  Frequency band min/max (Hz):",sngl(1.0/MAX_ATTENUATION_PERIOD),sngl(1.0/MIN_ATTENUATION_PERIOD)
      write(IMAIN,*) "  Period band min/max (s):",sngl(MIN_ATTENUATION_PERIOD),sngl(MAX_ATTENUATION_PERIOD)
      write(IMAIN,*) "  Logarithmic central frequency (Hz):",sngl(f_c_source)," period (s):",sngl(1.0/f_c_source)
      write(IMAIN,*) "  Using full attenuation with both Q_kappa and Q_mu."
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! clear memory variables if attenuation
    ! initialize memory variables for attenuation
    epsilondev_trace(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_xx(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_yy(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_xy(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_xz(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_yz(:,:,:,:) = 0._CUSTOM_REAL

    R_trace(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xx(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yy(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xy(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xz(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yz(:,:,:,:,:) = 0._CUSTOM_REAL

    if (FIX_UNDERFLOW_PROBLEM) then
      R_trace(:,:,:,:,:) = VERYSMALLVAL
      R_xx(:,:,:,:,:) = VERYSMALLVAL
      R_yy(:,:,:,:,:) = VERYSMALLVAL
      R_xy(:,:,:,:,:) = VERYSMALLVAL
      R_xz(:,:,:,:,:) = VERYSMALLVAL
      R_yz(:,:,:,:,:) = VERYSMALLVAL
    endif
  endif

  end subroutine prepare_timerun_attenuation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_gravity()

! precomputes gravity factors

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  double precision RICB,RCMB,RTOPDDOUBLEPRIME, &
    R80,R220,R400,R600,R670,R771,RMOHO,RMIDDLE_CRUST,ROCEAN
  double precision :: rspl_gravity(NR),gspl(NR),gspl2(NR)
  double precision :: radius,g,dg ! radius_km
  !double precision :: g_cmb_dble,g_icb_dble
  double precision :: rho,drhodr,vp,vs,Qkappa,Qmu
  integer :: nspl_gravity !int_radius
  integer :: i,j,k,iglob,ier

  ! sets up weights needed for integration of gravity
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = sngl( wxgll(i)*wygll(j)*wzgll(k) )
      enddo
    enddo
  enddo

  ! store g, rho and dg/dr=dg using normalized radius in lookup table every 100 m
  ! get density and velocity from PREM model using dummy doubling flag
  ! this assumes that the gravity perturbations are small and smooth
  ! and that we can neglect the 3D model and use PREM every 100 m in all cases
  ! this is probably a rather reasonable assumption
  if (GRAVITY) then

    ! allocates gravity arrays
    allocate( minus_deriv_gravity(NGLOB_AB), &
             minus_g(NGLOB_AB), stat=ier)
    if (ier /= 0) stop 'error allocating gravity arrays'

    ! sets up spline table
    call make_gravity(nspl_gravity,rspl_gravity,gspl,gspl2, &
                          ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                          R771,RTOPDDOUBLEPRIME,RCMB,RICB)

    ! pre-calculates gravity terms for all global points
    do iglob = 1,NGLOB_AB

      ! normalized radius ( zstore values given in m, negative values for depth)
      radius = ( R_EARTH + zstore(iglob) ) / R_EARTH
      call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g)

      ! use PREM density profile to calculate gravity (fine for other 1D models)
      call model_prem_iso(radius,rho,drhodr,vp,vs,Qkappa,Qmu, &
                        RICB,RCMB,RTOPDDOUBLEPRIME, &
                        R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

      dg = 4.0d0*rho - 2.0d0*g/radius

      ! re-dimensionalize
      g = g * R_EARTH*(PI*GRAV*RHOAV) ! in m / s^2 ( should be around 10 m/s^2)
      dg = dg * R_EARTH*(PI*GRAV*RHOAV) / R_EARTH ! gradient d/dz g , in 1/s^2

      minus_deriv_gravity(iglob) = - dg
      minus_g(iglob) = - g ! in negative z-direction

      ! debug
      !if (iglob == 1 .or. iglob == 1000 .or. iglob == 10000) then
      !  ! re-dimensionalize
      !  radius = radius * R_EARTH ! in m
      !  vp = vp * R_EARTH*dsqrt(PI*GRAV*RHOAV)  ! in m / s
      !  rho = rho  * RHOAV  ! in kg / m^3
      !  print *,'gravity: radius=',radius,'g=',g,'depth=',radius-R_EARTH
      !  print *,'vp=',vp,'rho=',rho,'kappa=',(vp**2) * rho
      !  print *,'minus_g..=',minus_g(iglob)
      !endif
    enddo

  else

    ! allocates dummy gravity arrays
    allocate( minus_deriv_gravity(0), &
             minus_g(0), stat=ier)
    if (ier /= 0) stop 'error allocating gravity arrays'

  endif

  end subroutine prepare_timerun_gravity

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
    allocate(R_xx_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
             R_yy_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
             R_xy_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
             R_xz_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
             R_yz_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS),stat=ier)
    if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

    allocate(R_trace_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK,N_SLS))
    if (ier /= 0) stop 'Error allocating array R_trace_lddrk etc.'

    if (SIMULATION_TYPE == 3) then
      allocate(b_R_xx_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
               b_R_yy_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
               b_R_xy_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
               b_R_xz_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS), &
               b_R_yz_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK ,N_SLS),stat=ier)
      if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

      allocate(b_R_trace_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK,N_SLS))
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
  use specfem_par, only: NGNOD,myrank,SIMULATION_TYPE,GPU_MODE
  use constants, only: IMAIN,NGNOD_EIGHT_CORNERS

  implicit none

  ! local parameters
  integer :: ispec,ispec_CPML,NSPEC_CPML_GLOBAL

  ! safety stops
  if (SIMULATION_TYPE /= 1) &
    stop 'Error C-PML for adjoint simulations not supported yet'
  if (GPU_MODE) &
    stop 'Error C-PML only supported in CPU mode for now'

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

  ! checks that 8-node mesh elements are used (27-node elements are not supported)
  if (NGNOD /= NGNOD_EIGHT_CORNERS) &
    stop 'error: the C-PML code works for 8-node bricks only; should be made more general'

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

  ! seismograms
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

    allocate(seismograms_eps(NDIM,NDIM,nrec_local,NSTEP),stat=ier)
    if (ier /= 0) stop 'error allocating array seismograms_eps'
    seismograms_eps(:,:,:,:) = 0._CUSTOM_REAL
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

      if (APPROXIMATE_HESS_KL) &
        hess_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL

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
          print *,'bit size fortran: ',bit_size(b_reclen_field)
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
        b_reclen_potential = CUSTOM_REAL * NGLLSQUARE * num_abs_boundary_faces

        ! check integer size limit: size of b_reclen_potential must fit onto an 4-byte integer
        if (num_abs_boundary_faces > 2147483646 / (CUSTOM_REAL * NGLLSQUARE)) then
          print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_potential
          print *,'  ',CUSTOM_REAL, NGLLSQUARE, num_abs_boundary_faces
          print *,'bit size fortran: ',bit_size(b_reclen_potential)
          call exit_MPI(myrank,"error b_reclen_potential integer limit")
        endif

        ! total file size (two lines to implicitly convert to 8-byte integers)
        filesize = b_reclen_potential
        filesize = filesize*NSTEP

        ! debug check size limit
        !if (NSTEP > 2147483646 / b_reclen_potential) then
        !  print *,'file size needed exceeds integer 4-byte limit: ',b_reclen_potential,NSTEP
        !  print *,'  ',CUSTOM_REAL, NGLLSQUARE, num_abs_boundary_faces,NSTEP
        !  print *,'file size fortran: ',filesize
        !  print *,'file bit size fortran: ',bit_size(filesize)
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
          print *,'bit size fortran: ',bit_size(b_reclen_field_poro)
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

  subroutine prepare_timerun_noise()

! prepares noise simulations

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: ier

  ! for noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then

    ! checks if free surface is defined
    if (num_free_surface_faces == 0) then
      write(*,*) myrank, " doesn't have a free_surface_face"
      ! stop 'error: noise simulations need a free surface'
    endif

    ! allocates arrays
    allocate(noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP),stat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error allocating noise source array')

    allocate(normal_x_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) stop 'error allocating array normal_x_noise'
    allocate(normal_y_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) stop 'error allocating array normal_y_noise'
    allocate(normal_z_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) stop 'error allocating array normal_z_noise'
    allocate(mask_noise(NGLLSQUARE*num_free_surface_faces),stat=ier)
    if (ier /= 0) stop 'error allocating array mask_noise'
    allocate(noise_surface_movie(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
    if (ier /= 0) stop 'error allocating array noise_surface_movie'

    ! initializes
    noise_sourcearray(:,:,:,:,:) = 0._CUSTOM_REAL
    normal_x_noise(:)            = 0._CUSTOM_REAL
    normal_y_noise(:)            = 0._CUSTOM_REAL
    normal_z_noise(:)            = 0._CUSTOM_REAL
    mask_noise(:)                = 0._CUSTOM_REAL
    noise_surface_movie(:,:,:) = 0._CUSTOM_REAL

    ! sets up noise source for master receiver station
    call read_parameters_noise(myrank,nrec,NSTEP,NGLLSQUARE*num_free_surface_faces, &
                               islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver,nu, &
                               noise_sourcearray,xigll,yigll,zigll, &
                               ibool, &
                               xstore,ystore,zstore, &
                               irec_master_noise,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                               NSPEC_AB,NGLOB_AB, &
                               num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                               ispec_is_acoustic)

    ! checks flags for noise simulation
    call check_parameters_noise(myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                LOCAL_PATH, &
                                num_free_surface_faces,NSTEP)
  endif

  end subroutine prepare_timerun_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_GPU()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use fault_solver_dynamic

  implicit none

  ! local parameters
  real :: free_mb,used_mb,total_mb


  ! GPU_MODE now defined in Par_file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "GPU Preparing Fields and Constants on Device."
    call flush_IMAIN()
  endif

  ! evaluates memory required
  call memory_eval_gpu()

  ! prepares general fields on GPU
  ! add GPU support for the C-PML routines
  call prepare_constants_device(Mesh_pointer, &
                                NGLLX, NSPEC_AB, NGLOB_AB, &
                                xix, xiy, xiz, etax,etay,etaz, gammax, gammay, gammaz, &
                                kappastore, mustore, &
                                ibool, &
                                num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
                                nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
                                hprime_xx,hprimewgll_xx, &
                                wgllwgll_xy, wgllwgll_xz, wgllwgll_yz, &
                                STACEY_ABSORBING_CONDITIONS, &
                                abs_boundary_ispec, abs_boundary_ijk, &
                                abs_boundary_normal, &
                                abs_boundary_jacobian2Dw, &
                                num_abs_boundary_faces, &
                                ispec_is_inner, &
                                NSOURCES, nsources_local, &
                                sourcearrays, islice_selected_source, ispec_selected_source, &
                                number_receiver_global, ispec_selected_rec, &
                                nrec, nrec_local, &
                                SIMULATION_TYPE, &
                                USE_MESH_COLORING_GPU, &
                                nspec_acoustic,nspec_elastic, &
                                myrank,SAVE_FORWARD)


  ! prepares fields on GPU for acoustic simulations
  if (ACOUSTIC_SIMULATION) then
    call prepare_fields_acoustic_device(Mesh_pointer, &
                                rmass_acoustic,rhostore,kappastore, &
                                num_phase_ispec_acoustic,phase_ispec_inner_acoustic, &
                                ispec_is_acoustic, &
                                NOISE_TOMOGRAPHY,num_free_surface_faces, &
                                free_surface_ispec,free_surface_ijk, &
                                b_reclen_potential,b_absorb_potential, &
                                ELASTIC_SIMULATION, num_coupling_ac_el_faces, &
                                coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                coupling_ac_el_normal,coupling_ac_el_jacobian2Dw, &
                                num_colors_outer_acoustic,num_colors_inner_acoustic, &
                                num_elem_colors_acoustic)

    if (SIMULATION_TYPE == 3) &
      call prepare_fields_acoustic_adj_dev(Mesh_pointer, &
                                APPROXIMATE_HESS_KL)

  endif

  ! prepares fields on GPU for elastic simulations
  ! add GPU support for the C-PML routines
  if (ELASTIC_SIMULATION) then
    call prepare_fields_elastic_device(Mesh_pointer, &
                                rmassx,rmassy,rmassz, &
                                rho_vp,rho_vs, &
                                num_phase_ispec_elastic,phase_ispec_inner_elastic, &
                                ispec_is_elastic, &
                                b_absorb_field,b_reclen_field, &
                                COMPUTE_AND_STORE_STRAIN, &
                                epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                epsilondev_xz,epsilondev_yz, &
                                ATTENUATION, &
                                size(R_xx), &
                                R_xx,R_yy,R_xy,R_xz,R_yz, &
                                one_minus_sum_beta,factor_common, &
                                alphaval,betaval,gammaval, &
                                APPROXIMATE_OCEAN_LOAD,rmass_ocean_load, &
                                NOISE_TOMOGRAPHY, &
                                free_surface_normal,free_surface_ispec,free_surface_ijk, &
                                num_free_surface_faces, &
                                ACOUSTIC_SIMULATION, &
                                num_colors_outer_elastic,num_colors_inner_elastic, &
                                num_elem_colors_elastic, &
                                ANISOTROPY, &
                                c11store,c12store,c13store,c14store,c15store,c16store, &
                                c22store,c23store,c24store,c25store,c26store, &
                                c33store,c34store,c35store,c36store, &
                                c44store,c45store,c46store,c55store,c56store,c66store)

    if (SIMULATION_TYPE == 3) &
      call prepare_fields_elastic_adj_dev(Mesh_pointer, &
                                NDIM*NGLOB_AB, &
                                COMPUTE_AND_STORE_STRAIN, &
                                epsilon_trace_over_3, &
                                b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                                b_epsilondev_xz,b_epsilondev_yz, &
                                b_epsilon_trace_over_3, &
                                ATTENUATION,size(R_xx), &
                                b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                                b_alphaval,b_betaval,b_gammaval, &
                                ANISOTROPIC_KL, &
                                APPROXIMATE_HESS_KL)

  endif

  ! prepares fields on GPU for poroelastic simulations
  if (POROELASTIC_SIMULATION) then
    stop 'todo poroelastic simulations on GPU'
  endif

  ! synchronizes processes
  !call synchronize_all()

  ! prepares needed receiver array for adjoint runs
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) &
    call prepare_sim2_or_3_const_device(Mesh_pointer, &
                                islice_selected_rec,size(islice_selected_rec), &
                                nadj_rec_local,nrec)

  ! prepares fields on GPU for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    ! note: noise tomography is only supported for elastic domains so far.

    ! copies noise  arrays to GPU
    call prepare_fields_noise_device(Mesh_pointer, &
                                NSPEC_AB, NGLOB_AB, &
                                free_surface_ispec, &
                                free_surface_ijk, &
                                num_free_surface_faces, &
                                NOISE_TOMOGRAPHY, &
                                NSTEP,noise_sourcearray, &
                                normal_x_noise,normal_y_noise,normal_z_noise, &
                                mask_noise,free_surface_jacobian2Dw)

  endif ! NOISE_TOMOGRAPHY

  ! prepares gravity arrays
  if (GRAVITY) then
    call prepare_fields_gravity_device(Mesh_pointer,GRAVITY, &
                                minus_deriv_gravity,minus_g,wgll_cube, &
                                ACOUSTIC_SIMULATION,rhostore)
  endif

  ! prepares Kelvin-Voigt damping around the fault
  if (SIMULATION_TYPE_DYN) then
    call Transfer_faultdata_GPU()
    call prepare_fault_device(Mesh_pointer,allocated(Kelvin_Voigt_eta),Kelvin_Voigt_eta)
    call rsf_GPU_init()
  endif

  ! synchronizes processes
  call synchronize_all()

  ! sends initial data to device

  ! puts acoustic initial fields onto GPU
  if (ACOUSTIC_SIMULATION) then
    call transfer_fields_ac_to_device(NGLOB_AB,potential_acoustic, &
                                      potential_dot_acoustic,potential_dot_dot_acoustic,Mesh_pointer)
    if (SIMULATION_TYPE == 3) &
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,b_potential_dot_dot_acoustic,Mesh_pointer)
  endif

  ! puts elastic initial fields onto GPU
  if (ELASTIC_SIMULATION) then
    ! transfer forward and backward fields to device with initial values
    call transfer_fields_el_to_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)
    if (SIMULATION_TYPE == 3) &
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
  endif

  ! synchronizes processes
  call synchronize_all()

  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if (myrank == 0) then
    call get_free_device_memory(free_mb,used_mb,total_mb)
    write(IMAIN,*) "GPU usage: free  =",free_mb," MB",nint(free_mb/total_mb*100.0),"%"
    write(IMAIN,*) "           used  =",used_mb," MB",nint(used_mb/total_mb*100.0),"%"
    write(IMAIN,*) "           total =",total_mb," MB",nint(total_mb/total_mb*100.0),"%"
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  contains

    subroutine memory_eval_gpu()

    implicit none

    ! local parameters
    double precision :: memory_size
    integer,parameter :: NGLL2 = 25
    integer,parameter :: NGLL3 = 125
    integer,parameter :: NGLL3_PADDED = 128

    memory_size = 0.d0

    ! add size of each set of arrays multiplied by the number of such arrays
    ! d_hprime_xx,d_hprimewgll_xx
    memory_size = memory_size + 2.d0 * NGLL2 * dble(CUSTOM_REAL)
    ! padded xix,..gammaz
    memory_size = memory_size + 9.d0 * NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
    ! padded kappav,muv
    memory_size = memory_size + 2.d0 * NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
    ! ibool
    memory_size = memory_size + NGLL3 * NSPEC_AB * dble(SIZE_INTEGER)
    ! d_ibool_interfaces_ext_mesh
    memory_size = memory_size + num_interfaces_ext_mesh * max_nibool_interfaces_ext_mesh * dble(SIZE_INTEGER)
    ! ispec_is_inner
    memory_size = memory_size + NSPEC_AB * dble(SIZE_INTEGER)

    if (STACEY_ABSORBING_CONDITIONS) then
      ! d_abs_boundary_ispec
      memory_size = memory_size + num_abs_boundary_faces * dble(SIZE_INTEGER)
      ! d_abs_boundary_ijk
      memory_size = memory_size + num_abs_boundary_faces * 3.d0 * NGLL2 * dble(SIZE_INTEGER)
      ! d_abs_boundary_normal
      memory_size = memory_size + num_abs_boundary_faces * NDIM * NGLL2 * dble(CUSTOM_REAL)
      ! d_abs_boundary_jacobian2Dw
      memory_size = memory_size + num_abs_boundary_faces * NGLL2 * dble(CUSTOM_REAL)
    endif

    ! sources
    ! d_sourcearrays
    memory_size = memory_size + NGLL3 * NSOURCES * NDIM * dble(CUSTOM_REAL)
    ! d_islice_selected_source,d_ispec_selected_source
    memory_size = memory_size + 2.0 * NSOURCES * dble(SIZE_INTEGER)

    ! receivers
    !d_number_receiver_global
    memory_size = memory_size + nrec_local * dble(SIZE_INTEGER)
    ! d_ispec_selected_rec
    memory_size = memory_size + nrec * dble(SIZE_INTEGER)

    ! acoustic simulations
    if (ACOUSTIC_SIMULATION) then
      ! d_potential_acoustic,d_potential_dot_acoustic,d_potential_dot_dot_acoustic
      memory_size = memory_size + 3.d0 * NGLOB_AB * dble(CUSTOM_REAL)
      ! d_rmass_acoustic
      memory_size = memory_size + NGLOB_AB * dble(CUSTOM_REAL)
      ! padded d_rhostore
      memory_size = memory_size + NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
      ! d_kappastore
      memory_size = memory_size + NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      ! d_phase_ispec_inner_acoustic
      memory_size = memory_size + 2.d0 * num_phase_ispec_acoustic * dble(SIZE_INTEGER)
      ! d_ispec_is_acoustic
      memory_size = memory_size + NSPEC_AB * dble(SIZE_INTEGER)
    endif

    ! elastic simulations
    if (ELASTIC_SIMULATION) then
      ! d_displ,..
      memory_size = memory_size + 3.d0 * NDIM * NGLOB_AB * dble(CUSTOM_REAL)
      ! d_send_accel_buffer
      memory_size = memory_size + 3.d0 * num_interfaces_ext_mesh * max_nibool_interfaces_ext_mesh * dble(CUSTOM_REAL)
      ! d_rmassx,..
      memory_size = memory_size + 3.d0 * NGLOB_AB * dble(CUSTOM_REAL)
      ! d_ispec_is_elastic
      memory_size = memory_size + NSPEC_AB * dble(SIZE_INTEGER)
      ! d_phase_ispec_inner_elastic
      memory_size = memory_size + 2.d0 * num_phase_ispec_elastic * dble(SIZE_INTEGER)
      ! d_station_seismo_field
      memory_size = memory_size + 3.d0 * NGLL3 * nrec_local * dble(CUSTOM_REAL)

      if (STACEY_ABSORBING_CONDITIONS .or. PML_CONDITIONS) then
        ! d_rho_vp,..
        memory_size = memory_size + 2.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      endif
      if (COMPUTE_AND_STORE_STRAIN) then
        ! d_epsilondev_xx,..
        memory_size = memory_size + 5.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      endif
      if (ATTENUATION) then
        ! d_R_xx,..
        memory_size = memory_size + 5.d0 * size(R_xx) * dble(CUSTOM_REAL)
        ! d_one_minus_sum_beta
        memory_size = memory_size + NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
        ! d_factor_common
        memory_size = memory_size + N_SLS * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
        ! alphaval,..
        memory_size = memory_size + 3.d0 * N_SLS * dble(CUSTOM_REAL)
      endif
      if (ANISOTROPY) then
        ! padded d_c11store,..
        memory_size = memory_size + 21.d0 * NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
      endif
      if (APPROXIMATE_OCEAN_LOAD) then
        ! d_rmass_ocean_load
        memory_size = memory_size + NGLOB_AB * dble(CUSTOM_REAL)
        ! d_free_surface_normal
        memory_size = memory_size + 3.d0 * NGLL2 * num_free_surface_faces * dble(CUSTOM_REAL)
      endif
    endif

    ! noise simulations
    if (NOISE_TOMOGRAPHY > 0) then
      ! d_free_surface_ispec
      memory_size = memory_size + num_free_surface_faces * dble(SIZE_INTEGER)
      ! d_free_surface_ijk
      memory_size = memory_size + 3.d0 * NGLL2 * num_free_surface_faces * dble(SIZE_INTEGER)
      ! d_noise_surface_movie
      memory_size = memory_size + 3.d0 * NGLL2 * num_free_surface_faces * dble(CUSTOM_REAL)
      if (NOISE_TOMOGRAPHY == 1) then
        ! d_noise_sourcearray
        memory_size = memory_size + 3.d0 * NGLL3 * NSTEP * dble(CUSTOM_REAL)
      endif
      if (NOISE_TOMOGRAPHY > 1) then
        ! d_normal_x_noise,..
        memory_size = memory_size + 5.d0 * NGLL2 * num_free_surface_faces * dble(CUSTOM_REAL)
      endif
      if (NOISE_TOMOGRAPHY == 3) then
        ! d_sigma_kl
        memory_size = memory_size + NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      endif
    endif

    if (GRAVITY) then
      ! d_minus_deriv_gravity,d_minus_g
      memory_size = memory_size + 2.d0 * NGLOB_AB * dble(CUSTOM_REAL)
    endif

    ! poor estimate for kernel simulations...
    if (SIMULATION_TYPE == 3) memory_size = 2.d0 * memory_size

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  minimum memory requested     : ', &
                    memory_size / 1024. / 1024., &
                     'MB per process'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    end subroutine memory_eval_gpu

  end subroutine prepare_timerun_GPU

!
!-------------------------------------------------------------------------------------------------
!

! OpenMP version uses "special" compute_forces_viscoelastic routine
! we need to set num_elem_colors_elastic arrays

#ifdef OPENMP_MODE
  subroutine prepare_timerun_OpenMP()

  use specfem_par
  use specfem_par_elastic

  implicit none

  ! local parameters
  integer :: ier
  integer :: NUM_THREADS
  integer :: OMP_GET_MAX_THREADS

  ! OpenMP for elastic simulation only supported yet
  if (ELASTIC_SIMULATION) then

    NUM_THREADS = OMP_GET_MAX_THREADS()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Using:',NUM_THREADS, ' OpenMP threads'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! allocate cfe_openmp local arrays for OpenMP version
!! DK DK July 2014: I do not know who wrote the OpenMP version, but it is currently broken
!! DK DK July 2014: because the arrays below are undeclared; I therefore need to comment them out
!! DK DK July 2014: for now and put a stop statement instead
    stop 'from DK DK, July 2014: the OpenMP version is currently broken here, not sure who wrote it, please fix it if possible'
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

!
!-------------------------------------------------------------------------------------------------
!

  ! compute gravity integrals of that slice, and then total integrals for the whole mesh

  subroutine compute_gravity_integrals()

  use constants

  use specfem_par

  implicit none

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  double precision :: jacobianl
  integer :: i,j,k,ispec,iglob,ier
  double precision :: xval,yval,zval
  double precision :: xval_squared,yval_squared,zval_squared
  double precision :: x_meshpoint,y_meshpoint,z_meshpoint
  double precision :: distance_squared,distance_cubed, &
                      three_over_distance_squared,one_over_distance_cubed,three_over_distance_fifth_power
  double precision :: common_multiplying_factor,common_mult_times_one_over,common_mult_times_three_over

  integer :: iobservation
  character(len=MAX_STRING_LEN) :: outputname

! read the observation surface
  x_observation(:) = 0.d0
  y_observation(:) = 0.d0
  z_observation(:) = 0.d0
  if (myrank == 0) then
    open(unit=IIN,file=OBSERVATION_GRID_FILE(1:len_trim(OBSERVATION_GRID_FILE)),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file observation_grid_to_use_for_gravity.txt')
    do iobservation = 1,NTOTAL_OBSERVATION
      read(IIN,*) x_observation(iobservation),y_observation(iobservation),z_observation(iobservation)
    enddo
    close(unit=IIN)
  endif

! broadcast the observation surface read
  call bcast_all_dp(x_observation, NTOTAL_OBSERVATION)
  call bcast_all_dp(y_observation, NTOTAL_OBSERVATION)
  call bcast_all_dp(z_observation, NTOTAL_OBSERVATION)

! initialize the gravity arrays
  g_x(:) = 0.d0
  g_y(:) = 0.d0
  g_z(:) = 0.d0

  G_xx(:) = 0.d0
  G_yy(:) = 0.d0
  G_zz(:) = 0.d0
  G_xy(:) = 0.d0
  G_xz(:) = 0.d0
  G_yz(:) = 0.d0

  ! calculates volume of all elements in mesh
  do ispec = 1,NSPEC_AB

    ! print information about number of elements done so far
    if (myrank == 0 .and. (mod(ispec,NSPEC_DISPLAY_INTERVAL) == 0 .or. ispec == 1 .or. ispec == NSPEC_AB)) then
       write(IMAIN,*) 'for gravity integrals, ',ispec,' elements computed out of ',NSPEC_AB
       ! write time stamp file to give information about progression of the calculation of gravity integrals
       write(outputname,"('/timestamp_gravity_calculations_ispec',i7.7,'_out_of_',i7.7)") ispec,NSPEC_AB
       ! timestamp file output
       open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write')
       write(IOUT,*) ispec,' elements done for gravity calculations out of ',NSPEC_AB
       close(unit=IOUT)
    endif

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

          ! compute the Jacobian
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! do this in double precision for accuracy
          jacobianl = 1.d0 / dble(xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianl <= ZERO) stop 'error: negative Jacobian found in integral calculation'

          iglob = ibool(i,j,k,ispec)
          x_meshpoint = xstore(iglob)
          y_meshpoint = ystore(iglob)
          z_meshpoint = zstore(iglob)

          common_multiplying_factor = jacobianl * weight * rhostore(i,j,k,ispec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! beginning of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! loop on all the points in the observation surface
  do iobservation = 1,NTOTAL_OBSERVATION

    xval = x_meshpoint - x_observation(iobservation)
    yval = y_meshpoint - y_observation(iobservation)
    zval = z_meshpoint - z_observation(iobservation)

    xval_squared = xval**2
    yval_squared = yval**2
    zval_squared = zval**2

    distance_squared = xval_squared + yval_squared + zval_squared
    distance_cubed = distance_squared * sqrt(distance_squared)

    three_over_distance_squared = 3.d0 / distance_squared
    one_over_distance_cubed = 1.d0 / distance_cubed
    three_over_distance_fifth_power = three_over_distance_squared * one_over_distance_cubed

    common_mult_times_one_over = common_multiplying_factor * one_over_distance_cubed
    common_mult_times_three_over = common_multiplying_factor * three_over_distance_fifth_power

    g_x(iobservation) = g_x(iobservation) + common_mult_times_one_over * xval
    g_y(iobservation) = g_y(iobservation) + common_mult_times_one_over * yval
    g_z(iobservation) = g_z(iobservation) + common_mult_times_one_over * zval

    G_xx(iobservation) = G_xx(iobservation) + common_mult_times_one_over * (xval_squared * three_over_distance_squared - 1.d0)
    G_yy(iobservation) = G_yy(iobservation) + common_mult_times_one_over * (yval_squared * three_over_distance_squared - 1.d0)
    G_zz(iobservation) = G_zz(iobservation) + common_mult_times_one_over * (zval_squared * three_over_distance_squared - 1.d0)

    G_xy(iobservation) = G_xy(iobservation) + common_mult_times_three_over * xval*yval
    G_xz(iobservation) = G_xz(iobservation) + common_mult_times_three_over * xval*zval
    G_yz(iobservation) = G_yz(iobservation) + common_mult_times_three_over * yval*zval

  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! end of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
      enddo
    enddo
  enddo

    ! the result is displayed in Eotvos = 1.e+9 s-2
    G_xx(:) = G_xx(:) * SI_UNITS_TO_EOTVOS
    G_yy(:) = G_yy(:) * SI_UNITS_TO_EOTVOS
    G_zz(:) = G_zz(:) * SI_UNITS_TO_EOTVOS
    G_xy(:) = G_xy(:) * SI_UNITS_TO_EOTVOS
    G_xz(:) = G_xz(:) * SI_UNITS_TO_EOTVOS
    G_yz(:) = G_yz(:) * SI_UNITS_TO_EOTVOS

    ! use an MPI reduction to compute the total value of the integral into a temporary array
    ! and then copy it back into the original array
    call sum_all_1Darray_dp(g_x,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) g_x(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(g_y,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) g_y(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(g_z,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) g_z(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_xx,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_xx(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_yy,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_yy(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_zz,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_zz(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_xy,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_xy(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_xz,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_xz(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_yz,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_yz(:) = temporary_array_for_sum(:)

  !--- print number of points and elements in the mesh for each region
  if (myrank == 0) then

      temporary_array_for_sum(:) = sqrt(g_x(:)**2 + g_y(:)**2 + g_z(:)**2)
      write(IMAIN,*)
      write(IMAIN,*) 'minval of norm of g vector on whole observation surface = ',minval(temporary_array_for_sum),' m.s-2'
      write(IMAIN,*) 'maxval of norm of g vector on whole observation surface = ',maxval(temporary_array_for_sum),' m.s-2'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xx on whole observation surface = ',minval(G_xx),' Eotvos'
      write(IMAIN,*) 'maxval of G_xx on whole observation surface = ',maxval(G_xx),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_yy on whole observation surface = ',minval(G_yy),' Eotvos'
      write(IMAIN,*) 'maxval of G_yy on whole observation surface = ',maxval(G_yy),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_zz on whole observation surface = ',minval(G_zz),' Eotvos'
      write(IMAIN,*) 'maxval of G_zz on whole observation surface = ',maxval(G_zz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xy on whole observation surface = ',minval(G_xy),' Eotvos'
      write(IMAIN,*) 'maxval of G_xy on whole observation surface = ',maxval(G_xy),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xz on whole observation surface = ',minval(G_xz),' Eotvos'
      write(IMAIN,*) 'maxval of G_xz on whole observation surface = ',maxval(G_xz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_yz on whole observation surface = ',minval(G_yz),' Eotvos'
      write(IMAIN,*) 'maxval of G_yz on whole observation surface = ',maxval(G_yz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'Minval and maxval of trace of G, which in principle should be zero:'
      write(IMAIN,*)
      temporary_array_for_sum(:) = abs(G_xx(:) + G_yy(:) + G_zz(:))
      write(IMAIN,*) 'minval of abs(G_xx + G_yy + G_zz) on whole observation surface = ',minval(temporary_array_for_sum),' Eotvos'
      write(IMAIN,*) 'maxval of abs(G_xx + G_yy + G_zz) on whole observation surface = ',maxval(temporary_array_for_sum),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) '-----------------------------'
      write(IMAIN,*)
      write(IMAIN,*) 'displaying the fields computed at observation point = ',iobs_receiver,' out of ',NTOTAL_OBSERVATION
      write(IMAIN,*)
      write(IMAIN,*) 'computed g_x  = ',g_x(iobs_receiver),' m.s-2'
      write(IMAIN,*) 'computed g_y  = ',g_y(iobs_receiver),' m.s-2'
      write(IMAIN,*) 'computed g_z  = ',g_z(iobs_receiver),' m.s-2'
      write(IMAIN,*)
      write(IMAIN,*) 'computed norm of g vector = ',sqrt(g_x(iobs_receiver)**2 + g_y(iobs_receiver)**2 + &
                                                                 g_z(iobs_receiver)**2),' m.s-2'

      write(IMAIN,*)
      write(IMAIN,*) 'computed G_xx = ',G_xx(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_yy = ',G_yy(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_zz = ',G_zz(iobs_receiver),' Eotvos'
      write(IMAIN,*)
      write(IMAIN,*) 'G tensor should be traceless, G_xx + G_yy + G_zz = 0.'
      write(IMAIN,*) 'Actual sum obtained = ',G_xx(iobs_receiver) + G_yy(iobs_receiver) + G_zz(iobs_receiver)
      if (max(abs(G_xx(iobs_receiver)),abs(G_yy(iobs_receiver)),abs(G_zz(iobs_receiver))) > TINYVAL) &
           write(IMAIN,*) ' i.e., ',sngl(100.d0*abs(G_xx(iobs_receiver) + G_yy(iobs_receiver) + G_zz(iobs_receiver)) / &
                                     max(abs(G_xx(iobs_receiver)),abs(G_yy(iobs_receiver)),abs(G_zz(iobs_receiver)))), &
                                     '% of max(abs(G_xx),abs(G_yy),abs(G_zz))'
      write(IMAIN,*)
      write(IMAIN,*) 'computed G_xy = ',G_xy(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_xz = ',G_xz(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_yz = ',G_yz(iobs_receiver),' Eotvos'

      ! save the results
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_x_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) g_x(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_y_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) g_y(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_z_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) g_z(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_norm_of_g_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) sqrt(g_x(iobservation)**2 + g_y(iobservation)**2 + g_z(iobservation)**2)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xx_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_xx(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_yy_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_yy(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_zz_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_zz(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xy_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_xy(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xz_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_xz(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_yz_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_yz(iobservation)
      enddo
      close(unit=IOUT)

  endif

  end subroutine compute_gravity_integrals

