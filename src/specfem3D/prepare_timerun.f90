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

  subroutine prepare_timerun()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie

  implicit none
  character(len=256) :: plot_file

  ! flag for any movie simulation
  if( EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP .or. &
     MOVIE_SURFACE .or. CREATE_SHAKEMAP .or. MOVIE_VOLUME .or. PNM_GIF_IMAGE ) then
    MOVIE_SIMULATION = .true.
  else
    MOVIE_SIMULATION = .false.
  endif

  ! user info
  if(myrank == 0) then

    write(IMAIN,*)
    if(ATTENUATION) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
      if(USE_OLSEN_ATTENUATION) then
        write(IMAIN,*) 'using Olsen''s attenuation'
      else
        write(IMAIN,*) 'not using Olsen''s attenuation'
      endif
    else
      write(IMAIN,*) 'no attenuation'
    endif

    write(IMAIN,*)
    if(ANISOTROPY) then
      write(IMAIN,*) 'incorporating anisotropy'
    else
      write(IMAIN,*) 'no anisotropy'
    endif

    write(IMAIN,*)
    if(OCEANS) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) 'no oceans'
    endif

    write(IMAIN,*)
    if(ACOUSTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating acoustic simulation'
    else
      write(IMAIN,*) 'no acoustic simulation'
    endif

    write(IMAIN,*)
    if(ELASTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating elastic simulation'
    else
      write(IMAIN,*) 'no elastic simulation'
    endif

    write(IMAIN,*)
    if(POROELASTIC_SIMULATION) then
      write(IMAIN,*) 'incorporating poroelastic simulation'
    else
      write(IMAIN,*) 'no poroelastic simulation'
    endif
    write(IMAIN,*)

    write(IMAIN,*)
    if(MOVIE_SIMULATION) then
      write(IMAIN,*) 'incorporating movie simulation'
    else
      write(IMAIN,*) 'no movie simulation'
    endif
    write(IMAIN,*)

  endif

  ! synchronize all the processes before assembling the mass matrix
  ! to make sure all the nodes have finished to read their databases
  call sync_all()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()


  ! initialize acoustic arrays to zero
  if( ACOUSTIC_SIMULATION ) then
    potential_acoustic(:) = 0._CUSTOM_REAL
    potential_dot_acoustic(:) = 0._CUSTOM_REAL
    potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if(FIX_UNDERFLOW_PROBLEM) potential_acoustic(:) = VERYSMALLVAL
  endif

  ! initialize elastic arrays to zero/verysmallvall
  if( ELASTIC_SIMULATION ) then
    displ(:,:) = 0._CUSTOM_REAL
    veloc(:,:) = 0._CUSTOM_REAL
    accel(:,:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if(FIX_UNDERFLOW_PROBLEM) displ(:,:) = VERYSMALLVAL
  endif


  ! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    deltat = sngl(DT)
  else
    deltat = DT
  endif
  deltatover2 = deltat/2._CUSTOM_REAL
  deltatsqover2 = deltat*deltat/2._CUSTOM_REAL

  ! seismograms
  if (nrec_local > 0) then
    ! allocate seismogram array
    allocate(seismograms_d(NDIM,nrec_local,NSTEP))
    allocate(seismograms_v(NDIM,nrec_local,NSTEP))
    allocate(seismograms_a(NDIM,nrec_local,NSTEP))

    ! initialize seismograms
    seismograms_d(:,:,:) = 0._CUSTOM_REAL
    seismograms_v(:,:,:) = 0._CUSTOM_REAL
    seismograms_a(:,:,:) = 0._CUSTOM_REAL
  endif

  ! prepares attenuation arrays
  call prepare_timerun_attenuation()

  ! initializes PML arrays
  if( ABSORBING_CONDITIONS  ) then
    if (SIMULATION_TYPE /= 1 .and. ABSORB_USE_PML )  then
      write(IMAIN,*) 'NOTE: adjoint simulations and PML not supported yet...'
    else
      if( ABSORB_USE_PML ) then
        call PML_initialize()
      endif
    endif
  endif

  ! opens source time function file
  if(PRINT_SOURCE_TIME_FUNCTION .and. myrank == 0) then
    ! print the source-time function
    if(NSOURCES == 1) then
      plot_file = '/plot_source_time_function.txt'
    else
     if(NSOURCES < 10) then
        write(plot_file,"('/plot_source_time_function',i1,'.txt')") NSOURCES
      else
        write(plot_file,"('/plot_source_time_function',i2,'.txt')") NSOURCES
      endif
    endif
    open(unit=IOSTF,file=trim(OUTPUT_FILES)//plot_file,status='unknown')
  endif

  ! user output
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(NSTEP*DT),' seconds'
    write(IMAIN,*)
  endif

  ! prepares ADJOINT simulations
  call prepare_timerun_adjoint()

  end subroutine prepare_timerun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_mass_matrices()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  implicit none

! the mass matrix needs to be assembled with MPI here once and for all
  if(ACOUSTIC_SIMULATION) then
    call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,rmass_acoustic, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh)

    ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmass_acoustic <= 0._CUSTOM_REAL) rmass_acoustic = 1._CUSTOM_REAL
    rmass_acoustic(:) = 1._CUSTOM_REAL / rmass_acoustic(:)

  endif ! ACOUSTIC_SIMULATION

  if(ELASTIC_SIMULATION) then
    call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,rmass, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

    ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmass <= 0._CUSTOM_REAL) rmass = 1._CUSTOM_REAL
    rmass(:) = 1._CUSTOM_REAL / rmass(:)

    if(OCEANS ) then
      if( minval(rmass_ocean_load(:)) <= 0._CUSTOM_REAL) &
        call exit_MPI(myrank,'negative ocean load mass matrix term')
      rmass_ocean_load(:) = 1. / rmass_ocean_load(:)
    endif

  endif ! ELASTIC_SIMULATION

  if(POROELASTIC_SIMULATION) then

    stop 'poroelastic simulation not implemented yet'
    ! but would be something like this...
    call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,rmass_solid_poroelastic, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

    call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,rmass_fluid_poroelastic, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

    ! fills mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmass_solid_poroelastic <= 0._CUSTOM_REAL) rmass_solid_poroelastic = 1._CUSTOM_REAL
    where(rmass_fluid_poroelastic <= 0._CUSTOM_REAL) rmass_fluid_poroelastic = 1._CUSTOM_REAL
    rmass_solid_poroelastic(:) = 1._CUSTOM_REAL / rmass_solid_poroelastic(:)
    rmass_fluid_poroelastic(:) = 1._CUSTOM_REAL / rmass_fluid_poroelastic(:)

  endif ! POROELASTIC_SIMULATION

  if(myrank == 0) write(IMAIN,*) 'end assembling MPI mass matrix'


  end subroutine prepare_timerun_mass_matrices

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
  double precision, dimension(N_SLS) :: tau_sigma_dble,beta_dble
  double precision factor_scale_dble,one_minus_sum_beta_dble
  double precision :: Q_mu
  double precision :: f_c_source
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tauinv
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: beta
  real(kind=CUSTOM_REAL):: scale_factor
  real(kind=CUSTOM_REAL):: vs_val
  integer :: i,j,k,ispec
  integer :: iselected ! iattenuation


  ! if attenuation is on, shift shear moduli to right frequency
  ! rescale mu to average (central) frequency for attenuation
  if(ATTENUATION) then

    ! initializes arrays
    one_minus_sum_beta(:,:,:,:) = 1._CUSTOM_REAL
    factor_common(:,:,:,:,:) = 1._CUSTOM_REAL

    ! precalculates tau_sigma depending on period band (constant for all Q_mu)
    call get_attenuation_constants(myrank,min_resolved_period,&
              tau_sigma_dble,f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

    ! determines alphaval,betaval,gammaval for runge-kutta
    if(CUSTOM_REAL == SIZE_REAL) then
      tau_sigma(:) = sngl(tau_sigma_dble(:))
    else
      tau_sigma(:) = tau_sigma_dble(:)
    endif
    call get_attenuation_memory_values(tau_sigma,deltat,alphaval,betaval,gammaval)

    ! precalculates the negative inverse of tau_sigma
    tauinv(:) = - 1._CUSTOM_REAL / tau_sigma(:)

    ! rescale shear modulus according to attenuation model
    do ispec = 1,NSPEC_AB

      ! skips non elastic elements
      if( ispec_is_elastic(ispec) .eqv. .false. ) cycle

      ! determines attenuation factors for each GLL point
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

            ! use scaling rule similar to Olsen et al. (2003)
            !! We might need to fix the attenuation part for the anisotropy case
            !! At this stage, we turn the ATTENUATION flag off always, and still keep mustore
            if(USE_OLSEN_ATTENUATION) then
              vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
              call get_attenuation_model_olsen( vs_val, Q_mu )
            else
              ! takes iflag set in (CUBIT) mesh
              iselected = iflag_attenuation_store(i,j,k,ispec)

              ! attenuation zero flag
              if( iselected == 0 ) cycle

              ! determines Q_mu
              select case(iselected)
              case(IATTENUATION_SEDIMENTS_40)
                Q_mu = 40.000000d0
              case(IATTENUATION_SEDIMENTS_50)
                Q_mu = 50.000000d0
              case(IATTENUATION_SEDIMENTS_60)
                Q_mu = 60.000000d0
              case(IATTENUATION_SEDIMENTS_70)
                Q_mu = 70.000000d0
              case(IATTENUATION_SEDIMENTS_80)
                Q_mu = 80.000000d0
              case(IATTENUATION_SEDIMENTS_90)
                Q_mu = 90.000000d0
              case(IATTENUATION_SEDIMENTS_100)
                Q_mu = 100.000000d0
              case(IATTENUATION_SEDIMENTS_110)
                Q_mu = 110.000000d0
              case(IATTENUATION_SEDIMENTS_120)
                Q_mu = 120.000000d0
              case(IATTENUATION_SEDIMENTS_130)
                Q_mu = 130.000000d0
              case(IATTENUATION_SEDIMENTS_140)
                Q_mu = 140.000000d0
              case(IATTENUATION_SEDIMENTS_150)
                Q_mu = 150.000000d0
              case(IATTENUATION_BEDROCK)
                Q_mu = 9000.d0
              case default
                call exit_MPI(myrank,'wrong attenuation flag in mesh')
              end select

            endif

            ! gets beta, on_minus_sum_beta and factor_scale
            call get_attenuation_factors(myrank,Q_mu,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                              f_c_source,tau_sigma_dble, &
                              beta_dble,one_minus_sum_beta_dble,factor_scale_dble)

            ! stores factor for unrelaxed parameter
            one_minus_sum_beta(i,j,k,ispec) = one_minus_sum_beta_dble

            ! stores scale factor for runge-kutta
            beta(:) = beta_dble(:)
            factor_common(:,i,j,k,ispec) = 2._CUSTOM_REAL * beta(:) * tauinv(:)

            ! scales only mu moduli
            scale_factor = factor_scale_dble
            mustore(i,j,k,ispec) = mustore(i,j,k,ispec) * scale_factor

          enddo
        enddo
      enddo
    enddo

    ! clear memory variables if attenuation
    ! initialize memory variables for attenuation
    epsilondev_xx(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_yy(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_xy(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_xz(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_yz(:,:,:,:) = 0._CUSTOM_REAL

    R_xx(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yy(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xy(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xz(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yz(:,:,:,:,:) = 0._CUSTOM_REAL

    if(FIX_UNDERFLOW_PROBLEM) then
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

  subroutine prepare_timerun_adjoint()

! prepares adjoint simulations

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  implicit none

  integer :: ier

! seismograms
  if (nrec_local > 0 .and. SIMULATION_TYPE == 2 ) then
    ! allocate Frechet derivatives array
    allocate(Mxx_der(nrec_local),Myy_der(nrec_local), &
            Mzz_der(nrec_local),Mxy_der(nrec_local), &
            Mxz_der(nrec_local),Myz_der(nrec_local), &
            sloc_der(NDIM,nrec_local))
    Mxx_der = 0._CUSTOM_REAL
    Myy_der = 0._CUSTOM_REAL
    Mzz_der = 0._CUSTOM_REAL
    Mxy_der = 0._CUSTOM_REAL
    Mxz_der = 0._CUSTOM_REAL
    Myz_der = 0._CUSTOM_REAL
    sloc_der = 0._CUSTOM_REAL

    allocate(seismograms_eps(NDIM,NDIM,nrec_local,NSTEP))
    seismograms_eps(:,:,:,:) = 0._CUSTOM_REAL
  endif

! timing
  if (SIMULATION_TYPE == 3) then

    ! backward/reconstructed wavefields: time stepping is in time-reversed sense
    ! (negative time increments)
    if(CUSTOM_REAL == SIZE_REAL) then
      b_deltat = - sngl(DT)
    else
      b_deltat = - DT
    endif
    b_deltatover2 = b_deltat/2._CUSTOM_REAL
    b_deltatsqover2 = b_deltat*b_deltat/2._CUSTOM_REAL

  endif

! attenuation backward memories
  if( ATTENUATION .and. SIMULATION_TYPE == 3 ) then

    ! precompute Runge-Kutta coefficients if attenuation
    call get_attenuation_memory_values(tau_sigma,b_deltat,b_alphaval,b_betaval,b_gammaval)

  endif

! initializes adjoint kernels and reconstructed/backward wavefields
  if (SIMULATION_TYPE == 3)  then
    ! elastic domain
    if( ELASTIC_SIMULATION ) then
      rho_kl(:,:,:,:)   = 0._CUSTOM_REAL
      mu_kl(:,:,:,:)    = 0._CUSTOM_REAL
      kappa_kl(:,:,:,:) = 0._CUSTOM_REAL

      ! reconstructed/backward elastic wavefields
      b_displ = 0._CUSTOM_REAL
      b_veloc = 0._CUSTOM_REAL
      b_accel = 0._CUSTOM_REAL
      if(FIX_UNDERFLOW_PROBLEM) b_displ = VERYSMALLVAL

      ! memory variables if attenuation
      if( ATTENUATION ) then
         b_R_xx = 0._CUSTOM_REAL
         b_R_yy = 0._CUSTOM_REAL
         b_R_xy = 0._CUSTOM_REAL
         b_R_xz = 0._CUSTOM_REAL
         b_R_yz = 0._CUSTOM_REAL
         b_epsilondev_xx = 0._CUSTOM_REAL
         b_epsilondev_yy = 0._CUSTOM_REAL
         b_epsilondev_xy = 0._CUSTOM_REAL
         b_epsilondev_xz = 0._CUSTOM_REAL
         b_epsilondev_yz = 0._CUSTOM_REAL
      endif

    endif

    ! acoustic domain
    if( ACOUSTIC_SIMULATION ) then
      rho_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
      kappa_ac_kl(:,:,:,:) = 0._CUSTOM_REAL

      ! reconstructed/backward acoustic potentials
      b_potential_acoustic = 0._CUSTOM_REAL
      b_potential_dot_acoustic = 0._CUSTOM_REAL
      b_potential_dot_dot_acoustic = 0._CUSTOM_REAL
      if(FIX_UNDERFLOW_PROBLEM) b_potential_acoustic = VERYSMALLVAL

    endif
  endif

! initialize Moho boundary index
  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    ispec2D_moho_top = 0
    ispec2D_moho_bot = 0
  endif

! stacey absorbing fields will be reconstructed for adjoint simulations
! using snapshot files of wavefields
  if( ABSORBING_CONDITIONS ) then

    ! opens absorbing wavefield saved/to-be-saved by forward simulations
    if( num_abs_boundary_faces > 0 .and. (SIMULATION_TYPE == 3 .or. &
          (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) ) then

      b_num_abs_boundary_faces = num_abs_boundary_faces

      ! elastic domains
      if( ELASTIC_SIMULATION) then
        ! allocates wavefield
        allocate(b_absorb_field(NDIM,NGLLSQUARE,b_num_abs_boundary_faces))

        b_reclen_field = CUSTOM_REAL * (NDIM * NGLLSQUARE * num_abs_boundary_faces)

        if (SIMULATION_TYPE == 3) then
          ! opens existing files
          open(unit=IOABS,file=trim(prname)//'absorb_field.bin',status='old',&
                action='read',form='unformatted',access='direct', &
                recl=b_reclen_field+2*4,iostat=ier )
          if( ier /= 0 ) call exit_mpi(myrank,'error opening proc***_absorb_field.bin file')
        else
          ! opens new file
          open(unit=IOABS,file=trim(prname)//'absorb_field.bin',status='unknown',&
                form='unformatted',access='direct',&
                recl=b_reclen_field+2*4,iostat=ier )
          if( ier /= 0 ) call exit_mpi(myrank,'error opening proc***_absorb_field.bin file')
        endif
      endif

      ! acoustic domains
      if( ACOUSTIC_SIMULATION) then
        ! allocates wavefield
        allocate(b_absorb_potential(NGLLSQUARE,b_num_abs_boundary_faces))

        b_reclen_potential = CUSTOM_REAL * (NGLLSQUARE * num_abs_boundary_faces)

        if (SIMULATION_TYPE == 3) then
          ! opens existing files
          open(unit=IOABS_AC,file=trim(prname)//'absorb_potential.bin',status='old',&
                action='read',form='unformatted',access='direct', &
                recl=b_reclen_potential+2*4,iostat=ier )
          if( ier /= 0 ) call exit_mpi(myrank,'error opening proc***_absorb_potential.bin file')
        else
          ! opens new file
          open(unit=IOABS_AC,file=trim(prname)//'absorb_potential.bin',status='unknown',&
                form='unformatted',access='direct',&
                recl=b_reclen_potential+2*4,iostat=ier )
          if( ier /= 0 ) call exit_mpi(myrank,'error opening proc***_absorb_potential.bin file')
        endif
      endif

    else
      ! dummy array
      b_num_abs_boundary_faces = 1
      if( ELASTIC_SIMULATION ) &
        allocate(b_absorb_field(NDIM,NGLLSQUARE,b_num_abs_boundary_faces))

      if( ACOUSTIC_SIMULATION ) &
        allocate(b_absorb_potential(NGLLSQUARE,b_num_abs_boundary_faces))

    endif
  endif

  end subroutine prepare_timerun_adjoint
