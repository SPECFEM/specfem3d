!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

  subroutine iterate_time()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  implicit none

!
!   s t a r t   t i m e   i t e r a t i o n s
!

! synchronize all processes to make sure everybody is ready to start time loop
  call sync_all()
  if(myrank == 0) write(IMAIN,*) 'All processes are synchronized before time loop'

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop...'
    write(IMAIN,*)
  endif

! create an empty file to monitor the start of the simulation
  if(myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown')
    write(IOUT,*) 'starting time loop'
    close(IOUT)
  endif

! get MPI starting time
  time_start = wtime()

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

  do it = 1,NSTEP

    ! simulation status output and stability check
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      call it_check_stability()
    endif

    ! update displacement using Newmark time scheme
    call it_update_displacement_scheme()

    ! acoustic solver
    ! (needs to be done first, before elastic one)
    if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic()

    ! elastic solver
    ! (needs to be done first, before poroelastic one)
    if( ELASTIC_SIMULATION ) call compute_forces_elastic()

    ! poroelastic solver
    if( POROELASTIC_SIMULATION ) call compute_forces_poroelastic()

    ! restores last time snapshot saved for backward/reconstruction of wavefields
    ! note: this must be read in after the Newmark time scheme
    if( SIMULATION_TYPE == 3 .and. it == 1 ) then
     call it_read_forward_arrays()
    endif

    ! write the seismograms with time shift (GPU_MODE transfer included)
    if (nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0 ) ) then
      call write_seismograms()
    endif

    ! resetting d/v/a/R/eps for the backward reconstruction with attenuation
    if (ATTENUATION ) then
      call it_store_attenuation_arrays()
    endif

    ! adjoint simulations: kernels
    if( SIMULATION_TYPE == 3 ) then
      call compute_kernels()
    endif

    ! outputs movie files
    if( MOVIE_SIMULATION ) then
      call write_movie_output()
    endif

    ! first step of noise tomography, i.e., save a surface movie at every time step
    if ( NOISE_TOMOGRAPHY == 1) then
      if( num_free_surface_faces > 0) then
        call noise_save_surface_movie(displ,ibool, &
                            noise_surface_movie,it, &
                            NSPEC_AB,NGLOB_AB, &
                            num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                            Mesh_pointer,GPU_MODE)
      endif
    endif

!
!---- end of time iteration loop
!
  enddo   ! end of main time loop

  call it_print_elapsed_time()

  ! Transfer fields from GPU card to host for further analysis
  if(GPU_MODE) call it_transfer_from_GPU()

  end subroutine iterate_time


!=====================================================================

  subroutine it_check_stability()

! computes the maximum of the norm of the displacement
! in all the slices using an MPI reduction
! and output timestamp file to check that simulation is running fine

  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic
  implicit none

  double precision :: tCPU,t_remain,t_total
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total

  ! maximum of the norm of the displacement
  real(kind=CUSTOM_REAL) Usolidnorm,Usolidnorm_all ! elastic
  real(kind=CUSTOM_REAL) Usolidnormp,Usolidnormp_all ! acoustic
  real(kind=CUSTOM_REAL) Usolidnorms,Usolidnorms_all ! solid poroelastic
  real(kind=CUSTOM_REAL) Usolidnormw,Usolidnormw_all ! fluid (w.r.t.s) poroelastic

  ! norm of the backward displacement
  real(kind=CUSTOM_REAL) b_Usolidnorm, b_Usolidnorm_all

  ! initializes
  Usolidnorm_all = 0.0_CUSTOM_REAL
  Usolidnormp_all = 0.0_CUSTOM_REAL
  Usolidnorms_all = 0.0_CUSTOM_REAL
  Usolidnormw_all = 0.0_CUSTOM_REAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!chris: Rewrite to get norm for each material when coupled simulations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! compute maximum of norm of displacement in each slice
  if( ELASTIC_SIMULATION ) then
    if( GPU_MODE) then
      ! way 2: just get maximum of field from GPU
      call get_norm_elastic_from_device(Usolidnorm,Mesh_pointer,1)
    else
      Usolidnorm = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2 + displ(3,:)**2))
    endif

    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    !if(Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0.0_CUSTOM_REAL) &
    !  call exit_MPI(myrank,'single forward simulation became unstable and blew up')

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnorm,Usolidnorm_all)
  endif

  if( ACOUSTIC_SIMULATION ) then
    if(GPU_MODE) then
      ! way 2: just get maximum of field from GPU
      call get_norm_acoustic_from_device(Usolidnormp,Mesh_pointer,1)
    else
      Usolidnormp = maxval(abs(potential_dot_dot_acoustic(:)))
    endif

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnormp,Usolidnormp_all)
  endif

  if( POROELASTIC_SIMULATION ) then
    Usolidnorms = maxval(sqrt(displs_poroelastic(1,:)**2 + displs_poroelastic(2,:)**2 + &
                             displs_poroelastic(3,:)**2))
    Usolidnormw = maxval(sqrt(displw_poroelastic(1,:)**2 + displw_poroelastic(2,:)**2 + &
                             displw_poroelastic(3,:)**2))

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnorms,Usolidnorms_all)
    call max_all_cr(Usolidnormw,Usolidnormw_all)
  endif


  ! adjoint simulations
  if( SIMULATION_TYPE == 3 ) then
    if( ELASTIC_SIMULATION ) then
      ! way 2
      if(GPU_MODE) then
        call get_norm_elastic_from_device(b_Usolidnorm,Mesh_pointer,3)
      else
        b_Usolidnorm = maxval(sqrt(b_displ(1,:)**2 + b_displ(2,:)**2 + b_displ(3,:)**2))
      endif
    endif
    if( ACOUSTIC_SIMULATION ) then
        ! way 2
        if(GPU_MODE) then
          call get_norm_acoustic_from_device(b_Usolidnorm,Mesh_pointer,3)
        else
          b_Usolidnorm = maxval(abs(b_potential_dot_dot_acoustic(:)))
        endif
    endif
    if( POROELASTIC_SIMULATION ) then
        b_Usolidnorm = maxval(sqrt(b_displs_poroelastic(1,:)**2 + b_displs_poroelastic(2,:)**2 + &
                                 b_displs_poroelastic(3,:)**2))
    endif
    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    !if(b_Usolidnorm > STABILITY_THRESHOLD .or. b_Usolidnorm < 0.0_CUSTOM_REAL) &
    !  call exit_MPI(myrank,'single backward simulation became unstable and blew up')

    ! compute max of all slices
    call max_all_cr(b_Usolidnorm,b_Usolidnorm_all)
  endif

  ! user output
  if(myrank == 0) then

    write(IMAIN,*) 'Time step # ',it
    write(IMAIN,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'

    ! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',sngl(tCPU/dble(it))
    if( ELASTIC_SIMULATION ) then
      write(IMAIN,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all
    endif
    if( ACOUSTIC_SIMULATION ) then
        write(IMAIN,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all
    endif
    if( POROELASTIC_SIMULATION ) then
        write(IMAIN,*) 'Max norm displacement vector Us in all slices (m) = ',Usolidnorms_all
        write(IMAIN,*) 'Max norm displacement vector W in all slices (m) = ',Usolidnormw_all
    endif
    ! adjoint simulations
    if (SIMULATION_TYPE == 3) write(IMAIN,*) &
           'Max norm U (backward) in all slices = ',b_Usolidnorm_all

    ! compute estimated remaining simulation time
    t_remain = (NSTEP - it) * (tCPU/dble(it))
    int_t_remain = int(t_remain)
    ihours_remain = int_t_remain / 3600
    iminutes_remain = (int_t_remain - 3600*ihours_remain) / 60
    iseconds_remain = int_t_remain - 3600*ihours_remain - 60*iminutes_remain
    write(IMAIN,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IMAIN,*) 'Time steps remaining = ',NSTEP - it
    write(IMAIN,*) 'Estimated remaining time in seconds = ',sngl(t_remain)
    write(IMAIN,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain

    ! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
    write(IMAIN,*) 'Estimated total run time in seconds = ',sngl(t_total)
    write(IMAIN,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'

    if(it < 100) then
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '**** BEWARE: the above time estimates are not reliable'
      write(IMAIN,*) '**** because fewer than 100 iterations have been performed'
      write(IMAIN,*) '************************************************************'
    endif
    write(IMAIN,*)

    ! flushes file buffer for main output file (IMAIN)
    call flush_IMAIN()

    ! write time stamp file to give information about progression of simulation
    write(outputname,"('/timestamp',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown')
    write(IOUT,*) 'Time step # ',it
    write(IOUT,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'
    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    if( ELASTIC_SIMULATION ) then
      write(IOUT,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all
    endif
    if( ACOUSTIC_SIMULATION ) then
        write(IOUT,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all
    endif
    if( POROELASTIC_SIMULATION ) then
        write(IOUT,*) 'Max norm displacement vector Us in all slices (m) = ',Usolidnorms_all
        write(IOUT,*) 'Max norm displacement vector W in all slices (m) = ',Usolidnormw_all
    endif
    ! adjoint simulations
    if (SIMULATION_TYPE == 3) write(IOUT,*) &
           'Max norm U (backward) in all slices = ',b_Usolidnorm_all
    ! estimation
    write(IOUT,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IOUT,*) 'Time steps remaining = ',NSTEP - it
    write(IOUT,*) 'Estimated remaining time in seconds = ',t_remain
    write(IOUT,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain
    write(IOUT,*) 'Estimated total run time in seconds = ',t_total
    write(IOUT,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IOUT,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'
    close(IOUT)

    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if(Usolidnorm_all > STABILITY_THRESHOLD .or. Usolidnorm_all < 0.0_CUSTOM_REAL &
     .or. Usolidnormp_all > STABILITY_THRESHOLD .or. Usolidnormp_all < 0.0_CUSTOM_REAL &
     .or. Usolidnorms_all > STABILITY_THRESHOLD .or. Usolidnorms_all < 0.0_CUSTOM_REAL &
     .or. Usolidnormw_all > STABILITY_THRESHOLD .or. Usolidnormw_all < 0.0_CUSTOM_REAL) &
        call exit_MPI(myrank,'forward simulation became unstable and blew up')
    ! adjoint simulations
    if(SIMULATION_TYPE == 3 .and. (b_Usolidnorm_all > STABILITY_THRESHOLD &
      .or. b_Usolidnorm_all < 0.0)) &
        call exit_MPI(myrank,'backward simulation became unstable and blew up')

  endif ! myrank

  end subroutine it_check_stability


!=====================================================================

  subroutine it_update_displacement_scheme()

! explicit Newmark time scheme with acoustic & elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 delta_t chi_dot_dot(t+delta_t)
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!
! note that this stage calculates the predictor terms
!
!   for
!   potential chi_dot(t+delta) requires + 1/2 delta_t chi_dot_dot(t+delta_t)
!                                   at a later stage (corrector) once where chi_dot_dot(t+delta) is calculated
!   and similar,
!   velocity v(t+delta_t) requires  + 1/2 delta_t a(t+delta_t)
!                                   at a later stage once where a(t+delta) is calculated
! also:
!   boundary term B_elastic requires chi_dot_dot(t+delta)
!                                   thus chi_dot_dot has to be updated first before the elastic boundary term is considered

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use PML_par
  use PML_par_acoustic
  implicit none

! updates acoustic potentials
  if( ACOUSTIC_SIMULATION ) then

    if(.NOT. GPU_MODE) then
      ! on CPU
      potential_acoustic(:) = potential_acoustic(:) &
                            + deltat * potential_dot_acoustic(:) &
                            + deltatsqover2 * potential_dot_dot_acoustic(:)
      potential_dot_acoustic(:) = potential_dot_acoustic(:) &
                                + deltatover2 * potential_dot_dot_acoustic(:)
      potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
    else
      ! on GPU
      call it_update_displacement_ac_cuda(Mesh_pointer, NGLOB_AB, &
                                          deltat, deltatsqover2, deltatover2, &
                                          b_deltat, b_deltatsqover2, b_deltatover2)
    endif

    ! time marching potentials
    if(ABSORB_USE_PML .and. ABSORBING_CONDITIONS) then
      if( GPU_MODE ) call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                              potential_dot_acoustic, potential_dot_dot_acoustic, Mesh_pointer)

      call PML_acoustic_time_march(NSPEC_AB,NGLOB_AB,ibool,&
                        potential_acoustic,potential_dot_acoustic,&
                        deltat,deltatsqover2,deltatover2,&
                        num_PML_ispec,PML_ispec,PML_damping_d,&
                        chi1,chi2,chi2_t,chi3,chi4,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot,&
                        iglob_is_PML_interface,PML_mask_ibool,&
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,&
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh,NPROC,&
                        ispec_is_acoustic)

      if( GPU_MODE ) call transfer_fields_ac_to_device(NGLOB_AB,potential_acoustic, &
                              potential_dot_acoustic, potential_dot_dot_acoustic, Mesh_pointer)
    endif

  endif ! ACOUSTIC_SIMULATION

! updates elastic displacement and velocity
  if( ELASTIC_SIMULATION ) then

    if(.NOT. GPU_MODE) then
      ! on CPU
      displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
      veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
      if( SIMULATION_TYPE /= 1 ) accel_adj_coupling(:,:) = accel(:,:)
      accel(:,:) = 0._CUSTOM_REAL
    else
      ! on GPU
      ! Includes SIM_TYPE 1 & 3 (for noise tomography)
      call it_update_displacement_cuda(Mesh_pointer, size(displ), deltat, deltatsqover2,&
                                       deltatover2, b_deltat, b_deltatsqover2, b_deltatover2)
    endif
  endif

! updates poroelastic displacements and velocities
  if( POROELASTIC_SIMULATION ) then
    ! solid phase
    displs_poroelastic(:,:) = displs_poroelastic(:,:) + deltat*velocs_poroelastic(:,:) + &
                              deltatsqover2*accels_poroelastic(:,:)
    velocs_poroelastic(:,:) = velocs_poroelastic(:,:) + deltatover2*accels_poroelastic(:,:)
    accels_poroelastic(:,:) = 0._CUSTOM_REAL

    ! fluid phase
    displw_poroelastic(:,:) = displw_poroelastic(:,:) + deltat*velocw_poroelastic(:,:) + &
                              deltatsqover2*accelw_poroelastic(:,:)
    velocw_poroelastic(:,:) = velocw_poroelastic(:,:) + deltatover2*accelw_poroelastic(:,:)
    accelw_poroelastic(:,:) = 0._CUSTOM_REAL
  endif

! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. .NOT. GPU_MODE) then
    ! acoustic backward fields
    if( ACOUSTIC_SIMULATION ) then
      b_potential_acoustic(:) = b_potential_acoustic(:) &
                              + b_deltat * b_potential_dot_acoustic(:) &
                              + b_deltatsqover2 * b_potential_dot_dot_acoustic(:)
      b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) &
                                  + b_deltatover2 * b_potential_dot_dot_acoustic(:)
      b_potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
    endif

    ! elastic backward fields
    if( ELASTIC_SIMULATION ) then
      b_displ(:,:) = b_displ(:,:) + b_deltat*b_veloc(:,:) + b_deltatsqover2*b_accel(:,:)
      b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)
      b_accel(:,:) = 0._CUSTOM_REAL
    endif
    ! poroelastic backward fields
    if( POROELASTIC_SIMULATION ) then
    ! solid phase
    b_displs_poroelastic(:,:) = b_displs_poroelastic(:,:) + b_deltat*b_velocs_poroelastic(:,:) + &
                              b_deltatsqover2*b_accels_poroelastic(:,:)
    b_velocs_poroelastic(:,:) = b_velocs_poroelastic(:,:) + b_deltatover2*b_accels_poroelastic(:,:)
    b_accels_poroelastic(:,:) = 0._CUSTOM_REAL

    ! fluid phase
    b_displw_poroelastic(:,:) = b_displw_poroelastic(:,:) + b_deltat*b_velocw_poroelastic(:,:) + &
                              b_deltatsqover2*b_accelw_poroelastic(:,:)
    b_velocw_poroelastic(:,:) = b_velocw_poroelastic(:,:) + b_deltatover2*b_accelw_poroelastic(:,:)
    b_accelw_poroelastic(:,:) = 0._CUSTOM_REAL
    endif
  endif

! adjoint simulations: moho kernel
  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    ispec2D_moho_top = 0
    ispec2D_moho_bot = 0
  endif


  end subroutine it_update_displacement_scheme

!=====================================================================

  subroutine it_read_forward_arrays()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  implicit none

  integer :: ier

! restores last time snapshot saved for backward/reconstruction of wavefields
! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
!          and adjoint sources will become more complicated
!          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields

  ! reads in wavefields
  open(unit=27,file=trim(prname)//'save_forward_arrays.bin',status='old',&
        action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error: opening save_forward_arrays'
    print*,'path: ',trim(prname)//'save_forward_arrays.bin'
    call exit_mpi(myrank,'error open file save_forward_arrays.bin')
  endif

  if( ACOUSTIC_SIMULATION ) then
    read(27) b_potential_acoustic
    read(27) b_potential_dot_acoustic
    read(27) b_potential_dot_dot_acoustic

    ! transfers fields onto GPU
    if(GPU_MODE) &
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                          b_potential_dot_acoustic, b_potential_dot_dot_acoustic, Mesh_pointer)
  endif

  ! elastic wavefields
  if( ELASTIC_SIMULATION ) then
    read(27) b_displ
    read(27) b_veloc
    read(27) b_accel

    ! puts elastic wavefield to GPU
    if(GPU_MODE) &
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)

    ! memory variables if attenuation
    if( ATTENUATION ) then
      read(27) b_R_xx
      read(27) b_R_yy
      read(27) b_R_xy
      read(27) b_R_xz
      read(27) b_R_yz
      read(27) b_epsilondev_xx
      read(27) b_epsilondev_yy
      read(27) b_epsilondev_xy
      read(27) b_epsilondev_xz
      read(27) b_epsilondev_yz

      ! puts elastic attenuation arrays to GPU
      if(GPU_MODE) &
          call transfer_b_fields_att_to_device(Mesh_pointer, &
                                            b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,size(b_R_xx), &
                                            b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz, &
                                            size(b_epsilondev_xx))
    endif

  endif

  ! poroelastic wavefields
  if( POROELASTIC_SIMULATION ) then
    read(27) b_displs_poroelastic
    read(27) b_velocs_poroelastic
    read(27) b_accels_poroelastic
    read(27) b_displw_poroelastic
    read(27) b_velocw_poroelastic
    read(27) b_accelw_poroelastic
  endif

  close(27)

  end subroutine it_read_forward_arrays

!=====================================================================

  subroutine it_store_attenuation_arrays()

! resetting d/v/a/R/eps for the backward reconstruction with attenuation

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  if( it > 1 .and. it < NSTEP) then
    ! adjoint simulations

! note backward/reconstructed wavefields:
!       storing wavefield displ() at time step it, corresponds to time (it-1)*DT - t0 (see routine write_seismograms_to_file )
!       reconstucted wavefield b_displ() at it corresponds to time (NSTEP-it-1)*DT - t0
!       we read in the reconstructed wavefield at the end of the time iteration loop, i.e. after the Newmark scheme,
!       thus, indexing is NSTEP-it (rather than something like NSTEP-(it-1) )
    if (SIMULATION_TYPE == 3 .and. mod(NSTEP-it,NSTEP_Q_SAVE) == 0) then
      ! reads files content
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") NSTEP-it
      open(unit=27,file=trim(prname_Q)//trim(outputname),status='old',&
            action='read',form='unformatted')
      if( ELASTIC_SIMULATION ) then
        ! reads arrays from disk files
        read(27) b_displ
        read(27) b_veloc
        read(27) b_accel

        ! puts elastic fields onto GPU
        if(GPU_MODE) then
          ! wavefields
          call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel, Mesh_pointer)
        endif

        read(27) b_R_xx
        read(27) b_R_yy
        read(27) b_R_xy
        read(27) b_R_xz
        read(27) b_R_yz
        read(27) b_epsilondev_xx
        read(27) b_epsilondev_yy
        read(27) b_epsilondev_xy
        read(27) b_epsilondev_xz
        read(27) b_epsilondev_yz

        ! puts elastic fields onto GPU
        if(GPU_MODE) then
          ! attenuation arrays
          call transfer_b_fields_att_to_device(Mesh_pointer, &
                                            b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,size(b_R_xx), &
                                            b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz, &
                                            size(b_epsilondev_xx))
        endif
      endif

      if( ACOUSTIC_SIMULATION ) then
        ! reads arrays from disk files
        read(27) b_potential_acoustic
        read(27) b_potential_dot_acoustic
        read(27) b_potential_dot_dot_acoustic

        ! puts acoustic fields onto GPU
        if(GPU_MODE) &
          call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                              b_potential_dot_acoustic, b_potential_dot_dot_acoustic, Mesh_pointer)

      endif
      close(27)
    else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. mod(it,NSTEP_Q_SAVE) == 0) then
      ! stores files content
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") it
      open(unit=27,file=trim(prname_Q)//trim(outputname),status='unknown',&
           action='write',form='unformatted')
      if( ELASTIC_SIMULATION ) then
        ! gets elastic fields from GPU onto CPU
        if(GPU_MODE) then
          call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc, accel, Mesh_pointer)
        endif

        ! writes to disk file
        write(27) displ
        write(27) veloc
        write(27) accel

        if(GPU_MODE) then
          ! attenuation arrays
          call transfer_fields_att_from_device(Mesh_pointer, &
                                            R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                                            epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                            size(epsilondev_xx))
        endif

        write(27) R_xx
        write(27) R_yy
        write(27) R_xy
        write(27) R_xz
        write(27) R_yz
        write(27) epsilondev_xx
        write(27) epsilondev_yy
        write(27) epsilondev_xy
        write(27) epsilondev_xz
        write(27) epsilondev_yz
      endif
      if( ACOUSTIC_SIMULATION ) then
       ! gets acoustic fields from GPU onto CPU
        if(GPU_MODE) &
          call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                              potential_dot_acoustic, potential_dot_dot_acoustic, Mesh_pointer)

        ! writes to disk file
        write(27) potential_acoustic
        write(27) potential_dot_acoustic
        write(27) potential_dot_dot_acoustic
      endif
      close(27)
    endif ! SIMULATION_TYPE
  endif ! it

  end subroutine it_store_attenuation_arrays

!=====================================================================

  subroutine it_print_elapsed_time()

    use specfem_par
    use specfem_par_elastic
    use specfem_par_acoustic
    implicit none

    ! local parameters
    double precision :: tCPU
    integer :: ihours,iminutes,iseconds,int_tCPU

    if(myrank == 0) then
       ! elapsed time since beginning of the simulation
       tCPU = wtime() - time_start
       int_tCPU = int(tCPU)
       ihours = int_tCPU / 3600
       iminutes = (int_tCPU - 3600*ihours) / 60
       iseconds = int_tCPU - 3600*ihours - 60*iminutes
       write(IMAIN,*) 'Time-Loop Complete. Timing info:'
       write(IMAIN,*) 'Total elapsed time in seconds = ',tCPU
       write(IMAIN,"(' Total elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    endif

  end subroutine it_print_elapsed_time

!=====================================================================

  subroutine it_transfer_from_GPU()

! transfers fields on GPU back onto CPU

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  ! to store forward wave fields
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then

    ! acoustic potentials
    if( ACOUSTIC_SIMULATION ) &
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                            potential_dot_acoustic, potential_dot_dot_acoustic, Mesh_pointer)

    ! elastic wavefield
    if( ELASTIC_SIMULATION ) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc, accel, Mesh_pointer)

      if (ATTENUATION) &
        call transfer_fields_att_from_device(Mesh_pointer, &
                                            R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                                            epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                            size(epsilondev_xx))

    endif
  else if (SIMULATION_TYPE == 3) then

    ! to store kernels
    ! acoustic domains
    if( ACOUSTIC_SIMULATION ) then
      ! only in case needed...
      !call transfer_b_fields_ac_from_device(NGLOB_AB,b_potential_acoustic, &
      !                      b_potential_dot_acoustic, b_potential_dot_dot_acoustic, Mesh_pointer)

      ! acoustic kernels
      call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
    endif

    ! elastic domains
    if( ELASTIC_SIMULATION ) then
      ! only in case needed...
      !call transfer_b_fields_from_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel, Mesh_pointer)

      ! elastic kernels
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,NSPEC_AB)
    endif

    ! specific noise strength kernel
    if( NOISE_TOMOGRAPHY == 3 ) then
      call transfer_kernels_noise_to_host(Mesh_pointer,Sigma_kl,NSPEC_AB)
    endif

    ! approximative hessian for preconditioning kernels
    if ( APPROXIMATE_HESS_KL ) then
      if( ELASTIC_SIMULATION ) call transfer_kernels_hess_el_tohost(Mesh_pointer,hess_kl,NSPEC_AB)
      if( ACOUSTIC_SIMULATION ) call transfer_kernels_hess_ac_tohost(Mesh_pointer,hess_ac_kl,NSPEC_AB)
    endif

  endif

  ! frees allocated memory on GPU
  call prepare_cleanup_device(Mesh_pointer, &
                              SAVE_FORWARD, &
                              ACOUSTIC_SIMULATION,ELASTIC_SIMULATION, &
                              ABSORBING_CONDITIONS,NOISE_TOMOGRAPHY,COMPUTE_AND_STORE_STRAIN, &
                              ATTENUATION,ANISOTROPY,OCEANS, &
                              APPROXIMATE_HESS_KL)

  end subroutine it_transfer_from_GPU
