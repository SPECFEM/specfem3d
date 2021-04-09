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


  subroutine iterate_time_undoatt()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use gravity_perturbation, only: gravity_timeseries, GRAVITY_SIMULATION

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_potential_acoustic_buffer,b_potential_dot_dot_acoustic_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_elastic_buffer,b_veloc_elastic_buffer,b_accel_elastic_buffer
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: b_noise_surface_movie_buffer

  double precision :: sizeval
  integer :: it_temp,seismo_current_temp,ier

  ! timing
  double precision, external :: wtime
#ifdef VTK_VIS
  logical :: do_restart = .false.
#endif

  ! checks if anything to do
  if (.not. UNDO_ATTENUATION_AND_OR_PML) return

  ! safety checks
  if (NOISE_TOMOGRAPHY /= 0) &
    call exit_MPI(myrank,'for undo_attenuation, NOISE_TOMOGRAPHY is not supported yet')
  if (SIMULATION_TYPE == 3 .and. POROELASTIC_SIMULATION) &
    call exit_MPI(myrank,'for undo_attenuation, POROELASTIC kernel simulation is not supported yet') ! not working yet...
  if (SIMULATION_TYPE == 3 .and. NOISE_TOMOGRAPHY /= 0) &
    call exit_MPI(myrank,'for undo_attenuation, noise kernel simulation is not supported yet') ! not working yet...

  ! note: NSTEP must not be a multiple of NT_DUMP_ATTENUATION, but should be equal or larger
  !
  !       also, we don't make sure if buffer size is not too big for system memory to hold in RAM.
  !       allocating too much memory will resort to memory paging into virtual memory and slow down performance.
  if (NT_DUMP_ATTENUATION > NSTEP) NT_DUMP_ATTENUATION = NSTEP

  ! number of time subsets for time loop
  NSUBSET_ITERATIONS = ceiling( dble(NSTEP)/dble(NT_DUMP_ATTENUATION) )

  ! checks
  if (NSUBSET_ITERATIONS <= 0) call exit_MPI(myrank,'Error invalid number of time subsets for undoing attenuation')

  ! user output
  if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) then
    if (myrank == 0) then
      write(IMAIN,*) 'undoing attenuation:'
      write(IMAIN,*) '  total number of time subsets                     = ',NSUBSET_ITERATIONS
      write(IMAIN,*) '  wavefield snapshots at every NT_DUMP_ATTENUATION = ',NT_DUMP_ATTENUATION
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! allocates buffers
  if (SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      if (ACOUSTIC_SIMULATION) then
        ! buffer(NGLOB_AB,NT_DUMP_ATTENUATION) in MB
        sizeval = 2.0 * dble(NGLOB_AB) * dble(NT_DUMP_ATTENUATION) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        write(IMAIN,*) '  size of acoustic wavefield buffer per slice      = ', sngl(sizeval),'MB'
      endif
      if (ELASTIC_SIMULATION) then
        ! buffer(3,NGLOB_AB,NT_DUMP_ATTENUATION) in MB
        sizeval = dble(NDIM) * dble(NGLOB_AB) * dble(NT_DUMP_ATTENUATION) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        if (APPROXIMATE_HESS_KL) sizeval = 2 * sizeval
        write(IMAIN,*) '  size of elastic wavefield buffer per slice       = ', sngl(sizeval),'MB'
      endif
    endif
    ! buffer arrays
    if (ACOUSTIC_SIMULATION) then
      allocate(b_potential_acoustic_buffer(NGLOB_AB,NT_DUMP_ATTENUATION), &
               b_potential_dot_dot_acoustic_buffer(NGLOB_AB,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_potential_acoustic arrays')
    endif
    if (ELASTIC_SIMULATION) then
      allocate(b_displ_elastic_buffer(NDIM,NGLOB_AB,NT_DUMP_ATTENUATION),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_elastic')
      if (APPROXIMATE_HESS_KL) then
        allocate(b_veloc_elastic_buffer(NDIM,NGLOB_AB,NT_DUMP_ATTENUATION), &
                 b_accel_elastic_buffer(NDIM,NGLOB_AB,NT_DUMP_ATTENUATION),stat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'error allocating b_accel_elastic arrays')
      endif
      ! noise kernel for source strength (sigma_kernel) needs buffer for reconstructed noise_surface_movie array,
      ! otherwise we need file i/o which will considerably slow down performance
      !if (NOISE_TOMOGRAPHY == 3) then
      !  allocate(b_noise_surface_movie_buffer(NDIM,NGLLX,NGLLY,NSPEC_TOP,NT_DUMP_ATTENUATION),stat=ier)
      !  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_noise_surface_movie_buffer')
      !endif
    endif

    ! for faster GPU mem copies
    if (GPU_MODE) then
      if (ACOUSTIC_SIMULATION) call register_host_array(NGLOB_AB,b_potential_acoustic)
      if (ELASTIC_SIMULATION) then
        call register_host_array(NDIM*NGLOB_AB,b_displ)
        if (APPROXIMATE_HESS_KL) call register_host_array(NDIM*NGLOB_AB,b_accel)
      endif
    endif
  endif

  !----  create a Gnuplot script to display the energy curve in log scale
  if (OUTPUT_ENERGY .and. myrank == 0) then
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set terminal x11'
    write(IOUT_ENERGY,*) '#set terminal wxt'
    write(IOUT_ENERGY,*) '#set terminal postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time step number"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(a152)') 'plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1, "energy.dat" us 1:3 &
                         &t "Potential Energy" w l lc 2, "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:3 t "Potential Energy" w l lc 2'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

 !! CD CD adds this (temporary)
  if (RECIPROCITY_AND_KH_INTEGRAL) open(unit=158,file='KH_integral',status='unknown')

  ! open the file in which we will store the energy curve
  if (OUTPUT_ENERGY .and. myrank == 0) &
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')


#ifdef VTK_VIS
  ! restart: goto starting point
123 continue
  if (VTK_MODE) then
    if (do_restart) then
      if (myrank == 0) print *,'VTK_VIS: restarting simulation'
    endif
  endif
#endif

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()
  if (myrank == 0) write(IMAIN,*) 'All processes are synchronized before the time loop'

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop in undoing attenuation...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! create an empty file to monitor the start of the simulation
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop'
    close(IOUT)
  endif

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  ! time loop increments begin/end
  it_begin = 1
  it_end = NSTEP

  ! initialize variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! initializes time increments
  it = 0

  ! get MPI starting
  time_start = wtime()

  ! loops over time subsets
  do iteration_on_subset = 1, NSUBSET_ITERATIONS

    ! wavefield storage
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! saves forward wavefields
      call save_forward_arrays_undoatt()

    else if (SIMULATION_TYPE == 3) then
      ! reads in last stored forward wavefield
      call read_forward_arrays_undoatt()

      ! note: after reading the restart files of displacement back from disk, recompute the strain from displacement;
      !       this is better than storing the strain to disk as well, which would drastically increase I/O volume
      ! computes strain based on current backward/reconstructed wavefield
      !if (COMPUTE_AND_STORE_STRAIN) call compute_strain_att_backward()   ! for future use to reduce file sizes...
    endif

    ! time loop within this iteration subset
    select case (SIMULATION_TYPE)
    case (1, 2)
      ! forward and adjoint simulations

      ! increment end of this subset
      if (iteration_on_subset < NSUBSET_ITERATIONS) then
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      else
        ! loops over remaining steps in last subset
        it_subset_end = NSTEP - (iteration_on_subset-1)*NT_DUMP_ATTENUATION
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability()
        endif

        ! simulation status output and stability check
        if (OUTPUT_ENERGY) then
          if (mod(it,NTSTEP_BETWEEN_OUTPUT_ENERGY) == 0 .or. it == 5 .or. it == NSTEP) call compute_energy()
        endif

        ! forward simulations
        do istage = 1, NSTAGE_TIME_SCHEME  ! is equal to 1 if Newmark because only one stage then
          if (USE_LDDRK) then
            ! updates wavefields using LDDRK time scheme
            call update_displ_lddrk()
          else
            ! updates wavefields using Newmark time scheme
            call update_displ_Newmark()
          endif

          ! calculates stiffness term
          ! acoustic solver (needs to be done first, before elastic one)
          if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_forward_calling()
          ! elastic solver
          if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_calling()
          ! poroelastic solver
          if (POROELASTIC_SIMULATION) call compute_forces_poroelastic_calling()
        enddo ! istage

        ! calculating gravity field at current timestep
        if (GRAVITY_SIMULATION) call gravity_timeseries()

        ! write the seismograms with time shift (GPU_MODE transfer included)
        call write_seismograms()

        ! outputs movie files
        if (MOVIE_SIMULATION) call write_movie_output()

        ! first step of noise tomography, i.e., save a surface movie at every time step
        ! modified from the subroutine 'write_movie_surface'
        if (NOISE_TOMOGRAPHY == 1) then
          call noise_save_surface_movie()
        endif

#ifdef VTK_VIS
        if (VTK_MODE) then
          ! updates VTK window
          call vtk_window_update()
        endif
#endif

        !! CD CD add this : under validation option
        if (RECIPROCITY_AND_KH_INTEGRAL) then
          if (.not. SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
            call surface_or_volume_integral_on_whole_domain()
            write(158,*) it*DT, integral_boun(1), integral_boun(2), integral_boun(3)
          endif
        endif

      enddo ! subset loop

    case (3)
      ! kernel simulations

      ! intermediate storage of it and seismo_current positions
      it_temp = it
      seismo_current_temp = seismo_current

      ! increment end of this subset
      if (iteration_on_subset == 1) then
        ! loops over remaining steps in last forward subset
        it_subset_end = NSTEP - (NSUBSET_ITERATIONS-1)*NT_DUMP_ATTENUATION
      else
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! reconstructs forward wavefields based on last stored wavefield data

      ! note: we step forward in time here, starting from last snapshot.
      !       the newly computed, reconstructed forward wavefields (b_displ_..) get stored in buffers.

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability_backward()
        endif

        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then
          if (USE_LDDRK) then
            ! update displacement using Runge-Kutta time scheme
            call update_displ_lddrk_backward()
          else
            ! update displacement using Newmark time scheme
            call update_displ_Newmark_backward()
          endif

          ! calculates stiffness term
          ! backward/reconstructed wavefields (propagates forward in time)
          ! acoustic solver (needs to be done first in forward direction)
          if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_backward_calling()
          ! elastic solver
          if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_backward_calling()
          ! poroelastic solver
          if (POROELASTIC_SIMULATION) stop 'POROELASTIC simulations in undo_att iterations not supported yet'
        enddo ! istage

        ! transfers wavefields from GPU to CPU for buffering
        if (GPU_MODE) then
          ! note these transfers are blocking at the moment, might be done async in future...
          if (ACOUSTIC_SIMULATION) then
            call transfer_b_potential_ac_from_device(NGLOB_AB,b_potential_acoustic,Mesh_pointer)
            call transfer_b_potential_dot_dot_ac_from_device(NGLOB_AB,b_potential_dot_dot_acoustic,Mesh_pointer)
          endif
          if (ELASTIC_SIMULATION) then
            call transfer_b_displ_from_device(NDIM*NGLOB_AB,b_displ,Mesh_pointer)
            if (APPROXIMATE_HESS_KL) then
              call transfer_b_veloc_from_device(NDIM*NGLOB_AB,b_veloc,Mesh_pointer)
              call transfer_b_accel_from_device(NDIM*NGLOB_AB,b_accel,Mesh_pointer)
            endif
          endif
        endif

        ! stores wavefield in buffers
        if (ACOUSTIC_SIMULATION) then
          b_potential_acoustic_buffer(:,it_of_this_subset) = b_potential_acoustic(:)
          b_potential_dot_dot_acoustic_buffer(:,it_of_this_subset) = b_potential_dot_dot_acoustic(:)
        endif
        if (ELASTIC_SIMULATION) then
          b_displ_elastic_buffer(:,:,it_of_this_subset) = b_displ(:,:)
          if (APPROXIMATE_HESS_KL) then
            b_veloc_elastic_buffer(:,:,it_of_this_subset) = b_veloc(:,:)
            b_accel_elastic_buffer(:,:,it_of_this_subset) = b_accel(:,:)
          endif
        endif

        ! for noise kernel
        !if (NOISE_TOMOGRAPHY == 3) then
        !  b_noise_surface_movie_buffer(:,:,:,:,it_of_this_subset) = noise_surface_movie(:,:,:,:)
        !endif
      enddo ! subset loop

      ! resets current it and seismo_current positions
      it = it_temp
      seismo_current = seismo_current_temp

      ! computes strain based on current adjoint wavefield
      !if (COMPUTE_AND_STORE_STRAIN) call compute_strain_att()  ! for future use to reduce file sizes...

      ! adjoint wavefield simulation
      do it_of_this_subset = 1, it_subset_end

        ! reads backward/reconstructed wavefield from buffers
        ! note: uses wavefield at corresponding time (NSTEP - it + 1 ), i.e. we have now time-reversed wavefields

        ! only the displacement needs to be stored in memory buffers in order to compute the sensitivity kernels,
        ! not the memory variables R_ij, because the sensitivity kernel calculations only involve the displacement
        ! and the strain, not the stress, and the strain can be recomputed on the fly by computing the gradient
        ! of the displacement read back from the memory buffers
        if (ACOUSTIC_SIMULATION) then
          b_potential_acoustic(:) = b_potential_acoustic_buffer(:,it_subset_end-it_of_this_subset+1)
          b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic_buffer(:,it_subset_end-it_of_this_subset+1)
        endif
        if (ELASTIC_SIMULATION) then
          b_displ(:,:) = b_displ_elastic_buffer(:,:,it_subset_end-it_of_this_subset+1)
          if (APPROXIMATE_HESS_KL) then
            b_veloc(:,:) = b_veloc_elastic_buffer(:,:,it_subset_end-it_of_this_subset+1)
            b_accel(:,:) = b_accel_elastic_buffer(:,:,it_subset_end-it_of_this_subset+1)
          endif
          ! for noise kernel
          !if (NOISE_TOMOGRAPHY == 3) then
          !  noise_surface_movie(:,:,:,:) = b_noise_surface_movie_buffer(:,:,:,:,it_subset_end-it_of_this_subset+1)
          !endif
        endif

        ! transfers wavefields from CPU to GPU
        if (GPU_MODE) then
          ! daniel debug: check if these transfers could be made async to overlap
          if (ACOUSTIC_SIMULATION) then
            call transfer_b_potential_ac_to_device(NGLOB_AB,b_potential_acoustic,Mesh_pointer)
            call transfer_b_potential_dot_dot_ac_to_device(NGLOB_AB,b_potential_dot_dot_acoustic,Mesh_pointer)
          endif
          if (ELASTIC_SIMULATION) then
            call transfer_b_displ_to_device(NDIM*NGLOB_AB,b_displ,Mesh_pointer)
            if (APPROXIMATE_HESS_KL) then
              call transfer_b_veloc_to_device(NDIM*NGLOB_AB,b_veloc,Mesh_pointer)
              call transfer_b_accel_to_device(NDIM*NGLOB_AB,b_accel,Mesh_pointer)
            endif
          endif
        endif

        it = it + 1

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability()
        endif

        ! computes adjoint wavefield
        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then
          if (USE_LDDRK) then
            ! updates wavefields using LDDRK time scheme
            call update_displ_lddrk()
          else
            ! updates wavefields using Newmark time scheme
            call update_displ_Newmark()
          endif

          ! calculates stiffness term
          ! acoustic solver
          ! (needs to be done first, before elastic one)
          if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_forward_calling()
          ! elastic solver
          if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_calling()
          ! poroelastic solver
          if (POROELASTIC_SIMULATION) call compute_forces_poroelastic_calling()
        enddo ! istage

        ! write the seismograms with time shift (GPU_MODE transfer included)
        call write_seismograms()

        ! adjoint simulations: kernels
        call compute_kernels()

      enddo ! subset loop

    end select ! SIMULATION_TYPE

  !
  !---- end of time iteration loop
  !
  enddo   ! end of main time loop

  ! frees undo_attenuation buffers
  if (SIMULATION_TYPE == 3) then
    if (ACOUSTIC_SIMULATION) deallocate(b_potential_acoustic_buffer,b_potential_dot_dot_acoustic_buffer)
    if (ELASTIC_SIMULATION) then
      deallocate(b_displ_elastic_buffer)
      if (APPROXIMATE_HESS_KL) deallocate(b_veloc_elastic_buffer,b_accel_elastic_buffer)
    endif
    ! GPU un-maps memory lock
    if (GPU_MODE) then
      if (ACOUSTIC_SIMULATION) call unregister_host_array(b_potential_acoustic)
      if (ELASTIC_SIMULATION) then
        call unregister_host_array(b_displ)
        if (APPROXIMATE_HESS_KL) call unregister_host_array(b_accel)
      endif
    endif
  endif

  ! user output
  call it_print_elapsed_time()

  ! safety check of last time loop increment
  if (it /= it_end) then
    print *,'Error time increments: it_end = ',it_end,' and last it = ',it,' do not match!'
    call exit_MPI(myrank,'Error invalid time increment ending')
  endif

  !! CD CD added this
  if (RECIPROCITY_AND_KH_INTEGRAL) then
    close(158)
    close(237)
    close(238)
  endif

  ! Transfer fields from GPU card to host for further analysis
  if (GPU_MODE) call it_transfer_from_GPU()

  ! closes energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)

#ifdef VTK_VIS
  if (VTK_MODE) then
    ! frees memory
    call vtk_window_cleanup(do_restart)

    ! check if restart time iterations
    if (do_restart) then
      goto 123
    endif
  endif
#endif

  ! cleanup GPU arrays
  if (GPU_MODE) call it_cleanup_GPU()

  end subroutine iterate_time_undoatt

