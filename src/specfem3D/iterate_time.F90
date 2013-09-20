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
  use gravity_perturbation, only : gravity_timeseries, GRAVITY_SIMULATION

  implicit none

  !----  create a Gnuplot script to display the energy curve in log scale
  if( OUTPUT_ENERGY .and. myrank == 0) then
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time step number"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(a152)') '#plot "energy.dat" us 1:2 t ''Kinetic Energy'' w l lc 1, "energy.dat" us 1:3 &
                         &t ''Potential Energy'' w l lc 2, "energy.dat" us 1:4 t ''Total Energy'' w l lc 4'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:2 t ''Kinetic Energy'' w l lc 1'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:3 t ''Potential Energy'' w l lc 2'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:4 t ''Total Energy'' w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

  ! open the file in which we will store the energy curve
  if( OUTPUT_ENERGY .and. myrank == 0 ) &
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')

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
    call flush_IMAIN()
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
    if( mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP ) &
      call check_stability()

    ! simulation status output and stability check
    if( OUTPUT_ENERGY ) then
      if( mod(it,NTSTEP_BETWEEN_OUTPUT_ENERGY) == 0 .or. it == 5 .or. it == NSTEP ) &
        call compute_total_energy()
    endif

    ! updates wavefields using Newmark time scheme
    call update_displacement_scheme()

    ! calculates stiffness term
    if( .not. GPU_MODE )then
      ! wavefields on CPU

      ! note: the order of the computations for acoustic and elastic domains is crucial for coupled simulations
      if( SIMULATION_TYPE == 3 ) then
        ! kernel/adjoint simulations

        ! adjoint wavefields
        if( ELASTIC_SIMULATION .and. ACOUSTIC_SIMULATION )then
          ! coupled acoustic-elastic simulations
          ! 1. elastic domain w/ adjoint wavefields
          call compute_forces_viscoelastic()
          ! 2. acoustic domain w/ adjoint wavefields
          call compute_forces_acoustic()
        else
          ! non-coupled simulations
          ! (purely acoustic or elastic)
          if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic()
          if( ELASTIC_SIMULATION ) call compute_forces_viscoelastic()
        endif

        ! backward/reconstructed wavefields
        ! acoustic solver
        ! (needs to be done after elastic one)
        if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic_bpwf()
        ! elastic solver
        ! (needs to be done first, before poroelastic one)
        if( ELASTIC_SIMULATION ) call compute_forces_viscoelastic_bpwf()

      else
        ! forward simulations

        ! 1. acoustic domain
        if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic()
        ! 2. elastic domain
        if( ELASTIC_SIMULATION ) call compute_forces_viscoelastic()
      endif

      ! poroelastic solver
      if( POROELASTIC_SIMULATION ) call compute_forces_poroelastic()

    else
      ! wavefields on GPU
      ! acoustic solver
      if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic_GPU()
      ! elastic solver
      ! (needs to be done first, before poroelastic one)
      if( ELASTIC_SIMULATION ) call compute_forces_viscoelastic_GPU()
    endif

    ! restores last time snapshot saved for backward/reconstruction of wavefields
    ! note: this must be read in after the Newmark time scheme
    if( SIMULATION_TYPE == 3 .and. it == 1 ) then
      call it_read_forward_arrays()
    endif

    ! write the seismograms with time shift (GPU_MODE transfer included)
    if( nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0 ) ) then
      call write_seismograms()
    endif

    ! calculating gravity field at current timestep
    if( GRAVITY_SIMULATION ) call gravity_timeseries()

    ! resetting d/v/a/R/eps for the backward reconstruction with attenuation
    if( ATTENUATION ) then
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
    if( NOISE_TOMOGRAPHY == 1 ) then
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
  if( GPU_MODE ) call it_transfer_from_GPU()

!----  close energy file
  if( OUTPUT_ENERGY .and. myrank == 0 ) close(IOUT_ENERGY)

  end subroutine iterate_time


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
  if (ADIOS_FOR_FORWARD_ARRAYS) then
    call read_forward_arrays_adios()
  else
    ! reads in wavefields
    open(unit=IIN,file=trim(prname)//'save_forward_arrays.bin',status='old',&
          action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error: opening save_forward_arrays'
      print*,'path: ',trim(prname)//'save_forward_arrays.bin'
      call exit_mpi(myrank,'error open file save_forward_arrays.bin')
    endif

    if( ACOUSTIC_SIMULATION ) then
      read(IIN) b_potential_acoustic
      read(IIN) b_potential_dot_acoustic
      read(IIN) b_potential_dot_dot_acoustic
    endif

    ! elastic wavefields
    if( ELASTIC_SIMULATION ) then
      read(IIN) b_displ
      read(IIN) b_veloc
      read(IIN) b_accel
      ! memory variables if attenuation
      if( ATTENUATION ) then
        if(FULL_ATTENUATION_SOLID) read(IIN) b_R_trace
        read(IIN) b_R_xx
        read(IIN) b_R_yy
        read(IIN) b_R_xy
        read(IIN) b_R_xz
        read(IIN) b_R_yz
        if(FULL_ATTENUATION_SOLID) read(IIN) b_epsilondev_trace
        read(IIN) b_epsilondev_xx
        read(IIN) b_epsilondev_yy
        read(IIN) b_epsilondev_xy
        read(IIN) b_epsilondev_xz
        read(IIN) b_epsilondev_yz
      endif ! ATTENUATION
    endif

    ! poroelastic wavefields
    if( POROELASTIC_SIMULATION ) then
      read(IIN) b_displs_poroelastic
      read(IIN) b_velocs_poroelastic
      read(IIN) b_accels_poroelastic
      read(IIN) b_displw_poroelastic
      read(IIN) b_velocw_poroelastic
      read(IIN) b_accelw_poroelastic
    endif

    close(IIN)
  endif

  if(GPU_MODE) then
    if( ACOUSTIC_SIMULATION ) then
    ! transfers fields onto GPU
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,      &
                                          b_potential_dot_dot_acoustic,  &
                                          Mesh_pointer)
    endif
    ! elastic wavefields
    if( ELASTIC_SIMULATION ) then
      ! puts elastic wavefield to GPU
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
      ! memory variables if attenuation
      if( ATTENUATION ) then
        call transfer_b_fields_att_to_device(Mesh_pointer,                    &
                           b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,                &
                           size(b_R_xx),                                      &
                           b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,   &
                           b_epsilondev_xz,b_epsilondev_yz,                   &
                           size(b_epsilondev_xx))
      endif
    endif
  endif

  end subroutine it_read_forward_arrays

!=====================================================================

  subroutine it_store_attenuation_arrays()

! resetting d/v/a/R/eps for the backward reconstruction with attenuation

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  integer :: ier

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
      open(unit=IIN,file=trim(prname_Q)//trim(outputname),status='old',&
            action='read',form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: opening save_Q_arrays'
        print*,'path: ',trim(prname_Q)//trim(outputname)
        call exit_mpi(myrank,'error open file save_Q_arrays_***.bin for reading')
      endif

      if( ELASTIC_SIMULATION ) then
        ! reads arrays from disk files
        read(IIN) b_displ
        read(IIN) b_veloc
        read(IIN) b_accel

        ! puts elastic fields onto GPU
        if(GPU_MODE) then
          ! wavefields
          call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel, Mesh_pointer)
        endif

        if(FULL_ATTENUATION_SOLID) read(IIN) b_R_trace
        read(IIN) b_R_xx
        read(IIN) b_R_yy
        read(IIN) b_R_xy
        read(IIN) b_R_xz
        read(IIN) b_R_yz
        if(FULL_ATTENUATION_SOLID) read(IIN) b_epsilondev_trace
        read(IIN) b_epsilondev_xx
        read(IIN) b_epsilondev_yy
        read(IIN) b_epsilondev_xy
        read(IIN) b_epsilondev_xz
        read(IIN) b_epsilondev_yz

        ! puts elastic fields onto GPU
        if(GPU_MODE) then
          ! attenuation arrays
          call transfer_b_fields_att_to_device(Mesh_pointer, &
                  b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,size(b_R_xx), &
!!!               b_R_trace,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,size(b_R_xx), &
                  b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz, &
!!!               b_epsilondev_trace,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz, &
                  size(b_epsilondev_xx))
        endif
      endif

      if( ACOUSTIC_SIMULATION ) then
        ! reads arrays from disk files
        read(IIN) b_potential_acoustic
        read(IIN) b_potential_dot_acoustic
        read(IIN) b_potential_dot_dot_acoustic

        ! puts acoustic fields onto GPU
        if(GPU_MODE) &
          call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                              b_potential_dot_acoustic, b_potential_dot_dot_acoustic, Mesh_pointer)

      endif
      close(IIN)

    else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. mod(it,NSTEP_Q_SAVE) == 0) then
      ! stores files content
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") it
      open(unit=IOUT,file=trim(prname_Q)//trim(outputname),status='unknown',&
           action='write',form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: opening save_Q_arrays'
        print*,'path: ',trim(prname_Q)//trim(outputname)
        call exit_mpi(myrank,'error open file save_Q_arrays_***.bin for writing')
      endif

      if( ELASTIC_SIMULATION ) then
        ! gets elastic fields from GPU onto CPU
        if(GPU_MODE) then
          call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc, accel, Mesh_pointer)
        endif

        ! writes to disk file
        write(IOUT) displ
        write(IOUT) veloc
        write(IOUT) accel

        if(GPU_MODE) then
          ! attenuation arrays
          call transfer_fields_att_from_device(Mesh_pointer, &
                     R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
!!!                  R_trace,R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                     epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
!!!                  epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                     size(epsilondev_xx))
        endif

        if(FULL_ATTENUATION_SOLID) write(IOUT) R_trace
        write(IOUT) R_xx
        write(IOUT) R_yy
        write(IOUT) R_xy
        write(IOUT) R_xz
        write(IOUT) R_yz
        if(FULL_ATTENUATION_SOLID) write(IOUT) epsilondev_trace
        write(IOUT) epsilondev_xx
        write(IOUT) epsilondev_yy
        write(IOUT) epsilondev_xy
        write(IOUT) epsilondev_xz
        write(IOUT) epsilondev_yz
      endif

      if( ACOUSTIC_SIMULATION ) then
        ! gets acoustic fields from GPU onto CPU
        if(GPU_MODE) then
          call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                                              potential_dot_acoustic, potential_dot_dot_acoustic, Mesh_pointer)
        endif

        ! writes to disk file
        write(IOUT) potential_acoustic
        write(IOUT) potential_dot_acoustic
        write(IOUT) potential_dot_dot_acoustic
      endif
      close(IOUT)

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

  if( myrank == 0 ) then
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
  if( SIMULATION_TYPE == 1 .and. SAVE_FORWARD ) then

    ! acoustic potentials
    if( ACOUSTIC_SIMULATION ) &
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                                          potential_dot_acoustic, potential_dot_dot_acoustic, &
                                          Mesh_pointer)

    ! elastic wavefield
    if( ELASTIC_SIMULATION ) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)

      if (ATTENUATION) &
        call transfer_fields_att_from_device(Mesh_pointer, &
                    R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
!!!                 R_trace,R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
!!!                 epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                    size(epsilondev_xx))

    endif
  else if( SIMULATION_TYPE == 3 ) then

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
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,cijkl_kl,NSPEC_AB)
    endif

    ! specific noise strength kernel
    if( NOISE_TOMOGRAPHY == 3 ) then
      call transfer_kernels_noise_to_host(Mesh_pointer,Sigma_kl,NSPEC_AB)
    endif

    ! approximative hessian for preconditioning kernels
    if( APPROXIMATE_HESS_KL ) then
      if( ELASTIC_SIMULATION ) call transfer_kernels_hess_el_tohost(Mesh_pointer,hess_kl,NSPEC_AB)
      if( ACOUSTIC_SIMULATION ) call transfer_kernels_hess_ac_tohost(Mesh_pointer,hess_ac_kl,NSPEC_AB)
    endif

  endif

  ! frees allocated memory on GPU
  call prepare_cleanup_device(Mesh_pointer, &
                              ACOUSTIC_SIMULATION,ELASTIC_SIMULATION, &
                              STACEY_ABSORBING_CONDITIONS,NOISE_TOMOGRAPHY,COMPUTE_AND_STORE_STRAIN, &
                              ATTENUATION,ANISOTROPY,APPROXIMATE_OCEAN_LOAD, &
                              APPROXIMATE_HESS_KL)

  end subroutine it_transfer_from_GPU
