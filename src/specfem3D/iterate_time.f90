!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5) then
      call it_check_stability()
    endif

    ! update displacement using Newark time scheme
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
    ! note: this must be read in after the Newark time scheme
    if( SIMULATION_TYPE == 3 .and. it == 1 ) then
     call it_read_foward_arrays()
    endif

    ! write the seismograms with time shift
    if (nrec_local > 0) then
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
    if ( NOISE_TOMOGRAPHY == 1 ) then
      call noise_save_surface_movie(displ, &
                              ibool, &
                              noise_surface_movie,it, &
                              NSPEC_AB,NGLOB_AB, &
                              num_free_surface_faces,free_surface_ispec,free_surface_ijk)
    endif

!
!---- end of time iteration loop
!
  enddo   ! end of main time loop

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!chris: Rewrite to get norm for each material when coupled simulations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! compute maximum of norm of displacement in each slice
  if( ELASTIC_SIMULATION ) &
    Usolidnorm = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2 + displ(3,:)**2))
  if( ACOUSTIC_SIMULATION ) &
      Usolidnormp = maxval(abs(potential_dot_dot_acoustic(:)))
  if( POROELASTIC_SIMULATION ) then
    Usolidnorms = maxval(sqrt(displs_poroelastic(1,:)**2 + displs_poroelastic(2,:)**2 + &
                             displs_poroelastic(3,:)**2))
    Usolidnormw = maxval(sqrt(displw_poroelastic(1,:)**2 + displw_poroelastic(2,:)**2 + &
                             displw_poroelastic(3,:)**2))
  endif

! compute the maximum of the maxima for all the slices using an MPI reduction
  call max_all_cr(Usolidnorm,Usolidnorm_all)
  call max_all_cr(Usolidnormp,Usolidnormp_all)
  call max_all_cr(Usolidnorms,Usolidnorms_all)
  call max_all_cr(Usolidnormw,Usolidnormw_all)

! adjoint simulations
  if( SIMULATION_TYPE == 3 ) then
    if( ELASTIC_SIMULATION ) then
      b_Usolidnorm = maxval(sqrt(b_displ(1,:)**2 + b_displ(2,:)**2 + b_displ(3,:)**2))
    else
      if( ACOUSTIC_SIMULATION ) then
        b_Usolidnorm = maxval(abs(b_potential_dot_dot_acoustic(:)))
      endif
    endif
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
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    if( ELASTIC_SIMULATION ) &
      write(IMAIN,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all
    if( ACOUSTIC_SIMULATION ) &
        write(IMAIN,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all
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
    write(IMAIN,*) 'Estimated remaining time in seconds = ',t_remain
    write(IMAIN,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain

! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
    write(IMAIN,*) 'Estimated total run time in seconds = ',t_total
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

! write time stamp file to give information about progression of simulation
    write(outputname,"('/timestamp',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown')
    write(IOUT,*) 'Time step # ',it
    write(IOUT,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'
    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    if( ELASTIC_SIMULATION ) &
      write(IOUT,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all
    if( ACOUSTIC_SIMULATION ) &
        write(IOUT,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all
    if( POROELASTIC_SIMULATION ) then
        write(IOUT,*) 'Max norm displacement vector Us in all slices (m) = ',Usolidnorms_all
        write(IOUT,*) 'Max norm displacement vector W in all slices (m) = ',Usolidnormw_all
    endif
    ! adjoint simulations
    if (SIMULATION_TYPE == 3) write(IOUT,*) &
           'Max norm U (backward) in all slices = ',b_Usolidnorm_all
    close(IOUT)


! check stability of the code, exit if unstable
! negative values can occur with some compilers when the unstable value is greater
! than the greatest possible floating-point number of the machine
    if(Usolidnorm_all > STABILITY_THRESHOLD .or. Usolidnorm_all < 0 &
     .or. Usolidnormp_all > STABILITY_THRESHOLD .or. Usolidnormp_all < 0 &
     .or. Usolidnorms_all > STABILITY_THRESHOLD .or. Usolidnorms_all < 0 &
     .or. Usolidnormw_all > STABILITY_THRESHOLD .or. Usolidnormw_all < 0) &
        call exit_MPI(myrank,'forward simulation became unstable and blew up')
    ! adjoint simulations
    if(SIMULATION_TYPE == 3 .and. (b_Usolidnorm_all > STABILITY_THRESHOLD &
      .or. b_Usolidnorm_all < 0)) &
        call exit_MPI(myrank,'backward simulation became unstable and blew up')

  endif ! myrank

  end subroutine it_check_stability


!=====================================================================

  subroutine it_update_displacement_scheme()

! explicit Newark time scheme with acoustic & elastic domains:
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
    potential_acoustic(:) = potential_acoustic(:) &
                            + deltat * potential_dot_acoustic(:) &
                            + deltatsqover2 * potential_dot_dot_acoustic(:)
    potential_dot_acoustic(:) = potential_dot_acoustic(:) &
                                + deltatover2 * potential_dot_dot_acoustic(:)
    potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL

    ! time marching potentials
    if(PML) call PML_acoustic_time_march(NSPEC_AB,NGLOB_AB,ibool,&
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
  endif

! updates elastic displacement and velocity
  if( ELASTIC_SIMULATION ) then
    displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
    veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
    if( SIMULATION_TYPE /= 1 ) accel_adj_coupling(:,:) = accel(:,:)
    accel(:,:) = 0._CUSTOM_REAL
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
  if (SIMULATION_TYPE == 3) then
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
  endif

! adjoint simulations: moho kernel
  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    ispec2D_moho_top = 0
    ispec2D_moho_bot = 0
  endif


  end subroutine it_update_displacement_scheme

!=====================================================================

  subroutine it_read_foward_arrays()

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
  endif

  ! elastic wavefields
  if( ELASTIC_SIMULATION ) then
    read(27) b_displ
    read(27) b_veloc
    read(27) b_accel

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
    endif

  endif

  close(27)

  end subroutine it_read_foward_arrays

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
!       we read in the reconstructed wavefield at the end of the time iteration loop, i.e. after the Newark scheme,
!       thus, indexing is NSTEP-it (rather than something like NSTEP-(it-1) )
    if (SIMULATION_TYPE == 3 .and. mod(NSTEP-it,NSTEP_Q_SAVE) == 0) then
      ! reads files content
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") NSTEP-it
      open(unit=27,file=trim(prname_Q)//trim(outputname),status='old',&
            action='read',form='unformatted')
      if( ELASTIC_SIMULATION ) then
        read(27) b_displ
        read(27) b_veloc
        read(27) b_accel
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
      endif
      if( ACOUSTIC_SIMULATION ) then
        read(27) b_potential_acoustic
        read(27) b_potential_dot_acoustic
        read(27) b_potential_dot_dot_acoustic
      endif
      close(27)
    else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. mod(it,NSTEP_Q_SAVE) == 0) then
      ! stores files content
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") it
      open(unit=27,file=trim(prname_Q)//trim(outputname),status='unknown',&
           action='write',form='unformatted')
      if( ELASTIC_SIMULATION ) then
        write(27) displ
        write(27) veloc
        write(27) accel
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
        write(27) potential_acoustic
        write(27) potential_dot_acoustic
        write(27) potential_dot_dot_acoustic
      endif
      close(27)
    endif ! SIMULATION_TYPE
  endif ! it

  end subroutine it_store_attenuation_arrays
