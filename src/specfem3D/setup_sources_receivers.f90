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

  subroutine setup_sources_receivers()

  use specfem_par

  use kdtree_search, only: kdtree_delete,kdtree_nodes_location,kdtree_nodes_index

  implicit none

  ! builds search tree
  if (.not. DO_BRUTE_FORCE_POINT_SEARCH) call setup_search_kdtree()

  ! sets number of timestep parameters
  call setup_timesteps()

  ! locates sources and determines simulation start time t0
  call setup_sources()

  ! reads in stations file and locates receivers
  call setup_receivers()

  ! pre-compute source arrays
  call setup_sources_precompute_arrays()

  ! pre-compute receiver interpolation factors
  call setup_receivers_precompute_intp()

  ! write source and receiver VTK files for Paraview
  call setup_sources_receivers_VTKfile()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all the slices'
    write(IMAIN,*)
    if (NSOURCES > 1) then
      write(IMAIN,*) 'Using ',NSOURCES,' point sources'
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! frees tree memory
  if (.not. DO_BRUTE_FORCE_POINT_SEARCH) then
    ! deletes tree arrays
    deallocate(kdtree_nodes_location)
    deallocate(kdtree_nodes_index)
    ! deletes search tree nodes
    call kdtree_delete()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine setup_sources_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_timesteps()

  use specfem_par, only: myrank,NSTEP,NSOURCES,SIMULATION_TYPE, &
    NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC, &
    USE_EXTERNAL_SOURCE_FILE,NSTEP_STF,NSOURCES_STF, &
    SAVE_ALL_SEISMOS_IN_ONE_FILE

  implicit none

  ! subsets used to save seismograms must not be larger than the whole time series
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP

  ! subsets used to save adjoint sources must not be larger than the whole time series,
  ! otherwise we waste memory
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    if (NTSTEP_BETWEEN_READ_ADJSRC > NSTEP) NTSTEP_BETWEEN_READ_ADJSRC = NSTEP
  endif

  !! VM VM set the size of user_source_time_function
  if (USE_EXTERNAL_SOURCE_FILE) then
     NSTEP_STF = NSTEP
     NSOURCES_STF = NSOURCES
  else !! We don't need the array user_source_time_function : use a small dummy array
     NSTEP_STF = 1
     NSOURCES_STF = 1
  endif

  ! check
  if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
    ! requires to have full length of seismograms
    if (NTSTEP_BETWEEN_OUTPUT_SEISMOS < NSTEP) then
      if (myrank == 0) then
        print *, 'Error: Setting SAVE_ALL_SEISMOS_IN_ONE_FILE to .true. requires NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP'
      endif
      call exit_MPI(myrank, &
                    'Setting SAVE_ALL_SEISMOS_IN_ONE_FILE not supported for current NTSTEP_BETWEEN_OUTPUT_SEISMOS in Par_file')
    endif
  endif

  end subroutine setup_timesteps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  implicit none

  double precision :: min_tshift_src_original
  integer :: isource,ier
  character(len=MAX_STRING_LEN) :: SOURCE_FILE,path_to_add

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources:',NSOURCES
    call flush_IMAIN()
  endif

  ! allocate arrays for source
  allocate(islice_selected_source(NSOURCES), &
           ispec_selected_source(NSOURCES), &
           Mxx(NSOURCES), &
           Myy(NSOURCES), &
           Mzz(NSOURCES), &
           Mxy(NSOURCES), &
           Mxz(NSOURCES), &
           Myz(NSOURCES), &
           xi_source(NSOURCES), &
           eta_source(NSOURCES), &
           gamma_source(NSOURCES), &
           tshift_src(NSOURCES), &
           hdur(NSOURCES), &
           hdur_Gaussian(NSOURCES), &
           utm_x_source(NSOURCES), &
           utm_y_source(NSOURCES), &
           nu_source(NDIM,NDIM,NSOURCES), &
           stat=ier)
  if (ier /= 0) stop 'error allocating arrays for sources'

  if (USE_FORCE_POINT_SOURCE) then
    allocate(factor_force_source(NSOURCES), &
             comp_dir_vect_source_E(NSOURCES), &
             comp_dir_vect_source_N(NSOURCES), &
             comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
  else
    allocate(factor_force_source(1), &
             comp_dir_vect_source_E(1), &
             comp_dir_vect_source_N(1), &
             comp_dir_vect_source_Z_UP(1),stat=ier)
  endif
  if (ier /= 0) stop 'error allocating arrays for force point sources'

  !! allocate the array contains the user defined source time function
  allocate(user_source_time_function(NSTEP_STF, NSOURCES_STF),stat=ier)
  if (ier /= 0) stop 'error allocating arrays for user sources time function'

  if (USE_FORCE_POINT_SOURCE) then
    SOURCE_FILE = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FORCESOLUTION'
  else
    SOURCE_FILE = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'CMTSOLUTION'
  endif

  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run
  ! (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    SOURCE_FILE = path_to_add(1:len_trim(path_to_add))//SOURCE_FILE(1:len_trim(SOURCE_FILE))
  endif

  ! locate sources in the mesh
  call locate_source(SOURCE_FILE,tshift_src,min_tshift_src_original,utm_x_source,utm_y_source, &
                     hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                     islice_selected_source,ispec_selected_source, &
                     factor_force_source,comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP, &
                     xi_source,eta_source,gamma_source,nu_source,user_source_time_function)

  ! determines onset time
  call setup_stf_constants(hdur,hdur_Gaussian,tshift_src,min_tshift_src_original, &
                           islice_selected_source,ispec_selected_source,t0)

  ! count number of sources located in this slice
  nsources_local = 0
  do isource = 1, NSOURCES
    if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
  enddo

  ! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there
  call setup_sources_check_acoustic()

  ! prints source time functions to output files
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  call get_run_number_of_the_source()

  end subroutine setup_sources
!
!-------------------------------------------------------------------------------------------------
!
!EB EB When NB_RUNS_ACOUSTIC_GPU > 1, the file SOURCE_FILE actually contains the sources for all the runs.
!This routine is intended to get the array that contains the run number of each source described in SOURCE_FILE.
!The line i of the file run_number_of_the_source contains the run number \in [ 0;NB_RUNS_ACOUSTIC_GPU-1] of the source i
  subroutine get_run_number_of_the_source()

  use constants
  use specfem_par, only: run_number_of_the_source,NSOURCES
  character(len=MAX_STRING_LEN) :: filename,string
  integer :: ier,isource,icounter

  allocate(run_number_of_the_source(NSOURCES))

  if (NB_RUNS_ACOUSTIC_GPU == 1) then
    run_number_of_the_source(:) = 0
  else

    filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'run_number_of_the_source'
    open(unit=IIN,file=filename,status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',filename
      stop 'Error opening run_number_of_the_source file'
    endif

    ! Checks if the number of lines is correct
    icounter = 0
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) string
      if (ier == 0) icounter = icounter + 1
    enddo
    close(IIN)
    if (icounter /= NSOURCES) stop 'Error total number of lines in run_number_of_the_source file is not equal to NSOURCES'

    ! Fills the array run_number_of_the_source
    open(unit=IIN,file=filename,status='old',action='read')
    ! reads run number for each source
    do isource = 1,NSOURCES
      read(IIN,"(a)") string
      read(string,*) run_number_of_the_source(isource)
    enddo
  endif

  end subroutine get_run_number_of_the_source

!
!-------------------------------------------------------------------------------------------------
!
  subroutine setup_stf_constants(hdur,hdur_Gaussian,tshift_src,min_tshift_src_original, &
                                  islice_selected_source,ispec_selected_source,t0)

  use constants
  use specfem_par, only: myrank,NSOURCES,MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,HDUR_MOVIE, &
                          USE_RICKER_TIME_FUNCTION,USE_EXTERNAL_SOURCE_FILE
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_movie

  implicit none

  double precision, dimension(NSOURCES) :: tshift_src,hdur,hdur_Gaussian
  double precision :: min_tshift_src_original,t0
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source

  ! local parameters
  double precision :: t0_acoustic
  integer :: isource,ispec

  if (abs(minval(tshift_src)) > TINYVAL) &
    call exit_MPI(myrank,'one tshift_src must be zero, others must be positive')

  ! filter source time function by Gaussian with hdur = HDUR_MOVIE when outputing movies or shakemaps
  if (MOVIE_SURFACE .or. MOVIE_VOLUME .or. CREATE_SHAKEMAP) then
    hdur = sqrt(hdur**2 + HDUR_MOVIE**2)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Each source is being convolved with HDUR_MOVIE = ',HDUR_MOVIE
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! convert the half duration for triangle STF to the one for Gaussian STF
  hdur_Gaussian(:) = hdur(:)/SOURCE_DECAY_MIMIC_TRIANGLE

  ! define t0 as the earliest start time
  ! note: an earlier start time also reduces numerical noise due to a
  !          non-zero offset at the beginning of the source time function
  t0 = - 2.0d0 * minval(tshift_src(:) - hdur(:))   ! - 1.5d0 * minval(tshift_src-hdur)

  ! uses an earlier start time if source is acoustic with a Gaussian source time function
  t0_acoustic = 0.0d0
  do isource = 1,NSOURCES
    if (myrank == islice_selected_source(isource)) then
      ispec = ispec_selected_source(isource)
      if (ispec_is_acoustic(ispec)) then
        ! uses an earlier start time
        t0_acoustic = - 3.0d0 * ( tshift_src(isource) - hdur(isource) )
        if (t0_acoustic > t0 ) t0 = t0_acoustic
      endif
    endif
  enddo
  ! passes maximum value to all processes
  ! note: t0 is defined positive and will be subtracted from simulation time (it-1)*DT
  t0_acoustic = t0
  call max_all_all_dp(t0_acoustic,t0)

  ! point force sources will start depending on the frequency given by hdur
  !if (USE_FORCE_POINT_SOURCE .or. USE_RICKER_TIME_FUNCTION) then
  !! COMMENTED BY FS FS -> account for the case USE_FORCE_POINT_SOURCE but NOT using a ricker (i.e. using a Gaussian),
  ! in this case the above defined t0 = - 2.0d0 * minval(tshift_src(:) - hdur(:)) is correct
  ! (analogous to using error function in case of moment tensor sources). You only need to be aware that hdur=0
  ! then has a different behaviour for point forces (compared with moment tensor sources):
  ! then hdur is set to TINYVAL and NOT to 5*DT as in case of moment tensor source (Heaviside)
  if (USE_RICKER_TIME_FUNCTION) then !! ADDED BY FS FS
    ! note: point force sources will give the dominant frequency in hdur,
    !       thus the main period is 1/hdur.
    !       also, these sources use a Ricker source time function instead of a Gaussian.
    !       for a Ricker source time function, a start time ~1.2 * dominant_period is a good choice
    t0 = - 1.2d0 * minval(tshift_src(:) - 1.0d0/hdur(:))
  endif

  !! VM VM for external source the time will begin with simulation
  if (USE_EXTERNAL_SOURCE_FILE) then
    t0 = 0.d0
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'External STF:'
      write(IMAIN,*) '  simulation start time set to zero: ', t0
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! checks if user set USER_T0 to fix simulation start time
  ! note: USER_T0 has to be positive
  if (USER_T0 > 0.d0) then
    ! user cares about origin time and time shifts of the CMTSOLUTION
    ! and wants to fix simulation start time to a constant start time
    ! time 0 on time axis will correspond to given origin time

    ! notifies user
    if (myrank == 0) then
      write(IMAIN,*) 'USER_T0: ',USER_T0
      write(IMAIN,*) 't0: ',t0,'min_tshift_src_original: ',min_tshift_src_original
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! checks if automatically set t0 is too small
    ! note: min_tshift_src_original can be a positive or negative time shift (minimum from all tshift)
    if (t0 <= USER_T0 + min_tshift_src_original) then
      ! by default, tshift_src(:) holds relative time shifts with a minimum time shift set to zero
      ! re-adds (minimum) original time shift such that sources will kick in
      ! according to their absolute time shift
      tshift_src(:) = tshift_src(:) + min_tshift_src_original

      ! sets new simulation start time such that
      ! simulation starts at t = - t0 = - USER_T0
      t0 = USER_T0

      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) '  set new simulation start time: ', - t0
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    else
      ! start time needs to be at least t0 for numerical stability
      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) 'error: USER_T0 is too small'
        write(IMAIN,*) '       must make one of three adjustements:'
        write(IMAIN,*) '       - increase USER_T0 to be at least: ',t0-min_tshift_src_original
        write(IMAIN,*) '       - decrease time shift in CMTSOLUTION file'
        write(IMAIN,*) '       - decrease hdur in CMTSOLUTION file'
        call flush_IMAIN()
      endif
      call exit_mpi(myrank,'error USER_T0 is set but too small')
    endif
  else if (USER_T0 < 0.d0) then
    if (myrank == 0) then
      write(IMAIN,*) 'error: USER_T0 is negative, must be set zero or positive!'
      call flush_IMAIN()
    endif
    call exit_mpi(myrank,'error negative USER_T0 parameter in constants.h')
  endif

  end subroutine setup_stf_constants

!
!-------------------------------------------------------------------------------------------------
!


  subroutine setup_sources_check_acoustic()

! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there

  use specfem_par
  use specfem_par_acoustic
  implicit none

  integer :: isource,ixmin,ixmax,iymin,iymax,izmin,izmax,iface,ispec
  logical :: is_on,is_on_all

! outputs a warning in case of an acoustic source lying on the free surface
  do isource = 1,NSOURCES
    ! checks if source is close to face
    is_on = .false.

    ! only receivers in this process
    if (myrank == islice_selected_source(isource)) then

      ispec = ispec_selected_source(isource)
      ! only if receiver is in an acoustic element
      if (ispec_is_acoustic(ispec)) then

        ! checks with free surface face
        do iface = 1,num_free_surface_faces

          if (ispec == free_surface_ispec(iface)) then

            ! determine face
            ixmin = minval( free_surface_ijk(1,:,iface) )
            ixmax = maxval( free_surface_ijk(1,:,iface) )

            iymin = minval( free_surface_ijk(2,:,iface) )
            iymax = maxval( free_surface_ijk(2,:,iface) )

            izmin = minval( free_surface_ijk(3,:,iface) )
            izmax = maxval( free_surface_ijk(3,:,iface) )

            ! xmin face
            if (ixmin == 1 .and. ixmax == 1) then
              if (xi_source(isource) < -0.99d0) is_on = .true.
            ! xmax face
            else if (ixmin == NGLLX .and. ixmax == NGLLX) then
              if (xi_source(isource) > 0.99d0) is_on = .true.
            ! ymin face
            else if (iymin == 1 .and. iymax == 1) then
              if (eta_source(isource) < -0.99d0) is_on = .true.
            ! ymax face
            else if (iymin == NGLLY .and. iymax == NGLLY) then
              if (eta_source(isource) > 0.99d0) is_on = .true.
            ! zmin face
            else if (izmin == 1 .and. izmax == 1) then
              if (gamma_source(isource) < -0.99d0) is_on = .true.
            ! zmax face
            else if (izmin == NGLLZ .and. izmax == NGLLZ) then
              if (gamma_source(isource) > 0.99d0) is_on = .true.
            endif

          endif ! free_surface_ispec
        enddo ! iface
      endif ! ispec_is_acoustic
    endif ! islice_selected_rec

    ! user output
    call any_all_l( is_on, is_on_all )
    if (myrank == 0 .and. is_on_all) then
      write(IMAIN,*) '********************************************************************'
      write(IMAIN,*) '*** source: ',isource,'                                          ***'
      write(IMAIN,*) '*** Warning: acoustic source located exactly on the free surface ***'
      write(IMAIN,*) '*** will be zeroed                                               ***'
      write(IMAIN,*) '********************************************************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  enddo ! num_free_surface_faces

  end subroutine setup_sources_check_acoustic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers()

  use specfem_par
  use specfem_par_acoustic

  implicit none

  ! local parameters
  integer :: nrec_simulation
  integer :: irec,isource,ier
  character(len=MAX_STRING_LEN) :: path_to_add

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! reads in station file
  if (SIMULATION_TYPE == 1) then
    rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS'
    filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_FILTERED'
  else
    rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT'
    filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT_FILTERED'
  endif

  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run
  ! (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    rec_filename = path_to_add(1:len_trim(path_to_add))//rec_filename(1:len_trim(rec_filename))
    filtered_rec_filename = path_to_add(1:len_trim(path_to_add))//filtered_rec_filename(1:len_trim(filtered_rec_filename))
  endif

  call station_filter(rec_filename,filtered_rec_filename,nrec)

  if (nrec < 1) call exit_MPI(myrank,'need at least one receiver')
  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*)
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) 'Total number of receivers = ', nrec
    else
      write(IMAIN,*) 'Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate memory for receiver arrays, i.e. stations given in STATIONS file
  allocate(islice_selected_rec(nrec), &
           ispec_selected_rec(nrec), &
           xi_receiver(nrec), &
           eta_receiver(nrec), &
           gamma_receiver(nrec), &
           station_name(nrec), &
           network_name(nrec), &
           nu(NDIM,NDIM,nrec), &
           stat=ier)
  if (ier /= 0) stop 'error allocating arrays for receivers'

  ! locate receivers in the mesh
  call locate_receivers(filtered_rec_filename,nrec,islice_selected_rec,ispec_selected_rec, &
                        xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                        utm_x_source(1),utm_y_source(1))

  ! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! number of receivers are given by stations
    ! in STATIONS (forward runs) or STATIONS_ADJOINT (kernel runs) file
    nrec_simulation = nrec
    do irec = 1,nrec
      if (myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    ! adjoint simulation:
    ! station locations (in STATIONS_ADJOINT file) become adjoint sources
    ! and source locations (in CMTSOLUTION file) become adjoint "receivers"
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if (myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

! check that the sum of the number of receivers in each slice is nrec
  call sum_all_i(nrec_local,nrec_tot_found)
  if (myrank == 0) then
    if (nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    endif
  endif
  call synchronize_all()

! checks if acoustic receiver is exactly on the free surface because pressure is zero there
  call setup_receivers_check_acoustic()

  end subroutine setup_receivers


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_check_acoustic()

! checks if acoustic receiver is exactly on the free surface because pressure is zero there

  use specfem_par
  use specfem_par_acoustic
  implicit none

  integer :: irec,ixmin,ixmax,iymin,iymax,izmin,izmax,iface,ispec
  logical :: is_on,is_on_all

! outputs a warning in case the receiver is lying on the free surface
  do irec = 1,nrec

    ! checks if receiver is close to face
    is_on = .false.

    ! only receivers in this process
    if (myrank == islice_selected_rec(irec)) then

      ispec = ispec_selected_rec(irec)
      ! only if receiver is in an acoustic element
      if (ispec_is_acoustic(ispec)) then

        ! checks with free surface face
        do iface = 1,num_free_surface_faces

          if (ispec == free_surface_ispec(iface)) then

            ! determine face
            ixmin = minval( free_surface_ijk(1,:,iface) )
            ixmax = maxval( free_surface_ijk(1,:,iface) )

            iymin = minval( free_surface_ijk(2,:,iface) )
            iymax = maxval( free_surface_ijk(2,:,iface) )

            izmin = minval( free_surface_ijk(3,:,iface) )
            izmax = maxval( free_surface_ijk(3,:,iface) )

            ! xmin face
            if (ixmin == 1 .and. ixmax == 1) then
              if (xi_receiver(irec) < -0.99d0) is_on = .true.
            ! xmax face
            else if (ixmin == NGLLX .and. ixmax == NGLLX) then
              if (xi_receiver(irec) > 0.99d0) is_on = .true.
            ! ymin face
            else if (iymin == 1 .and. iymax == 1) then
              if (eta_receiver(irec) < -0.99d0) is_on = .true.
            ! ymax face
            else if (iymin == NGLLY .and. iymax == NGLLY) then
              if (eta_receiver(irec) > 0.99d0) is_on = .true.
            ! zmin face
            else if (izmin == 1 .and. izmax == 1) then
              if (gamma_receiver(irec) < -0.99d0) is_on = .true.
            ! zmax face
            else if (izmin == NGLLZ .and. izmax == NGLLZ) then
              if (gamma_receiver(irec) > 0.99d0) is_on = .true.
            endif

          endif ! free_surface_ispec
        enddo ! iface
      endif ! ispec_is_acoustic
    endif ! islice_selected_rec

    ! user output
    call any_all_l( is_on, is_on_all )
    if (myrank == 0 .and. is_on_all) then
      ! limits user output if too many receivers
      if (nrec < 1000 .and. (.not. SU_FORMAT )) then
        write(IMAIN,*) '**********************************************************************'
        write(IMAIN,*) '*** station:',irec,'                                               ***'
        write(IMAIN,*) '*** Warning: acoustic receiver located exactly on the free surface ***'
        write(IMAIN,*) '*** Warning: tangential component will be zero there               ***'
        write(IMAIN,*) '**********************************************************************'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    endif

  enddo ! num_free_surface_faces

  end subroutine setup_receivers_check_acoustic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_precompute_arrays()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic

  implicit none

  real(kind=CUSTOM_REAL) :: factor_source
  real(kind=CUSTOM_REAL) :: junk

  integer :: isource,ispec
  integer :: irec
  integer :: icomp,itime,nadj_files_found,nadj_files_found_tot,ier

  character(len=3),dimension(NDIM) :: comp
  character(len=MAX_STRING_LEN) :: filename

  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas
  double precision :: norm,comp_x,comp_y,comp_z
  double precision :: sizeval
  logical :: does_source_encoding

  ! user output info
  ! sources
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      ! note: all process allocate the full sourcearrays array
      ! sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES)
      sizeval = dble(NSOURCES) * dble(NDIM * NGLLX * NGLLY * NGLLZ * CUSTOM_REAL / 1024. / 1024. )
      ! outputs info
      write(IMAIN,*) 'source arrays:'
      write(IMAIN,*) '  number of sources is ',NSOURCES
      write(IMAIN,*) '  size of source array                 = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! for source encoding (acoustic sources only so far)
  if (USE_SOURCE_ENCODING) then
    allocate(pm1_source_encoding(NSOURCES),stat=ier)
  else
    allocate(pm1_source_encoding(1),stat=ier)
  endif
  if (ier /= 0) stop 'error allocating arrays for sources'
  pm1_source_encoding(:) = 1._CUSTOM_REAL
  does_source_encoding = .false.

  ! forward simulations
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    allocate(sourcearray(NDIM,NGLLX,NGLLY,NGLLZ), &
             sourcearrays(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating array sourcearray'

    ! compute source arrays
    do isource = 1,NSOURCES

      ! initializes
      sourcearray(:,:,:,:) = ZERO

      !   check that the source slice number is okay
      if (islice_selected_source(isource) < 0 .or. islice_selected_source(isource) > NPROC-1) &
        call exit_MPI(myrank,'something is wrong with the source slice number')

      !   compute source arrays in source slice
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        ! compute Lagrange polynomials at the source location
        call lagrange_any(xi_source(isource),NGLLX,xigll,hxis,hpxis)
        call lagrange_any(eta_source(isource),NGLLY,yigll,hetas,hpetas)
        call lagrange_any(gamma_source(isource),NGLLZ,zigll,hgammas,hpgammas)

        if (USE_FORCE_POINT_SOURCE) then
          ! use of FORCESOLUTION files
          !
          ! note: for use_force_point_source xi/eta/gamma are also in the range [-1,1], for exact positioning
          factor_source = factor_force_source(isource)

          ! elastic sources
          ! note: point forces in poroelastic element act on both, solid and fluid parts (see compute_add_sources_poroelastic)
          !       we allow thus the force to have a given direction.
          if (ispec_is_elastic(ispec) .or. ispec_is_poroelastic(ispec)) then
            ! length of component vector
            norm = dsqrt( comp_dir_vect_source_E(isource)**2 &
                        + comp_dir_vect_source_N(isource)**2 &
                        + comp_dir_vect_source_Z_UP(isource)**2 )
            ! checks norm of component vector
            if (norm < TINYVAL) then
              call exit_MPI(myrank,'error force point source: component vector has (almost) zero norm')
            endif
            ! normalizes given component vector
            comp_x = comp_dir_vect_source_E(isource)/norm
            comp_y = comp_dir_vect_source_N(isource)/norm
            comp_z = comp_dir_vect_source_Z_UP(isource)/norm
            call compute_arrays_source_forcesolution(sourcearray,hxis,hetas,hgammas,factor_source, &
                                                     comp_x,comp_y,comp_z,nu_source(:,:,isource))
          endif

          ! acoustic sources
          if (ispec_is_acoustic(ispec)) then
            ! dipole
            if (DIPOLE_SOURCE_IN_FLUID) then
              ! length of component vector
              norm = dsqrt( comp_dir_vect_source_E(isource)**2 &
                          + comp_dir_vect_source_N(isource)**2 &
                          + comp_dir_vect_source_Z_UP(isource)**2 )
              comp_x = comp_dir_vect_source_E(isource)/norm
              comp_y = comp_dir_vect_source_N(isource)/norm
              comp_z = comp_dir_vect_source_Z_UP(isource)/norm
              write(*,*) " DIPOLE FLUID ", comp_x, comp_y, comp_z
              call compute_arrays_source_forcesolution_fluid(ispec, sourcearray, hxis,hetas,hgammas,hpxis,hpetas,hpgammas, &
                                                             factor_source,comp_x,comp_y,comp_z,nu_source(:,:,isource))
            else
              ! point force
              ! identical source array components in x,y,z-direction
              ! source array interpolated on all element GLL points
              ! we use an tilted force defined by its magnitude and the projections
              ! of an arbitrary (non-unitary) direction vector on the E/N/Z_UP basis
              comp_x = 1.0d0
              comp_y = 1.0d0
              comp_z = 1.0d0
              call compute_arrays_source_forcesolution(sourcearray,hxis,hetas,hgammas,factor_source, &
                                                       comp_x,comp_y,comp_z,nu_source(:,:,isource))
            endif
          endif ! acoustic

        else
          ! use of CMTSOLUTION files

          ! elastic or poroelastic moment tensor source
          if (ispec_is_elastic(ispec) .or. ispec_is_poroelastic(ispec)) then
            call compute_arrays_source_cmt(ispec,sourcearray, &
                                           hxis,hetas,hgammas,hpxis,hpetas,hpgammas, &
                                           Mxx(isource),Myy(isource),Mzz(isource), &
                                           Mxy(isource),Mxz(isource),Myz(isource))
          endif

          ! acoustic case
          if (ispec_is_acoustic(ispec)) then
            ! scalar moment of moment tensor values read in from CMTSOLUTION
            ! note: M0 by Dahlen and Tromp, eq. 5.91
            factor_source = 1.0/sqrt(2.0) * sqrt( Mxx(isource)**2 + Myy(isource)**2 + Mzz(isource)**2 &
                                                  + 2*( Myz(isource)**2 + Mxz(isource)**2 + Mxy(isource)**2) )

            ! scales source such that it would be equivalent to explosion source moment tensor,
            ! where Mxx=Myy=Mzz, others Mxy,.. = zero, in equivalent elastic media
            ! (and getting rid of 1/sqrt(2) factor from scalar moment tensor definition above)
            factor_source = factor_source * sqrt(2.0) / sqrt(3.0)

            ! source encoding
            ! determines factor +/-1 depending on sign of moment tensor
            ! (see e.g. Krebs et al., 2009. Fast full-wavefield seismic inversion using encoded sources,
            !   Geophysics, 74 (6), WCC177-WCC188.)
            if (USE_SOURCE_ENCODING) then
              pm1_source_encoding(isource) = sign(1.0d0,Mxx(isource))
              does_source_encoding = .true.
            endif

            comp_x = 1.0d0
            comp_y = 1.0d0
            comp_z = 1.0d0

            ! source array interpolated on all element GLL points
            call compute_arrays_source_forcesolution(sourcearray,hxis,hetas,hgammas,factor_source, &
                                                     comp_x,comp_y,comp_z,nu_source(:,:,isource))
          endif

        endif

        ! stores source excitations
        sourcearrays(isource,:,:,:,:) = sourcearray(:,:,:,:)

      endif
    enddo
  else
    ! SIMULATION_TYPE == 2
    ! allocate dummy array (needed for subroutine calls)
    allocate(sourcearrays(1,1,1,1,1),stat=ier)
    if (ier /= 0) stop 'error allocating dummy sourcearrays'
  endif

  ! ADJOINT simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! counts local receivers which become adjoint sources
    nadj_rec_local = 0

    ! gets number of local adjoint sources, i.e. located in this slice (nadj_rec_local)
    if (.not. SU_FORMAT) then
      ! prepares channel names
      do icomp=1,NDIM
        call write_channel_name(icomp,comp(icomp))
      enddo

      ! temporary counter to check if any files are found at all
      nadj_files_found = 0
      do irec = 1,nrec
        if (myrank == islice_selected_rec(irec)) then
          ! checks that the source slice number is okay
          if (islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROC-1) &
            call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')

          ! updates counter
          nadj_rec_local = nadj_rec_local + 1

          if (.not. READ_ADJSRC_ASDF) then
            ! checks **net**.**sta**.**BH**.adj files for correct number of time steps
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))

            ! loops over file components E/N/Z
            do icomp = 1,NDIM
              filename = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))// &
                       '/../SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
              open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
              if (ier == 0) then
                ! checks length of file
                itime = 0
                do while (ier == 0)
                  read(IIN,*,iostat=ier) junk,junk
                  if (ier == 0) itime = itime + 1
                enddo

                ! checks length
                if (itime /= NSTEP) then
                  print *,'adjoint source error: ',trim(filename),' has length',itime,' but should be',NSTEP
                  call exit_MPI(myrank, &
                    'file '//trim(filename)//' has wrong length, please check your adjoint sources and your simulation duration')
                endif

                ! updates counter for found files
                nadj_files_found = nadj_files_found + 1
              endif
              ! closes file
              close(IIN)
            enddo
          else
            adj_source_file = trim(network_name(irec))//'_'//trim(station_name(irec))
            do icomp = 1,NDIM
              filename = trim(adj_source_file) // '_'// comp(icomp)
              call check_adjoint_sources_asdf(filename)
            enddo
            nadj_files_found = nadj_files_found + 1
          endif
        endif
      enddo

    else
      ! SU_FORMAT file
      ! adjoint sources given in single SU_FORMAT file;
      ! skip counting, because only one file per component per proc in SU_FORMAT
      nadj_rec_local = nrec_local
      nadj_files_found = nrec_local
    endif !if (.not. SU_FORMAT)

    ! checks if any adjoint source files found at all
    call sum_all_i(nadj_files_found,nadj_files_found_tot)
    if (myrank == 0) then
      ! user output
      write(IMAIN,*)
      write(IMAIN,*) '    ',nadj_files_found_tot,' adjoint component trace files found in all slices'
      call flush_IMAIN()

      ! master process checks if any adjoint files found
      if (nadj_files_found_tot == 0) then
        print *,'Error no adjoint traces found: ',nadj_files_found_tot
        print *,'in directory : ',OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'/../SEM/'
        if (.not. SU_FORMAT) then
          print *,'with endings : ', '**.'//comp(1)//'.adj',' ','**.'//comp(2)//'.adj',' ','**.'//comp(3)//'.adj'
        endif
        print *
        call exit_MPI(myrank,'no adjoint traces found, please check adjoint sources in directory SEM/')
      endif
    endif

    ! initializes adjoint sources

    allocate(source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC),stat=ier)
    if (ier /= 0) stop 'error allocating array source_adjoint'

    ! note:
    ! computes adjoint sources in chunks/blocks during time iterations.
    ! we moved it to compute_add_sources_viscoelastic.f90 & compute_add_sources_acoustic.f90,
    ! because we may need to read in adjoint sources block by block

  endif

  ! user info
  if (USE_SOURCE_ENCODING) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'using source encoding:'
      if (does_source_encoding) then
        write(IMAIN,*) '  sources have been encoded'
      else
        write(IMAIN,*) '  source encoding has no effect (only supported for acoustic sources)'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! frees dynamically allocated memory
  deallocate(factor_force_source)
  deallocate(comp_dir_vect_source_E)
  deallocate(comp_dir_vect_source_N)
  deallocate(comp_dir_vect_source_Z_UP)

  end subroutine setup_sources_precompute_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_precompute_intp()

  use specfem_par
  implicit none

  ! local parameters
  integer :: irec,irec_local,isource,ier
  integer,dimension(0:NPROC-1) :: tmp_rec_local_all
  integer :: maxrec,maxproc(1)
  double precision :: sizeval

  ! statistics about allocation memory for seismograms & source_adjoint
  ! seismograms
  ! gather from slaves on master
  call gather_all_singlei(nrec_local,tmp_rec_local_all,NPROC)
  ! user output
  if (myrank == 0) then
    ! determines maximum number of local receivers and corresponding rank
    maxrec = maxval(tmp_rec_local_all(:))
    ! note: MAXLOC will determine the lower bound index as '1'.
    maxproc = maxloc(tmp_rec_local_all(:)) - 1
    ! seismograms array size in MB
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! seismograms need seismograms(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      sizeval = dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_OUTPUT_SEISMOS * CUSTOM_REAL / 1024. / 1024. )
    else
      ! adjoint seismograms need seismograms(NDIM*NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      sizeval = dble(maxrec) * dble(NDIM * NDIM * NTSTEP_BETWEEN_OUTPUT_SEISMOS * CUSTOM_REAL / 1024. / 1024. )
    endif
    ! outputs info
    write(IMAIN,*) 'seismograms:'
    if (WRITE_SEISMOGRAMS_BY_MASTER) then
      write(IMAIN,*) '  seismograms written by master process only'
    else
      write(IMAIN,*) '  seismograms written by all processes'
    endif
    write(IMAIN,*) '  writing out seismograms at every NTSTEP_BETWEEN_OUTPUT_SEISMOS = ',NTSTEP_BETWEEN_OUTPUT_SEISMOS
    write(IMAIN,*) '  maximum number of local receivers is ',maxrec,' in slice ',maxproc(1)
    write(IMAIN,*) '  size of maximum seismogram array       = ', sngl(sizeval),'MB'
    write(IMAIN,*) '                                         = ', sngl(sizeval/1024.d0),'GB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! gather from slaves on master
    call gather_all_singlei(nadj_rec_local,tmp_rec_local_all,NPROC)
    ! user output
    if (myrank == 0) then
      ! determines maximum number of local receivers and corresponding rank
      maxrec = maxval(tmp_rec_local_all(:))
      ! note: MAXLOC will determine the lower bound index as '1'.
      maxproc = maxloc(tmp_rec_local_all(:)) - 1
      ! source_adjoint size in MB
      !source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC)
      sizeval = dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_READ_ADJSRC * CUSTOM_REAL / 1024. / 1024. )
      ! outputs info
      write(IMAIN,*) 'adjoint source arrays:'
      write(IMAIN,*) '  reading adjoint sources at every NTSTEP_BETWEEN_READ_ADJSRC = ',NTSTEP_BETWEEN_READ_ADJSRC
      write(IMAIN,*) '  maximum number of local adjoint sources is ',maxrec,' in slice ',maxproc(1)
      write(IMAIN,*) '  size of maximum adjoint source array = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif


  ! needs to be allocate for subroutine calls (even if nrec_local == 0)
  allocate(number_receiver_global(nrec_local),stat=ier)
  if (ier /= 0) stop 'error allocating array number_receiver_global'

  ! stores local receivers interpolation factors
  if (nrec_local > 0) then
    ! allocate Lagrange interpolators for receivers
    allocate(hxir_store(nrec_local,NGLLX), &
             hetar_store(nrec_local,NGLLY), &
             hgammar_store(nrec_local,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating array hxir_store etc.'

    ! allocates derivatives
    if (SIMULATION_TYPE == 2) then
      allocate(hpxir_store(nrec_local,NGLLX), &
               hpetar_store(nrec_local,NGLLY), &
               hpgammar_store(nrec_local,NGLLZ),stat=ier)
      if (ier /= 0) stop 'error allocating array hpxir_store'
    endif

    ! define local to global receiver numbering mapping
    irec_local = 0
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      do irec = 1,nrec
        if (myrank == islice_selected_rec(irec)) then
          irec_local = irec_local + 1
          number_receiver_global(irec_local) = irec
        endif
      enddo
    else
      ! note: pure adjoint runs (SIMULATION_TYPE == 2):
      !       uses CMT source locations to store "receiver" seismograms
      do isource = 1,NSOURCES
        if (myrank == islice_selected_source(isource)) then
          irec_local = irec_local + 1
          number_receiver_global(irec_local) = isource
        endif
      enddo
    endif

    ! define and store Lagrange interpolators at all the receivers
    do irec_local = 1,nrec_local
      irec = number_receiver_global(irec_local)

      if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
        ! receiver positions
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
        call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
        call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
      else
        ! source positions
        ! note: pure adjoint runs (SIMULATION_TYPE == 2):
        !       uses CMT source locations to store "receiver" seismograms
        call lagrange_any(xi_source(irec),NGLLX,xigll,hxir,hpxir)
        call lagrange_any(eta_source(irec),NGLLY,yigll,hetar,hpetar)
        call lagrange_any(gamma_source(irec),NGLLZ,zigll,hgammar,hpgammar)
      endif

      ! stores interpolators
      hxir_store(irec_local,:) = hxir(:)
      hetar_store(irec_local,:) = hetar(:)
      hgammar_store(irec_local,:) = hgammar(:)

      ! stores derivatives
      if (SIMULATION_TYPE == 2) then
        hpxir_store(irec_local,:) = hpxir(:)
        hpetar_store(irec_local,:) = hpetar(:)
        hpgammar_store(irec_local,:) = hpgammar(:)
      endif
    enddo
  else
    ! VM VM need to allocate Lagrange interpolators for receivers with 0 because it is used
    ! in calling subroutines parmeters. (otherwise it can be crash at runtime).
    allocate(hxir_store(nrec_local,NGLLX), &
             hetar_store(nrec_local,NGLLY), &
             hgammar_store(nrec_local,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating array hxir_store etc.'
    ! allocates derivatives
    if (SIMULATION_TYPE == 2) then
      allocate(hpxir_store(nrec_local,NGLLX), &
               hpetar_store(nrec_local,NGLLY), &
               hpgammar_store(nrec_local,NGLLZ),stat=ier)
      if (ier /= 0) stop 'error allocating array hpxir_store'
    endif
  endif ! nrec_local > 0

  ! seismograms
  if (nrec_local > 0) then
    ! allocate seismogram array
    if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
      allocate(seismograms_d(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
    else
      allocate(seismograms_d(1,1,1),stat=ier)
    endif
    if (ier /= 0) stop 'error allocating array seismograms_d'

    if (SAVE_SEISMOGRAMS_VELOCITY) then
      allocate(seismograms_v(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
    else
      allocate(seismograms_v(1,1,1),stat=ier)
    endif
    if (ier /= 0) stop 'error allocating array seismograms_v'

    if (SAVE_SEISMOGRAMS_ACCELERATION) then
      allocate(seismograms_a(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
    else
      allocate(seismograms_a(1,1,1),stat=ier)
    endif
    if (ier /= 0) stop 'error allocating array seismograms_a'

    if (SAVE_SEISMOGRAMS_PRESSURE) then
      !NB_RUNS_ACOUSTIC_GPU is set to 1 by default in constants.h
      allocate(seismograms_p(NDIM,nrec_local*NB_RUNS_ACOUSTIC_GPU,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
    else
      allocate(seismograms_p(1,1,1),stat=ier)
    endif
    if (ier /= 0) stop 'error allocating array seismograms_p'

    ! initialize seismograms
    seismograms_d(:,:,:) = 0._CUSTOM_REAL
    seismograms_v(:,:,:) = 0._CUSTOM_REAL
    seismograms_a(:,:,:) = 0._CUSTOM_REAL
    seismograms_p(:,:,:) = 0._CUSTOM_REAL
  endif

  ! seismograms
  if (SIMULATION_TYPE == 2) then
    if (nrec_local > 0) then
      allocate(seismograms_eps(NDIM,NDIM,nrec_local,NSTEP),stat=ier)
      if (ier /= 0) stop 'error allocating array seismograms_eps'
      seismograms_eps(:,:,:,:) = 0._CUSTOM_REAL
    endif
  endif
  ! adjoint seismograms writing
  it_adj_written = 0


  end subroutine setup_receivers_precompute_intp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_VTKfile()

  use specfem_par
  implicit none

  double precision :: shape3D(NGNOD)
  double precision :: xil,etal,gammal
  double precision :: xmesh,ymesh,zmesh
  real(kind=CUSTOM_REAL),dimension(NGNOD) :: xelm,yelm,zelm
  integer :: ia,ispec,isource,irec,ier,totalpoints

  character(len=MAX_STRING_LEN) :: filename,filename_new,system_command,system_command1,system_command2

  ! determines number of points for vtk file
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    totalpoints = NSOURCES + nrec
  else
    ! pure adjoint simulation only needs receivers
    totalpoints = nrec
  endif

  if (myrank == 0) then
    ! vtk file
    open(IOVTK,file=trim(OUTPUT_FILES)//'/sr.vtk',status='unknown',iostat=ier)
    if (ier /= 0) stop 'error opening sr.vtk file'
    ! vtk header
    write(IOVTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOVTK,'(a)') 'Source and Receiver VTK file'
    write(IOVTK,'(a)') 'ASCII'
    write(IOVTK,'(a)') 'DATASET POLYDATA'
    write(IOVTK, '(a,i6,a)') 'POINTS ', totalpoints, ' float'
  endif

  ! sources
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do isource=1,NSOURCES
      ! spectral element id
      ispec = ispec_selected_source(isource)

      ! gets element ancor nodes
      if (myrank == islice_selected_source(isource)) then
        ! find the coordinates of the eight corner nodes of the element
        call eval_shape3D_element_corners(xelm,yelm,zelm,ispec, &
                        ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)

      endif
      ! master collects corner locations
      if (islice_selected_source(isource) /= 0) then
        if (myrank == 0) then
          call recvv_cr(xelm,NGNOD,islice_selected_source(isource),0)
          call recvv_cr(yelm,NGNOD,islice_selected_source(isource),0)
          call recvv_cr(zelm,NGNOD,islice_selected_source(isource),0)
        else if (myrank == islice_selected_source(isource)) then
          call sendv_cr(xelm,NGNOD,0,0)
          call sendv_cr(yelm,NGNOD,0,0)
          call sendv_cr(zelm,NGNOD,0,0)
        endif
      endif

      if (myrank == 0) then
        ! get the 3-D shape functions
        xil = xi_source(isource)
        etal = eta_source(isource)
        gammal = gamma_source(isource)
        call eval_shape3D_single(myrank,shape3D,xil,etal,gammal,NGNOD)

        ! interpolates source locations
        xmesh = 0.0
        ymesh = 0.0
        zmesh = 0.0
        do ia=1,NGNOD
          xmesh = xmesh + shape3D(ia)*xelm(ia)
          ymesh = ymesh + shape3D(ia)*yelm(ia)
          zmesh = zmesh + shape3D(ia)*zelm(ia)
        enddo

        ! writes out to VTK file
        write(IOVTK,'(3e18.6)') xmesh,ymesh,zmesh
      endif
    enddo ! NSOURCES
  endif

  ! receivers
  do irec=1,nrec
    ispec = ispec_selected_rec(irec)

    ! find the coordinates of the eight corner nodes of the element
    if (myrank == islice_selected_rec(irec)) then
      call eval_shape3D_element_corners(xelm,yelm,zelm,ispec, &
                      ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)
    endif
    ! master collects corner locations
    if (islice_selected_rec(irec) /= 0) then
      if (myrank == 0) then
        call recvv_cr(xelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(yelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(zelm,NGNOD,islice_selected_rec(irec),0)
      else if (myrank == islice_selected_rec(irec)) then
        call sendv_cr(xelm,NGNOD,0,0)
        call sendv_cr(yelm,NGNOD,0,0)
        call sendv_cr(zelm,NGNOD,0,0)
      endif
    endif

    if (myrank == 0) then
      ! get the 3-D shape functions
      xil = xi_receiver(irec)
      etal = eta_receiver(irec)
      gammal = gamma_receiver(irec)
      call eval_shape3D_single(myrank,shape3D,xil,etal,gammal,NGNOD)

      ! interpolates receiver locations
      xmesh = 0.0
      ymesh = 0.0
      zmesh = 0.0
      do ia=1,NGNOD
        xmesh = xmesh + shape3D(ia)*xelm(ia)
        ymesh = ymesh + shape3D(ia)*yelm(ia)
        zmesh = zmesh + shape3D(ia)*zelm(ia)
      enddo

      ! writes out to VTK file
      write(IOVTK,'(3e18.6)') xmesh,ymesh,zmesh
    endif
  enddo

  ! closes vtk file
  if (myrank == 0) then
    write(IOVTK,*)
    close(IOVTK)

    ! creates additional receiver and source files
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! extracts receiver locations
      filename = trim(OUTPUT_FILES)//'/sr.vtk'
      filename_new = trim(OUTPUT_FILES)//'/receiver.vtk'

      ! vtk file for receivers only
      write(system_command, &
  "('awk ',a1,'{if (NR < 5) print $0;if (NR == 5)print ',a1,'POINTS',i6,' float',a1,';if &
      &(NR > 5+',i6,')print $0}',a1,' < ',a,' > ',a)")&
      "'",'"',nrec,'"',NSOURCES,"'",trim(filename),trim(filename_new)

      ! extracts source locations
      filename_new = trim(OUTPUT_FILES)//'/source.vtk'

      write(system_command1, &
  "('awk ',a1,'{if (NR < 5) print $0;if (NR == 5)print ',a1,'POINTS',i6,' float',a1,';')") &
        "'",'"',NSOURCES,'"'

      write(system_command2, &
  "('if (NR > 5 && NR < 6+',i6,')print $0}END{print ',a,'}',a1,' < ',a,' > ',a)") &
        NSOURCES,'" "',"'",trim(filename),trim(filename_new)

      system_command = trim(system_command1)//trim(system_command2)

    endif
  endif

  end subroutine setup_sources_receivers_VTKfile

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_search_kdtree()

  use specfem_par

  use kdtree_search, only: kdtree_setup,kdtree_set_verbose,kdtree_delete,kdtree_find_nearest_neighbor, &
    kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

  implicit none
  integer :: ier,ispec,iglob

  ! kdtree search

  ! set number of tree nodes
  kdtree_num_nodes = NSPEC_AB

  ! allocates tree arrays
  allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
  if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
  allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
  if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

  ! tree verbosity
  if (myrank == 0) call kdtree_set_verbose()

  ! prepares search arrays, each element takes its internal GLL points for tree search
  kdtree_nodes_index(:) = 0
  kdtree_nodes_location(:,:) = 0.0

  ! fills kd-tree arrays
  do ispec = 1,NSPEC_AB
    ! adds node index ( index points to same ispec for all internal GLL points)
    kdtree_nodes_index(ispec) = ispec

    ! adds node location (of midpoint)
    iglob = ibool(MIDX,MIDY,MIDZ,ispec)
    kdtree_nodes_location(1,ispec) = xstore(iglob)
    kdtree_nodes_location(2,ispec) = ystore(iglob)
    kdtree_nodes_location(3,ispec) = zstore(iglob)
  enddo

  ! creates kd-tree for searching
  ! serial way
  !do i = 0,NPROCTOT_VAL-1
  !  if (myrank == i) then
  !    print *,'kd-tree setup for process: ',myrank
  !    call kdtree_setup()
  !  endif
  !  call synchronize_all()
  !enddo
  ! parallel way
  call kdtree_setup()

  ! synchronizes all mpi-processes
  call synchronize_all()

  end subroutine setup_search_kdtree
