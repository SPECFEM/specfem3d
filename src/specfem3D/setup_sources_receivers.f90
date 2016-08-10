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

  subroutine setup_sources_receivers()

  use specfem_par

  implicit none

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
    if (NSOURCES > 1) write(IMAIN,*) 'Using ',NSOURCES,' point sources'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine setup_sources_receivers

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

  double precision :: t0_acoustic,min_tshift_src_original
  integer :: yr,jda,ho,mi
  integer :: isource,ispec,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources:'
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
           nu_source(3,3,NSOURCES), stat=ier)
  if (ier /= 0) stop 'error allocating arrays for sources'

  if (USE_FORCE_POINT_SOURCE) then
    allocate(factor_force_source(NSOURCES), &
             comp_dir_vect_source_E(NSOURCES), &
             comp_dir_vect_source_N(NSOURCES), &
             comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
    if (ier /= 0) stop 'error allocating arrays for force point sources'
  endif

  !! VM VM set the size of user_source_time_function
  if (EXTERNAL_STF) then
     NSTEP_STF = NSTEP
     NSOURCES_STF = NSOURCES
  else !! We don't need the array user_source_time_function : use a small dummy array
     NSTEP_STF = 1
     NSOURCES_STF = 1
  endif
  !! allocate the array contains the user defined source time function
  allocate(user_source_time_function(NSTEP_STF, NSOURCES_STF),stat=ier)
  if (ier /= 0) stop 'error allocating arrays for user sources time function'

! locate sources in the mesh
!
! returns:  islice_selected_source & ispec_selected_source,
!                xi_source, eta_source & gamma_source
  call locate_source(ibool,NSOURCES,myrank,NSPEC_AB,NGLOB_AB,NGNOD, &
          xstore,ystore,zstore,xigll,yigll,zigll,NPROC, &
          tshift_src,min_tshift_src_original,yr,jda,ho,mi,utm_x_source,utm_y_source, &
          DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
          islice_selected_source,ispec_selected_source, &
          xi_source,eta_source,gamma_source, &
          nu_source,iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
          ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
          num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  if (abs(minval(tshift_src)) > TINYVAL) call exit_MPI(myrank,'one tshift_src must be zero, others must be positive')

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
  if (EXTERNAL_STF) then
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

  ! count number of sources located in this slice
  nsources_local = 0
  do isource = 1, NSOURCES
    if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
  enddo

  ! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there
  call setup_sources_check_acoustic()

  ! prints source time functions to output files
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  end subroutine setup_sources

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

  real(kind=CUSTOM_REAL):: minl,maxl,min_all,max_all
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! gets model dimensions
  minl = minval( xstore )
  maxl = maxval( xstore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LONGITUDE_MIN = min_all
  LONGITUDE_MAX = max_all

  minl = minval( ystore )
  maxl = maxval( ystore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LATITUDE_MIN = min_all
  LATITUDE_MAX = max_all

  ! reads in station file
  if (SIMULATION_TYPE == 1) then
    rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS'
    filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_FILTERED'

    ! see if we are running several independent runs in parallel
    ! if so, add the right directory for that run
    ! (group numbers start at zero, but directory names start at run0001, thus we add one)
    ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      rec_filename = path_to_add(1:len_trim(path_to_add))//rec_filename(1:len_trim(rec_filename))
      filtered_rec_filename = path_to_add(1:len_trim(path_to_add))//filtered_rec_filename(1:len_trim(filtered_rec_filename))
    endif

    call station_filter(SUPPRESS_UTM_PROJECTION,UTM_PROJECTION_ZONE,myrank,rec_filename,filtered_rec_filename,nrec, &
                        LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)

    if (nrec < 1) call exit_MPI(myrank,'need at least one receiver')
    call synchronize_all()

  else
    rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT'
    filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT_FILTERED'

    ! see if we are running several independent runs in parallel
    ! if so, add the right directory for that run
    ! (group numbers start at zero, but directory names start at run0001, thus we add one)
    ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      rec_filename = path_to_add(1:len_trim(path_to_add))//rec_filename(1:len_trim(rec_filename))
      filtered_rec_filename = path_to_add(1:len_trim(path_to_add))//filtered_rec_filename(1:len_trim(filtered_rec_filename))
    endif

    call station_filter(SUPPRESS_UTM_PROJECTION,UTM_PROJECTION_ZONE,myrank,rec_filename,filtered_rec_filename,nrec, &
                        LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)
    if (nrec < 1) call exit_MPI(myrank, 'adjoint simulation needs at least one receiver')
    call synchronize_all()
  endif

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

  if (nrec < 1) call exit_MPI(myrank,'need at least one receiver')

  ! allocate memory for receiver arrays, i.e. stations given in STATIONS file
  allocate(islice_selected_rec(nrec), &
           ispec_selected_rec(nrec), &
           xi_receiver(nrec), &
           eta_receiver(nrec), &
           gamma_receiver(nrec), &
           station_name(nrec), &
           network_name(nrec), &
           nu(NDIM,NDIM,nrec),stat=ier)
  if (ier /= 0) stop 'error allocating arrays for receivers'

  ! locate receivers in the mesh
  call locate_receivers(ibool,myrank,NSPEC_AB,NGLOB_AB,NGNOD, &
                        xstore,ystore,zstore,xigll,yigll,zigll,filtered_rec_filename, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                        NPROC,utm_x_source(1),utm_y_source(1), &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                        iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
                        num_free_surface_faces,free_surface_ispec,free_surface_ijk,SU_FORMAT)

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
  integer :: i,j,k
  integer :: icomp,itime,nadj_files_found,nadj_files_found_tot,ier

  character(len=3),dimension(NDIM) :: comp
  character(len=MAX_STRING_LEN) :: filename

  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision :: hlagrange
  double precision :: norm,comp_x,comp_y,comp_z

  logical :: does_source_encoding

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

        if (USE_FORCE_POINT_SOURCE) then ! use of FORCESOLUTION files

          ! note: for use_force_point_source xi/eta/gamma are also in the range [-1,1], for exact positioning

          ! initializes source array
          sourcearrayd(:,:,:,:) = 0.0d0

          ! calculates source array for interpolated location
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                hlagrange = hxis(i) * hetas(j) * hgammas(k)

                ! acoustic source
                ! identical source array components in x,y,z-direction
                if (ispec_is_acoustic(ispec)) then
                  sourcearrayd(:,i,j,k) = factor_force_source(isource) * hlagrange
                endif

                ! elastic source
                if (ispec_is_elastic(ispec)) then
                  ! length of component vector
                  norm = sqrt( comp_dir_vect_source_E(isource)**2 &
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

                  ! we use an tilted force defined by its magnitude and the projections
                  ! of an arbitrary (non-unitary) direction vector on the E/N/Z_UP basis
                  sourcearrayd(:,i,j,k) = factor_force_source(isource) * hlagrange * &
                                          ( nu_source(1,:,isource) * comp_x + &
                                            nu_source(2,:,isource) * comp_y + &
                                            nu_source(3,:,isource) * comp_z )

                endif
              enddo
            enddo
          enddo

          ! distinguish between single and double precision for reals
          sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:),kind=CUSTOM_REAL)

        else ! use of CMTSOLUTION files

          ! elastic or poroelastic moment tensor source
          if (ispec_is_elastic(ispec) .or. ispec_is_poroelastic(ispec)) then
            call compute_arrays_source(ispec,sourcearray, &
                                       Mxx(isource),Myy(isource),Mzz(isource), &
                                       Mxy(isource),Mxz(isource),Myz(isource), &
                                       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                       hxis,hpxis,hetas,hpetas,hgammas,hpgammas,NSPEC_AB)
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

            ! source array interpolated on all element GLL points
            call compute_arrays_source_acoustic(sourcearray,hxis,hetas,hgammas,factor_source)
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
    allocate(adj_sourcearrays(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating array adj_sourcearrays'
    adj_sourcearrays(:,:,:,:,:,:) = 0._CUSTOM_REAL

    ! note:
    ! computes adjoint sources in chunks/blocks during time iterations.
    ! we moved it to compute_add_sources_viscoelastic.f90 & compute_add_sources_acoustic.f90,
    ! because we may need to read in adjoint sources block by block

  else
    ! (SIMULATION_TYPE == 1)
    ! allocate dummy array in order to be able to use it as a subroutine argument, even if unused
    nadj_rec_local = 0
    NTSTEP_BETWEEN_READ_ADJSRC = 0
    allocate(adj_sourcearrays(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating dummy array adj_sourcearrays'
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

  end subroutine setup_sources_precompute_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_precompute_intp()

  use specfem_par
  implicit none

  integer :: irec,irec_local,isource,ier

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
  endif ! nrec_local > 0

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

  subroutine print_stf_file()

  use specfem_par
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: stf_used,time_source
  real(kind=CUSTOM_REAL),dimension(NSTEP) :: source_time_function

  double precision :: stf,time_source_dble
  double precision,external :: get_stf_acoustic,get_stf_viscoelastic,get_stf_poroelastic

  integer :: isource,ispec,ier
  character(len=MAX_STRING_LEN) :: plot_file

  ! check
  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    stop 'unrecognized SIMULATION_TYPE value in printing stf file'

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'printing the source-time function'
    call flush_IMAIN()
  endif

  ! note: the source time function will be output for each source separately,
  !       instead of summing up all stf from single source contributions

  ! source contributions
  do isource = 1,NSOURCES

    ! initializes
    source_time_function(:) = 0._CUSTOM_REAL

    ! compute the source contribution (only if this proc carries the source)
    if (myrank == islice_selected_source(isource)) then

      ! time loop
      do it = 1,NSTEP
        ! current source time
        if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 2) then
          if (USE_LDDRK) then
            time_source_dble = dble(it-1)*DT + dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(it-1)*DT - t0 - tshift_src(isource)
          endif
        else
          ! backward simulation (SIMULATION_TYPE == 3)
          if (USE_LDDRK) then
            time_source_dble = dble(NSTEP-1)*DT - dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(NSTEP-1)*DT - t0 - tshift_src(isource)
          endif
        endif

        ispec = ispec_selected_source(isource)

        ! determines source time function value
        if (ispec_is_acoustic(ispec)) then
          stf = get_stf_acoustic(time_source_dble,isource)
        else if (ispec_is_elastic(ispec)) then
          stf = get_stf_viscoelastic(time_source_dble,isource)
        else if (ispec_is_poroelastic(ispec)) then
          stf = get_stf_poroelastic(time_source_dble,isource)
        else
          call exit_MPI(myrank,'Invalid source element type, please check your mesh...')
        endif

        !! VM VM add external source time function
        if (EXTERNAL_STF) stf = user_source_time_function(it, isource)

        ! distinguishes between single and double precision for reals
        stf_used = real(stf,kind=CUSTOM_REAL)

        ! for file output
        source_time_function(it) = stf_used
      enddo
    endif

    ! master collects stf (if it does not already have it, i.e. if this source is not on the master)
    if (islice_selected_source(isource) /= 0) then
      if (myrank == 0) then
        ! master collects
        call recvv_cr(source_time_function,NSTEP,islice_selected_source(isource),0)
      else if (myrank == islice_selected_source(isource)) then
        ! slave sends to master
        call sendv_cr(source_time_function,NSTEP,0,0)
      endif
    endif

    ! master prints out to file
    if (myrank == 0) then
      ! opens source time function file
      if (NSOURCES == 1) then
        plot_file = '/plot_source_time_function.txt'
      else if (isource < 10) then
        write(plot_file,"('/plot_source_time_function',i1,'.txt')") isource
      else if (isource < 100) then
        write(plot_file,"('/plot_source_time_function',i2,'.txt')") isource
      else
        write(plot_file,"('/plot_source_time_functionA',i7.7,'.txt')") isource
      endif
      open(unit=IOSTF,file=trim(OUTPUT_FILES)//trim(plot_file),status='unknown',iostat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error opening plot_source_time_function file')

      do it = 1,NSTEP
        ! overall time, note that t_shift_src will start at zero for simulation
        if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 2) then
          if (USE_LDDRK) then
            time_source = dble(it-1)*DT + dble(C_LDDRK(istage))*DT - t0
          else
            time_source = dble(it-1)*DT - t0
          endif
        else
          ! backward simulation (SIMULATION_TYPE == 3)
          if (USE_LDDRK) then
            time_source = dble(NSTEP-1)*DT - dble(C_LDDRK(istage))*DT - t0
          else
            time_source = dble(NSTEP-1)*DT - t0
          endif
        endif

        ! file output
        write(IOSTF,*) time_source,source_time_function(it)
      enddo

      close(IOSTF)

    endif ! of if (myrank == 0) then

  enddo ! NSOURCES

  end subroutine print_stf_file

