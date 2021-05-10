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

!----
!----  locate_source finds the correct position of the source
!----

  subroutine locate_source(filename,tshift_src,min_tshift_src_original,utm_x_source,utm_y_source, &
                           hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                           islice_selected_source,ispec_selected_source, &
                           factor_force_source,comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP, &
                           xi_source,eta_source,gamma_source,nu_source,user_source_time_function)

  use constants

  use specfem_par, only: USE_FORCE_POINT_SOURCE,USE_RICKER_TIME_FUNCTION, &
      UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,USE_SOURCES_RECEIVERS_Z, &
      NSTEP_STF,NSOURCES_STF,USE_EXTERNAL_SOURCE_FILE, &
      USE_TRICK_FOR_BETTER_PRESSURE,COUPLE_WITH_INJECTION_TECHNIQUE, &
      myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,DT, &
      NSOURCES

  ! PML
  use pml_par, only: is_CPML

  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: filename
  double precision,dimension(NSOURCES),intent(inout) :: tshift_src
  double precision,intent(inout) :: min_tshift_src_original
  double precision, dimension(NSOURCES),intent(inout) :: utm_x_source,utm_y_source

  ! CMTS
  double precision, dimension(NSOURCES),intent(inout) :: hdur
  double precision, dimension(NSOURCES),intent(inout) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  integer, dimension(NSOURCES), intent(inout) :: islice_selected_source,ispec_selected_source

  ! force
  double precision, dimension(NSOURCES) :: factor_force_source
  double precision, dimension(NSOURCES) :: comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP

  double precision, dimension(NSOURCES),intent(inout) :: xi_source,eta_source,gamma_source
  double precision, dimension(NDIM,NDIM,NSOURCES),intent(out) :: nu_source
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF,NSOURCES_STF) :: user_source_time_function

  ! local parameters
  ! sources
  integer :: isource
  double precision :: f0,t0_ricker
  ! CMTs
  double precision, dimension(:), allocatable :: lat,long,depth
  double precision, dimension(:,:), allocatable ::  moment_tensor
  ! positioning
  double precision, dimension(:), allocatable :: x_found,y_found,z_found
  double precision, dimension(:), allocatable :: elevation
  double precision, dimension(:), allocatable :: final_distance
  double precision, dimension(:), allocatable :: x_target,y_target,z_target
  integer, dimension(:), allocatable :: idomain

  double precision, external :: get_cmt_scalar_moment
  double precision, external :: get_cmt_moment_magnitude

  ! location search
  integer :: ispec_found,idomain_found
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  double precision :: x,y,z,x_new,y_new,z_new
  double precision :: xi,eta,gamma,final_distance_squared
  double precision, dimension(NDIM,NDIM) :: nu_found

  ! subset arrays
  integer :: nsources_subset_current_size,isource_in_this_subset,isources_already_done
  integer, dimension(NSOURCES_SUBSET_MAX) :: ispec_selected_source_subset,idomain_subset
  double precision, dimension(NSOURCES_SUBSET_MAX) :: xi_source_subset,eta_source_subset,gamma_source_subset
  double precision, dimension(NSOURCES_SUBSET_MAX) :: x_found_subset,y_found_subset,z_found_subset
  double precision, dimension(NSOURCES_SUBSET_MAX) :: final_distance_subset
  double precision, dimension(NDIM,NDIM,NSOURCES_SUBSET_MAX) :: nu_subset

  logical,dimension(:),allocatable :: is_CPML_source,is_CPML_source_all

  integer :: ispec,ier
  logical :: is_done_sources

  ! timer MPI
  double precision :: tstart,tCPU
  double precision, external :: wtime

  !-----------------------------------------------------------------------------------

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating sources'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,'(1x,a,a,a)') 'reading source information from ', trim(filename), ' file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocates temporary arrays
  allocate(moment_tensor(6,NSOURCES), &
           lat(NSOURCES), &
           long(NSOURCES), &
           depth(NSOURCES), &
           x_found(NSOURCES), &
           y_found(NSOURCES), &
           z_found(NSOURCES), &
           elevation(NSOURCES), &
           final_distance(NSOURCES), &
           x_target(NSOURCES), &
           y_target(NSOURCES), &
           z_target(NSOURCES), &
           idomain(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating source arrays')
  moment_tensor(:,:) = 0.d0
  lat(:) = 0.d0; long(:) = 0.d0; depth(:) = 0.d0
  x_found(:) = 0.d0; y_found(:) = 0.d0; z_found(:) = 0.d0
  x_target(:) = 0.d0; y_target(:) = 0.d0; z_target(:) = 0.d0
  elevation(:) = 0.d0
  final_distance(:) = HUGEVAL
  idomain(:) = 0

  ! clear the arrays
  Mxx(:) = 0.d0
  Myy(:) = 0.d0
  Mzz(:) = 0.d0
  Mxy(:) = 0.d0
  Mxz(:) = 0.d0
  Myz(:) = 0.d0

  ! read all the sources
  if (USE_FORCE_POINT_SOURCE) then
    ! point forces
    if (myrank == 0) then
      ! only main process reads in FORCESOLUTION file
      call get_force(filename,tshift_src,hdur,lat,long,depth,NSOURCES,min_tshift_src_original,factor_force_source, &
                     comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP, &
                     user_source_time_function)
    endif
    ! broadcasts specific point force infos
    call bcast_all_dp(factor_force_source,NSOURCES)
    call bcast_all_dp(comp_dir_vect_source_E,NSOURCES)
    call bcast_all_dp(comp_dir_vect_source_N,NSOURCES)
    call bcast_all_dp(comp_dir_vect_source_Z_UP,NSOURCES)
  else
    ! CMT moment tensors
    if (myrank == 0) then
      ! only main process reads in CMTSOLUTION file
      call get_cmt(filename,tshift_src,hdur,lat,long,depth,moment_tensor, &
                   DT,NSOURCES,min_tshift_src_original,user_source_time_function)
    endif
    ! broadcasts specific moment tensor infos
    call bcast_all_dp(moment_tensor,6*NSOURCES)
  endif

  ! broadcasts general source information read on the main to the nodes
  call bcast_all_dp(tshift_src,NSOURCES)
  call bcast_all_dp(hdur,NSOURCES)
  call bcast_all_dp(lat,NSOURCES)
  call bcast_all_dp(long,NSOURCES)
  call bcast_all_dp(depth,NSOURCES)
  call bcast_all_singledp(min_tshift_src_original)
  call bcast_all_cr(user_source_time_function,NSOURCES_STF*NSTEP_STF)

  ! compute typical size of elements
  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  ! get MPI starting time
  tstart = wtime()

  ! user output
  if (myrank == 0) then
    if (SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,*) 'no UTM projection'
    else
      write(IMAIN,*) 'UTM projection:'
      write(IMAIN,*) '  UTM zone: ',UTM_PROJECTION_ZONE
    endif
    if (USE_SOURCES_RECEIVERS_Z) then
      write(IMAIN,*) 'using sources/receivers Z:'
      write(IMAIN,*) '  (depth) becomes directly (z) coordinate'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! determines target point locations (need to locate z coordinate of all sources)
  ! note: we first locate all the target positions in the mesh to reduces the need of MPI communication
  if (.not. USE_SOURCES_RECEIVERS_Z) then
    ! converts km to m
    do isource = 1,NSOURCES
      depth(isource) = depth(isource)*1000.0d0
    enddo
  endif
  call get_elevation_and_z_coordinate_all(NSOURCES,long,lat,depth,utm_x_source,utm_y_source,elevation, &
                                          x_target,y_target,z_target)
  !debug
  !print *,'source elevations:',elevation

  ! reference frame convertion:
  !   Harvard CMT convention: r is up, t is south, and p is east
  !                           (x = South, y = East, z = Up)
  !
  !   SPECFEM reference frame convention: x = East, y = North, z = Up
  !
  !   to convert from CMT to SPECFEM: r -> z, theta -> -y, phi -> x
  !
  ! moment_tensor(1,:) = Mrr =  Mzz
  ! moment_tensor(2,:) = Mtt =  Myy
  ! moment_tensor(3,:) = Mpp =  Mxx
  ! moment_tensor(4,:) = Mrt = -Myz
  ! moment_tensor(5,:) = Mrp =  Mxz
  ! moment_tensor(6,:) = Mtp = -Mxy

  ! get the moment tensor
  Mzz(:) = + moment_tensor(1,:) ! Mrr
  Mxx(:) = + moment_tensor(3,:) ! Mpp
  Myy(:) = + moment_tensor(2,:) ! Mtt
  Mxz(:) = + moment_tensor(5,:) ! Mrp
  Myz(:) = - moment_tensor(4,:) ! - Mrt
  Mxy(:) = - moment_tensor(6,:) ! - Mtp

  ! loop on all the sources
  do isources_already_done = 0, NSOURCES, NSOURCES_SUBSET_MAX

    ! the size of the subset can be the maximum size, or less (if we are in the last subset,
    ! or if there are fewer sources than the maximum size of a subset)
    nsources_subset_current_size = min(NSOURCES_SUBSET_MAX, NSOURCES - isources_already_done)

    ! loop over sources within this subset
    do isource_in_this_subset = 1,nsources_subset_current_size

      ! mapping from source number in current subset to real source number in all the subsets
      isource = isource_in_this_subset + isources_already_done

      x = x_target(isource)
      y = y_target(isource)
      z = z_target(isource)

      ! locates point in mesh
      call locate_point_in_mesh(x, y, z, &
                                SOURCES_CAN_BE_BURIED, elemsize_max_glob, &
                                ispec_found, xi, eta, gamma, &
                                x_new, y_new, z_new, &
                                idomain_found, nu_found, final_distance_squared)

      ispec_selected_source_subset(isource_in_this_subset) = ispec_found

      x_found_subset(isource_in_this_subset) = x_new
      y_found_subset(isource_in_this_subset) = y_new
      z_found_subset(isource_in_this_subset) = z_new

      xi_source_subset(isource_in_this_subset) = xi
      eta_source_subset(isource_in_this_subset) = eta
      gamma_source_subset(isource_in_this_subset) = gamma

      idomain_subset(isource_in_this_subset) = idomain_found
      nu_subset(:,:,isource_in_this_subset) = nu_found(:,:)
      final_distance_subset(isource_in_this_subset) = sqrt(final_distance_squared)

      ! user output progress
      if (myrank == 0 .and. NSOURCES > 1000) then
        if (mod(isource,500) == 0) then
          write(IMAIN,*) '  located source ',isource,'out of',NSOURCES
          call flush_IMAIN()
        endif
      endif
    enddo ! loop over subset

    ! main process locates best location in all slices
    call locate_MPI_slice(nsources_subset_current_size,isources_already_done, &
                          ispec_selected_source_subset, &
                          x_found_subset, y_found_subset, z_found_subset, &
                          xi_source_subset,eta_source_subset,gamma_source_subset, &
                          idomain_subset,nu_subset,final_distance_subset, &
                          NSOURCES,ispec_selected_source, islice_selected_source, &
                          x_found,y_found,z_found, &
                          xi_source, eta_source, gamma_source, &
                          idomain,nu_source,final_distance)

  enddo ! end of loop on all the sources

  ! bcast from main process
  call bcast_all_i(islice_selected_source,NSOURCES)
  call bcast_all_i(idomain,NSOURCES)
  call bcast_all_i(ispec_selected_source,NSOURCES)

  call bcast_all_dp(xi_source,NSOURCES)
  call bcast_all_dp(eta_source,NSOURCES)
  call bcast_all_dp(gamma_source,NSOURCES)

  call bcast_all_dp(x_found,NSOURCES)
  call bcast_all_dp(y_found,NSOURCES)
  call bcast_all_dp(z_found,NSOURCES)

  call bcast_all_dp(nu_source,NDIM*NDIM*NSOURCES)
  call bcast_all_dp(final_distance,NSOURCES)

  ! warning if source in C-PML region
  allocate(is_CPML_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating is_CPML_source array')
  if (myrank == 0) then
    ! only main collects
    allocate(is_CPML_source_all(NSOURCES),stat=ier)
  else
    ! dummy
    allocate(is_CPML_source_all(1),stat=ier)
  endif
  ! sets flag if source element in PML
  is_CPML_source(:) = .false.
  do isource = 1,NSOURCES
    if (islice_selected_source(isource) == myrank) then
      ispec = ispec_selected_source(isource)
      if (is_CPML(ispec)) then
        is_CPML_source(isource) = .true.
        ! debug
        !print *,'Warning: rank ',myrank,' has source ', isource,' in C-PML region'
      endif
    endif
  enddo
  call any_all_1Darray_l(is_CPML_source,is_CPML_source_all,NSOURCES)


  ! sets new utm coordinates for best locations
  utm_x_source(:) = x_found(:)
  utm_y_source(:) = y_found(:)

  ! output source information to a file so that we can load it and write to SU headers later
  if (myrank == 0) then
    open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_sources.txt',status='unknown')
    do isource=1,NSOURCES
      write(IOUT_SU,*) x_found(isource),y_found(isource),z_found(isource)
    enddo
    close(IOUT_SU)
  endif

  is_done_sources = .false.

  ! user output
  if (myrank == 0) then

    ! skip source info for wavefield injection simulations
    if (COUPLE_WITH_INJECTION_TECHNIQUE) then
      write(IMAIN,*)
      write(IMAIN,*) 'note on using injection technique:'
      write(IMAIN,*) '  source (inside the mesh) will be ignored if we are coupling with DSM/AxiSEM/FK,'
      write(IMAIN,*) '  because the source is precisely the wavefield coming from the injection boundary'
      write(IMAIN,*)
      ! nothing to display anymore
      is_done_sources = .true.
    endif
  endif

  ! user output with source details
  if (myrank == 0 .and. (.not. is_done_sources)) then
    ! loops over all sources
    do isource = 1,NSOURCES

      if (SHOW_DETAILS_LOCATE_SOURCE .or. NSOURCES < 10) then
        ! source info
        write(IMAIN,*)
        write(IMAIN,*) 'source # ',isource

        ! source type
        write(IMAIN,*) '  source located in slice ',islice_selected_source(isource)
        write(IMAIN,*) '                 in element ',ispec_selected_source(isource)
        if (idomain(isource) == IDOMAIN_ACOUSTIC) then
          write(IMAIN,*) '                 in acoustic domain'
        else if (idomain(isource) == IDOMAIN_ELASTIC) then
          write(IMAIN,*) '                 in elastic domain'
        else if (idomain(isource) == IDOMAIN_POROELASTIC) then
          write(IMAIN,*) '                 in poroelastic domain'
        else
          write(IMAIN,*) '                 in unknown domain'
        endif
        write(IMAIN,*)

        ! source location (reference element)
        if (USE_FORCE_POINT_SOURCE) then
          ! single point force
          write(IMAIN,*) '  using force point source: '
          write(IMAIN,*) '    xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '    eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) '    gamma coordinate of source in that element: ',gamma_source(isource)

          write(IMAIN,*)
          write(IMAIN,*) '    component of direction vector in East direction: ',sngl(comp_dir_vect_source_E(isource))
          write(IMAIN,*) '    component of direction vector in North direction: ',sngl(comp_dir_vect_source_N(isource))
          write(IMAIN,*) '    component of direction vector in Vertical direction: ',sngl(comp_dir_vect_source_Z_UP(isource))
          write(IMAIN,*)
          write(IMAIN,*) '    nu1 = ',nu_source(1,:,isource)
          write(IMAIN,*) '    nu2 = ',nu_source(2,:,isource)
          write(IMAIN,*) '    nu3 = ',nu_source(3,:,isource)
          write(IMAIN,*)
          write(IMAIN,*) '    at (x,y,z) coordinates = ',x_found(isource),y_found(isource),z_found(isource)
        else
          ! moment tensor
          write(IMAIN,*) '  using moment tensor source: '
          write(IMAIN,*) '    xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '    eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) '    gamma coordinate of source in that element: ',gamma_source(isource)
        endif
        write(IMAIN,*)

        ! source time function info
        write(IMAIN,*) '  source time function:'
        if (USE_EXTERNAL_SOURCE_FILE) then
          ! external STF
          write(IMAIN,*) '    using external source time function'
          write(IMAIN,*)
        else
          ! STF details
          if (USE_FORCE_POINT_SOURCE) then
            ! force sources
            if (USE_RICKER_TIME_FUNCTION) then
              write(IMAIN,*) '    using Ricker source time function'
            else
              ! acoustic/elastic/.. domain by default uses a Gaussian STF
              write(IMAIN,*) '    using Gaussian source time function'
            endif
          else
            ! CMT sources
            if (USE_RICKER_TIME_FUNCTION) then
              write(IMAIN,*) '    using Ricker source time function'
            else
              if (idomain(isource) == IDOMAIN_ACOUSTIC .or. idomain(isource) == IDOMAIN_POROELASTIC) then
                ! acoustic/poroelastic domain by default uses a Gaussian STF
                write(IMAIN,*) '    using Gaussian source time function'
              else
                ! elastic domain by default uses a Heaviside STF
                write(IMAIN,*) '    using Heaviside source time function'
              endif
            endif
          endif
          ! acoustic pressure trick
          if (idomain(isource) == IDOMAIN_ACOUSTIC) then
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              write(IMAIN,*) '    using trick for better pressure (second derivatives)'
            endif
          endif

          ! frequency/half-duration
          if (USE_FORCE_POINT_SOURCE) then
            ! single point force
            ! prints frequency content for point forces
            f0 = hdur(isource)
            if (USE_RICKER_TIME_FUNCTION) then
              write(IMAIN,*) '    using a source of dominant frequency ',f0

              t0_ricker = 1.2d0/f0
              write(IMAIN,*) '    t0_ricker = ',t0_ricker
              write(IMAIN,*) '    Ricker frequency: ',hdur(isource),' Hz'
            else
              if (idomain(isource) == IDOMAIN_ACOUSTIC) then
                write(IMAIN,*) '    Gaussian half duration: ',5.d0*DT,' seconds'
              else if (idomain(isource) == IDOMAIN_ELASTIC) then
                write(IMAIN,*) '    Gaussian half duration: ',hdur(isource)/SOURCE_DECAY_MIMIC_TRIANGLE,' seconds'
              else if (idomain(isource) == IDOMAIN_POROELASTIC) then
                write(IMAIN,*) '    Gaussian half duration: ',5.d0*DT,' seconds'
              endif
            endif
            write(IMAIN,*)
          else
            ! moment-tensor
            if (USE_RICKER_TIME_FUNCTION) then
              write(IMAIN,*) '    Ricker frequency: ',hdur(isource),' Hz'
            else
              ! add message if source is a Heaviside
              if (hdur(isource) <= 5.*DT) then
                write(IMAIN,*)
                write(IMAIN,*) '    Source time function is a Heaviside, convolve later'
                write(IMAIN,*)
              endif
              write(IMAIN,*) '    half duration: ',hdur(isource),' seconds'
            endif
          endif
        endif
        write(IMAIN,*) '    time shift: ',tshift_src(isource),' seconds'
        write(IMAIN,*)

        ! magnitude
#ifdef DEBUG_COUPLED
        write(IMAIN,*)
        write(IMAIN,*) 'Coupled activated, thus not including any internal source'
        write(IMAIN,*)
#else
        write(IMAIN,*) '  magnitude of the source:'
        if (USE_FORCE_POINT_SOURCE) then
          ! single point force
          write(IMAIN,*) '    factor = ', factor_force_source(isource)
        else
          ! moment-tensor
          write(IMAIN,*) '       scalar moment M0 = ', &
            get_cmt_scalar_moment(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource)),' dyne-cm'
          write(IMAIN,*) '    moment magnitude Mw = ', &
            get_cmt_moment_magnitude(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource))
        endif
#endif
        write(IMAIN,*)

        ! location accuracy
        write(IMAIN,*) '  original (requested) position of the source:'
        write(IMAIN,*)
        write(IMAIN,*) '            latitude: ',lat(isource)
        write(IMAIN,*) '           longitude: ',long(isource)
        write(IMAIN,*)
        if (SUPPRESS_UTM_PROJECTION) then
          write(IMAIN,*) '               x: ',utm_x_source(isource)
          write(IMAIN,*) '               y: ',utm_y_source(isource)
        else
          write(IMAIN,*) '           UTM x: ',utm_x_source(isource)
          write(IMAIN,*) '           UTM y: ',utm_y_source(isource)
        endif
        if (USE_SOURCES_RECEIVERS_Z) then
          write(IMAIN,*) '           z: ',depth(isource),' m'
        else
          write(IMAIN,*) '           depth: ',depth(isource)/1000.0,' km'
          write(IMAIN,*) '  topo elevation: ',elevation(isource)
        endif

        write(IMAIN,*)
        write(IMAIN,*) '  position of the source that will be used:'
        write(IMAIN,*)
        if (SUPPRESS_UTM_PROJECTION) then
          write(IMAIN,*) '               x: ',x_found(isource)
          write(IMAIN,*) '               y: ',y_found(isource)
        else
          write(IMAIN,*) '           UTM x: ',x_found(isource)
          write(IMAIN,*) '           UTM y: ',y_found(isource)
        endif
        if (USE_SOURCES_RECEIVERS_Z) then
          write(IMAIN,*) '               z: ',z_found(isource)
        else
          write(IMAIN,*) '           depth: ',dabs(z_found(isource) - elevation(isource))/1000.,' km'
          write(IMAIN,*) '               z: ',z_found(isource)
        endif
        write(IMAIN,*)

        ! display error in location estimate
        write(IMAIN,*) '  error in location of the source: ',sngl(final_distance(isource)),' m'
        write(IMAIN,*)

        ! add warning if estimate is poor
        ! (usually means source outside the mesh given by the user)
        if (final_distance(isource) > elemsize_max_glob) then
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
        endif
        ! add warning if located in PML
        if (is_CPML_source_all(isource)) then
          write(IMAIN,*) '*******************************************************'
          write(IMAIN,*) '***** WARNING: source located in C-PML region *********'
          write(IMAIN,*) '*******************************************************'
          write(IMAIN,*)
        endif
        write(IMAIN,*)
      endif  ! end of detailed output to locate source

      ! checks CMTSOLUTION format for acoustic case
      if (idomain(isource) == IDOMAIN_ACOUSTIC .and. .not. USE_FORCE_POINT_SOURCE) then
        if (Mxx(isource) /= Myy(isource) .or. Myy(isource) /= Mzz(isource) .or. &
           Mxy(isource) > TINYVAL .or. Mxz(isource) > TINYVAL .or. Myz(isource) > TINYVAL) then
          write(IMAIN,*)
          write(IMAIN,*) 'Error CMTSOLUTION format for acoustic source:'
          write(IMAIN,*) '   acoustic source needs explosive moment tensor with'
          write(IMAIN,*) '      Mrr = Mtt = Mpp '
          write(IMAIN,*) '   and '
          write(IMAIN,*) '      Mrt = Mrp = Mtp = zero'
          write(IMAIN,*)
          ! stop
          print *,'Error: invalid acoustic source, must have: Mrr == Mtt == Mpp and Mrt == Mrp == Mtp == zero'
          print *,'       current Mrr = ',Mzz(isource),' Mtt = ',Myy(isource),' Mpp = ',Mxx(isource)
          print *,'               Mrt = ',-Myz(isource),' Mrp = ',Mxz(isource),' Mtp = ',-Mxy(isource)
          call exit_mpi(myrank,'Invalid acoustic source')
        endif
      endif

      ! checks source domain
      if (idomain(isource) /= IDOMAIN_ACOUSTIC .and. &
          idomain(isource) /= IDOMAIN_ELASTIC .and. &
          idomain(isource) /= IDOMAIN_POROELASTIC) then
        call exit_MPI(myrank,'source located in unknown domain')
      endif

    ! end of loop on all the sources
    enddo

    if (.not. SHOW_DETAILS_LOCATE_SOURCE .and. NSOURCES > 1) then
      write(IMAIN,*)
      write(IMAIN,*) '*************************************'
      write(IMAIN,*) ' using sources ',NSOURCES
      write(IMAIN,*) '*************************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! display maximum error in location estimate
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(maxval(final_distance)),' m'
    write(IMAIN,*)
    call flush_IMAIN()

  endif     ! end of section executed by main process only

  ! elapsed time since beginning of source detection
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for detection of sources in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of source detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! frees temporary arrays
  deallocate(moment_tensor,lat,long,depth)
  deallocate(x_found,y_found,z_found,elevation,final_distance)
  deallocate(x_target,y_target,z_target,idomain)
  deallocate(is_CPML_source,is_CPML_source_all)

  end subroutine locate_source

