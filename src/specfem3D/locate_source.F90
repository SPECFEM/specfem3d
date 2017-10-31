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
      NSTEP_STF,NSOURCES_STF,USE_EXTERNAL_SOURCE_FILE,USE_TRICK_FOR_BETTER_PRESSURE, &
      ibool,myrank,NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,DT,NSOURCES


  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: filename
  double precision,dimension(NSOURCES),intent(inout) :: tshift_src
  double precision,intent(inout) :: min_tshift_src_original
  double precision, dimension(NDIM,NDIM,NSOURCES),intent(out) :: nu_source
  double precision, dimension(NSOURCES) :: factor_force_source
  double precision, dimension(NSOURCES) :: comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF,NSOURCES_STF) :: user_source_time_function

  integer isource

  double precision, dimension(NSOURCES),intent(inout) :: utm_x_source,utm_y_source
  double precision final_distance_source(NSOURCES)
  double precision x_target_source,y_target_source,z_target_source

  integer,intent(inout) :: islice_selected_source(NSOURCES)

  ! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! sources
  integer,intent(inout) :: ispec_selected_source(NSOURCES)
  double precision :: f0,t0_ricker

  ! CMTs
  double precision, dimension(NSOURCES),intent(inout) :: hdur
  double precision, dimension(NSOURCES),intent(inout) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(NSOURCES) :: lat,long,depth
  double precision, dimension(6,NSOURCES) ::  moment_tensor

  ! positioning
  double precision, dimension(NSOURCES),intent(inout) :: xi_source,eta_source,gamma_source
  double precision, dimension(NSOURCES) :: x_found_source,y_found_source,z_found_source
  double precision, dimension(NSOURCES) :: elevation

  integer, dimension(NSOURCES) :: idomain

  double precision, external :: get_cmt_scalar_moment
  double precision, external :: get_cmt_moment_magnitude

  ! location search
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  !-----------------------------------------------------------------------------------

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
      ! only master process reads in FORCESOLUTION file
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
      ! only master process reads in CMTSOLUTION file
      call get_cmt(filename,tshift_src,hdur,lat,long,depth,moment_tensor, &
                   DT,NSOURCES,min_tshift_src_original,user_source_time_function)
    endif
    ! broadcasts specific moment tensor infos
    call bcast_all_dp(moment_tensor,6*NSOURCES)
  endif

  ! broadcasts general source information read on the master to the nodes
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
  time_start = wtime()

  ! user output
  if (myrank == 0) then
    if (SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,*) 'no UTM projection'
    else
      write(IMAIN,*) 'UTM projection:'
      write(IMAIN,*) '  UTM zone: ',UTM_PROJECTION_ZONE
    endif
    if (USE_SOURCES_RECEIVERS_Z) then
      write(IMAIN,*) '  (depth) becomes directly (z) coordinate'
    endif
  endif

  !
  ! r -> z, theta -> -y, phi -> x
  !
  !  Mrr =  Mzz
  !  Mtt =  Myy
  !  Mpp =  Mxx
  !  Mrt = -Myz
  !  Mrp =  Mxz
  !  Mtp = -Mxy

  ! get the moment tensor
  Mzz(:) = + moment_tensor(1,:)
  Mxx(:) = + moment_tensor(3,:)
  Myy(:) = + moment_tensor(2,:)
  Mxz(:) = + moment_tensor(5,:)
  Myz(:) = - moment_tensor(4,:)
  Mxy(:) = - moment_tensor(6,:)

  ! loop on all the sources
  do isource = 1,NSOURCES

    ! get z target coordinate, depending on the topography
    if (.not. USE_SOURCES_RECEIVERS_Z) depth(isource) = depth(isource)*1000.0d0
    call get_elevation_and_z_coordinate(long(isource),lat(isource),utm_x_source(isource),utm_y_source(isource),z_target_source, &
                                        elevation(isource),depth(isource))
    x_target_source = utm_x_source(isource)
    y_target_source = utm_y_source(isource)

    call locate_point_in_mesh(x_target_source, y_target_source, z_target_source, SOURCES_CAN_BE_BURIED, elemsize_max_glob, &
            ispec_selected_source(isource), xi_source(isource), eta_source(isource), gamma_source(isource), &
            x_found_source(isource), y_found_source(isource), z_found_source(isource), idomain(isource),nu_source(:,:,isource))

    ! synchronize all the processes to make sure all the estimates are available
    call synchronize_all()

    call locate_MPI_slice_and_bcast_to_all(x_target_source, y_target_source, z_target_source, &
                                           x_found_source(isource), y_found_source(isource), z_found_source(isource), &
                                           xi_source(isource), eta_source(isource), gamma_source(isource), &
                                           ispec_selected_source(isource), islice_selected_source(isource), &
                                           final_distance_source(isource), idomain(isource),nu_source(:,:,isource))

  enddo ! end of loop on all the sources

  if (myrank == 0) then

    do isource = 1,NSOURCES

      if (SHOW_DETAILS_LOCATE_SOURCE .or. NSOURCES == 1) then
        ! source info
        write(IMAIN,*)
        write(IMAIN,*) '*************************************'
        write(IMAIN,*) ' locating source ',isource
        write(IMAIN,*) '*************************************'
        write(IMAIN,*)

        ! source type
        write(IMAIN,*) 'source located in slice ',islice_selected_source(isource)
        write(IMAIN,*) '               in element ',ispec_selected_source(isource)
        if (idomain(isource) == IDOMAIN_ACOUSTIC) then
          write(IMAIN,*) '               in acoustic domain'
        else if (idomain(isource) == IDOMAIN_ELASTIC) then
          write(IMAIN,*) '               in elastic domain'
        else if (idomain(isource) == IDOMAIN_POROELASTIC) then
          write(IMAIN,*) '               in poroelastic domain'
        else
          write(IMAIN,*) '               in unknown domain'
        endif
        write(IMAIN,*)

        ! source location (reference element)
        if (USE_FORCE_POINT_SOURCE) then
          ! single point force
          write(IMAIN,*) 'using force point source: '
          write(IMAIN,*) '  xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '  eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) '  gamma coordinate of source in that element: ',gamma_source(isource)

          write(IMAIN,*)
          write(IMAIN,*) '  component of direction vector in East direction: ',comp_dir_vect_source_E(isource)
          write(IMAIN,*) '  component of direction vector in North direction: ',comp_dir_vect_source_N(isource)
          write(IMAIN,*) '  component of direction vector in Vertical direction: ',comp_dir_vect_source_Z_UP(isource)
          write(IMAIN,*)
          write(IMAIN,*) '  nu1 = ',nu_source(1,:,isource)
          write(IMAIN,*) '  nu2 = ',nu_source(2,:,isource)
          write(IMAIN,*) '  nu3 = ',nu_source(3,:,isource)
          write(IMAIN,*)
          write(IMAIN,*) '  at (x,y,z) coordinates = ',x_found_source(isource),y_found_source(isource),z_found_source(isource)
        else
          ! moment tensor
          write(IMAIN,*) 'using moment tensor source: '
          write(IMAIN,*) '  xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '  eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) '  gamma coordinate of source in that element: ',gamma_source(isource)
        endif
        write(IMAIN,*)

        ! source time function info
        write(IMAIN,*) 'source time function:'
        if (USE_EXTERNAL_SOURCE_FILE) then
          ! external STF
          write(IMAIN,*) '  using external source time function'
          write(IMAIN,*)
        else
          ! STF details
          if (USE_RICKER_TIME_FUNCTION) then
            write(IMAIN,*) '  using Ricker source time function'
          else
            write(IMAIN,*) '  using Gaussian source time function'
          endif
          if (idomain(isource) == IDOMAIN_ACOUSTIC) then
            if (USE_TRICK_FOR_BETTER_PRESSURE) then
              write(IMAIN,*) '  using trick for better pressure (second derivatives)'
            endif
          endif

          ! frequency/half-duration
          if (USE_FORCE_POINT_SOURCE) then
            ! single point force
            ! prints frequency content for point forces
            f0 = hdur(isource)
            if (USE_RICKER_TIME_FUNCTION) then
              write(IMAIN,*) '  using a source of dominant frequency ',f0

              t0_ricker = 1.2d0/f0
              write(IMAIN,*) '  t0_ricker = ',t0_ricker
              write(IMAIN,*) '  Ricker frequency: ',hdur(isource),' Hz'
            else
              if (idomain(isource) == IDOMAIN_ACOUSTIC) then
                write(IMAIN,*) '  Gaussian half duration: ',5.d0*DT,' seconds'
              else if (idomain(isource) == IDOMAIN_ELASTIC) then
                write(IMAIN,*) '  Gaussian half duration: ',hdur(isource)/SOURCE_DECAY_MIMIC_TRIANGLE,' seconds'
              else if (idomain(isource) == IDOMAIN_POROELASTIC) then
                write(IMAIN,*) '  Gaussian half duration: ',5.d0*DT,' seconds'
              endif
            endif
            write(IMAIN,*)
          else
            ! moment-tensor
            if (USE_RICKER_TIME_FUNCTION) then
              write(IMAIN,*) '  Ricker frequency: ',hdur(isource),' Hz'
            else
              ! add message if source is a Heaviside
              if (hdur(isource) <= 5.*DT) then
                write(IMAIN,*)
                write(IMAIN,*) '  Source time function is a Heaviside, convolve later'
                write(IMAIN,*)
              endif
              write(IMAIN,*) '  half duration: ',hdur(isource),' seconds'
            endif
          endif
        endif
        write(IMAIN,*) '  time shift: ',tshift_src(isource),' seconds'
        write(IMAIN,*)

        ! magnitude
#ifdef DEBUG_COUPLED
        write(IMAIN,*)
        write(IMAIN,*) 'Coupled activated, thus not including any internal source'
        write(IMAIN,*)
#else
        write(IMAIN,*) 'magnitude of the source:'
        if (USE_FORCE_POINT_SOURCE) then
          ! single point force
          write(IMAIN,*) '  factor = ', factor_force_source(isource)
        else
          ! moment-tensor
          write(IMAIN,*) '     scalar moment M0 = ', &
            get_cmt_scalar_moment(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource)),' dyne-cm'
          write(IMAIN,*) '  moment magnitude Mw = ', &
            get_cmt_moment_magnitude(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource))
        endif
#endif
        write(IMAIN,*)

        ! location accuracy
        write(IMAIN,*) 'original (requested) position of the source:'
        write(IMAIN,*)
        write(IMAIN,*) '          latitude: ',lat(isource)
        write(IMAIN,*) '         longitude: ',long(isource)
        write(IMAIN,*)
        if (SUPPRESS_UTM_PROJECTION) then
          write(IMAIN,*) '             x: ',utm_x_source(isource)
          write(IMAIN,*) '             y: ',utm_y_source(isource)
        else
          write(IMAIN,*) '         UTM x: ',utm_x_source(isource)
          write(IMAIN,*) '         UTM y: ',utm_y_source(isource)
        endif
        if (USE_SOURCES_RECEIVERS_Z) then
          write(IMAIN,*) '         z: ',depth(isource),' m'
        else
          write(IMAIN,*) '         depth: ',depth(isource)/1000.0,' km'
          write(IMAIN,*) 'topo elevation: ',elevation(isource)
        endif

        write(IMAIN,*)
        write(IMAIN,*) 'position of the source that will be used:'
        write(IMAIN,*)
        if (SUPPRESS_UTM_PROJECTION) then
          write(IMAIN,*) '             x: ',x_found_source(isource)
          write(IMAIN,*) '             y: ',y_found_source(isource)
        else
          write(IMAIN,*) '         UTM x: ',x_found_source(isource)
          write(IMAIN,*) '         UTM y: ',y_found_source(isource)
        endif
        if (USE_SOURCES_RECEIVERS_Z) then
          write(IMAIN,*) '             z: ',z_found_source(isource)
        else
          write(IMAIN,*) '         depth: ',dabs(z_found_source(isource) - elevation(isource))/1000.,' km'
          write(IMAIN,*) '             z: ',z_found_source(isource)
        endif
        write(IMAIN,*)

        ! display error in location estimate
        write(IMAIN,*) 'error in location of the source: ',sngl(final_distance_source(isource)),' m'

        ! add warning if estimate is poor
        ! (usually means source outside the mesh given by the user)
        if (final_distance_source(isource) > elemsize_max_glob) then
          write(IMAIN,*)
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
        endif

      endif  ! end of detailed output to locate source

      ! checks CMTSOLUTION format for acoustic case
      if (idomain(isource) == IDOMAIN_ACOUSTIC .and. .not. USE_FORCE_POINT_SOURCE) then
        if (Mxx(isource) /= Myy(isource) .or. Myy(isource) /= Mzz(isource) .or. &
           Mxy(isource) > TINYVAL .or. Mxz(isource) > TINYVAL .or. Myz(isource) > TINYVAL) then
          write(IMAIN,*)
          write(IMAIN,*) ' error CMTSOLUTION format for acoustic source:'
          write(IMAIN,*) '   acoustic source needs explosive moment tensor with'
          write(IMAIN,*) '      Mrr = Mtt = Mpp '
          write(IMAIN,*) '   and '
          write(IMAIN,*) '      Mrt = Mrp = Mtp = zero'
          write(IMAIN,*)
          call exit_mpi(myrank,'error acoustic source')
        endif
      endif

      ! checks source domain
      if (idomain(isource) /= IDOMAIN_ACOUSTIC .and. idomain(isource) /= IDOMAIN_ELASTIC .and. &
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
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(maxval(final_distance_source)),' m'
    write(IMAIN,*)
    call flush_IMAIN()

  endif     ! end of section executed by main process only

  ! sets new utm coordinates for best locations
  utm_x_source(:) = x_found_source(:)
  utm_y_source(:) = y_found_source(:)

  ! elapsed time since beginning of source detection
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for detection of sources in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of source detection - done'
    write(IMAIN,*)
    call flush_IMAIN()

    ! output source information to a file so that we can load it and write to SU headers later
    open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_sources.txt',status='unknown')
    do isource=1,NSOURCES
      write(IOUT_SU,*) x_found_source(isource),y_found_source(isource),z_found_source(isource)
    enddo
    close(IOUT_SU)
  endif

  end subroutine locate_source

