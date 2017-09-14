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

  subroutine locate_source(NSOURCES,tshift_src,min_tshift_src_original,yr,jda,ho,mi,utm_x_source,utm_y_source, &
                           hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                           islice_selected_source,ispec_selected_source, &
                           xi_source,eta_source,gamma_source)

  use constants

  use specfem_par, only: USE_FORCE_POINT_SOURCE,USE_RICKER_TIME_FUNCTION, &
      UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,USE_SOURCES_RECEIVERS_Z, &
      factor_force_source,comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP, &
      user_source_time_function,NSTEP_STF,NSOURCES_STF,USE_EXTERNAL_SOURCE_FILE,USE_TRICK_FOR_BETTER_PRESSURE,&
      ibool,myrank,NSPEC_AB,NGLOB_AB,NGNOD,xstore,ystore,zstore,xigll,yigll,zigll,NPROC, &
      DT,iglob_is_surface_external_mesh,ispec_is_surface_external_mesh,num_free_surface_faces,free_surface_ispec,free_surface_ijk

  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_poroelastic, only: ispec_is_poroelastic


  implicit none

  integer,intent(in) :: NSOURCES
  integer,intent(inout) :: yr,jda,ho,mi
  double precision,dimension(NSOURCES),intent(inout) :: tshift_src
  double precision,intent(inout) :: min_tshift_src_original
  double precision :: sec

  integer iprocloop

  integer i,j,k,ispec,iglob,isource
  integer iproc(1)

  double precision, dimension(NSOURCES),intent(inout) :: utm_x_source,utm_y_source
  double precision dist_squared
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta

  ! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

  ! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz
  double precision dgamma

  double precision final_distance_source(NSOURCES)

  double precision x_target_source,y_target_source,z_target_source

  double precision,dimension(1) :: altitude_source,distmin_ele

  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all

  real(kind=CUSTOM_REAL) :: xloc,yloc,loc_ele,loc_distmin

  integer,intent(inout) :: islice_selected_source(NSOURCES)

  ! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! sources
  integer,intent(inout) :: ispec_selected_source(NSOURCES)
  integer :: ngather, ns, ne, ig, is, ng

  integer, dimension(NGATHER_SOURCES,0:NPROC-1) :: ispec_selected_source_all
  double precision, dimension(NGATHER_SOURCES,0:NPROC-1) :: xi_source_all,eta_source_all,gamma_source_all
  double precision, dimension(NGATHER_SOURCES,0:NPROC-1) :: x_found_source_all,y_found_source_all,z_found_source_all
  double precision, dimension(NGATHER_SOURCES,0:NPROC-1) :: final_distance_source_all

  double precision, dimension(:), allocatable :: tmp_local
  double precision, dimension(:,:),allocatable :: tmp_all_local

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
  double precision :: distmin_squared,distmin_not_squared

  integer, dimension(:), allocatable :: tmp_i_local
  integer, dimension(:,:),allocatable :: tmp_i_all_local

  integer ix_initial_guess_source,iy_initial_guess_source,iz_initial_guess_source
  integer ier

  integer, dimension(NSOURCES) :: idomain
  integer, dimension(NGATHER_SOURCES,0:NPROC-1) :: idomain_all

  double precision, external :: get_cmt_scalar_moment
  double precision, external :: get_cmt_moment_magnitude

  ! location search
  logical :: located_target
  double precision :: typical_size_squared
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
      call get_force(tshift_src,hdur,lat,long,depth,NSOURCES,min_tshift_src_original,factor_force_source, &
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
      call get_cmt(yr,jda,ho,mi,sec,tshift_src,hdur,lat,long,depth,moment_tensor, &
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

  ! define topology of the control element
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

  ! compute typical size of elements
  if (USE_DISTANCE_CRITERION_SOURCES) then
    ! gets mesh dimensions
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)

    ! sets typical element size for search
    typical_size_squared =  elemsize_max_glob

    ! use 10 times the distance as a criterion for source detection
    typical_size_squared = (10. * typical_size_squared)**2
  endif

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

  ! loop on all the sources
  do isource = 1,NSOURCES

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
    Mzz(isource) = + moment_tensor(1,isource)
    Mxx(isource) = + moment_tensor(3,isource)
    Myy(isource) = + moment_tensor(2,isource)
    Mxz(isource) = + moment_tensor(5,isource)
    Myz(isource) = - moment_tensor(4,isource)
    Mxy(isource) = - moment_tensor(6,isource)

    ! gets UTM x,y
    call utm_geo(long(isource),lat(isource),utm_x_source(isource),utm_y_source(isource), &
                   UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

    ! get approximate topography elevation at source long/lat coordinates
    xloc = utm_x_source(isource)
    yloc = utm_y_source(isource)
    call get_topo_elevation_free(xloc,yloc,loc_ele,loc_distmin, &
                              NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              num_free_surface_faces,free_surface_ispec,free_surface_ijk)
    altitude_source(1) = loc_ele
    distmin_ele(1) = loc_distmin

    !  MPI communications to determine the best slice
    call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
    call gather_all_dp(altitude_source,1,elevation_all,1,NPROC)

    if (myrank == 0) then
      iproc = minloc(distmin_ele_all)
      altitude_source(1) = elevation_all(iproc(1))
    endif
    call bcast_all_dp(altitude_source,1)
    elevation(isource) = altitude_source(1)

    x_target_source = utm_x_source(isource)
    y_target_source = utm_y_source(isource)

    ! source Z coordinate
    if (USE_SOURCES_RECEIVERS_Z) then
      ! alternative: depth is given as z value directly
      z_target_source = depth(isource)
    else
      ! depth in CMTSOLUTION given in km
      z_target_source =  - depth(isource)*1000.0d0 + elevation(isource)
    endif

    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    ! flag to check that we located at least one target element
    located_target = .false.

    ispec_selected_source(isource) = 0
    ix_initial_guess_source = 1
    iy_initial_guess_source = 1
    iz_initial_guess_source = 1

    do ispec = 1,NSPEC_AB

      ! exclude elements that are too far from target
      if (USE_DISTANCE_CRITERION_SOURCES) then
        iglob = ibool(MIDX,MIDY,MIDZ,ispec)
        dist_squared = (x_target_source - dble(xstore(iglob)))**2 &
                     + (y_target_source - dble(ystore(iglob)))**2 &
                     + (z_target_source - dble(zstore(iglob)))**2
        if (dist_squared > typical_size_squared) cycle
      endif

      do k = 2,NGLLZ - 1
        do j = 2,NGLLY - 1
          do i = 2,NGLLX - 1

            iglob = ibool(i,j,k,ispec)

            if (.not. SOURCES_CAN_BE_BURIED) then
              if ((.not. iglob_is_surface_external_mesh(iglob)) .or. (.not. ispec_is_surface_external_mesh(ispec))) then
                cycle
              endif
            endif

            ! keep this point if it is closer to the source
            !  we compare squared distances instead of distances themselves to significantly speed up calculations
            dist_squared = (x_target_source - dble(xstore(iglob)))**2 &
                         + (y_target_source - dble(ystore(iglob)))**2 &
                         + (z_target_source - dble(zstore(iglob)))**2
            if (dist_squared < distmin_squared) then
              distmin_squared = dist_squared
              ispec_selected_source(isource) = ispec
              ix_initial_guess_source = i
              iy_initial_guess_source = j
              iz_initial_guess_source = k
              located_target = .true.

              x_found_source(isource) = xstore(iglob)
              y_found_source(isource) = ystore(iglob)
              z_found_source(isource) = zstore(iglob)

            endif

          enddo
        enddo
      enddo

    ! end of loop on all the elements in current slice
    enddo

    ! if we have not located a target element, the source is not in this slice
    ! therefore use first element only for fictitious iterative search
    if (.not. located_target) then
      ispec_selected_source(isource) = 1
      ix_initial_guess_source = 1
      iy_initial_guess_source = 1
      iz_initial_guess_source = 1
      final_distance_source(isource) = HUGEVAL
    endif

    ! sets whether acoustic (1) or elastic (2)
    if (ispec_is_acoustic( ispec_selected_source(isource) )) then
      idomain(isource) = IDOMAIN_ACOUSTIC
    else if (ispec_is_elastic( ispec_selected_source(isource) )) then
      idomain(isource) = IDOMAIN_ELASTIC
    else if (ispec_is_poroelastic( ispec_selected_source(isource) )) then
      idomain(isource) = IDOMAIN_POROELASTIC
    else
      idomain(isource) = 0
    endif

    ! store xi,eta,gamma and x,y,z of point found
    ! note: here they have range [1.0d0,NGLLX/Y/Z]
    xi_source(isource) = dble(ix_initial_guess_source)
    eta_source(isource) = dble(iy_initial_guess_source)
    gamma_source(isource) = dble(iz_initial_guess_source)


! *******************************************
! find the best (xi,eta,gamma) for the source
! *******************************************

    ! this tries to find best location
    if (USE_BEST_LOCATION_FOR_SOURCE) then

      ! uses actual location interpolators, in range [-1,1]
      xi = xigll(ix_initial_guess_source)
      eta = yigll(iy_initial_guess_source)
      gamma = zigll(iz_initial_guess_source)

      ! define coordinates of the control points of the element
      do ia=1,NGNOD
        iax = 0
        iay = 0
        iaz = 0
        if (iaddx(ia) == 0) then
          iax = 1
        else if (iaddx(ia) == 1) then
          iax = (NGLLX+1)/2
        else if (iaddx(ia) == 2) then
          iax = NGLLX
        else
          call exit_MPI(myrank,'incorrect value of iaddx')
        endif

        if (iaddy(ia) == 0) then
          iay = 1
        else if (iaddy(ia) == 1) then
          iay = (NGLLY+1)/2
        else if (iaddy(ia) == 2) then
          iay = NGLLY
        else
          call exit_MPI(myrank,'incorrect value of iaddy')
        endif

        if (iaddz(ia) == 0) then
          iaz = 1
        else if (iaddz(ia) == 1) then
          iaz = (NGLLZ+1)/2
        else if (iaddz(ia) == 2) then
          iaz = NGLLZ
        else
          call exit_MPI(myrank,'incorrect value of iaddz')
        endif

        iglob = ibool(iax,iay,iaz,ispec_selected_source(isource))
        xelm(ia) = dble(xstore(iglob))
        yelm(ia) = dble(ystore(iglob))
        zelm(ia) = dble(zstore(iglob))

      enddo

      ! iterate to solve the non linear system
      do iter_loop = 1,NUM_ITER

        ! recompute jacobian for the new point
        call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

        ! compute distance to target location
        dx = - (x - x_target_source)
        dy = - (y - y_target_source)
        dz = - (z - z_target_source)

        ! compute increments
        dxi  = xix*dx + xiy*dy + xiz*dz
        deta = etax*dx + etay*dy + etaz*dz
        dgamma = gammax*dx + gammay*dy + gammaz*dz

        ! update values
        xi = xi + dxi
        eta = eta + deta
        gamma = gamma + dgamma

        ! impose that we stay in that element
        ! (useful if user gives a receiver outside the mesh for instance)
        ! we can go slightly outside the [1,1] segment since with finite elements
        ! the polynomial solution is defined everywhere
        ! this can be useful for convergence of itertive scheme with distorted elements.
        ! this is purposely set to 1.10, do *NOT* set it back to 1.00 because 1.10 gives more accurate locations
        if (xi > 1.10d0) xi = 1.10d0
        if (xi < -1.10d0) xi = -1.10d0
        if (eta > 1.10d0) eta = 1.10d0
        if (eta < -1.10d0) eta = -1.10d0
        if (gamma > 1.10d0) gamma = 1.10d0
        if (gamma < -1.10d0) gamma = -1.10d0

      enddo

      ! compute final coordinates of point found
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

      ! store xi,eta,gamma and x,y,z of point found
      ! note: xi/eta/gamma will be in range [-1,1]
      xi_source(isource) = xi
      eta_source(isource) = eta
      gamma_source(isource) = gamma

      x_found_source(isource) = x
      y_found_source(isource) = y
      z_found_source(isource) = z

    else

      ! takes initial GLL point guess and uses actual location interpolators, in range [-1,1]
      ! note: xi/eta/gamma will be in range [-1,1]
      xi_source(isource) = xigll(ix_initial_guess_source)
      eta_source(isource) = yigll(iy_initial_guess_source)
      gamma_source(isource) = zigll(iz_initial_guess_source)

    endif ! USE_BEST_LOCATION_FOR_SOURCE

    ! compute final distance between asked and found (converted to km)
    final_distance_source(isource) = dsqrt((x_target_source-x_found_source(isource))**2 + &
                (y_target_source-y_found_source(isource))**2 + (z_target_source-z_found_source(isource))**2)

  ! end of loop on all the sources
  enddo

  ! now gather information from all the nodes
  !
  ! note: gathering source information is done here for bins of NGATHER_SOURCES
  !       to save memory for very large number of sources

  ! number of gather bins
  ngather = NSOURCES/NGATHER_SOURCES
  if (mod(NSOURCES,NGATHER_SOURCES) /= 0) ngather = ngather+1

  ! loops over single bin
  do ig = 1, ngather
    ! source index start
    ns = (ig-1) * NGATHER_SOURCES + 1
    ! source index end
    ne = min(ig*NGATHER_SOURCES, NSOURCES)
    ! number of gathering sources
    ng = ne - ns + 1

    ispec_selected_source_all(:,:) = -1

    ! avoids warnings about temporary creations of arrays for function call by compiler
    allocate(tmp_i_local(ng),tmp_i_all_local(ng,0:NPROC-1),stat=ier)
    if (ier /= 0) stop 'error allocating array tmp_i_local'
    tmp_i_local(:) = ispec_selected_source(ns:ne)
    call gather_all_i(tmp_i_local,ng,tmp_i_all_local,ng,NPROC)
    ispec_selected_source_all(1:ng,:) = tmp_i_all_local(:,:)

    ! acoustic/elastic domain
    tmp_i_local(:) = idomain(ns:ne)
    call gather_all_i(tmp_i_local,ng,tmp_i_all_local,ng,NPROC)
    idomain_all(1:ng,:) = tmp_i_all_local(:,:)

    deallocate(tmp_i_local,tmp_i_all_local)

    ! avoids warnings about temporary creations of arrays for function call by compiler
    allocate(tmp_local(ng),tmp_all_local(ng,0:NPROC-1),stat=ier)
    if (ier /= 0) stop 'error allocating array tmp_local'
    tmp_local(:) = xi_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    xi_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = eta_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    eta_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = gamma_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    gamma_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = final_distance_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    final_distance_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = x_found_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    x_found_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = y_found_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    y_found_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = z_found_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    z_found_source_all(1:ng,:) = tmp_all_local(:,:)

    deallocate(tmp_local,tmp_all_local)

    ! this is executed by main process only
    if (myrank == 0) then

      ! check that the gather operation went well
      if (any(ispec_selected_source_all(1:ng,:) == -1)) call exit_MPI(myrank,'gather operation failed for source')

      ! loop on all the sources
      do is = 1,ng
        ! global source index (between 1 and NSOURCES)
        isource = ns + is - 1

        ! loop on all the results to determine the best slice
        distmin_not_squared = HUGEVAL
        do iprocloop = 0,NPROC-1
          if (final_distance_source_all(is,iprocloop) < distmin_not_squared) then
            ! minimum distance
            distmin_not_squared = final_distance_source_all(is,iprocloop)

            ! stores best source information
            ! slice and element
            islice_selected_source(isource) = iprocloop
            ispec_selected_source(isource) = ispec_selected_source_all(is,iprocloop)
            ! source positioning (inside element)
            xi_source(isource) = xi_source_all(is,iprocloop)
            eta_source(isource) = eta_source_all(is,iprocloop)
            gamma_source(isource) = gamma_source_all(is,iprocloop)
            ! location
            x_found_source(isource) = x_found_source_all(is,iprocloop)
            y_found_source(isource) = y_found_source_all(is,iprocloop)
            z_found_source(isource) = z_found_source_all(is,iprocloop)
            ! domain (elastic/acoustic/poroelastic)
            idomain(isource) = idomain_all(is,iprocloop)
          endif
        enddo
        ! stores minimum distance to desired location
        final_distance_source(isource) = distmin_not_squared
      enddo
    endif !myrank
  enddo ! ngather

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
          write(IMAIN,*) '         z: ',depth(isource),' km'
        else
          write(IMAIN,*) '         depth: ',depth(isource),' km'
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
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
        if (final_distance_source(isource) > 3000.d0) then
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

    ! sets new utm coordinates for best locations
    utm_x_source(:) = x_found_source(:)
    utm_y_source(:) = y_found_source(:)

  endif     ! end of section executed by main process only

  ! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_source,NSOURCES)
  call bcast_all_i(ispec_selected_source,NSOURCES)
  call bcast_all_dp(xi_source,NSOURCES)
  call bcast_all_dp(eta_source,NSOURCES)
  call bcast_all_dp(gamma_source,NSOURCES)
  call bcast_all_dp(utm_x_source,NSOURCES)
  call bcast_all_dp(utm_y_source,NSOURCES)

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

