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

!----
!----  locate_source finds the correct position of the source
!----

  subroutine locate_source(ibool,NSOURCES,myrank,NSPEC_AB,NGLOB_AB,NGNOD, &
                           xstore,ystore,zstore,xigll,yigll,zigll,NPROC, &
                           tshift_src,min_tshift_src_original,yr,jda,ho,mi,utm_x_source,utm_y_source, &
                           DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                           islice_selected_source,ispec_selected_source, &
                           xi_source,eta_source,gamma_source, &
                           nu_source,iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
                           ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                           num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  use constants

  use specfem_par, only: USE_FORCE_POINT_SOURCE,USE_RICKER_TIME_FUNCTION, &
      UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
      factor_force_source,comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP, &
      user_source_time_function,NSTEP_STF,NSOURCES_STF,EXTERNAL_STF,USE_TRICK_FOR_BETTER_PRESSURE

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NSPEC_AB,NGLOB_AB,NSOURCES,NGNOD

  double precision,intent(in) :: DT

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  integer,intent(in) :: myrank

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: xstore,ystore,zstore

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic

  integer,intent(inout) :: yr,jda,ho,mi

  double precision,dimension(NSOURCES),intent(inout) :: tshift_src
  double precision,intent(inout) :: min_tshift_src_original
  double precision :: sec

  integer iprocloop

  integer i,j,k,ispec,iglob,isource
  integer imin,imax,jmin,jmax,kmin,kmax
  integer iproc(1)

  double precision, dimension(NSOURCES),intent(inout) :: utm_x_source,utm_y_source
  double precision dist_squared
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta

  ! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

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

  double precision, dimension(3,3,NGATHER_SOURCES,0:NPROC-1) :: nu_source_all

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
  double precision, dimension(3,3,NSOURCES),intent(inout) :: nu_source

  double precision, dimension(NSOURCES) :: x_found_source,y_found_source,z_found_source
  double precision, dimension(NSOURCES) :: elevation
  double precision :: distmin_squared,distmin_not_squared

  integer, dimension(:), allocatable :: tmp_i_local
  integer, dimension(:,:),allocatable :: tmp_i_all_local

  ! for surface locating and normal computing with external mesh
  integer :: pt0_ix,pt0_iy,pt0_iz,pt1_ix,pt1_iy,pt1_iz,pt2_ix,pt2_iy,pt2_iz
  real(kind=CUSTOM_REAL), dimension(3) :: u_vector,v_vector,w_vector

  logical, dimension(NGLOB_AB),intent(in) :: iglob_is_surface_external_mesh
  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_surface_external_mesh

  integer,intent(in) :: num_free_surface_faces
  integer, dimension(num_free_surface_faces),intent(in) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces),intent(in) :: free_surface_ijk

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
  if (USE_DISTANCE_CRITERION) then
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

    ! orientation consistent with the UTM projection
    ! East
    nu_source(1,1,isource) = 1.d0
    nu_source(1,2,isource) = 0.d0
    nu_source(1,3,isource) = 0.d0
    ! North
    nu_source(2,1,isource) = 0.d0
    nu_source(2,2,isource) = 1.d0
    nu_source(2,3,isource) = 0.d0
    ! Vertical
    nu_source(3,1,isource) = 0.d0
    nu_source(3,2,isource) = 0.d0
    nu_source(3,3,isource) = 1.d0

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

    do ispec=1,NSPEC_AB

      ! exclude elements that are too far from target
      if (USE_DISTANCE_CRITERION) then
        iglob = ibool(MIDX,MIDY,MIDZ,ispec)
        dist_squared = (x_target_source - dble(xstore(iglob)))**2 &
                     + (y_target_source - dble(ystore(iglob)))**2 &
                     + (z_target_source - dble(zstore(iglob)))**2
        if (dist_squared > typical_size_squared) cycle
      endif

      ! define the interval in which we look for points
      if (USE_FORCE_POINT_SOURCE) then
        !imin = 1
        !imax = NGLLX

        !jmin = 1
        !jmax = NGLLY

        !kmin = 1
        !kmax = NGLLZ
        !! VM VM exclude edges to ensure this point is not shared with other elements
        !! otherwise a error location on source can occur in the case of a force point source (FORCESOLUTION)
        imin = 2
        imax = NGLLX - 1

        jmin = 2
        jmax = NGLLY - 1

        kmin = 2
        kmax = NGLLZ - 1
      else
        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        imin = 2
        imax = NGLLX - 1

        jmin = 2
        jmax = NGLLY - 1

        kmin = 2
        kmax = NGLLZ - 1
      endif

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

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

              ! compute final distance between asked and found (converted to km)
              final_distance_source(isource) = dsqrt((x_target_source-x_found_source(isource))**2 + &
                (y_target_source-y_found_source(isource))**2 + (z_target_source-z_found_source(isource))**2)

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

    ! get normal to the face of the hexahedra if receiver is on the surface
    if ((.not. SOURCES_CAN_BE_BURIED) .and. .not. (ispec_selected_source(isource) == 0)) then

      ! note: at this point, xi_source,.. are in range [1.0d0,NGLLX/Y/Z] for point sources only,
      !            for non-point sources the range is limited to [2.0d0,NGLLX/Y/Z - 1]
      if (.not. USE_FORCE_POINT_SOURCE) call exit_MPI(myrank,'error locate source: no point source at surface')

      ! initialize indices
      pt0_ix = -1
      pt0_iy = -1
      pt0_iz = -1
      pt1_ix = -1
      pt1_iy = -1
      pt1_iz = -1
      pt2_ix = -1
      pt2_iy = -1
      pt2_iz = -1

      ! we get two vectors of the face (three points) to compute the normal
      if (nint(xi_source(isource)) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(1,2,2,ispec_selected_source(isource)))) then
        pt0_ix = 1
        pt0_iy = NGLLY
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = 1
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif
      if (nint(xi_source(isource)) == NGLLX .and. &
         iglob_is_surface_external_mesh(ibool(NGLLX,2,2,ispec_selected_source(isource)))) then
        pt0_ix = NGLLX
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = NGLLX
        pt1_iy = NGLLY
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = 1
        pt2_iz = NGLLZ
      endif
      if (nint(eta_source(isource)) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,1,2,ispec_selected_source(isource)))) then
        pt0_ix = 1
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = NGLLX
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = 1
        pt2_iy = 1
        pt2_iz = NGLLZ
      endif
      if (nint(eta_source(isource)) == NGLLY .and. &
         iglob_is_surface_external_mesh(ibool(2,NGLLY,2,ispec_selected_source(isource)))) then
        pt0_ix = NGLLX
        pt0_iy = NGLLY
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = NGLLY
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif
      if (nint(gamma_source(isource)) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,2,1,ispec_selected_source(isource)))) then
        pt0_ix = NGLLX
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = NGLLY
        pt2_iz = 1
      endif
      if (nint(gamma_source(isource)) == NGLLZ .and. &
         iglob_is_surface_external_mesh(ibool(2,2,NGLLZ,ispec_selected_source(isource)))) then
        pt0_ix = 1
        pt0_iy = 1
        pt0_iz = NGLLZ
        pt1_ix = NGLLX
        pt1_iy = 1
        pt1_iz = NGLLZ
        pt2_ix = 1
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif

      if (pt0_ix < 0 .or. pt0_iy < 0 .or. pt0_iz < 0 .or. &
         pt1_ix < 0 .or. pt1_iy < 0 .or. pt1_iz < 0 .or. &
         pt2_ix < 0 .or. pt2_iy < 0 .or. pt2_iz < 0) then
        call exit_mpi(myrank,'error in computing normal for sources.')
      endif

      u_vector(1) = xstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_source(isource))) &
         - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      u_vector(2) = ystore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_source(isource))) &
         - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      u_vector(3) = zstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_source(isource))) &
         - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      v_vector(1) = xstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_source(isource))) &
         - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      v_vector(2) = ystore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_source(isource))) &
         - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      v_vector(3) = zstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_source(isource))) &
         - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))

      ! cross product
      w_vector(1) = u_vector(2)*v_vector(3) - u_vector(3)*v_vector(2)
      w_vector(2) = u_vector(3)*v_vector(1) - u_vector(1)*v_vector(3)
      w_vector(3) = u_vector(1)*v_vector(2) - u_vector(2)*v_vector(1)

      ! normalize vector w
      w_vector(:) = w_vector(:)/sqrt(w_vector(1)**2+w_vector(2)**2+w_vector(3)**2)

      ! build the two other vectors for a direct base: we normalize u, and v=w^u
      u_vector(:) = u_vector(:)/sqrt(u_vector(1)**2+u_vector(2)**2+u_vector(3)**2)
      v_vector(1) = w_vector(2)*u_vector(3) - w_vector(3)*u_vector(2)
      v_vector(2) = w_vector(3)*u_vector(1) - w_vector(1)*u_vector(3)
      v_vector(3) = w_vector(1)*u_vector(2) - w_vector(2)*u_vector(1)

      ! build rotation matrice nu for seismograms
      ! East (u)
      nu_source(1,1,isource) = u_vector(1)
      nu_source(1,2,isource) = v_vector(1)
      nu_source(1,3,isource) = w_vector(1)
      ! North (v)
      nu_source(2,1,isource) = u_vector(2)
      nu_source(2,2,isource) = v_vector(2)
      nu_source(2,3,isource) = w_vector(2)
      ! Vertical (w)
      nu_source(3,1,isource) = u_vector(3)
      nu_source(3,2,isource) = v_vector(3)
      nu_source(3,3,isource) = w_vector(3)

    endif ! of if (.not. SOURCES_CAN_BE_BURIED)

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
        ! (useful if user gives a source outside the mesh for instance)
        if (xi > 1.d0) xi = 1.d0
        if (xi < -1.d0) xi = -1.d0
        if (eta > 1.d0) eta = 1.d0
        if (eta < -1.d0) eta = -1.d0
        if (gamma > 1.d0) gamma = 1.d0
        if (gamma < -1.d0) gamma = -1.d0

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

      ! compute final distance between asked and found (converted to km)
      final_distance_source(isource) = dsqrt((x_target_source-x_found_source(isource))**2 + &
        (y_target_source-y_found_source(isource))**2 + (z_target_source-z_found_source(isource))**2)
    else

      ! takes initial GLL point guess and uses actual location interpolators, in range [-1,1]
      ! note: xi/eta/gamma will be in range [-1,1]
      xi_source(isource) = xigll(ix_initial_guess_source)
      eta_source(isource) = yigll(iy_initial_guess_source)
      gamma_source(isource) = zigll(iz_initial_guess_source)

    endif ! USE_BEST_LOCATION_FOR_SOURCE

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

    do i=1,3
      do j=1,3
        tmp_local(:) = nu_source(i,j,ns:ne)
        call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
        nu_source_all(i,j,1:ng,:) = tmp_all_local(:,:)
      enddo
    enddo
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
            ! orientation
            nu_source(:,:,isource) = nu_source_all(:,:,is,iprocloop)
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
        if (EXTERNAL_STF) then
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
              write(IMAIN,*) '  lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
              write(IMAIN,*) '  lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)

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

