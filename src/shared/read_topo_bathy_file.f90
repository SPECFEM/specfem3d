!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

  subroutine read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO)

! reads topography and bathymetry file

  use constants

  implicit none

  ! use integer array to store topography values
  integer :: NX_TOPO,NY_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

  ! local parameters
  integer :: ix,iy,ier

  ! initializes
  itopo_bathy(:,:) = 0

  ! opens file
  open(unit=13,file=trim(TOPO_FILE),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'error opening topography file: ',trim(TOPO_FILE)
    stop 'error opening topography file'
  endif

  ! reads in values
  do iy = 1,NY_TOPO
    do ix = 1,NX_TOPO
      read(13,*) itopo_bathy(ix,iy)
    enddo
  enddo
  close(13)

  end subroutine read_topo_bathy_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_bathy_elevation(x_target,y_target,target_elevation, &
                                      itopo_bathy,NX_TOPO,NY_TOPO)

! finds elevation from topography file

  use constants

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target

  real(kind=CUSTOM_REAL),intent(out) :: target_elevation

  integer,intent(in) :: NX_TOPO,NY_TOPO
  integer, dimension(NX_TOPO,NY_TOPO),intent(in) :: itopo_bathy

  ! local parameters
  double precision :: xval,yval,long,lat
  double precision :: long_corner,lat_corner,ratio_xi,ratio_eta
  integer :: icornerlong,icornerlat

  ! get coordinates of current point
  xval = dble(x_target)
  yval = dble(y_target)

  ! project x and y in UTM back to long/lat since topo file is in long/lat
  call utm_geo(long,lat,xval,yval,IUTM2LONGLAT)

  ! get coordinate of corner in bathy/topo model
  icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
  icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

  ! avoid edge effects and extend with identical point if outside model
  if (icornerlong < 1) icornerlong = 1
  if (icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
  if (icornerlat < 1) icornerlat = 1
  if (icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

  ! compute coordinates of corner
  long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
  lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

  ! compute ratio for interpolation
  ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
  ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

  ! avoid edge effects
  if (ratio_xi < 0.d0) ratio_xi = 0.d0
  if (ratio_xi > 1.d0) ratio_xi = 1.d0
  if (ratio_eta < 0.d0) ratio_eta = 0.d0
  if (ratio_eta > 1.d0) ratio_eta = 1.d0

  ! interpolate elevation at current point
  target_elevation = &
        itopo_bathy(icornerlong,icornerlat) * (1.d0-ratio_xi)*(1.d0-ratio_eta) + &
        itopo_bathy(icornerlong+1,icornerlat) * ratio_xi*(1.d0-ratio_eta) + &
        itopo_bathy(icornerlong+1,icornerlat+1) * ratio_xi*ratio_eta + &
        itopo_bathy(icornerlong,icornerlat+1) * (1.d0-ratio_xi)*ratio_eta

  end subroutine get_topo_bathy_elevation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_elevation_free(x_target,y_target,target_elevation,target_distmin, &
                                     NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                     num_free_surface_faces,free_surface_ispec,free_surface_ijk)

! get approximate topography elevation at source long/lat coordinates

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,HUGEVAL,TINYVAL,MIDX,MIDY,MIDZ,USE_DISTANCE_CRITERION_TOPO

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target

  real(kind=CUSTOM_REAL),intent(out) :: target_elevation
  real(kind=CUSTOM_REAL),intent(out) :: target_distmin

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: xstore,ystore,zstore

  ! free surface
  integer,intent(in) :: num_free_surface_faces
  integer, dimension(num_free_surface_faces),intent(in) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces),intent(in) :: free_surface_ijk

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(4) :: elevation_node,dist_node
  real(kind=CUSTOM_REAL),dimension(4) :: weight
  real(kind=CUSTOM_REAL) :: distmin,dist
  real(kind=CUSTOM_REAL) :: dist_node_min,norm

  integer :: iface,i,j,k,ispec,iglob,igll,jgll,kgll,ijk
  integer :: iselected,jselected,iface_selected
  integer :: inode,iadjust,jadjust
  integer :: ilocmin(1)

  real(kind=CUSTOM_REAL) :: typical_size
  logical :: located_target

  ! initialize
  target_elevation = 0.0_CUSTOM_REAL
  target_distmin = HUGEVAL

  if (num_free_surface_faces > 0) then

    ! computes typical size of elements at the surface (uses first element for estimation)
    if (USE_DISTANCE_CRITERION_TOPO) then
      ispec = free_surface_ispec(1)
      typical_size =  (xstore(ibool(1,1,1,ispec)) - xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2 &
                    + (ystore(ibool(1,1,1,ispec)) - ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2
      ! use 10 times the distance as a criterion for point detection
      typical_size = 10.0_CUSTOM_REAL * typical_size
    endif

    ! flag to check that we located at least one target element
    located_target = .false.

    !   set distance to huge initial value
    distmin = HUGEVAL
    iselected = 2
    jselected = 2
    iface_selected = 1

    ! loops over all free surface faces
    do iface = 1,num_free_surface_faces
      ispec = free_surface_ispec(iface)

      ! exclude elements that are too far from target
      if (USE_DISTANCE_CRITERION_TOPO) then
        iglob = ibool(MIDX,MIDY,MIDZ,ispec)
        dist = (x_target - xstore(iglob))**2 + (y_target - ystore(iglob))**2
        if (dist > typical_size) cycle
      endif

      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do j = 2,NGLLY - 1
        do i = 2,NGLLX - 1

          igll = free_surface_ijk(1,(j-1)*NGLLY+i,iface)
          jgll = free_surface_ijk(2,(j-1)*NGLLY+i,iface)
          kgll = free_surface_ijk(3,(j-1)*NGLLY+i,iface)

          iglob = ibool(igll,jgll,kgll,ispec)

          ! distance (squared) to target
          dist = ( x_target - xstore(iglob) )**2 + &
                 ( y_target - ystore(iglob) )**2

          ! keep this point if it is closer to the receiver
          if (dist < distmin) then
            distmin = dist
            iface_selected = iface
            iselected = i
            jselected = j
            ! elevation (given in z - coordinate)
            target_elevation = zstore(iglob)
            located_target = .true.
          endif
        enddo
      enddo
    enddo

    ! if we have not located a target element, the point is not in this slice
    ! therefore use first element only for fictitious iterative search
    if (.not. located_target) then
      iselected = 2
      jselected = 2
      iface_selected = 1
    endif

    !  weighted mean at current point of topography elevation of the four closest nodes
    !  set distance to huge initial value
    distmin = HUGEVAL

    ! for example: jselected = 2, iselected = 2
    !     -> j = 2,3 ; i = 2,3
    !
    !     .     .       .
    !     |     |       |
    !     x --- x ----- x -- ..
    !     |     |       |
    !     |     |       |
    !     x --- o ----- x -- ..
    !     |     |(2,2)  |
    !     x --- x ----- x -- ..
    !
    ! evaluates distances & elevation in 4 quads around selected point (iselected,jselected)
    do j = jselected,jselected+1
      do i = iselected,iselected+1
        ! distances to target
        dist_node(:) = HUGEVAL
        inode = 0
        do jadjust = 0,1
          do iadjust = 0,1
            ! for example: j = 2, i = 2
            !   -> ijk_00 = 1 * NGLLY + 2  = 2 + 5
            !   -> ijk_01 = 1 * NGLLY + 1  = 1 + 5
            !   -> ijk_10 = 0 * NGLLY + 2  = 2
            !   -> ijk_11 = 0 * NGLLY + 1  = 1
            ijk = (j-jadjust-1)*NGLLY + i - iadjust

            ispec = free_surface_ispec(iface_selected)

            igll = free_surface_ijk(1,ijk,iface_selected)
            jgll = free_surface_ijk(2,ijk,iface_selected)
            kgll = free_surface_ijk(3,ijk,iface_selected)

            iglob = ibool(igll,jgll,kgll,ispec)

            ! stores node infos
            inode = inode + 1
            elevation_node(inode) = zstore(iglob)
            dist_node(inode) = sqrt( (x_target - xstore(iglob))**2 + (y_target - ystore(iglob))**2)
          enddo
        enddo

        ! weighted elevation
        dist = sum( dist_node(:) )
        if (dist < distmin) then

          ! sets new minimum distance (of all 4 closest nodes)
          distmin = dist
          target_distmin = distmin

          ! interpolates elevation
          if (dist > TINYVAL) then
            !original weighting: by distance, which means larger distance -> stronger weight
            !                    but we would prefer closer points having more influence.
            !                    therefore, will use an inverse distance weighting instead...
            !
            !target_elevation =  (dist_node(1)/dist) * elevation_node(1) + &
            !                    (dist_node(2)/dist) * elevation_node(2) + &
            !                    (dist_node(3)/dist) * elevation_node(3) + &
            !                    (dist_node(4)/dist) * elevation_node(4)

            ! gets minimum distance value & index
            ilocmin = minloc(dist_node(:))
            dist_node_min = dist_node(ilocmin(1))

            ! checks if a node is almost exactly at target location
            if (dist_node_min < TINYVAL) then
              ! takes elevation of node
              target_elevation = elevation_node(ilocmin(1))
            else
              ! inverse distance weighting
              ! (closer points have higher weight)
              do k = 1,4
                if (dist_node(k) > 0.0_CUSTOM_REAL) then
                  weight(k) = 1.0_CUSTOM_REAL / dist_node(k)
                else
                  weight(k) = 1.e10  ! very large weight for point on target, will dominate interpolation
                endif
              enddo

              ! normalize weights: w_i = w_i / sum(w_i)
              norm = sum(weight(:))
              if (norm > TINYVAL) then
                weight(:) = weight(:) / norm
              else
                ! all 4 weights almost zero, meaning distances are all very large; uses equal weighting
                weight(:) = 1.0 / 4.0
              endif
              ! interpolation
              target_elevation = weight(1)*elevation_node(1) &
                               + weight(2)*elevation_node(2) &
                               + weight(3)*elevation_node(3) &
                               + weight(4)*elevation_node(4)
            endif
          else
            stop 'Error summed distance to node is zero'
          endif
        endif

      enddo
    enddo

  endif

  end subroutine get_topo_elevation_free

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_elevation_free_closest(x_target,y_target,target_elevation,target_distmin, &
                                             NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                             num_free_surface_faces,free_surface_ispec,free_surface_ijk)

! get approximate topography elevation at long/lat coordinates from closest point

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,HUGEVAL,MIDX,MIDY,MIDZ,USE_DISTANCE_CRITERION_TOPO

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target

  real(kind=CUSTOM_REAL),intent(out) :: target_elevation
  real(kind=CUSTOM_REAL),intent(out) :: target_distmin

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: xstore,ystore,zstore

  ! free surface
  integer,intent(in) :: num_free_surface_faces
  integer, dimension(num_free_surface_faces),intent(in) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces),intent(in) :: free_surface_ijk

  ! local parameters
  real(kind=CUSTOM_REAL) :: distmin,dist

  integer :: iface,i,ispec,iglob,igll,jgll,kgll

  real(kind=CUSTOM_REAL) :: typical_size
  logical :: located_target

  ! initialize
  target_elevation = 0.0_CUSTOM_REAL
  target_distmin = HUGEVAL

  if (num_free_surface_faces > 0) then

    ! computes typical size of elements at the surface (uses first element for estimation)
    if (USE_DISTANCE_CRITERION_TOPO) then
      ispec = free_surface_ispec(1)
      typical_size =  (xstore(ibool(1,1,1,ispec)) - xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2 &
                    + (ystore(ibool(1,1,1,ispec)) - ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2
      ! use 10 times the distance as a criterion for point detection
      typical_size = 10.0_CUSTOM_REAL * typical_size
    endif

    ! flag to check that we located at least one target element
    located_target = .false.

    !   set distance to huge initial value
    distmin = HUGEVAL

    ! loops over all free surface faces
    do iface = 1,num_free_surface_faces
      ispec = free_surface_ispec(iface)

      ! excludes elements that are too far from target
      if (USE_DISTANCE_CRITERION_TOPO) then
        iglob = ibool(MIDX,MIDY,MIDZ,ispec)
        dist = (x_target - xstore(iglob))**2 + (y_target - ystore(iglob))**2
        if (dist > typical_size) cycle
      endif

      ! loop only on points inside the element
      do i = 1,NGLLSQUARE
        igll = free_surface_ijk(1,i,iface)
        jgll = free_surface_ijk(2,i,iface)
        kgll = free_surface_ijk(3,i,iface)

        iglob = ibool(igll,jgll,kgll,ispec)

        ! distance (squared) to target
        dist = ( x_target - xstore(iglob) )**2 + &
               ( y_target - ystore(iglob) )**2

        ! keep this point if it is closer to the receiver
        if (dist < distmin) then
          distmin = dist

          ! elevation (given in z - coordinate)
          target_elevation = zstore(iglob)
          target_distmin = dist
          located_target = .true.
        endif
      enddo
    enddo

    ! if we have not located a target element, the point is not in this slice
    ! therefore use first element only for fictitious iterative search
    if (.not. located_target) then
      !stop 'error: point was not located in get_elevation_closest()'
      ! takes first point for estimation
      iglob = ibool(1,1,1,ispec)
      ! elevation (given in z - coordinate)
      target_elevation = zstore(iglob)
      target_distmin = ( x_target - xstore(iglob) )**2 + ( y_target - ystore(iglob) )**2
      located_target = .true.
    endif

  endif

  end subroutine get_topo_elevation_free_closest

