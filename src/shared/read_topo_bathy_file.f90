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

  subroutine read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO)

! reads topography and bathymetry file

  implicit none

  include "constants.h"

  ! use integer array to store topography values
  integer :: NX_TOPO,NY_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

  ! local parameters
  integer :: ix,iy,ier

  ! initializes
  itopo_bathy(:,:) = 0

  ! opens file
  open(unit=13,file=trim(TOPO_FILE),status='old',action='read',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening topography file: ',trim(TOPO_FILE)
    stop 'error opening topography file'
  endif

  ! reads in values
  do iy=1,NY_TOPO
    do ix=1,NX_TOPO
      read(13,*) itopo_bathy(ix,iy)
    enddo
  enddo
  close(13)

  end subroutine read_topo_bathy_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_bathy_elevation(x_target,y_target,target_elevation, &
                                itopo_bathy,NX_TOPO,NY_TOPO, &
                                UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION)

! finds elevation from topography file

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target

  real(kind=CUSTOM_REAL),intent(out) :: target_elevation

  integer :: NX_TOPO,NY_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION

  ! local parameters
  double precision :: xval,yval,long,lat
  double precision :: long_corner,lat_corner,ratio_xi,ratio_eta
  integer :: icornerlong,icornerlat

  ! get coordinates of current point
  xval = dble(x_target)
  yval = dble(y_target)

  ! project x and y in UTM back to long/lat since topo file is in long/lat
  call utm_geo(long,lat,xval,yval,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

  ! get coordinate of corner in bathy/topo model
  icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
  icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

  ! avoid edge effects and extend with identical point if outside model
  if(icornerlong < 1) icornerlong = 1
  if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
  if(icornerlat < 1) icornerlat = 1
  if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

  ! compute coordinates of corner
  long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
  lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

  ! compute ratio for interpolation
  ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
  ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

  ! avoid edge effects
  if(ratio_xi < 0.) ratio_xi = 0.
  if(ratio_xi > 1.) ratio_xi = 1.
  if(ratio_eta < 0.) ratio_eta = 0.
  if(ratio_eta > 1.) ratio_eta = 1.

  ! interpolate elevation at current point
  target_elevation = &
        itopo_bathy(icornerlong,icornerlat) * (1.0_CUSTOM_REAL-ratio_xi)*(1.0_CUSTOM_REAL-ratio_eta) + &
        itopo_bathy(icornerlong+1,icornerlat) * ratio_xi*(1.0_CUSTOM_REAL-ratio_eta) + &
        itopo_bathy(icornerlong+1,icornerlat+1) * ratio_xi*ratio_eta + &
        itopo_bathy(icornerlong,icornerlat+1) * (1.0_CUSTOM_REAL-ratio_xi)*ratio_eta

  end subroutine get_topo_bathy_elevation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_elevation_free(x_target,y_target,target_elevation,target_distmin, &
                                    NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    num_free_surface_faces,free_surface_ispec,free_surface_ijk)

! get approximate topography elevation at source long/lat coordinates

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target

  real(kind=CUSTOM_REAL),intent(out) :: target_elevation
  real(kind=CUSTOM_REAL),intent(out) :: target_distmin

  integer :: NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! free surface
  integer :: num_free_surface_faces
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(4) :: elevation_node,dist_node
  real(kind=CUSTOM_REAL) :: distmin,dist

  integer :: iface,i,j,ispec,iglob,igll,jgll,kgll
  integer :: iselected,jselected,iface_selected
  integer :: inode,iadjust,jadjust

  ! faster element search
  logical,parameter :: USE_DISTANCE_CRITERION = .true.
  integer,parameter :: MIDX = (NGLLX+1)/2
  integer,parameter :: MIDY = (NGLLY+1)/2
  integer,parameter :: MIDZ = (NGLLZ+1)/2

  real(kind=CUSTOM_REAL) :: typical_size
  logical :: located_target

  ! initialize
  target_elevation = 0.0_CUSTOM_REAL
  target_distmin = HUGEVAL


  if(num_free_surface_faces > 0) then

    ! computes typical size of elements at the surface (uses first element for estimation)
    if( USE_DISTANCE_CRITERION ) then
      ispec = free_surface_ispec(1)
      typical_size =  (xstore(ibool(1,1,1,ispec)) - xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2 &
                    + (ystore(ibool(1,1,1,ispec)) - ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2
      ! use 10 times the distance as a criterion for point detection
      typical_size = 10. * typical_size
    endif

    ! flag to check that we located at least one target element
    located_target = .false.

    !   set distance to huge initial value
    distmin = HUGEVAL
    iselected = 2
    jselected = 2
    iface_selected = 1

    ! loops over all free surface faces
    do iface=1,num_free_surface_faces
      ispec = free_surface_ispec(iface)

      ! exclude elements that are too far from target
      if( USE_DISTANCE_CRITERION ) then
        iglob = ibool(MIDX,MIDY,MIDZ,ispec)
        dist = (x_target - xstore(iglob))**2 + (y_target - ystore(iglob))**2
        if( dist > typical_size ) cycle
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
          if(dist < distmin) then
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
    if(.not. located_target) then
      iselected = 2
      jselected = 2
      iface_selected = 1
    endif

    !  weighted mean at current point of topography elevation of the four closest nodes
    !  set distance to huge initial value
    distmin = HUGEVAL
    do j=jselected,jselected+1
      do i=iselected,iselected+1
        ! distances to target
        dist_node(:) = HUGEVAL
        inode = 0
        do jadjust=0,1
          do iadjust= 0,1
            ispec = free_surface_ispec(iface_selected)
            igll = free_surface_ijk(1,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
            jgll = free_surface_ijk(2,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
            kgll = free_surface_ijk(3,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
            iglob = ibool(igll,jgll,kgll,ispec)

            ! stores node infos
            inode = inode + 1
            elevation_node(inode) = zstore(iglob)
            dist_node(inode) = sqrt( (x_target - xstore(iglob))**2 + (y_target - ystore(iglob))**2 )
          enddo
        enddo

        ! weighted elevation
        dist = sum( dist_node(:) )
        if(dist < distmin) then

          ! sets new minimum distance (of all 4 closest nodes)
          distmin = dist
          target_distmin = distmin

          ! interpolates elevation
          if( dist > TINYVAL ) then
            target_elevation =  (dist_node(1)/dist)*elevation_node(1) + &
                                (dist_node(2)/dist)*elevation_node(2) + &
                                (dist_node(3)/dist)*elevation_node(3) + &
                                (dist_node(4)/dist)*elevation_node(4)
          else
            stop 'error summed distance to node is zero'
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

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target

  real(kind=CUSTOM_REAL),intent(out) :: target_elevation
  real(kind=CUSTOM_REAL),intent(out) :: target_distmin

  integer :: NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! free surface
  integer :: num_free_surface_faces
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  ! local parameters
  real(kind=CUSTOM_REAL) :: distmin,dist

  integer :: iface,i,ispec,iglob,igll,jgll,kgll

  ! faster element search
  logical,parameter :: USE_DISTANCE_CRITERION = .true.
  integer,parameter :: MIDX = (NGLLX+1)/2
  integer,parameter :: MIDY = (NGLLY+1)/2
  integer,parameter :: MIDZ = (NGLLZ+1)/2

  real(kind=CUSTOM_REAL) :: typical_size
  logical :: located_target

  ! initialize
  target_elevation = 0.0_CUSTOM_REAL
  target_distmin = HUGEVAL


  if(num_free_surface_faces > 0) then

    ! computes typical size of elements at the surface (uses first element for estimation)
    if( USE_DISTANCE_CRITERION ) then
      ispec = free_surface_ispec(1)
      typical_size =  (xstore(ibool(1,1,1,ispec)) - xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2 &
                    + (ystore(ibool(1,1,1,ispec)) - ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec)))**2
      ! use 10 times the distance as a criterion for point detection
      typical_size = 10. * typical_size
    endif

    ! flag to check that we located at least one target element
    located_target = .false.

    !   set distance to huge initial value
    distmin = HUGEVAL

    ! loops over all free surface faces
    do iface=1,num_free_surface_faces
      ispec = free_surface_ispec(iface)

      ! excludes elements that are too far from target
      if( USE_DISTANCE_CRITERION ) then
        iglob = ibool(MIDX,MIDY,MIDZ,ispec)
        dist = (x_target - xstore(iglob))**2 + (y_target - ystore(iglob))**2
        if( dist > typical_size ) cycle
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
        if(dist < distmin) then
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
    if(.not. located_target) then
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

