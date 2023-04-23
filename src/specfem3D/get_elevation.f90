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


!--------------------------------------------------------------------------------------------------------------------
! get z target coordinate, depending on the topography
!--------------------------------------------------------------------------------------------------------------------

  subroutine get_elevation_and_z_coordinate_all(npoints,plon,plat,pbur,utm_x,utm_y,elevation,x_target,y_target,z_target)

! routine for an array of points

  use constants
  use specfem_par, only: ibool,myrank,NSPEC_AB,NGLOB_AB,USE_SOURCES_RECEIVERS_Z, &
                         xstore,ystore,zstore,NPROC,num_free_surface_faces,free_surface_ispec,free_surface_ijk
  implicit none

  integer, intent(in) :: npoints
  double precision, dimension(npoints), intent(in)  :: plon,plat,pbur

  double precision, dimension(npoints), intent(inout) :: utm_x,utm_y,elevation
  double precision, dimension(npoints), intent(inout) :: x_target,y_target,z_target

  ! local parameters
  integer :: ipoin
  double precision :: z,distmin,lon,lat
  double precision, dimension(:), allocatable :: elevation_distmin
  double precision, dimension(:,:), allocatable :: elevation_all,elevation_distmin_all
  real(kind=CUSTOM_REAL) :: xloc,yloc,loc_ele,loc_distmin
  integer :: iproc,ier

  !debug
  !print *,'get elevation: ',npoints

  ! determine x/y from lat/lon
  do ipoin = 1, npoints
    ! point longitude/latitude
    lon = plon(ipoin)
    lat = plat(ipoin)

    ! convert station location to UTM
    call utm_geo(lon,lat,utm_x(ipoin),utm_y(ipoin),ILONGLAT2UTM)

    ! target location of point
    x_target(ipoin) = utm_x(ipoin)
    y_target(ipoin) = utm_y(ipoin)
  enddo

  ! check if z location specified
  if (USE_SOURCES_RECEIVERS_Z) then
    ! burial depth is given as z value directly
    do ipoin = 1, npoints
      ! target location of receiver
      z_target(ipoin) = pbur(ipoin)
    enddo
    ! all done
    return
  endif ! USE_SOURCES_RECEIVERS_Z


  ! needs to find elevation of mesh at point location lat/lon to use burial depth for z coordinate
  !
  ! allocates temporary arrays
  allocate(elevation_distmin(npoints),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2184')
  if (ier /= 0) stop 'Error allocating elevation arrays'
  elevation_distmin(:) = HUGEVAL

  if (myrank == 0) then
    ! only main gathers all
    allocate(elevation_all(npoints,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2185')
    allocate(elevation_distmin_all(npoints,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2186')
    if (ier /= 0) stop 'Error allocating elevation gather arrays'
  else
    allocate(elevation_all(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2187')
    allocate(elevation_distmin_all(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2188')
    if (ier /= 0) stop 'Error allocating elevation gather arrays'
  endif
  elevation_all(:,:) = 0.d0
  elevation_distmin_all(:,:) = HUGEVAL

  ! gets elevation
  do ipoin = 1, npoints
    xloc = x_target(ipoin)
    yloc = y_target(ipoin)

    ! get approximate topography elevation at point x/y coordinates
    call get_topo_elevation_free(xloc,yloc,loc_ele,loc_distmin, &
                                 NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                 num_free_surface_faces,free_surface_ispec,free_surface_ijk)

    elevation(ipoin) = loc_ele
    elevation_distmin(ipoin) = loc_distmin
  enddo
  ! MPI communications to determine the best slice
  call gather_all_dp(elevation,npoints,elevation_all,npoints,NPROC)
  call gather_all_dp(elevation_distmin,npoints,elevation_distmin_all,npoints,NPROC)

  ! main process selects closest
  if (myrank == 0) then
    ! loops over all receivers
    do ipoin = 1,npoints
      ! find closest estimate in all slices
      distmin = HUGEVAL
      do iproc = 0,NPROC-1
        ! closest point
        if (elevation_distmin_all(ipoin,iproc) < distmin) then
          distmin = elevation_distmin_all(ipoin,iproc)
          elevation(ipoin) = elevation_all(ipoin,iproc)

          ! point's Z coordinate
          ! burial depth read in file assumed to be given in m
          z = elevation(ipoin) - pbur(ipoin)

          ! target coordinate, depending on the topography
          z_target(ipoin) = z
        endif
      enddo
    enddo
  endif

  ! free temporary arrays
  deallocate(elevation_distmin,elevation_all,elevation_distmin_all)

  ! main broadcasts to all slices
  call bcast_all_dp(z_target,npoints)

  end subroutine get_elevation_and_z_coordinate_all

!
!-----------------------------------------------------------------------------------------
!

  subroutine get_elevation_and_z_coordinate(lon,lat,utm_x,utm_y,z_target,elevation,bury)

! routine for a single point

  use constants
  use specfem_par, only: USE_SOURCES_RECEIVERS_Z,ibool,myrank,NSPEC_AB,NGLOB_AB, &
                         xstore,ystore,zstore,NPROC,num_free_surface_faces,free_surface_ispec,free_surface_ijk

  implicit none

  double precision,intent(inout) :: lon,lat,bury
  double precision,intent(inout) :: utm_x,utm_y
  double precision,intent(inout) :: z_target,elevation

  ! local parameters
  integer,dimension(1)              :: iproc
  double precision,dimension(1)     :: altitude_rec,distmin_ele
  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all
  real(kind=CUSTOM_REAL)            :: xloc,yloc,loc_ele,loc_distmin

  ! convert station location to UTM
  call utm_geo(lon,lat,utm_x,utm_y,ILONGLAT2UTM)

  xloc = utm_x
  yloc = utm_y

  ! point's Z coordinate
  if (USE_SOURCES_RECEIVERS_Z) then
    ! alternative: burial depth is given as z value directly
    z_target = bury
    ! all done
    return
  endif

  ! get approximate topography elevation at point long/lat coordinates
  call get_topo_elevation_free(xloc,yloc,loc_ele,loc_distmin, &
                               NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                               num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  altitude_rec(1) = loc_ele
  distmin_ele(1)  = loc_distmin

  !  MPI communications to determine the best slice
  call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
  call gather_all_dp(altitude_rec,1,elevation_all,1,NPROC)

  if (myrank == 0) then
    iproc = minloc(distmin_ele_all)
    altitude_rec(1) = elevation_all(iproc(1))
  endif

  call bcast_all_dp(altitude_rec,1)
  elevation = altitude_rec(1)

  ! point's Z coordinate
  ! burial depth read in file given in m
  z_target = elevation - bury

  end subroutine get_elevation_and_z_coordinate

