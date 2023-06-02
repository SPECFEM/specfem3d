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

  subroutine detect_surface(NPROC,nglob,nspec,ibool, &
                            ispec_is_surface, &
                            iglob_is_surface, &
                            nfaces_surface, &
                            num_interfaces_ext_mesh, &
                            max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh, &
                            my_neighbors_ext_mesh, &
                            ibool_interfaces_ext_mesh)

! detects surface (points/elements) of model based upon valence
!
! returns: ispec_is_surface, iglob_is_surface and nfaces_surface

  use constants

  implicit none

  ! global indexing
  integer,intent(in) :: NPROC,nglob,nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  ! surface
  logical, dimension(nspec),intent(inout) :: ispec_is_surface
  logical, dimension(nglob),intent(inout) :: iglob_is_surface
  integer,intent(inout) :: nfaces_surface

  ! MPI partitions
  integer,intent(in) :: num_interfaces_ext_mesh
  integer,intent(in) :: max_nibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer,dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh

  ! local parameters
  integer, dimension(:), allocatable :: valence
  integer :: ispec,i,j,k,iglob,ier

  ! detecting surface points/elements (based on valence check on NGLL points) for external mesh
  allocate(valence(nglob),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1223')
  if (ier /= 0) stop 'error allocate valence array'

  ! initialize surface indices
  ispec_is_surface(:) = .false.
  iglob_is_surface(:) = .false.
  valence(:) = 0

  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if (iglob < 1 .or. iglob > nglob) then
            print *,'error valence iglob:',iglob,i,j,k,ispec
            stop 'error valence'
          endif
          valence(iglob) = valence(iglob) + 1
        enddo
      enddo
    enddo
  enddo

  ! adds contributions from different partitions to valence
  call assemble_MPI_scalar_i_blocking(NPROC,nglob,valence, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

  ! determines spectral elements containing surface points
  do ispec = 1, nspec

    ! loops over GLL points not on edges or corners
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          if ((k == 1 .or. k == NGLLZ) .and. (j /= 1 .and. j /= NGLLY) .and. (i /= 1 .and. i /= NGLLX) .or. &
              (j == 1 .or. j == NGLLY) .and. (k /= 1 .and. k /= NGLLZ) .and. (i /= 1 .and. i /= NGLLX) .or. &
              (i == 1 .or. i == NGLLX) .and. (k /= 1 .and. k /= NGLLZ) .and. (j /= 1 .and. j /= NGLLY)) then
            iglob = ibool(i,j,k,ispec)
            if (valence(iglob) == 1) then
              ! sets surface flags for element and global points
              call ds_set_surface_flags(nspec,ispec_is_surface, &
                                        nglob,iglob_is_surface, &
                                        i,j,k,ispec,ibool)
            endif

          endif
        enddo
      enddo
    enddo

  enddo ! nspec

  ! safety check
  if ((NGLLX <= 2) .or. (NGLLY <= 2) .or. (NGLLZ <= 2)) &
    stop 'Error: invalid NGLL number for surface detection, must be > 2'

  ! counts faces for movies and shakemaps
  nfaces_surface = 0
  do ispec = 1, nspec
    ! takes an inner point of the surface, assuming assuming (2,2,*) is not a corner, i.e., NGLLX > 2, NGLLY > 2
    iglob = ibool(2,2,1,ispec)
    if (iglob_is_surface(iglob)) then
      nfaces_surface = nfaces_surface + 1
    endif
    iglob = ibool(2,2,NGLLZ,ispec)
    if (iglob_is_surface(iglob)) then
      nfaces_surface = nfaces_surface + 1
    endif
    iglob = ibool(2,1,2,ispec)
    if (iglob_is_surface(iglob)) then
      nfaces_surface = nfaces_surface + 1
    endif
    iglob = ibool(2,NGLLY,2,ispec)
    if (iglob_is_surface(iglob)) then
      nfaces_surface = nfaces_surface + 1
    endif
    iglob = ibool(1,2,2,ispec)
    if (iglob_is_surface(iglob)) then
      nfaces_surface = nfaces_surface + 1
    endif
    iglob = ibool(NGLLX,2,2,ispec)
    if (iglob_is_surface(iglob)) then
      nfaces_surface = nfaces_surface + 1
    endif
  enddo

  end subroutine detect_surface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine ds_set_surface_flags(nspec,ispec_is_surface, &
                                  nglob,iglob_is_surface, &
                                  i,j,k,ispec,ibool)

  ! put this into separate subroutine to compile faster, otherwise compilers will try to unroll all do loops

  use constants

  implicit none

  ! global indexing
  integer,intent(in) :: nglob,nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  integer,intent(in) :: i,j,k,ispec

  !   surface flags
  logical, dimension(nspec),intent(inout) :: ispec_is_surface
  logical, dimension(nglob),intent(inout) :: iglob_is_surface

  ! local parameters
  integer :: kk,jj,ii

  ! sets surface flag for element
  ispec_is_surface(ispec) = .true.

  ! sets flags for all GLL points on this face
  if (k == 1 .or. k == NGLLZ) then
    do jj = 1, NGLLY
      do ii = 1, NGLLX
        iglob_is_surface(ibool(ii,jj,k,ispec)) = .true.
      enddo
    enddo
  endif
  if (j == 1 .or. j == NGLLY) then
    do kk = 1, NGLLZ
      do ii = 1, NGLLX
        iglob_is_surface(ibool(ii,j,kk,ispec)) = .true.
      enddo
    enddo
  endif
  if (i == 1 .or. i == NGLLX) then
    do kk = 1, NGLLZ
      do jj = 1, NGLLY
        iglob_is_surface(ibool(i,jj,kk,ispec)) = .true.
      enddo
    enddo
  endif

  end subroutine ds_set_surface_flags

!
!-------------------------------------------------------------------------------------------------
!

  subroutine detect_surface_cross_section(NPROC,nglob,nspec,ibool, &
                                          ispec_is_surface_cross_section, &
                                          iglob_is_surface_cross_section, &
                                          nfaces_surface, &
                                          num_interfaces_ext_mesh, &
                                          max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh, &
                                          ibool_interfaces_ext_mesh, &
                                          x_section,y_section,z_section, &
                                          xstore,ystore,zstore)

! instead of surface of model, this returns cross-section surfaces through model
! at specified x,y,z - coordinates
!
! note: x,y,z coordinates must coincide with the element (outer-)faces, no planes inside elements are taken
!         (this is only a quick & dirty cross-section implementation, no sophisticated interpolation of points considered...)
!
! returns: ispec_is_surface_cross_section, iglob_is_surface_cross_section and nfaces_surface

  use constants

  implicit none

  ! global indexing
  integer,intent(in) :: NPROC,nglob,nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  ! surface
  logical, dimension(nspec),intent(inout) :: ispec_is_surface_cross_section
  logical, dimension(nglob),intent(inout) :: iglob_is_surface_cross_section
  integer,intent(inout) :: nfaces_surface

  ! MPI partitions
  integer,intent(in) :: num_interfaces_ext_mesh
  integer,intent(in) :: max_nibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer,dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh

  ! specified x,y,z - coordinates
  real(kind=CUSTOM_REAL),intent(in) :: x_section,y_section,z_section

  ! mesh global point coordinates
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(6) :: midpoint_faces_x,midpoint_faces_y, &
                                         midpoint_faces_z
  real(kind=CUSTOM_REAL),dimension(6) :: midpoint_dist_x,midpoint_dist_y,midpoint_dist_z
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord_face,ycoord_face,zcoord_face
  real(kind=CUSTOM_REAL) :: min_dist,normal(NDIM)
  integer, dimension(:), allocatable :: valence

  integer :: ispec,i,j,k,iglob,ier,countval
  integer :: iface,icorner
  logical, dimension(:),allocatable :: ispec_has_points

  ! corners indices of reference cube faces
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/)) ! xmin
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
       reshape((/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/)) !xmax
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
       reshape((/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/)) !ymin
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
       reshape((/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/)) ! ymax
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/)) ! bottom
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
       reshape((/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/)) !top
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
       reshape((/ iface1_corner_ijk,iface2_corner_ijk, &
                  iface3_corner_ijk,iface4_corner_ijk, &
                  iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/)) ! all faces
  integer,dimension(3,6),parameter :: iface_midpoint_ijk = &
             reshape( (/ 1,3,3, NGLLX,3,3, 3,1,3, 3,NGLLY,3, 3,3,1, 3,3,NGLLZ  /),(/3,6/))   ! top

  ! detecting surface points/elements (based on valence check on NGLL points) for external mesh
  allocate(valence(nglob),ispec_has_points(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1224')
  if (ier /= 0) stop 'error allocate valence array'

  ! an estimation of the minimum distance between global points (for an element width)
  min_dist = minval( (xstore(ibool(1,3,3,:)) - xstore(ibool(NGLLX,3,3,:)))**2 &
                  + (ystore(ibool(1,3,3,:)) - ystore(ibool(NGLLX,3,3,:)))**2 &
                  + (zstore(ibool(1,3,3,:)) - zstore(ibool(NGLLX,3,3,:)))**2)
  min_dist = sqrt(min_dist)

  ! initialize surface indices
  ispec_is_surface_cross_section(:) = .false.
  iglob_is_surface_cross_section(:) = .false.
  nfaces_surface  = 0
  valence(:) = 0

  ! sets valence value to one corresponding to process rank  for points on cross-sections
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)

          ! x cross-section
          if (abs( xstore(iglob) - x_section ) < 0.2*min_dist) then
            ! sets valence to 1 for points on cross-sections
            valence(iglob) = myrank+1
          endif

          ! y cross-section
          if (abs( ystore(iglob) - y_section ) < 0.2*min_dist) then
            ! sets valence to 1 for points on cross-sections
            valence(iglob) = myrank+1
          endif

          ! z cross-section
          if (abs( zstore(iglob) - z_section ) < 0.2*min_dist) then
            ! sets valence to 1 for points on cross-sections
            valence(iglob) = myrank+1
          endif

        enddo
      enddo
    enddo
  enddo

  ! adds contributions from different partitions to valence
  call assemble_MPI_scalar_i_blocking(NPROC,nglob,valence, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)


  ! determines spectral elements containing surface points
  ! (only counts element outer faces, no planes inside element)
  ispec_has_points(:) = .false.
  countval = 0
  do ispec = 1, nspec

    ! loops over GLL points not on edges or corners, but inside faces
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          iglob = ibool(i,j,k,ispec)

          ! sets flag if element has points
          if (valence(iglob) > 0) ispec_has_points(ispec) = .true.

          ! checks element surfaces for valence points
          if (((k == 1 .or. k == NGLLZ) .and. (j == 2 .and. i == 2)) .or. &
              ((j == 1 .or. j == NGLLY) .and. (k == 2 .and. i == 2)) .or. &
              ((i == 1 .or. i == NGLLX) .and. (k == 2 .and. j == 2))) then

            iglob = ibool(i,j,k,ispec)

            ! considers only points in same process or, if point is shared between two processes,
            ! only with higher process ranks than itself
            if (valence(iglob) == myrank+1 .or. valence(iglob) > 2*(myrank+1)) then
              ! sets surface flags for cross section
              call ds_set_cross_section_flags(nspec,ispec_is_surface_cross_section, &
                                              nglob,iglob_is_surface_cross_section, &
                                              i,j,k,ispec,ibool, &
                                              valence,countval)
            endif

          endif
        enddo
      enddo
    enddo

  enddo ! nspec

  ! tries to find closest face if points are inside
  do ispec = 1,nspec

    ! in case element has still unresolved points in interior,
    ! we take closest element face to cross-section plane
    if (ispec_has_points(ispec)) then

      ! an estimation of the element width
      min_dist = sqrt((xstore(ibool(1,3,3,ispec)) - xstore(ibool(NGLLX,3,3,ispec)))**2 &
                   + (ystore(ibool(1,3,3,ispec)) - ystore(ibool(NGLLX,3,3,ispec)))**2 &
                   + (zstore(ibool(1,3,3,ispec)) - zstore(ibool(NGLLX,3,3,ispec)))**2)

      ! determines element face by minimum distance of midpoints
      midpoint_faces_x(:) = 0.0
      midpoint_faces_y(:) = 0.0
      midpoint_faces_z(:) = 0.0

      do iface=1,6

        ! face corners
        do icorner = 1,NGNOD2D_FOUR_CORNERS
          i = iface_all_corner_ijk(1,icorner,iface)
          j = iface_all_corner_ijk(2,icorner,iface)
          k = iface_all_corner_ijk(3,icorner,iface)

          ! coordinates
          iglob = ibool(i,j,k,ispec)
          xcoord_face(icorner) = xstore(iglob)
          ycoord_face(icorner) = ystore(iglob)
          zcoord_face(icorner) = zstore(iglob)

          ! face midpoint coordinates
          midpoint_faces_x(iface) =  midpoint_faces_x(iface) + xcoord_face(icorner)
          midpoint_faces_y(iface) =  midpoint_faces_y(iface) + ycoord_face(icorner)
          midpoint_faces_z(iface) =  midpoint_faces_z(iface) + zcoord_face(icorner)
        enddo

        midpoint_faces_x(iface) = midpoint_faces_x(iface) / 4.0
        midpoint_faces_y(iface) = midpoint_faces_y(iface) / 4.0
        midpoint_faces_z(iface) = midpoint_faces_z(iface) / 4.0

        ! gets face normal
        normal(:) = 0._CUSTOM_REAL
        call get_element_face_normal(ispec,iface,xcoord_face,ycoord_face,zcoord_face, &
                                     ibool,nspec,nglob,xstore,ystore,zstore, &
                                     normal)

        ! distance to cross-section planes
        midpoint_dist_x(iface) = abs(midpoint_faces_x(iface) - x_section)
        midpoint_dist_y(iface) = abs(midpoint_faces_y(iface) - y_section)
        midpoint_dist_z(iface) = abs(midpoint_faces_z(iface) - z_section)


        ! x cross-section plane
        i = iface_midpoint_ijk(1,iface)
        j = iface_midpoint_ijk(2,iface)
        k = iface_midpoint_ijk(3,iface)
        if (midpoint_dist_x(iface) < 0.5*min_dist .and. &
           valence(ibool(i,j,k,ispec)) /= -1) then
          ! checks face normal points in similar direction as cross-section normal
          if (abs(normal(1)) > 0.6) then
            ! sets surfaces flags
            call ds_set_plane_flags(iface,ispec, &
                                    nspec,ispec_is_surface_cross_section, &
                                    nglob,iglob_is_surface_cross_section, &
                                    ibool,valence)
          endif
        endif

        ! y cross-section plane
        i = iface_midpoint_ijk(1,iface)
        j = iface_midpoint_ijk(2,iface)
        k = iface_midpoint_ijk(3,iface)
        if (midpoint_dist_y(iface) < 0.5*min_dist .and. &
           valence(ibool(i,j,k,ispec)) /= -1) then
          ! checks face normal points in similar direction as cross-section normal
          if (abs(normal(2)) > 0.6) then
            ! sets surfaces flags
            call ds_set_plane_flags(iface,ispec, &
                                    nspec,ispec_is_surface_cross_section, &
                                    nglob,iglob_is_surface_cross_section, &
                                    ibool,valence)
          endif
        endif

        ! z cross-section plane
        i = iface_midpoint_ijk(1,iface)
        j = iface_midpoint_ijk(2,iface)
        k = iface_midpoint_ijk(3,iface)
        if (midpoint_dist_z(iface) < 0.5*min_dist .and. &
           valence(ibool(i,j,k,ispec)) /= -1) then
          ! checks face normal points in similar direction as cross-section normal
          if (abs(normal(3)) > 0.6) then
            ! sets surfaces flags
            call ds_set_plane_flags(iface,ispec, &
                                    nspec,ispec_is_surface_cross_section, &
                                    nglob,iglob_is_surface_cross_section, &
                                    ibool,valence)
          endif
        endif

      enddo ! iface

    endif
  enddo

  ! counts faces for external-mesh movies and shakemaps
  nfaces_surface = 0
  do ispec = 1, nspec
    if (ispec_is_surface_cross_section(ispec)) then
      ! zmin face
      if (iglob_is_surface_cross_section(ibool(2,2,1,ispec))) then
        nfaces_surface = nfaces_surface + 1
      endif
      ! zmax
      if (iglob_is_surface_cross_section(ibool(2,2,NGLLZ,ispec))) then
        nfaces_surface = nfaces_surface + 1
      endif
      ! ymin
      if (iglob_is_surface_cross_section(ibool(2,1,2,ispec))) then
        nfaces_surface = nfaces_surface + 1
      endif
      ! ymax
      if (iglob_is_surface_cross_section(ibool(2,NGLLY,2,ispec))) then
        nfaces_surface = nfaces_surface + 1
      endif
      !xmin
      if (iglob_is_surface_cross_section(ibool(1,2,2,ispec))) then
        nfaces_surface = nfaces_surface + 1
      endif
      !xmax
      if (iglob_is_surface_cross_section(ibool(NGLLX,2,2,ispec))) then
        nfaces_surface = nfaces_surface + 1
      endif
    endif
  enddo

  end subroutine detect_surface_cross_section


!
!-------------------------------------------------------------------------------------------------
!

  subroutine ds_set_cross_section_flags(nspec,ispec_is_surface, &
                                        nglob,iglob_is_surface, &
                                        i,j,k,ispec,ibool, &
                                        valence,countval)

  ! put this into separate subroutine to compile faster, otherwise compilers will try to unroll all do loops

  use constants

  implicit none

  ! global indexing
  integer,intent(in) :: nglob,nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  integer,intent(in) :: i,j,k,ispec

  integer,intent(inout) :: countval
  integer, dimension(nglob),intent(inout) :: valence

  !   surface flags
  logical, dimension(nspec),intent(inout) :: ispec_is_surface
  logical, dimension(nglob),intent(inout) :: iglob_is_surface

  ! local parameters
  integer :: kk,jj,ii
  logical :: has_face

  ! initialize element flag
  has_face = .false.

  ! sets flags for all GLL points on a face and makes sure it's not inside the element
  ! zmin & zmax face
  if ((k == 1 .or. k == NGLLZ) .and. valence(ibool(3,3,k,ispec)) >= 1) then
    has_face = .true.
    do jj = 1, NGLLY
      do ii = 1, NGLLX
        iglob_is_surface(ibool(ii,jj,k,ispec)) = .true.
        ! resets valence to count face only once
        valence(ibool(ii,jj,k,ispec)) = -1
      enddo
    enddo
  endif

  ! ymin & ymax
  if ((j == 1 .or. j == NGLLY) .and. valence(ibool(3,j,3,ispec)) >= 1) then
    has_face = .true.
    do kk = 1, NGLLZ
      do ii = 1, NGLLX
        iglob_is_surface(ibool(ii,j,kk,ispec)) = .true.
        ! resets valence to count face only once
        valence(ibool(ii,j,kk,ispec)) = -1
      enddo
    enddo
  endif

  ! xmin & xmax
  if ((i == 1 .or. i == NGLLX) .and. valence(ibool(i,3,3,ispec)) >= 1) then
    has_face = .true.
    do kk = 1, NGLLZ
      do jj = 1, NGLLY
        iglob_is_surface(ibool(i,jj,kk,ispec)) = .true.
        ! resets valence to count face only once
        valence(ibool(i,jj,kk,ispec)) = -1
      enddo
    enddo
  endif

  ! sets flag for element to indicate that it has a face on surface
  if (has_face) then
    ispec_is_surface(ispec) = .true.
    countval = countval + 1
  endif

  end subroutine ds_set_cross_section_flags

!
!-------------------------------------------------------------------------------------------------
!

  subroutine ds_set_plane_flags(iface,ispec, &
                                nspec,ispec_is_surface, &
                                nglob,iglob_is_surface, &
                                ibool,valence)

  ! put this into separate subroutine to compile faster, otherwise compilers will try to unroll all do loops

  use constants

  implicit none

  integer,intent(in) :: iface,ispec

  ! global indexing
  integer,intent(in) :: nglob,nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  integer, dimension(nglob),intent(inout) :: valence

  !   surface flags
  logical, dimension(nspec),intent(inout) :: ispec_is_surface
  logical, dimension(nglob),intent(inout) :: iglob_is_surface


  ! local parameters
  integer :: jj,ii,i,j,k
  integer,dimension(3,NGLLX,NGLLX) :: face_ijk

  call get_element_face_gll_indices(iface,face_ijk,NGLLX,NGLLX)

  do jj = 1, NGLLY
    do ii = 1, NGLLX
      i = face_ijk(1,ii,jj)
      j = face_ijk(2,ii,jj)
      k = face_ijk(3,ii,jj)
      ! sets iglob flag on face points
      iglob_is_surface(ibool(i,j,k,ispec)) = .true.
      ! sets ispec flag
      ispec_is_surface(ispec) = .true.
      ! resets valence
      valence(ibool(i,j,k,ispec)) = -1
    enddo
  enddo

  end subroutine ds_set_plane_flags

