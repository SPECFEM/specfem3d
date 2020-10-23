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

  subroutine store_coords(xstore,ystore,zstore,xelm,yelm,zelm,ispec,nspec,shape3D)

  use constants, only: NGNOD_EIGHT_CORNERS
  use constants_meshfem3D, only: NGLLX_M,NGLLY_M,NGLLZ_M
  use shared_parameters, only: NGNOD

  implicit none

  integer,intent(in) :: ispec,nspec
  double precision, dimension(NGNOD_EIGHT_CORNERS),intent(in) :: xelm,yelm,zelm
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),intent(inout) :: xstore,ystore,zstore
  double precision, dimension(NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M),intent(in) :: shape3D

  ! local parameters
  double precision, dimension(NGNOD) :: offset_x,offset_y,offset_z
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M) :: x_elem,y_elem,z_elem
  integer :: i,j,k

  ! element nodes
  if (NGLLX_M == 2 .and. NGLLY_M == 2 .and. NGLLZ_M == 2) then
    ! safety check
    if (NGNOD /= 8) stop 'need NGNOD == 8 for store_coords() routine when NGLLX_M == 2'

    ! only corners
    xstore(1,1,1,ispec) = xelm(1)
    ystore(1,1,1,ispec) = yelm(1)
    zstore(1,1,1,ispec) = zelm(1)

    xstore(NGLLX_M,1,1,ispec) = xelm(2)
    ystore(NGLLX_M,1,1,ispec) = yelm(2)
    zstore(NGLLX_M,1,1,ispec) = zelm(2)

    xstore(NGLLX_M,NGLLY_M,1,ispec) = xelm(3)
    ystore(NGLLX_M,NGLLY_M,1,ispec) = yelm(3)
    zstore(NGLLX_M,NGLLY_M,1,ispec) = zelm(3)

    xstore(1,NGLLY_M,1,ispec) = xelm(4)
    ystore(1,NGLLY_M,1,ispec) = yelm(4)
    zstore(1,NGLLY_M,1,ispec) = zelm(4)

    xstore(1,1,NGLLZ_M,ispec) = xelm(5)
    ystore(1,1,NGLLZ_M,ispec) = yelm(5)
    zstore(1,1,NGLLZ_M,ispec) = zelm(5)

    xstore(NGLLX_M,1,NGLLZ_M,ispec) = xelm(6)
    ystore(NGLLX_M,1,NGLLZ_M,ispec) = yelm(6)
    zstore(NGLLX_M,1,NGLLZ_M,ispec) = zelm(6)

    xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) = xelm(7)
    ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) = yelm(7)
    zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) = zelm(7)

    xstore(1,NGLLY_M,NGLLZ_M,ispec) = xelm(8)
    ystore(1,NGLLY_M,NGLLZ_M,ispec) = yelm(8)
    zstore(1,NGLLY_M,NGLLZ_M,ispec) = zelm(8)

    ! all done
    return

  else
    ! computes all GLL node positions
    ! anchor positions at corners
    offset_x(1:NGNOD_EIGHT_CORNERS) = xelm(1:NGNOD_EIGHT_CORNERS)
    offset_y(1:NGNOD_EIGHT_CORNERS) = yelm(1:NGNOD_EIGHT_CORNERS)
    offset_z(1:NGNOD_EIGHT_CORNERS) = zelm(1:NGNOD_EIGHT_CORNERS)

    ! adds missing anchor points
    if (NGNOD == 27) then
      call add_missing_nodes(offset_x,offset_y,offset_z)
    endif

    ! calculates positions
    call calc_coords_mesh(x_elem,y_elem,z_elem,offset_x,offset_y,offset_z,shape3D)

    ! stores all element GLL node positions
    do k = 1,NGLLZ_M
      do j = 1,NGLLY_M
        do i = 1,NGLLX_M
          xstore(i,j,k,ispec) = x_elem(i,j,k)
          ystore(i,j,k,ispec) = y_elem(i,j,k)
          zstore(i,j,k,ispec) = z_elem(i,j,k)
        enddo
      enddo
    enddo
  endif

  end subroutine store_coords

!
!-------------------------------------------------------------------------------------------------
!

  subroutine add_missing_nodes(offset_x,offset_y,offset_z)

! compute the missing nodes of a 27-node element when only the 8 corners have been given

! the topology of the nodes is described in file hex_nodes.f90 as well as in
! UTILS/chunk_notes_scanned/numbering_convention_27_nodes.*

  use constants, only: NGNOD_EIGHT_CORNERS
  use shared_parameters, only: NGNOD

  implicit none

  double precision, dimension(NGNOD) :: offset_x,offset_y,offset_z

! list of corners defining the edges and the faces
  integer, parameter :: NEDGES = 12, NFACES = 6
  integer, dimension(NEDGES,2) :: list_corners_edge
  integer, dimension(NFACES,4) :: list_corners_face

  integer :: iedge,iface,ignod

  ! check if anything to do
  if (NGNOD == 8) return

! list of corners defining the edges
! the edge number is sorted according to the numbering convention defined in file hex_nodes.f90
! as well as in DATA/util/YYYYYYYYYYYYYYYYYYYYYYYYYYY DK DK UGLY YYYYYYYYYYYYYYYYYYY

  list_corners_edge( 1,1) = 1
  list_corners_edge( 1,2) = 2

  list_corners_edge( 2,1) = 2
  list_corners_edge( 2,2) = 3

  list_corners_edge( 3,1) = 3
  list_corners_edge( 3,2) = 4

  list_corners_edge( 4,1) = 4
  list_corners_edge( 4,2) = 1

  list_corners_edge( 5,1) = 1
  list_corners_edge( 5,2) = 5

  list_corners_edge( 6,1) = 2
  list_corners_edge( 6,2) = 6

  list_corners_edge( 7,1) = 3
  list_corners_edge( 7,2) = 7

  list_corners_edge( 8,1) = 4
  list_corners_edge( 8,2) = 8

  list_corners_edge( 9,1) = 5
  list_corners_edge( 9,2) = 6

  list_corners_edge(10,1) = 6
  list_corners_edge(10,2) = 7

  list_corners_edge(11,1) = 7
  list_corners_edge(11,2) = 8

  list_corners_edge(12,1) = 8
  list_corners_edge(12,2) = 5

! list of corners defining the faces
! the face number is sorted according to the numbering convention defined in file hex_nodes.f90
! as well as in DATA/util/YYYYYYYYYYYYYYYYYYYYYYYYYYY DK DK UGLY YYYYYYYYYYYYYYYYYYY

  list_corners_face(1,1) = 1
  list_corners_face(1,2) = 2
  list_corners_face(1,3) = 3
  list_corners_face(1,4) = 4

  list_corners_face(2,1) = 1
  list_corners_face(2,2) = 2
  list_corners_face(2,3) = 6
  list_corners_face(2,4) = 5

  list_corners_face(3,1) = 2
  list_corners_face(3,2) = 3
  list_corners_face(3,3) = 7
  list_corners_face(3,4) = 6

  list_corners_face(4,1) = 4
  list_corners_face(4,2) = 3
  list_corners_face(4,3) = 7
  list_corners_face(4,4) = 8

  list_corners_face(5,1) = 1
  list_corners_face(5,2) = 4
  list_corners_face(5,3) = 8
  list_corners_face(5,4) = 5

  list_corners_face(6,1) = 5
  list_corners_face(6,2) = 6
  list_corners_face(6,3) = 7
  list_corners_face(6,4) = 8

! midside nodes (nodes located in the middle of an edge)
  do iedge = 1,NEDGES

! node numbers for edge centers start at 9
    ignod = (iedge - 1) + 9

    offset_x(ignod) = (offset_x(list_corners_edge(iedge,1)) + offset_x(list_corners_edge(iedge,2))) / 2.d0

    offset_y(ignod) = (offset_y(list_corners_edge(iedge,1)) + offset_y(list_corners_edge(iedge,2))) / 2.d0

    offset_z(ignod) = (offset_z(list_corners_edge(iedge,1)) + offset_z(list_corners_edge(iedge,2))) / 2.d0

  enddo

! side center nodes (nodes located in the middle of a face)
  do iface = 1,NFACES

! node numbers for face centers start at 21
    ignod = (iface - 1) + 21

    offset_x(ignod) = (offset_x(list_corners_face(iface,1)) + &
                       offset_x(list_corners_face(iface,2)) + &
                       offset_x(list_corners_face(iface,3)) + &
                       offset_x(list_corners_face(iface,4))) / 4.d0

    offset_y(ignod) = (offset_y(list_corners_face(iface,1)) + &
                       offset_y(list_corners_face(iface,2)) + &
                       offset_y(list_corners_face(iface,3)) + &
                       offset_y(list_corners_face(iface,4))) / 4.d0

    offset_z(ignod) = (offset_z(list_corners_face(iface,1)) + &
                       offset_z(list_corners_face(iface,2)) + &
                       offset_z(list_corners_face(iface,3)) + &
                       offset_z(list_corners_face(iface,4))) / 4.d0

  enddo

! center node (barycenter of the eight corners)
  offset_x(27) = sum(offset_x(1:NGNOD_EIGHT_CORNERS)) / dble(NGNOD_EIGHT_CORNERS)
  offset_y(27) = sum(offset_y(1:NGNOD_EIGHT_CORNERS)) / dble(NGNOD_EIGHT_CORNERS)
  offset_z(27) = sum(offset_z(1:NGNOD_EIGHT_CORNERS)) / dble(NGNOD_EIGHT_CORNERS)

  end subroutine add_missing_nodes

!
!-------------------------------------------------------------------------------------------------
!

  subroutine calc_coords_mesh(x_elem,y_elem,z_elem,offset_x,offset_y,offset_z,shape3D)

  use constants_meshfem3D, only: NGLLX_M,NGLLY_M,NGLLZ_M
  use shared_parameters, only: NGNOD

  implicit none

  double precision, dimension(NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M),intent(in) :: shape3D
  double precision, dimension(NGNOD),intent(in) :: offset_x,offset_y,offset_z

  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M),intent(out) :: x_elem,y_elem,z_elem

  ! local parameters
  integer :: i,j,k,ia
  double precision :: xmesh,ymesh,zmesh

  do k = 1,NGLLZ_M
    do j = 1,NGLLY_M
      do i = 1,NGLLX_M
        xmesh = 0.d0
        ymesh = 0.d0
        zmesh = 0.d0

        do ia = 1,NGNOD
          xmesh = xmesh + shape3D(ia,i,j,k)*offset_x(ia)
          ymesh = ymesh + shape3D(ia,i,j,k)*offset_y(ia)
          zmesh = zmesh + shape3D(ia,i,j,k)*offset_z(ia)
        enddo

        x_elem(i,j,k) = xmesh
        y_elem(i,j,k) = ymesh
        z_elem(i,j,k) = zmesh
      enddo
    enddo
  enddo

  end subroutine calc_coords_mesh
