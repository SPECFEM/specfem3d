!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! compute the missing nodes of a 27-node element when only the 8 corners have been given

! the topology of the nodes is described in file hex_nodes.f90 as well as in
! UTILS/chunk_notes_scanned/numbering_convention_27_nodes.*

  subroutine add_missing_nodes(offset_x,offset_y,offset_z)

  implicit none

  include "constants.h"

  double precision, dimension(NGNOD) :: offset_x,offset_y,offset_z

! list of corners defining the edges and the faces
  integer, parameter :: NEDGES = 12, NFACES = 6
  integer, dimension(NEDGES,2) :: list_corners_edge
  integer, dimension(NFACES,4) :: list_corners_face

  integer :: iedge,iface,ignod

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

