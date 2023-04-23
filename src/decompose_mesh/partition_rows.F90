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

! ROWS partitioning:
!
!  If you don't want to use an external graph partitioner like SCOTCH, METIS, PATOH ..
!  then you can use this ROWS_PART partitioning as a simple partitioning scheme.
!
!  This partioning schemes make slices along the x-direction.
!
!  These is no need to add compiler flags like -DUSE_SCOTCH, -DUSE_METIS, .. and libraries in Makefile in the SPECFEM3D/ directory.
!
!  To select ROWS_PART instead of SCOTCH, in Par_file choose:
!    PARTITIONING_TYPE = 4


  subroutine partition_rows()

! uses a simple row partioning scheme, i.e., only partitions along x-direction.
! the resulting partitions will look like bread slices across the model.
!
! note: this scheme is mostly meant for testing and comparisons
!       (as an example of an inefficient partitioning scheme),
!       thus better use SCOTCH/METIS for "serious" simulations.

  use decompose_mesh_par, only: nparts,nodes_coords,nnodes,nspec,NGNOD,elmnts,part

  implicit none

  ! local parameters
  integer :: ispec,inode
  double precision :: x_max, x_min
  double precision :: element_center,range,slice_width,dloc
  integer :: sliceid,id

  ! nodes_coords(:,inode) = (x, y, z) of node inode
  ! elmnts(:,ielem) = (1,2,3,4,5,6,7,8) corresponding to nodes_coords (if NGNOD==8)

  ! checks if partitioning needed
  if (nparts == 1) return

  ! ROWS partitioning
  print *,'ROWS partitioning'

  ! x-direction
  x_min = minval(nodes_coords(1,:))
  x_max = maxval(nodes_coords(1,:))
  range = x_max - x_min
  slice_width = range/nparts

  ! user output
  print *,'  x-direction: min = ',sngl(x_min),' max = ',sngl(x_max)
  print *,'  x-direction: range = ',sngl(range)
  print *,'  x-direction: estimated slice width = ',sngl(slice_width)

  ! check (for non-zero division)
  if (range < 1.d-16) stop 'Error zero range in x-direction'

  ! assignes elements to partitions according to location
  do ispec = 1,nspec
    ! debug: element size
    !element_center = 1.d24
    !do inode = 2,NGNOD
    !  ! coordinates
    !  id = elmnts(inode,ispec) + 1
    !  ! distances between points
    !  dloc = (nodes_coords(1,elmnts(1,ispec) + 1) - nodes_coords(1,id))**2
    !  if (dloc > 1.d-9 .and. dloc < element_center) element_center = dloc
    !enddo
    !dloc = sqrt(element_center)

    ! center location
    element_center = 0.d0
    do inode = 1,NGNOD
      ! node index (note that elmnts array starts indexing from 0 here)
      id = elmnts(inode,ispec) + 1
      if (id < 1 .or. id > nnodes) stop 'Error node index exceeds bounds in partition_rows'
      ! debug
      !if (ispec == 200) print *,'ispec: ',ispec,'nodes_coords',id,nodes_coords(1,id),element_center,sngl(dloc)
      element_center = element_center + nodes_coords(1,id)
    enddo
    element_center = element_center/NGNOD

    ! location with respect to x-range
    dloc = (element_center - x_min)/range

    ! slice id for this element
    sliceid = floor(nparts * dloc)

    ! checks limit
    if (sliceid < 0) sliceid = 0
    if (sliceid > nparts-1) sliceid = nparts-1

    ! assignes partition number for this element
    part(ispec) = sliceid

    ! debug
    !if (ispec == 200) print *,'ispec: ',ispec,'element_center = ',element_center,'dloc = ',dloc,'slice id = ',sliceid
  enddo

  end subroutine partition_rows
