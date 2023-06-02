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

! read an external mesh file (list of points and list of elements) composed of tetrahedra
! and create a mesh of hexahedra by cutting each tetrahedron into four hexahedra
! using the middle of each edge, each face and the barycenter.
! For a picture of what this gives, see e.g. http://www.tetgen.org/figs/Delaunay-Voronoi-3D.gif
! (or a copy of it in utils/Cubit_or_Gmsh/logo_of_TetGen_showing_Delaunay_Voronoi_3D_how_to_cut_a_tetra_into_four_hexas.gif)

! Dimitri Komatitsch, CNRS, Marseille, France, June 2015.

  program convert_tetra_mesh_to_hexa_mesh

  implicit none

!
! work in single or in double precision (4 or 8 bytes)
!
  integer, parameter :: CUSTOM_REAL = 4 ! 8

! read list of elements stored in new Gmsh 2.9.3 format or in old Gmsh 2.4.2 format (the old one has one extra dummy value)
  logical, parameter :: USE_OLD_GMSH_MESH_FORMAT = .true. ! .false.

  real(kind=CUSTOM_REAL), parameter :: ONE_THIRD = 1._CUSTOM_REAL / 3._CUSTOM_REAL

  integer :: nglob,ntet,nelem_in_file

  integer :: i,k,ihexa,iread,itype,idummy1,idummy2,idummy3,ivolume

  real(kind=CUSTOM_REAL) :: xread,yread,zread

  real(kind=CUSTOM_REAL), allocatable, dimension(:) :: x,y,z

  integer, dimension(15) :: number_of_points_per_element_type
  integer, dimension(0:3) :: inode_read

! coordinates of nodes in the middle of the tetrahedron edges
  real(kind=CUSTOM_REAL) :: mid_0_1_x,mid_0_2_x,mid_0_3_x,mid_1_2_x,mid_1_3_x,mid_2_3_x
  real(kind=CUSTOM_REAL) :: mid_0_1_y,mid_0_2_y,mid_0_3_y,mid_1_2_y,mid_1_3_y,mid_2_3_y
  real(kind=CUSTOM_REAL) :: mid_0_1_z,mid_0_2_z,mid_0_3_z,mid_1_2_z,mid_1_3_z,mid_2_3_z

! coordinates of nodes in the middle of the tetrahedron faces
  real(kind=CUSTOM_REAL) :: mid_face_0_1_2_x,mid_face_0_1_3_x,mid_face_0_2_3_x,mid_face_1_2_3_x
  real(kind=CUSTOM_REAL) :: mid_face_0_1_2_y,mid_face_0_1_3_y,mid_face_0_2_3_y,mid_face_1_2_3_y
  real(kind=CUSTOM_REAL) :: mid_face_0_1_2_z,mid_face_0_1_3_z,mid_face_0_2_3_z,mid_face_1_2_3_z

! coordinates of the barycenter of the tetrahedron
  real(kind=CUSTOM_REAL) :: barycenter_x,barycenter_y,barycenter_z

  integer :: iglob_0,iglob_1,iglob_2,iglob_3
  integer :: iglob_mid_0_1,iglob_mid_0_2,iglob_mid_0_3,iglob_mid_1_2,iglob_mid_1_3,iglob_mid_2_3,iglob_barycenter
  integer :: iglob_mid_face_0_1_2,iglob_mid_face_0_1_3,iglob_mid_face_0_2_3,iglob_mid_face_1_2_3

  number_of_points_per_element_type(:) = 0

! for Gmsh element types (from http://gmsh.info/doc/texinfo/gmsh.html#Low-order-elements )
  number_of_points_per_element_type(1)  = 2  ! 2-node line
  number_of_points_per_element_type(2)  = 3  ! 3-node triangle
  number_of_points_per_element_type(4)  = 4  ! 4-node tetrahedron
  number_of_points_per_element_type(15) = 1  ! 1-node point

! point numbering convention for tetrahedra in Gmsh (from http://gmsh.info/doc/texinfo/gmsh.html#Low-order-elements )
!
! Tetrahedron4:                               Tetrahedron10:
!                    v
!                  .
!                ,/
!               /
!            2                                     2
!          ,/|`\                                 ,/|`\
!        ,/  |  `\                             ,/  |  `\
!      ,/    '.   `\                         ,6    '.   `5
!    ,/       |     `\                     ,/       8     `\
!  ,/         |       `\                 ,/         |       `\
! 0-----------'.--------1 --> u         0--------4--'.--------1
!  `\.         |      ,/                 `\.         |      ,/
!     `\.      |    ,/                      `\.      |    ,9
!        `\.   '. ,/                           `7.   '. ,/
!           `\. |/                                `\. |/
!              `3                                    `3
!                 `\.
!                    ` w
!

! --------- read mesh points ---------

open(unit=10,file='points.txt',status='old',action='read')

read(10,*) nglob

print *,'reading ',nglob,' mesh points from the Gmsh database'

allocate(x(nglob))
allocate(y(nglob))
allocate(z(nglob))

do i = 1,nglob
  read(10,*) iread,xread,yread,zread
  if (iread < 1 .or. iread > nglob) stop 'incorrect point read'
  x(iread) = xread
  y(iread) = yread
  z(iread) = zread
enddo

close(10)

print *
print *,'x min max read = ',minval(x),maxval(x)
print *,'y min max read = ',minval(y),maxval(y)
print *,'z min max read = ',minval(z),maxval(z)
print *

! --------- read mesh points ---------

open(unit=10,file='elements.txt',status='old',action='read')

read(10,*) nelem_in_file

print *,'reading ',nelem_in_file,' mesh elements of any geometrical kind from the Gmsh database'

ntet = 0

do i = 1,nelem_in_file
  inode_read(:) = 0
! read list of elements stored in new Gmsh 2.9.3 format or in old Gmsh 2.4.2 format (the old one has one extra dummy value)
  if (USE_OLD_GMSH_MESH_FORMAT) then
    read(10,*) iread,itype,idummy1,idummy2,ivolume,idummy3,(inode_read(k), k=0,number_of_points_per_element_type(itype)-1)
  else
    read(10,*) iread,itype,idummy1,idummy2,ivolume,(inode_read(k), k=0,number_of_points_per_element_type(itype)-1)
  endif
  if (number_of_points_per_element_type(itype) <= 0) stop 'incorrect element type read'
  if (iread < 1 .or. iread > nelem_in_file) stop 'incorrect element read'
  if (itype == 4) then
    ntet = ntet + 1
  endif
enddo

close(10)

print *
print *,'number of tetrahedra read = ',ntet

! writing the database for the hexahedral mesh

print *
print *,'writing the database for the new mesh consisting of hexahedra...'

! create the subdirectory to store the mesh if it does not exist
call system('mkdir -p MESH')

open(unit=9,file='points.txt',status='old',action='read')

open(unit=10,file='elements.txt',status='old',action='read')

! file with hexahedra points (four hexahedra created out of each tetrahedron)
open(unit=14,file='MESH/nodes_coords_file',status='unknown',action='write')

! file with hexahedra (four hexahedra created out of each tetrahedron)
open(unit=15,file='MESH/mesh_file',status='unknown',action='write')

read(10,*) nelem_in_file

! we need to copy the existing list of points and then add 11 new points in each tetrahedron
write(14,*) nglob + 11*ntet

! out of each tetrahedron we create four hexahedra
write(15,*) 4*ntet

print *
print *,'the mesh of hexahedra will have ',nglob + 11*ntet,' points'
print *,'and ',4*ntet,' elements'
print *

read(9,*) nglob
do i = 1,nglob
  read(9,*) iread,xread,yread,zread
  write(14,*) iread,xread,yread,zread
enddo

ntet = 0
ihexa = 0

do i = 1,nelem_in_file

  inode_read(:) = 0

! read list of elements stored in new Gmsh 2.9.3 format or in old Gmsh 2.4.2 format (the old one has one extra dummy value)
  if (USE_OLD_GMSH_MESH_FORMAT) then
    read(10,*) iread,itype,idummy1,idummy2,ivolume,idummy3,(inode_read(k), k=0,number_of_points_per_element_type(itype)-1)
  else
    read(10,*) iread,itype,idummy1,idummy2,ivolume,(inode_read(k), k=0,number_of_points_per_element_type(itype)-1)
  endif
  if (number_of_points_per_element_type(itype) <= 0) stop 'incorrect element type read'
  if (iread < 1 .or. iread > nelem_in_file) stop 'incorrect element read'

  ! processing only the elements that are tetrahedra
  if (itype == 4) then
    ntet = ntet + 1

! now let us cut each tetrahedron into four hexahedra using the middle of each edge, each face and the barycenter.
! For a picture of what this gives, see e.g. http://www.tetgen.org/figs/Delaunay-Voronoi-3D.gif
! (or a copy of it in utils/Cubit_or_Gmsh/logo_of_TetGen_showing_Delaunay_Voronoi_3D_how_to_cut_a_tetra_into_four_hexas.gif)

  ! new points located in the middle of the edges
  mid_0_1_x = 0.5_CUSTOM_REAL * (x(inode_read(0)) + x(inode_read(1)))
  mid_0_2_x = 0.5_CUSTOM_REAL * (x(inode_read(0)) + x(inode_read(2)))
  mid_0_3_x = 0.5_CUSTOM_REAL * (x(inode_read(0)) + x(inode_read(3)))
  mid_1_2_x = 0.5_CUSTOM_REAL * (x(inode_read(1)) + x(inode_read(2)))
  mid_1_3_x = 0.5_CUSTOM_REAL * (x(inode_read(1)) + x(inode_read(3)))
  mid_2_3_x = 0.5_CUSTOM_REAL * (x(inode_read(2)) + x(inode_read(3)))

  mid_0_1_y = 0.5_CUSTOM_REAL * (y(inode_read(0)) + y(inode_read(1)))
  mid_0_2_y = 0.5_CUSTOM_REAL * (y(inode_read(0)) + y(inode_read(2)))
  mid_0_3_y = 0.5_CUSTOM_REAL * (y(inode_read(0)) + y(inode_read(3)))
  mid_1_2_y = 0.5_CUSTOM_REAL * (y(inode_read(1)) + y(inode_read(2)))
  mid_1_3_y = 0.5_CUSTOM_REAL * (y(inode_read(1)) + y(inode_read(3)))
  mid_2_3_y = 0.5_CUSTOM_REAL * (y(inode_read(2)) + y(inode_read(3)))

  mid_0_1_z = 0.5_CUSTOM_REAL * (z(inode_read(0)) + z(inode_read(1)))
  mid_0_2_z = 0.5_CUSTOM_REAL * (z(inode_read(0)) + z(inode_read(2)))
  mid_0_3_z = 0.5_CUSTOM_REAL * (z(inode_read(0)) + z(inode_read(3)))
  mid_1_2_z = 0.5_CUSTOM_REAL * (z(inode_read(1)) + z(inode_read(2)))
  mid_1_3_z = 0.5_CUSTOM_REAL * (z(inode_read(1)) + z(inode_read(3)))
  mid_2_3_z = 0.5_CUSTOM_REAL * (z(inode_read(2)) + z(inode_read(3)))

  ! new points located in the middle of the faces
  mid_face_0_1_2_x = ONE_THIRD * (x(inode_read(0)) + x(inode_read(1)) + x(inode_read(2)))
  mid_face_0_1_3_x = ONE_THIRD * (x(inode_read(0)) + x(inode_read(1)) + x(inode_read(3)))
  mid_face_0_2_3_x = ONE_THIRD * (x(inode_read(0)) + x(inode_read(2)) + x(inode_read(3)))
  mid_face_1_2_3_x = ONE_THIRD * (x(inode_read(1)) + x(inode_read(2)) + x(inode_read(3)))

  mid_face_0_1_2_y = ONE_THIRD * (y(inode_read(0)) + y(inode_read(1)) + y(inode_read(2)))
  mid_face_0_1_3_y = ONE_THIRD * (y(inode_read(0)) + y(inode_read(1)) + y(inode_read(3)))
  mid_face_0_2_3_y = ONE_THIRD * (y(inode_read(0)) + y(inode_read(2)) + y(inode_read(3)))
  mid_face_1_2_3_y = ONE_THIRD * (y(inode_read(1)) + y(inode_read(2)) + y(inode_read(3)))

  mid_face_0_1_2_z = ONE_THIRD * (z(inode_read(0)) + z(inode_read(1)) + z(inode_read(2)))
  mid_face_0_1_3_z = ONE_THIRD * (z(inode_read(0)) + z(inode_read(1)) + z(inode_read(3)))
  mid_face_0_2_3_z = ONE_THIRD * (z(inode_read(0)) + z(inode_read(2)) + z(inode_read(3)))
  mid_face_1_2_3_z = ONE_THIRD * (z(inode_read(1)) + z(inode_read(2)) + z(inode_read(3)))

  ! new point located in the middle of the element (i.e. at its barycenter)
  barycenter_x = 0.25_CUSTOM_REAL * (x(inode_read(0)) + x(inode_read(1)) + x(inode_read(2)) + x(inode_read(3)))
  barycenter_y = 0.25_CUSTOM_REAL * (y(inode_read(0)) + y(inode_read(1)) + y(inode_read(2)) + y(inode_read(3)))
  barycenter_z = 0.25_CUSTOM_REAL * (z(inode_read(0)) + z(inode_read(1)) + z(inode_read(2)) + z(inode_read(3)))

  ! write these new points in the database of points

  iglob_mid_0_1 = nglob + 1
  iglob_mid_0_2 = nglob + 2
  iglob_mid_0_3 = nglob + 3
  iglob_mid_1_2 = nglob + 4
  iglob_mid_1_3 = nglob + 5
  iglob_mid_2_3 = nglob + 6

  iglob_mid_face_0_1_2 = nglob + 7
  iglob_mid_face_0_1_3 = nglob + 8
  iglob_mid_face_0_2_3 = nglob + 9
  iglob_mid_face_1_2_3 = nglob + 10

  iglob_barycenter = nglob + 11

  write(14,*) iglob_mid_0_1,mid_0_1_x,mid_0_1_y,mid_0_1_z
  write(14,*) iglob_mid_0_2,mid_0_2_x,mid_0_2_y,mid_0_2_z
  write(14,*) iglob_mid_0_3,mid_0_3_x,mid_0_3_y,mid_0_3_z
  write(14,*) iglob_mid_1_2,mid_1_2_x,mid_1_2_y,mid_1_2_z
  write(14,*) iglob_mid_1_3,mid_1_3_x,mid_1_3_y,mid_1_3_z
  write(14,*) iglob_mid_2_3,mid_2_3_x,mid_2_3_y,mid_2_3_z

  write(14,*) iglob_mid_face_0_1_2,mid_face_0_1_2_x,mid_face_0_1_2_y,mid_face_0_1_2_z
  write(14,*) iglob_mid_face_0_1_3,mid_face_0_1_3_x,mid_face_0_1_3_y,mid_face_0_1_3_z
  write(14,*) iglob_mid_face_0_2_3,mid_face_0_2_3_x,mid_face_0_2_3_y,mid_face_0_2_3_z
  write(14,*) iglob_mid_face_1_2_3,mid_face_1_2_3_x,mid_face_1_2_3_y,mid_face_1_2_3_z

  write(14,*) iglob_barycenter,barycenter_x,barycenter_y,barycenter_z

  nglob = nglob + 11 ! because we have added eleven points

  ! define some useful aliases for the four existing points
  iglob_0 = inode_read(0)
  iglob_1 = inode_read(1)
  iglob_2 = inode_read(2)
  iglob_3 = inode_read(3)

  ! write the four new elements in the database of elements
  ihexa = ihexa + 1
  write(15,"(9i12)") ihexa,iglob_1,iglob_mid_1_2,iglob_mid_face_0_1_2,iglob_mid_0_1, &
                    iglob_mid_1_3,iglob_mid_face_1_2_3,iglob_barycenter,iglob_mid_face_0_1_3

  ihexa = ihexa + 1
  write(15,"(9i12)") ihexa,iglob_mid_1_2,iglob_2,iglob_mid_0_2,iglob_mid_face_0_1_2, &
                    iglob_mid_face_1_2_3,iglob_mid_2_3,iglob_mid_face_0_2_3,iglob_barycenter

  ihexa = ihexa + 1
  write(15,"(9i12)") ihexa,iglob_0,iglob_mid_0_1,iglob_mid_face_0_1_2,iglob_mid_0_2, &
                    iglob_mid_0_3,iglob_mid_face_0_1_3,iglob_barycenter,iglob_mid_face_0_2_3

  ihexa = ihexa + 1
  write(15,"(9i12)") ihexa,iglob_mid_face_0_1_3,iglob_mid_1_3,iglob_mid_face_1_2_3,iglob_barycenter, &
                    iglob_mid_0_3,iglob_3,iglob_mid_2_3,iglob_mid_face_0_2_3

  endif

enddo

  close(9)
  close(10)
  close(14)
  close(15)

  end program convert_tetra_mesh_to_hexa_mesh

