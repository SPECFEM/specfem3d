
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

  program create_OpenDX_file_to_view_the_mesh

! Dimitri Komatitsch, CNRS Marseille, France, September 2018

  implicit none

  integer :: nspec,nglob
  integer :: ispec,iglob
  integer :: iglob_read,ispec_loop,iformat,NGNOD,ia

  double precision, dimension(:), allocatable :: x,y,z

  double precision :: xread,yread,zread,xmin,xmax,ymin,ymax,zmin,zmax

  integer, dimension(:,:), allocatable :: ibool
  integer, dimension(:), allocatable :: imaterial

! read the mesh files in ASCII format (that is the standard case)
  iformat = 1

! the mesh contains HEX8 elements
  NGNOD = 8

! open SPECFEM3D_Cartesian mesh file to read the points
  if (iformat == 1) then
    open(unit=23,file='DATA_3D/nodes_coords_file',status='old',action='read')
    read(23,*) nglob
  else
    open(unit=23,file='DATA_3D/nodes_coords_file.bin',form='unformatted',status='old',action='read')
    read(23) nglob
  endif
  allocate(x(nglob))
  allocate(y(nglob))
  allocate(z(nglob))
  if (iformat == 1) then
    do iglob = 1,nglob
      read(23,*) iglob_read,xread,yread,zread
      x(iglob_read) = xread
      y(iglob_read) = yread
      z(iglob_read) = zread
    enddo
  else
    read(23) x
    read(23) y
    read(23) z
  endif
  close(23)

! compute the min and max values of each coordinate
   xmin = minval(x)
   xmax = maxval(x)

   ymin = minval(y)
   ymax = maxval(y)

   zmin = minval(z)
   zmax = maxval(z)

  print *,'Xmin and Xmax of the mesh read = ',xmin,xmax
  print *,'Ymin and Ymax of the mesh read = ',ymin,ymax
  print *,'Zmin and Zmax of the mesh read = ',zmin,zmax
  print *

! ************* read mesh elements *************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  if (iformat == 1) then
    open(unit=23,file='DATA_3D/mesh_file',status='old',action='read')
    read(23,*) nspec
  else
    open(unit=23,file='DATA_3D/mesh_file.bin',form='unformatted',status='old',action='read')
    read(23) nspec
  endif

  allocate(ibool(NGNOD,nspec))
  allocate(imaterial(nspec))

! loop on the whole mesh
  if (iformat == 1) then
    do ispec_loop = 1,nspec
      read(23,*) ispec,(ibool(ia,ispec), ia = 1,NGNOD)
    enddo
  else
    read(23) ibool
  endif

  close(23)

  print *,'Total number of elements in the mesh read = ',nspec
  print *

! read the material property numbers
  open(unit=23,file='DATA_3D/materials_file',status='old',action='read')
  do ispec_loop = 1,nspec
    read(23,*) ispec,imaterial(ispec)
  enddo
  close(23)

!! DK DK create an OpenDX file showing all the elements of the mesh

    open(unit=11,file='DX_fullmesh.dx',status='unknown')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'

    do iglob = 1,nglob
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") x(iglob),y(iglob),z(iglob)
    enddo

! ************* generate elements ******************

  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',nspec,' data follows'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
! in the case of OpenDX, node numbers start at zero
      write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
            ibool(4,ispec)-1, ibool(1,ispec)-1, ibool(8,ispec)-1, ibool(5,ispec)-1, &
            ibool(3,ispec)-1, ibool(2,ispec)-1, ibool(7,ispec)-1, ibool(6,ispec)-1
  enddo

! ************* generate element data values ******************

! output AVS or DX header for data
    write(11,*) 'attribute "element type" string "cubes"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
        write(11,*) imaterial(ispec)
  enddo

! define OpenDX field
    write(11,*) 'attribute "dep" string "connections"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'

  close(11)

  end program create_OpenDX_file_to_view_the_mesh

