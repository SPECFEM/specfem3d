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


  subroutine create_visual_files(CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES,&
                                nspec,nglob,prname,nodes_coords,ibool,true_material_num)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

! Mesh files for visualization
  logical CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES

! number of spectral elements in each block
  integer nspec

! number of vertices in each block
  integer nglob

! name of the database files
  character(len=256) prname

! arrays with the mesh
  integer true_material_num(nspec)
  integer ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision :: nodes_coords(nglob,3)

!  ---------------
  integer i,ipoin,ispec

  if(CREATE_ABAQUS_FILES) then

     open(unit=64,file=prname(1:len_trim(prname))//'mesh.INP',status='unknown',action='write',form='formatted')
     write(64,'(a8)') '*HEADING'
!     write(64,'(a52)') 'cubit(mesh): 04/17/2009: 18:11:24'
     write(64,'(a27)') 'SPECFEM3D meshfem3D(mesh): '
     write(64,'(a5)') '*NODE'
     do i=1,nglob
!        write(64,*) i,',',nodes_coords(i,1),',',nodes_coords(i,2),',',nodes_coords(i,3)
         write(64,'(i10,3(a,e15.6))') i,',',nodes_coords(i,1),',',nodes_coords(i,2),',',nodes_coords(i,3)
     enddo
     write(64,'(a31)') '*ELEMENT, TYPE=C3D8R, ELSET=EB1'
     do ispec=1,nspec
        write(64,'(i10,8(a,i10))') ispec,',',ibool(1,1,2,ispec),',',ibool(1,1,1,ispec),',',ibool(1,2,1,ispec), &
             ',',ibool(1,2,2,ispec),',',ibool(2,1,2,ispec),',',ibool(2,1,1,ispec),',',ibool(2,2,1,ispec),',', &
             ibool(2,2,2,ispec)
     enddo
     close(64)

  endif


  if(CREATE_DX_FILES) then

     open(unit=66,file=prname(1:len_trim(prname))//'mesh.dx',status='unknown')

     ! write OpenDX header
     write(66,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'

     do ipoin = 1,nglob
        write(66,*) sngl(nodes_coords(ipoin,1)),sngl(nodes_coords(ipoin,2)),sngl(nodes_coords(ipoin,3))
     enddo

     ! ************* generate elements ******************

     write(66,*) 'object 2 class array type int rank 1 shape ',8,' items ',nspec,' data follows'

     do ispec=1,nspec

        ! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
        ! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
        ! in the case of OpenDX, node numbers start at zero
        write(66,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
             ibool(1,1,2,ispec)-1,ibool(2,1,2,ispec)-1,ibool(1,2,2,ispec)-1,ibool(2,2,2,ispec)-1,&
             ibool(1,1,1,ispec)-1,ibool(2,1,1,ispec)-1,ibool(1,2,1,ispec)-1,ibool(2,2,1,ispec)-1
     enddo



     ! ************* generate element data values ******************

     ! output OpenDX header for data
     write(66,*) 'attribute "element type" string "cubes"'

     write(66,*) 'attribute "ref" string "positions"'
     write(66,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

     ! loop on all the elements
     do ispec = 1,nspec
        write(66,*)  true_material_num(ispec)
     enddo

     write(66,*) 'attribute "dep" string "connections"'
     write(66,*) 'object "irregular positions irregular connections" class field'
     write(66,*) 'component "positions" value 1'
     write(66,*) 'component "connections" value 2'
     write(66,*) 'component "data" value 3'
     write(66,*) 'end'

     close(66)

  endif

  if( CREATE_VTK_FILES ) then
    ! vtk file output
    open(66,file=prname(1:len_trim(prname))//'mesh.vtk',status='unknown')
    write(66,'(a)') '# vtk DataFile Version 3.1'
    write(66,'(a)') 'material model VTK file'
    write(66,'(a)') 'ASCII'
    write(66,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(66, '(a,i12,a)') 'POINTS ', nglob, ' float'
    do ipoin = 1,nglob
      write(66,*) sngl(nodes_coords(ipoin,1)),sngl(nodes_coords(ipoin,2)),sngl(nodes_coords(ipoin,3))
    enddo
    write(66,*) ""

    ! note: indices for vtk start at 0
    write(66,'(a,i12,i12)') "CELLS ",nspec,nspec*9
    do ispec=1,nspec
      write(66,'(9i12)') 8, &
            ibool(1,1,1,ispec)-1,ibool(2,1,1,ispec)-1,ibool(2,2,1,ispec)-1,ibool(1,2,1,ispec)-1,&
            ibool(1,1,2,ispec)-1,ibool(2,1,2,ispec)-1,ibool(2,2,2,ispec)-1,ibool(1,2,2,ispec)-1
    enddo
    write(66,*) ""

    ! type: hexahedrons
    write(66,'(a,i12)') "CELL_TYPES ",nspec
    write(66,'(6i12)') (12,ispec=1,nspec)
    write(66,*) ""

    write(66,'(a,i12)') "CELL_DATA ",nspec
    write(66,'(a)') "SCALARS elem_val float"
    write(66,'(a)') "LOOKUP_TABLE default"
    do ispec = 1,nspec
      write(66,*) true_material_num(ispec)
    enddo
    write(66,*) ""
    close(66)

  endif

  call sync_all()

  end subroutine create_visual_files

