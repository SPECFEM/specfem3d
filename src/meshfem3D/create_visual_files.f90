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


  subroutine create_visual_files(CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
                                 nspec,nglob,prname,nodes_coords,ibool,ispec_material_id)

  use constants
  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  implicit none

! Mesh files for visualization
  logical,intent(in) :: CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES

! number of spectral elements in each block
  integer,intent(in) :: nspec

! number of vertices in each block
  integer,intent(in) :: nglob

! name of the database files
  character(len=MAX_STRING_LEN),intent(in) :: prname

! arrays with the mesh
  integer,dimension(nspec),intent(in) :: ispec_material_id
  integer,dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),intent(in) :: ibool
  double precision,dimension(nglob,NDIM),intent(in) :: nodes_coords

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer :: i,ipoin,ispec,ier

  if (CREATE_ABAQUS_FILES) then

    filename = prname(1:len_trim(prname))//'mesh.INP'
    open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='formatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error opening file for CREATE_ABAQUS_FILES'
    endif

    write(IOUT,'(a8)') '*HEADING'
    write(IOUT,'(a27)') 'SPECFEM3D meshfem3D(mesh): '
    write(IOUT,'(a5)') '*NODE'
    do i=1,nglob
      write(IOUT,'(i10,3(a,e15.6))') i,',',nodes_coords(i,1),',',nodes_coords(i,2),',',nodes_coords(i,3)
    enddo
    write(IOUT,'(a31)') '*ELEMENT, TYPE=C3D8R, ELSET=EB1'
    do ispec=1,nspec
      write(IOUT,'(i10,8(a,i10))') ispec,',',ibool(1,1,NGLLZ_M,ispec),',',ibool(1,1,1,ispec), &
                                   ',',ibool(1,NGLLY_M,1,ispec), ',',ibool(1,NGLLY_M,NGLLZ_M,ispec), &
                                   ',',ibool(NGLLX_M,1,NGLLZ_M,ispec),',',ibool(NGLLX_M,1,1,ispec), &
                                   ',',ibool(NGLLX_M,NGLLY_M,1,ispec),',',ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    enddo
    close(IOUT)

  endif


  if (CREATE_DX_FILES) then

    filename = prname(1:len_trim(prname))//'mesh.dx'
    open(unit=IOUT,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error opening file for CREATE_DX_FILES'
    endif

    ! write OpenDX header
    write(IOUT,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'

    do ipoin = 1,nglob
      write(IOUT,*) sngl(nodes_coords(ipoin,1)),sngl(nodes_coords(ipoin,2)),sngl(nodes_coords(ipoin,3))
    enddo

    ! ************* generate elements ******************

    write(IOUT,*) 'object 2 class array type int rank 1 shape ',8,' items ',nspec,' data follows'

    do ispec=1,nspec

      ! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
      ! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
      ! in the case of OpenDX, node numbers start at zero
      write(IOUT,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
                 ibool(1,1,NGLLZ_M,ispec)-1,ibool(NGLLX_M,1,NGLLZ_M,ispec)-1, &
                 ibool(1,NGLLY_M,NGLLX_M,ispec)-1,ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)-1, &
                 ibool(1,1,1,ispec)-1,ibool(NGLLX_M,1,1,ispec)-1, &
                 ibool(1,NGLLY_M,1,ispec)-1,ibool(NGLLX_M,NGLLY_M,1,ispec)-1
    enddo

    ! ************* generate element data values ******************

    ! output OpenDX header for data
    write(IOUT,*) 'attribute "element type" string "cubes"'

    write(IOUT,*) 'attribute "ref" string "positions"'
    write(IOUT,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

    ! loop on all the elements
    do ispec = 1,nspec
      write(IOUT,*)  ispec_material_id(ispec)
    enddo

    write(IOUT,*) 'attribute "dep" string "connections"'
    write(IOUT,*) 'object "irregular positions irregular connections" class field'
    write(IOUT,*) 'component "positions" value 1'
    write(IOUT,*) 'component "connections" value 2'
    write(IOUT,*) 'component "data" value 3'
    write(IOUT,*) 'end'

    close(IOUT)

  endif

  if (CREATE_VTK_FILES) then

    ! vtk file output
    filename = prname(1:len_trim(prname))//'mesh.vtk'
    open(IOUT,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error opening file mesh.INP for CREATE_ABAQUS_FILES'
    endif

    write(IOUT,'(a)') '# vtk DataFile Version 3.1'
    write(IOUT,'(a)') 'material model VTK file'
    write(IOUT,'(a)') 'ASCII'
    write(IOUT,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(IOUT, '(a,i12,a)') 'POINTS ', nglob, ' float'
    do ipoin = 1,nglob
      write(IOUT,*) sngl(nodes_coords(ipoin,1)),sngl(nodes_coords(ipoin,2)),sngl(nodes_coords(ipoin,3))
    enddo
    write(IOUT,*)

    ! note: indices for vtk start at 0
    write(IOUT,'(a,i12,i12)') "CELLS ",nspec,nspec*9
    do ispec=1,nspec
      write(IOUT,'(9i12)') 8, &
            ibool(1,1,1,ispec)-1,ibool(NGLLX_M,1,1,ispec)-1, &
            ibool(NGLLX_M,NGLLY_M,1,ispec)-1,ibool(1,NGLLY_M,1,ispec)-1, &
            ibool(1,1,NGLLZ_M,ispec)-1,ibool(NGLLX_M,1,NGLLZ_M,ispec)-1, &
            ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)-1,ibool(1,NGLLY_M,NGLLZ_M,ispec)-1
    enddo
    write(IOUT,*)

    ! type: hexahedrons
    write(IOUT,'(a,i12)') "CELL_TYPES ",nspec
    write(IOUT,'(6i12)') (12,ispec=1,nspec)
    write(IOUT,*)

    write(IOUT,'(a,i12)') "CELL_DATA ",nspec
    write(IOUT,'(a)') "SCALARS material_id float"
    write(IOUT,'(a)') "LOOKUP_TABLE default"
    do ispec = 1,nspec
      write(IOUT,*) ispec_material_id(ispec)
    enddo
    write(IOUT,*)
    close(IOUT)

  endif

  call synchronize_all()

  end subroutine create_visual_files

