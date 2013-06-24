!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  program create_slice_VTK

! this programs creates a vtk file that specifies the (velocity) model values on each element point,
! rather than on global points. this will lead to differences especially where a (velocity) discontinuity
! from one element to the other exists. in such cases, this output file should be more accurate
! in how it is visualized.
!
! creates for each slice and each region a new vtk-file, doesn't combine different slices into one file
! (only outputs values on corner points of each element)
!
! for compilation, this file 'create_slice_VTK.f90' needs the constants.h file located in SPECFEM3D/src/shared
!
! in this current directory SPECFEM3D/utils/Visualization/Paraview:
! compile with:
!   > f90 -o xcreate_slice_VTK create_slice_VTK.f90
! to run :
!   > ./xcreate_slice_VTK my_slices.txt vp ~/DATABASES_MPI/ ~/DATABASES_MPI ~/OUTPUT_FILES/
!
! (or see usage below)
!
  implicit none

  include "../../../src/shared/constants.h"

  character(len=256) :: sline, arg(5), filename, indir, outdir
  character(len=256) :: prname, prname_lp
  character(len=256) :: mesh_file,local_data_file

  integer, dimension(300) :: node_list
  integer :: iproc, num_node, i,ios, it, ier
  integer :: njunk

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: data
  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  integer, dimension(:,:,:,:),allocatable :: ibool

  integer :: NSPEC_AB, NGLOB_AB


  ! starts here---------------------------------------------------------------------------------------
  do i = 1, 4
    call getarg(i,arg(i))
    if (i < 4 .and. trim(arg(i)) == '') then
      print *
      print *, ' Usage: xcreate_slice_VTK slice_list filename input_dir output_dir'
      print *, ' '
      print *, '   - slice_list:    file containing slice/proc ids '
      print *, '   - filename:    looks for filename.bin must be array of (NGLLX,NGLLY,NGLLZ,nspec)'
      print *, '   - input_dir:    includes proc***_external_mesh.bin and proc****filename.bin '
      print *, '   - output_dir:    output mesh files go to here '
      print *
      stop ' Reenter command line options'
    endif
  enddo


  ! get slices id
  num_node = 0
  open(unit = 20, file = trim(arg(1)), status = 'old',iostat = ios)
  if (ios /= 0) then
    print*,'no file: ',trim(arg(1))
    stop 'Error opening slices file'
  endif

  do while (1 == 1)
    read(20,'(a)',iostat=ios) sline
    if (ios /= 0) exit
    read(sline,*,iostat=ios) njunk
    if (ios /= 0) exit
    num_node = num_node + 1
    node_list(num_node) = njunk
  enddo
  close(20)
  print *, 'slice list: '
  print *, node_list(1:num_node)
  print *, ' '

  ! file to collect
  filename = arg(2)

  ! input and output dir
  indir= arg(3)
  outdir = arg(4)


  ! write points information
  do it = 1, num_node

    iproc = node_list(it)

    print *, ' '
    print *, 'Reading slice ', iproc

    ! gets number of elements and global points for this partition
    write(prname_lp,'(a,i6.6,a)') trim(indir)//'/proc',iproc,'_'
    open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin',&
          status='old',action='read',form='unformatted')
    read(27) NSPEC_AB
    read(27) NGLOB_AB

    ! ibool file
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibool'
    read(27) ibool

    ! global point arrays
    allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array xstore etc.'
    read(27) xstore
    read(27) ystore
    read(27) zstore
    close(27)


    ! filename.bin
    ! data file
    write(prname,'(a,i6.6,a)') trim(indir)//'proc',iproc,'_'
    local_data_file = trim(prname) // trim(filename) // '.bin'
    open(unit = 28,file = trim(local_data_file),status='old',&
          action='read', iostat = ios,form ='unformatted')
    if (ios /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array data'
    read(28) data
    close(28)

    print *, trim(local_data_file)
    print *,'  min/max value: ',minval(data),maxval(data)
    print *


    write(mesh_file,'(a,i6.6,a)') trim(outdir)//'/' // 'slice_',iproc,'_'//trim(filename)
    print *, trim(mesh_file)


    ! writes out vtk file
    call write_VTK_data_gll_cr(NSPEC_AB,NGLOB_AB, &
            xstore,ystore,zstore,ibool, &
            data,mesh_file)

    deallocate(ibool,xstore,ystore,zstore,data)

  enddo  ! all slices for points

  print *, 'Done writing slice files'
  print *, ' '


  end program create_slice_VTK

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_gll_cr(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            gll_data,prname_file)

! external mesh routine for saving vtk files for custom_real values on all gll points

  implicit none

  include "../../../src/shared/constants.h"

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! gll data values array
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: gll_data

! file name
  character(len=256) prname_file

  integer :: ispec,i

! write source and receiver VTK files for Paraview
  write(IMAIN,*) '  vtk file: '
  write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'

  ! writes out all points for each element, not just global ones
  write(IOVTK, '(a,i12,a)') 'POINTS ', nspec*8, ' float'
  do ispec=1,nspec
    i = ibool(1,1,1,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

    i = ibool(NGLLX,1,1,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

    i = ibool(NGLLX,NGLLY,1,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

    i = ibool(1,NGLLY,1,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

    i = ibool(1,1,NGLLZ,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

    i = ibool(NGLLX,1,NGLLZ,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

    i = ibool(NGLLX,NGLLY,NGLLZ,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

    i = ibool(1,NGLLY,NGLLZ,ispec)
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,(ispec-1)*8,(ispec-1)*8+1,(ispec-1)*8+2,(ispec-1)*8+3,&
          (ispec-1)*8+4,(ispec-1)*8+5,(ispec-1)*8+6,(ispec-1)*8+7
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,'(6i12)') (12,ispec=1,nspec)
  write(IOVTK,*) ""

  ! writes out gll-data (velocity) for each element point
  write(IOVTK,'(a,i12)') "POINT_DATA ",nspec*8
  write(IOVTK,'(a)') "SCALARS gll_data float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    i = ibool(1,1,1,ispec)
    write(IOVTK,'(3e18.6)') gll_data(1,1,1,ispec)

    i = ibool(NGLLX,1,1,ispec)
    write(IOVTK,'(3e18.6)') gll_data(NGLLX,1,1,ispec)

    i = ibool(NGLLX,NGLLY,1,ispec)
    write(IOVTK,'(3e18.6)') gll_data(NGLLX,NGLLY,1,ispec)

    i = ibool(1,NGLLY,1,ispec)
    write(IOVTK,'(3e18.6)') gll_data(1,NGLLY,1,ispec)

    i = ibool(1,1,NGLLZ,ispec)
    write(IOVTK,'(3e18.6)') gll_data(1,1,NGLLZ,ispec)

    i = ibool(NGLLX,1,NGLLZ,ispec)
    write(IOVTK,'(3e18.6)') gll_data(NGLLX,1,NGLLZ,ispec)

    i = ibool(NGLLX,NGLLY,NGLLZ,ispec)
    write(IOVTK,'(3e18.6)') gll_data(NGLLX,NGLLY,NGLLZ,ispec)

    i = ibool(1,NGLLY,NGLLZ,ispec)-1
    write(IOVTK,'(3e18.6)') gll_data(1,NGLLY,NGLLZ,ispec)
  enddo
  write(IOVTK,*) ""

  close(IOVTK)


  end subroutine write_VTK_data_gll_cr

