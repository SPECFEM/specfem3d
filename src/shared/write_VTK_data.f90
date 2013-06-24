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


! routine for saving vtk file holding integer flag on each spectral element

  subroutine write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)


  implicit none

  include "constants.h"

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element flag array
  integer, dimension(nspec) :: elem_flag

  ! file name
  character(len=256) :: prname_file
  !character(len=2), optional, intent(in) :: str_id

  ! local parameters
  integer :: ispec,i

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,'(6i12)') (12,ispec=1,nspec)
  write(IOVTK,*) ""

  write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
  !if( present( str_id ) ) then
  !  write(IOVTK,'(a)') "SCALARS elem_flag_"//str_id//" integer"
  !else
  !  write(IOVTK,'(a)') "SCALARS elem_flag integer"
  !endif
  write(IOVTK,'(a)') "SCALARS elem_flag integer"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOVTK,*) elem_flag(ispec)
  enddo
  write(IOVTK,*) ""
  close(IOVTK)


  end subroutine write_VTK_data_elem_i

!=====================================================================


! routine for saving vtk file holding logical flag on each spectral element

  subroutine write_VTK_data_elem_l(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)


  implicit none

  include "constants.h"

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! element flag array
  logical, dimension(nspec) :: elem_flag
  integer :: ispec,i

! file name
  character(len=256) prname_file

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,'(6i12)') (12,ispec=1,nspec)
  write(IOVTK,*) ""

  write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOVTK,'(a)') "SCALARS elem_flag integer"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    if( elem_flag(ispec) .eqv. .true. ) then
      write(IOVTK,*) 1
    else
      write(IOVTK,*) 0
    endif
  enddo
  write(IOVTK,*) ""
  close(IOVTK)


  end subroutine write_VTK_data_elem_l


!=============================================================

! external mesh routine for saving vtk files for custom_real values on all gll points

  subroutine write_VTK_data_gll_cr(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            gll_data,prname_file)

  implicit none

  include "constants.h"

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! gll data values array
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: gll_data

! masking arrays (takes first data value assigned on a global point, ignores any data values later on for the same global point)
  real, dimension(:),allocatable :: flag_val
  logical, dimension(:),allocatable :: mask_ibool

! file name
  character(len=256) prname_file

  integer :: ispec,i,j,k,ier,iglob

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,'(6i12)') (12,ispec=1,nspec)
  write(IOVTK,*) ""

  ! iflag field on global nodeset
  allocate(mask_ibool(nglob),flag_val(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating mask'

  mask_ibool = .false.
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if( .not. mask_ibool(iglob) ) then
            flag_val(iglob) = gll_data(i,j,k,ispec)
            mask_ibool(iglob) = .true.
          endif
        enddo
      enddo
    enddo
  enddo

  write(IOVTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOVTK,'(a)') "SCALARS gll_data float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOVTK,*) flag_val(i)
  enddo
  write(IOVTK,*) ""

  close(IOVTK)


  end subroutine write_VTK_data_gll_cr

!=============================================================

! external mesh routine for saving vtk files for integer values on all gll points

  subroutine write_VTK_data_gll_i(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            gll_data,prname_file)

  implicit none

  include "constants.h"

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! gll data values array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: gll_data

! masking arrays (takes first data value assigned on a global point, ignores any data values later on for the same global point)
  real, dimension(:),allocatable :: flag_val
  logical, dimension(:),allocatable :: mask_ibool

! file name
  character(len=256) prname_file

  integer :: ispec,i,j,k,ier,iglob

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,'(6i12)') (12,ispec=1,nspec)
  write(IOVTK,*) ""

  ! iflag field on global nodeset
  allocate(mask_ibool(nglob),flag_val(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating mask'

  mask_ibool = .false.
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if( .not. mask_ibool(iglob) ) then
            flag_val(iglob) = gll_data(i,j,k,ispec)
            mask_ibool(iglob) = .true.
          endif
        enddo
      enddo
    enddo
  enddo

  write(IOVTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOVTK,'(a)') "SCALARS gll_data float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOVTK,*) flag_val(i)
  enddo
  write(IOVTK,*) ""

  close(IOVTK)


  end subroutine write_VTK_data_gll_i

!=============================================================

! external mesh routine for saving vtk files for points locations

  subroutine write_VTK_data_points(nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy, &
            points_globalindices,num_points_globalindices, &
            prname_file)

  implicit none

  include "constants.h"

  integer :: nglob

! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! gll data values array
  integer :: num_points_globalindices
  integer, dimension(num_points_globalindices) :: points_globalindices

! file name
  character(len=256) prname_file

  integer :: i,iglob

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', num_points_globalindices, ' float'
  do i=1,num_points_globalindices
    iglob = points_globalindices(i)
    if( iglob <= 0 .or. iglob > nglob ) then
      print*,'error: '//prname_file(1:len_trim(prname_file))//'.vtk'
      print*,'error global index: ',iglob,i
      close(IOVTK)
      stop 'error vtk points file'
    endif

    write(IOVTK,'(3e18.6)') xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
  enddo
  write(IOVTK,*) ""

  close(IOVTK)


  end subroutine write_VTK_data_points


!=============================================================

! external mesh routine for saving vtk files for points locations

  subroutine write_VTK_data_elem_vectors(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_vector,prname_file)


  implicit none

  include "constants.h"

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element flag array
  real(kind=CUSTOM_REAL), dimension(3,nspec) :: elem_vector
  integer :: ispec,i

  ! file name
  character(len=256) prname_file

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,'(6i12)') (12,ispec=1,nspec)
  write(IOVTK,*) ""

  ! vector data for each cell
  write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOVTK,'(a)') "VECTORS _vectors_ float"
  do i=1,nspec
    write(IOVTK,*) elem_vector(1,i),elem_vector(2,i),elem_vector(3,i)
  enddo

  write(IOVTK,*) ""
  close(IOVTK)


  end subroutine write_VTK_data_elem_vectors

!=============================================================

  subroutine write_VTK_data_elem_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)


  implicit none

  include "constants.h"

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element flag array
  real(kind=CUSTOM_REAL), dimension(nspec) :: elem_flag

  ! file name
  character(len=256) :: prname_file
  !character(len=2), optional, intent(in) :: str_id

  ! local parameters
  integer :: ispec,i

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,'(6i12)') (12,ispec=1,nspec)
  write(IOVTK,*) ""

  write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
  !if( present( str_id ) ) then
  !  write(IOVTK,'(a)') "SCALARS elem_flag_"//str_id//" integer"
  !else
  !  write(IOVTK,'(a)') "SCALARS elem_flag integer"
  !endif
  write(IOVTK,'(a)') "SCALARS elem_val float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOVTK,*) elem_flag(ispec)
  enddo
  write(IOVTK,*) ""
  close(IOVTK)


  end subroutine write_VTK_data_elem_cr
