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

! routine for saving vtk file holding integer flag on each spectral element

! note: due to some compilers output, we usually put formatting specifiers here (e.g., (3e20.12)).
!
!       example: some compilers tend to output a line like:
!                   -1500.  -1500.   1312.5
!                with same numbers as:
!                   2*-1500., 1312.5
!
!                this will corrupt the VTK files and become unreadable.
!                thus, by putting
!
!                   write(IOUT_VTK,'(3e20.12)') ..
!
!                the lines becomes fine again. unfortunately, we need to assume a "good" number range
!                for this. since we could have SPECFEM simulations for microscales on micrometer up to
!                macroscales for kilometers, this range specifier might not always be optimal.
!                maybe we find a better way in future...

  subroutine write_VTK_data_elem_i(nspec,nglob, &
                                   xstore,ystore,zstore,ibool, &
                                   elem_flag,prname_file)

  use constants, only: MAX_STRING_LEN,IOUT_VTK,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

  ! element flag array
  integer, dimension(nspec) :: elem_flag

  ! file name
  character(len=MAX_STRING_LEN) :: prname_file

  ! local parameters
  integer :: ispec,i

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOUT_VTK,'(3e20.12)') xstore(i),ystore(i),zstore(i)
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_flag integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_flag(ispec)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_i

!
!-------------------------------------------------------------------------------------------------
!

! routine for saving vtk file holding logical flag on each spectral element

  subroutine write_VTK_data_elem_l(nspec,nglob, &
                                   xstore,ystore,zstore,ibool, &
                                   elem_flag,prname_file)

  use constants, only: MAX_STRING_LEN,IOUT_VTK,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

! element flag array
  logical, dimension(nspec) :: elem_flag
  integer :: ispec,i

! file name
  character(len=MAX_STRING_LEN) :: prname_file

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOUT_VTK,'(3e20.12)') xstore(i),ystore(i),zstore(i)
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_flag integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    if (elem_flag(ispec) .eqv. .true.) then
      write(IOUT_VTK,*) 1
    else
      write(IOUT_VTK,*) 0
    endif
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_l


!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving vtk files for CUSTOM_REAL values on all GLL points

  subroutine write_VTK_data_gll_cr(nspec,nglob, &
                                   xstore,ystore,zstore,ibool, &
                                   gll_data,prname_file)

  use constants, only: MAX_STRING_LEN,IOUT_VTK,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

! GLL data values array
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: gll_data

! masking arrays (takes first data value assigned on a global point, ignores any data values later on for the same global point)
  real, dimension(:),allocatable :: flag_val
  logical, dimension(:),allocatable :: mask_ibool

! file name
  character(len=MAX_STRING_LEN) :: prname_file

  integer :: ispec,i,j,k,ier,iglob

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOUT_VTK,'(3e20.12)') xstore(i),ystore(i),zstore(i)
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  ! iflag field on global nodeset
  if (.not. allocated(mask_ibool)) then
     allocate(mask_ibool(nglob),flag_val(nglob),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 1234')
     if (ier /= 0) stop 'error allocating mask'
  endif

  mask_ibool = .false.
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
            flag_val(iglob) = gll_data(i,j,k,ispec)
            mask_ibool(iglob) = .true.
          endif
        enddo
      enddo
    enddo
  enddo

  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOUT_VTK,'(a)') "SCALARS gll_data float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) flag_val(i)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_gll_cr

!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving vtk files for integer values on all GLL points

  subroutine write_VTK_data_gll_i(nspec,nglob, &
                                  xstore,ystore,zstore,ibool, &
                                  gll_data,prname_file)

  use constants, only: MAX_STRING_LEN,IOUT_VTK,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

! GLL data values array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: gll_data

! masking arrays (takes first data value assigned on a global point, ignores any data values later on for the same global point)
  real, dimension(:),allocatable :: flag_val
  logical, dimension(:),allocatable :: mask_ibool

! file name
  character(len=MAX_STRING_LEN) :: prname_file

  integer :: ispec,i,j,k,ier,iglob

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOUT_VTK,'(3e20.12)') xstore(i),ystore(i),zstore(i)
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  ! iflag field on global nodeset
  allocate(mask_ibool(nglob),flag_val(nglob),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1235')
  if (ier /= 0) stop 'error allocating mask'

  mask_ibool = .false.
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
            flag_val(iglob) = gll_data(i,j,k,ispec)
            mask_ibool(iglob) = .true.
          endif
        enddo
      enddo
    enddo
  enddo

  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOUT_VTK,'(a)') "SCALARS gll_data float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) flag_val(i)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_gll_i

!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving vtk files for points locations

  subroutine write_VTK_data_points(nglob, &
                                   xstore,ystore,zstore, &
                                   points_globalindices,num_points_globalindices, &
                                   prname_file)

  use constants, only: MAX_STRING_LEN,IOUT_VTK,CUSTOM_REAL

  implicit none

  integer :: nglob

! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

! GLL data values array
  integer :: num_points_globalindices
  integer, dimension(num_points_globalindices) :: points_globalindices

! file name
  character(len=MAX_STRING_LEN) :: prname_file

  integer :: i,iglob

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', num_points_globalindices, ' float'
  do i=1,num_points_globalindices
    iglob = points_globalindices(i)
    if (iglob <= 0 .or. iglob > nglob) then
      print *,'error: '//prname_file(1:len_trim(prname_file))//'.vtk'
      print *,'error global index: ',iglob,i
      close(IOUT_VTK)
      stop 'error vtk points file'
    endif

    write(IOUT_VTK,'(3e20.12)') xstore(iglob),ystore(iglob),zstore(iglob)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_points


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_points_elem(NGNOD,xelm,yelm,zelm,val,prname_file)

! single element output as points

  use constants, only: MAX_STRING_LEN,IOUT_VTK

  implicit none

  integer,intent(in) :: NGNOD

  ! coordinates
  double precision, dimension(NGNOD),intent(in) :: xelm,yelm,zelm
  double precision, intent(in) :: val

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: i,itype

  ! write source and receiver VTK files for Paraview
  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', NGNOD, ' float'
  do i = 1,NGNOD
    ! note: the format specifiers is automatic to have a more readable file, but might have issues on some compiler outputs.
    !       see the remark at the top of the file.
    write(IOUT_VTK,*) xelm(i),yelm(i),zelm(i)
  enddo
  write(IOUT_VTK,*) ''

  ! single element
  !
  ! cell type: hexahedrons
  !
  ! VTK type definition numbers:
  !   https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
  ! ordering images see e.g.:
  !   https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/
  if (NGNOD == 8) then
    ! hex with 8 nodes
    !
    ! VTK_HEXAHEDRON == 12
    !
    ! The hexahedron is defined by an ordered list of eight points.
    ! ordering: 0,1,2,3  -> bottom (left,front),(right,front),(right,back),(left,back)
    !           4,5,6,7  -> top    (left,front),(right,front),(right,back),(left,back)
    itype = 12
  else if (NGNOD == 27) then
    ! hex with 27 nodes
    !
    ! VTK_TRIQUADRATIC_HEXAHEDRON == 29
    !
    ! The tri-quadratic hexahedron is a primary three-dimensional cell. It is defined by twenty-seven points.
    ! The first eight points are located at the vertices of the hexahedron;
    ! the next twelve are located in the middle of each of the twelves edges;
    ! the next six are located in the center of each faces; and the last one is located in the center of the hexahedron
    ! ordering: 0,1,2,3  -> bottom (left,front),(right,front),(right,back),(left,back)
    !           4,5,6,7  -> top    (left,front),(right,front),(right,back),(left,back)
    !           8,9,10,11   -> bottom edges mid-points (bottom,front),(bottom,right),(bottom,back),(bottom,left)
    !           12,13,14,15 -> top edges mid-points (top,front),(top,right),(top,back),(top,left)
    !           16,17,18,19 -> bottom-to-top edges mid-points (front,left),(front,right),(back,right),(back,left)
    !           20,21,22,23,24,25 -> face centers (left),(right),(front),(back),(bottom),(top)
    !           27          -> center point of hexahedron
    itype = 29
  else
    print *,'Invalid NGNOD ',NGNOD,'in routine write_VTK_data_points_elem() routine'
    stop 'Invalid NGNOD in write_VTK_data_points_elem()'
  endif

  ! SPECFEM order:
  ! uses ordering:
  !   corners      bottom (left,front),(right,front),(right,back),(left,back)  -> vtk 0,1,2,3
  !                top    (left,front),(right,front),(right,back),(left,back)  ->     4,5,6,7     same as vtk
  !
  !   edges    (bottom edges: front,right,back,left)
  !            (bottom-to-top edges: front left,right, back right,left)
  !            (top edges: front,right,back,left)
  !            -> vtk needs bottom,top,mid edges   -> 8,9,10,11, 16,17,18,19, 12,13,14,15
  !   faces    (bottom),(front),(right),(back),(left),(top)
  !            -> vtk needs (left),(right),(front),(back),(bottom),(top) -> 24,22,21,23,20,25
  !   center node -> 26
  !
  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",1,1*(NGNOD + 1)
  if (NGNOD == 27) then
    ! hex27 single element order
    write(IOUT_VTK,'(28i4)') NGNOD,0,1,2,3, 4,5,6,7, 8,9,10,11, 16,17,18,19, 12,13,14,15, 24,22,21,23,20,25, 26
  else
    ! hex8 single element
    write(IOUT_VTK,'(9i12)') NGNOD,(i-1,i=1,NGNOD)
  endif
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",1
  write(IOUT_VTK,'(i12)') itype
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",1
  write(IOUT_VTK,'(a)') "SCALARS elem_val float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  write(IOUT_VTK,*) val
  write(IOUT_VTK,*) ''

  close(IOUT_VTK)

  end subroutine write_VTK_data_points_elem

!
!-------------------------------------------------------------------------------------------------
!


! external mesh routine for saving vtk files for points locations

  subroutine write_VTK_data_elem_vectors(nspec,nglob, &
                                         xstore,ystore,zstore,ibool, &
                                         elem_vector,prname_file)

  use constants, only: MAX_STRING_LEN,IOUT_VTK,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

  ! element flag array
  real(kind=CUSTOM_REAL), dimension(3,nspec) :: elem_vector
  integer :: ispec,i

  ! file name
  character(len=MAX_STRING_LEN) :: prname_file

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOUT_VTK,'(3e20.12)') xstore(i),ystore(i),zstore(i)
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  ! vector data for each cell
  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "VECTORS _vectors_ float"
  do i=1,nspec
    write(IOUT_VTK,*) elem_vector(1,i),elem_vector(2,i),elem_vector(3,i)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_vectors

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_elem_cr(nspec,nglob, &
                                    xstore,ystore,zstore,ibool, &
                                    elem_data,prname_file)

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

  ! element flag array
  real(kind=CUSTOM_REAL), dimension(nspec) :: elem_data

  ! file name
  character(len=MAX_STRING_LEN) :: prname_file
  !character(len=2), optional, intent(in) :: str_id

  ! local parameters
  integer :: ispec,i

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e20.12)') xstore(i),ystore(i),zstore(i)
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_val float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_data(ispec)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_cr


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTU_data_elem_cr_binary(nspec,nglob, &
                                           xstore,ystore,zstore,ibool, &
                                           elem_data,prname_file)

! saves vtu files in binary format, for CUSTOM_REAL values on all spectral elements

  use constants, only: CUSTOM_REAL,SIZE_REAL,SIZE_INTEGER,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer,intent(in) :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  ! element data values array
  real(kind=CUSTOM_REAL), dimension(nspec),intent(in) :: elem_data

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,it,ier

  ! local parameters
  integer :: len_bytes,offset
  character(len=MAX_STRING_LEN) :: var_name
  character(len=16) :: str1,str_endian
  character(len=12) :: str2,str3
  character(len=1),parameter :: LF = achar(10) ! line feed character

  ! endianness - determine endianness at run time:
  ! https://www.ibm.com/developerworks/aix/library/au-endianc/index.html
  !
  ! for the Fortran version:
  ! given the integer value of 1, big endian uses 00 00 00 01 and little endian uses 01 00 00 00 as bit representation.
  ! using transfer(1,'a') takes the bit-representation of integer value 1 (multi-byte, usually 32-bit)
  ! and interprets it as a character type (of 1-byte length).
  ! thus, looking at the first byte, it is either 0 or 1, respectively.
  ! finally, ichar(c) takes a character and returns its position number.
  logical, parameter :: is_big_endian = (ichar(transfer(1,'a')) == 0)

  ! note: VTK by default assumes binary data is big endian for "legacy" VTK format,
  !       we use here the new .vtu file with XML format. in this case, we can specify the endianness by byte_order=".."

  if (is_big_endian) then
    str_endian = 'BigEndian'
  else
    str_endian = 'LittleEndian'
  endif

  ! variable name
  var_name = 'elem_val'

  ! numbers as strings
  write(str1,'(i16)') nglob
  write(str2,'(i12)') nspec

  ! data offset for appended datasets
  offset = 0
  write(str3,'(i12)') offset

  ! opens unstructured grid file as binary file
  ! convert='BIG_ENDIAN' not needed, will be specified in XML format
  open(IOUT_VTK,file=trim(prname_file)//'.vtu',access='stream',form='unformatted', &
         status='unknown',action='write',iostat=ier)
  if (ier /= 0 ) then
    print *, 'Error opening VTU output file: ',trim(prname_file)//'.vtu'
    stop 'Error opening VTU output file'
  endif

  ! XML file header
  ! see: https://www.vtk.org/Wiki/VTK_XML_Formats
  write(IOUT_VTK) '<?xml version="1.0"?>' // LF
  write(IOUT_VTK) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="'// trim(str_endian) // '">' // LF
  write(IOUT_VTK) '<UnstructuredGrid>' // LF
  write(IOUT_VTK) '<Piece NumberOfPoints="'// str1 // '" NumberOfCells="' // str2 // '">' // LF

  ! points
  write(IOUT_VTK) '<Points>' // LF
  ! binary format: not working properly yet - using appended instead
  !write(IOUT_VTK) '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary" encoding="raw">' // LF
  !! number of bytes to follow
  !! see format description: https://www.vtk.org/Wiki/VTK_XML_Formats#Uncompressed_Data
  !len_bytes = 3 * nglob * SIZE_REAL
  !write(IOUT_VTK) len_bytes
  !do i = 1,nglob
  !  write(IOUT_VTK) real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
  !enddo
  !write(IOUT_VTK) '</DataArray>' // LF
  !
  ! appended format:
  write(IOUT_VTK) '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="' &
                   // str3 // '"/>' // LF
  ! updates offset
  ! array length in bytes
  len_bytes = 3 * nglob * SIZE_REAL
  ! new offset: data array size (len_bytes) + 4 bytes for length specifier at the beginning
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset
  write(IOUT_VTK) '</Points>'//LF

  ! cells
  write(IOUT_VTK) '<Cells>' // LF
  ! connectivity
  write(IOUT_VTK) '<DataArray type="Int32" Name="connectivity" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = 8 * nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset

  ! offsets
  write(IOUT_VTK) '<DataArray type="Int32" Name="offsets" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset

  ! type: hexahedrons
  write(IOUT_VTK) '<DataArray type="Int32" Name="types" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset
  write(IOUT_VTK) '</Cells>' // LF

  ! cell data
  write(IOUT_VTK) '<CellData Scalars="Scalars_">' // LF
  write(IOUT_VTK) '<DataArray type="Float32" Name="' // trim(var_name) // '" format="appended" offset="' &
                  // str3 // '"/>' // LF
  write(IOUT_VTK) '</CellData>' // LF

  ! empty point data values
  write(IOUT_VTK) '<PointData>' // '</PointData>' // LF

  ! finishes XML file
  write(IOUT_VTK) '</Piece>' // LF
  write(IOUT_VTK) '</UnstructuredGrid>' // LF

  ! in case of appended data format, with offsets:
  !write(IOUT_VTK) '<AppendedData encoding="base64">' // LF
  write(IOUT_VTK) '<AppendedData encoding="raw">' // LF
  write(IOUT_VTK) '_'
  ! points
  len_bytes = 3 * nglob * SIZE_REAL
  write(IOUT_VTK) len_bytes
  do i = 1,nglob
    write(IOUT_VTK) real(xstore(i),kind=4),real(ystore(i),kind=4),real(zstore(i),kind=4)
  enddo
  ! cells
  ! connectivity
  ! number of bytes to follow
  len_bytes = 8 * nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! note: indices for VTK start at 0
  do ispec = 1,nspec
    write(IOUT_VTK) ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
      ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  ! offsets
  ! number of bytes to follow
  len_bytes = nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! offsets (8 points each)
  write(IOUT_VTK) ((it*8),it = 1,nspec)
  ! types
  ! number of bytes to follow
  len_bytes = nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! type for hex elements is 12
  write(IOUT_VTK) (12,it = 1,nspec)

  ! cell data values
  ! number of bytes to follow
  len_bytes = nspec * SIZE_REAL
  write(IOUT_VTK) len_bytes
  ! data values
  do ispec = 1,nspec
    write(IOUT_VTK) real(elem_data(ispec),kind=4)
  enddo

  write(IOUT_VTK) '</AppendedData>' // LF
  write(IOUT_VTK) '</VTKFile>' // LF
  close(IOUT_VTK)

  end subroutine write_VTU_data_elem_cr_binary

!
!------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_ngnod_elem_i(nspec,nglob,NGNOD, &
                                         xstore_db,ystore_db,zstore_db,elmnts, &
                                         elem_flag,filename)

  use constants, only: IOUT_VTK,MAX_STRING_LEN

  implicit none

  integer :: nspec,nglob,NGNOD

  ! global coordinates
  integer, dimension(NGNOD,nspec) :: elmnts
  double precision, dimension(nglob) :: xstore_db,ystore_db,zstore_db

  ! element flag array
  integer, dimension(nspec) :: elem_flag

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: ispec,i

  ! safety check
  if (NGNOD < 8) stop 'Error invalid NGNOD in write_VTK_data_ngnod_elem_i routine'

  ! write source and receiver VTK files for Paraview
  open(IOUT_VTK,file=filename(1:len_trim(filename))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e20.12)') xstore_db(i),ystore_db(i),zstore_db(i)
  enddo
  write(IOUT_VTK,*) ""

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
                             elmnts(1,ispec)-1,elmnts(2,ispec)-1,elmnts(3,ispec)-1,elmnts(4,ispec)-1, &
                             elmnts(5,ispec)-1,elmnts(6,ispec)-1,elmnts(7,ispec)-1,elmnts(8,ispec)-1
  enddo
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ""

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_flag integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_flag(ispec)
  enddo
  write(IOUT_VTK,*) ""
  close(IOUT_VTK)

  end subroutine write_VTK_data_ngnod_elem_i


!------------------------------------------------------------------------------------
!
! VTK routines for meshfem3D
!
!------------------------------------------------------------------------------------

  subroutine write_VTK_data_elem_cr_meshfem(nspec,nglob,xstore_db,ystore_db,zstore_db,ibool, &
                                            elem_data,filename)

! special routine for meshfem3D with simpler mesh arrays

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGNOD_EIGHT_CORNERS,IOUT_VTK

  implicit none

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC) :: ibool
  double precision, dimension(nglob) :: xstore_db,ystore_db,zstore_db

  ! element flag array
  real(kind=CUSTOM_REAL), dimension(nspec) :: elem_data

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: ispec,i

  open(IOUT_VTK,file=trim(filename),status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,*) sngl(xstore_db(i)),sngl(ystore_db(i)),sngl(zstore_db(i))
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
                             ibool(1,ispec)-1,ibool(2,ispec)-1,ibool(4,ispec)-1,ibool(3,ispec)-1, &
                             ibool(5,ispec)-1,ibool(6,ispec)-1,ibool(8,ispec)-1,ibool(7,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedra
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS skewness float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_data(ispec)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_cr_meshfem

!
!------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_elem_i_earthmesh(nspec,nglob,x_mesh,y_mesh,z_mesh, &
                                             iglob,EtoV,elem_data,filename)

! special routine for meshfem3D with simpler earth_mesh arrays

  use constants, only: MAX_STRING_LEN,IOUT_VTK

  implicit none

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(nglob) :: iglob
  integer, dimension(8,nspec) :: EtoV
  double precision, dimension(nglob) :: x_mesh,y_mesh,z_mesh

  ! element flag array
  integer, dimension(nspec) :: elem_data

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: i

  open(IOUT_VTK,file=trim(filename),status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
     !!if (ifseg(i)) then
        write(IOUT_VTK,*) x_mesh(i), y_mesh(i), z_mesh(i)
     !!endif
  enddo
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do i = 1,nspec
     write(27,'(9i15)') 8, &
                        iglob((EtoV(1, i)))-1, iglob((EtoV(2, i)))-1, iglob((EtoV(3, i)))-1, iglob((EtoV(4, i)))-1, &
                        iglob((EtoV(5, i)))-1, iglob((EtoV(6, i)))-1, iglob((EtoV(7, i)))-1, iglob((EtoV(8, i)))-1
  enddo
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,i=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_val float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nspec
      write(IOUT_VTK,*) elem_data(i)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_i_earthmesh

!
!------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_elem_i_meshfemCPML(nglob,nspec,NGLL,nodes_coords,ibool, &
                                               nspec_CPML,CPML_regions,is_CPML,filename)

! special routine for meshfem3D with simpler earth_mesh arrays

  use constants, only: MAX_STRING_LEN,NDIM,IOUT_VTK

  implicit none

  integer :: nspec,nglob,NGLL

  ! global coordinates
  integer, dimension(NGLL,NGLL,NGLL,nspec) :: ibool
  double precision, dimension(nglob,NDIM) :: nodes_coords

  ! element flag array
  integer :: nspec_CPML
  integer, dimension(nspec_CPML) :: CPML_regions
  logical, dimension(nspec) :: is_CPML

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: ispec,i,ispec_CPML

  open(IOUT_VTK,file=trim(filename),status='unknown')

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,*) sngl(nodes_coords(i,1)),sngl(nodes_coords(i,2)),sngl(nodes_coords(i,3))
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
                             ibool(1,1,1,ispec)-1,ibool(NGLL,1,1,ispec)-1, &
                             ibool(NGLL,NGLL,1,ispec)-1,ibool(1,NGLL,1,ispec)-1, &
                             ibool(1,1,NGLL,ispec)-1,ibool(NGLL,1,NGLL,ispec)-1, &
                             ibool(NGLL,NGLL,NGLL,ispec)-1,ibool(1,NGLL,NGLL,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_flag integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  ispec_CPML = 0
  do ispec = 1,nspec
    if (is_CPML(ispec) .eqv. .true.) then
      ispec_CPML = ispec_CPML + 1
      write(IOUT_VTK,*) CPML_regions(ispec_CPML)
    else
      write(IOUT_VTK,*) 0
    endif
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_i_meshfemCPML

!
!------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_elem_i_meshfem(nspec,NGLL,xstore_mesh,ystore_mesh,zstore_mesh, &
                                           ibool,elem_data,filename)

! special routine for meshfem3D with simpler earth_mesh arrays

  use constants, only: MAX_STRING_LEN,IOUT_VTK

  implicit none

  integer :: nspec,NGLL

  ! global coordinates
  integer, dimension(NGLL,NGLL,NGLL,nspec) :: ibool
  double precision, dimension(NGLL,NGLL,NGLL,nspec) :: xstore_mesh,ystore_mesh,zstore_mesh

  ! element flag array
  integer, dimension(nspec) :: elem_data

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: ispec,i,j,k

  open(IOUT_VTK,file=trim(filename),status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nspec*NGLL*NGLL*NGLL, ' float'
  do ispec = 1,nspec
    do k = 1,NGLL
      do j = 1,NGLL
        do i = 1,NGLL
          write(IOUT_VTK,*) xstore_mesh(i,j,k,ispec),ystore_mesh(i,j,k,ispec),zstore_mesh(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
                             ibool(1,1,1,ispec)-1,ibool(NGLL,1,1,ispec)-1, &
                             ibool(NGLL,NGLL,1,ispec)-1,ibool(1,NGLL,1,ispec)-1, &
                             ibool(1,1,NGLL,ispec)-1,ibool(NGLL,1,NGLL,ispec)-1, &
                             ibool(NGLL,NGLL,NGLL,ispec)-1,ibool(1,NGLL,NGLL,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_val float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_data(ispec)
  enddo
  write(IOUT_VTK,*) ''

  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_i_meshfem
