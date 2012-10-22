!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
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

! utility to locate partition which is closest to given point location
!
! compile with one of these (use your default):
! > gfortran -I src/shared -o bin/xlocate_partition utils/locate_partition.f90
! > ifort -assume byterecl -I src/shared -o bin/xlocate_partition utils/locate_partition.f90
! > pgf90 -I src/shared -o bin/xlocate_partition utils/locate_partition.f90
!
! specify a target (x,y) may be in UTM, not lon-lat, then run with:
! > ./bin/xlocate_partition 70000.0 11000.0 -3000.0 ./OUTPUT_FILES/DATABASES_MPI/
!
! this will generate the output file OUTPUT_FILES/DATABASES_MPI/partition_bounds.dat

  program locate_partition

! works for external, unregular meshes

  implicit none

  include 'constants.h'

  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  integer, dimension(:,:,:,:),allocatable :: ibool
  integer :: NSPEC_AB, NGLOB_AB

  integer :: i,ios,ier
  integer :: iproc
  character(len=256) :: arg(4)
  character(len=256) :: LOCAL_PATH
  character(len=256) :: prname_lp

  real(kind=CUSTOM_REAL) :: x_found,y_found,z_found,distance
  real(kind=CUSTOM_REAL) :: target_x,target_y,target_z
  real(kind=CUSTOM_REAL) :: total_x,total_y,total_z
  real(kind=CUSTOM_REAL) :: total_distance
  integer :: total_partition

! checks given arguments
  print *
  print *,'locate partition'
  print *,'----------------------------'

  do i = 1, 4
    call getarg(i,arg(i))
    if (i <= 4 .and. trim(arg(i)) == '') then
      print *, 'Usage: '
      print *, '        xlocate_partition x y z Databases_directory'
      stop ' Reenter command line options'
    endif
  enddo
  read(arg(1),*) target_x
  read(arg(2),*) target_y
  read(arg(3),*) target_z
  LOCAL_PATH = arg(4)

  print *,'search location: '
  print *,'  x = ',target_x
  print *,'  y = ',target_y
  print *,'  z = ',target_z
  print *,'in directory: ',trim(LOCAL_PATH)
  print *,'----------------------------'
  print *

  ! open a text file to list the maximal bounds of each partition
  open(11,file=trim(LOCAL_PATH)//'partition_bounds.dat',status='unknown')

  ! loops over slices (process partitions)
  total_distance = HUGEVAL
  total_partition = -1
  total_x = 0.0
  total_y = 0.0
  total_z = 0.0
  iproc = -1
  do while( iproc < 10000000 )
    ! starts with 0
    iproc = iproc + 1

    ! gets number of elements and global points for this partition
    write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'
    open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin',&
          status='old',action='read',form='unformatted',iostat=ios)
    if( ios /= 0 ) exit

    read(27,iostat=ier) NSPEC_AB
    if( ier /= 0 ) stop 'please check your compilation, use the same compiler & flags as for SPECFEM3D'
    read(27,iostat=ier) NGLOB_AB
    if( ier /= 0 ) stop 'please check your compilation, use the same compiler & flags as for SPECFEM3D'

    ! ibool file
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibool'
    read(27,iostat=ier) ibool
    if( ier /= 0 ) stop 'please check your compilation, use the same compiler & flags as for SPECFEM3D'

    ! global point arrays
    allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array xstore etc.'
    read(27,iostat=ier) xstore
    if( ier /= 0 ) stop 'please check your compilation, use the same compiler & flags as for SPECFEM3D'
    read(27,iostat=ier) ystore
    if( ier /= 0 ) stop 'please check your compilation, use the same compiler & flags as for SPECFEM3D'
    read(27,iostat=ier) zstore
    if( ier /= 0 ) stop 'please check your compilation, use the same compiler & flags as for SPECFEM3D'
    close(27)

    print*,'partition: ',iproc
    print*,'  min/max x = ',minval(xstore),maxval(xstore)
    print*,'  min/max y = ',minval(ystore),maxval(ystore)
    print*,'  min/max z = ',minval(zstore),maxval(zstore)
    print*

    write(11,'(i10,6e18.6)') iproc,minval(xstore),maxval(xstore),minval(ystore),maxval(ystore),minval(zstore),maxval(zstore)

    ! gets distance to target location
    call get_closest_point(target_x,target_y,target_z, &
                         NGLOB_AB,NSPEC_AB,xstore,ystore,zstore,ibool, &
                         distance,x_found,y_found,z_found)

    if( distance < total_distance ) then
      total_distance = distance
      total_partition = iproc
      total_x = x_found
      total_y = y_found
      total_z = z_found
    endif

    ! cleans up memory allocations
    deallocate(ibool,xstore,ystore,zstore)

  enddo  ! all slices for points

  close(11)

  ! checks
  if (total_partition < 0 ) then
    print*,'Error: partition not found among ',iproc,'partitions searched'
    stop 'Error: partition not found'
  endif

  ! output
  print*,'number of partitions searched: ',iproc
  print*
  print*,'closest grid point location found:'
  print*,'  x = ',total_x
  print*,'  y = ',total_y
  print*,'  z = ',total_z
  print*,'  distance to search location = ',sqrt(total_distance)
  print*,'closest partition: '
  print*,'  partition = ',total_partition
  print*

  end program locate_partition

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_closest_point(target_x,target_y,target_z, &
                         NGLOB_AB,NSPEC_AB,xstore,ystore,zstore,ibool, &
                         distance,x_found,y_found,z_found)

  implicit none
  include 'constants.h'

  real(kind=CUSTOM_REAL),intent(in) :: target_x,target_y,target_z

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: xstore, ystore, zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  real(kind=CUSTOM_REAL),intent(out) :: distance
  real(kind=CUSTOM_REAL),intent(out) :: x_found,y_found,z_found

  ! local parameters
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: dist

  distance = HUGEVAL
  x_found = 0.0
  y_found = 0.0
  z_found = 0.0

  do ispec=1,NSPEC_AB
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dist =  (target_x - xstore(iglob))*(target_x - xstore(iglob)) &
                + (target_y - ystore(iglob))*(target_y - ystore(iglob)) &
                + (target_z - zstore(iglob))*(target_z - zstore(iglob))

          if( dist < distance ) then
            distance = dist
            x_found = xstore(iglob)
            y_found = ystore(iglob)
            z_found = zstore(iglob)
          endif
        enddo
      enddo
    enddo
  enddo

  end subroutine
