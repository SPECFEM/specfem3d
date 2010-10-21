!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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
! Last revision: Hom Nath Gharti, April 06,2010

subroutine cvd_count_totals_ext_mesh(nproc,proc_list,proc_width,inp_dir,&
                        node_count,elmt_count,out_res)

! counts total number of points and elements for external meshes in given slice list
! returns: total number of elements (elmt_count) and number of points (node_count)
use visualize_constants
!use visualize_par
implicit none
!include 'constants.h'

integer,intent(in) :: nproc, proc_width
integer,dimension(nproc),intent(in) :: proc_list
character(len=60),intent(in) :: inp_dir
integer,intent(out) :: node_count,elmt_count
integer,intent(in) :: out_res

! local parameters
integer, dimension(:,:,:,:),allocatable :: ibool
logical, dimension(:),allocatable :: mask_ibool
integer :: NSPEC_AB, NGLOB_AB
integer :: iproc,npoint,nelement,ios,i_spec
integer,dimension(NENOD_OUT) :: iglob
character(len=256) :: mesh_fname
character(len=20) :: tmp_str,proc_width_str,format_str
integer :: i_proc

! loops over all slices (process partitions)
node_count = 0
elmt_count = 0

write(tmp_str,*)proc_width
write(proc_width_str,*)trim(adjustl(tmp_str))//'.'//trim(adjustl(tmp_str)) 
format_str='(a,i'//trim(adjustl(proc_width_str))//',a)'

do i_proc = 1, nproc

  ! gets number of elements and points for this slice
  iproc = proc_list(i_proc)
  
  write(mesh_fname,fmt=format_str) trim(inp_dir)//'/proc',iproc,'_external_mesh.bin'    
  open(unit=27,file=trim(mesh_fname),status='old',action='read',form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(*,'(/,a)')'ERROR: file '//trim(mesh_fname)//' cannot be opend!'
    stop
  endif

  read(27) NSPEC_AB
  read(27) NGLOB_AB 
  ! gets ibool
  if(out_res/=2) then
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    read(27) ibool
  endif
  close(27)   
      
  ! calculates totals
  if(out_res==2) then
    ! total number of global points
    node_count = node_count + NGLOB_AB

    ! total number of elements
    ! each spectral elements gets subdivided by GLL points, 
    ! which form (NGLLX-1)*(NGLLY-1)*(NGLLZ-1) sub-elements
    nelement = NSPEC_AB * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1) 
    elmt_count = elmt_count + nelement

  elseif (out_res==1)then ! Medium resolution

    ! mark element corners (global AVS or DX points)
    allocate(mask_ibool(NGLOB_AB))      
    mask_ibool = .false.
    do i_spec=1,NSPEC_AB
      ! Bottom corners
      iglob(1)=ibool(1,1,1,i_spec)
      iglob(2)=ibool(NGLLX,1,1,i_spec)
      iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
      iglob(4)=ibool(1,NGLLY,1,i_spec)

      ! Top corners
      iglob(5)=ibool(1,1,NGLLZ,i_spec)
      iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
      iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
      iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)

      ! Bottom mid points
      iglob(9)=ibool(NGLLX_MID,1,1,i_spec)
      iglob(10)=ibool(NGLLX,NGLLY_MID,1,i_spec)
      iglob(11)=ibool(NGLLX_MID,NGLLY,1,i_spec)
      iglob(12)=ibool(1,NGLLY_MID,1,i_spec)

      ! Top mid points
      iglob(13)=ibool(NGLLX_MID,1,NGLLZ,i_spec)
      iglob(14)=ibool(NGLLX,NGLLY_MID,NGLLZ,i_spec)
      iglob(15)=ibool(NGLLX_MID,NGLLY,NGLLZ,i_spec)
      iglob(16)=ibool(1,NGLLY_MID,NGLLZ,i_spec)

      ! Center mid points
      iglob(17)=ibool(1,1,NGLLZ_MID,i_spec)
      iglob(18)=ibool(NGLLX,1,NGLLZ_MID,i_spec)
      iglob(19)=ibool(NGLLX,NGLLY,NGLLZ_MID,i_spec)
      iglob(20)=ibool(1,NGLLY,NGLLZ_MID,i_spec)

      mask_ibool(iglob) = .true.
      !mask_ibool(iglob2) = .true.
      !mask_ibool(iglob3) = .true.
      !mask_ibool(iglob4) = .true.
      !mask_ibool(iglob5) = .true.
      !mask_ibool(iglob6) = .true.
      !mask_ibool(iglob7) = .true.
      !mask_ibool(iglob8) = .true.
      !mask_ibool(iglob9) = .true.
      !mask_ibool(iglob10) = .true.
      !mask_ibool(iglob11) = .true.
      !mask_ibool(iglob12) = .true.
      !mask_ibool(iglob13) = .true.
      !mask_ibool(iglob14) = .true.
      !mask_ibool(iglob15) = .true.
      !mask_ibool(iglob16) = .true.
      !mask_ibool(iglob17) = .true.
      !mask_ibool(iglob18) = .true.
      !mask_ibool(iglob19) = .true.
      !mask_ibool(iglob20) = .true.
    enddo        

    ! count global number of AVS or DX points
    npoint = count(mask_ibool(:))      

    node_count = node_count + npoint
    
    ! total number of spectral elements
    elmt_count = elmt_count + NSPEC_AB
    deallocate(mask_ibool)
  elseif (out_res==0)then ! 
    ! mark element corners (global AVS or DX points)
    allocate(mask_ibool(NGLOB_AB))      
    mask_ibool = .false.
    do i_spec=1,NSPEC_AB
      iglob(1)=ibool(1,1,1,i_spec)
      iglob(2)=ibool(NGLLX,1,1,i_spec)
      iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
      iglob(4)=ibool(1,NGLLY,1,i_spec)
      iglob(5)=ibool(1,1,NGLLZ,i_spec)
      iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
      iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
      iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)
      mask_ibool(iglob) = .true.
      !mask_ibool(iglob2) = .true.
      !mask_ibool(iglob3) = .true.
      !mask_ibool(iglob4) = .true.
      !mask_ibool(iglob5) = .true.
      !mask_ibool(iglob6) = .true.
      !mask_ibool(iglob7) = .true.
      !mask_ibool(iglob8) = .true.
    enddo        

    ! count global number of AVS or DX points
    npoint = count(mask_ibool(:))      

    node_count = node_count + npoint
    
    ! total number of spectral elements
    elmt_count = elmt_count + NSPEC_AB
    deallocate(mask_ibool)
  else
    write(*,'(/,a)')'ERROR: wrong out_res value!'
    stop
  endif ! out_res
  if(out_res/=2) then
    deallocate(ibool)
  endif
        
enddo
  
end subroutine cvd_count_totals_ext_mesh  
!=============================================================


subroutine cvd_write_corners_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,&
                              numpoin,fd_x,fd_y,fd_z)

! writes out locations of spectral element corners only
use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore, ystore, zstore


integer:: i_node,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
real :: x, y, z
integer :: i_spec
integer,dimension(NENOD_OUT) :: iglob
integer :: fd_x, fd_y, fd_z

! writes our corner point locations  
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)
  
  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      x = xstore(iglob(i_node))
      y = ystore(iglob(i_node))
      z = zstore(iglob(i_node))
      call write_real(x,fd_x)
      call write_real(y,fd_y)
      call write_real(z,fd_z)
    
      mask_ibool(iglob(i_node)) = .true.
    endif
  enddo

enddo ! i_spec
  !stop 
end subroutine cvd_write_corners_only
!=============================================================

subroutine cvd_write_hexa20_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,&
                              numpoin,fd_x,fd_y,fd_z)

! writes out locations of spectral element corners only
use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore, ystore, zstore

integer:: i_node,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
real :: x, y, z
integer :: i_spec
integer,dimension(NENOD_OUT) :: iglob
integer :: fd_x, fd_y, fd_z  

! writes our corner point locations  
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  ! Bottom corners
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)

  ! Top corners
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)

  ! Bottom mid points
  iglob(9)=ibool(NGLLX_MID,1,1,i_spec)
  iglob(10)=ibool(NGLLX,NGLLY_MID,1,i_spec)
  iglob(11)=ibool(NGLLX_MID,NGLLY,1,i_spec)
  iglob(12)=ibool(1,NGLLY_MID,1,i_spec)

  ! Top mid points
  iglob(13)=ibool(NGLLX_MID,1,NGLLZ,i_spec)
  iglob(14)=ibool(NGLLX,NGLLY_MID,NGLLZ,i_spec)
  iglob(15)=ibool(NGLLX_MID,NGLLY,NGLLZ,i_spec)
  iglob(16)=ibool(1,NGLLY_MID,NGLLZ,i_spec)

  ! Center mid points
  iglob(17)=ibool(1,1,NGLLZ_MID,i_spec)
  iglob(18)=ibool(NGLLX,1,NGLLZ_MID,i_spec)
  iglob(19)=ibool(NGLLX,NGLLY,NGLLZ_MID,i_spec)
  iglob(20)=ibool(1,NGLLY,NGLLZ_MID,i_spec)

  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      x = xstore(iglob(i_node))
      y = ystore(iglob(i_node))
      z = zstore(iglob(i_node))
      call write_real(x,fd_x)
      call write_real(y,fd_y)
      call write_real(z,fd_z)
    
      mask_ibool(iglob(i_node)) = .true.
    endif
  enddo   

enddo ! i_spec
  !stop 
end subroutine cvd_write_hexa20_only
!=============================================================

subroutine cvd_write_GLL_points_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,&
                                numpoin,fd_x,fd_y,fd_z)

! writes out locations of all GLL points of spectral elements

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore, ystore, zstore

integer:: numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
real :: x, y, z
integer :: i_spec,i,j,k,iglob1
integer :: fd_x, fd_y, fd_z

! writes out point locations and values
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        iglob1 = ibool(i,j,k,i_spec)
        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          x = xstore(iglob1)
          y = ystore(iglob1)
          z = zstore(iglob1)
          call write_real(x,fd_x)
          call write_real(y,fd_y)
          call write_real(z,fd_z)            
          mask_ibool(iglob1) = .true.
        endif
      enddo ! i
    enddo ! j
  enddo ! k
enddo !i_spec

end subroutine cvd_write_GLL_points_only
!=============================================================

subroutine cvd_write_corners_data(NSPEC_AB,NGLOB_AB,ibool,dat,&
                              numpoin,fd)

! writes out locations of spectral element corners only

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool  
real,dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: dat
integer:: i_node,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
!real :: x, y, z
integer :: i_spec
integer,dimension(NENOD_OUT) :: iglob
real,dimension(NENOD_OUT) :: tmp_dat
integer :: fd

! writes out corner point locations  
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)
  
  ! Nodal data in a element
  tmp_dat(1)=dat(1,1,1,i_spec)
  tmp_dat(2)=dat(NGLLX,1,1,i_spec)
  tmp_dat(3)=dat(NGLLX,NGLLY,1,i_spec)
  tmp_dat(4)=dat(1,NGLLY,1,i_spec)
  tmp_dat(5)=dat(1,1,NGLLZ,i_spec)
  tmp_dat(6)=dat(NGLLX,1,NGLLZ,i_spec)
  tmp_dat(7)=dat(NGLLX,NGLLY,NGLLZ,i_spec)
  tmp_dat(8)=dat(1,NGLLY,NGLLZ,i_spec)
  
  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      call write_real(tmp_dat(i_node),fd)
      mask_ibool(iglob(i_node)) = .true.
    endif
  enddo
  
enddo ! i_spec
  
end subroutine cvd_write_corners_data
!=============================================================

subroutine cvd_write_corners_data_glob(NSPEC_AB,NGLOB_AB,ibool,dat,&
                              numpoin,fd)

! writes out locations of spectral element corners only

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool  
real,dimension(NGLOB_AB),intent(in) :: dat
integer:: i_node,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
!real :: x, y, z
integer :: i_spec
integer,dimension(NENOD_OUT) :: iglob
integer :: fd

! writes out corner point locations  
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)
  
  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      call write_real(dat(iglob(i_node)),fd)
      mask_ibool(iglob(i_node)) = .true.
    endif
  enddo

enddo ! i_spec
  
end subroutine cvd_write_corners_data_glob
!=============================================================

subroutine cvd_write_hexa20_data(NSPEC_AB,NGLOB_AB,ibool,dat,&
                              numpoin,fd)

! writes out locations of spectral element corners only

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool  
real,dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: dat
integer:: i_node,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
!real :: x, y, z
integer :: i_spec
integer,dimension(NENOD_OUT) :: iglob
real,dimension(NENOD_OUT) :: tmp_dat
integer :: fd

! writes out corner point locations  
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  ! Bottom corners
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)

  ! Top corners
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)

  ! Bottom mid points
  iglob(9)=ibool(NGLLX_MID,1,1,i_spec)
  iglob(10)=ibool(NGLLX,NGLLY_MID,1,i_spec)
  iglob(11)=ibool(NGLLX_MID,NGLLY,1,i_spec)
  iglob(12)=ibool(1,NGLLY_MID,1,i_spec)

  ! Top mid points
  iglob(13)=ibool(NGLLX_MID,1,NGLLZ,i_spec)
  iglob(14)=ibool(NGLLX,NGLLY_MID,NGLLZ,i_spec)
  iglob(15)=ibool(NGLLX_MID,NGLLY,NGLLZ,i_spec)
  iglob(16)=ibool(1,NGLLY_MID,NGLLZ,i_spec)

  ! Center mid points
  iglob(17)=ibool(1,1,NGLLZ_MID,i_spec)
  iglob(18)=ibool(NGLLX,1,NGLLZ_MID,i_spec)
  iglob(19)=ibool(NGLLX,NGLLY,NGLLZ_MID,i_spec)
  iglob(20)=ibool(1,NGLLY,NGLLZ_MID,i_spec)
  
  ! Nodal data in a element
  ! Bottom corners
  tmp_dat(1)=dat(1,1,1,i_spec)
  tmp_dat(2)=dat(NGLLX,1,1,i_spec)
  tmp_dat(3)=dat(NGLLX,NGLLY,1,i_spec)
  tmp_dat(4)=dat(1,NGLLY,1,i_spec)

  ! Top corners
  tmp_dat(5)=dat(1,1,NGLLZ,i_spec)
  tmp_dat(6)=dat(NGLLX,1,NGLLZ,i_spec)
  tmp_dat(7)=dat(NGLLX,NGLLY,NGLLZ,i_spec)
  tmp_dat(8)=dat(1,NGLLY,NGLLZ,i_spec)

  ! Bottom mid points
  tmp_dat(9)=dat(NGLLX_MID,1,1,i_spec)
  tmp_dat(10)=dat(NGLLX,NGLLY_MID,1,i_spec)
  tmp_dat(11)=dat(NGLLX_MID,NGLLY,1,i_spec)
  tmp_dat(12)=dat(1,NGLLY_MID,1,i_spec)

  ! Top mid points
  tmp_dat(13)=dat(NGLLX_MID,1,NGLLZ,i_spec)
  tmp_dat(14)=dat(NGLLX,NGLLY_MID,NGLLZ,i_spec)
  tmp_dat(15)=dat(NGLLX_MID,NGLLY,NGLLZ,i_spec)
  tmp_dat(16)=dat(1,NGLLY_MID,NGLLZ,i_spec)

  ! Center mid points
  tmp_dat(17)=dat(1,1,NGLLZ_MID,i_spec)
  tmp_dat(18)=dat(NGLLX,1,NGLLZ_MID,i_spec)
  tmp_dat(19)=dat(NGLLX,NGLLY,NGLLZ_MID,i_spec)
  tmp_dat(20)=dat(1,NGLLY,NGLLZ_MID,i_spec)
  
  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      call write_real(tmp_dat(i_node),fd)
      mask_ibool(iglob(i_node)) = .true.
    endif
  enddo    
    
enddo ! i_spec
  
end subroutine cvd_write_hexa20_data
!=============================================================

subroutine cvd_write_hexa20_data_glob(NSPEC_AB,NGLOB_AB,ibool,dat,&
                              numpoin,fd)

! writes out locations of spectral element corners only

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool  
real,dimension(NGLOB_AB),intent(in) :: dat
integer:: i_node,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
!real :: x, y, z
integer :: i_spec
integer,dimension(NENOD_OUT) :: iglob
integer :: fd

! writes out corner point locations  
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  ! Bottom corners
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)

  ! Top corners
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)

  ! Bottom mid points
  iglob(9)=ibool(NGLLX_MID,1,1,i_spec)
  iglob(10)=ibool(NGLLX,NGLLY_MID,1,i_spec)
  iglob(11)=ibool(NGLLX_MID,NGLLY,1,i_spec)
  iglob(12)=ibool(1,NGLLY_MID,1,i_spec)

  ! Top mid points
  iglob(13)=ibool(NGLLX_MID,1,NGLLZ,i_spec)
  iglob(14)=ibool(NGLLX,NGLLY_MID,NGLLZ,i_spec)
  iglob(15)=ibool(NGLLX_MID,NGLLY,NGLLZ,i_spec)
  iglob(16)=ibool(1,NGLLY_MID,NGLLZ,i_spec)

  ! Center mid points
  iglob(17)=ibool(1,1,NGLLZ_MID,i_spec)
  iglob(18)=ibool(NGLLX,1,NGLLZ_MID,i_spec)
  iglob(19)=ibool(NGLLX,NGLLY,NGLLZ_MID,i_spec)
  iglob(20)=ibool(1,NGLLY,NGLLZ_MID,i_spec)

  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      call write_real(dat(iglob(i_node)),fd)
      mask_ibool(iglob(i_node)) = .true.
    endif
  enddo    

enddo ! i_spec
  
end subroutine cvd_write_hexa20_data_glob
!=============================================================  

subroutine cvd_write_GLL_points_data(NSPEC_AB,NGLOB_AB,ibool,dat,&
                                numpoin,fd)

! writes out locations of all GLL points of spectral elements

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool  
real,dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: dat
integer:: numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
!real :: x, y, z
integer :: i_spec,i,j,k,iglob1
integer :: fd

! writes out point locations and values
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        iglob1 = ibool(i,j,k,i_spec)
        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1            
          call write_real(dat(i,j,k,i_spec),fd)
          mask_ibool(iglob1) = .true.
        endif
      enddo ! i
    enddo ! j
  enddo ! k
enddo !i_spec

end subroutine cvd_write_GLL_points_data  
!=============================================================

subroutine cvd_write_GLL_points_data_glob(NSPEC_AB,NGLOB_AB,ibool,dat,&
                                numpoin,fd)

! writes out locations of all GLL points of spectral elements

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool  
real,dimension(NGLOB_AB),intent(in) :: dat
integer:: numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
!real :: x, y, z
integer :: i_spec,i,j,k,iglob1
integer :: fd

! writes out point locations and values
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        iglob1 = ibool(i,j,k,i_spec)
        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1            
          call write_real(dat(iglob1),fd)
          mask_ibool(iglob1) = .true.
        endif
      enddo ! i
    enddo ! j
  enddo ! k
enddo !i_spec

end subroutine cvd_write_GLL_points_data_glob  
!=============================================================

! writes out locations of spectral element corners only

subroutine cvd_write_corner_elements(NSPEC_AB,NGLOB_AB,ibool,&
                                    np,nelement,numpoin,fd)

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
integer:: i_node,np,nelement,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
integer,dimension(:),allocatable :: num_ibool  
integer :: i_spec 
integer,dimension(NENOD_OUT) :: iglob
integer :: inode
integer :: fd

! writes out element indices
allocate(mask_ibool(NGLOB_AB))
allocate(num_ibool(NGLOB_AB))
mask_ibool(:) = .false.
num_ibool(:) = 0
numpoin = 0    
do i_spec=1,NSPEC_AB  
  ! gets corner indices
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)

  ! sets increasing numbering
  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      num_ibool(iglob(i_node)) = numpoin
      mask_ibool(iglob(i_node)) = .true.          
    endif
  enddo
  do i_node=1,NENOD_OUT
    inode = num_ibool(iglob(i_node)) + np ! -1
    call write_integer(inode,fd)
  enddo
enddo

! elements written
nelement = NSPEC_AB

! updates points written
np = np + numpoin
  
end subroutine cvd_write_corner_elements  
!=============================================================

! writes out locations of spectral element corners only

subroutine cvd_write_hexa20_elements(NSPEC_AB,NGLOB_AB,ibool,&
                                    np,nelement,numpoin,fd)

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
integer:: i_node,np,nelement,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
integer,dimension(:),allocatable :: num_ibool  
integer :: i_spec 
integer,dimension(NENOD_OUT) :: iglob
integer :: inode  
integer :: fd

! writes out element indices
allocate(mask_ibool(NGLOB_AB))
allocate(num_ibool(NGLOB_AB))
mask_ibool(:) = .false.
num_ibool(:) = 0
numpoin = 0    
do i_spec=1,NSPEC_AB  
  ! Bottom corners
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)

  ! Top corners
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)

  ! Bottom mid points
  iglob(9)=ibool(NGLLX_MID,1,1,i_spec)
  iglob(10)=ibool(NGLLX,NGLLY_MID,1,i_spec)
  iglob(11)=ibool(NGLLX_MID,NGLLY,1,i_spec)
  iglob(12)=ibool(1,NGLLY_MID,1,i_spec)

  ! Top mid points
  iglob(13)=ibool(NGLLX_MID,1,NGLLZ,i_spec)
  iglob(14)=ibool(NGLLX,NGLLY_MID,NGLLZ,i_spec)
  iglob(15)=ibool(NGLLX_MID,NGLLY,NGLLZ,i_spec)
  iglob(16)=ibool(1,NGLLY_MID,NGLLZ,i_spec)

  ! Center mid points
  iglob(17)=ibool(1,1,NGLLZ_MID,i_spec)
  iglob(18)=ibool(NGLLX,1,NGLLZ_MID,i_spec)
  iglob(19)=ibool(NGLLX,NGLLY,NGLLZ_MID,i_spec)
  iglob(20)=ibool(1,NGLLY,NGLLZ_MID,i_spec)

  ! sets increasing numbering
  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      num_ibool(iglob(i_node)) = numpoin
      mask_ibool(iglob(i_node)) = .true.          
    endif
  enddo
  do i_node=1,NENOD_OUT
    inode = num_ibool(iglob(i_node)) + np ! -1
    call write_integer(inode,fd)
  enddo
  
enddo

! elements written
nelement = NSPEC_AB

! updates points written
np = np + numpoin
  
end subroutine cvd_write_hexa20_elements  
!=============================================================


subroutine cvd_write_GLL_elements(NSPEC_AB,NGLOB_AB,ibool, &
                                  np,nelement,numpoin,fd)

! writes out indices of elements given by GLL points 

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in):: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
integer:: i_node,np,numpoin,nelement

! local parameters
logical,dimension(:),allocatable :: mask_ibool
integer,dimension(:),allocatable :: num_ibool    
integer :: i_spec,i,j,k
integer,dimension(NENOD_OUT) :: iglob
integer :: iglob1,inode
integer :: fd

! sets numbering num_ibool respecting mask
allocate(mask_ibool(NGLOB_AB))
allocate(num_ibool(NGLOB_AB))
mask_ibool(:) = .false.
num_ibool(:) = 0
numpoin = 0  
do i_spec=1,NSPEC_AB
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        iglob1 = ibool(i,j,k,i_spec)
        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          num_ibool(iglob1) = numpoin
          mask_ibool(iglob1) = .true.
        endif
      enddo ! i
    enddo ! j
  enddo ! k
enddo !i_spec

! outputs GLL subelement
do i_spec = 1, NSPEC_AB
  do k = 1, NGLLZ-1
    do j = 1, NGLLY-1
      do i = 1, NGLLX-1
        iglob(1) = ibool(i,j,k,i_spec)
        iglob(2) = ibool(i+1,j,k,i_spec)
        iglob(3) = ibool(i+1,j+1,k,i_spec)
        iglob(4) = ibool(i,j+1,k,i_spec)
        iglob(5) = ibool(i,j,k+1,i_spec)
        iglob(6) = ibool(i+1,j,k+1,i_spec)
        iglob(7) = ibool(i+1,j+1,k+1,i_spec)
        iglob(8) = ibool(i,j+1,k+1,i_spec)
        do i_node=1,NENOD_OUT
          inode = num_ibool(iglob(i_node)) + np ! -1
          call write_integer(inode,fd)
        enddo
      enddo
    enddo
  enddo
enddo
! elements written
nelement = NSPEC_AB * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1) 

! updates points written
np = np + numpoin

end subroutine cvd_write_GLL_elements
!=============================================================

! subroutines below may not be used

subroutine cvd_write_corners(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,dat,&
                              i_proc,npp,numpoin,fd)

! writes out locations of spectral element corners only

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore, ystore, zstore
real,dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: dat
integer:: i_node,i_proc
integer :: npp,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
real :: x, y, z
integer :: i_spec
integer,dimension(NENOD_OUT) :: iglob
real,dimension(NENOD_OUT) :: tmp_dat
integer :: fd

!writes out total number of points
if (i_proc == 1) then
  call write_integer(npp)
endif

! writes our corner point locations  
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  iglob(1)=ibool(1,1,1,i_spec)
  iglob(2)=ibool(NGLLX,1,1,i_spec)
  iglob(3)=ibool(NGLLX,NGLLY,1,i_spec)
  iglob(4)=ibool(1,NGLLY,1,i_spec)
  iglob(5)=ibool(1,1,NGLLZ,i_spec)
  iglob(6)=ibool(NGLLX,1,NGLLZ,i_spec)
  iglob(7)=ibool(NGLLX,NGLLY,NGLLZ,i_spec)
  iglob(8)=ibool(1,NGLLY,NGLLZ,i_spec)
  
  tmp_dat(1)=dat(1,1,1,i_spec)
  tmp_dat(2)=dat(NGLLX,1,1,i_spec)
  tmp_dat(3)=dat(NGLLX,NGLLY,1,i_spec)
  tmp_dat(4)=dat(1,NGLLY,1,i_spec)
  tmp_dat(5)=dat(1,1,NGLLZ,i_spec)
  tmp_dat(6)=dat(NGLLX,1,NGLLZ,i_spec)
  tmp_dat(7)=dat(NGLLX,NGLLY,NGLLZ,i_spec)
  tmp_dat(8)=dat(1,NGLLY,NGLLZ,i_spec)
  
  do i_node=1,NENOD_OUT
    if(.not. mask_ibool(iglob(i_node))) then
      numpoin = numpoin + 1
      x = xstore(iglob(i_node))
      y = ystore(iglob(i_node))
      z = zstore(iglob(i_node))
      call write_real(x,fd)
      call write_real(y,fd)
      call write_real(z,fd)
      call write_real(tmp_dat(i_node),fd)
      mask_ibool(iglob(i_node)) = .true.
    endif
  enddo

enddo ! i_spec
  
end subroutine cvd_write_corners
!=============================================================

subroutine cvd_write_GLL_points(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,dat,&
                                i_proc,npp,numpoin,fd)

! writes out locations of all GLL points of spectral elements

use visualize_constants
implicit none
!include 'constants.h'

integer,intent(in) :: NSPEC_AB,NGLOB_AB
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore, ystore, zstore
real,dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: dat
integer:: i_proc,npp,numpoin

! local parameters
logical,dimension(:),allocatable :: mask_ibool
real :: x, y, z
integer :: i_spec,i,j,k,iglob1
integer :: fd

!writes out total number of points
if (i_proc == 1) then
  call write_integer(npp)
endif

! writes out point locations and values
allocate(mask_ibool(NGLOB_AB))
mask_ibool(:) = .false.
numpoin = 0
do i_spec=1,NSPEC_AB
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        iglob1 = ibool(i,j,k,i_spec)
        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          x = xstore(iglob1)
          y = ystore(iglob1)
          z = zstore(iglob1)
          call write_real(x,fd)
          call write_real(y,fd)
          call write_real(z,fd)
          call write_real(dat(i,j,k,i_spec),fd)
          mask_ibool(iglob1) = .true.
        endif
      enddo ! i
    enddo ! j
  enddo ! k
enddo !i_spec

end subroutine cvd_write_GLL_points
!=============================================================
  
