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

module visualize_constants

implicit none
!----------------------------------------------
! constants section (imported from constants.h)
integer, parameter :: SIZE_REAL = 4
integer, parameter :: SIZE_DOUBLE = 8

! usually the size of integer and logical variables is the same as regular single-precision real variable
integer, parameter :: SIZE_INTEGER = SIZE_REAL
integer, parameter :: SIZE_LOGICAL = SIZE_REAL

! set to SIZE_REAL to run in single precision
! set to SIZE_DOUBLE to run in double precision (increases memory size by 2)
integer, parameter :: CUSTOM_REAL = SIZE_REAL

! number of GLL points in each direction of an element (degree plus one)
integer, parameter :: NGLLX = 5
integer, parameter :: NGLLY = NGLLX
integer, parameter :: NGLLZ = NGLLX
integer,parameter :: NGLLX_MID = 3 ! = NGLLX/2 + 1
integer,parameter :: NGLLY_MID = NGLLX_MID
integer,parameter :: NGLLZ_MID = NGLLX_MID

! 3-D simulation
integer, parameter :: NDIM = 3
integer, parameter :: NGNOD = 8
!----------------------------------------------
  
integer, parameter :: SIZE_INT = 4
integer, parameter :: SIZE_FLOAT = 4
character(len=1),parameter :: CR = achar(13) ! Carriage-return for overwriting
integer,save :: NENOD_OUT ! number of elemental nodes (nodes per element) for output

end module visualize_constants
!=====================================================================

module visualize_par
use visualize_constants

implicit none  

! mesh coordinates
real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
integer, dimension(:,:,:,:),allocatable :: ibool

integer :: NSPEC_AB, NGLOB_AB    

! input information
integer :: inp_ncomp
integer :: t_nstep,t_width,t_inc,t_start
double precision :: dt
character(len=80),dimension(:),allocatable :: inp_head
character(len=20) :: inp_ext,t_form,t_width_str
character(len=80) :: inp_path

! output information
integer :: out_format,out_ncomp,out_res,out_nslice
character(len=20) :: out_ext
character(len=80) :: out_path,out_fname,out_head,out_vname    

! File descriptors
integer :: fd,fd_con,fd_x,fd_y,fd_z
integer,dimension(:),allocatable :: fd_array
integer :: dat_topo

! processor and slice infromation
integer :: proc_width
character(len=20) :: proc_form,proc_width_str
character(len=60) :: proc_head  
integer,dimension(:),allocatable :: slice_nnode,slice_nelmt ! Number of nodes and element in output slice 
integer,allocatable,dimension(:,:) :: slice_proc_list ! Procesor list in ouput slice
integer,allocatable,dimension(:) :: slice_nproc ! Number of input procesors in output slice
character(len=60),allocatable,dimension(:) :: server_name, server_exec

character(len=80) :: format_str1,format_str2,format_str3,num_str1,num_str2,tmp_str     
real :: tmp_real  
integer :: tmp_int

end module visualize_par
!=====================================================================
  
