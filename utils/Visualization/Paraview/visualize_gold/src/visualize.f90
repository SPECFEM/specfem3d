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
program visualize

! This programs writes Ensight Gold binary file collecting binary mesh and volume data
! files produced by SPECFEM3D. The Ensight Gold file can then be visualized in 
! VTK/ParaView. See http://www.vtk.org and http://www.paraview.org for details.
!------------------------------------------
! DEPENDENCY:
!   cfunc4fortran.c, visualize_par.f90, visualize_collect.f90, write_ensight.f90,
!   write_vtu.f90 
! COMPILE
!   >> make
! USAGE
!   ./xvisualize [input filename]
!   e.g., ./xvisualize visualize.in
!   see visualize.in for input detail
! HISTORY:
!   Hom Nath Gharti, NORSAR
!   Mar 10,2010 (NORSAR)
! FEEDBACK:
!   homnath_AT_norsar_DOT_no
!------------------------------------------  
  
use visualize_par
implicit none

integer :: i_slice
character(len=80) :: inp_fname  
  
!write(*,*)iargc()
if (iargc() <= 0) then
  write(*,'(/,a)')'ERROR: no input file!'
  stop! counts total number of nodes and elementsp
endif
  
write(*,'(a)',advance='no')'reading main input file...'
call getarg(1,inp_fname)
   
call read_input(inp_fname)
write(*,'(a)')'complete!'
  
! check and display key input parameters...'
if (out_res == 0) then
  NENOD_OUT = 8
elseif (out_res == 1) then
  NENOD_OUT = 20
elseif (out_res == 2) then
  NENOD_OUT = 8
else
  write(*,'(/,a)')'ERROR: wrong out_res value!'
  stop
endif  

if (out_ncomp > inp_ncomp)then
  write(*,'(/,a)')'ERROR: number of components for output cannot be greater than for input!'
  stop
elseif (out_ncomp>1 .and. out_ncomp /= inp_ncomp)then
  write(*,'(/,a)')'ERROR: not supported components transformation!'
  stop
endif

! Display key information
write(*,'(a)')'-------------------------------'
!write(*,*)'Number of data slices: ',nproc
write(*,*)'number of image frames: ',t_nstep
write(*,*)'input directory:',inp_path
if (inp_ncomp==1) then
  write(*,*) 'input data type: SCALAR'
elseif (inp_ncomp==3) then
  write(*,*) 'input data type: VECTOR'
elseif (inp_ncomp==6) then
  write(*,*) 'input data type: 9C SYMMETRIC TENSOR'
else
  write(*,'(/,a)')'ERROR: unsupported data type!'
  stop
endif
write(*,*)'output directory:',out_path
if (out_ncomp==1) then
  write(*,*) 'output data type: SCALAR'
elseif (out_ncomp==3) then
  write(*,*) 'output data type: VECTOR'
elseif (out_ncomp==6) then
  write(*,*) 'output data type: 9C SYMMETRIC TENSOR'
else
  write(*,'(/,a)')'ERROR: unsupported data type!'
  stop
endif
if (out_format==0) then
  write(*,*)'output format: VTK'
elseif (out_format==1) then
  write(*,*)'output format: Ensight Gold'
else
  write(*,'(/,a)')'ERROR: unsupported output format!'
  stop
endif
if (out_res==0) then
  write(*,*)'resolution: LOW'
elseif (out_res==1) then
  write(*,*)'resolution: MEDIUM'
elseif (out_res==2) then
  write(*,*)'resolution: HIGH'
else
  write(*,'(/,a)')'ERROR: unsupported resolution!'
  stop
endif

write(*,*)'number of output slices: ',out_nslice
  
write(*,'(a)')'-------------------------------'
write(*,'(a)',advance='no')'counting meshes...'
  
! count total number of nodes and elements in all slices
allocate(slice_nnode(out_nslice))
allocate(slice_nelmt(out_nslice))
! Loop over output slices
do i_slice=1,out_nslice    
  
  slice_nnode(i_slice) = 0
  slice_nelmt(i_slice) = 0
  
  call cvd_count_totals_ext_mesh(slice_nproc(i_slice), &
  slice_proc_list(i_slice,1:slice_nproc(i_slice)),proc_width, &
  inp_path,slice_nnode(i_slice),slice_nelmt(i_slice),out_res)    
  !write(*,*)i_slice,slice_nnode(i_slice),slice_nelmt(i_slice)   
enddo
write(*,'(a)')'complete!'
   
if (out_format==0)then
  ! VTK files
  call write_vtu()
elseif (out_format==1)then
  ! Ensight Gold files
  call write_ensight()
else
  write(*,'(/,a)')'ERROR: unsupported ouput format!'
  stop
endif
write(*,'(a)')'-------------------------------'
 
end program visualize
