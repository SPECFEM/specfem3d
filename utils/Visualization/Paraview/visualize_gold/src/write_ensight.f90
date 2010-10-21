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
! Revision: April 06,2010 HNG
subroutine write_ensight ()
  
use visualize_par
implicit none
  
integer :: ios,i_slice
character(len=80) :: file_head
!character(len=80),dimension(out_nslice),optional :: server_name, server_exec
  
! Write a Ensight Gold SOS file
if (out_nslice>1)then  
  open(unit=101, file=trim(out_path)// '/' // trim(out_head)//'.sos', status='replace', action='write', iostat=ios)
    
  write(101,'(a)')'FORMAT'
  write(101,'(a,/)')'type:  master_server gold'

  write(101,'(a)')'SERVERS'
  write(101,'(a,i2,/)')'number of servers:    ',out_nslice
    
  ! Loop over output slices
  do i_slice=1,out_nslice
    write(101,'(a,i2)')'#Server ',i_slice
    write(101,'(a)')'machine id: '//trim(server_name(i_slice))
    write(101,'(a)')'executable: '//trim(server_exec(i_slice))
    write(101,'(a)')'#login id: '
    write(101,'(a)')'#data_path: '
    write(tmp_str,*)i_slice
    file_head=trim(out_head)//'_server'//trim(adjustl(tmp_str))
    write(101,'(a,/)')'casefile: '//trim(file_head)//'.case'
  enddo
  close(101)
endif
  
write(*,'(a)')'writing ensight gold files...'
! Loop over output slices
do i_slice=1,out_nslice
  !write(*,*)'slice: ',i_slice   
  ! counts total number of points
  !nnode = 0
  !nelmt = 0
    
  !call cvd_count_totals_ext_mesh(slice_nproc(i_slice), &
  !slice_proc_list(i_slice,1:slice_nproc(i_slice)),proc_width, &
  !inp_path,nnode,nelmt,out_res)
  !write(*,'(a)')'complete!'
  !write(*,*)'  Total number of nodes: ',nnode
  !write(*,*)'  Total number of elements: ',nelmt

  ! Ensight Gold files
  call write_ensight_serial(i_slice,slice_nproc(i_slice),slice_proc_list(i_slice,1:slice_nproc(i_slice)))    
enddo 
write(*,*)
write(*,'(a)')'complete'
end subroutine write_ensight
!=====================================================================

subroutine write_ensight_serial (i_slice,nproc,proc_list)  

use visualize_par
implicit none
  
integer,intent(in) :: i_slice
integer,intent(in) :: nproc
integer,dimension(nproc),intent(in) :: proc_list
  
character(len=20), parameter :: wild_char='********************'
!integer,dimension(6),parameter :: file_unit = (/ 111, 222, 333, 444, 555, 666 /)
  
integer :: i,i_t,i_proc,iproc
integer :: ios
character(len=80) :: geo_file
character(len=80) :: buffer ! this must be 80 characters long
integer :: npart,ts,tstep
integer :: node_count,elmt_count
integer :: i_comp

! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: dat ! Data points in local nodes
real(kind=CUSTOM_REAL),dimension(:),allocatable :: dat_glob ! data points in global points
real(kind=4),dimension(:,:,:,:),allocatable :: tmp_dat
real(kind=4),dimension(:),allocatable :: tmp_dat_glob
character(len=80) :: file_head
character(len=256) :: mesh_file
character(len=20) :: ensight_etype
character(len=80) :: inp_fname
integer :: nnode,nelmt
  
! Ensight element type 
if (out_res==1)then
  ! Medium resolution
  ! 20-noded hexahedra
  ensight_etype='hexa20'
else
  ! Low or high resolution
  ! 8-noded hexahedra
  ensight_etype='hexa8'
endif

write(tmp_str,*)i_slice
file_head=trim(out_head)//'_server'//trim(adjustl(tmp_str))

! write Ensight Gold case file
ts=1 ! Time set

!write(*,'(a)',advance='no')'writing Ensight case file...'
open(unit=11, file=trim(out_path)// '/' // trim(file_head)//'.case', status='replace', action='write', iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: output file "'//trim(file_head)//'.case'//'" cannot be opened!'
  stop
endif  

write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  ensight gold'

write(11,'(a)')'GEOMETRY'
write(11,'(a,a,/)')'model:    ',trim(file_head)//'.geo'

write(11,'(a)')'VARIABLE'
if (out_ncomp == 1)then
  write(11,'(a,i10,a,a,a,a,/)')'scalar per node: ',ts,' ',trim(out_vname),' ',trim(file_head)//'_'//wild_char(1:t_width)//'.scl'
elseif (out_ncomp == 3)then
  write(11,'(a,i10,a,a,a,a,/)')'vector per node: ',ts,' ',trim(out_vname),' ',trim(file_head)//'_'//wild_char(1:t_width)//'.vec'
elseif (out_ncomp == 6)then
  write(11,'(a,i10,a,a,a,a,/)')'tensor symm per node: ',ts,' ',trim(out_vname),' ', &
  trim(file_head)//wild_char(1:t_width)//'.tns'
else
  write(*,'(/,a,i5,a)')'ERROR: number of components ',out_ncomp,' not supported!'
  stop
endif

write(11,'(a)')'TIME'
write(11,'(a,i10)')'time set:',ts
write(11,'(a,i10)')'number of steps:',t_nstep
write(11,'(a,i10)')'filename start number:',t_start
write(11,'(a,i10)')'filename increment:',t_inc
write(11,'(a)',advance='no')'time values: '
  
do i=1,t_nstep
  write(11,'(e12.5)',advance='yes')(t_start+(i-1)*t_inc)*DT
enddo
close(11)
!write(*,'(a)')'complete!'

!write(*,'(a)',advance='no')'writing Ensight mesh file...'
  
!write(*,*) 'Slice list: '
!write(*,*) proc_list(1:nproc)

! open Ensight Gold geo file to store mesh data
geo_file = trim(out_path) // '/' // trim(file_head)//'.geo' !; write(*,*)geo_file
call open_file2write(trim(geo_file)//char(0),fd)    
   
npart=1        
!write(*,*)nnode,nelmt
buffer='C Binary'
call write_string(buffer//char(0),fd)
buffer='Created by write_ensight Routine'
call write_string(buffer//char(0),fd)
buffer='specfem3d_sesame'   
call write_string(buffer//char(0),fd)
buffer='node id off'     
call write_string(buffer//char(0),fd)
buffer='element id off'
call write_string(buffer//char(0),fd)
!call write_string(extent_str//char(0),fd)
!do j=1,3
!  do i=1,2
!    call write_float(real(extent(i,j)),fd)
!  enddo
!enddo
buffer='part' 
call write_string(buffer//char(0),fd)
call write_integer(npart,fd)
buffer='unstructured meshes'
call write_string(buffer//char(0),fd)
buffer='coordinates' 
call write_string(buffer//char(0),fd)
call write_integer(slice_nnode(i_slice),fd)
  
! writes point and scalar information  
! loops over slices (process partitions)
node_count = 0 

! Open temporary files to store coordinates
call open_file2write('../tmp/tmp_x'//char(0),fd_x)
call open_file2write('../tmp/tmp_y'//char(0),fd_y)
call open_file2write('../tmp/tmp_z'//char(0),fd_z)

write(tmp_str,*)proc_width
write(proc_width_str,*)trim(adjustl(tmp_str))//'.'//trim(adjustl(tmp_str)) 
format_str1='(a,i'//trim(adjustl(proc_width_str))//',a)'
write(tmp_str,*)t_width
write(t_width_str,*)trim(adjustl(tmp_str))//'.'//trim(adjustl(tmp_str))  
format_str2='(a,i'//trim(adjustl(proc_width_str))//',a,i'//trim(adjustl(t_width_str))//',a)'
do i_proc = 1, nproc

  iproc = proc_list(i_proc)

  !write(*,*) ' '
  !write(*,*) 'Reading slice ', iproc

  ! gets number of elements and global points for this partition
  write(mesh_file,fmt=format_str1) trim(inp_path)//'/proc',iproc,'_external_mesh.bin'    
  open(unit=27,file=trim(mesh_file),status='old',action='read',form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(*,'(/,a)')'ERROR: file '//trim(mesh_file)//' cannot be opened!'
    stop
  endif    
  read(27) NSPEC_AB
  read(27) NGLOB_AB    
  ! ibool file
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  read(27) ibool    

  ! global point arrays
  allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB)) 
  read(27) xstore
  read(27) ystore
  read(27) zstore
  close(27)    
    
  ! writes point coordinates and scalar value to mesh file
  if (out_res==0)then
    ! writes out element corners only
    call cvd_write_corners_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
    nnode,fd_x,fd_y,fd_z)
  elseif (out_res==1)then
    ! writes out element corners only
    call cvd_write_hexa20_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
    nnode,fd_x,fd_y,fd_z)
  elseif (out_res==2)then  
    ! high resolution, all GLL points
    call cvd_write_GLL_points_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,&
    nnode,fd_x,fd_y,fd_z)
  else
    write(*,'(/,a)')'ERROR: wrong out_res value!'
    stop
  endif
    
  !write(*,*)'  points:',node_count,nnode
    
  ! stores total number of points written
  node_count = node_count + nnode

  ! cleans up memory allocations
  deallocate(ibool,xstore,ystore,zstore)
    
enddo  ! all slices for points
call close_file(fd_x)
call close_file(fd_y)
call close_file(fd_z)

if (node_count /=  slice_nnode(i_slice))then
  write(*,'(/,a)')'Error: number of total points are not consistent!'
  stop
endif
!write(*,*) 'Total number of points: ', node_count
!write(*,*) ' '

!write coordinates to file
!call open_file2read('tmp_x'//char(0),fd_x)
call open_file2read('../tmp/tmp_x'//char(0),fd_x)
do i=1,slice_nnode(i_slice)
  call read_float(tmp_real,fd_x)
  !write(*,*)'new:',tmp_real
  !stop
  call write_float(tmp_real,fd)   
enddo

call close_delete_file('../tmp/tmp_x'//char(0),fd_x)


call open_file2read('../tmp/tmp_y'//char(0),fd_y)
do i=1,slice_nnode(i_slice)
  call read_float(tmp_real,fd_y)
  call write_float(tmp_real,fd)   
enddo
call close_delete_file('../tmp/tmp_y'//char(0),fd_y)

call open_file2read('../tmp/tmp_z'//char(0),fd_z)
do i=1,slice_nnode(i_slice)
  call read_float(tmp_real,fd_z)
  call write_float(tmp_real,fd)   
enddo
call close_delete_file('../tmp/tmp_z'//char(0),fd_z)


! writes element information
buffer=ensight_etype
call write_string(buffer//char(0),fd)
call write_integer(slice_nelmt(i_slice),fd)
elmt_count = 0
node_count = 0
do i_proc = 1, nproc

  iproc = proc_list(i_proc)

  ! gets number of elements and global points for this partition
  write(mesh_file,fmt=format_str1) trim(inp_path)//'/proc',iproc,'_external_mesh.bin'    
  open(unit=27,file=trim(mesh_file),status='old',action='read',form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(*,'(/,a)')'ERROR: file '//trim(mesh_file)//' cannot be opened!'
    stop
  endif    
  read(27) NSPEC_AB
  read(27) NGLOB_AB    
  ! ibool file
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  read(27) ibool
  close(27)

  ! writes out element corner indices
  if (out_res==0) then
    ! spectral elements
    call cvd_write_corner_elements(NSPEC_AB,NGLOB_AB,ibool, &
    node_count,nelmt,nnode,fd)  
  elseif (out_res==1) then
    ! spectral elements
    call cvd_write_hexa20_elements(NSPEC_AB,NGLOB_AB,ibool, &
    node_count,nelmt,nnode,fd)
  elseif (out_res==2)then 
    ! subdivided spectral elements
    call cvd_write_GLL_elements(NSPEC_AB,NGLOB_AB,ibool, &
    node_count,nelmt,nnode,fd)  
  else
    write(*,'(/,a)')'ERROR: wrong out_res value!'
    stop
  endif
  
  !write(*,*)'  elements:',elmt_count,nelmt
  !write(*,*)'  points : ',node_count,nnode
  
  elmt_count = elmt_count + nelmt

  deallocate(ibool)

enddo ! nproc

! close mesh file
call close_file(fd)

! checks with total number of elements
if (elmt_count /= slice_nelmt(i_slice)) then 
  !write(*,'(/,a)')'ERROR: number of elements counted:',elmt_count,'total:',slice_nelmt(i_slice)
  write(*,'(/,a)')'ERROR: number of total elements are not consistent!'
  stop
endif
!write(*,*) 'Total number of elements: ', elmt_count

!write(*,'(a)')'complete!'  

!write(*,*) 'Done writing '//trim(geo_file)

! Write data files

! Format for progress display
write(tmp_str,*)ceiling(log10(real(t_nstep)+1))
format_str3='(a,a,i'//trim(adjustl(tmp_str))//',a,i'//trim(adjustl(tmp_str))
write(tmp_str,*)ceiling(log10(real(out_nslice)+1))
format_str3=trim(format_str3)//',a,i'//trim(adjustl(tmp_str))//',a,i'//trim(adjustl(tmp_str))//')'

if (out_ncomp==1) then
  out_ext='.scl'
elseif (out_ncomp==3) then
  out_ext='.vec'
elseif (out_ncomp==6) then
  out_ext='.tns'
else
  write(*,'(/,a,i5,a)')'ERROR: number of components ',out_ncomp,' not supported!'
  stop
endif 

!write(*,'(a)',advance='no')'time step: '
do i_t=1,t_nstep
  tstep=t_start + (i_t-1)*t_inc

  ! Open Ensight Gold data file to store data
  write(out_fname,fmt=format_str1)trim(out_path) // '/'//trim(file_head)//'_',tstep,trim(out_ext)  
  npart=1;
  call open_file2write(trim(out_fname)//char(0),fd)
  buffer='Scalar data'
  call write_string(buffer//char(0),fd)
  buffer='part'
  call write_string(buffer//char(0),fd)
  call write_integer(npart,fd)
  buffer='coordinates'

  call write_string(buffer//char(0),fd)

  !Open temporary files to store data
  !call open_file2write('tmp_val'//char(0),fd_x)
  

  if (out_ncomp>1) then
    do i_comp=1,out_ncomp
      node_count = 0
      do i_proc = 1, nproc

        iproc = proc_list(i_proc)

        write(mesh_file,fmt=format_str1) trim(inp_path)//'/proc',iproc,'_external_mesh.bin'           
        open(unit=27,file=trim(mesh_file),status='old',action='read',form='unformatted')
        read(27) NSPEC_AB
        read(27) NGLOB_AB
  
        ! ibool file
        allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))    
        read(27) ibool
        close(27)
        if (dat_topo == 0)then
          allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB)) 

          ! data file  
          write(inp_fname,fmt=format_str2)trim(out_path)//'/'//trim(proc_head), &
          iproc,trim(inp_head(i_comp)),tstep,trim(inp_ext)
          open(unit = 11,file = trim(inp_fname),status='old',&
            action='read', iostat = ios,form ='unformatted')
          if (ios /= 0) then
            write(*,*)'Error opening '//trim(inp_fname)
            stop
          endif
            
          read(11) dat     
      
          ! writes point coordinates and scalar value to mesh file
          if (out_res==0) then 
            call cvd_write_corners_data(NSPEC_AB,NGLOB_AB,ibool,real(dat), &
            nnode,fd)  
          elseif (out_res==1) then 
            call cvd_write_hexa20_data(NSPEC_AB,NGLOB_AB,ibool,real(dat), &
            nnode,fd)
          elseif (out_res==2) then  
            call cvd_write_GLL_points_data(NSPEC_AB,NGLOB_AB,ibool,real(dat), &
            nnode,fd)
          else
            write(*,'(/,a)')'ERROR: wrong out_res value!'
            stop
          endif
          ! cleans up memory allocations
          deallocate(dat)
        elseif (dat_topo==1)then
          allocate(dat_glob(NGLOB_AB)) 

          ! data file  
          write(inp_fname,fmt=format_str2)trim(out_path)//'/'//trim(proc_head), &
          iproc,trim(inp_head(i_comp)),tstep,trim(inp_ext)
          open(unit = 11,file = trim(inp_fname),status='old',&
            action='read', iostat = ios,form ='unformatted')
          if (ios /= 0) then
            write(*,*)'Error opening ',trim(inp_fname)
            stop
          endif
            
          read(11) dat_glob     
      
          ! writes point coordinates and scalar value to mesh file
          if (out_res==0) then 
            call cvd_write_corners_data_glob(NSPEC_AB,NGLOB_AB,ibool,real(dat_glob), &
            nnode,fd)  
          elseif (out_res==1) then 
            call cvd_write_hexa20_data_glob(NSPEC_AB,NGLOB_AB,ibool,real(dat_glob), &
            nnode,fd)
          elseif (out_res==2) then  
            call cvd_write_GLL_points_data_glob(NSPEC_AB,NGLOB_AB,ibool,real(dat_glob), &
            nnode,fd)
          else
            write(*,'(/,a)')'ERROR: wrong out_res value!'
            stop
          endif
          ! cleans up memory allocations
          deallocate(dat_glob)
        endif ! if dat_topo == 0
        deallocate(ibool)
  
        !write(*,*)'  points:',node_count,nnode
  
        ! stores total number of points written
        node_count = node_count + nnode          
  
      enddo  ! i_proc = 1, nproc
      if (node_count /=  slice_nnode(i_slice))then
        write(*,'(/,a)')'Error: Number of total points are not consistent'
        stop
      endif
    enddo
  else ! if (out_ncomp>1)
    node_count = 0
    do i_proc = 1, nproc

      iproc = proc_list(i_proc)

      write(mesh_file,fmt=format_str1) trim(inp_path)//'/proc',iproc,'_external_mesh.bin'           
      open(unit=27,file=trim(mesh_file),status='old',action='read',form='unformatted')
      read(27) NSPEC_AB
      read(27) NGLOB_AB
  
      ! ibool file
      allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))    
      read(27) ibool
      close(27)
      
      if (dat_topo==0)then
        allocate(tmp_dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
        allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB)) 
        tmp_dat=0.0
        do i_comp=1,inp_ncomp 
          ! data file  
          write(inp_fname,fmt=format_str2)trim(out_path)//'/'//trim(proc_head), &
            iproc,trim(inp_head(i_comp)),tstep,trim(inp_ext)
          open(unit = 222,file = trim(inp_fname),status='old',&
            action='read', iostat = ios,form ='unformatted')
          if (ios /= 0) then
            write(*,*)'Error opening '//trim(inp_fname)
            stop
          endif
            
          read(222) dat
          close(222)
          tmp_dat=tmp_dat+real(dat)
          !write(*,*)inp_fname
        enddo
        if (inp_ncomp==3 .and. out_ncomp==1)then
          tmp_dat=0.5*tmp_dat ! Equivalent to S-wave potential
        endif       
      
        ! writes point coordinates and scalar value to mesh file
        if (out_res==0) then 
          call cvd_write_corners_data(NSPEC_AB,NGLOB_AB,ibool,tmp_dat, &
          nnode,fd)  
        elseif (out_res==1) then 
          call cvd_write_hexa20_data(NSPEC_AB,NGLOB_AB,ibool,tmp_dat, &
          nnode,fd)
        elseif (out_res==2) then  
          call cvd_write_GLL_points_data(NSPEC_AB,NGLOB_AB,ibool,tmp_dat, &
          nnode,fd)
        else
          write(*,'(/,a)')'ERROR: wrong out_res value!'
          stop
        endif
        ! cleans up memory allocations
        deallocate(ibool,dat,tmp_dat)
      elseif (dat_topo==1)then
        allocate(tmp_dat_glob(NGLOB_AB))
        allocate(dat_glob(NGLOB_AB)) 
        tmp_dat=0.0
        do i_comp=1,inp_ncomp 
          ! data file  
          write(inp_fname,fmt=format_str2)trim(out_path)//'/'//trim(proc_head), &
            iproc,trim(inp_head(i_comp)),tstep,trim(inp_ext)
          open(unit = 11,file = trim(inp_fname),status='old',&
            action='read', iostat = ios,form ='unformatted')
          if (ios /= 0) then
            write(*,'(/,a)')'ERROR: opening '//trim(inp_fname)
            stop
          endif
            
          read(11) dat_glob
          tmp_dat_glob=tmp_dat_glob+real(dat_glob)
          !write(*,*)inp_fname
        enddo
        if (inp_ncomp==3 .and. out_ncomp==1)then
          tmp_dat_glob=0.5*tmp_dat_glob ! Equivalent to S-wave potential
        endif       
      
        ! writes point coordinates and scalar value to mesh file
        if (out_res==0) then 
          call cvd_write_corners_data_glob(NSPEC_AB,NGLOB_AB,ibool,tmp_dat_glob, &
          nnode,fd)  
        elseif (out_res==1) then 
          call cvd_write_hexa20_data_glob(NSPEC_AB,NGLOB_AB,ibool,tmp_dat_glob, &
          nnode,fd)
        elseif (out_res==2) then
          call cvd_write_GLL_points_data_glob(NSPEC_AB,NGLOB_AB,ibool,tmp_dat_glob, &
          nnode,fd)
        else
          write(*,'(/,a)')'ERROR: wrong out_res value!'
          stop
        endif
        ! cleans up memory allocations
        deallocate(ibool,dat_glob,tmp_dat_glob)
      endif ! if dat_topo == 0
  
      !write(*,*)'  points:',node_count,nnode
  
      ! stores total number of points written
      node_count = node_count + nnode
      
  
    enddo  ! i_proc = 1, nproc

    if (node_count /=  slice_nnode(i_slice)) then
      write(*,'(/,a)')'Error: Number of total points are not consistent'
      stop
    endif
    !write(*,*) 'Total number of points: ', node_count
    !write(*,*) ' '

    call close_file(fd)
  endif ! if (out_ncomp>1)

  ! Display progress
  write(*,fmt=format_str3,advance='no')CR,' slice: ',i_slice,'/',out_nslice,', time step: ',i_t,'/',t_nstep
  
enddo ! do i_t=1,t_nstep

!write(*,'(a)')' complete!'

end subroutine write_ensight_serial
!=============================================================
 
