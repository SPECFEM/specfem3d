subroutine write_vtu

! This programs writes vtu (VTK binary unstructed data file)
! files produced by SPECFEM3D. The vtu file can then be visualized in 
! VTK/ParaView. See http://www.vtk.org and http://www.paraview.org for details.
!------------------------------------------
! DEPENDENCY:
!   cfunc4fortran.c, visualize_par.f90, visualize_collect.f90  
! COMPILE
!   >> gfortran -c write_vtu
! USAGE
!   write_vtu
! HISTORY:
!   Hom Nath Gharti, NORSAR
!   April 06,2010; Mar 10,2010 (NORSAR)
! FEEDBACK:
!   homnath_AT_norsar_DOT_no
!------------------------------------------
use visualize_par
implicit none

!integer,dimension(nproc),intent(in) :: proc_list

integer :: vtk_etype
!integer,dimension(6),parameter :: file_unit = (/ 111, 222, 333, 444, 555, 666 /)
integer,parameter :: plot_nvar = 5
integer :: endian
integer,dimension(10) :: bytes,off
integer,parameter :: LE=0,BE=1
integer :: i,j,i_t,i_proc,iproc,i_slice 
integer :: ios
integer :: tstep
integer :: node_count,elmt_count
integer :: i_comp

! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: dat
real(kind=CUSTOM_REAL),dimension(:),allocatable :: dat_glob 
real(kind=4),dimension(:,:,:,:),allocatable :: tmp_dat
real(kind=4),dimension(:),allocatable :: tmp_rvect,tmp_dat_glob

character(len=12) :: byte_order
!character(len=60) :: tmp_str,
character(len=256) :: buffer,pvd_file,pvtu_file,vtu_file,mesh_file
integer,parameter :: pvd_unit=11,pvtu_unit=22,vtu_unit=33
character(len=80) :: file_head,inp_fname
integer :: nnode,nelmt,tmp_nnode !,tmp_nelmt  

write(*,'(a)')'writing VTK files...'
! Determine the Endianness of the Architecture
call get_endian(endian)
if(endian == LE)then
  byte_order='LittleEndian'
elseif(endian == BE)then
  byte_order='BigEndian'
else
  write(*,'(/,a)')'ERROR: illegal endianness!'
  stop
endif

! vtk element type
if (out_res==1)then
  ! Medium resolution
  ! 20-noded hexahedra
  vtk_etype=25
else
  ! Low and high resolution
  ! 12-noded hexahedra
  vtk_etype=12
endif

write(tmp_str,*)proc_width
write(proc_width_str,*)trim(adjustl(tmp_str))//'.'//trim(adjustl(tmp_str))  
write(tmp_str,*)t_width
write(t_width_str,*)trim(adjustl(tmp_str))//'.'//trim(adjustl(tmp_str))
format_str1='(a,i'//trim(adjustl(proc_width_str))//',a)'  
format_str2='(a,i'//trim(adjustl(proc_width_str))//',a,i'//trim(adjustl(t_width_str))//',a)'

! Format for progress display
write(tmp_str,*)ceiling(log10(real(t_nstep)+1))
format_str3='(a,a,i'//trim(adjustl(tmp_str))//',a,i'//trim(adjustl(tmp_str))
write(tmp_str,*)ceiling(log10(real(out_nslice)+1))
format_str3=trim(format_str3)//',a,i'//trim(adjustl(tmp_str))//',a,i'//trim(adjustl(tmp_str))//')'
out_ext='.vtu'

!write(*,*)file_head
!stop
! Open pvd file
pvd_file=trim(out_path)//'/'// trim(out_head)//'.pvd'
open(unit=pvd_unit, file=trim(pvd_file), status='replace', action='write', iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: output file "'//trim(pvd_file)//'" cannot be opened!'
  stop
endif
buffer='<?xml version="1.0"?>'
write(pvd_unit,'(a)')trim(buffer)
buffer='<VTKFile type="Collection" version="0.1" byte_order="'//trim(byte_order)//'">'
write(pvd_unit,'(a)')trim(buffer)
buffer='<Collection>'
write(pvd_unit,'(a)')trim(buffer)

!slice_nnode=0
!slice_nelmt=0
do i_t=1,t_nstep

  tstep=t_start + (i_t-1)*t_inc
  
  ! create pvtu file name   
  write(pvtu_file,fmt=format_str1)trim(out_head)//'_',tstep,'.pvtu'
  
  ! collect pvtu file name in pvd file   
  write(num_str1,'(f16.6)')(t_start+(i_t-1)*t_inc)*dt !  Change format here if time is so big and so tiny
  buffer='<DataSet timestep="'//trim(adjustl(num_str1))//'" part="001" file="'//trim(pvtu_file)//'"/>'
  write(pvd_unit,'(a)')trim(buffer)
  
  ! open pvtu file
  pvtu_file = trim(out_path) // '/' // trim(pvtu_file)    
  open(unit=pvtu_unit, file=trim(pvtu_file), action='write', status='replace')
  ! write headers
  
  buffer='<?xml version="1.0"?>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="'//trim(byte_order)//'">';
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='<PUnstructuredGrid GhostLevel="0">'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='<PPoints>'
  write(pvtu_unit,'(a)')trim(buffer)
  write(num_str1,*)off(1);
  buffer='<PDataArray type="Float32" NumberOfComponents="3" format="appended"/>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='</PPoints>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='<PCells>'
  write(pvtu_unit,'(a)')trim(buffer)
  write(num_str1,*)off(3);
  buffer='<PDataArray type="Int32" Name="connectivity" format="appended"/>'
  write(pvtu_unit,'(a)')trim(buffer)
  write(num_str1,*)off(4);
  buffer='<PDataArray type="Int32" Name="offsets" format="appended"/>'
  write(pvtu_unit,'(a)')trim(buffer)
  write(num_str1,*)off(5);
  buffer='<PDataArray type="Int32" Name="types" format="appended"/>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='</PCells>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='<PPointData>'
  write(pvtu_unit,'(a)')trim(buffer)
  write(num_str1,*)off(2);
  buffer='<PDataArray type="Float32" NumberOfComponents="1" Name="'//trim(out_vname)// &
          '" format="appended"/>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='</PPointData>'
  write(pvtu_unit,'(a)')trim(buffer)
  
  buffer='<PCellData>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='</PCellData>'
  write(pvtu_unit,'(a)')trim(buffer)    
  
  
  do i_slice=1,out_nslice ! out processors
    
    ! counts total number of points in each slice
    !slice_nnode(i_slice) = 0
    !slice_nelmt(i_slice) = 0
    
    !call cvd_count_totals_ext_mesh(slice_nproc(i_slice), &
    !slice_proc_list(i_slice,1:slice_nproc(i_slice)),proc_width, &
    !inp_path,nnode,nelmt,out_res)
    !write(*,'(a)')'complete!'
    !write(*,*)'  Total number of nodes: ',nnode
    !write(*,*)'  Total number of elements: ',nelmt
    
    !slice_nnode(i_slice)=nnode
    !slice_nelmt(i_slice)=nelmt
    
    ! Compute bytes and offsets
    bytes(1) = (ndim*slice_nnode(i_slice))*size_float ! Coordinates
    bytes(2) = (NENOD_OUT*slice_nelmt(i_slice))*size_int ! Connectivity
    bytes(3) = (slice_nelmt(i_slice))*size_int ! Offsets
    bytes(4) = (slice_nelmt(i_slice))*size_int ! Types
    bytes(5) = (out_ncomp*slice_nnode(i_slice))*size_float ! Nodal values        
  
    off(1)=0; ! 1st offset
    do i=1,plot_nvar
      if(i < plot_nvar)then
        off(i+1)=off(i)+size_int+bytes(i)
      endif
      bytes(i)=bytes(i)+size_int
    enddo
    
    ! create vtu file
    write(tmp_str,*)i_slice
    file_head=trim(out_head)//'_server'//trim(adjustl(tmp_str))      
    write(out_fname,fmt=format_str1)trim(file_head)//'_',tstep,trim(out_ext)
    !write(*,*)tstep,trim(out_fname)
    ! write vtu file to pvtu file      
    buffer='<Piece Source="'//trim(out_fname)//'"/>'
    write(pvtu_unit,'(a)')trim(buffer)        

    ! open vtu file
    vtu_file = trim(out_path) // '/' // trim(out_fname)    
    open(unit=vtu_unit, file=trim(vtu_file), action='write', status='replace',iostat=ios)
    !write(*,*)trim(vtu_file),vtu_unit
    if (ios/=0)then
      write(*,'(/,a)')'ERROR: file '//trim(vtu_file)//' cannot be opened!'
      stop
    endif
        
    ! write header
    buffer='<?xml version="1.0"?>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="'//trim(byte_order)//'">';
    write(vtu_unit,'(a)')trim(buffer)
    buffer='<UnstructuredGrid>'
    write(vtu_unit,'(a)')trim(buffer)
    write(num_str1,*)slice_nnode(i_slice); write(num_str2,*)slice_nelmt(i_slice); 
    buffer='<Piece NumberOfPoints="'//trim(adjustl(num_str1))//'" NumberOfCells="'//trim(adjustl(num_str2))//'">'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='<Points>'
    write(vtu_unit,'(a)')trim(buffer)
    write(num_str1,*)off(1);
    buffer='<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="'//trim(adjustl(num_str1))//'"/>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='</Points>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='<Cells>'
    write(vtu_unit,'(a)')trim(buffer)
    write(num_str1,*)off(2);
    buffer='<DataArray type="Int32" Name="connectivity" format="appended" offset="'//trim(adjustl(num_str1))//'"/>'
    write(vtu_unit,'(a)')trim(buffer)
    write(num_str1,*)off(3);
    buffer='<DataArray type="Int32" Name="offsets" format="appended" offset="'//trim(adjustl(num_str1))//'"/>'
    write(vtu_unit,'(a)')trim(buffer)
    write(num_str1,*)off(4);
    buffer='<DataArray type="Int32" Name="types" format="appended" offset="'//trim(adjustl(num_str1))//'"/>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='</Cells>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='<PointData>'
    write(vtu_unit,'(a)')trim(buffer)
    write(num_str1,*)off(5);
    buffer='<DataArray type="Float32" NumberOfComponents="1" Name="'//trim(out_vname)// &
            '" format="appended" offset="'//trim(adjustl(num_str1))//'"/>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='</PointData>'
    write(vtu_unit,'(a)')trim(buffer)      
    buffer='<CellData>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='</CellData>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='</Piece>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='</UnstructuredGrid>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='<AppendedData encoding="raw">'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='_'
    write(vtu_unit,'(a)',advance='no')trim(buffer)
    close(vtu_unit)
      
    if (i_t==1)then
      ! open temporary files to store coordinates
      write(tmp_str,*)i_slice
      call open_file2write('../tmp/tmp_x_slice'//trim(adjustl(tmp_str))//char(0),fd_x)
      call open_file2write('../tmp/tmp_y_slice'//trim(adjustl(tmp_str))//char(0),fd_y)
      call open_file2write('../tmp/tmp_z_slice'//trim(adjustl(tmp_str))//char(0),fd_z)
      ! open temporary files to store connectivity list
      call open_file2write('../tmp/tmp_con_slice'//trim(adjustl(tmp_str))//char(0),fd_con)
      
      ! open mesh files and store mesh data
      tmp_nnode = 0
      elmt_count = 0
      node_count = 0
      do i_proc = 1, slice_nproc(i_slice)

        iproc = slice_proc_list(i_slice,i_proc)          

        ! Read and store mesh data
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
        if (out_res==0) then
          ! writes out element corners only
          call cvd_write_corners_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
          nnode,fd_x,fd_y,fd_z)
        elseif (out_res==1)then
          ! writes out element corners only
          call cvd_write_hexa20_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
          nnode,fd_x,fd_y,fd_z)
        elseif(out_res==2)then
          ! high resolution, all GLL points
          call cvd_write_GLL_points_only(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,&
          nnode,fd_x,fd_y,fd_z)
        else
          write(*,'(/,a)')'ERROR: wrong out_res value!'
          stop
        endif
        
        !write(*,*)'  points:',node_count,nnode
        
        ! stores total number of points written
        tmp_nnode = tmp_nnode + nnode

        ! cleans up memory allocations
        deallocate(xstore,ystore,zstore)
        
        ! Read and store connectivity list 

        ! writes out element corner indices
        if(out_res==0) then
          ! spectral elements
          call cvd_write_corner_elements(NSPEC_AB,NGLOB_AB,ibool, &
          node_count,nelmt,nnode,fd_con)
        elseif (out_res==1) then
          ! spectral elements
          call cvd_write_hexa20_elements(NSPEC_AB,NGLOB_AB,ibool, &
          node_count,nelmt,nnode,fd_con)
        elseif(out_res==2)then
          ! subdivided spectral elements
          call cvd_write_GLL_elements(NSPEC_AB,NGLOB_AB,ibool, &
          node_count,nelmt,nnode,fd_con)  
        else
          write(*,'(/,a)')'ERROR: wrong out_res value!'
          stop
        endif
        deallocate(ibool)
        
        !write(*,*)'  elements:',elmt_count,nelmt
        !write(*,*)'  points : ',node_count,nnode
        !write(*,*)tmp_nnode,node_count
        if (tmp_nnode/=node_count)then
          write(*,'(/,a)')'ERROR: inconsistent number of nodes!'
          stop
        endif
        
        !write(*,*)node_count
        !node_count = node_count + nnode          
        elmt_count = elmt_count + nelmt          
        
      enddo ! i_proc
      ! close temporary files
      call close_file(fd_x)
      call close_file(fd_y)
      call close_file(fd_z)
      call close_file(fd_con)
        
      !write(*,*)node_count,slice_nnode(i_slice)
      if (node_count /=  slice_nnode(i_slice)) stop 'Error: Number of total points are not consistent'
      ! checks with total number of elements
      if (elmt_count /= slice_nelmt(i_slice)) then 
        !write(*,'(/,a)')'ERROR: number of elements counted:',elmt_count,'total:',slice_nelmt(i_slice)
        write(*,*)'Number of total elements are not consistent!'
        stop
      endif
    endif
    
    !write(*,*)trim(vtu_file)
    ! write coordinates to vtu file          
    call open_file2append(trim(vtu_file)//char(0),fd) 
    !write(*,*)'hi'
    call write_integer(bytes(1),fd) 
    write(tmp_str,*)i_slice
    call open_file2read('../tmp/tmp_x_slice'//trim(adjustl(tmp_str))//char(0),fd_x)
    call open_file2read('../tmp/tmp_y_slice'//trim(adjustl(tmp_str))//char(0),fd_y)
    call open_file2read('../tmp/tmp_z_slice'//trim(adjustl(tmp_str))//char(0),fd_z)
    do i=1,slice_nnode(i_slice)
      call read_float(tmp_real,fd_x); call write_float(tmp_real,fd)
      call read_float(tmp_real,fd_y); call write_float(tmp_real,fd)
      call read_float(tmp_real,fd_z); call write_float(tmp_real,fd)         
    enddo 
    call close_file(fd_x)
    call close_file(fd_y)
    call close_file(fd_z)  
      
    ! write connectivity to vtu file
    call write_integer(bytes(2),fd) 
    call open_file2read('../tmp/tmp_con_slice'//trim(adjustl(tmp_str))//char(0),fd_con) 
    do i=1,slice_nelmt(i_slice)
      do j=1,NENOD_OUT
        call read_integer(tmp_int,fd_con)
        call write_integer(tmp_int-1,fd)
        !write(*,*)tmp_int-1
      enddo
    enddo
    call close_file(fd_con)    

    ! write offsets
    call write_integer(bytes(3),fd)        
    tmp_int=0
    do i=1,slice_nelmt(i_slice)
      tmp_int=tmp_int+NENOD_OUT
      call write_integer(tmp_int,fd)
    enddo

    ! Write element types
    call write_integer(bytes(4),fd) 
    do i=1,slice_nelmt(i_slice) 
      call write_integer(vtk_etype,fd)
    enddo        
    
    ! write data to vtu file
    call write_integer(bytes(5),fd)
    
    if (out_ncomp>1)then ! vector or tensor data
      if (i_t==1)then
        allocate(fd_array(out_ncomp))
        allocate(tmp_rvect(out_ncomp))
      endif
      ! open temporary files
      do i_comp=1,out_ncomp
        write(tmp_str,*)i_comp
        call open_file2write('../tmp/tmp_data_comp'//trim(adjustl(tmp_str))//char(0),fd_array(i_comp))
      enddo
      node_count=0
      do i_proc = 1, slice_nproc(i_slice)

        iproc=slice_proc_list(i_slice,i_proc)               

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
        
        if (dat_topo==0)then ! Data from local points            
          allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB)) 
          
          do i_comp=1,inp_ncomp 
            ! data file  
            write(inp_fname,fmt=format_str2)trim(out_path)//'/'//trim(proc_head), &
              iproc,trim(inp_head(i_comp)),tstep,trim(inp_ext)
            open(unit = 222,file = trim(inp_fname),status='old',&
              action='read', iostat = ios,form ='unformatted')
            if (ios /= 0) then
              write(*,'(/,a)')'ERROR: opening '//trim(inp_fname)
              stop
            endif
              
            read(222) dat
            close(222)
            
            ! writes point data to file
            if (out_res==0) then 
              call cvd_write_corners_data(NSPEC_AB,NGLOB_AB,ibool,real(dat), &
              nnode,fd_array(i_comp))  
            elseif (out_res==1) then 
              call cvd_write_hexa20_data(NSPEC_AB,NGLOB_AB,ibool,real(dat), &
              nnode,fd_array(i_comp))
            elseif (out_res==2) then  
              call cvd_write_GLL_points_data(NSPEC_AB,NGLOB_AB,ibool,real(dat), &
              nnode,fd_array(i_comp))
            else
              write(*,'(/,a)')'ERROR: wrong out_res value!'
              stop
            endif          
            
          enddo ! i_comp
          ! cleans up memory allocations
          deallocate(dat)          
          
        elseif (dat_topo==1)then ! Data from global points            
          allocate(dat_glob(NGLOB_AB)) 
          
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
            
            ! writes point coordinates and scalar value to mesh file
            if (out_res==0) then 
              call cvd_write_corners_data_glob(NSPEC_AB,NGLOB_AB,ibool,real(dat_glob), &
              nnode,fd_array(i_comp))  
            elseif (out_res==1) then 
              call cvd_write_hexa20_data_glob(NSPEC_AB,NGLOB_AB,ibool,real(dat_glob), &
              nnode,fd_array(i_comp))
            elseif (out_res==2) then
              call cvd_write_GLL_points_data_glob(NSPEC_AB,NGLOB_AB,ibool,real(dat_glob), &
              nnode,fd_array(i_comp))
            else
              write(*,'(/,a)')'ERROR: wrong out_res value!'
              stop
            endif
          enddo ! i_comp
        
          
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
      
      ! close temporary files
      do i_comp=1,out_ncomp
        call close_file(fd_array(i_comp))
      enddo
      
      ! open temporary files to read
      do i_comp=1,out_ncomp
        write(tmp_str,*)i_comp
        call open_file2read('../tmp/tmp_data_comp'//trim(adjustl(tmp_str))//char(0),fd_array(i_comp))
      enddo
      if (out_ncomp==3)then
      ! vector
        do i=1,slice_nnode(i_slice)
          do i_comp=1,out_ncomp
            call read_float(tmp_real,fd_array(i_slice))
            call write_float(tmp_real,fd)
          enddo
        enddo
      elseif (out_ncomp==6)then
      ! 9-component symmetric tensor
        do i=1,slice_nnode(i_slice)
          do i_comp=1,out_ncomp
            call read_float(tmp_rvect(i_comp),fd_array(i_slice))              
          enddo
          call write_float(tmp_rvect(1),fd); call write_float(tmp_rvect(4),fd); call write_float(tmp_rvect(5),fd)
          call write_float(tmp_rvect(4),fd); call write_float(tmp_rvect(2),fd); call write_float(tmp_rvect(6),fd)
          call write_float(tmp_rvect(5),fd); call write_float(tmp_rvect(6),fd); call write_float(tmp_rvect(3),fd)            
        enddo
      else
        write(*,'(/,a)')'ERROR: wrong out_ncomp value!'
        stop
      endif
      ! close temporary data files
      do i_comp=1,out_ncomp
        call close_file(fd_array(i_comp))
      enddo
    else ! out_ncomp        
    
      node_count=0
      do i_proc = 1, slice_nproc(i_slice)

        iproc=slice_proc_list(i_slice,i_proc)               

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
        
        if (dat_topo==0)then ! Data from local points
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
              write(*,'(/,a)')'ERROR: opening '//trim(inp_fname)
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
        
          ! writes point data to file
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
        elseif (dat_topo==1)then ! Data from global points
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
      call close_file(fd)
      if (node_count /=  slice_nnode(i_slice))then
        write(*,'(/,a)')'Error: Number of total points are not consistent'
        stop
        endif
    endif ! out_ncomp

    ! Write post header for vtu file
    open(unit=vtu_unit, file=trim(vtu_file), action='write', status='old',position='append')    
    write(vtu_unit,*) ! Write new line    
    buffer='</AppendedData>'
    write(vtu_unit,'(a)')trim(buffer)
    buffer='</VTKFile>'
    write(vtu_unit,'(a)')trim(buffer)
    close(vtu_unit);
    
    ! Display progress
    write(*,fmt=format_str3,advance='no')CR,' slice: ',i_slice,'/',out_nslice,', time step: ',i_t,'/',t_nstep
      
  enddo  ! do i_slice
  
  ! Write post header for pvtu file        
  buffer='</PUnstructuredGrid>'
  write(pvtu_unit,'(a)')trim(buffer)
  buffer='</VTKFile>'
  write(pvtu_unit,'(a)')trim(buffer)
  close(pvtu_unit);
  
  ! Display progress
  !write(*,fmt=format_str3,advance='no')CR,'writing VTK files... ',i_t,'/',t_nstep
  
enddo ! i_t  

! write post header for pvd file
buffer='</Collection>'
write(pvd_unit,'(a)')trim(buffer)
buffer='</VTKFile>'
write(pvd_unit,'(a)')trim(buffer)
close(pvd_unit);

! delete temporary mesh files
do i_slice=1,out_nslice    
  write(tmp_str,*)i_slice
  call delete_file('../tmp/tmp_x_slice'//trim(adjustl(tmp_str))//char(0))
  call delete_file('../tmp/tmp_y_slice'//trim(adjustl(tmp_str))//char(0))
  call delete_file('../tmp/tmp_z_slice'//trim(adjustl(tmp_str))//char(0))
  call delete_file('../tmp/tmp_con_slice'//trim(adjustl(tmp_str))//char(0))
enddo

if (out_ncomp>1)then
  ! free memory
  deallocate(fd_array,tmp_rvect)
  ! delete temporary data files
  do i_comp=1,out_ncomp
    write(tmp_str,*)i_comp
    call delete_file('../tmp/tmp_data_comp'//trim(adjustl(tmp_str))//char(0))
  enddo
endif
write(*,*)
write(*,'(a)')' complete!'

end subroutine write_vtu
!=============================================================
