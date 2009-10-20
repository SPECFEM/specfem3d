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

  program combine_paraview_data_ext_mesh

! puts the output of SPECFEM3D in ParaView format.
! see http://www.paraview.org for details

! combines the database files on several slices.
! the local database file needs to have been collected onto the frontend (copy_local_database.pl)

! works for external, unregular meshes

  implicit none

  include 'constants.h'
  include 'OUTPUT_FILES/values_from_mesher.h'

! comment next line if using old basin version
  integer :: NSPEC_AB, NGLOB_AB
  
! parameters  
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: data
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore

  integer, dimension(:,:,:,:),allocatable :: ibool
  logical, dimension(:),allocatable :: mask_ibool
  integer,dimension(:),allocatable :: num_ibool
  real,dimension(:,:,:,:),allocatable :: dat
    
  integer :: numpoin
  integer :: i, ios, it
  integer :: iproc, proc1, proc2, num_node, node_list(300), nspec, nglob
  integer :: np, ne, npp, nee, npoint, nelement, njunk
    
  character(len=150) :: sline, arg(6), filename, indir, outdir, prname
  character(len=150) :: mesh_file,local_data_file, local_ibool_file
  logical :: HIGH_RESOLUTION_MESH
  integer :: ires

! checks given arguments
  print *
  print *,'Recombining ParaView data for slices'
  print *

  do i = 1, 6
    call getarg(i,arg(i))
    if (i < 6 .and. trim(arg(i)) == '') then
      print *, 'Usage: '
      print *, '        xcombine_data start_slice end_slice filename input_dir output_dir high/low-resolution'
      print *, '    or '
      print *, '        xcombine_data slice_list filename input_dir output_dir high/low-resolution'
      print *, ' possible filenames are '
      print *, '   rho_vp, rho_vs, kappastore, mustore etc'
      print *, '   that are stored in the local directory as real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,nspec)  '
      print *, '   in filename.bin'
      print *, ' files have been collected in input_dir, output mesh file goes to output_dir '
      print *, ' give 0 for low resolution and 1 for high resolution'
      stop ' Reenter command line options'
    endif
  enddo

! get slice list
  if (trim(arg(6)) == '') then
    num_node = 0
    open(unit = 20, file = trim(arg(1)), status = 'unknown',iostat = ios)
    if (ios /= 0) then
      print *,'Error opening ',trim(arg(1))
      stop
    endif
    do while ( 1 == 1)
      read(20,'(a)',iostat=ios) sline
      if (ios /= 0) exit
      read(sline,*,iostat=ios) njunk
      if (ios /= 0) exit
      num_node = num_node + 1
      node_list(num_node) = njunk
    enddo
    close(20)
    filename = arg(2)
    indir= arg(3)
    outdir = arg(4)
    read(arg(5),*) ires
  else
    read(arg(1),*) proc1
    read(arg(2),*) proc2
    do iproc = proc1, proc2
      node_list(iproc - proc1 + 1) = iproc
    enddo
    num_node = proc2 - proc1 + 1
    filename = arg(3)
    indir = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
  endif

  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
  else
    HIGH_RESOLUTION_MESH = .true.
  endif

  print *, 'Slice list: '
  print *, node_list(1:num_node)

  ! open paraview output mesh file
  mesh_file = trim(outdir) // '/' // trim(filename)//'.mesh'
  call open_file(trim(mesh_file)//char(0))
    
  ! counts total number of points (all slices)
  npp = 0
  nee = 0
  call combine_vol_data_count_totals_ext_mesh(num_node,node_list,indir,npp,nee,HIGH_RESOLUTION_MESH)    


  ! write point and scalar information  
  np = 0 
  do it = 1, num_node

    iproc = node_list(it)

    print *, ' '
    print *, 'Reading slice ', iproc
    write(prname,'(a,i6.6,a)') trim(indir)//'/proc',iproc,'_'

    open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',action='read',form='unformatted')
    read(27) NSPEC_AB
    read(27) NGLOB_AB 
    close(27)   
    nspec = NSPEC_AB
    nglob = NGLOB_AB

    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(mask_ibool(NGLOB_AB))
    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB)) 
    
    ! data file  
    local_data_file = trim(prname) // trim(filename) // '.bin'
    open(unit = 27,file = trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
    if (ios /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    read(27) data
    close(27)
    print *, trim(local_data_file)

    ! uses implicit conversion to real values
    dat = data

    ! ibool file
    local_ibool_file = trim(prname) // 'ibool' // '.bin'
    open(unit = 28,file = trim(local_ibool_file),status='old',action='read', iostat = ios, form='unformatted')
    if (ios /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    read(28) ibool
    close(28)
    print *, trim(local_ibool_file)

    ! writes point coordinates and scalar value to mesh file
    if (.not. HIGH_RESOLUTION_MESH) then
      ! writes out element corners only
      call combine_vol_data_write_corners(nspec,nglob,ibool,mask_ibool,&
                                            xstore,ystore,zstore,dat,npoint,&
                                            it,npp,num_node,prname,numpoin)
    else  
      ! high resolution, all GLL points
      call combine_vol_data_write_GLL_points(nspec,nglob,ibool,mask_ibool,&
                                            xstore,ystore,zstore,dat,&
                                            it,npp,prname,numpoin)
    endif
    
    print*,'  points:',np,numpoin
    
    ! stores total number of points written
    np = np + numpoin

    deallocate(ibool,mask_ibool,data,dat,xstore,ystore,zstore)
    
  enddo  ! all slices for points

  if (np /=  npp) stop 'Error: Number of total points are not consistent'
  print *, 'Total number of points: ', np
  print *, ' '


! write element information
  ne = 0
  np = 0
  do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc
    write(prname,'(a,i6.6,a)') trim(indir)//'/proc',iproc,'_'

    open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',action='read',form='unformatted')
    read(27) NSPEC_AB
    read(27) NGLOB_AB 
    close(27)   
    nspec = NSPEC_AB
    nglob = NGLOB_AB
    
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(mask_ibool(NGLOB_AB))
    allocate(num_ibool(NGLOB_AB))

    ! ibool file
    local_ibool_file = trim(prname) // 'ibool' // '.bin'
    open(unit = 28,file = trim(local_ibool_file),status='old',action='read', iostat = ios, form='unformatted')
    if (ios /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    read(28) ibool
    close(28)

    if (.not. HIGH_RESOLUTION_MESH) then
      ! spectral elements
      call combine_vol_data_write_corner_elements(nspec,nglob,ibool,mask_ibool,num_ibool, &
                                            np,nelement, &
                                            it,nee,num_node,prname,numpoin)  
    else 
      ! subdivided spectral elements
      call combine_vol_data_write_GLL_elements(nspec,nglob,ibool,mask_ibool,num_ibool, &
                                            np,nelement, &
                                            it,nee,numpoin)  
    endif
    
    print*,'  elements:',ne,nelement
    print*,'  points : ',np,numpoin
    
    ne = ne + nelement

    deallocate(ibool,mask_ibool,num_ibool)

  enddo ! num_node
  
  ! checks with total number of elements
  if (ne /= nee) then 
    print*,'error: number of elements counted:',ne,'total:',nee
    stop 'Number of total elements are not consistent'
  endif
  print *, 'Total number of elements: ', ne

  ! close mesh file
  call close_file()

  print *, 'Done writing '//trim(mesh_file)

  end program combine_paraview_data


!=============================================================

! counts total number of points and elements for external meshes in given slice list

  subroutine combine_vol_data_count_totals_ext_mesh(num_node,node_list,indir,npp,nee,HIGH_RESOLUTION_MESH)

  implicit none
  include 'constants.h'
  
  integer,intent(in) :: num_node,node_list(300)
  character(len=150),intent(in) :: indir
  integer,intent(out) :: npp,nee
  logical,intent(in) :: HIGH_RESOLUTION_MESH
  
  ! local parameters
  integer, dimension(:,:,:,:),allocatable :: ibool
  logical, dimension(:),allocatable :: mask_ibool
  integer :: NSPEC_AB, NGLOB_AB
  integer :: it,iproc,npoint,nelement,ios,ispec
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  character(len=150) :: prname
  
  npp = 0
  nee = 0
  do it = 1, num_node
    ! gets number of elements and points for this slice
    iproc = node_list(it)
    write(prname,'(a,i6.6,a)') trim(indir)//'/proc',iproc,'_'
    open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',action='read',form='unformatted')
    read(27) NSPEC_AB
    read(27) NGLOB_AB 
    close(27)   
    
    if( HIGH_RESOLUTION_MESH ) then
      ! total number of global points
      npp = npp + NGLOB_AB

      ! total number of elements
      ! each spectral elements gets subdivided by GLL points, which form (NGLLX-1)*(NGLLY-1)*(NGLLZ-1) sub-elements
      nelement = NSPEC_AB * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1) 
      nee = nee + nelement

    else
      ! counts element corners only
      allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(mask_ibool(NGLOB_AB))

      ! ibool file
      open(unit = 28,file = prname(1:len_trim(prname))//'ibool'//'.bin',status='old',action='read',&
            iostat = ios,form='unformatted')
      if (ios /= 0) then
        print *,'Error opening: ',prname(1:len_trim(prname))//'ibool'//'.bin'
        stop
      endif
      read(28) ibool
      close(28)

      ! mark element corners (global AVS or DX points)
      mask_ibool = .false.
      do ispec=1,NSPEC_AB
        iglob1=ibool(1,1,1,ispec)
        iglob2=ibool(NGLLX,1,1,ispec)
        iglob3=ibool(NGLLX,NGLLY,1,ispec)
        iglob4=ibool(1,NGLLY,1,ispec)
        iglob5=ibool(1,1,NGLLZ,ispec)
        iglob6=ibool(NGLLX,1,NGLLZ,ispec)
        iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
        iglob8=ibool(1,NGLLY,NGLLZ,ispec)
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob5) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob8) = .true.
      enddo        

      ! count global number of AVS or DX points
      npoint = count(mask_ibool(:))
      npp = npp + npoint
      
      ! total number of spectral elements
      nee = nee + NSPEC_AB

      deallocate(ibool,mask_ibool)
    endif ! HIGH_RESOLUTION_MESH      
  enddo
  
  
  end subroutine combine_vol_data_count_totals_ext_mesh
  
!=============================================================

! writes out locations of spectral element corners only

  subroutine combine_vol_data_write_corners(nspec,nglob,ibool,mask_ibool,&
                                            xstore,ystore,zstore,dat,npoint,&
                                            it,npp,num_node,prname,numpoin)

  implicit none
  include 'constants.h'
  
  integer,intent(in) :: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  logical,dimension(nglob) :: mask_ibool
  real(kind=CUSTOM_REAL),dimension(nglob) :: xstore, ystore, zstore
  real,dimension(NGLLY,NGLLY,NGLLZ,nspec),intent(in) :: dat
  integer:: it  
  integer :: npp,num_node,npoint,numpoin
  character(len=150) :: prname

  ! local parameters
  real :: x, y, z
  integer :: ios,ispec,njunk
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  character(len=150) :: local_file

! corner locations  
  ! reads in coordinate files
  local_file = trim(prname)//'x.bin'
  open(unit = 27,file = trim(prname)//'x.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) then
    print *,'Error opening ',trim(local_file)
    stop
  endif
  read(27) xstore
  close(27)
  local_file = trim(prname)//'y.bin'
  open(unit = 27,file = trim(prname)//'y.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) then
    print *,'Error opening ',trim(local_file)
    stop
  endif
  read(27) ystore
  close(27)
  local_file = trim(prname)//'z.bin'
  open(unit = 27,file = trim(prname)//'z.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) then
    print *,'Error opening ',trim(local_file)
    stop
  endif
  read(27) zstore
  close(27)  

  ! writes out total number of points
  if (it == 1) then
    call write_integer(npp)
  endif

  mask_ibool(:) = .false.
  numpoin = 0

  do ispec=1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)

    if(.not. mask_ibool(iglob1)) then
      numpoin = numpoin + 1
      x = xstore(iglob1)
      y = ystore(iglob1)
      z = zstore(iglob1)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(1,1,1,ispec))
      mask_ibool(iglob1) = .true.
    endif
    if(.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      x = xstore(iglob2)
      y = ystore(iglob2)
      z = zstore(iglob2)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(NGLLX,1,1,ispec))
      mask_ibool(iglob2) = .true.
    endif
    if(.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      x = xstore(iglob3)
      y = ystore(iglob3)
      z = zstore(iglob3)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(NGLLX,NGLLY,1,ispec))
      mask_ibool(iglob3) = .true.
    endif
    if(.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      x = xstore(iglob4)
      y = ystore(iglob4)
      z = zstore(iglob4)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(1,NGLLY,1,ispec))
      mask_ibool(iglob4) = .true.
    endif
    if(.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      x = xstore(iglob5)
      y = ystore(iglob5)
      z = zstore(iglob5)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(1,1,NGLLZ,ispec))
      mask_ibool(iglob5) = .true.
    endif
    if(.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      x = xstore(iglob6)
      y = ystore(iglob6)
      z = zstore(iglob6)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(NGLLX,1,NGLLZ,ispec))
      mask_ibool(iglob6) = .true.
    endif
    if(.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      x = xstore(iglob7)
      y = ystore(iglob7)
      z = zstore(iglob7)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(NGLLX,NGLLY,NGLLZ,ispec))
      mask_ibool(iglob7) = .true.
    endif
    if(.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      x = xstore(iglob8)
      y = ystore(iglob8)
      z = zstore(iglob8)
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(dat(1,NGLLY,NGLLZ,ispec))
      mask_ibool(iglob8) = .true.
    endif
  enddo ! ispec
    
  end subroutine combine_vol_data_write_corners


!=============================================================

! writes out locations of all GLL points of spectral elements

  subroutine combine_vol_data_write_GLL_points(nspec,nglob,ibool,mask_ibool,&
                                            xstore,ystore,zstore,dat,&
                                            it,npp,prname,numpoin)

  implicit none
  include 'constants.h'
  
  integer,intent(in) :: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  logical,dimension(nglob) :: mask_ibool
  real(kind=CUSTOM_REAL),dimension(nglob) :: xstore, ystore, zstore
  real,dimension(NGLLY,NGLLY,NGLLZ,nspec),intent(in) :: dat
  integer:: it,npp,numpoin
  character(len=150) :: prname

  ! local parameters
  real :: x, y, z
  integer :: ios,ispec,i,j,k,iglob
  character(len=150) :: local_file

  ! writes out total number of points
  if (it == 1) then
    !npoint = nglob
    call write_integer(npp)
  endif

  ! reads in coordinate files
  local_file = trim(prname)//'x.bin'
  open(unit = 27,file = trim(prname)//'x.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) then
    print *,'Error opening ',trim(local_file)
    stop
  endif
  read(27) xstore
  close(27)
  local_file = trim(prname)//'y.bin'
  open(unit = 27,file = trim(prname)//'y.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) then
    print *,'Error opening ',trim(local_file)
    stop
  endif
  read(27) ystore
  close(27)
  local_file = trim(prname)//'z.bin'
  open(unit = 27,file = trim(prname)//'z.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) then
    print *,'Error opening ',trim(local_file)
    stop
  endif
  read(27) zstore
  close(27)

  mask_ibool(:) = .false.
  numpoin = 0

  do ispec=1,nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if(.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            x = xstore(iglob)
            y = ystore(iglob)
            z = zstore(iglob)
            call write_real(x)
            call write_real(y)
            call write_real(z)
            call write_real(dat(i,j,k,ispec))
            mask_ibool(iglob) = .true.
          endif
        enddo ! i
      enddo ! j
    enddo ! k
  enddo !ispec

  end subroutine combine_vol_data_write_GLL_points
  
!=============================================================

! writes out locations of spectral element corners only

  subroutine combine_vol_data_write_corner_elements(nspec,nglob,ibool,mask_ibool,num_ibool,&
                                            np,nelement,&
                                            it,nee,num_node,prname,numpoin)

  implicit none
  include 'constants.h'
  
  integer,intent(in) :: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  logical,dimension(nglob) :: mask_ibool
  integer,dimension(nglob) :: num_ibool
  integer:: it,nee,num_node,np,nelement,numpoin
  character(len=150) :: prname

  ! local parameters
  integer :: ios,ispec,i
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  integer :: njunk, njunk2, n1, n2, n3, n4, n5, n6, n7, n8  
  character(len=150) :: local_element_file


  ! outputs total number of elements for all slices
  if (it == 1) then
    call write_integer(nee)
  end if

  num_ibool(:) = 0
  mask_ibool(:) = .false.
  numpoin = 0
  
  do ispec=1,nspec
    ! gets corner indices
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)

    ! sets increasing numbering
    if(.not. mask_ibool(iglob1)) then
      numpoin = numpoin + 1
      num_ibool(iglob1) = numpoin
      mask_ibool(iglob1) = .true.          
    endif
    if(.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      num_ibool(iglob2) = numpoin
      mask_ibool(iglob2) = .true.        
    endif
    if(.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      num_ibool(iglob3) = numpoin
      mask_ibool(iglob3) = .true.        
    endif
    if(.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      num_ibool(iglob4) = numpoin
      mask_ibool(iglob4) = .true.        
    endif
    if(.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      num_ibool(iglob5) = numpoin
      mask_ibool(iglob5) = .true.        
    endif
    if(.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      num_ibool(iglob6) = numpoin
      mask_ibool(iglob6) = .true.        
    endif
    if(.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      num_ibool(iglob7) = numpoin
      mask_ibool(iglob7) = .true.        
    endif
    if(.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      num_ibool(iglob8) = numpoin
      mask_ibool(iglob8) = .true.        
    endif
  
    ! outputs corner indices (starting with 0 )
    n1 = num_ibool(iglob1) -1 + np 
    n2 = num_ibool(iglob2) -1 + np 
    n3 = num_ibool(iglob3) -1 + np 
    n4 = num_ibool(iglob4) -1 + np 
    n5 = num_ibool(iglob5) -1 + np 
    n6 = num_ibool(iglob6) -1 + np 
    n7 = num_ibool(iglob7) -1 + np 
    n8 = num_ibool(iglob8) -1 + np 
    
    call write_integer(n1)
    call write_integer(n2)
    call write_integer(n3)
    call write_integer(n4)
    call write_integer(n5)
    call write_integer(n6)
    call write_integer(n7)
    call write_integer(n8)

  enddo

  ! elements written
  nelement = nspec
  
  ! updates points written
  np = np + numpoin
    
  end subroutine combine_vol_data_write_corner_elements
  
  
!=============================================================

! writes out locations of spectral element corners only

  subroutine combine_vol_data_write_GLL_elements(nspec,nglob,ibool,mask_ibool,num_ibool,&
                                            np,nelement,&
                                            it,nee,numpoin)

  implicit none
  include 'constants.h'
  
  integer,intent(in):: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  logical,dimension(nglob) :: mask_ibool
  integer,dimension(nglob) :: num_ibool
  integer:: it,nee,np,numpoin,nelement

  ! local parameters
  integer :: ispec,i,j,k
  integer :: iglob,iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  integer :: n1, n2, n3, n4, n5, n6, n7, n8
  
  ! outputs total number of elements for all slices
  if (it == 1) then
    !nee = nelement * num_node
    call write_integer(nee)
  endif

  numpoin = 0
  mask_ibool = .false.
  
  ! sets numbering num_ibool respecting mask
  do ispec=1,nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if(.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            num_ibool(iglob) = numpoin
            mask_ibool(iglob) = .true.
          endif
        enddo ! i
      enddo ! j
    enddo ! k
  enddo !ispec

  ! outputs GLL subelement
  do ispec = 1, nspec
    do k = 1, NGLLZ-1
      do j = 1, NGLLY-1
        do i = 1, NGLLX-1
          iglob1 = ibool(i,j,k,ispec)
          iglob2 = ibool(i+1,j,k,ispec)
          iglob3 = ibool(i+1,j+1,k,ispec)
          iglob4 = ibool(i,j+1,k,ispec)
          iglob5 = ibool(i,j,k+1,ispec)
          iglob6 = ibool(i+1,j,k+1,ispec)
          iglob7 = ibool(i+1,j+1,k+1,ispec)
          iglob8 = ibool(i,j+1,k+1,ispec)
          n1 = num_ibool(iglob1)+np-1
          n2 = num_ibool(iglob2)+np-1
          n3 = num_ibool(iglob3)+np-1
          n4 = num_ibool(iglob4)+np-1
          n5 = num_ibool(iglob5)+np-1
          n6 = num_ibool(iglob6)+np-1
          n7 = num_ibool(iglob7)+np-1
          n8 = num_ibool(iglob8)+np-1
          call write_integer(n1)
          call write_integer(n2)
          call write_integer(n3)
          call write_integer(n4)
          call write_integer(n5)
          call write_integer(n6)
          call write_integer(n7)
          call write_integer(n8)
        enddo
      enddo
    enddo
  enddo
  ! elements written
  nelement = nspec * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1) 
  
  ! updates points written
  np = np + numpoin

  end subroutine combine_vol_data_write_GLL_elements