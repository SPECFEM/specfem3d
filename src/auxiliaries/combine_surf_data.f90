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

program combine_surf_data

  ! puts the output of SPECFEM3D in ParaView format.
  ! see http://www.paraview.org for details

  ! combines the database files on several slices.
  ! the local database file needs to have been collected onto the frontend (copy_local_database.pl)

  use constants

  implicit none

  integer :: i,j,k,ispec, it
  integer :: iproc, proc1, proc2, num_node, node_list(300), nspec, nglob
  integer :: np, ne, npp, nee, npoint, nelement, njunk, n1, n2, n3, n4
  integer :: numpoin, iglob1, iglob2, iglob3, iglob4, iglob
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: data_3D
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: data_2D
  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  integer, dimension(:,:,:,:),allocatable :: ibool
  logical, dimension(:),allocatable :: mask_ibool
  integer :: NSPEC_AB,NGLOB_AB

  real :: x, y, z
  real, dimension(:,:,:,:), allocatable :: dat3D
  real, dimension(:,:,:), allocatable :: dat2D
  character(len=MAX_STRING_LEN) :: sline, arg(8), filename, indir, outdir, prname, surfname
  character(len=MAX_STRING_LEN*2) :: mesh_file, local_file, local_data_file, local_ibool_file
  character(len=MAX_STRING_LEN*2) :: local_ibool_surf_file

!  integer :: num_ibool(NGLOB_AB)
  integer,dimension(:),allocatable :: num_ibool

  logical :: HIGH_RESOLUTION_MESH,  FILE_ARRAY_IS_3D
  integer :: ires, nspec_surf, npoint1, npoint2, ispec_surf, inx, iny, idimval, ier
  integer,dimension(:), allocatable ::  ibelm_surf


  do i = 1, 8
    call get_command_argument(i,arg(i))
    if (i < 6 .and. trim(arg(i)) == '') then
      print *, 'Usage: xcombine_surface start_slice end_slice filename surfacename input_dir output_dir high/low-resolution 3D/2D'
      print *, '    or xcombine_surface slice_list filename surfacename input_dir output_dir high/low-resolution 3D/2D'
      print *, ' possible filenames are kappastore(NGLLX,NGLLY,NGLLZ,nspec), alpha_kernel(NGLLX,NGLLY,nspec_surf)'
      print *, ' possible surface name: moho   as in ibelm_moho.bin'
      print *, ' files have been collected in input_dir, output mesh file goes to output_dir '
      print *, ' give 0 for low resolution and 1 for high resolution'
      print *, ' give 0 for 2D and 1 for 3D filenames'
      stop ' Reenter command line options'
    endif
  enddo

  ! get slice list
  if (trim(arg(8)) == '') then
    num_node = 0
    open(unit = 20, file = trim(arg(1)), status = 'unknown',iostat = ier)
    if (ier /= 0) then
      print *,'Error opening ',trim(arg(1))
      stop
    endif
    do while (1 == 1)
      read(20,'(a)',iostat=ier) sline
      if (ier /= 0) exit
      read(sline,*,iostat=ier) njunk
      if (ier /= 0) exit
      num_node = num_node + 1
      node_list(num_node) = njunk
    enddo
    close(20)
    filename = arg(2)
    surfname = arg(3)
    indir= arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
    read(arg(7),*) idimval
  else
    read(arg(1),*) proc1
    read(arg(2),*) proc2
    do iproc = proc1, proc2
      node_list(iproc - proc1 + 1) = iproc
    enddo
    num_node = proc2 - proc1 + 1
    filename = arg(3)
    surfname = arg(4)
    indir = arg(5)
    outdir = arg(6)
    read(arg(7),*) ires
    read(arg(8),*) idimval
  endif

  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
    inx = NGLLX-1
    iny = NGLLY-1
  else
    HIGH_RESOLUTION_MESH = .true.
    inx = 1
    iny = 1
  endif

  if (idimval == 0) then
    FILE_ARRAY_IS_3D = .false.
  else
    FILE_ARRAY_IS_3D = .true.
  endif

  print *, 'Slice list: '
  print *, node_list(1:num_node)

  ! open paraview output mesh file
  mesh_file = trim(outdir) // '/' // trim(filename)//'.surf'
  call open_file_create(trim(mesh_file)//char(0))

!  nspec = NSPEC_AB
!  nglob = NGLOB_AB

  np = 0

  ! =======  loop over all slices, write point and scalar information ======

  do it = 1, num_node

    iproc = node_list(it)

    print *, ' '
    print *, 'Reading slice ', iproc
    write(prname,'(a,i6.6,a)') trim(indir)//'/proc',iproc,'_'

    ! gets number of elements and global points for this partition
    open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    read(27) NSPEC_AB
    read(27) NGLOB_AB
    close(27)
    nspec = NSPEC_AB
    nglob = NGLOB_AB

    ! allocates arrays
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1107')
    allocate(mask_ibool(NGLOB_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1108')
    allocate(num_ibool(NGLOB_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1109')
    allocate(xstore(NGLOB_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1110')
    allocate(ystore(NGLOB_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1111')
    allocate(zstore(NGLOB_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1112')
    if (ier /= 0) stop 'error allocating array ibool etc.'

    ! surface file
    local_ibool_surf_file = trim(prname) // 'ibelm_' //trim(surfname)// '.bin'
    open(unit = 28,file = trim(local_ibool_surf_file),status='old', iostat = ier, form='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(local_ibool_surf_file)
      stop
    endif
    read(28) nspec_surf
    read(28) npoint1
    read(28) npoint2

    if (it == 1) then
      allocate(ibelm_surf(nspec_surf),stat=ier)
      if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1113')
      if (ier /= 0) stop 'error allocating array ibelm_surf'
    endif
    read(28) ibelm_surf
    close(28)
    print *, trim(local_ibool_surf_file)

    if (it == 1) then
      if (FILE_ARRAY_IS_3D) then
        allocate(data_3D(NGLLX,NGLLY,NGLLZ,NSPEC_AB),dat3D(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
        if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1114')
        if (ier /= 0) stop 'error allocating array data_3D'
      else
        allocate(data_2D(NGLLX,NGLLY,nspec_surf),dat2D(NGLLX,NGLLY,nspec_surf),stat=ier)
        if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1115')
        if (ier /= 0) stop 'error allocating array data_2D'
      endif
    endif

    ! data file
    local_data_file = trim(prname) // trim(filename) // '.bin'
    open(unit = 27,file = trim(local_data_file),status='old', iostat = ier,form ='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    if (FILE_ARRAY_IS_3D) then
      read(27) data_3D
      dat3D = data_3D
    else
      read(27) data_2D
      dat2D = data_2D
    endif
    close(27)
    print *, trim(local_data_file)

    ! ibool file
    local_ibool_file = trim(prname) // 'ibool' // '.bin'
    open(unit = 28,file = trim(local_ibool_file),status='old', iostat = ier, form='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    read(28) ibool
    close(28)
    print *, trim(local_ibool_file)


    mask_ibool(:) = .false.
    numpoin = 0

    if (it == 1) then
      if (HIGH_RESOLUTION_MESH) then
        npoint = npoint2
      else
        npoint = npoint1
      endif
      npp = npoint * num_node
      call write_integer(npp)
    endif

    local_file = trim(prname)//'x.bin'
    open(unit = 27,file = trim(prname)//'x.bin',status='old', iostat = ier,form ='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(local_file)
      stop
    endif
    read(27) xstore
    close(27)
    local_file = trim(prname)//'y.bin'
    open(unit = 27,file = trim(prname)//'y.bin',status='old', iostat = ier,form ='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(local_file)
      stop
    endif
    read(27) ystore
    close(27)
    local_file = trim(prname)//'z.bin'
    open(unit = 27,file = trim(prname)//'z.bin',status='old', iostat = ier,form ='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(local_file)
      stop
    endif
    read(27) zstore
    close(27)

    do ispec_surf=1,nspec_surf
      ispec = ibelm_surf(ispec_surf)
      k = 1
      do j = 1, NGLLY, iny
        do i = 1, NGLLX, inx
          iglob = ibool(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            x = xstore(iglob)
            y = ystore(iglob)
            z = zstore(iglob)
            call write_real(x)
            call write_real(y)
            call write_real(z)
            if (FILE_ARRAY_IS_3D) then
              call write_real(dat3D(i,j,k,ispec))
            else
              call write_real(dat2D(i,j,ispec_surf))
            endif
            mask_ibool(iglob) = .true.
          endif
        enddo ! i
      enddo ! j
    enddo !ispec

    if (numpoin /= npoint) stop 'Error: number of points are not consistent'
    np = np + npoint

    ! frees arrays
    deallocate(ibool,mask_ibool,num_ibool,xstore,ystore,zstore)

  enddo  ! all slices for points

  if (np /= npp) stop 'Error: Number of total points are not consistent'
  print *, 'Total number of points: ', np
  print *, ' '


  ne = 0
  ! ============  write element information =====================
  do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc
    write(prname,'(a,i6.6,a)') trim(indir)//'/proc',iproc,'_'

    ! gets number of elements and global points for this partition
    open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    read(27) NSPEC_AB
    read(27) NGLOB_AB
    close(27)
    nspec = NSPEC_AB
    nglob = NGLOB_AB

    ! allocates arrays
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1116')
    allocate(mask_ibool(NGLOB_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1117')
    allocate(num_ibool(NGLOB_AB),stat=ier)
    if (ier /= 0) call my_local_exit_MPI_without_rank('error allocating array 1118')
    if (ier /= 0) stop 'error allocating array ibool etc.'


    np = npoint * (it-1)

! surface file
    local_ibool_surf_file = trim(prname) // 'ibelm_' //trim(surfname)// '.bin'
    open(unit = 28,file = trim(local_ibool_surf_file),status='old', iostat = ier, form='unformatted')
    read(28) nspec_surf
    read(28) njunk
    read(28) njunk
    read(28) ibelm_surf
    close(28)

! ibool file
    local_ibool_file = trim(prname) // 'ibool' // '.bin'
    open(unit = 28,file = trim(local_ibool_file),status='old', iostat = ier, form='unformatted')
    read(28) ibool
    close(28)

    if (it == 1) then
      if (HIGH_RESOLUTION_MESH) then
        nelement = nspec_surf  * (NGLLX-1) * (NGLLY-1)
      else
        nelement = nspec_surf
      endif
      nee = nelement * num_node
      call write_integer(nee)
    endif

    numpoin = 0
    mask_ibool = .false.
    do ispec_surf=1,nspec_surf
      ispec = ibelm_surf(ispec_surf)
      k = 1
      do j = 1, NGLLY, iny
        do i = 1, NGLLX, inx
          iglob = ibool(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            num_ibool(iglob) = numpoin
            mask_ibool(iglob) = .true.
          endif
        enddo ! i
      enddo ! j
    enddo !ispec

    do ispec_surf = 1, nspec_surf
      ispec = ibelm_surf(ispec_surf)
      k = 1
      do j = 1, NGLLY-1, iny
        do i = 1, NGLLX-1, inx
          iglob1 = ibool(i,j,k,ispec)
          iglob2 = ibool(i+inx,j,k,ispec)
          iglob3 = ibool(i+inx,j+iny,k,ispec)
          iglob4 = ibool(i,j+iny,k,ispec)

          n1 = num_ibool(iglob1)+np-1
          n2 = num_ibool(iglob2)+np-1
          n3 = num_ibool(iglob3)+np-1
          n4 = num_ibool(iglob4)+np-1

          call write_integer(n1)
          call write_integer(n2)
          call write_integer(n3)
          call write_integer(n4)

        enddo
      enddo
    enddo
    ne = ne + nelement

    ! frees arrays
    deallocate(ibool,mask_ibool,num_ibool)

  enddo ! num_node
  if (ne /= nee) stop 'Number of total elements are not consistent'
  print *, 'Total number of elements: ', ne

  call close_file()

  print *, 'Done writing '//trim(mesh_file)

end program combine_surf_data

!
!-------------------------------------------------------------------------------------------------
!

! version without rank number printed in the error message

  subroutine my_local_exit_MPI_without_rank(error_msg)

  implicit none

  character(len=*) error_msg

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

  stop 'Fatal error'

  end subroutine my_local_exit_MPI_without_rank

