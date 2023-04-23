!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

  program combine_vol_data_vtk_binary

! puts the output of SPECFEM3D into '***.mesh' format,
! which can be converted via mesh2vtu into ParaView format.
!
! for Paraview, see http://www.paraview.org for details
!
! combines the database files on several slices.
! the local database file needs to have been collected onto the frontend (copy_local_database.pl)
!
! works for external, unregular meshes

  use constants

  use shared_parameters

  use combine_vol_data_mod

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! not supported yet with adios
  !use combine_vol_data_adios_mod
#endif

  use combine_vtk_par
  implicit none

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  double precision,dimension(:,:,:,:),allocatable :: data_dp
  ! real array for data
  real,dimension(:,:,:,:),allocatable :: data_sp

  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  real,dimension(:),allocatable :: pts
  integer, dimension(:,:,:,:),allocatable :: ibool

  integer :: NSPEC_AB, NGLOB_AB, NSPEC_IRREGULAR
  integer :: numpoin

  integer :: i, it, ier
  integer :: iproc, num_node

  integer,dimension(MAX_NUM_NODES) :: node_list
  integer,dimension(:), allocatable :: celltype
  integer,dimension(:,:), allocatable :: conn
  integer :: np, ne, npp, nee, nelement

  character(len=MAX_STRING_LEN) :: arg(9), filename, indir, outdir
  character(len=MAX_STRING_LEN) :: prname, prname_lp
  character(len=MAX_STRING_LEN*2) :: mesh_file,local_data_file
  logical :: HIGH_RESOLUTION_MESH,BROADCAST_AFTER_READ
  integer :: ires
  integer :: sizeprocs

  ! MPI initialization
  call init_mpi()
  call world_size(sizeprocs)
  if (sizeprocs /= 1) then
    print *, "sequential program. Only mpirun -np 1 ..."
    call abort_mpi()
  endif

  print *
  print *,'Recombining ParaView data for slices'
  print *

  ! needs local_path for mesh files
  myrank = 0
  BROADCAST_AFTER_READ = .false.
  call read_parameter_file(BROADCAST_AFTER_READ)

  ! reads in arguments
  do i = 1, command_argument_count()
    call get_command_argument(i,arg(i))
  enddo

  call read_args(arg, MAX_NUM_NODES, node_list, num_node, filename, indir, outdir, ires, NPROC)

  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
  else
    HIGH_RESOLUTION_MESH = .true.
  endif

  print *, 'Slice list: '
  print *, node_list(1:num_node)

  mesh_file = trim(outdir) // '/' // trim(filename)//'.vtk'

  ! counts total number of points (all slices)
  npp = 0
  nee = 0

  call cvd_count_totals_ext_mesh(num_node,node_list,LOCAL_PATH, &
                                 npp,nee,HIGH_RESOLUTION_MESH)

  ! writes point and scalar information
  ! loops over slices (process partitions)
  np = 0
  allocate(pts(3*npp),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1156')
  do it = 1, num_node

    iproc = node_list(it)
    print *, ' '
    print *, 'Reading slice ', iproc

    write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'
    open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    read(27) NSPEC_AB
    read(27) NGLOB_AB
    read(27) NSPEC_IRREGULAR

    ! ibool and global point arrays file
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1157')
    if (ier /= 0) stop 'error allocating array ibool'
    allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1158')
    if (ier /= 0) stop 'error allocating array xstore etc.'

    read(27) ibool
    read(27) xstore
    read(27) ystore
    read(27) zstore
    close(27)

    allocate(data_sp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1159')
    if (ier /= 0) stop 'error allocating single precision data array'

    if (CUSTOM_REAL == SIZE_DOUBLE) then
      allocate(data_dp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1160')
      if (ier /= 0) stop 'error allocating double precision data array'
    endif

    ! data file
    write(prname,'(a,i6.6,a)') trim(indir)//'proc',iproc,'_'
    local_data_file = trim(prname) // trim(filename) // '.bin'
    open(unit = 28,file = trim(local_data_file),status='old', &
          action='read',form ='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif

    ! Read either SP or DP floating point numbers.
    if (CUSTOM_REAL == SIZE_DOUBLE) then
      read(28) data_dp
    else
      read(28) data_sp
    endif
    close(28)
    print *, trim(local_data_file)

    ! uses conversion to real values
    if (CUSTOM_REAL == SIZE_DOUBLE) then
      data_sp(:,:,:,:) = sngl(data_dp(:,:,:,:))
      deallocate(data_dp)
    endif

    print *, trim(local_data_file)

    ! writes point coordinates and scalar value to mesh file
    if (.not. HIGH_RESOLUTION_MESH) then
      ! writes out element corners only
      call cvd_write_corners_binary(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,data_sp, &
                                    it,npp,numpoin,np,pts)
    else
      ! high resolution, all GLL points
      call cvd_write_GLL_points_binary(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,data_sp, &
                                       it,npp,numpoin,np,pts)
    endif

    print *,'  points:',np,numpoin

    ! stores total number of points written
    np = np + numpoin

    ! cleans up memory allocations
    deallocate(ibool,xstore,ystore,zstore)
    deallocate(data_sp)

  enddo  ! all slices for points

  if (np /= npp) stop 'Error: Number of total points are not consistent'
  print *, 'Total number of points: ', np
  print *, ' '

! writes element information
  ne = 0
  np = 0

  allocate(conn(8,nee),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1161')

  do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc

    write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'

    ! gets number of elements and global points for this partition
    open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
          status='old',action='read',form='unformatted')
    read(27) NSPEC_AB
    read(27) NGLOB_AB
    read(27) NSPEC_IRREGULAR

    ! ibool file
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1162')
    if (ier /= 0) stop 'error allocating array ibool'
    read(27) ibool
    close(27)

    ! writes out element corner indices
    if (.not. HIGH_RESOLUTION_MESH) then
      ! spectral elements
      call cvd_write_corner_elements_binary(NSPEC_AB,NGLOB_AB,ibool, &
                                            np,ne,nelement,it,nee,numpoin,conn)
    else
      ! subdivided spectral elements
      call cvd_write_GLL_elements_binary(NSPEC_AB,NGLOB_AB,ibool, &
                                         np,ne,nelement,it,nee,numpoin,conn)
    endif

    print *,'  elements:',ne,nelement
    print *,'  points : ',np,numpoin

    ne = ne + nelement

    deallocate(ibool)

  enddo ! num_node

  ! checks with total number of elements
  if (ne /= nee) then
    print *,'error: number of elements counted:',ne,'total:',nee
    stop 'Number of total elements are not consistent'
  endif
  print *, 'Total number of elements: ', ne

  allocate(celltype(nee),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1163')
  celltype=12

  call write_unstructured_mesh(mesh_file,len_trim(mesh_file), 1, npp, pts, nee, celltype, conn, &
                               filename,len_trim(filename),total_dat)

  call finalize_mpi()

  print *, 'Done writing '//trim(mesh_file)

  end program combine_vol_data_vtk_binary


!=============================================================

  subroutine cvd_count_totals_ext_mesh(num_node,node_list,LOCAL_PATH,npp,nee,HIGH_RESOLUTION_MESH)

! counts total number of points and elements for external meshes in given slice list
! returns: total number of elements (nee) and number of points (npp)

  use constants

  use combine_vtk_par

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! not supported yet with adios
  !use combine_vol_data_adios_mod
#endif

  implicit none

  integer,intent(in) :: num_node
  integer,dimension(MAX_NUM_NODES),intent(in) :: node_list

  integer,intent(out) :: npp,nee
  logical,intent(in) :: HIGH_RESOLUTION_MESH
  character(len=MAX_STRING_LEN),intent(in) :: LOCAL_PATH

  ! local parameters
  integer, dimension(:,:,:,:),allocatable :: ibool
  logical, dimension(:),allocatable :: mask_ibool
  integer :: NSPEC_AB, NGLOB_AB, NSPEC_IRREGULAR
  integer :: it,iproc,npoint,nelement,ispec,ier
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  character(len=MAX_STRING_LEN) :: prname_lp

  ! loops over all slices (process partitions)
  npp = 0
  nee = 0

  do it = 1, num_node

    ! gets number of elements and points for this slice
    iproc = node_list(it)
    write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'
    open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',prname_lp(1:len_trim(prname_lp))//'external_mesh.bin'
      stop
    endif

    read(27) NSPEC_AB
    read(27) NGLOB_AB
    read(27) NSPEC_IRREGULAR

    ! gets ibool
    if (.not. HIGH_RESOLUTION_MESH) then
      allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1164')
      if (ier /= 0) stop 'error allocating array ibool'
      read(27) ibool
    endif

    close(27)

    ! calculates totals
    if (HIGH_RESOLUTION_MESH) then
      ! total number of global points
      npp = npp + NGLOB_AB

      ! total number of elements
      ! each spectral elements gets subdivided by GLL points,
      ! which form (NGLLX-1)*(NGLLY-1)*(NGLLZ-1) sub-elements
      nelement = NSPEC_AB * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
      nee = nee + nelement

    else

      ! mark element corners (global AVS or DX points)
      allocate(mask_ibool(NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1165')
      if (ier /= 0) stop 'error allocating array mask_ibool'
      mask_ibool = .false.
      do ispec = 1,NSPEC_AB
        iglob1 = ibool(1,1,1,ispec)
        iglob2 = ibool(NGLLX,1,1,ispec)
        iglob3 = ibool(NGLLX,NGLLY,1,ispec)
        iglob4 = ibool(1,NGLLY,1,ispec)
        iglob5 = ibool(1,1,NGLLZ,ispec)
        iglob6 = ibool(NGLLX,1,NGLLZ,ispec)
        iglob7 = ibool(NGLLX,NGLLY,NGLLZ,ispec)
        iglob8 = ibool(1,NGLLY,NGLLZ,ispec)
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

    endif ! HIGH_RESOLUTION_MESH

    ! frees arrays
    if (allocated(mask_ibool)) deallocate( mask_ibool)
    if (allocated(ibool)) deallocate(ibool)

  enddo

  end subroutine cvd_count_totals_ext_mesh

!=============================================================


  subroutine cvd_write_corners_binary(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,dat, &
                                      it,npp,numpoin,np,pts)

! writes out locations of spectral element corners only

  use constants
  use combine_vtk_par

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore, ystore, zstore
  real,dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: dat
  integer,intent(in):: it
  integer,intent(in) :: npp,np
  integer,intent(inout) :: numpoin
  real,dimension(3*npp),intent(inout) :: pts

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  real :: x, y, z
  integer :: ispec,ier
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8

  ! writes out total number of points
  if (it == 1) then
    allocate(total_dat(npp),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1166')
    if (ier /= 0) stop 'error allocating total dat array'
    total_dat(:) = 0.0
  endif

  ! writes our corner point locations
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1167')
  if (ier /= 0) stop 'error allocating array mask_ibool'
  mask_ibool(:) = .false.
  numpoin = 0
  do ispec = 1,NSPEC_AB
    iglob1 = ibool(1,1,1,ispec)
    iglob2 = ibool(NGLLX,1,1,ispec)
    iglob3 = ibool(NGLLX,NGLLY,1,ispec)
    iglob4 = ibool(1,NGLLY,1,ispec)
    iglob5 = ibool(1,1,NGLLZ,ispec)
    iglob6 = ibool(NGLLX,1,NGLLZ,ispec)
    iglob7 = ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8 = ibool(1,NGLLY,NGLLZ,ispec)

    if (.not. mask_ibool(iglob1)) then
      numpoin = numpoin + 1
      x = xstore(iglob1)
      y = ystore(iglob1)
      z = zstore(iglob1)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(1,1,1,ispec)
      mask_ibool(iglob1) = .true.
    endif
    if (.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      x = xstore(iglob2)
      y = ystore(iglob2)
      z = zstore(iglob2)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(NGLLX,1,1,ispec)
      mask_ibool(iglob2) = .true.
    endif
    if (.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      x = xstore(iglob3)
      y = ystore(iglob3)
      z = zstore(iglob3)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(NGLLX,NGLLY,1,ispec)
      mask_ibool(iglob3) = .true.
    endif
    if (.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      x = xstore(iglob4)
      y = ystore(iglob4)
      z = zstore(iglob4)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(1,NGLLY,1,ispec)
      mask_ibool(iglob4) = .true.
    endif
    if (.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      x = xstore(iglob5)
      y = ystore(iglob5)
      z = zstore(iglob5)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(1,1,NGLLZ,ispec)
      mask_ibool(iglob5) = .true.
    endif
    if (.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      x = xstore(iglob6)
      y = ystore(iglob6)
      z = zstore(iglob6)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(NGLLX,1,NGLLZ,ispec)
      mask_ibool(iglob6) = .true.
    endif
    if (.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      x = xstore(iglob7)
      y = ystore(iglob7)
      z = zstore(iglob7)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(NGLLX,NGLLY,NGLLZ,ispec)
      mask_ibool(iglob7) = .true.
    endif
    if (.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      x = xstore(iglob8)
      y = ystore(iglob8)
      pts(3*(np+numpoin-1)+1) = x
      pts(3*(np+numpoin-1)+2) = y
      pts(3*(np+numpoin-1)+3) = z
      total_dat(np+numpoin) = dat(1,NGLLY,NGLLZ,ispec)
      mask_ibool(iglob8) = .true.
    endif
  enddo ! ispec

  end subroutine cvd_write_corners_binary

!=============================================================

  subroutine cvd_write_GLL_points_binary(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,dat, &
                                         it,npp,numpoin,np,pts)

! writes out locations of all GLL points of spectral elements

  use constants
  use combine_vtk_par

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore, ystore, zstore
  real,dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: dat
  integer,intent(in) :: it,npp,np
  integer,intent(inout) :: numpoin
  real,dimension(3*npp),intent(inout) :: pts

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  real :: x, y, z
  integer :: ispec,i,j,k,iglob,ier

  ! writes out total number of points
  if (it == 1) then
    allocate(total_dat(npp),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1168')
    if (ier /= 0) stop 'error allocating total dat array'
    total_dat(:) = 0.0
  endif

  ! writes out point locations and values
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1169')
  if (ier /= 0) stop 'error allocating array mask_ibool'

  mask_ibool(:) = .false.
  numpoin = 0
  do ispec=1,NSPEC_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            x = xstore(iglob)
            y = ystore(iglob)
            z = zstore(iglob)
            pts(3*(np+numpoin-1)+1) = x
            pts(3*(np+numpoin-1)+2) = y
            pts(3*(np+numpoin-1)+3) = z
            total_dat(np+numpoin) = dat(i,j,k,ispec)
            mask_ibool(iglob) = .true.
          endif
        enddo ! i
      enddo ! j
    enddo ! k
  enddo !ispec

  end subroutine cvd_write_GLL_points_binary

!=============================================================

! writes out locations of spectral element corners only

  subroutine cvd_write_corner_elements_binary(NSPEC_AB,NGLOB_AB,ibool, &
                                              np,ne,nelement,it,nee,numpoin,conn)

  use constants
  use combine_vtk_par

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  integer,intent(in) :: it,nee,ne
  integer,intent(inout) :: np,nelement,numpoin
  integer,dimension(8,nee),intent(inout) :: conn

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  integer,dimension(:),allocatable :: num_ibool
  integer :: ispec,ier,n
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8

  ! outputs total number of elements for all slices
  if (it == 1) then
    ! VTK
    ! note: indices for vtk start at 0
    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nee,nee*9
  endif

  ! writes out element indices
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1170')
  allocate(num_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1171')
  if (ier /= 0) stop 'error allocating array mask_ibool'

  mask_ibool(:) = .false.
  num_ibool(:) = 0
  numpoin = 0
  n = 0
  do ispec = 1,NSPEC_AB
    ! gets corner indices
    iglob1 = ibool(1,1,1,ispec)
    iglob2 = ibool(NGLLX,1,1,ispec)
    iglob3 = ibool(NGLLX,NGLLY,1,ispec)
    iglob4 = ibool(1,NGLLY,1,ispec)
    iglob5 = ibool(1,1,NGLLZ,ispec)
    iglob6 = ibool(NGLLX,1,NGLLZ,ispec)
    iglob7 = ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8 = ibool(1,NGLLY,NGLLZ,ispec)

    ! sets increasing numbering
    if (.not. mask_ibool(iglob1)) then
      numpoin = numpoin + 1
      num_ibool(iglob1) = numpoin
      mask_ibool(iglob1) = .true.
    endif
    if (.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      num_ibool(iglob2) = numpoin
      mask_ibool(iglob2) = .true.
    endif
    if (.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      num_ibool(iglob3) = numpoin
      mask_ibool(iglob3) = .true.
    endif
    if (.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      num_ibool(iglob4) = numpoin
      mask_ibool(iglob4) = .true.
    endif
    if (.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      num_ibool(iglob5) = numpoin
      mask_ibool(iglob5) = .true.
    endif
    if (.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      num_ibool(iglob6) = numpoin
      mask_ibool(iglob6) = .true.
    endif
    if (.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      num_ibool(iglob7) = numpoin
      mask_ibool(iglob7) = .true.
    endif
    if (.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      num_ibool(iglob8) = numpoin
      mask_ibool(iglob8) = .true.
    endif
    n = n + 1

    ! outputs corner indices (starting with 0)
    conn(1,ne +n) = num_ibool(iglob1) -1 + np
    conn(2,ne +n) = num_ibool(iglob2) -1 + np
    conn(3,ne +n) = num_ibool(iglob3) -1 + np
    conn(4,ne +n) = num_ibool(iglob4) -1 + np
    conn(5,ne +n) = num_ibool(iglob5) -1 + np
    conn(6,ne +n) = num_ibool(iglob6) -1 + np
    conn(7,ne +n) = num_ibool(iglob7) -1 + np
    conn(8,ne +n) = num_ibool(iglob8) -1 + np

  enddo

  ! elements written
  nelement = NSPEC_AB

  ! updates points written
  np = np + numpoin

  end subroutine cvd_write_corner_elements_binary


!=============================================================


  subroutine cvd_write_GLL_elements_binary(NSPEC_AB,NGLOB_AB,ibool, &
                                           np,ne,nelement,it,nee,numpoin,conn)

! writes out indices of elements given by GLL points

  use constants
  use combine_vtk_par

  implicit none

  integer,intent(in):: NSPEC_AB,NGLOB_AB
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  integer,intent(in) :: it,nee,ne
  integer,intent(inout) :: np,numpoin,nelement
  integer,dimension(8,nee),intent(inout) :: conn

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  integer,dimension(:),allocatable :: num_ibool
  integer :: ispec,i,j,k,ier,n
  integer :: iglob,iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8

  ! outputs total number of elements for all slices
  if (it == 1) then
    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nee,nee*9
  endif

  ! sets numbering num_ibool respecting mask
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1172')
  allocate(num_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1173')
  if (ier /= 0) stop 'error allocating array mask_ibool'

  mask_ibool(:) = .false.
  num_ibool(:) = 0
  numpoin = 0
  do ispec = 1,NSPEC_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            num_ibool(iglob) = numpoin
            mask_ibool(iglob) = .true.
          endif
        enddo ! i
      enddo ! j
    enddo ! k
  enddo !ispec

  n = 0
  ! outputs GLL subelement
  do ispec = 1, NSPEC_AB
    do k = 1, NGLLZ-1
      do j = 1, NGLLY-1
        do i = 1, NGLLX-1
          n = n + 1
          iglob1 = ibool(i,j,k,ispec)
          iglob2 = ibool(i+1,j,k,ispec)
          iglob3 = ibool(i+1,j+1,k,ispec)
          iglob4 = ibool(i,j+1,k,ispec)
          iglob5 = ibool(i,j,k+1,ispec)
          iglob6 = ibool(i+1,j,k+1,ispec)
          iglob7 = ibool(i+1,j+1,k+1,ispec)
          iglob8 = ibool(i,j+1,k+1,ispec)
          conn(1,ne +n) = num_ibool(iglob1)+np-1
          conn(2,ne +n) = num_ibool(iglob2)+np-1
          conn(3,ne +n) = num_ibool(iglob3)+np-1
          conn(4,ne +n) = num_ibool(iglob4)+np-1
          conn(5,ne +n) = num_ibool(iglob5)+np-1
          conn(6,ne +n) = num_ibool(iglob6)+np-1
          conn(7,ne +n) = num_ibool(iglob7)+np-1
          conn(8,ne +n) = num_ibool(iglob8)+np-1
        enddo
      enddo
    enddo
  enddo

  ! elements written
  nelement = NSPEC_AB * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
  ! updates points written
  np = np + numpoin

  end subroutine cvd_write_GLL_elements_binary

