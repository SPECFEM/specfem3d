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

  program combine_vol_data

! combines mesh array data (NGLLX,NGLLY,NGLLZ,NSPEC) from different processes/ranks into a single file
!
! different format are possible:
!  xcombine_vol_data    : .mesh format - will create *.mesh files, which then can be converted to *.vtu files with the helper
!                                        mesh2vtu utility (in utils/Visualization/Paraview/mesh2vtu/ folder)
!
!  xcombine_vol_data_vtk: *.vtk format - will create *.vtk files (w/out a need of a VTK library installation)
!
!  xcombine_vol_data_vtu: *.vtu format - will create *.vtu binary files (w/out a need of a VTK library installation),
!                                        which have smaller file sizes as the ascii *.vtk files.

! puts the output of SPECFEM3D either into '***.mesh' format (which can further be converted via mesh2vtu)
! or directly into '***.vtk/***.vtu' format which can be visualized, e.g., by ParaView.
!
! for Paraview, see http://www.paraview.org for details
!
! combines the database files on several slices.
! the local database file needs to have been collected onto the frontend (copy_local_database.pl)
!
! works for external, irregular meshes

  use constants

  use shared_parameters

  use combine_vol_data_mod

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use combine_vol_data_adios_mod
#endif

  use combine_vtk_par

  implicit none

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: data

  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  integer, dimension(:,:,:,:),allocatable :: ibool

  integer :: NSPEC_AB, NGLOB_AB

  integer :: numpoin

  integer :: i, it, ier
  integer :: iproc, num_node

  integer,dimension(MAX_NUM_NODES) :: node_list

  integer :: np, ne, npp, nee, nelement

  character(len=MAX_STRING_LEN) :: arg(9)
  character(len=MAX_STRING_LEN) :: filename, outdir
  character(len=MAX_STRING_LEN*2) :: mesh_file,local_data_file
  character(len=MAX_STRING_LEN) :: data_array_name
  logical :: HIGH_RESOLUTION_MESH,BROADCAST_AFTER_READ
  integer :: ires

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! Variables to read ADIOS files
  integer :: sizeprocs
  character(len=MAX_STRING_LEN) :: var_name, value_file_name, mesh_file_name
#else
  integer :: NSPEC_IRREGULAR
  character(len=MAX_STRING_LEN) :: prname,prname_lp
  character(len=MAX_STRING_LEN) :: indir
#endif

  ! starts here---------------------------------------------------------------
  ier = 0 ! avoids compiler warning in case of ADIOS and VTK output

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS needs MPI initialization
  ! MPI initialization
  call init_mpi()
  call world_size(sizeprocs)
  ! checks number of processes
  ! note: must run as a single process with: mpirun -np 1 ..
  if (sizeprocs /= 1) then
    call world_rank(myrank)
    ! usage info
    if (myrank == 0) then
      print *, "ADIOS requires MPI functionality. However, this program executes as sequential program."
      print *, "Invalid number of processes used: ", sizeprocs, " procs"
      print *
      print *, "Please run: mpirun -np 1 ./bin/xcombine_vol_data_**_adios"
    endif
    call abort_mpi()
  endif
#endif

  ! user output
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

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! input
  call read_args_adios(arg, MAX_NUM_NODES, node_list, num_node, &
                       var_name, value_file_name, mesh_file_name, &
                       outdir, ires, NPROC)
  filename = var_name
#else
  ! default
  call read_args(arg, MAX_NUM_NODES, node_list, num_node, filename, indir, outdir, ires, NPROC)
#endif

  ! resolution
  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
  else
    HIGH_RESOLUTION_MESH = .true.
  endif

  ! output info
  print *, 'combine volumetric data'
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    print *, '  using ADIOS'
    print *, '  mesh topology file: ',trim(mesh_file_name)
    print *, '  input         file: ',trim(value_file_name)
    print *, '  variable name     : ',trim(var_name)
#else
    print *, '  mesh topology dir : ',trim(indir)
    print *, '  input file    dir : ',trim(filename)
#endif
  print *
  print *, '  output directory  : ',trim(outdir)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
#ifdef USE_VTK_INSTEAD_OF_MESH
  print *, '  using VTK format for file output'
#endif
#ifdef USE_VTU_INSTEAD_OF_MESH
  print *, '  using VTU format for file output'
#endif
#else
  ! .mesh format
  print *, '  using .mesh format for file output'
#endif
  print *
  print *, 'Slice list: '
  print *, node_list(1:num_node)
  print *

  ! checks
  if (num_node < 1) stop 'Error need at least one slice for combining data arrays, please check your input arguments...'

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! opens mesh and data/value file
  call init_adios(value_file_name, mesh_file_name)
#endif

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU format
  ! no open file required, will store data into arrays first, then write out to file
  !continue
#else
  ! .mesh format
  ! open paraview output mesh file
  mesh_file = trim(outdir) // '/' // trim(filename)//'.mesh'
  call open_file_create(trim(mesh_file)//char(0))
#endif

  ! counts total number of points (all slices)
  npp = 0
  nee = 0

  call cvd_count_totals_ext_mesh(num_node,node_list,LOCAL_PATH, &
                                 npp,nee,HIGH_RESOLUTION_MESH)

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
  print *
  print *,'VTK initial total points: ',npp
  print *,'VTK initial total elements: ',nee
  print *
  ! creates arrays to hold point data
  allocate(total_dat(npp), &
           total_dat_xyz(3,npp), &
           total_dat_con(8,nee),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat_con array'
  total_dat(:) = 0.0_CUSTOM_REAL
  total_dat_xyz(:,:) = 0.0_CUSTOM_REAL
  total_dat_con(:,:) = 0
#else
  ! .mesh
  ! writes out total number of points
  call write_integer(npp)
#endif

  ! writes point and scalar information
  ! loops over slices (process partitions)
  np = 0
  do it = 1, num_node

    iproc = node_list(it)
    print *, ' '
    print *, 'Reading slice ', iproc

    ! gets number of elements and global points for this partition
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    ! reads mesh nglob & nspec
    call read_scalars_adios_mesh(iproc, NGLOB_AB, NSPEC_AB)
#else
    ! default binary
    write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'
    open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    read(27) NSPEC_AB
    read(27) NGLOB_AB
    read(27) NSPEC_IRREGULAR
#endif

    ! ibool and global point arrays file
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1141')
    if (ier /= 0) stop 'error allocating array ibool'
    allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1142')
    if (ier /= 0) stop 'error allocating array xstore etc.'

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    call read_ibool_adios_mesh(iproc, NGLLX, NGLLY, NGLLZ, NSPEC_AB, ibool)
    call read_coordinates_adios_mesh(iproc, NGLOB_AB, xstore, ystore, zstore)
#else
    ! default binary
    read(27) ibool
    read(27) xstore
    read(27) ystore
    read(27) zstore
    close(27)
#endif

    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1143')
    if (ier /= 0) stop 'error allocating single precision data array'
    data(:,:,:,:) = 0.0_CUSTOM_REAL

    ! reads in kernel/data values
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    call read_values_adios(var_name, iproc, NSPEC_AB, data)
    local_data_file = trim(var_name)
#else
    ! default binary data file
    write(prname,'(a,i6.6,a)') trim(indir)// '/proc',iproc,'_'
    local_data_file = trim(prname) // trim(filename) // '.bin'
    open(unit = 28,file = trim(local_data_file),status='old',action='read',form ='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    ! Read in data
    read(28) data
    close(28)
#endif

    ! output info
    print *,'  file: ',trim(local_data_file)
    print *,'  min/max value: ',minval(data(:,:,:,:)),maxval(data(:,:,:,:))
    print *

    ! writes point coordinates and scalar value to mesh file
    if (.not. HIGH_RESOLUTION_MESH) then
      ! writes out element corners only
      call cvd_write_corners(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,data,numpoin,np)
    else
      ! high resolution, all GLL points
      call cvd_write_GLL_points(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,data,numpoin,np)
    endif

    print *,'  points:',np,numpoin

    ! checks integer overflow
    if (numpoin > 2147483646.d0 - np ) stop 'Error number of points might exceed integer limit'

    ! stores total number of points written
    np = np + numpoin

    ! cleans up memory allocations
    deallocate(ibool,xstore,ystore,zstore)
    deallocate(data)

  enddo  ! all slices for points

  if (np /= npp) stop 'Error: Number of total points are not consistent'

  print *
  print *, 'Total number of points: ', np
  print *, ' '

  ! outputs total number of elements for all slices
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
  ! continue
#else
  ! .mesh
  ! writes out total number of elements
  call write_integer(nee)
#endif

  ! writes element information
  ne = 0
  np = 0
  do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    call read_scalars_adios_mesh(iproc, NGLOB_AB, NSPEC_AB)
#else
    ! default binary
    write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'

    ! gets number of elements and global points for this partition
    open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin',status='old',action='read',form='unformatted')
    read(27) NSPEC_AB
    read(27) NGLOB_AB
    read(27) NSPEC_IRREGULAR
#endif

    ! ibool file
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1145')

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    call read_ibool_adios_mesh(iproc, NGLLX, NGLLY, NGLLZ, NSPEC_AB, ibool)
#else
    ! default binary
    if (ier /= 0) stop 'error allocating array ibool'
    read(27) ibool
    close(27)
#endif

    ! writes out element corner indices
    if (.not. HIGH_RESOLUTION_MESH) then
      ! spectral elements
      call cvd_write_corner_elements(NSPEC_AB,NGLOB_AB,ibool,np,ne,numpoin,nelement)
    else
      ! subdivided spectral elements
      call cvd_write_GLL_elements(NSPEC_AB,NGLOB_AB,ibool,np,ne,numpoin,nelement)
    endif

    print *,'  elements: ',ne,nelement
    print *,'  points  : ',np,numpoin

    deallocate(ibool)

  enddo ! num_node

  ! checks with total number of elements
  if (ne /= nee) then
    print *,'error: number of elements counted:',ne,'total:',nee
    stop 'Number of total elements are not consistent'
  endif

  print *
  print *, 'Total number of elements: ', ne
  print *

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU format
  ! point data values
  ! data array name
  data_array_name = trim(filename)
  ! movie snapshot files: uses a common short name (like "velocity_X","div","curl_X",..)
  ! note: changing the data array name for each timestep creates problems in paraview
  !       when reading in a bunch of these vtk-files and animating the wavefield
  ! velocity_X_it002000.bin,..,velocity_Y_it**.bin,velocity_Z_it**.bin
  if (filename(1:13) == 'velocity_X_it' .or. filename(1:13) == 'velocity_Y_it' .or. filename(1:13) == 'velocity_Z_it') then
    data_array_name = trim(filename(1:10)) ! "velocity_X",..
  endif
  if (filename(1:10) == 'displ_X_it' .or. filename(1:10) == 'displ_Y_it' .or. filename(1:10) == 'displ_Z_it') then
    data_array_name = trim(filename(1:10)) ! "displ_X",..
  endif
  if (filename(1:6) == 'div_it') then
    data_array_name = trim(filename(1:3)) ! "div"
  endif
  if (filename(1:9) == 'curl_X_it' .or. filename(1:9) == 'curl_Y_it' .or. filename(1:9) == 'curl_Z_it') then
    data_array_name = trim(filename(1:6)) ! "curl_X",..
  endif
  if (filename(1:11) == 'pressure_it' ) then
    data_array_name = trim(filename(1:8)) ! "pressure"
  endif

#ifdef USE_VTK_INSTEAD_OF_MESH
  ! default VTK
  ! outputs unstructured grid file
  mesh_file = trim(outdir) // '/' // trim(filename)//'.vtk'
  call write_VTK_movie_data(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,data_array_name)
#endif
#ifdef USE_VTU_INSTEAD_OF_MESH
  ! VTU binary format
  mesh_file = trim(outdir) // '/' // trim(filename)//'.vtu'
  call write_VTU_movie_data_binary(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,data_array_name)
#endif

  ! free arrays for this region
  deallocate(total_dat,total_dat_xyz,total_dat_con)
#else
  ! .mesh
  ! close mesh file
  call close_file()
  ! to avoid compiler warning
  data_array_name = ''
#endif

  ! cleanup
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  call clean_adios()
  ! shuts down mpi
  call finalize_mpi()
#endif

  print *
  print *, 'Done writing '//trim(mesh_file)
  print *

  end program combine_vol_data


!=============================================================

  subroutine cvd_count_totals_ext_mesh(num_node,node_list,LOCAL_PATH, &
                                       npp,nee,HIGH_RESOLUTION_MESH)

! counts total number of points and elements for external meshes in given slice list
! returns: total number of elements (nee) and number of points (npp)

  use constants

  use combine_vtk_par

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use combine_vol_data_adios_mod
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
  integer :: NSPEC_AB, NGLOB_AB
  integer :: it,iproc,npoint,nelement,ispec,ier
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! Variables for ADIOS
#else
  integer :: NSPEC_IRREGULAR
  character(len=MAX_STRING_LEN) :: prname_lp
#endif

  ! loops over all slices (process partitions)
  npp = 0
  nee = 0

  do it = 1, num_node

    iproc = node_list(it)

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    call read_scalars_adios_mesh(iproc, NGLOB_AB, NSPEC_AB)
    ! to avoid compiler warning
    ier = len_trim(LOCAL_PATH)
#else
    ! default
    ! gets number of elements and points for this slice
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
#endif

    ! gets ibool
    if (.not. HIGH_RESOLUTION_MESH) then
      allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1146')
      if (ier /= 0) stop 'error allocating array ibool'
#ifdef USE_ADIOS_INSTEAD_OF_MESH
      ! ADIOS
      call read_ibool_adios_mesh(iproc, NGLLX, NGLLY, NGLLZ, NSPEC_AB, ibool)
#else
      read(27) ibool
#endif
    endif

    ! closes file
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! nothing to do, continue
#else
    close(27)
#endif

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
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1147')
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


  subroutine cvd_write_corners(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,data,numpoin,np)

! writes out locations of spectral element corners only
  use constants
  use combine_vtk_par
  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: xstore, ystore, zstore
  real(kind=CUSTOM_REAL),dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: data

  integer, intent(inout) :: numpoin,np

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  real :: x, y, z, val
  integer :: ispec,ier
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8

  ! writes our corner point locations
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1149')
  if (ier /= 0) stop 'error allocating array mask_ibool'
  mask_ibool(:) = .false.

  ! dummy statement to avoid compiler warning
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
  ! nothing to do
#else
  ! dummy statement
  ispec = np
#endif

  ! adds corner points
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
      val = data(1,1,1,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      ! .mesh
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob1) = .true.
    endif

    if (.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      x = xstore(iglob2)
      y = ystore(iglob2)
      z = zstore(iglob2)
      val = data(NGLLX,1,1,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob2) = .true.
    endif

    if (.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      x = xstore(iglob3)
      y = ystore(iglob3)
      z = zstore(iglob3)
      val = data(NGLLX,NGLLY,1,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob3) = .true.
    endif

    if (.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      x = xstore(iglob4)
      y = ystore(iglob4)
      z = zstore(iglob4)
      val = data(1,NGLLY,1,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob4) = .true.
    endif

    if (.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      x = xstore(iglob5)
      y = ystore(iglob5)
      z = zstore(iglob5)
      val = data(1,1,NGLLZ,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob5) = .true.
    endif

    if (.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      x = xstore(iglob6)
      y = ystore(iglob6)
      z = zstore(iglob6)
      val = data(NGLLX,1,NGLLZ,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob6) = .true.
    endif

    if (.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      x = xstore(iglob7)
      y = ystore(iglob7)
      z = zstore(iglob7)
      val = data(NGLLX,NGLLY,NGLLZ,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob7) = .true.
    endif

    if (.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      x = xstore(iglob8)
      y = ystore(iglob8)
      z = zstore(iglob8)
      val = data(1,NGLLY,NGLLZ,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
      ! VTK/VTU
      total_dat_xyz(1,np+numpoin) = x
      total_dat_xyz(2,np+numpoin) = y
      total_dat_xyz(3,np+numpoin) = z
      !write(IOUT_VTK,'(3e18.6)') x,y,z
      total_dat(np+numpoin) = val
#else
      call write_real(x)
      call write_real(y)
      call write_real(z)
      call write_real(val)
#endif
      mask_ibool(iglob8) = .true.
    endif
  enddo ! ispec

  end subroutine cvd_write_corners

!=============================================================

  subroutine cvd_write_GLL_points(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,data,numpoin,np)

! writes out locations of all GLL points of spectral elements
  use constants
  use combine_vtk_par
  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: xstore, ystore, zstore
  real(kind=CUSTOM_REAL),dimension(NGLLY,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: data

  integer,intent(inout) :: numpoin,np

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  real :: x, y, z, val
  integer :: ispec,i,j,k,iglob,ier

  ! writes out point locations and values
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1151')
  if (ier /= 0) stop 'error allocating array mask_ibool'
  mask_ibool(:) = .false.

  ! dummy statement to avoid compiler warning
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
  ! nothing to do
#else
  ! dummy statement
  ispec = np
#endif

  ! adds GLL points
  numpoin = 0
  do ispec = 1,NSPEC_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            x = xstore(iglob)
            y = ystore(iglob)
            z = zstore(iglob)
            val = data(i,j,k,ispec)
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
            ! VTK/VTU
            total_dat_xyz(1,np+numpoin) = x
            total_dat_xyz(2,np+numpoin) = y
            total_dat_xyz(3,np+numpoin) = z
            !write(IOUT_VTK,'(3e18.6)') x,y,z
            total_dat(np+numpoin) = val
#else
            call write_real(x)
            call write_real(y)
            call write_real(z)
            call write_real(val)
#endif
            mask_ibool(iglob) = .true.
          endif
        enddo ! i
      enddo ! j
    enddo ! k
  enddo !ispec

  end subroutine cvd_write_GLL_points

!=============================================================

! writes out locations of spectral element corners only

  subroutine cvd_write_corner_elements(NSPEC_AB,NGLOB_AB,ibool,np,ne,numpoin,nelement)

  use constants
  use combine_vtk_par
  implicit none

  integer, intent(in) :: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB), intent(in) :: ibool
  integer, intent(inout) :: np,ne
  integer, intent(out) :: numpoin,nelement

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  integer,dimension(:),allocatable :: num_ibool
  integer :: ispec,ier
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  integer :: n1, n2, n3, n4, n5, n6, n7, n8

!  ! outputs total number of elements for all slices
!  if (it == 1) then
!#ifdef USE_VTK_INSTEAD_OF_MESH
!    ! VTK
!    ! note: indices for vtk start at 0
!    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nee,nee*9
!#else
!    call write_integer(nee)
!#endif
!  endif

  ! writes out element indices
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1152')
  allocate(num_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1153')
  if (ier /= 0) stop 'error allocating array mask_ibool'
  mask_ibool(:) = .false.
  num_ibool(:) = 0

  numpoin = 0
  nelement = 0
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

    ! outputs corner indices (starting with 0)
    n1 = num_ibool(iglob1) -1 + np
    n2 = num_ibool(iglob2) -1 + np
    n3 = num_ibool(iglob3) -1 + np
    n4 = num_ibool(iglob4) -1 + np
    n5 = num_ibool(iglob5) -1 + np
    n6 = num_ibool(iglob6) -1 + np
    n7 = num_ibool(iglob7) -1 + np
    n8 = num_ibool(iglob8) -1 + np

    ! counts elements
    nelement = nelement + 1

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
    ! VTK/VTU
    ! note: indices for VTK start at 0
    total_dat_con(1,nelement + ne) = n1
    total_dat_con(2,nelement + ne) = n2
    total_dat_con(3,nelement + ne) = n3
    total_dat_con(4,nelement + ne) = n4
    total_dat_con(5,nelement + ne) = n5
    total_dat_con(6,nelement + ne) = n6
    total_dat_con(7,nelement + ne) = n7
    total_dat_con(8,nelement + ne) = n8
    !write(IOUT_VTK,'(9i12)') 8,n1,n2,n3,n4,n5,n6,n7,n8
#else
    ! .mesh
    call write_integer(n1)
    call write_integer(n2)
    call write_integer(n3)
    call write_integer(n4)
    call write_integer(n5)
    call write_integer(n6)
    call write_integer(n7)
    call write_integer(n8)
#endif

  enddo

  ! elements written
  !nelement = NSPEC_AB
  ne = ne + nelement

  ! updates points written
  np = np + numpoin

  end subroutine cvd_write_corner_elements


!=============================================================


  subroutine cvd_write_GLL_elements(NSPEC_AB,NGLOB_AB,ibool,np,ne,numpoin,nelement)

! writes out indices of elements given by GLL points
  use constants
  use combine_vtk_par
  implicit none

  integer, intent(in):: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB), intent(in) :: ibool
  integer, intent(inout) :: np,ne
  integer, intent(out) :: numpoin,nelement

  ! local parameters
  logical,dimension(:),allocatable :: mask_ibool
  integer,dimension(:),allocatable :: num_ibool
  integer :: ispec,i,j,k,ier
  integer :: iglob,iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  integer :: n1, n2, n3, n4, n5, n6, n7, n8

  ! outputs total number of elements for all slices
!  if (it == 1) then
!#ifdef USE_VTK_INSTEAD_OF_MESH
!    ! VTK
!    ! note: indices for vtk start at 0
!    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nee,nee*9
!#else
!    !nee = nelement * num_node
!    call write_integer(nee)
!#endif
!  endif

  ! sets numbering num_ibool respecting mask
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1154')
  allocate(num_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1155')
  if (ier /= 0) stop 'error allocating array mask_ibool'
  mask_ibool(:) = .false.
  num_ibool(:) = 0

  ! sets mask
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

  ! outputs GLL subelement
  nelement = 0
  do ispec = 1, NSPEC_AB
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

          ! counts elements
          nelement = nelement + 1

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
          ! VTK/VTU
          ! note: indices for VTK start at 0
          total_dat_con(1,nelement + ne) = n1
          total_dat_con(2,nelement + ne) = n2
          total_dat_con(3,nelement + ne) = n3
          total_dat_con(4,nelement + ne) = n4
          total_dat_con(5,nelement + ne) = n5
          total_dat_con(6,nelement + ne) = n6
          total_dat_con(7,nelement + ne) = n7
          total_dat_con(8,nelement + ne) = n8
          !write(IOUT_VTK,'(9i12)') 8,n1,n2,n3,n4,n5,n6,n7,n8
#else
          ! .mesh
          call write_integer(n1)
          call write_integer(n2)
          call write_integer(n3)
          call write_integer(n4)
          call write_integer(n5)
          call write_integer(n6)
          call write_integer(n7)
          call write_integer(n8)
#endif
        enddo
      enddo
    enddo
  enddo

  ! elements written
  !nelement = NSPEC_AB * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
  ne = ne + nelement

  ! updates points written
  np = np + numpoin

  end subroutine cvd_write_GLL_elements

