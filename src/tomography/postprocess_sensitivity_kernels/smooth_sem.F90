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

! XSMOOTH_SEM
!
! USAGE
!   mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR [USE_GPU]
!
!
! COMMAND LINE ARGUMENTS
!   SIGMA_H                - horizontal smoothing radius
!   SIGMA_V                - vertical smoothing radius
!   KERNEL_NAME            - kernel name, e.g. alpha_kernel
!   INPUT_DIR              - directory from which kernels are read
!   OUTPUT_DIR             - directory to which smoothed kernels are written
!   USE_GPU                - (optional) use GPUs for computation
!
! DESCRIPTION
!   Smooths kernels by convolution with a Gaussian. Writes the resulting
!   smoothed kernels to OUTPUT_DIR.
!
!   Files written to OUTPUT_DIR have the suffix 'smooth' appended,
!   e.g. proc***alpha_kernel.bin becomes proc***alpha_kernel_smooth.bin
!
!   This program's primary use case is to smooth kernels. It can be used though on
!   any scalar field of dimension (NGLLX,NGLLY,NGLLZ,NSPEC).
!
!   This is a parrallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarassingly-parallel
!   fashion.

#include "config.fh"

program smooth_sem

  use constants, only: USE_QUADRATURE_RULE_FOR_SMOOTHING

  use postprocess_par, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
    MAX_STRING_LEN,IIN,IOUT,GAUSSALPHA,GAUSSBETA,PI,MAX_KERNEL_NAMES

  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic,min_resolved_period
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
    phistore,tortstore,rhoarraystore
  use specfem_par_movie

  implicit none

  integer, parameter :: NARGS = 6

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat,dat_smooth
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dummy ! for jacobian read
  integer :: NSPEC_N, NGLOB_N, NSPEC_IRREGULAR_N

  integer :: i,j,k,iglob,ier,ispec2,ispec,ispec_irreg,inum
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif
  integer :: icounter,num_slices
  integer :: iproc,ncuda_devices

  integer(kind=8) :: Container

  integer,parameter :: MAX_NODE_LIST = 300
  integer :: node_list(MAX_NODE_LIST)
  logical :: do_include_slice

  character(len=MAX_STRING_LEN) :: arg(6)
  character(len=MAX_STRING_LEN) :: kernel_name, input_dir, output_dir
  character(len=MAX_STRING_LEN) :: prname_lp
  character(len=MAX_STRING_LEN*2) :: local_data_file


  character(len=MAX_STRING_LEN),dimension(:),allocatable :: kernel_names
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  integer :: nker

  ! smoothing parameters
  character(len=MAX_STRING_LEN*2) :: ks_file

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3
  real(kind=CUSTOM_REAL) :: sigma_h2_inv,sigma_v2_inv
  real(kind=CUSTOM_REAL) :: sigma_h3_sq,sigma_v3_sq

  real(kind=CUSTOM_REAL) :: x0, y0, z0, norm, norm_h, norm_v
  real(kind=CUSTOM_REAL) :: center_x0, center_y0, center_z0
  real(kind=CUSTOM_REAL) :: center_x, center_y, center_z

  real(kind=CUSTOM_REAL) :: max_old,max_new,max_old_all,max_new_all
  real(kind=CUSTOM_REAL) :: min_old,min_new,min_old_all,min_new_all

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: exp_val,factor
  real(kind=CUSTOM_REAL) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: tk, bk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: xl, yl, zl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: xx, yy, zz

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx0, cy0, cz0
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx, cy, cz

  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: element_size

  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  ! reference slice
  real(kind=CUSTOM_REAL) :: x_min_ref,x_max_ref
  real(kind=CUSTOM_REAL) :: y_min_ref,y_max_ref
  real(kind=CUSTOM_REAL) :: z_min_ref,z_max_ref
  real(kind=CUSTOM_REAL) :: x_min,x_max
  real(kind=CUSTOM_REAL) :: y_min,y_max
  real(kind=CUSTOM_REAL) :: z_min,z_max
  real(kind=CUSTOM_REAL) :: dim_x,dim_y,dim_z

  logical :: BROADCAST_AFTER_READ, USE_GPU

  ! timing
  double precision, external :: wtime
  real :: t1,t2

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running XSMOOTH_SEM"
  call synchronize_all()

  ! timing
  t1 = wtime()

  ! parse command line arguments
  if (command_argument_count() < NARGS-1) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR [GPU_MODE]'
      print *,'  with'
      print *,'   SIGMA_H SIGMA_V  - horizontal & vertical smoothing lengths'
      print *,'   KERNEL_NAME      - sensitivity kernel name (e.g., alpha_kernel for proc***_alpha_kernel.bin files)'
      print *,'   INPUT_DIR        - input directory holding kernel files'
      print *,'   OUPUT_DIR        - output directory for smoothed kernel'
      print *,'   GPU_MODE         - (optional) set to .true. to use GPU, otherwise set to .false. for CPU run (default off)'
      print *
      stop 'Please check command line arguments'
    endif
  endif
  call synchronize_all()

  ! allocates array
  allocate(kernel_names(MAX_KERNEL_NAMES),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel_names array'
  kernel_names(:) = ''

  ! parse command line arguments
  do i = 1, NARGS
    if (command_argument_count() >= i) then
      call get_command_argument(i,arg(i), status=ier)
    else
      ! optional argument
      arg(i) = ''
    endif
  enddo

  read(arg(1),*) sigma_h
  read(arg(2),*) sigma_v
  kernel_names_comma_delimited = arg(3)
  input_dir= arg(4)
  output_dir = arg(5)
  if (command_argument_count() == 6) then
    read(arg(6),*) USE_GPU
  else
    USE_GPU = .false.
  endif

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  kernel_name = kernel_names(1)

  if (USE_GPU) call initialize_gpu_device(myrank,ncuda_devices)

  if (nker > 1) then
    if (myrank == 0) then
      ! The machinery for reading multiple names from the command line is in place,
      ! but the smoothing routines themselves have not yet been modified to work on multiple arrays.
      if (myrank == 0) then
        print *,'Smoothing only first name in list: ',trim(kernel_name)
        print *
      endif
    endif
  endif
  call synchronize_all()

  ! check smoothing radii
  sigma_h2 = 2.0 * sigma_h ** 2  ! factor two for Gaussian distribution with standard variance sigma
  sigma_v2 = 2.0 * sigma_v ** 2

  if (sigma_h2 < 1.e-18) stop 'Error sigma_h2 zero, must non-zero'
  if (sigma_v2 < 1.e-18) stop 'Error sigma_v2 zero, must non-zero'

  ! adds margin to search radius
  element_size = max(sigma_h,sigma_v) * 0.5

  ! search radius
  sigma_h3 = 3.0  * sigma_h + element_size
  sigma_v3 = 3.0  * sigma_v + element_size

  ! helper variables
  sigma_h2_inv = 1.0_CUSTOM_REAL / sigma_h2
  sigma_v2_inv = 1.0_CUSTOM_REAL / sigma_v2

  sigma_h3_sq = sigma_h3 * sigma_h3
  sigma_v3_sq = sigma_v3 * sigma_v3

  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
  ! note: smoothing is using a Gaussian (ellipsoid for sigma_h /= sigma_v),
  norm_h = real(2.0*PI*sigma_h**2,kind=CUSTOM_REAL)
  norm_v = real(sqrt(2.0*PI) * sigma_v,kind=CUSTOM_REAL)
  norm   = norm_h * norm_v

  ! user output
  if (myrank == 0) then
    print *,'command line arguments:'
    print *,'  smoothing sigma_h , sigma_v                : ',sigma_h,sigma_v
    ! scale length: approximately S ~ sigma * sqrt(8.0) for a Gaussian smoothing
    print *,'  smoothing scalelengths horizontal, vertical: ',sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *,'  input dir : ',trim(input_dir)
    print *,'  output dir: ',trim(output_dir)
    print *,"  GPU_MODE: ", USE_GPU
    print *
  endif

  ! reads the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(BROADCAST_AFTER_READ)

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED not supported yet for smoothing, please rerun program...'

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *,'Error number of processors supposed to run on: ',NPROC
      print *,'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *,'Please rerun with: mpirun -np ',NPROC,' bin/xsmooth_sem .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  call read_mesh_for_init()

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),irregular_element_number(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 980')

  if (NSPEC_IRREGULAR > 0) then
    ! allocate arrays for storing the databases
    allocate(xixstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 981')
    allocate(xiystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 982')
    allocate(xizstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 983')
    allocate(etaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 984')
    allocate(etaystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 985')
    allocate(etazstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 986')
    allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 987')
    allocate(gammaystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 988')
    allocate(gammazstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 989')
    allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 990')
   else
    ! allocate arrays for storing the databases
    allocate(xixstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 991')
    allocate(xiystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 992')
    allocate(xizstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 993')
    allocate(etaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 994')
    allocate(etaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 995')
    allocate(etazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 996')
    allocate(gammaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 997')
    allocate(gammaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 998')
    allocate(gammazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 999')
    allocate(jacobianstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1000')
  endif
  if (ier /= 0) stop 'Error allocating arrays for databases'

  ! mesh node locations
  allocate(xstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1001')
  allocate(ystore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1002')
  allocate(zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1003')
  if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1004')
  allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1005')
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating rho array 1005')
  if (ier /= 0) stop 'Error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1006')
  allocate(ispec_is_elastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1007')
  allocate(ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1008')
  if (ier /= 0) stop 'Error allocating arrays for material flags'
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! reads in external mesh
  call read_mesh_databases()

  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  ! outputs infos
  if (myrank == 0) then
    print *,'mesh dimensions:'
    print *,'  Xmin and Xmax of the model = ',x_min_glob,x_max_glob
    print *,'  Ymin and Ymax of the model = ',y_min_glob,y_max_glob
    print *,'  Zmin and Zmax of the model = ',z_min_glob,z_max_glob
    print *
    print *,'  Max GLL point distance = ',distance_max_glob
    print *,'  Min GLL point distance = ',distance_min_glob
    print *,'  Max/min ratio = ',distance_max_glob/distance_min_glob
    print *
    print *,'  Max element size = ',elemsize_max_glob
    print *,'  Min element size = ',elemsize_min_glob
    print *,'  Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
    print *
  endif

  ! mesh resolution
  call check_mesh_resolution(NSPEC_AB,NGLOB_AB, &
                             ibool,xstore,ystore,zstore, &
                             ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                             kappastore,mustore,rhostore, &
                             phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                             DT,model_speed_max,min_resolved_period)

  ! for smoothing, we use cell centers to find and locate nearby elements
  !
  ! sets the location of the center of the elements and local points
  allocate(xl(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1013')
  allocate(yl(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1014')
  allocate(zl(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1015')
  allocate(cx0(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1016')
  allocate(cy0(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1017')
  allocate(cz0(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1018')
  if (ier /= 0) stop 'Error allocating array xl etc.'

  ! sets element center location
  do ispec = 1, nspec_AB

    DO_LOOP_IJK
      iglob = ibool(INDEX_IJK,ispec)
      xl(INDEX_IJK,ispec) = xstore(iglob)
      yl(INDEX_IJK,ispec) = ystore(iglob)
      zl(INDEX_IJK,ispec) = zstore(iglob)
    ENDDO_LOOP_IJK

    cx0(ispec) = (xl(1,1,1,ispec) + xl(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cy0(ispec) = (yl(1,1,1,ispec) + yl(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cz0(ispec) = (zl(1,1,1,ispec) + zl(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
  enddo

  ! reference slice dimension
  x_min_ref = minval(xstore)
  x_max_ref = maxval(xstore)
  y_min_ref = minval(ystore)
  y_max_ref = maxval(ystore)
  z_min_ref = minval(zstore)
  z_max_ref = maxval(zstore)

  dim_x = abs(x_max_ref - x_min_ref)
  dim_y = abs(y_max_ref - y_min_ref)
  dim_z = abs(z_max_ref - z_min_ref)

  ! user output
  if (myrank == 0) then
    print *,'smoothing:'
    print *,'  single slice dimensions in x/y/z-direction: ',dim_x,' / ',dim_y,' / ',dim_z
    print *,'  Gaussian search radius horizontal =',sigma_h3,' vertical =',sigma_v3
    print *
  endif

  ! checks if Gaussian support exceeds slice dimensions
  if (sigma_h3 >= dim_x .or. sigma_h3 >= dim_y .or. sigma_v3 >= dim_z) then
    ! Gaussian support is likely larger than the direct neighbor and has support in a much wider area
    ! user output
    if (myrank == 0) then
      print *,'  using large Gaussian with respect to slice dimension for smoothing'
      print *
    endif
  endif

  ! frees memory
  deallocate(xstore,ystore,zstore)
  deallocate(xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore)
  deallocate(ibool,irregular_element_number)
  deallocate(jacobianstore)

  call synchronize_all()

  ! initializes node list
  node_list(:) = -1

  ! adds this partition itself first
  icounter = 1
  node_list(icounter) = myrank

  ! sets up slices to process
  do iproc = 0,NPROC-1
    ! skip own process slice, has already been added
    if (iproc == myrank) cycle

    ! checks if slice is a direct neighbor
    do_include_slice = .false.
    do i = 1,num_interfaces_ext_mesh
      if (iproc == my_neighbors_ext_mesh(i)) then
        ! found a neighbor slice
        do_include_slice = .true.
        exit
      endif
    enddo

    if (.not. do_include_slice) then
      ! note: Gaussian support might be larger than closest neighbor slices
      !       we add all slices close enough to still have an influence

      ! checks distances to this slice
      ! reads in slice mesh
      ! neighbor database file
      call create_name_database(prname,iproc,LOCAL_PATH)
      prname_lp = prname(1:len_trim(prname))//'external_mesh.bin'

      ! gets number of elements and global points for this partition
      open(unit=IIN,file=trim(prname_lp),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error could not open database file: ',trim(prname_lp)
        call exit_mpi(myrank, 'Error reading neighbors external mesh file')
      endif
      read(IIN) NSPEC_N
      read(IIN) NGLOB_N
      read(IIN) NSPEC_IRREGULAR_N

      ! allocates mesh arrays
      allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1019')
      if (ier /= 0) stop 'Error allocating array ibool'
      allocate(xstore(NGLOB_N),ystore(NGLOB_N),zstore(NGLOB_N),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1020')
      if (ier /= 0) stop 'Error allocating array xstore etc.'

      ! ibool file
      read(IIN) ibool

      ! global point arrays
      read(IIN) xstore
      read(IIN) ystore
      read(IIN) zstore
      close(IIN)

      ! determines min/max values of slice mesh
      x_min = minval(xstore)
      x_max = maxval(xstore)
      y_min = minval(ystore)
      y_max = maxval(ystore)
      z_min = minval(zstore)
      z_max = maxval(zstore)

      ! slice dimensions
      dim_x = abs(x_max - x_min)
      dim_y = abs(y_max - y_min)
      dim_z = abs(z_max - z_min)

      ! frees memory
      deallocate(ibool)
      deallocate(xstore,ystore,zstore)

      ! re-evalutes distances between slices
      ! checks if edges are within search radius from reference slice
      ! note: for slices with irregular shapes (e.g. from scotch decomposition), this assumes a box slice shape
      !       and might add slices even if further separated; a conservative approach
      if ((x_min < x_max_ref + sigma_h3 .and. x_max > x_min_ref - sigma_h3) &
          .and. (y_min < y_max_ref + sigma_h3 .and. y_max > y_min_ref - sigma_h3) &
          .and. (z_min < z_max_ref + sigma_v3 .and. z_max > z_min_ref - sigma_v3)) then
        do_include_slice = .true.
      endif

    endif

    ! adds to smoothing list neighbors
    if (do_include_slice) then
      icounter = icounter + 1

      ! checks bounds
      if (icounter > MAX_NODE_LIST) stop 'Error number of interfaces exceeds MAX_NODE_LIST'

      ! adds slice to smoothing list
      node_list(icounter) = iproc
    endif
  enddo
  ! set total number of slices to loop over
  num_slices = icounter

  ! synchronizes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *,'slices:',num_slices
    print *,'  rank ',myrank,'  has smoothing slices:'
    print *,node_list(1:num_slices)
    print *
  endif

  !do i=0,sizeprocs-1
  !  if (myrank == i) then
  !    print *,'rank:',myrank,'  smoothing slices'
  !    print *,node_list(1:num_slices)
  !    print *
  !  endif
  !enddo

  ! for jacobian and weights
  ! GLL points weights
  if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
        enddo
      enddo
    enddo
  endif

  ! loops over slices
  ! each process reads in his own neighbor slices and Gaussian filters the values
  allocate(tk(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1021')
  allocate(bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1022')
  if (ier /= 0) stop 'Error allocating array tk and bk'

  tk = 0.0_CUSTOM_REAL
  bk = 0.0_CUSTOM_REAL

  ! GPU setup
  if (USE_GPU) then
    call prepare_GPU_smooth(Container,xl,yl,zl,sigma_h2_inv,sigma_v2_inv,sigma_h3_sq,sigma_v3_sq,NSPEC_AB,nker,wgll_cube)

    ! synchronizes all processes
    call synchronize_all()
  endif


  do inum = 1,num_slices

    iproc = node_list(inum)

    if (myrank == 0) print *,'  reading slice:',iproc

    ! neighbor database file
    call create_name_database(prname,iproc,LOCAL_PATH)
    prname_lp = prname(1:len_trim(prname))//'external_mesh.bin'

    ! gets number of point locations
    open(unit=IIN,file=trim(prname_lp),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname_lp)
      call exit_mpi(myrank, 'Error reading neighbors external mesh file')
    endif
    read(IIN) NSPEC_N
    read(IIN) NGLOB_N
    read(IIN) NSPEC_IRREGULAR_N

    ! allocates arrays
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1023')
    if (ier /= 0) stop 'Error allocating array ibool'
    allocate(xstore(NGLOB_N),ystore(NGLOB_N),zstore(NGLOB_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1024')
    if (ier /= 0) stop 'Error allocating array xstore etc.'

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      if (NSPEC_IRREGULAR_N > 0) then
        allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR_N),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1025')
        if (ier /= 0) stop 'Error allocating array jacobian'
        allocate(dummy(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1026')
        if (ier /= 0) stop 'Error allocating array dummy'
      else
        allocate(jacobianstore(1,1,1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1027')
        if (ier /= 0) stop 'Error allocating array jacobian'
        allocate(dummy(1,1,1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1028')
        if (ier /= 0) stop 'Error allocating array dummy'
      endif
      allocate(irregular_element_number(NSPEC_N),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1029')
      if (ier /= 0) stop 'Error allocating array irregular_element_number'
    endif

    ! ibool file
    read(IIN) ibool

    ! global point arrays
    read(IIN) xstore
    read(IIN) ystore
    read(IIN) zstore

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      ! reads in jacobian
      read(IIN) irregular_element_number
      read(IIN) xix_regular
      read(IIN) jacobian_regular
      read(IIN) dummy ! xix
      read(IIN) dummy ! xiy
      read(IIN) dummy ! xiz
      read(IIN) dummy ! etax
      read(IIN) dummy ! etay
      read(IIN) dummy ! etaz
      read(IIN) dummy ! gammax
      read(IIN) dummy ! gammay
      read(IIN) dummy ! gammaz
      read(IIN) jacobianstore
    endif

    close(IIN)

    ! get the location of the center of the elements and local points
    allocate(xx(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1030')
    allocate(yy(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1031')
    allocate(zz(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1032')
    allocate(cx(NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1033')
    allocate(cy(NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1034')
    allocate(cz(NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1035')
    if (ier /= 0) stop 'Error allocating array xx etc.'

    ! sets element center location
    do ispec = 1, NSPEC_N

      DO_LOOP_IJK

        iglob = ibool(INDEX_IJK,ispec)
        xx(INDEX_IJK,ispec) = xstore(iglob)
        yy(INDEX_IJK,ispec) = ystore(iglob)
        zz(INDEX_IJK,ispec) = zstore(iglob)

      ENDDO_LOOP_IJK

      ! calculate element center location
      cx(ispec) = (xx(1,1,1,ispec) + xx(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
      cy(ispec) = (yy(1,1,1,ispec) + yy(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
      cz(ispec) = (zz(1,1,1,ispec) + zz(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    enddo

    deallocate(xstore,ystore,zstore)
    deallocate(ibool)

    ! data file
    write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_'
    local_data_file = trim(prname) // trim(kernel_name) // '.bin'

    open(unit = IIN,file = trim(local_data_file),status='old',action='read',form ='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening data file: ',trim(local_data_file)
      stop 'Error opening data file'
    endif

    allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1036')
    if (ier /= 0) stop 'Error allocating dat array'

    read(IIN) dat
    close(IIN)

    ! statistics
    if (iproc == myrank) then
      min_old = minval(dat(:,:,:,:))
      max_old = maxval(dat(:,:,:,:))
    endif

    ! finds closest elements for smoothing
    !if (myrank==0) print *, '  start looping over elements and points for smoothing ...'
    if (USE_GPU) then
      call compute_smooth(Container,jacobianstore,jacobian_regular,irregular_element_number, &
                          xx,yy,zz,dat,NSPEC_N,NSPEC_IRREGULAR_N)
    else
      ! loop over elements to be smoothed in the current slice
      do ispec = 1, NSPEC_AB

        ! element center position
        center_x0 = cx0(ispec)
        center_y0 = cy0(ispec)
        center_z0 = cz0(ispec)

        ! --- only double loop over the elements in the search radius ---
        do ispec2 = 1, NSPEC_N

          ! search element center position
          center_x = cx(ispec2)
          center_y = cy(ispec2)
          center_z = cz(ispec2)

          ! calculates horizontal and vertical distance between two element centers
          ! (squared distances)
          call get_distance_vec(dist_h,dist_v,center_x0,center_y0,center_z0,center_x,center_y,center_z)

          ! checks distance between centers of elements
          if (dist_h > sigma_h3_sq .or. dist_v > sigma_v3_sq) cycle

          ! integration factors
          if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
            ispec_irreg = irregular_element_number(ispec2)
            if (ispec_irreg == 0) jacobianl = jacobian_regular
            do k=1,NGLLZ
              do j=1,NGLLY
                do i=1,NGLLX
                  if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)
                  factor(i,j,k) = jacobianl * wgll_cube(i,j,k)
                enddo
              enddo
            enddo
          endif

          ! loop over GLL points of the elements in current slice (ispec)
          DO_LOOP_IJK

            ! reference location
            ! current point (i,j,k,ispec) location, Cartesian coordinates
            x0 = xl(INDEX_IJK,ispec)
            y0 = yl(INDEX_IJK,ispec)
            z0 = zl(INDEX_IJK,ispec)

            ! calculate weights based on Gaussian smoothing
            call smoothing_weights_vec(x0,y0,z0,sigma_h2_inv,sigma_v2_inv,exp_val, &
                                       xx(:,:,:,ispec2),yy(:,:,:,ispec2),zz(:,:,:,ispec2))

            ! adds GLL integration weights
            if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
              exp_val(:,:,:) = exp_val(:,:,:) * factor(:,:,:)
            endif

            ! adds contribution of element ispec2 to smoothed kernel values
            tk(INDEX_IJK,ispec) = tk(INDEX_IJK,ispec) + sum(exp_val(:,:,:) * dat(:,:,:,ispec2))

            ! normalization, integrated values of Gaussian smoothing function
            bk(INDEX_IJK,ispec) = bk(INDEX_IJK,ispec) + sum(exp_val(:,:,:))

          ENDDO_LOOP_IJK

        enddo ! ispec2
      enddo ! ispec

    endif ! GPU_MODE

    ! frees arrays
    deallocate(dat)
    deallocate(xx,yy,zz)
    deallocate(cx,cy,cz)

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      deallocate(irregular_element_number)
      deallocate(jacobianstore)
      deallocate(dummy)
    endif

  enddo ! iproc

  ! normalizes/scaling factor
  if (myrank == 0) then
    print *
    print *,'Scaling values: min/max = ',minval(bk),maxval(bk)
  endif

  allocate(dat_smooth(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1037')
  if (ier /= 0) stop 'Error allocating array dat_smooth'

  dat_smooth(:,:,:,:) = 0.0_CUSTOM_REAL

  if (USE_GPU) then
    call get_smooth(Container,dat_smooth)
  else
    do ispec = 1, NSPEC_AB

      DO_LOOP_IJK

        ! checks the normalization criterion
        !if (abs(bk(i,j,k,ispec) - norm) > 1.e-4) then
        !  print *, 'Problem norm here --- ', ispec, i, j, k, bk(i,j,k,ispec), norm
        !endif
        if (abs(bk(INDEX_IJK,ispec)) < 1.e-18) then
          print *, 'Problem norm here --- ', ispec, i, j, k, bk(i,j,k,ispec), norm
        endif

        ! normalizes smoothed kernel values by integral value of Gaussian weighting
        dat_smooth(INDEX_IJK,ispec) = tk(INDEX_IJK,ispec) / bk(INDEX_IJK,ispec)

      ENDDO_LOOP_IJK

    enddo !  ispec

  endif ! GPU_MODE

  ! frees memory
  deallocate(tk,bk)

  ! statistics
  min_new = minval(dat_smooth(:,:,:,:))
  max_new = maxval(dat_smooth(:,:,:,:))

  ! file output
  ! smoothed kernel file name
  write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name)//'_smooth.bin'

  open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening smoothed kernel file'
  write(IOUT) dat_smooth(:,:,:,:)
  close(IOUT)
  if (myrank == 0) print *,'written: ',trim(ks_file)

  ! frees memory
  deallocate(dat_smooth)

  ! synchronizes
  call synchronize_all()

  ! the min/maximum value for the smoothed kernel
  call min_all_cr(min_old, min_old_all)
  call min_all_cr(min_new, min_new_all)
  call max_all_cr(max_old, max_old_all)
  call max_all_cr(max_new,max_new_all)

  if (myrank == 0) then
    print *
    print *,'Minimum data value before smoothing = ', min_old_all
    print *,'Minimum data value after smoothing  = ', min_new_all
    print *
    print *,'Maximum data value before smoothing = ', max_old_all
    print *,'Maximum data value after smoothing  = ', max_new_all
    print *
    close(IMAIN)
  endif

  ! timing
  t2 = wtime()

  if (myrank == 0) then
    if (USE_GPU) then
      print *,'Computation time with GPU:',t2-t1
    else
      print *,'Computation time with CPU:',t2-t1
    endif
  endif

  ! stop all the processes and exit
  call finalize_mpi()

  contains
!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,y0,z0,sigma_h2_inv,sigma_v2_inv,exp_val, &
                                   xx_elem,yy_elem,zz_elem)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: xx_elem, yy_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,sigma_h2_inv,sigma_v2_inv

  ! local parameters
  real(kind=CUSTOM_REAL) :: dist_h_sq,dist_v_sq
  real(kind=CUSTOM_REAL) :: val
  real(kind=CUSTOM_REAL) :: x1,y1,z1
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  DO_LOOP_IJK

    ! point in second slice
    x1 = xx_elem(INDEX_IJK)
    y1 = yy_elem(INDEX_IJK)
    z1 = zz_elem(INDEX_IJK)

    ! gets vertical and horizontal distance
    ! vertical distance (squared)
    dist_v_sq = (z0-z1)*(z0-z1)

    ! horizontal distance (squared)
    dist_h_sq = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)

    ! Gaussian function
    val = exp(- dist_h_sq * sigma_h2_inv - dist_v_sq * sigma_v2_inv)

    ! stores values in element array
    exp_val(INDEX_IJK) = val

  ENDDO_LOOP_IJK

  end subroutine smoothing_weights_vec

!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_vec(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

! returns vector lengths as distances in radial and horizontal direction
! only for flat Earth with z in vertical direction

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,x1,y1,z1

  ! vertical distance
  !dist_v = sqrt( (z0-z1)*(z0-z1) )
  ! squared (to avoid costly square-root
  dist_v = (z0-z1)*(z0-z1)

  ! horizontal distance
  !dist_h = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) )
  ! squared
  dist_h = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)

  end subroutine get_distance_vec

end program smooth_sem
