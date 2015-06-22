!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

! XSMOOTH_SEM
!
! USAGE
!   mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR
!
!
! COMMAND LINE ARGUMENTS
!   SIGMA_H                - horizontal smoothing radius
!   SIGMA_V                - vertical smoothing radius
!   KERNEL_NAME            - kernel name, e.g. alpha_kernel
!   INPUT_DIR              - directory from which kernels are read
!   OUTPUT_DIR             - directory to which smoothed kernels are written
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


program smooth_sem

  use postprocess_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,NGLLSQUARE, &
    MAX_STRING_LEN,IIN,IOUT, &
    GAUSSALPHA,GAUSSBETA,PI,TWO_PI, &
    MAX_KERNEL_NAMES

  use specfem_par
  use specfem_par_elastic,only: ELASTIC_SIMULATION,ispec_is_elastic,rho_vp,rho_vs,min_resolved_period
  use specfem_par_acoustic,only: ACOUSTIC_SIMULATION,ispec_is_acoustic
  use specfem_par_poroelastic,only: POROELASTIC_SIMULATION,ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
    phistore,tortstore,rhoarraystore
  use specfem_par_movie

  implicit none

  integer, parameter :: NARGS = 5

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat,dat_smooth
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dummy
  integer :: NSPEC_N, NGLOB_N

  integer :: i,j,k,iglob,ier,ispec2,ispec,inum
  integer :: iproc

  integer,parameter :: MAX_NODE_LIST = 300
  integer :: node_list(MAX_NODE_LIST)

  character(len=MAX_STRING_LEN) :: arg(5)
  character(len=MAX_STRING_LEN) :: kernel_name, input_dir, output_dir
  character(len=MAX_STRING_LEN) :: prname_lp
  character(len=MAX_STRING_LEN*2) :: local_data_file


  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  integer :: nker

  ! smoothing parameters
  character(len=MAX_STRING_LEN*2) :: ks_file

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3
  real(kind=CUSTOM_REAL) :: x0, y0, z0, norm, norm_h, norm_v, max_old, max_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: exp_val !,factor

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

  logical :: BROADCAST_AFTER_READ

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running XSMOOTH_SEM"
  call synchronize_all()

  ! parse command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
        print *, 'USAGE:  mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR'
      stop ' Please check command line arguments'
    endif
  endif
  call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo

  read(arg(1),*) sigma_h
  read(arg(2),*) sigma_v
  kernel_names_comma_delimited = arg(3)
  input_dir= arg(4)
  output_dir = arg(5)

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  kernel_name = trim(kernel_names(1))

  if (nker > 1) then
    if (myrank == 0) then
      ! The machinery for reading multiple names from the command line is in place,
      ! but the smoothing routines themselves have not yet been modified to work
      !  on multiple arrays.
      if (myrank == 0) print *, 'Smoothing only first name in list: ', trim(kernel_name)
      if (myrank == 0) print *
    endif
  endif
  call synchronize_all()

  ! check smoothing radii
  sigma_h2 = 2.0 * sigma_h ** 2  ! factor two for gaussian distribution with standard variance sigma
  sigma_v2 = 2.0 * sigma_v ** 2

  if (sigma_h2 < 1.e-18) stop 'Error sigma_h2 zero, must non-zero'
  if (sigma_v2 < 1.e-18) stop 'Error sigma_v2 zero, must non-zero'

  ! adds margin to search radius
  element_size = max(sigma_h,sigma_v) * 0.5

  ! search radius
  sigma_h3 = 3.0  * sigma_h + element_size
  sigma_v3 = 3.0  * sigma_v + element_size

  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
  ! note: smoothing is using a gaussian (ellipsoid for sigma_h /= sigma_v),
  norm_h = 2.0*PI*sigma_h**2
  norm_v = sqrt(2.0*PI) * sigma_v
  norm   = norm_h * norm_v

  ! user output
  if (myrank == 0) then
    print *,"command line arguments:"
    print *,"  smoothing sigma_h , sigma_v                : ",sigma_h,sigma_v
    ! scalelength: approximately S ~ sigma * sqrt(8.0) for a gaussian smoothing
    print *,"  smoothing scalelengths horizontal, vertical: ",sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *,"  input dir : ",trim(input_dir)
    print *,"  output dir: ",trim(output_dir)
    print *
  endif

  ! reads the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED not supported yet for smoothing, please rerun program...'

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *, 'Error number of processors supposed to run on: ',NPROC
      print *, 'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *, 'please rerun with: mpirun -np ',NPROC,' bin/xsmooth_sem .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  call read_mesh_for_init()

  ! allocate arrays for storing the databases
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xix(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xiy(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xiz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for databases'

  ! mesh node locations
  allocate(xstore(NGLOB_AB), &
           ystore(NGLOB_AB), &
           zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB), &
           ispec_is_elastic(NSPEC_AB), &
           ispec_is_poroelastic(NSPEC_AB),stat=ier)
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

  if (ELASTIC_SIMULATION) then
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)

  else if (POROELASTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    rho_vp = 0.0_CUSTOM_REAL
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution_poro(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    DT,model_speed_max,min_resolved_period, &
                                    phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                    LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  else if (ACOUSTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vs'
    rho_vp = sqrt( kappastore / rhostore ) * rhostore
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  endif

  ! for smoothing, we use cell centers to find and locate nearby elements
  !
  ! sets the location of the center of the elements and local points
  allocate(xl(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           yl(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           zl(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           cx0(NSPEC_AB), &
           cy0(NSPEC_AB), &
           cz0(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array xl etc.'

  ! sets element center location
  do ispec = 1, nspec_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          xl(i,j,k,ispec) = xstore(iglob)
          yl(i,j,k,ispec) = ystore(iglob)
          zl(i,j,k,ispec) = zstore(iglob)
        enddo
      enddo
    enddo
    cx0(ispec) = (xl(1,1,1,ispec) + xl(NGLLX,NGLLY,NGLLZ,ispec))/2.0
    cy0(ispec) = (yl(1,1,1,ispec) + yl(NGLLX,NGLLY,NGLLZ,ispec))/2.0
    cz0(ispec) = (zl(1,1,1,ispec) + zl(NGLLX,NGLLY,NGLLZ,ispec))/2.0
  enddo

  ! frees memory
  deallocate(xstore,ystore,zstore)
  deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
  deallocate(ibool)
  deallocate(jacobian)

  ! sets up slices to process
  if (num_interfaces_ext_mesh+1 > MAX_NODE_LIST) stop 'Error number of neighbor interfaces exceeds MAX_NODE_LIST'
  node_list(:) = -1
  do i=1,num_interfaces_ext_mesh
    ! adds neighbors
    node_list(i) = my_neighbours_ext_mesh(i)
  enddo
  ! adds this partition itself
  node_list(num_interfaces_ext_mesh+1) = myrank

  ! user output
  if (myrank == 0) then
  print *,'slices:'
  print *,'  rank:',myrank,'  smoothing slices'
  print *,node_list(1:num_interfaces_ext_mesh+1)
  endif

  !do i=0,sizeprocs-1
  !  if (myrank == i) then
  !    print *,'rank:',myrank,'  smoothing slices'
  !    print *,node_list(1:num_interfaces_ext_mesh+1)
  !    print *
  !  endif
  !enddo

  ! GLL points weights
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

  ! synchronizes
  call synchronize_all()

! loops over slices
! each process reads in his own neighbor slices and gaussian filters the values
  allocate(tk(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array tk and bk'

  tk = 0.0_CUSTOM_REAL
  bk = 0.0_CUSTOM_REAL
  do inum = 1,num_interfaces_ext_mesh+1

    iproc = node_list(inum)

    if (myrank == 0) print *,'  reading slice:',iproc

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
    close(IIN)

    ! allocates arrays
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) stop 'Error allocating array ibool'
    allocate(xstore(NGLOB_N),ystore(NGLOB_N),zstore(NGLOB_N),stat=ier)
    if (ier /= 0) stop 'Error allocating array xstore etc.'
    allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) stop 'Error allocating array jacobian'
    allocate(dummy(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) stop 'Error allocating array dummy'

    ! gets number of point locations (and jacobian, but jacobian not used by default)
    open(unit=IIN,file=trim(prname_lp),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname_lp)
      call exit_mpi(myrank, 'Error reading neighbors external mesh file')
    endif
    read(IIN) NSPEC_N
    read(IIN) NGLOB_N

    ! ibool file
    read(IIN) ibool

    ! global point arrays
    read(IIN) xstore
    read(IIN) ystore
    read(IIN) zstore

    ! reads in jacobian
    read(IIN) dummy ! xix
    read(IIN) dummy ! xiy
    read(IIN) dummy ! xiz
    read(IIN) dummy ! etax
    read(IIN) dummy ! etay
    read(IIN) dummy ! etaz
    read(IIN) dummy ! gammax
    read(IIN) dummy ! gammay
    read(IIN) dummy ! gammaz
    read(IIN) jacobian
    close(IIN)

    deallocate(dummy)

    ! get the location of the center of the elements and local points
    allocate(xx(NGLLX,NGLLY,NGLLZ,NSPEC_N), &
             yy(NGLLX,NGLLY,NGLLZ,NSPEC_N), &
             zz(NGLLX,NGLLY,NGLLZ,NSPEC_N), &
             cx(NSPEC_N), &
             cy(NSPEC_N), &
             cz(NSPEC_N),stat=ier)
    if (ier /= 0) stop 'Error allocating array xx etc.'

    ! sets element center location
    do ispec = 1, nspec_N
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)
            xx(i,j,k,ispec) = xstore(iglob)
            yy(i,j,k,ispec) = ystore(iglob)
            zz(i,j,k,ispec) = zstore(iglob)
          enddo
        enddo
      enddo
      cx(ispec) = (xx(1,1,1,ispec) + xx(NGLLX,NGLLY,NGLLZ,ispec))/2.0
      cy(ispec) = (yy(1,1,1,ispec) + yy(NGLLX,NGLLY,NGLLZ,ispec))/2.0
      cz(ispec) = (zz(1,1,1,ispec) + zz(NGLLX,NGLLY,NGLLZ,ispec))/2.0
    enddo

    deallocate(xstore,ystore,zstore)
    deallocate(ibool)

    ! data file
    write(prname,'(a,i6.6,a)') trim(input_dir)//'proc',iproc,'_'
    local_data_file = trim(prname) // trim(kernel_name) // '.bin'

    open(unit = IIN,file = trim(local_data_file),status='old',action='read',form ='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening data file: ',trim(local_data_file)
      stop 'Error opening data file'
    endif

    allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) stop 'Error allocating dat array'

    read(IIN) dat
    close(IIN)

    if (iproc == myrank)  max_old = maxval(abs(dat(:,:,:,:)))

    ! finds closest elements for smoothing
    !if (myrank==0) print *, '  start looping over elements and points for smoothing ...'

    ! loop over elements to be smoothed in the current slice
    do ispec = 1, nspec_AB
      ! --- only double loop over the elements in the search radius ---
      do ispec2 = 1, nspec_N

        ! calculates horizontal and vertical distance between two element centers
        call get_distance_vec(dist_h,dist_v,cx0(ispec),cy0(ispec),cz0(ispec),&
                          cx(ispec2),cy(ispec2),cz(ispec2))

        ! checks distance between centers of elements
        if (dist_h > sigma_h3 .or. dist_v > sigma_v3) cycle

        ! integration factors
        !factor(:,:,:) = jacobian(:,:,:,ispec2) * wgll_cube(:,:,:)
        !factor(:,:,:) = 1.0_CUSTOM_REAL

        ! loop over GLL points of the elements in current slice (ispec)
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              ! reference location
              ! current point (i,j,k,ispec) location, cartesian coordinates
              x0 = xl(i,j,k,ispec)
              y0 = yl(i,j,k,ispec)
              z0 = zl(i,j,k,ispec)

              ! calculate weights based on gaussian smoothing
              exp_val = 0.0_CUSTOM_REAL
              call smoothing_weights_vec(x0,y0,z0,sigma_h2,sigma_v2,exp_val,&
                      xx(:,:,:,ispec2),yy(:,:,:,ispec2),zz(:,:,:,ispec2))

              ! adds GLL integration weights
              !exp_val(:,:,:) = exp_val(:,:,:) * factor(:,:,:)

              ! adds contribution of element ispec2 to smoothed kernel values
              tk(i,j,k,ispec) = tk(i,j,k,ispec) + sum(exp_val(:,:,:) * dat(:,:,:,ispec2))

              ! normalization, integrated values of gaussian smoothing function
              bk(i,j,k,ispec) = bk(i,j,k,ispec) + sum(exp_val(:,:,:))
            enddo
          enddo
        enddo ! i,j,k
      enddo ! ispec2
    enddo ! ispec

    ! frees arrays
    deallocate(xx,yy,zz)
    deallocate(cx,cy,cz)
    deallocate(jacobian)
    deallocate(dat)

  enddo ! iproc
  if (myrank == 0) print *

  ! normalizes/scaling factor
  if (myrank == 0) print *, 'Scaling values: min/max = ',minval(bk),maxval(bk)

  allocate(dat_smooth(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array dat_smooth'

  dat_smooth(:,:,:,:) = 0.0_CUSTOM_REAL
  do ispec = 1, nspec_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! checks the normalization criterion
          !if (abs(bk(i,j,k,ispec) - norm) > 1.e-4) then
          !  print *, 'Problem norm here --- ', ispec, i, j, k, bk(i,j,k,ispec), norm
          !endif
          if (abs(bk(i,j,k,ispec)) < 1.e-18) then
            print *, 'Problem norm here --- ', ispec, i, j, k, bk(i,j,k,ispec), norm
          endif

          ! normalizes smoothed kernel values by integral value of gaussian weighting
          dat_smooth(i,j,k,ispec) = tk(i,j,k,ispec) / bk(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo !  ispec
  deallocate(tk,bk)

  max_new = maxval(abs(dat_smooth(:,:,:,:)))

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

  ! the maximum value for the smoothed kernel
  norm = max_old
  call max_all_cr(norm, max_old)
  norm = max_new
  call max_all_cr(norm, max_new)
  if (myrank == 0) then
    print *
    print *,'  Maximum data value before smoothing = ', max_old
    print *,'  Maximum data value after smoothing  = ', max_new
    print *
    close(IMAIN)
  endif


  ! stop all the processes and exit
  call finalize_mpi()

end program smooth_sem

!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,y0,z0,sigma_h2,sigma_v2,exp_val,&
                              xx_elem,yy_elem,zz_elem)

  use constants
  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: xx_elem, yy_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,sigma_h2,sigma_v2

  ! local parameters
  integer :: ii,jj,kk
  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: sigma_h2_inv,sigma_v2_inv

  sigma_h2_inv = 1.0_CUSTOM_REAL / sigma_h2
  sigma_v2_inv = 1.0_CUSTOM_REAL / sigma_v2

  do kk = 1, NGLLZ
    do jj = 1, NGLLY
      do ii = 1, NGLLX
        ! point in second slice
        ! gets vertical and horizontal distance
        call get_distance_vec(dist_h,dist_v,x0,y0,z0, &
            xx_elem(ii,jj,kk),yy_elem(ii,jj,kk),zz_elem(ii,jj,kk))

        ! gaussian function
        exp_val(ii,jj,kk) = exp(- sigma_h2_inv*(dist_h*dist_h) - sigma_v2_inv*(dist_v*dist_v))
      enddo
    enddo
  enddo

  end subroutine smoothing_weights_vec


!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_vec(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

! returns vector lengths as distances in radial and horizontal direction
! only for flat earth with z in vertical direction

  use constants
  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,x1,y1,z1

  ! vertical distance
  dist_v = sqrt( (z0-z1)*(z0-z1) )

  ! horizontal distance
  dist_h = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) )

  end subroutine get_distance_vec


