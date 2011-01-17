!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

! this program can be used for smoothing a kernel,
! where it smooths files with a given input kernel name:
!
! compile with:
!     mpif90 -o bin/xsmooth_vol_data -Wall src/shared/smooth_vol_data.f90 \
!               obj/spec/gll_library.o obj/spec/read_parameter_file.o obj/spec/read_value_parameters.o \
!               obj/spec/get_value_parameters.o obj/spec/param_reader.o obj/spec/parallel.o obj/spec/exit_mpi.o
!     or
!     make xsmooth_vol_data
!
! Usage:
!   mpirun -np nprocs ./xsmooth_vol_data filename input_dir output_dir sigma_h sigma_v
!   e.g.
!   mpirun -np 8 ./smooth_vol_data alpha_kernel DATABASES_MPI/ OUTPUT_SUM/ 100 20
!
! where:
!   sigma_h                - gaussian width for horizontal smoothing
!   sigma_v                - gaussian width for vertical smoothing
!   kernel_file_name  - takes file with this kernel name,
!                                     e.g. "alpha_kernel"
!   scratch_file_dir     - directory containing kernel files,
!                                      e.g. proc***_alpha_kernel.bin
!   topo_dir                - directory containing mesh topo files:
!                                       proc***_solver_data1.bin, proc***_solver_data2.bin
! outputs:
!    puts the resulting, smoothed kernel files into the same directory as scratch_file_dir/
!    with a file ending "proc***_kernel_smooth.bin"

program smooth_vol_data

! this is the embarassingly-parallel program that smooths any specfem function (primarily the kernels)
! that has the dimension of (NGLLX,NGLLY,NGLLZ,NSPEC)
!
! NOTE:  smoothing can be different in vertical & horizontal directions; mesh is in Cartesian geometry.
!              algorithm uses vertical as Z, horizontal as X/Y direction

  implicit none
  include "constants.h"
  include 'mpif.h'
  include "precision.h"

 ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat,dat_smooth

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dummy
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: dummy_1
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: dummy_2
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: dummy_3
  integer, dimension(:),allocatable :: idummy
  logical, dimension(:),allocatable :: ldummy
  integer, dimension(:,:,:),allocatable :: idummy_3

  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  integer, dimension(:,:,:,:),allocatable :: ibool

  integer :: num_interfaces_ext_mesh
  integer, dimension(:),allocatable:: my_neighbours_ext_mesh

  integer :: NSPEC_AB, NGLOB_AB
  integer :: NSPEC_N, NGLOB_N

  integer :: i,j,k,ios,it,iglob,ier,ispec2,ispec
  integer :: iproc, node_list(300)
  integer :: np

  character(len=256) :: arg(7), filename, indir, outdir
  character(len=256) :: prname, prname_lp
  character(len=256) :: local_data_file

  ! for read_parameter_files
  double precision :: DT
  double precision :: HDUR_MOVIE
  integer :: NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP, &
            UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer :: NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  integer :: NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
            USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION
  logical :: ATTENUATION,USE_OLSEN_ATTENUATION, &
            OCEANS
  logical :: ABSORBING_CONDITIONS,SAVE_FORWARD
  logical :: ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION
  character(len=256) LOCAL_PATH

  ! smoothing parameters
  character(len=256) :: ks_file

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3
  real(kind=CUSTOM_REAL) :: x0, y0, z0, norm, norm_h, norm_v, max_old, max_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: exp_val !,factor

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: tk, bk, jacobian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: xl, yl, zl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: xx, yy, zz

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx0, cy0, cz0
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx, cy, cz

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll

  ! array with all the weights in the cube
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: element_size

  logical :: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION
  integer :: idummy_a
  integer :: myrank,sizeprocs

!------------------

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if (myrank == 0) print*,"smooth:"
  call mpi_barrier(MPI_COMM_WORLD,ier)

  ! reads arguments
  do i = 1, 5
    call getarg(i,arg(i))
    if (i <= 5 .and. trim(arg(i)) == '') then
      print *, 'Usage: '
      print *, '        xsmooth_data filename input_dir output_dir sigma_h sigma_v'
      print *, '    or '
      print *, '        xsmooth_data filename input_dir output_dir sigma_h sigma_v'
      print *
      print *, ' possible filenames are '
      print *, '   rho_vp, rho_vs, kappastore, mustore, alpha_kernel, etc'
      print *
      print *, '   that are stored in the local directory as real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,NSPEC_AB)  '
      print *, '   in filename.bin'
      print *
      print *, ' files have been collected in input_dir, output mesh file goes to output_dir '
      print *
      stop ' Reenter command line options'
    endif
  enddo

  ! gets arguments
  filename = arg(1)
  indir= arg(2)
  outdir = arg(3)
  read(arg(4),*) sigma_h
  read(arg(5),*) sigma_v

  ! initializes lengths
  sigma_h2 = 2.0 * sigma_h ** 2  ! factor two for gaussian distribution with standard variance sigma
  sigma_v2 = 2.0 * sigma_v ** 2

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
  if( myrank == 0 ) then
    print*,"defaults:"
    print*,"  element size   : ",element_size
    print*,"  smoothing sigma_h , sigma_v: ",sigma_h,sigma_v
    print*,"  in dir : ",trim(indir)
    print*,"  out dir: ",trim(outdir)
  endif

  ! needs local_path for mesh files
  call read_parameter_file( NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        OCEANS,ANISOTROPY,ABSORBING_CONDITIONS, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY)

  if (sizeprocs /= NPROC) call exit_mpi(myrank,'Error total number of slices')
  call mpi_barrier(MPI_COMM_WORLD,ier)

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

! ---------------------

! reads mesh file

  ACOUSTIC_SIMULATION = .false.
  ELASTIC_SIMULATION = .false.
  POROELASTIC_SIMULATION = .false.
  np = 0

  iproc = myrank
  if( myrank == 0 ) then
    print *, ' '
    !print *, 'Reading partition ', iproc
  endif

  ! gets number of elements and global points for this partition
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'//'external_mesh.bin'
  !if( myrank == 0 ) print*,trim(prname_lp)

  open(unit=27,file=trim(prname_lp),&
          status='old',action='read',form='unformatted',iostat=ios)
  if( ios /= 0 ) call exit_mpi(myrank, 'error reading external mesh file')

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

  ! reads in jacobian
  allocate(dummy(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  read(27) dummy ! xix
  read(27) dummy ! xiy
  read(27) dummy ! xiz
  read(27) dummy ! etax
  read(27) dummy ! etay
  read(27) dummy ! etaz
  read(27) dummy ! gammax
  read(27) dummy ! gammay
  read(27) dummy ! gammaz
  read(27) jacobian

  ! reads in partiton neighbors
  read(27) dummy ! kappastore
  read(27) dummy ! mustore

  allocate(ldummy(NSPEC_AB))
  read(27) ldummy ! ispec_is_acoustic
  call any_all_l( ANY(ldummy), ACOUSTIC_SIMULATION )

  read(27) ldummy ! ispec_is_elastic
  call any_all_l( ANY(ldummy), ELASTIC_SIMULATION )

  read(27) ldummy ! ispec_is_poroelastic
  call any_all_l( ANY(ldummy), POROELASTIC_SIMULATION )

  deallocate(ldummy)

  allocate(dummy_1(NGLOB_AB))
  if( ACOUSTIC_SIMULATION ) then
    read(27) dummy_1 ! rmass_acoustic
    read(27) dummy ! rhostore
  endif
  if( ELASTIC_SIMULATION ) then
    read(27) dummy_1 ! rmass
    if( OCEANS ) then
      read(27) dummy_1 ! rmass_ocean_load
    endif
    read(27) dummy ! rho_vp
    read(27) dummy ! rho_vs
  endif
  deallocate(dummy_1)
  deallocate(dummy)

  read(27) idummy_a ! num_abs_boundary_faces
  allocate(idummy(idummy_a))
  allocate(idummy_3(3,NGLLSQUARE,idummy_a))
  allocate(dummy_2(NGLLSQUARE,idummy_a))
  allocate(dummy_3(NDIM,NGLLSQUARE,idummy_a))
  read(27) idummy ! abs_boundary_ispec
  read(27) idummy_3 ! abs_boundary_ijk
  read(27) dummy_2 ! abs_boundary_jacobian2Dw
  read(27) dummy_3 ! abs_boundary_normal
  deallocate( idummy,idummy_3,dummy_2,dummy_3)

  read(27) idummy_a ! num_free_surface_faces
  allocate(idummy(idummy_a))
  allocate(idummy_3(3,NGLLSQUARE,idummy_a))
  allocate(dummy_2(NGLLSQUARE,idummy_a))
  allocate(dummy_3(NDIM,NGLLSQUARE,idummy_a))
  read(27) idummy ! free_surface_ispec
  read(27) idummy_3 ! free_surface_ijk
  read(27) dummy_2 ! free_surface_jacobian2Dw
  read(27) dummy_3 ! free_surface_normal
  deallocate( idummy,idummy_3,dummy_2,dummy_3)

  read(27) idummy_a ! num_coupling_ac_el_faces
  allocate(idummy(idummy_a))
  allocate(idummy_3(3,NGLLSQUARE,idummy_a))
  allocate(dummy_2(NGLLSQUARE,idummy_a))
  allocate(dummy_3(NDIM,NGLLSQUARE,idummy_a))
  read(27) idummy ! coupling_ac_el_ispec
  read(27) idummy_3 ! coupling_ac_el_ijk
  read(27) dummy_2 ! coupling_ac_el_jacobian2Dw
  read(27) dummy_3 ! coupling_ac_el_normal
  deallocate( idummy,idummy_3,dummy_2,dummy_3)

  read(27) num_interfaces_ext_mesh ! num_interfaces_ext_mesh
  read(27) idummy_a ! max_nibool_interfaces_ext_mesh
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh))
  read(27) my_neighbours_ext_mesh

  close(27)

  ! get the location of the center of the elements and local points
  allocate(xl(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(yl(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(zl(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(cx0(NSPEC_AB))
  allocate(cy0(NSPEC_AB))
  allocate(cz0(NSPEC_AB))
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
  deallocate(xstore,ystore,zstore)
  deallocate(ibool)
  deallocate(jacobian)

  ! sets up slices to process
  do i=1,num_interfaces_ext_mesh
    ! adds neighbors
    node_list(i) = my_neighbours_ext_mesh(i)
  enddo
  ! adds this partition itself
  node_list(num_interfaces_ext_mesh+1) = myrank

  ! user output
  if(myrank == 0) then
    print*
    print*,'  rank:',myrank,'  smoothing slices'
    print*,node_list(1:num_interfaces_ext_mesh+1)
  endif
  !do i=0,sizeprocs-1
  !  if( myrank == i ) then
  !    print*,'rank:',myrank,'  smoothing slices'
  !    print*,node_list(1:num_interfaces_ext_mesh+1)
  !    print*
  !  endif
  !enddo

  ! synchronizes
  call mpi_barrier(MPI_COMM_WORLD,ier)


!----------------------

! loops over slices
! each process reads in his own neighbor slices and gaussian filters the values

  allocate(tk(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  tk = 0.0_CUSTOM_REAL
  bk = 0.0_CUSTOM_REAL
  do it=1,num_interfaces_ext_mesh+1

    iproc = node_list(it)

    if( myrank == 0 ) print*,'  reading slice:',iproc

    ! gets number of elements and global points for this partition
    write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',iproc,'_'//'external_mesh.bin'
    !if( myrank == 0 ) print*,trim(prname_lp)

    open(unit=27,file=trim(prname_lp),&
            status='old',action='read',form='unformatted',iostat=ios)
    if( ios /= 0 ) call exit_mpi(myrank, 'error reading neighbors external mesh file')

    read(27) NSPEC_N
    read(27) NGLOB_N

    ! ibool file
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_N))
    read(27) ibool

    ! global point arrays
    allocate(xstore(NGLOB_N),ystore(NGLOB_N),zstore(NGLOB_N))
    read(27) xstore
    read(27) ystore
    read(27) zstore

    ! reads in jacobian
    allocate(dummy(NGLLX,NGLLY,NGLLZ,NSPEC_N), &
            jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_N))
    read(27) dummy ! xix
    read(27) dummy ! xiy
    read(27) dummy ! xiz
    read(27) dummy ! etax
    read(27) dummy ! etay
    read(27) dummy ! etaz
    read(27) dummy ! gammax
    read(27) dummy ! gammay
    read(27) dummy ! gammaz
    read(27) jacobian
    close(27)
    deallocate(dummy)

    ! get the location of the center of the elements and local points
    allocate(xx(NGLLX,NGLLY,NGLLZ,NSPEC_N))
    allocate(yy(NGLLX,NGLLY,NGLLZ,NSPEC_N))
    allocate(zz(NGLLX,NGLLY,NGLLZ,NSPEC_N))
    allocate(cx(NSPEC_N))
    allocate(cy(NSPEC_N))
    allocate(cz(NSPEC_N))
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
    write(prname,'(a,i6.6,a)') trim(indir)//'proc',iproc,'_'
    local_data_file = trim(prname) // trim(filename) // '.bin'
    !if( myrank == 0 ) print*,trim(local_data_file)

    open(unit = 28,file = trim(local_data_file),status='old',&
          action='read',form ='unformatted',iostat=ios)
    if (ios /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ios)
    if( ios /= 0 ) stop 'error allocating dat array'
    read(28) dat
    close(28)

    if( iproc == myrank )  max_old = maxval(abs(dat(:,:,:,:)))


    ! finds closest elements for smoothing
    !if(myrank==0) print*, '  start looping over elements and points for smoothing ...'


    ! loop over elements to be smoothed in the current slice
    do ispec = 1, nspec_AB
      ! --- only double loop over the elements in the search radius ---
      do ispec2 = 1, nspec_N

        ! calculates horizontal and vertical distance between two element centers
        call get_distance_vec(dist_h,dist_v,cx0(ispec),cy0(ispec),cz0(ispec),&
                          cx(ispec2),cy(ispec2),cz(ispec2))

        ! checks distance between centers of elements
        if ( dist_h > sigma_h3 .or. dist_v > sigma_v3 ) cycle

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

  ! normalizes
  !if(myrank==0) print*, 'normalizes values ...'
  allocate(dat_smooth(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  dat_smooth = 0.0_CUSTOM_REAL
  do ispec = 1, nspec_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! checks the normalization criterion
          !if (abs(bk(i,j,k,ispec) - norm) > 1.e-4 ) then
          !  print *, 'Problem norm here --- ', ispec, i, j, k, bk(i,j,k,ispec), norm
          !endif
          if (abs(bk(i,j,k,ispec)) < 1.e-18 ) then
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

!------------------

! file output

  ! smoothed kernel file name
  write(ks_file,'(a,i6.6,a)') trim(outdir)//'/proc',myrank,'_'//trim(filename)//'_smooth.bin'

  open(11,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  if( ier /= 0 ) stop 'error opening smoothed kernel file'
  write(11) dat_smooth(:,:,:,:)
  close(11)
  if( myrank == 0 ) print *,'written:',trim(ks_file)

  ! frees memory
  deallocate(dat_smooth)

  ! synchronizes
  call mpi_barrier(MPI_COMM_WORLD,ier)

  ! the maximum value for the smoothed kernel
  norm = max_old
  call mpi_reduce(norm,max_old,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  norm = max_new
  call mpi_reduce(norm,max_new,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  if( myrank == 0 ) then
    print *
    print *,'  Maximum data value before smoothing = ', max_old
    print *,'  Maximum data value after smoothing = ', max_new
    print *
  endif

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program smooth_vol_data

!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,y0,z0,sigma_h2,sigma_v2,exp_val,&
                              xx_elem,yy_elem,zz_elem)

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: xx_elem, yy_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,sigma_h2,sigma_v2

  ! local parameters
  integer :: ii,jj,kk
  real(kind=CUSTOM_REAL) :: dist_h,dist_v

  do kk = 1, NGLLZ
    do jj = 1, NGLLY
      do ii = 1, NGLLX
        ! point in second slice
        ! gets vertical and horizontal distance
        call get_distance_vec(dist_h,dist_v,x0,y0,z0, &
            xx_elem(ii,jj,kk),yy_elem(ii,jj,kk),zz_elem(ii,jj,kk))

        ! gaussian function
        exp_val(ii,jj,kk) = exp(- (dist_h*dist_h)/sigma_h2 - (dist_v*dist_v)/sigma_v2)
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

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,x1,y1,z1

  ! vertical distance
  dist_v = sqrt( (z0-z1)*(z0-z1) )

  ! horizontal distance
  dist_h = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) )

  end subroutine get_distance_vec

