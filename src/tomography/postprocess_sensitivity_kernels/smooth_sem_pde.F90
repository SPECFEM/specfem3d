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

!! Usage:
!! mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR GPU_MODE [SET_ZERO_IN_PML]
!! The last argument (optional, FALSE by default) only functions in case of PML.
!! If set to be TRUE, then the PML domain will be set to zero.
!! If the last argument is omitted, then the command-line arguments
!! are identical to smooth_sem.
!! Ideally, the result of this program should be identical to
!! smooth_sem, but will be tens of times faster

!! GPU is not supported yet, but should be straightforward to
!! implement.
!! Currently, only NGLL=5 is optimized with force inline

!! This implementation solves a diffusion PDE with explicit
!! time stepping. The parameter CFL_CONST is set to guarantee
!! stability. It can be increased in case of instability.

#include "config.fh"

program smooth_sem_pde

  use constants, only: m1, m2, CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, MAX_STRING_LEN,IIN,IOUT
  use specfem_par
  use specfem_par_elastic, only: &
      num_phase_ispec_elastic, &
      nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic
  !use specfem_par_acoustic, only: ispec_is_acoustic, nspec_acoustic
  !use specfem_par_poroelastic, only: ispec_is_poroelastic

  use pml_par, only: is_CPML, NSPEC_CPML, CPML_to_spec

  implicit none

  integer, parameter :: NARGS = 7
  integer, parameter :: PRINT_INFO_PER_STEP = 1000000
  real(kind=CUSTOM_REAL), parameter :: CFL_CONST = 9.0

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat
  !! spherical coordinate !!
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: rotate_r
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    stemp1,stemp2,stemp3, &
    snewtemp1,snewtemp2,snewtemp3, &
    dat_elem
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dat_glob, ddat_glob
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rvol
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: rvol_local
  integer :: i,j,k,l,iglob,ier,ispec,ispec_p,iphase,ispec_irreg

  character(len=MAX_STRING_LEN) :: arg(NARGS)
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  !character(len=MAX_STRING_LEN) :: prname_lp
  character(len=MAX_STRING_LEN*2) :: ks_file
  character(len=MAX_STRING_LEN*2) :: local_data_file

  !character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  !character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: kernel_name
  integer :: num_elements
  real t1,t2,tnow,tlast

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_v, ch, cv, cmax
  real(kind=CUSTOM_REAL) :: min_val, max_val, min_val_glob, max_val_glob

  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  real(kind=CUSTOM_REAL) :: &
    xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl, &
    stemp1l, stemp2l, stemp3l, ddxl, ddyl, ddzl
  integer :: ntstep, istep
  double precision :: weight
  !real(kind=CUSTOM_REAL), dimension(MAX_KERNEL_NAMES) :: max_old,max_new,max_old_all,max_new_all
  !real(kind=CUSTOM_REAL), dimension(MAX_KERNEL_NAMES) :: min_old,min_new,min_old_all,min_new_all
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_vector_ext_mesh_smooth
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_vector_ext_mesh_smooth
  logical :: USE_GPU, SET_ZERO_IN_PML

  integer(kind=8) :: Smooth_container

  ! initializes MPI
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running XSMOOTH_SEM"
  call synchronize_all()
  call cpu_time(t1)

  ! parse command line arguments
  if ((command_argument_count() /= NARGS) .and. &
      (command_argument_count() /= (NARGS - 1)) .and. &
      (command_argument_count() /= (NARGS - 2))) then
    if (myrank == 0) then
        print *,'USAGE:  mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V '
        print *,'KERNEL_NAME INPUT_DIR OUPUT_DIR GPU_MODE [SET_ZERO_IN_PML]'
      stop 'Please check command line arguments'
    endif
  endif
  call synchronize_all()

  if (command_argument_count() == NARGS) then
    do i = 1, NARGS
      call get_command_argument(i,arg(i), status=ier)
    enddo
  else if (command_argument_count() == (NARGS - 1)) then
    do i = 1, NARGS-1
      call get_command_argument(i,arg(i), status=ier)
    enddo
    arg(NARGS) = 'TRUE' ! SET_ZERO_IN_PML = .true. by default
  else if (command_argument_count() == (NARGS - 2)) then
    do i = 1, NARGS-2
      call get_command_argument(i,arg(i), status=ier)
    enddo
    arg(NARGS-1) = 'FALSE' ! USE_GPU = .false. by default
    arg(NARGS) = 'TRUE' ! SET_ZERO_IN_PML = .true. by default
  else
    if (myrank == 0) then
      stop 'Please check command line arguments'
    endif
  endif

  read(arg(1),*) sigma_h
  read(arg(2),*) sigma_v
  kernel_name = arg(3)
  input_dir= arg(4)
  output_dir = arg(5)
  read(arg(6),*) USE_GPU
  read(arg(7),*) SET_ZERO_IN_PML

  !call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  call synchronize_all()
  ! user output
  if (myrank == 0) then
    print *,'command line arguments:'
    print *,'  smoothing sigma_h , sigma_v                : ',sigma_h,sigma_v
    print *,'  input dir : ',trim(input_dir)
    print *,'  output dir: ',trim(output_dir)
    print *,"  GPU_MODE: ", USE_GPU
    print *,"  SET_ZERO_IN_PML: ", SET_ZERO_IN_PML
    print *
  endif

  !if (USE_GPU) then
  !  if (myrank == 0) stop 'GPU mode not available yet'
  !endif

  call initialize_simulation()

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

  !! broadcast distance_min_glob to other processors
  ! check before broadcast
  !if (myrank == 1) print *, 'distance_min_glob = ', distance_min_glob, 'myrank=', myrank
  call bcast_all_singlecr(distance_min_glob)
  ! check after broadcast
  !if (myrank == 1) print *, 'distance_min_glob = ', distance_min_glob, 'myrank=', myrank
  !! spherical coordinate !!
  !allocate(rotate_r(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1013')
  !
  !do ispec = 1, NSPEC_AB
  !  do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
  !    iglob = ibool(i,j,k,ispec)
  !    xl = xstore(iglob)
  !    yl = ystore(iglob)
  !    zl = zstore(iglob)
  !    rl = sqrt(xl*xl+yl*yl+zl*zl)
  !    rotate_r(1,i,j,k,ispec) = xl / rl
  !    rotate_r(2,i,j,k,ispec) = yl / rl
  !    rotate_r(3,i,j,k,ispec) = zl / rl
  !  enddo;enddo;enddo
  !enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!

  deallocate(xstore,ystore,zstore,kappastore,mustore)
  !deallocate(ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic)

  !! determine ch, cv, ntstep
  cmax = distance_min_glob ** 2 / CFL_CONST
  if (sigma_v >= sigma_h) then
    cv = cmax
    ch = cv * (sigma_h ** 2) / (sigma_v ** 2)
  else
    ch = cmax
    cv = ch * (sigma_v ** 2) / (sigma_h ** 2)
  endif
  ntstep = int(ceiling((max(sigma_h,sigma_v)**2)/(2.0*cmax)))

  if (myrank == 0) print *, 'cv=', cv, 'ch=', ch, 'ntstep=', ntstep
  print *, jacobian_regular, xix_regular

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! initialize time iteration
   ! set up GLL points, weights and derivation matrices for reference element
   ! (between -1,1)
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                                  hprime_xx,hprime_yy,hprime_zz, &
                                  hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                                  wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube)

  ! define transpose of derivation matrix
  do j = 1,NGLLY
    do i = 1,NGLLX
      hprime_xxT(j,i) = hprime_xx(i,j)
      hprime_yyT(j,i) = hprime_yy(i,j)
      hprime_zzT(j,i) = hprime_zz(i,j)

      hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
    enddo
  enddo

  ! define a 3D extension in order to be able to force vectorization in the compute_forces routines
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        wgllwgll_yz_3D(i,j,k) = wgllwgll_yz(j,k)
        wgllwgll_xz_3D(i,j,k) = wgllwgll_xz(i,k)
        wgllwgll_xy_3D(i,j,k) = wgllwgll_xy(i,j)
      enddo
    enddo
  enddo

  allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(dat_glob(NGLOB_AB))
  allocate(ddat_glob(NGLOB_AB))
  allocate(buffer_send_vector_ext_mesh_smooth(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_recv_vector_ext_mesh_smooth(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))

   ! prepare assemble array
  allocate(rvol(NGLOB_AB))
  rvol(:) = 0.0
  allocate(rvol_local(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  do ispec = 1, NSPEC_AB
    ispec_irreg = irregular_element_number(ispec)
    do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
      weight =  wxgll(i)*wygll(j)*wzgll(k)
      if (ispec_irreg /= 0) then
        jacobianl = jacobianstore(i,j,k,ispec_irreg)
      else
        jacobianl = jacobian_regular
      endif
      rvol_local(i,j,k,ispec) = real(dble(jacobianl)*weight,kind=CUSTOM_REAL)
      iglob = ibool(i,j,k,ispec)
      rvol(iglob) = rvol(iglob) + rvol_local(i,j,k,ispec)
    enddo;enddo;enddo
  enddo
  call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rvol, &
                                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                    my_neighbors_ext_mesh)
  rvol(:) = 1.0 / rvol(:)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read in data to be smoothed
  ! data file
  write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'
  local_data_file = trim(prname) // trim(kernel_name) // '.bin'

  open(unit = IIN,file = trim(local_data_file),status='old',action='read', &
       form ='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening data file: ',trim(local_data_file)
    stop 'Error opening data file'
  endif

  read(IIN) dat
  close(IIN)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! project
  dat_glob(:) = 0.0
  do ispec = 1, NSPEC_AB
    do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
      iglob = ibool(i,j,k,ispec)
      dat_glob(iglob) = dat_glob(iglob) + dat(i,j,k,ispec) * rvol_local(i,j,k,ispec)
    enddo;enddo;enddo
  enddo

  call assemble_MPI_send_smooth(NPROC,NGLOB_AB, &
                                dat_glob,buffer_send_vector_ext_mesh_smooth, &
                                buffer_recv_vector_ext_mesh_smooth, &
                                num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                my_neighbors_ext_mesh, &
                                request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
  call assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB, &
                                 dat_glob,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh, &
                                 max_nibool_interfaces_ext_mesh, &
                                 nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                 request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                 my_neighbors_ext_mesh,myrank)

  if (myrank == 0) print *, 'Before smoothing: '

  dat_glob(:) = dat_glob(:) * rvol(:)

  min_val = minval(dat_glob)
  max_val = maxval(dat_glob)
  call min_all_cr(min_val, min_val_glob)
  call max_all_cr(max_val, max_val_glob)
  if (myrank == 0) then
    print *, '  '//trim(kernel_name)
    print *, '    minval:', min_val_glob
    print *, '    maxval:', max_val_glob
    if (myrank == 0) call cpu_time(tlast)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! broadcast glob array back to local array
  do ispec = 1, NSPEC_AB
    do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
      iglob = ibool(i,j,k,ispec)
      dat(i,j,k,ispec) = dat_glob(iglob)
    enddo;enddo;enddo
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ddat_glob(:) = 0.0

  if (USE_GPU .and. GPU_MODE) then
  ! both GPU flags in command line and par file needs to be turned on
    ! user output
    call synchronize_all()
    if (myrank == 0) then
      print *, "preparing fields and constants on GPU devices"
    endif
    call synchronize_all()

    ! prepares general fields on GPU
    ! add GPU support for the C-PML routines
    call prepare_constants_smooth_device(Mesh_pointer, &
                                         NGLLX, NSPEC_AB, NGLOB_AB, NSPEC_IRREGULAR,irregular_element_number, &
                                         xixstore,xiystore,xizstore, &
                                         etaxstore,etaystore,etazstore, &
                                         gammaxstore,gammaystore,gammazstore, &
                                         xix_regular,jacobian_regular,ibool, &
                                         num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
                                         nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
                                         hprime_xx,hprimewgll_xx, &
                                         wgllwgll_xy, wgllwgll_xz, wgllwgll_yz, &
                                         myrank, &
                                         PML_CONDITIONS)


    if (myrank == 0) then
      print *, "preparing smoothing parameters on GPU"
    endif

    call prepare_smooth_pde_gpu(Mesh_pointer, &
                                Smooth_container, &
                                dat_glob, &
                                rvol, &
                                phase_ispec_inner_elastic, &
                                num_phase_ispec_elastic, &
                                CPML_to_spec, &
                                NSPEC_CPML, &
                                cv, &
                                ch)
  endif

  if (myrank == 0) then
    print *, "start stepping"
  endif

  do istep = 1, ntstep
    !if (myrank == 0) print *, istep
    do iphase = 1,2
      if (iphase == 1) then
        num_elements = nspec_outer_elastic
      else
        num_elements = nspec_inner_elastic
      endif

      if (USE_GPU .and. GPU_MODE) then
        call compute_update_element_smooth_pde_gpu(Mesh_pointer, &
                                                   Smooth_container, &
                                                   iphase, &
                                                   num_elements)
      else
        do ispec_p = 1,num_elements
          ispec = phase_ispec_inner_elastic(ispec_p,iphase)
          do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dat_elem(i,j,k) = dat_glob(iglob)
          enddo; enddo; enddo
          select case (NGLLX)
          case (5)
            call mxm5_single(hprime_xx,m1,dat_elem,stemp1,m2)
            call mxm5_3dmat_single(dat_elem,m1,hprime_xxT,m1, &
                                   stemp2,NGLLX)
            call mxm5_single(dat_elem,m2,hprime_xxT,stemp3,m1)
          case default
            !! derivative along xi, eta, zeta
            do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
              stemp1l = 0._CUSTOM_REAL
              stemp2l = 0._CUSTOM_REAL
              stemp3l = 0._CUSTOM_REAL
              ! we can merge the loops because NGLLX == NGLLY == NGLLZ
              do l = 1,NGLLX
                stemp1l = stemp1l + dat_elem(l,j,k)*hprime_xx(i,l)
                stemp2l = stemp2l + dat_elem(i,l,k)*hprime_yy(j,l)
                stemp3l = stemp3l + dat_elem(i,j,l)*hprime_zz(k,l)
              enddo
              stemp1(i,j,k) = stemp1l
              stemp2(i,j,k) = stemp2l
              stemp3(i,j,k) = stemp3l
            enddo;enddo;enddo
          end select

          ispec_irreg = irregular_element_number(ispec)
          if (ispec_irreg /= 0) then
            do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
              xixl = xixstore(i,j,k,ispec_irreg)
              xiyl = xiystore(i,j,k,ispec_irreg)
              xizl = xizstore(i,j,k,ispec_irreg)
              etaxl = etaxstore(i,j,k,ispec_irreg)
              etayl = etaystore(i,j,k,ispec_irreg)
              etazl = etazstore(i,j,k,ispec_irreg)
              gammaxl = gammaxstore(i,j,k,ispec_irreg)
              gammayl = gammaystore(i,j,k,ispec_irreg)
              gammazl = gammazstore(i,j,k,ispec_irreg)
              jacobianl = jacobianstore(i,j,k,ispec_irreg)

              ! derivatives along x, y, z
              ddxl = xixl*stemp1(i,j,k) + etaxl*stemp2(i,j,k) + &
                     gammaxl*stemp3(i,j,k)
              ddyl = xiyl*stemp1(i,j,k) + etayl*stemp2(i,j,k) + &
                     gammayl*stemp3(i,j,k)
              ddzl = xizl*stemp1(i,j,k) + etazl*stemp2(i,j,k) + &
                     gammazl*stemp3(i,j,k)
              !! spherical coordinate
              !rxl = rotate_r(1,i,j,k,ispec)
              !ryl = rotate_r(2,i,j,k,ispec)
              !rzl = rotate_r(3,i,j,k,ispec)
              !stemp1(i,j,k) = ((cv-ch) * (rxl*xixl+ryl*xiyl+rzl*xizl) * &
              !  (rxl*ddxl+ryl*ddyl+rzl*ddzl) +&
              !  ch * (xixl*ddxl+xiyl*ddyl+xizl*ddzl)) * jacobianl
              !stemp2(i,j,k) = ((cv-ch) * (rxl*etaxl+ryl*etayl+rzl*etazl) * &
              !  (rxl*ddxl+ryl*ddyl+rzl*ddzl) +&
              !  ch * (etaxl*ddxl+etayl*ddyl+etazl*dzl)) * jacobianl
              !stemp3(i,j,k) = ((cv-ch) * (rxl*gammaxl+ryl*gammayl+rzl*gammazl)* &
              !  (rxl*ddxl+ryl*ddyl+rzl*ddzl) +&
              !  ch * (gammaxl*ddxl+gammayl*ddyl+gammazl*ddzl)) * jacobianl
              !! Cartesian coordinate
              stemp1(i,j,k) = ((cv-ch) * xizl * ddzl +&
                ch * (xixl*ddxl+xiyl*ddyl+xizl*ddzl)) * jacobianl
              stemp2(i,j,k) = ((cv-ch) * etazl * ddzl +&
                ch * (etaxl*ddxl+etayl*ddyl+etazl*ddzl)) * jacobianl
              stemp3(i,j,k) = ((cv-ch) * gammazl * ddzl +&
                ch * (gammaxl*ddxl+gammayl*ddyl+gammazl*ddzl)) * jacobianl
            enddo;enddo;enddo
          else
            xixl = xix_regular
            jacobianl = jacobian_regular
            do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
              stemp1(i,j,k) = (&
                ch * xixl*xixl*stemp1(i,j,k)) * jacobianl
              stemp2(i,j,k) = (&
                ch * xixl*xixl*stemp2(i,j,k)) * jacobianl
              stemp3(i,j,k) = (&
                cv * xixl*xixl*stemp3(i,j,k)) * jacobianl
            enddo;enddo;enddo
          endif

          select case (NGLLX)
          case (5)
            call mxm5_single(hprimewgll_xxT,m1,stemp1,snewtemp1,m2)
            call mxm5_3dmat_single(stemp2,m1,hprimewgll_xx,m1, &
                                   snewtemp2,NGLLX)
            call mxm5_single(stemp3,m2,hprimewgll_xx,snewtemp3,m1)
          case default
            do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
              stemp1l = 0._CUSTOM_REAL
              stemp2l = 0._CUSTOM_REAL
              stemp3l = 0._CUSTOM_REAL
              ! we can merge these loops because NGLLX = NGLLY = NGLLZ
              do l = 1,NGLLX
                stemp1l = stemp1l + stemp1(l,j,k) * hprimewgll_xx(l,i)
                stemp2l = stemp2l + stemp2(i,l,k) * hprimewgll_yy(l,j)
                stemp3l = stemp3l + stemp3(i,j,l) * hprimewgll_zz(l,k)
              enddo
              ! to be compatible with matrix versions from above
              snewtemp1(i,j,k) = stemp1l
              snewtemp2(i,j,k) = stemp2l
              snewtemp3(i,j,k) = stemp3l
            enddo;enddo;enddo
          end select

          do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            ddat_glob(iglob) = ddat_glob(iglob) - (&
                            wgllwgll_yz_3D(i,j,k) * snewtemp1(i,j,k)+&
                            wgllwgll_xz_3D(i,j,k) * snewtemp2(i,j,k)+&
                            wgllwgll_xy_3D(i,j,k) * snewtemp3(i,j,k))
          enddo;enddo;enddo
        enddo  ! ispec_p = 1, num_elements
      endif ! GPU

      !! assemble MPI
      if (iphase == 1) then
        if (USE_GPU .and. GPU_MODE) then
          call transfer_boun_dat_smooth_pde_from_device(Mesh_pointer, &
                                                        Smooth_container, &
                                                        buffer_send_vector_ext_mesh_smooth)
          call assemble_MPI_send_smooth_cuda(NPROC, &
                                      buffer_send_vector_ext_mesh_smooth, &
                                      buffer_recv_vector_ext_mesh_smooth, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh, &
                                      request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
        else
          call assemble_MPI_send_smooth(NPROC,NGLOB_AB, &
                                      ddat_glob,buffer_send_vector_ext_mesh_smooth, &
                                      buffer_recv_vector_ext_mesh_smooth, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh, &
                                      request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
        endif
      else
        if (USE_GPU .and. GPU_MODE) then
          call assemble_MPI_write_smooth_cuda(NPROC, &
                                       Mesh_pointer, Smooth_container, &
                                       buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh, &
                                       max_nibool_interfaces_ext_mesh, &
                                       request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
        else
          call assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB, &
                                       ddat_glob,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh, &
                                       max_nibool_interfaces_ext_mesh, &
                                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                       request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                       my_neighbors_ext_mesh,myrank)
        endif
      endif
      !!!!!!!!!!!!!!!!!
    enddo !iphase = 1,2

    if (USE_GPU .and. GPU_MODE) then
      call kernel_3_smooth_pde_cuda(Mesh_pointer, Smooth_container)
    else
      ddat_glob(:) = ddat_glob(:) * rvol(:)
    endif
    !! update
    if (USE_GPU .and. GPU_MODE) then
      call update_dat_smooth_pde_cuda(Mesh_pointer, Smooth_container)
    else
      dat_glob(:) = dat_glob(:) + ddat_glob(:)
      ddat_glob(:) = 0.0
    endif

    !! info
    if (mod(istep, PRINT_INFO_PER_STEP) == 0) then
      if (myrank == 0) print *, 'Step:', istep
      if (USE_GPU .and. GPU_MODE) then
        call get_norm_smooth_pde_from_device(Mesh_pointer, Smooth_container, &
                                             max_val, 0)
      else
        max_val = maxval(abs(dat_glob))
      endif
      call max_all_cr(max_val, max_val_glob)
      if (myrank == 0) then
        print *, '  '//trim(kernel_name)
        print *, '    max abs val:', max_val_glob
        call cpu_time(tnow)
        print *, 'time since last message:', tnow-tlast
        call cpu_time(tlast)
      endif
    endif

    !!!!!!!!!!!!!
    !! mute everything in the PML domain
    if (SET_ZERO_IN_PML .and. PML_CONDITIONS) then
      if (USE_GPU .and. GPU_MODE) then
        call zero_pml_smooth_pde_cuda(Mesh_pointer, Smooth_container)
      else
        do ispec = 1, NSPEC_AB
          do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            if (is_CPML(ispec)) then
              dat_glob(iglob) = 0.0
            endif
          enddo;enddo;enddo
        enddo
      endif
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call synchronize_all()
  enddo

  if (USE_GPU .and. GPU_MODE) then
    call transfer_dat_smooth_pde_from_device(Mesh_pointer, Smooth_container, &
                                             dat_glob)
  endif

  call synchronize_all()
  call cpu_time(t2)
  if (myrank == 0) &
    print *, 'Computation time with PDE-based smoothing on CPU:', t2-t1

  !! statistics
  min_val = minval(dat_glob)
  max_val = maxval(dat_glob)
  call min_all_cr(min_val, min_val_glob)
  call max_all_cr(max_val, max_val_glob)
  if (myrank == 0) then
    print *, 'After smoothing:'
    print *, '  '//trim(kernel_name)
    print *, '    minval:', min_val_glob
    print *, '    maxval:', max_val_glob
  endif

  !! broadcast glob array back to local array
  do ispec = 1, NSPEC_AB
    do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
      iglob = ibool(i,j,k,ispec)
      dat(i,j,k,ispec) = dat_glob(iglob)
    enddo;enddo;enddo
  enddo

  !! output
  write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name)//'_smooth.bin'
  open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening smoothed kernel file'
  write(IOUT) dat
  close(IOUT)
  if (myrank == 0) print *,'written: ',trim(ks_file)

  deallocate(ibool,irregular_element_number)
  deallocate(xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
             gammaxstore,gammaystore,gammazstore,jacobianstore)
  !! spherical coordinate !!
  !deallocate(rotate_r)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  deallocate(dat, dat_glob, ddat_glob)
  deallocate(buffer_send_vector_ext_mesh_smooth, &
             buffer_recv_vector_ext_mesh_smooth)
  deallocate(rvol, rvol_local)

  call finalize_mpi()

  contains

  subroutine mxm5_single(A,n1,B,C,n3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_single
#else
! cray
!DIR$ INLINEALWAYS mxm5_single
#endif

! two-dimensional arrays (25,5)/(5,25)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_single

  subroutine mxm5_3dmat_single(A,n1,B,n2,C,n3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3dmat_single
#else
! cray
!DIR$ INLINEALWAYS mxm5_3dmat_single
#endif

! three-dimensional arrays (5,5,5) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
!DIR$ SIMD
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3dmat_single

end program smooth_sem_pde

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_send_smooth(NPROC,NGLOB_AB, &
                                      array_val,buffer_send_vector_ext_mesh_smooth, &
                                      buffer_recv_vector_ext_mesh_smooth, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh, &
                                      request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    ! sends data

  use constants, only: CUSTOM_REAL, itag

  implicit none

  integer :: NPROC
  integer :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh_smooth,buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: &
    nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh):: &
    ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer ipoin,iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_vector_ext_mesh_smooth(ipoin,iinterface) = &
             array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_vector_ext_mesh(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_send_smooth

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_send_smooth_cuda(NPROC, &
                                           buffer_send_vector_ext_mesh_smooth, &
                                           buffer_recv_vector_ext_mesh_smooth, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    ! sends data

  use constants, only: CUSTOM_REAL, itag

  implicit none

  integer :: NPROC

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh_smooth,buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: &
    nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer :: iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_vector_ext_mesh(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_send_smooth_cuda

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB, &
                                       array_val,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh, &
                                       max_nibool_interfaces_ext_mesh, &
                                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                       request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                       my_neighbors_ext_mesh,myrank)

! waits for data to receive and assembles

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: NPROC
  integer :: NGLOB_AB
! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,myrank

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh):: &
    ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbors_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
    mybuffer
  integer :: ipoin,iinterface,iglob
  logical :: need_add_my_contrib

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if (NPROC == 1) return

! move interface values of array_val to local buffers
  do iinterface = 1, num_interfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      mybuffer(ipoin,iinterface) = array_val(iglob)
     ! set them to zero right away to avoid counting it more than once during
     ! assembly:
     ! buffers of higher rank get zeros on nodes shared with current buffer
      array_val(iglob) = 0._CUSTOM_REAL
    enddo
  enddo

! wait for communications completion (recv)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_recv_vector_ext_mesh(iinterface))
  enddo

! adding all contributions in order of processor rank
  need_add_my_contrib = .true.
  do iinterface = 1, num_interfaces_ext_mesh
    if (need_add_my_contrib .and. myrank < my_neighbors_ext_mesh(iinterface)) &
      call add_my_contrib()
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      array_val(iglob) = array_val(iglob) + &
        buffer_recv_vector_ext_mesh_smooth(ipoin,iinterface)
    enddo
  enddo
  if (need_add_my_contrib) call add_my_contrib()

! wait for communications completion (send)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_send_vector_ext_mesh(iinterface))
  enddo

  contains

    subroutine add_my_contrib()

    integer :: my_iinterface,my_ipoin

    do my_iinterface = 1, num_interfaces_ext_mesh
      do my_ipoin = 1, nibool_interfaces_ext_mesh(my_iinterface)
        iglob = ibool_interfaces_ext_mesh(my_ipoin,my_iinterface)
        array_val(iglob) = array_val(iglob) + &
          mybuffer(my_ipoin,my_iinterface)
      enddo
    enddo
    need_add_my_contrib = .false.

    end subroutine add_my_contrib

  end subroutine assemble_MPI_w_ord_smooth

  subroutine assemble_MPI_write_smooth_cuda(NPROC,Mesh_pointer, &
                                            Smooth_container, &
                                            buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

! waits for data to receive and assembles

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: NPROC

  integer(kind=8), intent(in) :: Mesh_pointer, Smooth_container

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer :: iinterface

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if (NPROC == 1) return

! wait for communications completion (recv)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_recv_vector_ext_mesh(iinterface))
  enddo

  call transfer_asmbl_dat_smooth_pde_from_device(Mesh_pointer, &
                                           Smooth_container, &
                                           buffer_recv_vector_ext_mesh_smooth)
! wait for communications completion (send)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_send_vector_ext_mesh(iinterface))
  enddo

  end subroutine assemble_MPI_write_smooth_cuda
