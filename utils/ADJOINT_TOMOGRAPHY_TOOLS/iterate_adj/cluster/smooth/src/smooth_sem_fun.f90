program smooth_specfem_function

! this is the embarassingly-parallel program that smooth any specfem function (primarily
! the kernels) that has the dimension of (NGLLX,NGLLY,NGLLZ,NSPEC_MAX), notice that it
! uses the constants.h and precision.h files from the original SEM package, and the
! values_from_mesher.h file from the output of the mesher (or create_header_file),
! therefore, you need to compile it for your specific case

  implicit none
  include 'mpif.h'
  include 'constants.h'
  include 'precision.h'
  include 'values_from_mesher.h'

! ======================================================

  integer, parameter :: NSPEC_MAX=NSPEC_AB
  integer, parameter :: NGLOB_MAX=NGLOB_AB

! only include the neighboring 3 x 3 slices
  integer, parameter :: NSLICES = 3
  integer ,parameter :: NSLICES2 = NSLICES * NSLICES

  character(len=150) :: s_nproc_xi, s_nproc_eta, s_nchunks, s_element_size, s_sigma_h, s_sigma_v
  character(len=150) :: kernel_file_name, scratch_topo_dir, scratch_file_dir
  integer :: sizeprocs, ier,myrank, nproc_xi, nproc_eta, nchunks, ichunk, ixi, ieta, iglob
  integer :: islice(NSLICES2), islice0(NSLICES2), ns

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3

  real(kind=CUSTOM_REAL) :: x0, y0, z0, norm, norm_h, norm_v, element_size, max_old, max_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor, exp_val

  character(len=150) ::  ks_file, reg_name
  character(len=150), dimension(NSLICES2) :: x_file, y_file, z_file, j_file, k_file, i_file

  logical :: global_code

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: kernel, kernel_smooth
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: tk, bk, jacobian, xl, yl, zl, xx, yy, zz
  real(kind=CUSTOM_REAL), dimension(NGLOB_MAX) :: x, y, z
  real(kind=CUSTOM_REAL), dimension(NSPEC_MAX) :: cx0, cy0, cz0, cx, cy, cz

  integer :: i, j, k, ispec, ii, ispec2, nspec(NSLICES2), nglob(NSLICES2)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll

  ! array with all the weights in the cube
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  ! input arguments (nchunks=0 means basin code!)
  call get_command_argument(1,s_nproc_xi)
  call get_command_argument(2,s_nproc_eta)
  call get_command_argument(3,s_nchunks)
  call get_command_argument(4,s_element_size)
  call get_command_argument(5,s_sigma_h)
  call get_command_argument(6,s_sigma_v)
  call get_command_argument(7,kernel_file_name)
  call get_command_argument(8,scratch_file_dir)
  call get_command_argument(9,scratch_topo_dir)

  if (trim(s_nproc_xi) == '' .or. trim(s_nproc_eta) == '' .or. trim(s_nchunks) == '' &
             .or. trim(s_element_size) == '' .or. trim(kernel_file_name) == '' &
             .or. trim(scratch_file_dir) == '' .or. trim(scratch_topo_dir) == '') then
    call exit_MPI(myrank,'Usage: smooth_sem_fun nproc_xi nproc_eta nchunks element_size_on_surface(km) sigma_h(km) sigma_v(km) kernel_file_name scratch_file_dir scratch_topo_dir')
  endif

  ! read in parameter information
  read(s_nproc_xi,*) nproc_xi
  read(s_nproc_eta,*) nproc_eta
  read(s_nchunks,*) nchunks
  read(s_element_size,*) element_size
  read(s_sigma_h,*) sigma_h
  read(s_sigma_v,*) sigma_v

  if (nchunks == 0) then
    global_code = .false.
    reg_name='_'
    nchunks = 1
  else
    global_code = .true.
    reg_name='_reg1_'
  endif
  if (sizeprocs /= nproc_xi*nproc_eta*nchunks) call exit_mpi(myrank,'Error total number of slices')

  element_size = element_size * 1000  ! e.g. 9 km on the surface, 36 km at CMB

  sigma_h = sigma_h * 1000.0 ! m
  sigma_v = sigma_v * 1000.0 ! m

  sigma_h2 = sigma_h ** 2
  sigma_v2 = sigma_v ** 2

  ! search radius
  sigma_h3 = 3.0  * sigma_h + element_size
  sigma_v3 = 3.0  * sigma_v + element_size

  ! theoretic normal value
  norm_h = 2.0*PI*sigma_h**2
  norm_v = sqrt(2.0*PI) * sigma_v
  norm   = norm_h * norm_v
  !norm = (sqrt(2.0*PI) * sigma) ** 3

  ! GLL points
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

! ---- figure out the neighboring 8 or 7 slices: (ichunk,ixi,ieta) index start at 0------
  ichunk = myrank / (nproc_xi * nproc_eta)
  ieta = (myrank - ichunk * nproc_xi * nproc_eta) / nproc_xi
  ixi = myrank - ichunk * nproc_xi * nproc_eta - ieta * nproc_xi

  ! get the neighboring slices:
  call get_all_eight_slices(ichunk,ixi,ieta, &
             islice0(1),islice0(2),islice0(3),islice0(4),islice0(5),islice0(6),islice0(7),islice0(8), &
             nproc_xi,nproc_eta)

  ! remove the repeated slices (only 8 for corner slices in global case)
  islice(1) = myrank; j = 1
  do i = 1, 8
    if (.not. any(islice(1:i) == islice0(i)) .and. islice0(i) < sizeprocs) then
      j = j + 1
      islice(j) = islice0(i)
    endif
  enddo
  ns = j

  ! read in the topology files of the current and neighboring slices
  do i = 1, ns
    write(x_file(i),'(a,i6.6,a)') trim(scratch_topo_dir)//'/proc',islice(i),trim(reg_name)//'x.bin'
    write(y_file(i),'(a,i6.6,a)') trim(scratch_topo_dir)//'/proc',islice(i),trim(reg_name)//'y.bin'
    write(z_file(i),'(a,i6.6,a)') trim(scratch_topo_dir)//'/proc',islice(i),trim(reg_name)//'z.bin'
    write(j_file(i),'(a,i6.6,a)') trim(scratch_topo_dir)//'/proc',islice(i),trim(reg_name)//'jacobian.bin'
    write(i_file(i),'(a,i6.6,a)') trim(scratch_topo_dir)//'/proc',islice(i),trim(reg_name)//'ibool.bin'
    write(k_file(i),'(a,i6.6,a)') trim(scratch_file_dir)//'/proc',islice(i),trim(reg_name)//trim(kernel_file_name)//'.bin'

    nspec(i) = NSPEC_AB
    nglob(i) = NGLOB_AB
  enddo

  ! read in myrank slice
  open(11,file=x_file(1),status='old',form='unformatted')
  read(11) x(1:nglob(1))
  close(11)
  open(11,file=y_file(1),status='old',form='unformatted')
  read(11) y(1:nglob(1))
  close(11)
  open(11,file=z_file(1),status='old',form='unformatted')
  read(11) z(1:nglob(1))
  close(11)
  open(11,file=i_file(1),status='old',form='unformatted')
  read(11) ibool(:,:,:,1:nspec(1))
  close(11)

  ! get the location of the center of the elements
  do ispec = 1, nspec(1)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          xl(i,j,k,ispec) = x(iglob)
          yl(i,j,k,ispec) = y(iglob)
          zl(i,j,k,ispec) = z(iglob)
        enddo
      enddo
    enddo
    cx0(ispec) = (xl(1,1,1,ispec) + xl(NGLLX,NGLLY,NGLLZ,ispec))/2
    cy0(ispec) = (yl(1,1,1,ispec) + yl(NGLLX,NGLLY,NGLLZ,ispec))/2
    cz0(ispec) = (zl(1,1,1,ispec) + zl(NGLLX,NGLLY,NGLLZ,ispec))/2
  enddo

  if (myrank == 0) write(*,*) 'start looping over elements and points for smoothing ...'

  write(ks_file,'(a,i6.6,a)') trim(scratch_file_dir)//'/proc',myrank,trim(reg_name)//trim(kernel_file_name)//'_smooth.bin'

  tk = 0.0; bk = 0.0; kernel_smooth=0.0

  ! loop over all the slices
  do ii = 1, ns

    ! read in the topology, kernel files, calculate center of elements
    open(11,file=x_file(ii),status='old',form='unformatted')
    read(11) x(1:nglob(ii))
    close(11)
    open(11,file=y_file(ii),status='old',form='unformatted')
    read(11) y(1:nglob(ii))
    close(11)
    open(11,file=z_file(ii),status='old',form='unformatted')
    read(11) z(1:nglob(ii))
    close(11)
    open(11,file=i_file(ii),status='old',form='unformatted')
    read(11) ibool(:,:,:,1:nspec(ii))
    close(11)
    open(11,file=j_file(ii),status='old',form='unformatted')
    read(11) jacobian(:,:,:,1:nspec(ii))
    close(11)
    open(11,file=k_file(ii),status='old',form='unformatted')
    read(11) kernel(:,:,:,1:nspec(ii))
    close(11)

! get the global maximum value of the original kernel file
    if (ii == 1) then
      call mpi_reduce(maxval(abs(kernel(:,:,:,1:nspec(ii)))), max_old, 1, &
                 CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
    endif

    do ispec2 = 1, nspec(ii)
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec2)
            xx(i,j,k,ispec2) = x(iglob)
            yy(i,j,k,ispec2) = y(iglob)
            zz(i,j,k,ispec2) = z(iglob)
          enddo
        enddo
      enddo
      cx(ispec2) = (xx(1,1,1,ispec2) + xx(NGLLX,NGLLZ,NGLLY,ispec2))/2
      cy(ispec2) = (yy(1,1,1,ispec2) + yy(NGLLX,NGLLZ,NGLLY,ispec2))/2
      cz(ispec2) = (zz(1,1,1,ispec2) + zz(NGLLX,NGLLZ,NGLLY,ispec2))/2
    enddo

    if (myrank == 0) write(*,*) 'slice number = ', ii

    ! loop over elements to be smoothed in the current slice
    do ispec = 1, nspec(1)
      if (mod(ispec,100) == 0 .and. myrank == 0) write(*,*) 'ispec=', ispec

      ! --- only double loop over the elements in the search radius ---
      do ispec2 = 1, nspec(ii)

        !if ( sqrt( (cx(ispec2)-cx0(ispec)) **2 + (cy(ispec2)-cy0(ispec)) ** 2 + (cz(ispec2)-cz0(ispec)) ** 2) > sigma3) cycle
        if ( sqrt( (cx(ispec2)-cx0(ispec)) **2 + (cy(ispec2)-cy0(ispec)) ** 2 ) > sigma_h3 .or. sqrt( (cz(ispec2)-cz0(ispec)) ** 2) > sigma_v3 ) cycle

        factor(:,:,:) = jacobian(:,:,:,ispec2) * wgll_cube(:,:,:) ! integration factors

        ! loop over GLL points of the elements in current slice (ispec)
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX

              x0 = xl(i,j,k,ispec); y0 = yl(i,j,k,ispec); z0 = zl(i,j,k,ispec) ! current point (i,j,k,ispec)

              !exp_val(:,:,:) = exp( -((xx(:,:,:,ispec2)-x0)**2+(yy(:,:,:,ispec2)-y0)**2 &
              !          +(zz(:,:,:,ispec2)-z0)**2 )/(2*sigma2) )*factor(:,:,:)

              exp_val(:,:,:) = exp( -(xx(:,:,:,ispec2)-x0)**2/(2.0*sigma_h2) &
                                    -(yy(:,:,:,ispec2)-y0)**2/(2.0*sigma_h2) &
                                    -(zz(:,:,:,ispec2)-z0)**2/(2.0*sigma_v2) ) * factor(:,:,:)

              tk(i,j,k,ispec) = tk(i,j,k,ispec) + sum(exp_val(:,:,:) * kernel(:,:,:,ispec2))
              bk(i,j,k,ispec) = bk(i,j,k,ispec) + sum(exp_val(:,:,:))

            enddo
          enddo
        enddo ! (i,j,k)
      enddo ! (ispec2)
    enddo   ! (ispec)
  enddo     ! islice

  if (myrank == 0) write(*,*) 'Done with integration ...'

  ! compute the smoothed kernel values
  do ispec = 1, nspec(1)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          if (bk(i,j,k,ispec) < 0.01 * norm) then ! check the normalization criterion
            print *, 'Problem here --- ', myrank, ispec, i, j, k, bk(i,j,k,ispec), 0.01 * norm
            call exit_mpi(myrank, 'Error computing Gaussian function on the grid')
          endif

          kernel_smooth(i,j,k,ispec) = tk(i,j,k,ispec)/bk(i,j,k,ispec)

        enddo
      enddo
    enddo
  enddo

  open(11,file=trim(ks_file),status='unknown',form='unformatted')
  ! Note: output the following instead of kernel_smooth(:,:,:,1:nspec(1)) to create files of the same sizes
  write(11) kernel_smooth(:,:,:,:)
  close(11)

  ! the maximum value for the smoothed kernel
  call mpi_reduce(maxval(abs(kernel_smooth(:,:,:,1:nspec(1)))), max_new, 1, &
             CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)

  if (myrank == 0) then
    print *, 'Maximum data value before smoothing = ', max_old
    print *, 'Maximum data value after smoothing = ', max_new
  endif

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program smooth_specfem_function
