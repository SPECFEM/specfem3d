! smooth_sem_globe
!
! this program can be used for smoothing a (summed) event kernel,
! where it smooths files with a given input kernel name:
!
! Usage:
!   ./smooth_sem_globe sigma_h(km) sigma_v(km) kernel_file_name scratch_file_dir scratch_topo_dir
!   e.g.
!   ./smooth_sem_globe 160 10 bulk_c_kernel OUTPUT_SUM/ topo/
!
! where:
!   sigma_h                - Gaussian width for horizontal smoothing (in km)
!   sigma_v                - Gaussian width for vertical smoothing (in km)
!   kernel_file_name  - takes file with this kernel name,
!                                     e.g. "bulk_c_kernel"
!   scratch_file_dir     - directory containing kernel files,
!                                      e.g. proc***_reg1_bulk_c_kernel.bin
!   topo_dir                - directory containing mesh topo files:
!                                       proc***_solver_data1.bin, proc***_solver_data2.bin
! outputs:
!    puts the resulting, smoothed kernel files into the same directory as scratch_file_dir/
!    with a file ending "proc***_kernel_smooth.bin"

program smooth_sem_globe

! this is the embarassingly-parallel program that smooths any specfem function (primarily the kernels)
! that has the dimension of (NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT)
!
! notice that it uses the constants_globe.h and precision_globe.h files
! from the original SPECFEM3D_GLOBE package, and the
! values_from_mesher_globe.h file from the output of the mesher (or create_header_file),
! therefore, you need to compile it for your specific case
!
! NOTE:  smoothing can be different in radial & horizontal directions; mesh is in spherical geometry.
!              algorithm uses vector components in radial/horizontal direction

  implicit none
  include 'mpif.h'
  include '../../SHARE_FILES/HEADER_FILES/constants.h'
  include '../../SHARE_FILES/HEADER_FILES/precision.h'
  include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h'

! ======================================================
  ! USER PARAMETERS

  ! taken from values_from_mesher.h:
  !   average size of a spectral element in km = ...
  !   e.g. nproc 12x12, nex 192: element_size = 52.122262
!  real(kind=CUSTOM_REAL),parameter:: element_size = 52.12262
  real(kind=CUSTOM_REAL),parameter:: element_size = 41.69810

! ======================================================

  !takes region 1 kernels
  integer, parameter :: NSPEC_MAX = NSPEC_CRUST_MANTLE_ADJOINT
  integer, parameter :: NGLOB_MAX = NGLOB_CRUST_MANTLE

  ! only include the neighboring 3 x 3 slices
  integer, parameter :: NSLICES = 3
  integer ,parameter :: NSLICES2 = NSLICES * NSLICES

  character(len=256) :: s_sigma_h, s_sigma_v
  character(len=256) :: kernel_file_name, scratch_topo_dir, scratch_file_dir
  integer :: sizeprocs,ier,myrank,ichunk, ixi, ieta, iglob
  integer :: islice(NSLICES2), islice0(NSLICES2), nums

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3

  real(kind=CUSTOM_REAL) :: x0, y0, z0, norm, norm_h, norm_v, element_size_m, max_old, max_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor, exp_val

  character(len=256) ::  ks_file, reg_name
  character(len=256), dimension(NSLICES2) :: x_file, y_file, z_file, j_file, k_file, i_file
  character(len=256), dimension(NSLICES2) :: solver1_file,solver2_file

  logical :: global_code

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: kernel, kernel_smooth
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: tk, bk, jacobian, xl, yl, zl, xx, yy, zz

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_MAX) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NGLOB_MAX) :: x, y, z
  real(kind=CUSTOM_REAL), dimension(NSPEC_MAX) :: cx0, cy0, cz0, cx, cy, cz

  integer :: i,ii,j,jj,k,kk,ispec,iproc,ispec2,nspec(NSLICES2),nglob(NSLICES2)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll

  ! array with all the weights in the cube
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: r1,theta1

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if (myrank == 0) print *,"smooth:"
  call mpi_barrier(MPI_COMM_WORLD,ier)

  ! arguments
  call get_command_argument(1,s_sigma_h)
  call get_command_argument(2,s_sigma_v)
  call get_command_argument(3,kernel_file_name)
  call get_command_argument(4,scratch_file_dir)
  call get_command_argument(5,scratch_topo_dir)

  if ( trim(s_sigma_h) == '' .or. trim(s_sigma_v) == '' &
    .or. trim(kernel_file_name) == '' &
    .or. trim(scratch_file_dir) == '' &
    .or. trim(scratch_topo_dir) == '') then
    call exit_MPI(myrank,'Usage: smooth_sem_globe sigma_h(km) sigma_v(km) kernel_file_name scratch_file_dir scratch_topo_dir')
  endif

  ! read in parameter information
  read(s_sigma_h,*) sigma_h
  read(s_sigma_v,*) sigma_v

  ! checks if basin code or global code: global code uses nchunks /= 0
  if (NCHUNKS_VAL == 0) then
    global_code = .false.
    call exit_mpi(myrank,'Error nchunks')
  else
    global_code = .true.
    reg_name='_reg1_'
  endif
  if (sizeprocs /= NPROC_XI_VAL*NPROC_ETA_VAL*NCHUNKS_VAL) call exit_mpi(myrank,'Error total number of slices')

  ! user output
  if (myrank == 0) then
    print *,"defaults:"
    print *,"  NPROC_XI , NPROC_ETA: ",NPROC_XI_VAL,NPROC_ETA_VAL
    print *,"  NCHUNKS                         : ",NCHUNKS_VAL
    print *,"  element size on surface(km): ",element_size
    print *,"  smoothing sigma_h , sigma_v: ",sigma_h,sigma_v
  endif
  ! synchronizes
  call mpi_barrier(MPI_COMM_WORLD,ier)

  ! initializes lengths
  element_size_m = element_size * 1000  ! e.g. 9 km on the surface, 36 km at CMB
  if (global_code)  element_size_m = element_size_m/R_EARTH

  sigma_h = sigma_h * 1000.0 ! m
  if (global_code) sigma_h = sigma_h / R_EARTH ! scale
  sigma_v = sigma_v * 1000.0 ! m
  if (global_code) sigma_v = sigma_v / R_EARTH ! scale

  sigma_h2 = 2.0 * sigma_h ** 2  ! factor two for Gaussian distribution with standard variance sigma
  sigma_v2 = 2.0 * sigma_v ** 2

  ! search radius
  sigma_h3 = 3.0  * sigma_h + element_size_m
  sigma_v3 = 3.0  * sigma_v + element_size_m

  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )

! note: smoothing is using a Gaussian (ellipsoid for sigma_h /= sigma_v),
!          but in spherical coordinates, we use horizontal distance as epicentral distance
!          and vertical distance as radial distance?

! not squared since epicentral distance is taken? values from bk seem to be closer to squared ones...
  !norm_h = sqrt(2.0*PI) * sigma_h
  norm_h = 2.0*PI*sigma_h**2
  norm_v = sqrt(2.0*PI) * sigma_v
  norm   = norm_h * norm_v
  !norm = (sqrt(2.0*PI) * sigma) ** 3 ! for sigma_h = sigma_v = sigma


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

  ! ---- figure out the neighboring 8 or 7 slices: (ichunk,ixi,ieta) index start at 0------
  ichunk = myrank / (NPROC_XI_VAL * NPROC_ETA_VAL)
  ieta = (myrank - ichunk * NPROC_XI_VAL * NPROC_ETA_VAL) / NPROC_XI_VAL
  ixi = myrank - ichunk * NPROC_XI_VAL * NPROC_ETA_VAL - ieta * NPROC_XI_VAL

  ! get the neighboring slices:
  call get_all_eight_slices(ichunk,ixi,ieta, &
             islice0(1),islice0(2),islice0(3),islice0(4),islice0(5),islice0(6),islice0(7),islice0(8), &
             NPROC_XI_VAL,NPROC_ETA_VAL)

  ! remove the repeated slices (only 8 for corner slices in global case)
  islice(1) = myrank; j = 1
  do i = 1, 8
    if (.not. any(islice(1:i) == islice0(i)) .and. islice0(i) < sizeprocs) then
      j = j + 1
      islice(j) = islice0(i)
    endif
  enddo
  nums = j

  if ( myrank == 0 ) then
    print *,'slices:',nums
    print *,'  ',islice(:)
    print *
  endif

  ! read in the topology files of the current and neighboring slices
  do i = 1, nums
    write(k_file(i),'(a,i6.6,a)') &
      trim(scratch_file_dir)//'/proc',islice(i),trim(reg_name)//trim(kernel_file_name)//'.bin'

    write(solver1_file(i),'(a,i6.6,a)') &
      trim(scratch_topo_dir)//'/proc',islice(i),trim(reg_name)//'solver_data_1.bin'
    write(solver2_file(i),'(a,i6.6,a)') &
      trim(scratch_topo_dir)//'/proc',islice(i),trim(reg_name)//'solver_data_2.bin'

    nspec(i) = NSPEC_MAX
    nglob(i) = NGLOB_MAX
  enddo

  ! point locations
  open(11,file=solver2_file(1),status='old',form='unformatted',iostat=ier)
  if ( ier /= 0 ) call exit_mpi(myrank,'error opening solver2 file')

  read(11) x(1:nglob(1))
  read(11) y(1:nglob(1))
  read(11) z(1:nglob(1))
  read(11) ibool(:,:,:,1:nspec(1))
  close(11)

  ! jacobian
  open(11,file=solver1_file(1),status='old',form='unformatted',iostat=ier)
  if ( ier /= 0 ) call exit_mpi(myrank,'error opening solver1 file')

  read(11) xix
  read(11) xiy
  read(11) xiz
  read(11) etax
  read(11) etay
  read(11) etaz
  read(11) gammax
  read(11) gammay
  read(11) gammaz
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

          ! build jacobian
          ! get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          ! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          jacobian(i,j,k,ispec) = jacobianl



        enddo
      enddo
    enddo
    cx0(ispec) = (xl(1,1,1,ispec) + xl(NGLLX,NGLLY,NGLLZ,ispec))/2.0
    cy0(ispec) = (yl(1,1,1,ispec) + yl(NGLLX,NGLLY,NGLLZ,ispec))/2.0
    cz0(ispec) = (zl(1,1,1,ispec) + zl(NGLLX,NGLLY,NGLLZ,ispec))/2.0
  enddo

  if (myrank == 0) write(*,*) 'start looping over elements and points for smoothing ...'

  ! smoothed kernel file name
  write(ks_file,'(a,i6.6,a)') trim(scratch_file_dir)//'/proc',myrank, &
                          trim(reg_name)//trim(kernel_file_name)//'_smooth.bin'


  tk = 0.0
  bk = 0.0
  kernel_smooth=0.0

  ! loop over all the slices
  do iproc = 1, nums

    ! read in the topology, kernel files, calculate center of elements
    ! point locations
    ! given in Cartesian coordinates
    open(11,file=solver2_file(iproc),status='old',form='unformatted',iostat=ier)
    if ( ier /= 0 ) call exit_mpi(myrank,'error opening slices: solver2 file')

    read(11) x(1:nglob(iproc))
    read(11) y(1:nglob(iproc))
    read(11) z(1:nglob(iproc))
    read(11) ibool(:,:,:,1:nspec(iproc))
    close(11)

    open(11,file=solver1_file(iproc),status='old',form='unformatted',iostat=ier)
    if ( ier /= 0 ) call exit_mpi(myrank,'error opening slices: solver1 file')

    read(11) xix
    read(11) xiy
    read(11) xiz
    read(11) etax
    read(11) etay
    read(11) etaz
    read(11) gammax
    read(11) gammay
    read(11) gammaz
    close(11)

    ! get the location of the center of the elements
    do ispec = 1, nspec(iproc)
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! build jacobian
            ! get derivatives of ux, uy and uz with respect to x, y and z
            xixl = xix(i,j,k,ispec)
            xiyl = xiy(i,j,k,ispec)
            xizl = xiz(i,j,k,ispec)
            etaxl = etax(i,j,k,ispec)
            etayl = etay(i,j,k,ispec)
            etazl = etaz(i,j,k,ispec)
            gammaxl = gammax(i,j,k,ispec)
            gammayl = gammay(i,j,k,ispec)
            gammazl = gammaz(i,j,k,ispec)
            ! compute the jacobian
            jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                          - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                          + xizl*(etaxl*gammayl-etayl*gammaxl))
            jacobian(i,j,k,ispec) = jacobianl
          enddo
        enddo
      enddo
    enddo

    ! kernel file
    open(11,file=k_file(iproc),status='old',form='unformatted',iostat=ier)
    if ( ier /= 0 ) call exit_mpi(myrank,'error opening kernel file')

    read(11) kernel(:,:,:,1:nspec(iproc))
    close(11)


    ! get the global maximum value of the original kernel file
    if (iproc == 1) then
      call mpi_reduce(maxval(abs(kernel(:,:,:,1:nspec(iproc)))), max_old, 1, &
                 CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
    endif

    ! calculate element center location
    do ispec2 = 1, nspec(iproc)
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
      cx(ispec2) = (xx(1,1,1,ispec2) + xx(NGLLX,NGLLZ,NGLLY,ispec2))/2.0
      cy(ispec2) = (yy(1,1,1,ispec2) + yy(NGLLX,NGLLZ,NGLLY,ispec2))/2.0
      cz(ispec2) = (zz(1,1,1,ispec2) + zz(NGLLX,NGLLZ,NGLLY,ispec2))/2.0
    enddo

    ! loop over elements to be smoothed in the current slice
    do ispec = 1, nspec(1)

      ! --- only double loop over the elements in the search radius ---
      do ispec2 = 1, nspec(iproc)

        ! calculates horizontal and vertical distance between two element centers

        ! vector approximation
        call get_distance_vec(dist_h,dist_v,cx0(ispec),cy0(ispec),cz0(ispec), &
                          cx(ispec2),cy(ispec2),cz(ispec2))

        ! note: distances and sigmah, sigmav are normalized by R_EARTH

        ! checks distance between centers of elements
        if ( dist_h > sigma_h3 .or. abs(dist_v) > sigma_v3 ) cycle

        ! integration factors:
        ! uses volume assigned to GLL points
        factor(:,:,:) = jacobian(:,:,:,ispec2) * wgll_cube(:,:,:)
        ! no volume
        !factor(:,:,:) = 1.0_CUSTOM_REAL

        ! loop over GLL points of the elements in current slice (ispec)
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX

              ! reference location
              ! current point (i,j,k,ispec) location, Cartesian coordinates
              x0 = xl(i,j,k,ispec)
              y0 = yl(i,j,k,ispec)
              z0 = zl(i,j,k,ispec)

              ! calculate weights based on Gaussian smoothing
              call smoothing_weights_vec(x0,y0,z0,ispec2,sigma_h2,sigma_v2,exp_val, &
                      xx(:,:,:,ispec2),yy(:,:,:,ispec2),zz(:,:,:,ispec2))

              ! adds GLL integration weights
              exp_val(:,:,:) = exp_val(:,:,:) * factor(:,:,:)

              ! adds contribution of element ispec2 to smoothed kernel values
              tk(i,j,k,ispec) = tk(i,j,k,ispec) + sum(exp_val(:,:,:) * kernel(:,:,:,ispec2))

              ! normalization, integrated values of Gaussian smoothing function
              bk(i,j,k,ispec) = bk(i,j,k,ispec) + sum(exp_val(:,:,:))

              ! checks number
              !if ( isNaN(tk(i,j,k,ispec)) ) then
              !  print *,'error tk NaN: ',tk(i,j,k,ispec)
              !  print *,'rank:',myrank
              !  print *,'i,j,k,ispec:',i,j,k,ispec
              !  print *,'tk: ',tk(i,j,k,ispec),'bk:',bk(i,j,k,ispec)
              !  print *,'sum exp_val: ',sum(exp_val(:,:,:)),'sum factor:',sum(factor(:,:,:))
              !  print *,'sum kernel:',sum(kernel(:,:,:,ispec2))
              !  call exit_MPI(myrank,'error NaN')
              !endif

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

          ! checks the normalization criterion
          ! e.g. sigma_h 160km, sigma_v 40km:
          !     norm (not squared sigma_h ) ~ 0.001
          !     norm ( squared sigma_h) ~ 6.23 * e-5
          if (abs(bk(i,j,k,ispec) - norm) > 1.e-4 ) then
            print *, 'Problem norm here --- ', myrank, ispec, i, j, k, bk(i,j,k,ispec), norm
            !call exit_mpi(myrank, 'Error computing Gaussian function on the grid')
          endif

          ! normalizes smoothed kernel values by integral value of Gaussian weighting
          kernel_smooth(i,j,k,ispec) = tk(i,j,k,ispec) / bk(i,j,k,ispec)


          ! checks number
          if ( isNaN(kernel_smooth(i,j,k,ispec)) ) then
            print *,'error kernel_smooth NaN: ',kernel_smooth(i,j,k,ispec)
            print *,'rank:',myrank
            print *,'i,j,k,ispec:',i,j,k,ispec
            print *,'tk: ',tk(i,j,k,ispec),'bk:',bk(i,j,k,ispec)
            call exit_MPI(myrank,'error NaN')
          endif

        enddo
      enddo
    enddo
  enddo
  if (myrank == 0) write(*,*) '  norm: ',norm

  ! file output
  open(11,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  if ( ier /= 0 ) call exit_mpi(myrank,'error opening smoothed kernel file')

  ! Note: output the following instead of kernel_smooth(:,:,:,1:nspec(1)) to create files of the same sizes
  write(11) kernel_smooth(:,:,:,:)
  close(11)

  if (myrank == 0) print *,'  written:',trim(ks_file)



  ! the maximum value for the smoothed kernel
  call mpi_reduce(maxval(abs(kernel_smooth(:,:,:,1:nspec(1)))), max_new, 1, &
             CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)

  if (myrank == 0) then
    print *
    print *, 'Maximum data value before smoothing = ', max_old
    print *, 'Maximum data value after smoothing = ', max_new
    print *
  endif

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program smooth_sem_globe

!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,y0,z0,ispec2,sigma_h2,sigma_v2,exp_val, &
                              xx_elem,yy_elem,zz_elem)

  implicit none
  include "../../SHARE_FILES/HEADER_FILES/constants.h"

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: xx_elem, yy_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,sigma_h2,sigma_v2
  integer,intent(in) :: ispec2

  ! local parameters
  integer :: ii,jj,kk
  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  !real(kind=CUSTOM_REAL) :: r0,r1,theta1

  ! >>>>>
  ! uniform sigma
  !exp_val(:,:,:) = exp( -((xx(:,:,:,ispec2)-x0)**2+(yy(:,:,:,ispec2)-y0)**2 &
  !          +(zz(:,:,:,ispec2)-z0)**2 )/(2*sigma2) )*factor(:,:,:)

  ! from basin code smoothing:
  ! Gaussian function
  !exp_val(:,:,:) = exp( -(xx(:,:,:,ispec2)-x0)**2/(sigma_h2) &
  !                      -(yy(:,:,:,ispec2)-y0)**2/(sigma_h2) &
  !                      -(zz(:,:,:,ispec2)-z0)**2/(sigma_v2) ) * factor(:,:,:)
  ! >>>>>

  do kk = 1, NGLLZ
    do jj = 1, NGLLY
      do ii = 1, NGLLX
        ! point in second slice

        ! vector approximation:
        call get_distance_vec(dist_h,dist_v,x0,y0,z0, &
            xx_elem(ii,jj,kk),yy_elem(ii,jj,kk),zz_elem(ii,jj,kk))

        ! Gaussian function
        exp_val(ii,jj,kk) = exp( - (dist_h*dist_h)/sigma_h2 &
                                  - (dist_v*dist_v)/sigma_v2 )    ! * factor(ii,jj,kk)


        ! checks number
        !if ( isNaN(exp_val(ii,jj,kk)) ) then
        !  print *,'error exp_val NaN: ',exp_val(ii,jj,kk)
        !  print *,'i,j,k:',ii,jj,kk
        !  print *,'dist_h: ',dist_h,'dist_v:',dist_v
        !  print *,'sigma_h2:',sigma_h2,'sigma_v2:',sigma_v2
        !  call exit_MPI(myrank,'error NaN')
        !endif

      enddo
    enddo
  enddo

  end subroutine smoothing_weights_vec


!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_vec(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

! returns vector lengths as distances in radial and horizontal direction

  implicit none
  include "../../SHARE_FILES/HEADER_FILES/constants.h"

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,x1,y1,z1

  ! local parameters
  real(kind=CUSTOM_REAL) :: r0,r1
  real(kind=CUSTOM_REAL) :: theta,ratio
  !real(kind=CUSTOM_REAL) :: vx,vy,vz,alpha

  ! vertical distance
  r0 = sqrt( x0*x0 + y0*y0 + z0*z0 ) ! length of first position vector
  r1 = sqrt( x1*x1 + y1*y1 + z1*z1 )
  dist_v = r1 - r0
  ! only for flat earth with z in depth: dist_v = sqrt( (cz(ispec2)-cz0(ispec))** 2)

  ! epicentral distance
  ! (accounting for spherical curvature)
  ! calculates distance of circular segment
  ! angle between r0 and r1 in radian
  ! given by dot-product of two vectors
  ratio = (x0*x1 + y0*y1 + z0*z1)/(r0 * r1)

  ! checks boundaries of ratio (due to numerical inaccuracies)
  if ( ratio > 1.0_CUSTOM_REAL ) ratio = 1.0_CUSTOM_REAL
  if ( ratio < -1.0_CUSTOM_REAL ) ratio = -1.0_CUSTOM_REAL

  theta = acos( ratio )

  ! segment length at heigth of r1
  dist_h = r1 * theta

  ! vector approximation (fast computation): neglects curvature
  ! horizontal distance
  ! length of vector from point 0 to point 1
  ! assuming small earth curvature  (since only for neighboring elements)

  ! scales r0 to have same length as r1
  !alpha = r1 / r0
  !vx = alpha * x0
  !vy = alpha * y0
  !vz = alpha * z0

  ! vector in horizontal between new r0 and r1
  !vx = x1 - vx
  !vy = y1 - vy
  !vz = z1 - vz

  ! distance is vector length
  !dist_h = sqrt( vx*vx + vy*vy + vz*vz )

  end subroutine get_distance_vec

