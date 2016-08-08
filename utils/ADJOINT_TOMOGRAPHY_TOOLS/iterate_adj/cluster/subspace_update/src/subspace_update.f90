program subspace_update

  implicit none
  include 'mpif.h'
  include 'constants.h'
  include 'precision.h'
  include 'values_from_mesher.h'

! ======================================================

  integer, parameter :: NSPEC=NSPEC_AB
  integer, parameter :: NGLOB=NGLOB_AB

  character(len=150) :: kernel_file_list, kernel_list(1000), sline, k_file, kernel_name, dm_file, filename
  character(len=150) :: dcov_tag, idir, kdir, odir, beta_dir, bulk_dir
  logical :: global_code
  integer :: nchunks
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel, dm
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mu_beta, mu_bulk, dnorm, dcov_fac
  real(kind=CUSTOM_REAL) :: lambdas(1000), lam
  real(kind=CUSTOM_REAL) :: dVfac, kfac, sigma_structure
  integer :: isrc, nsrc, myrank, sizeprocs,  ier
  integer :: i, j, k, p, ispec, ios

  integer :: pmax, ptemp, npmax, pmax_list(1000)

 ! ============ program starts here =====================
 ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  !-----------------------------------------------
  ! input parameters

  kernel_file_list = 'kernels_run'
  kdir = 'INPUT_KERNELS/'                           ! input kernels directory
  idir = 'INPUT/mu_run/'                            ! input directory
  odir = 'OUTPUT_MODEL/'                            ! output directory

  ! read in the tag designating the factors used in the data covariance (event or window)
  open(unit=18,file=trim(idir)//'dcov_tag',status='old',iostat=ios)
  read(18,*) dcov_tag
  close(18)

  ! read in the tag designating the kernels to use
  open(unit=20,file=trim(idir)//'ker_tag',status='old',iostat=ios)
  read(20,*) kernel_name
  close(20)

  if (myrank == 0) then
     write(*,*) 'dcov_tag : ', trim(dcov_tag)
     write(*,*) 'kernel_name : ', trim(kernel_name)

     write(*,*) 'kdir : ', trim(kdir)
     write(*,*) 'idir : ', trim(idir)
     write(*,*) 'odir : ', trim(odir)
  endif

!!$  nchunks = 0
!!$  scratch_topo_dir = 'topo'

!!$  if (nchunks == 0) then
!!$    global_code = .false.
!!$    reg_name='_'
!!$    nchunks = 1
!!$  else
!!$    global_code = .true.
!!$    reg_name='_reg1_'
!!$  endif

  !-----------------------------------------------
  ! READ IN FILES

  ! read in list of event kernels
  nsrc=0
  open(unit=19,file=trim(kernel_file_list),status='old',iostat=ios)
  if (ios /= 0) then
    print *,'Error opening ',trim(kernel_file_list)
    stop
  endif
  do while (1 == 1)
    read(19,'(a)',iostat=ios) sline
    if (ios /= 0) exit
    nsrc=nsrc+1
    kernel_list(nsrc) = sline
  enddo
  close(19)

  if (myrank == 0) then
     write(*,*) 'total number of events (and kernels) : ',nsrc
     write(*,*) (kernel_list(isrc),isrc=1,nsrc)
  endif

  ! allocate arrays
  allocate(mu_beta(nsrc),mu_bulk(nsrc),dnorm(nsrc),dcov_fac(nsrc))

  ! read in the data-norm vector for subspace method
  if (myrank == 0) write(*,*) 'reading in data_norm'
  open(unit=22,file=trim(idir)//'data_norm',status='old',iostat=ios)
  do i=1,nsrc
     read(22,*) dnorm(i)
     if (myrank == 0) write(*,*) dnorm(i)
  enddo
  close(22)

  ! read in the number of window picks per event
  if (myrank == 0) write(*,*) 'reading in window tallies'
  open(unit=23,file=trim(idir)//'dcov_fac',status='old',iostat=ios)
  do i=1,nsrc
     read(23,*) dcov_fac(i)
     if (myrank == 0) write(*,*) dcov_fac(i)
  enddo
  close(23)

  ! read in the model covariance normalization term
  if (myrank == 0) write(*,*) 'reading in model covariance normalization term'
  open(unit=24,file='INPUT/dVfac',status='old',iostat=ios)
     read(24,*) dVfac
  close(24)
  if (myrank == 0) write(*,*) dVfac

  ! read in the sigma value (fixed) used for the model covariance
  if (myrank == 0) write(*,*) 'reading in sigma_structure'
  open(unit=25,file='INPUT/sigma_structure',status='old',iostat=ios)
     read(25,*) sigma_structure
  close(25)
  if (myrank == 0) write(*,*) sigma_structure

!---------------------------------------------------

  ! read in the pmax values
  ! pmax : max singular value for TSVD construction of subspace Hessian
  npmax=0
  open(unit=19,file=trim(idir)//'pmax',status='old',iostat=ios)
  if (ios /= 0) then
    print *,'Error opening ',trim(filename)
    stop
  endif
  do while (1 == 1)
    read(19,*,iostat=ios) ptemp
    if (ios /= 0) exit
    npmax=npmax+1
    pmax_list(npmax) = ptemp
  enddo
  close(19)

  if (myrank == 0) then
     write(*,*) 'total number of pmax update models : ',npmax
     write(*,*) (pmax_list(i),i=1,npmax)
  endif

  ! read in the damping values
  if (myrank == 0) write(*,*) 'reading in lambda values'
  open(unit=26,file=trim(idir)//'lams',status='old',iostat=ios)
  read(26,*) (lambdas(p),p=1,npmax)
  close(26)

  if (myrank == 0) then
     write(*,*) (lambdas(p),p=1,npmax)
  endif

!---------------------------------------------------

  do p = 1,npmax

     pmax = pmax_list(p)
     lam  = lambdas(p)
     if (myrank == 0) write(*,*) 'pmax = ',pmax,' -- ',p,' out of ',npmax
     if (myrank == 0) write(*,*) 'lambda = ', lam

     ! read in mu vector for subspace method
     if (myrank == 0) write(*,*) 'reading in mu vectors'
     write(filename,'(a,i3.3)') trim(idir) // 'mu_p', pmax
     open(unit=20,file=trim(filename),status='old',iostat=ios)
     do i=1,nsrc
        read(20,*) mu_beta(i)
        if (myrank == 0) write(*,*) mu_beta(i)
     enddo
     close(20)

     !-----------------------------------------------
     ! OPERATION: dm = Cm G^T mu
     ! NOTES:
     ! (1) We can omit the volumetric integration V_i, because it appears in the denominator
     !     of the Cm diagonal elements and also in the M x 1 vector (G^T mu) elements.
     !     In this manner, we should think of the operation Cm G^T as being applied FIRST.
     ! (2) The factor of 2 (in kfac) is because K_bulk = 2 K_kappa and K_beta = 2 K_mu .
     ! (3) The factor dcov_fac(isrc) is a factor from the data covariance matrix that was
     !     left out in computing the adjoint sources, but it WAS included in computing
     !     the weight terms dnorm.
     ! (4) The MINUS sign is needed for each gradient (kfac).
     ! (5) To invert H = G^T G, we added lamda^2 to the diagonal. This is like multiplying
     !     every gradient (kernel) by lambda, which we must also do here.

     if (myrank == 0) write(*,*) 'performing the operation G-t mu -- nsrc = ', nsrc

     dm(:,:,:,:) = 0.      ! initialize model update (on this slice)

     do isrc = 1, nsrc
        if (myrank == 0) write(*,*) 'model index ', p, ' -- event index ', isrc

        ! read in event kernel
        write(k_file,'(a,i6.6,a)') trim(kdir) // trim(kernel_list(isrc)) // '/proc',myrank,'_' // trim(kernel_name) // '.bin'
        open(12,file=trim(k_file),status='old',form='unformatted')
        read(12) kernel(:,:,:,1:NSPEC)
        close(12)

        !kfac = -2.0 * dVfac * sigma_structure**2 * mu_beta(isrc) / (dcov_fac(isrc) * dnorm(isrc) )
        kfac = -2.0 * dVfac * sigma_structure**2 * mu_beta(isrc) / dcov_fac(isrc)

        ! BE CAREFUL ABOUT SINGLE-PRECISION ADDITION OF SMALL NUMBERS

        do ispec = 1, NSPEC
           do k = 1, NGLLZ
              do j = 1, NGLLY
                 do i = 1, NGLLX
                    dm(i,j,k,ispec) = dm(i,j,k,ispec) + kfac * kernel(i,j,k,ispec)
                 enddo
              enddo
           enddo
        enddo

     enddo

     !-----------------------------------------------
     ! write out the model update

     if (myrank == 0) write(*,*) 'writing out the model update...'
     write(dm_file,'(a,i6.6,a,i3.3,a)') trim(odir) // 'proc', myrank, '_dm_' // trim(kernel_name) // '_p', pmax, '.bin'
     open(12,file=trim(dm_file),form='unformatted')
     write(12) dm(:,:,:,1:nspec)
     close(12)

  enddo  ! for p

  if (myrank == 0) write(*,*) 'done writing out the model updates'

  !-----------------------------------------------

  deallocate(mu_beta,mu_bulk,dnorm,dcov_fac)

  !-----------------------------------------------

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program subspace_update


