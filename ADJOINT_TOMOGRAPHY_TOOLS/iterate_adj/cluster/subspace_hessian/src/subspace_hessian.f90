program subspace_hessian

  implicit none
  include 'mpif.h'
  include 'constants.h'
  include 'precision.h'
  include 'values_from_mesher.h'

! ======================================================

!  integer, parameter :: NSLICES=168
  integer, parameter :: NSPEC=NSPEC_AB
  integer, parameter :: NGLOB=NGLOB_AB

  character(len=150) :: kernel_file_list, kernel_list(1000), sline, k_file, kernel_name, kernel_name1, kernel_name2, kdir
  character(len=150) :: smodel
  logical :: global_code
  integer :: nchunks
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel, kernel_j, dV
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: kdot,kdot_total
  real(kind=CUSTOM_REAL) :: dV_sum, dV_sum_total, sigma_structure, dVfac, coverage_structure, dot_prod_fac, kfac

  integer :: iker, jker, nker, ivar, nops, myrank, sizeprocs,  ier, ismooth
  integer :: i, j, k, ispec, ios
  !integer, dimension(:), allocatable :: nwin_vec

  character(len=150) :: scratch_topo_dir, reg_name
  character(len=150) :: j_file, filename
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: jacobian

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

  !-----------------------------------------------
  ! input parameters

  kernel_file_list = 'kernels_run'
  kdir = 'INPUT_KERNELS/'
  sigma_structure = 0.1
  coverage_structure = 1.0
  dot_prod_fac = 1.0e10         ! should be on the order of 1/norm(kernel)

  !ismooth = 1                   ! smoothed event kernels (1) or not (0)
  scratch_topo_dir = 'topo'
  nchunks = 0

  ! read in the tag designating the model number
  open(unit=11,file='INPUT/smodel',status='old',iostat=ios)
  read(11,*) smodel
  close(11)

  ! read in the tag designating whether the kernels are smoothed (1) or not (0)
  open(unit=12,file='INPUT/ismooth',status='old',iostat=ios)
  read(12,*) ismooth
  close(12)

  if(ismooth == 1) then
     kernel_name1 = 'mu_kernel_smooth'
     kernel_name2 = 'kappa_kernel_smooth'
  else
     kernel_name1 = 'mu_kernel'
     kernel_name2 = 'kappa_kernel'
  endif

  ! read in list of event kernels
  nker=0
  open(unit = 20, file = trim(kernel_file_list), status = 'old',iostat = ios)
  if (ios /= 0) then
    print *,'Error opening ',trim(kernel_file_list)
    stop
  endif
  do while (1 == 1)
    read(20,'(a)',iostat=ios) sline
    if (ios /= 0) exit
    nker=nker+1
    kernel_list(nker) = sline
  enddo
  close(20)

  !-----------------------------------------------
  ! read in number of windows picked per event
  ! --> this is a factor used in the data covariance matrix

!!$  allocate(nwin_vec(nsrc))
!!$  nwin_vec = 0
!!$  if (myrank == 0) write(*,*) 'reading in number of windows per event'
!!$  open(unit=20,file=trim(win_file),status='old',iostat=ios)
!!$  do i=1,nsrc
!!$     read(20,*) nwin_vec(i)
!!$     if(myrank==0) write(*,*) nwin_vec(i)
!!$  enddo
!!$  close(20)

  !-----------------------------------------------

  if (nchunks == 0) then
    global_code = .false.
    reg_name='_'
    nchunks = 1
  else
    global_code = .true.
    reg_name='_reg1_'
  endif

  !-----------------------------------------------
  ! read in the topology files

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

  if (myrank == 0) write(*,*) 'reading in the topology files...'

  ! read in the topology files
  ! BUG FIX 10-14-08: i-1 --> myrank
  write(j_file,'(a,i6.6,a)') trim(scratch_topo_dir) // '/proc',myrank,trim(reg_name) // 'jacobian.bin'
  open(11,file=j_file,status='old',form='unformatted')
  read(11) jacobian(:,:,:,1:NSPEC)
  close(11)

  ! volumetric integration array for local GLL points
  do ispec = 1,NSPEC  ! local
     dV(:,:,:,ispec) =  wgll_cube(:,:,:) * jacobian(:,:,:,ispec)
  enddo
  dV_sum = sum( dV(:,:,:,:) )

  ! GOAL: Output each dV_sum to a file (presently, it goes to the .o file) -- 168 numbers in total
  ! write(*,*) 'myrank = ',myrank,' --> sum of dV array is :', dV_sum

  ! sum the volume pieces from all slices
  ! 1 --> 1 real number
  ! 0 --> all info to the 0 slice
  call mpi_reduce(dV_sum,dV_sum_total,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  dVfac = dV_sum_total * coverage_structure
  !dVfac = 1.0

  ! output the total volume and sigma_structure to a file
  if (myrank == 0) then

     open(18,file='coverage_structure',status='unknown')
     write(18,*) coverage_structure
     close(18)

     open(19,file='sigma_structure',status='unknown')
     write(19,*) sigma_structure
     close(19)

     open(20,file='volume_total',status='unknown')
     write(20,*) dV_sum_total
     close(20)

     open(20,file='dVfac',status='unknown')
     write(20,*) dVfac
     close(20)

     write(*,*) 'sigma_structure :', sigma_structure
     write(*,*) 'structure coverage :', sigma_structure
     write(*,*) 'sum of dV array is :', dV_sum_total
     write(*,*) 'dVfac :', dVfac
  endif

  !-----------------------------------------------
  ! perform dot-product operations on event kernels, then output a vector of values
  ! OPERATION: H = G Cm G^T
  ! NOTES:
  ! (1) Model covariance matrix is assumed to be diagonal: Cm = sigma_i^2 V / V_i
  ! (2) We assume sigma_i is constant everywhere.
  ! (3) The factor of 4 = 2*2 is for K_c = 2 K_kappa and K_beta = 2 K_mu
  ! (4) g Cm g = (K_i V_i) (V/V_i) (K_i V_i) = K_i K_i V_i V
  ! (5) The WEIGHTS from the data are NOT included here (see subspace_specfem.m)
  ! (6) The MINUS SIGN for each kernel will cancel in the dot product operation.
  ! (7) The volume integration factor is dVfac = dV_sum_total * coverage_structure.
  !     If we expected to recover the full model, then this would give a norm of 1.
  !     coverage_structure expresses the fraction of the full volume that you expect
  !     to recover.

  nops = nker*(nker+1)/2     ! number of operations
  if (myrank == 0) write(*,*) 'number of kernels, number of operations :',nker,nops
  allocate(kdot(nops),kdot_total(nops))

  ! GOAL 3: Output a file with nops rows, each row a dot-product between two event kernels

  !-------------

  do ivar = 1,2

     if (myrank == 0) write(*,*) 'performing the dot-product operations...'

     if (ivar == 1) kernel_name = kernel_name1
     if (ivar == 2) kernel_name = kernel_name2

     k = 0
     kdot(:) = 0.
     kdot_total(:) = 0.

     do iker = 1, nker

        ! read in event kernel i
        write(k_file,'(a,i6.6,a)') trim(kdir) // trim(kernel_list(iker)) // '/proc',myrank,'_' // trim(kernel_name) // '.bin'
        open(12,file=trim(k_file),status='old',form='unformatted')
        read(12) kernel(:,:,:,1:NSPEC)
        close(12)

        if (myrank == 0) write(*,'(3i8)') iker, k

        do jker = iker, nker
           k = k+1

           ! read in event kernel j
           write(k_file,'(a,i6.6,a)') trim(kdir) // trim(kernel_list(jker)) // '/proc',myrank,'_' // trim(kernel_name) // '.bin'
           open(12,file=trim(k_file),status='old',form='unformatted')
           read(12) kernel_j(:,:,:,1:NSPEC)
           close(12)

           ! dot product (see notes above)
           ! NOTE: The numbers are slightly different if dV_sum_total is taken OUTSIDE;
           !       this is probably due to single precision computations or storage.
           !       This is why we use the factor dot_prod_fac .

           !kdot(k) = 4.0 * sigma_structure**2 / dot_prod_fac**2 * sum( (kernel(:,:,:,:)*dot_prod_fac) &
           !     * (kernel_j(:,:,:,:)*dot_prod_fac) * (dV(:,:,:,:)*dVfac) )      ! version 1

           kdot(k) = sum( (kernel(:,:,:,:)*dot_prod_fac) * (kernel_j(:,:,:,:)*dot_prod_fac) * dV(:,:,:,:) )

           !kdot(k) = 4.0 * dVfac * sum( kernel(:,:,:,:) * kernel_j(:,:,:,:) * dV(:,:,:,:) )   ! version 2
           !kdot(k) = sum( kernel(:,:,:,:) * kernel_j(:,:,:,:) * dV(:,:,:,:) )                 ! version 3

           !if (myrank == 0) write(*,'(3i8)') k, iker, jker

        enddo
     enddo

     call mpi_barrier(MPI_COMM_WORLD,ier)
     call mpi_reduce(kdot,kdot_total,k,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

     ! final scaling factor
     kfac = 4.0 * sigma_structure**2 / dot_prod_fac**2 * dVfac
     kdot_total(:) = kdot_total(:) * kfac

     if (myrank == 0) then
        write(*,*) 'writing out files...'

        write(filename,'(a,i3.3)') 'hessian_index_' // trim(smodel) // '_' // trim(kernel_name) // '_', nker
        open(20,file=trim(filename),status='unknown')
        k = 0
        do iker = 1, nker
           do jker = iker, nker
              k = k+1
              write(20,'(3i8,1e24.12)') k,iker,jker,kdot_total(k)
           enddo
        enddo
        close(20)

        !write(filename,'(a,i3.3)') 'hessian_' // trim(kernel_name) // '_', nker
        !open(19,file=trim(filename),status='unknown')
        !write(19,'(i8,1e24.12)') (kdot_total(k),k=1,nops)
        !close(19)

        write(filename,'(a)') 'kfac_' // trim(smodel) // '_' // trim(kernel_name)
        open(18,file=trim(filename),status='unknown')
        write(18,'(1e24.12)') kfac
        close(18)

     endif

  enddo  ! ivar

  !-----------------------------------------------

  deallocate(kdot,kdot_total)

  !-----------------------------------------------

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program subspace_hessian


