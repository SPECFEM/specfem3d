program cg_test

  implicit none
  include 'mpif.h'
  include 'constants.h'
  include 'precision.h'
  include 'values_from_mesher.h'

! ======================================================

!  integer, parameter :: NSLICES=168
  integer, parameter :: NSPEC=NSPEC_AB
  integer, parameter :: NGLOB=NGLOB_AB

  character(len=150) :: sline, k_file, kernel_name
  logical :: global_code
  integer :: nchunks
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel_mu, kernel_kappa, dV
  real(kind=CUSTOM_REAL) :: nu_step, nu_step_mu, nu_step_kappa, dmisfit, kfac, kfac_mu, kfac_kappa, dmisfit_stop
  real(kind=CUSTOM_REAL) :: gdot_kappa, gdot_kappa_total, gdot_mu, gdot_mu_total, gdot_total
  real(kind=CUSTOM_REAL) :: dVfac, sigma_structure, dcov_fac, dot_prod_fac

  integer :: iker, jker, nker, nops, myrank, sizeprocs,  ier
  integer :: i, j, k, ispec, ios
  integer :: nwin_tot

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

  dot_prod_fac = 1.0e8         ! should be on the order of 1/norm(kernel)
  scratch_topo_dir = 'topo'
  nchunks = 0

  !-----------------------------------------------
  ! READ in constants for model covariance and data covariance

  ! read in the data-norm vector for subspace method
  if (myrank == 0) write(*,*) 'reading in total volume'
  open(unit=24,file='INPUT/volume_total',status='old',iostat=ios)
     read(24,*) dVfac
  close(24)
  if (myrank == 0) write(*,*) dVfac

  ! read in the sigma value (fixed) used for the model covariance
  if (myrank == 0) write(*,*) 'reading in sigma_structure'
  open(unit=25,file='INPUT/sigma_structure',status='old',iostat=ios)
     read(25,*) sigma_structure
  close(25)
  if (myrank == 0) write(*,*) sigma_structure

  ! read in the total number of measurement windows
  if (myrank == 0) write(*,*) 'reading in nwin_tot'
  open(unit=26,file='INPUT/nwin_tot',status='old',iostat=ios)
     read(26,*) nwin_tot
  close(26)
  if (myrank == 0) write(*,*) nwin_tot

  ! read in the misfit value
  if (myrank == 0) write(*,*) 'reading in misfit function value'
  open(unit=27,file='INPUT/dmisfit',status='old',iostat=ios)
     read(27,*) dmisfit
  close(27)
  if (myrank == 0) write(*,*) dmisfit

  !-----------------------------------------------

  !dmisfit_stop = 0.5
  dmisfit_stop = dmisfit / 20.0
  dcov_fac = 1.0

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
  write(j_file,'(a,i6.6,a)') trim(scratch_topo_dir)//'/proc',i-1,trim(reg_name)//'jacobian.bin'
  open(11,file=j_file,status='old',form='unformatted')
  read(11) jacobian(:,:,:,1:NSPEC)
  close(11)

  ! volumetric integration array for local GLL points
  do ispec = 1,NSPEC
     dV(:,:,:,ispec) =  wgll_cube(:,:,:) * jacobian(:,:,:,ispec)
  enddo

  !-----------------------------------------------
  ! perform dot-product operation
  ! OPERATION:
  ! NOTES:

  if (myrank == 0) write(*,*) 'performing the dot-product operation'

  kernel_name = 'mu_kernel_smooth'

  ! read in summed event kernel
  write(k_file,'(a,i6.6,a)') 'inout/proc', myrank, '_'//trim(kernel_name)//'.bin'
  open(12,file=trim(k_file),status='old',form='unformatted')
  read(12) kernel_mu(:,:,:,1:NSPEC)
  close(12)

  ! dot product
  ! NOTE: The numbers are slightly different if dVfac is taken OUTSIDE;
  !       this is probably due to single precision computations or storage.

  gdot_mu = 4.0 * sigma_structure**2 / dot_prod_fac**2 * sum( (kernel_mu(:,:,:,:)*dot_prod_fac) * (kernel_mu(:,:,:,:)*dot_prod_fac) * dV(:,:,:,:) * dVfac / dcov_fac**2 )

  call mpi_barrier(MPI_COMM_WORLD,ier)

  call mpi_reduce(gdot_mu,gdot_mu_total,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  !-----------------------------------------------

  kernel_name = 'kappa_kernel_smooth'
  write(k_file,'(a,i6.6,a)') 'inout/proc', myrank, '_'//trim(kernel_name)//'.bin'
  open(12,file=trim(k_file),status='old',form='unformatted')
  read(12) kernel_kappa(:,:,:,1:NSPEC)
  close(12)

  gdot_kappa = 4.0 * sigma_structure**2 / dot_prod_fac**2 * sum( (kernel_kappa(:,:,:,:)*dot_prod_fac) * (kernel_kappa(:,:,:,:)*dot_prod_fac) * dV(:,:,:,:) * dVfac / dcov_fac**2 )

  !-----------------------------------------------

  call mpi_barrier(MPI_COMM_WORLD,ier)
  call mpi_reduce(gdot_mu,gdot_mu_total,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(gdot_kappa,gdot_kappa_total,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  !-----------------------------------------------
  ! compute the test model for the CG algorithm
  ! NOTES:
  ! (1) The V_i factor appears in both Cm and in g and CANCEL.
  ! (2) sigma_structure**2 * dVfac is for Cm.
  ! (3) The 2.0 is for K_bulk = 2 K_kappa (or K_beta = 2 K_mu ).
  ! (4) The 1/nwin_tot factor was left out of the adjoint source.
  ! (5) 0.5*dmisfit is because dmisfit does not include the 0.5 factor.
  ! (6) dmisfit_stop indicates that you stop evaluating the misfit function at this value.

  gdot_total = gdot_mu_total + gdot_kappa_total

  nu_step       = 2.0 * (dmisfit_stop - 0.5*dmisfit) / gdot_total
  nu_step_mu    = 2.0 * (dmisfit_stop - 0.5*dmisfit) / gdot_mu_total
  nu_step_kappa = 2.0 * (dmisfit_stop - 0.5*dmisfit) / gdot_kappa_total

  !kfac       = 2.0 * sigma_structure**2 * dVfac * nu_step / dcov_fac
  !kfac_mu    = 2.0 * sigma_structure**2 * dVfac * nu_step_mu / dcov_fac
  !kfac_kappa = 2.0 * sigma_structure**2 * dVfac * nu_step_kappa / dcov_fac

  kfac       = -0.1495e8
  kfac_mu    = -0.1517e8
  kfac_kappa = -0.1011e10

  !-----------------------------------------------

  kernel_name = 'mu_kernel_smooth'
  write(k_file,'(a,i6.6,a)') 'inout'//'/proc',myrank,'_'//trim(kernel_name)//'_dm.bin'
  open(12,file=trim(k_file),form='unformatted')
  write(12) kernel_mu(:,:,:,1:NSPEC) * kfac
  close(12)

  kernel_name = 'kappa_kernel_smooth'
  write(k_file,'(a,i6.6,a)') 'inout'//'/proc',myrank,'_'//trim(kernel_name)//'_dm.bin'
  open(12,file=trim(k_file),form='unformatted')
  write(12) kernel_kappa(:,:,:,1:NSPEC) * kfac
  close(12)

  if (myrank == 0) then
     filename = 'gdot_all'
     open(19,file=trim(filename),status='unknown')
     write(19,'(3e24.12)') gdot_mu_total, gdot_kappa_total, gdot_total
     close(19)

     filename = 'nu_step'
     open(20,file=trim(filename),status='unknown')
     write(20,'(3e24.12)') nu_step, nu_step_mu, nu_step_kappa
     close(20)

     filename = 'kfac'
     open(21,file=trim(filename),status='unknown')
     write(21,'(3e24.12)') kfac, kfac_mu, kfac_kappa
     close(21)
  endif

  !-----------------------------------------------

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program cg_test


