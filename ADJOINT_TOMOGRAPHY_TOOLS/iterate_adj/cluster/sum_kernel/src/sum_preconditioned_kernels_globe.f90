! sum_preconditioned_kernels_globe
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! this version uses the approximate Hessian stored as
!   - proc***_reg1_hess_kernel.bin
!          to precondition the transverse isotropic kernel files
!
! input file: kernels_run.globe
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernel_list.globe")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory

program sum_preconditioned_kernels_globe

  implicit none
  include 'mpif.h'
  include 'constants_globe.h'
  include 'precision_globe.h'
  include 'values_from_mesher_globe.h'

! ======================================================
  ! USER PARAMETERS

  ! by default, this algorithm uses transverse isotropic (bulk,bulk_betav,bulk_betah,eta) kernels to sum up
  ! if you prefer using isotropic kernels, set flags below accordingly

  ! if you prefer using isotropic kernels (bulk,bulk_beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ISO_KERNELS = .false.

  ! if you prefer isotropic  (alpha,beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ALPHA_BETA_RHO = .false.


! ======================================================

  character(len=150) :: kernel_file_list, kernel_list(1000), sline, kernel_name
  integer :: nker, myrank, sizeprocs,  ier
  integer :: ios

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if(myrank==0) then
    write(*,*) 'sum_preconditioned_kernels_globe:'
    write(*,*)
    write(*,*) 'reading kernel list: '
  endif

  ! default values
  kernel_file_list='kernels_run.globe'

  ! reads in event list
  nker=0
  open(unit = 20, file = trim(kernel_file_list), status = 'old',iostat = ios)
  if (ios /= 0) then
     print *,'Error opening ',trim(kernel_file_list),myrank
     stop 1
  endif
  do while (1 == 1)
     read(20,'(a)',iostat=ios) sline
     if (ios /= 0) exit
     nker = nker+1
     kernel_list(nker) = sline
  enddo
  close(20)
  if( myrank == 0 ) then
    write(*,*) '  ',nker,' events'
    write(*,*)
  endif

  ! sums up kernels
  if( USE_ISO_KERNELS ) then

    if( myrank == 0 ) write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'

    !  isotropic kernels
    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_beta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  else if( USE_ALPHA_BETA_RHO ) then

    ! isotropic kernels
    if( myrank == 0 ) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'reg1_alpha_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_beta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  else

    ! transverse isotropic kernels
    if( myrank == 0 ) write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah,eta'

    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betav_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betah_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_eta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  endif

  if(myrank==0) write(*,*) 'done writing all kernels'

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program sum_preconditioned_kernels_globe

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  implicit none

  include 'mpif.h'
  include 'constants_globe.h'
  include 'precision_globe.h'
  include 'values_from_mesher_globe.h'

  !----------------------------------------------------------------------------------------
  ! USER PARAMETER
  ! 1 permille of maximum for inverting hessian
  real(kind=CUSTOM_REAL),parameter :: THRESHOLD_HESS = 1.e-3
  ! sums all hessians before inverting and preconditioning
  logical, parameter :: USE_HESS_SUM = .true.
  ! uses source mask to blend out source elements
  logical, parameter :: USE_SOURCE_MASK = .false.
  
  !----------------------------------------------------------------------------------------
  real(kind=CUSTOM_REAL) :: norm,norm_sum
  character(len=150) :: kernel_name,kernel_list(1000)
  integer :: nker,myrank,ier

  ! local parameters
  character(len=150) :: k_file
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    kernel_crust_mantle,hess_crust_mantle,total_kernel
  integer :: iker,ios
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: total_hess,mask_source

  ! initializes arrays
  if( USE_HESS_SUM ) then
  
    allocate( total_hess(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) )
    total_hess(:,:,:,:) = 0.0_CUSTOM_REAL
    
  endif
  
  if( USE_SOURCE_MASK ) then
  
    allocate( mask_source(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) )
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL  
    
  endif
  
  ! loops over all event kernels
  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if(myrank==0) then
      write(*,*) 'reading in event kernel for: ',trim(kernel_name)
      write(*,*) '  and preconditioner: ','reg1_hess_kernel'
      write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel / frechet derivative
    kernel_crust_mantle = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,'_'//trim(kernel_name)//'.bin'

    open(12,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
    if( ios /= 0 ) then
      write(*,*) '  kernel not found:',trim(k_file)
      cycle
    endif
    read(12) kernel_crust_mantle
    close(12)

    ! outputs norm of kernel 
    norm = sum( kernel_crust_mantle * kernel_crust_mantle )
    call mpi_reduce(norm,norm_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
    if( myrank == 0 ) then
      print*,'  norm kernel: ',sqrt(norm_sum)
    endif


    ! approximate Hessian
    hess_crust_mantle = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,'_reg1_hess_kernel.bin'

    open(12,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
    if( ios /= 0 ) then
      write(*,*) '  hess not found:',trim(k_file)
    else
      read(12) hess_crust_mantle
      close(12)
    endif

    ! outputs norm of preconditioner
    norm = sum( hess_crust_mantle * hess_crust_mantle )
    call mpi_reduce(norm,norm_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
    if( myrank == 0 ) then
      print*,'  norm preconditioner: ',sqrt(norm_sum)
      print*
    endif

    ! note: we take absolute values for hessian (as proposed by Yang)
    hess_crust_mantle = abs(hess_crust_mantle)
    
    ! source mask
    if( USE_SOURCE_MASK ) then
      ! reads in mask
      write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                            //'/proc',myrank,'_reg1_mask_source.bin'
      open(12,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
      if( ios /= 0 ) then
        write(*,*) '  file not found:',trim(k_file)
        cycle
      endif
      read(12) mask_source
      close(12)
      
      ! masks source elements
      kernel_crust_mantle = kernel_crust_mantle * mask_source
      
    endif

    ! precondition
    if( USE_HESS_SUM ) then

      ! sums up hessians first
      total_hess = total_hess + hess_crust_mantle

    else

      ! inverts hessian
      call invert_hess( myrank,hess_crust_mantle,THRESHOLD_HESS )
      
      ! preconditions each event kernel with its hessian
      kernel_crust_mantle = kernel_crust_mantle * hess_crust_mantle

    endif
    
    ! sums all kernels from each event
    total_kernel = total_kernel + kernel_crust_mantle
      
  enddo

  ! preconditions summed kernels with summed hessians
  if( USE_HESS_SUM ) then
  
      ! inverts hessian matrix
      call invert_hess( myrank,total_hess,THRESHOLD_HESS )
      
      ! preconditions kernel
      total_kernel = total_kernel * total_hess
  
  endif

  ! stores summed kernels
  if(myrank==0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)

  write(k_file,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,'_'//trim(kernel_name)//'.bin'

  open(12,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ios)
  if( ios /= 0 ) then
    write(*,*) 'ERROR kernel not written:',trim(k_file)
    stop 'error kernel write'
  endif
  write(12) total_kernel
  close(12)

  if(myrank==0) write(*,*)

end subroutine sum_kernel_pre

!
!-------------------------------------------------------------------------------------------------
!

subroutine invert_hess( myrank,hess_matrix,THRESHOLD_HESS )

! inverts the hessian matrix
! the approximate hessian is only defined for diagonal elements: like
! H_nn = \frac{ \partial^2 \chi }{ \partial \rho_n \partial \rho_n }
! on all GLL points, which are indexed (i,j,k,ispec)
  
  implicit none
  include 'mpif.h'
  include 'constants_globe.h'
  include 'precision_globe.h'
  include 'values_from_mesher_globe.h'

  integer :: myrank
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    hess_matrix
  real(kind=CUSTOM_REAL) :: THRESHOLD_HESS

  ! local parameters
  real(kind=CUSTOM_REAL) :: maxh,maxh_all
  integer :: ier
  
  ! maximum value of hessian
  maxh = maxval( abs(hess_matrix) )

  ! determines maximum from all slices on master
  call mpi_allreduce(maxh,maxh_all,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,ier)
  
  ! user output
  if( myrank == 0 ) then
    print*
    print*,'hessian maximum: ',maxh_all
    print*
  endif
  
  ! normalizes hessian 
  if( maxh_all < 1.e-18 ) then
    ! hessian is zero, re-initializes 
    hess_matrix = 1.0_CUSTOM_REAL
    !call exit_mpi(myrank,'error hessian too small')
  else
    ! since hessian has absolute values, this scales between [0,1]
    hess_matrix = hess_matrix / maxh_all  
  endif
  

  ! inverts hessian values
  where( abs(hess_matrix(:,:,:,:)) > THRESHOLD_HESS )
    hess_matrix = 1.0_CUSTOM_REAL / hess_matrix
  elsewhere
    hess_matrix = 1.0_CUSTOM_REAL / THRESHOLD_HESS
  endwhere

  ! rescales hessian
  !hess_matrix = hess_matrix * maxh_all
  
end subroutine invert_hess
