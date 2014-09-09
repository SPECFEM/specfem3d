! sum_kernels_globe
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! input file: kernels_run.globe       
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernel_list.globe")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory

program sum_kernels_globe

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

  ! crust mantle region 1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    kernel_crust_mantle,total_kernel

  character(len=150) :: kernel_file_list, kernel_list(1000), sline, k_file, kernel_name
  integer :: iker, nker, myrank, sizeprocs,  ier
  integer :: i, j, k,ispec, iglob, ishell, n, it, j1, ib, npts_sem, ios

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if(myrank==0) write(*,*) 'reading kernel list: '

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
  
    !  isotropic kernels
    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)
    
  else if( USE_ALPHA_BETA_RHO ) then

    ! isotropic kernels
  
    kernel_name = 'reg1_alpha_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

  else
  
    ! transverse isotropic kernels
    
    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betav_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betah_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_eta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)
  
  endif
  
  if(myrank==0) write(*,*) 'done writing all kernels'

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program sum_kernels_globe

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel(kernel_name,kernel_list,nker,myrank)

  implicit none

  include 'mpif.h'
  include 'constants_globe.h'
  include 'precision_globe.h'
  include 'values_from_mesher_globe.h'

  character(len=150) :: kernel_name,kernel_list(1000)
  integer :: nker,myrank,ier

  ! local parameters
  character(len=150) :: k_file
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    kernel_crust_mantle,total_kernel
  real(kind=CUSTOM_REAL) :: norm,norm_sum
  integer :: iker,ios

  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if(myrank==0) then 
    write(*,*) 'reading in event kernel for: ',trim(kernel_name)
    write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel
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
    print*
    endif

    ! sums all kernels from each event
    total_kernel(:,:,:,1:NSPEC_CRUST_MANTLE_ADJOINT) = &
        total_kernel(:,:,:,1:NSPEC_CRUST_MANTLE_ADJOINT) &
        + kernel_crust_mantle(:,:,:,1:NSPEC_CRUST_MANTLE_ADJOINT)
  enddo

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

end subroutine sum_kernel

