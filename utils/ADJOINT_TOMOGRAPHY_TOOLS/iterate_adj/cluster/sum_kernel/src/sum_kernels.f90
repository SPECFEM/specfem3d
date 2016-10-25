program sum_kernels

  implicit none
  include 'mpif.h'
  include 'constants.h'
  include 'precision.h'
  include 'values_from_mesher.h'

  ! ======================================================

  integer, parameter :: NSPEC=NSPEC_AB_VAL

  character(len=150) :: kernel_file_list, kernel_list(1000), sline, k_file, kernel_name, ftagfile, ftag
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: kernel, total_kernel
  integer :: iker, nker, myrank, sizeprocs,  ier
  integer :: i, j, k,ispec, iglob, ishell, n, it, j1, ib, npts_sem, ios

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  kernel_file_list='kernels_run_sub'
  ftagfile = 'ftags'

  ! read in the list of event IDs
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

  ! read in the name of the kernel
  open(unit = 20, file = trim(ftagfile), status = 'old',iostat = ios)
  if (ios /= 0) then
     print *,'Error opening ',trim(ftagfile)
     stop
  endif
  read(20,'(a)',iostat=ios) ftag
  if (ios /= 0) then
     print *,'Error reading ',trim(ftagfile)
     stop
  endif
  close(20)

  !kernel_name = 'mu_kernel'
  !kernel_name = 'mu_kernel_smooth'
  kernel_name = ftag

  total_kernel=0.
  do iker = 1, nker
     if (myrank == 1) write(*,*) 'reading in event kernel for mu: ', iker, ' out of ', nker
     write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker))//'/proc',myrank,'_'//trim(kernel_name)//'.bin'

     open(12,file=trim(k_file),status='old',form='unformatted')
     read(12) kernel(:,:,:,1:nspec)
     close(12)

     !total_kernel(:,:,:,1:nspec) = total_kernel(:,:,:,1:nspec) + abs( kernel(:,:,:,1:nspec) )
     total_kernel(:,:,:,1:nspec) = total_kernel(:,:,:,1:nspec) + kernel(:,:,:,1:nspec)

  enddo
  if (myrank == 1) write(*,*) 'writing out summed kernel for mu'
  write(k_file,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,'_'//trim(kernel_name)//'.bin'
  open(12,file=trim(k_file),form='unformatted')
  write(12) total_kernel(:,:,:,1:nspec)
  close(12)

!!$  !kernel_name = 'kappa_kernel'
!!$  kernel_name = 'kappa_kernel_smooth'
!!$
!!$  total_kernel=0.
!!$  do iker = 1, nker
!!$     if (myrank==1) write(*,*) 'reading in event kernel for kappa: ', iker, ' out of ', nker
!!$     write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker))//'/proc',myrank,'_'//trim(kernel_name)//'.bin'
!!$
!!$     open(12,file=trim(k_file),status='old',form='unformatted')
!!$     read(12) kernel(:,:,:,1:nspec)
!!$     close(12)
!!$
!!$     !total_kernel(:,:,:,1:nspec) = total_kernel(:,:,:,1:nspec) + abs( kernel(:,:,:,1:nspec) )
!!$     total_kernel(:,:,:,1:nspec) = total_kernel(:,:,:,1:nspec) + kernel(:,:,:,1:nspec)
!!$
!!$  enddo
!!$  if (myrank==1) write(*,*) 'writing out summed kernel for kappa'
!!$  write(k_file,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,'_'//trim(kernel_name)//'.bin'
!!$  open(12,file=trim(k_file),form='unformatted')
!!$  write(12) total_kernel(:,:,:,1:nspec)
!!$  close(12)

  if (myrank == 1) write(*,*) 'done writing all kernels, now finishing...'

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program sum_kernels


