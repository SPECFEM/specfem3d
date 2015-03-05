! This program is used to precondtion misfit kernels
! Last modified: Fri Sep 14 10:35:58 EDT 2012


program precond_kernels
  implicit none
  include 'mpif.h'
  include '../../SHARE_FILES/HEADER_FILES/constants.h'
  include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h'
  include '../../SHARE_FILES/HEADER_FILES/precision.h'

  integer,parameter:: NSPEC=NSPEC_CRUST_MANTLE
  integer,parameter:: NKERNEL=4    !bulk_betah, bulk_betav, bulk_c, eta
  real(kind=CUSTOM_REAL),parameter:: THRESHOLD_HESS=5.e-4

  integer:: myrank, sizeprocs,ier
  integer:: ios,iker
  character(len=350):: kernel_file,hess_file,input_dir
  character(len=150):: kernel_name(NKERNEL)
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: kernel,hess,kernel_precond
  real(kind=CUSTOM_REAL):: maxh,maxh_all

  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)


  call getarg(1,input_dir)

  if (trim(input_dir) == '' ) then
     call exit_MPI(myrank,'USAGE: xprecond_kernels input_dir')
  endif

  if (myrank == 0) then
     write(*,*) 'PRECONDITION MISFIT KERNELS'
  endif

  kernel_name=(/"reg1_bulk_betah_kernel","reg1_bulk_betav_kernel","reg1_bulk_c_kernel","reg1_eta_kernel"/)


  do iker=1,NKERNEL
     kernel_precond=0.
     hess=0.
     kernel=0.

     write(kernel_file,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
     write(hess_file,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_reg1_hess_kernel.bin'

     if (myrank==0) then
        write(*,*) 'READING IN KERNELS:',trim(kernel_name(iker))
        write(*,*) 'KERNEL FILES:',trim(kernel_file)
        write(*,*) 'HESSIAN FILES:',trim(hess_file)
     endif

     open(unit=1002,file=trim(kernel_file),status='old',form='unformatted')
     read(1002) kernel(:,:,:,1:NSPEC)
     close(1002)


     open(unit=1002,file=trim(hess_file),status='old',form='unformatted')
     read(1002) hess(:,:,:,1:NSPEC)
     close(1002)
     maxh=maxval(hess)

     call mpi_allreduce(maxh,maxh_all,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,ier)
     if ( maxh_all < 1.e-18 ) then
          call MPI_ABORT(MPI_COMM_WORLD,30,ier)
     endif

     if (myrank==0) write(*,*) 'MAX HESSIAN FOR ALL PROCESSORS:',maxh_all

     ! normalized hess
     hess=hess/maxh_all

     where(hess(:,:,:,:) > THRESHOLD_HESS )
          hess=1.0_CUSTOM_REAL / hess
     elsewhere
          hess=1.0_CUSTOM_REAL / THRESHOLD_HESS
     endwhere

     kernel_precond=kernel*hess


     if (myrank==0) write(*,*) 'WRITING OUT PRECONDITIONED KERNEL FOR:',trim(kernel_name(iker))
     write(kernel_file,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'_precond.bin'

     open(1002,file=trim(kernel_file),form='unformatted')
     write(1002) kernel_precond(:,:,:,1:NSPEC)
     close(1002)

     if (iker==1) then
        write(hess_file,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_reg1_hess_kernel_precond.bin'
        open(1002,file=trim(hess_file),form='unformatted')
        write(1002) hess(:,:,:,1:NSPEC)
        close(1002)
     endif
  enddo

  if (myrank==0) write(*,*) 'DONE PRECONDITION ALL MISFIT KERNELS'

  call MPI_FINALIZE(ier)

end program precond_kernels
