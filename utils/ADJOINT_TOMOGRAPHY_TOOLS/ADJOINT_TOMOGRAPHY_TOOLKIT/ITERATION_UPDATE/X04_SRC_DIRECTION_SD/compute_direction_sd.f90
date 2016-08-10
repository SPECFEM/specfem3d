program xcompute_direction_sd
  implicit none

  include 'mpif.h'
  include '../../SHARE_FILES/HEADER_FILES/constants.h'
  include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h'
  include '../../SHARE_FILES/HEADER_FILES/precision.h'

  integer,parameter:: NSPEC=NSPEC_CRUST_MANTLE
  integer,parameter:: NKERNEL=4
  integer:: myrank, sizeprocs,ier
  integer:: iker,ispec,i,j,k

  character(len=512):: direction_dir, gradient_dir
  character(len=512):: direction_file, gradient_file
  character(len=256):: kernel_name(NKERNEL)

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: direction,gradient

  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  call getarg(1,direction_dir)
  call getarg(2,gradient_dir)

  if (trim(direction_dir) == '' .or. trim(gradient_dir) == '') then
        call exit_MPI(myrank,'USAGE: xcompute_direction_sd direction_dir gradient_dir')
  endif


  kernel_name=(/"reg1_bulk_betah_kernel_precond_smooth","reg1_bulk_betav_kernel_precond_smooth","reg1_bulk_c_kernel_precond_smooth","reg1_eta_kernel_precond_smooth"/)

  do iker = 1,NKERNEL
     direction=0._CUSTOM_REAL

     write(direction_file,'(a,i6.6,a)') trim(direction_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
     write(gradient_file,'(a,i6.6,a)') trim(gradient_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'


     open(1001,file=trim(gradient_file),status='old',form='unformatted',iostat=ier)
     if (myrank == 0) print *, 'reading gradient1:',trim(gradient_file)
     if (ier /= 0) then
        print *, 'error opening:',trim(gradient_file)
        call exit_mpi(myrank,'file not found')
     endif
     read(1001) gradient(:,:,:,1:NSPEC)
     close(1001)

     direction=-gradient

     open(1001,file=trim(direction_file),form='unformatted',action='write')
     if (myrank == 0) print *, 'writing direction1:',direction_file
     write(1001) direction
     close(1001)

  enddo ! kernel type

call MPI_FINALIZE(ier)

end program xcompute_direction_sd
