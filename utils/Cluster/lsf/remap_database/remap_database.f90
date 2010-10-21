program remap_databases

! this program remaps the databases to resume a solver calculation after
! the mesher 
! Qinya Liu, May 2007, Caltech

implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer, parameter :: MAX_PROCS = 1000
  integer ier,sizeprocs,myrank,ios,i
  character(len=256) old_machine_file,junk,junk2,slice_to_old_machine(MAX_PROCS), &
             mymachine, local_data_base, scp_outfile, command_string

  integer num_slices, num_slices2,num

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

! run the main program
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  call getarg(1,old_machine_file)
  if (trim(old_machine_file) == '') call exit_mpi(myrank,'Usage: remap_database old-mach num-slice')
  call getarg(2,junk)
  if (trim(junk) == '') call exit_mpi(myrank,'Usage: remap_database old-mach num-slice')
  read(junk,*) num_slices
  if (num_slices /= sizeprocs) call exit_mpi(myrank,'number of slices does not match')

  num_slices2 = num_slices
  open(11,file=trim(old_machine_file),status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening old machine file'
  do while (1 == 1)
    read(11,'(a)',iostat=ios) junk2
    if (ios /= 0) exit
    read(11,'(i)',iostat=ios) num
    if (ios /= 0) exit
    if (num /= 2) stop 'Number of slices per node not 2'
    do i = 1, num
      slice_to_old_machine(num_slices2-i+1) = junk2
    enddo
    num_slices2 = num_slices2 - num
  enddo
  if (num_slices2 /= 0) stop 'Error counting number of slices'
  close(11)

  mymachine = slice_to_old_machine(myrank+1)
 
  local_data_base = '/scratch/lqy/DATABASES_MPI'

  write(scp_outfile,'(a,i4.4)') 'OUTPUT_FILES/scp_out.',myrank

  write(command_string,'(a,i4.4,a)') 'scp lqy@'//trim(mymachine)//':'//trim(local_data_base)//'/*', &
             myrank, '*  '//trim(local_data_base)

  call system('echo '//trim(command_string)//' > '//trim(scp_outfile))
  
  call system(trim(command_string)//' >> '//trim(scp_outfile))
           
! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program remap_databases
          
