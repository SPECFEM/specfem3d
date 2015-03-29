program xcopy_local

  implicit none
  include 'mpif.h'

  integer myrank,sizeprocs,ier
  integer iregion_code
  character(len=256) :: command,filename,procname
  character(len=256) :: databases,scratchdir
  character(len=8) cp

  ! forward mesh databases
  !databases = "/tigress-hsm/dpeter/SPECFEM3D_GLOBE/DATABASES_MPI.middle_east"

  ! forward mesh databases (better filesystem for parallel i/o)
  databases = "/scratch/lustre/hejunzhu/XEUROPE_MODEL/MESH"

  ! local node scratch directory
  scratchdir = "/scratch/hejunzhu/"

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

  ! run the main program
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  ! creates local scratch directory
  !command = 'rm -rf '//trim(scratchdir)
  !call system(trim(command))
  !command = 'mkdir -p '//trim(scratchdir)
  !call system(trim(command))

  ! copy command
  if( myrank == 0 ) then
    cp = "cp -v"
  else
    cp = "cp"
  endif

  do iregion_code=1,3
    ! process name
    write(procname,"('/proc',i6.6,'_reg',i1,'_')") myrank,iregion_code

    ! database files

    filename = trim(databases)//trim(procname)//'array_dims.txt'

    command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    call system(trim(command))

    filename = trim(databases)//trim(procname)//'attenuation.bin'

    command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    call system(trim(command))

    filename = trim(databases)//trim(procname)//'boundary.bin'

    command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    call system(trim(command))

    filename = trim(databases)//trim(procname)//'ibool*.txt'

    command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    call system(trim(command))

    filename = trim(databases)//trim(procname)//'solver_data_1.bin'

    command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    call system(trim(command))

    filename = trim(databases)//trim(procname)//'solver_data_2.bin'

    command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    call system(trim(command))

    filename = trim(databases)//trim(procname)//'stacey.bin'

    command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    call system(trim(command))


    !filename = trim(databases)//trim(procname)//'absorb*.bin'
    !
    !command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    !call system(trim(command))

    !filename = trim(databases)//trim(procname)//'boundary_disc.bin'
    !
    !command =  trim(cp)//' '//trim(filename)//' '//trim(scratchdir)
    !call system(trim(command))


  enddo

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program
