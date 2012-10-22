!
! program to create adjoint source files in SU (Seismic Unix) format
!
! this template program uses SU-formatted files.
! by default, it takes a pure waveform misfit to construct adjoint sources (adj = syn - dat)
! please modify routine create_adjoint_trace() according to your needs.
!
! for example, compile with:
!  > gfortran -o ../../example/homogeneous/halfspace/bin/xSU_adjoint SU_adjoint.f90
!
! call in example directory example/homogeneous_halfspace/ by:
!  > ./bin/xSU_adjoint 4
!
!-------------------------------------------------------------------------------------------------

  module SU_adjoint_var

  implicit none

  real(kind=4),dimension(:),allocatable :: syn,dat,adj

  character(len=128) :: DATA_PATH = "./DATA/"
  character(len=128) :: SYN_PATH = "./OUTPUT_FILES/"
  character(len=128) :: ADJ_PATH = "./OUTPUT_FILES/SEM/"

  end module

!-------------------------------------------------------------------------------------------------

  program SU_adjoint

  use SU_adjoint_var
  implicit none

  integer :: NPROC
  integer :: iproc,icomp,ios,irec,i

  character(len=512) :: filename
  character(len=512) :: arg(1)
  character(len=16) :: procname
  character(len=1),dimension(3) :: compstr = (/ "x" , "y" , "z" /)

  ! SU header
  integer(kind=4) :: r4head(60)
  integer(kind=4) :: header4(1)
  integer(kind=2) :: header2(2)
  equivalence(header2,header4)

  integer :: NSTEP
  double precision :: DT

  ! reads in file arguments
  i = 1
  do while (1 == 1)
    call getarg(i,arg(i))
    if (i <= 1 .and. trim(arg(i)) == '') then
      print*,'Usage: '
      print*,'  ./xSU_adjoint NPROC'
      print*,'with'
      print*,'  NPROC: total number of partitions'
      stop
    endif
    if (trim(arg(i)) == '') exit
    if (i == 1) then
      read(arg(i),*,iostat=ios) NPROC
      if (ios /= 0) stop 'Error reading NPROC'
    endif
    i = i+1
  enddo

  ! user output
  print*, "SU adjoint "
  print*
  print*, "number of partitions : ",NPROC
  print*, "data path            : ",trim(DATA_PATH)
  print*, "synthetics path      : ",trim(SYN_PATH)
  print*, "adjoint sources path : ",trim(ADJ_PATH)
  print*

  ! reads NSTEP info from seismograms
  ! only do this once
  iproc = 0
  write(procname,"(i4)") iproc
  filename = trim(DATA_PATH)//trim(adjustl(procname))//"_dx_SU"
  open(11,file=trim(filename),access='direct',status='old', &
        recl=240,iostat=ios)
  if( ios /= 0 ) then
    print*,'error opening file: ',trim(filename)
    stop 'error opening data file'
  endif

  ! reads in header
  read(11,rec=1,iostat=ios) r4head
  close(11)

  ! extracts header information
  header4 = r4head(29)
  NSTEP = header2(2)
  header4 = r4head(30)
  DT = header2(1) * 1.0d-6

  ! user output
  print*, "NSTEP = ",NSTEP
  print*, "DT    = ",DT
  print*

  ! allocates arrays
  allocate(syn(NSTEP),dat(NSTEP),adj(NSTEP))

  ! loops over all partitions
  do iproc = 0, NPROC - 1

    ! partition name
    write(procname,"(i4)") iproc
    print*,"...reading partition: ",iproc

    ! loops over components
    do icomp = 1,3
      ! user output
      print*,"  component ",compstr(icomp)

      ! opens files
      ! data
      filename = trim(DATA_PATH)//trim(adjustl(procname))//"_d"//compstr(icomp)//"_SU"
      open(11,file=trim(filename),access='direct',status='old', &
            recl=240+4*NSTEP,iostat=ios)
      if( ios /= 0 ) then
        print*,'error opening file: ',trim(filename)
        stop 'error opening input data file '
      endif

      ! synthetics
      filename = trim(SYN_PATH)//trim(adjustl(procname))//"_d"//compstr(icomp)//"_SU"
      open(22,file=trim(filename),access='direct',status='old', &
            recl=240+4*NSTEP,iostat=ios)
      if( ios /= 0 ) then
        print*,'error opening file: ',trim(filename)
        stop 'error opening input file '
      endif

      ! adjoint sources
      filename = trim(ADJ_PATH)//trim(adjustl(procname))//"_d"//compstr(icomp)//"_SU"//".adj"
      open(33,file=trim(filename),access='direct',status='unknown', &
            recl=240+4*NSTEP,iostat = ios)
      if( ios /= 0 ) then
        print*,'error opening file: ',trim(filename)
        stop 'error opening output file '
      endif

      ! loops over all records
      irec=1
      do while(ios==0)
        ! reads in data
        read(11,rec=irec,iostat=ios) r4head,dat
        if (ios /= 0) cycle

        ! reads in synthetics
        read(22,rec=irec,iostat=ios) r4head,syn
        if (ios /= 0) cycle

        ! creates adjoint trace
        call create_adjoint_trace()

        ! writes out adjoint source
        write(33,rec=irec,iostat=ios) r4head,adj
        if (ios /= 0) cycle

        irec=irec+1
      enddo
      close(11)
      close(22)
      close(33)
    enddo

    ! user output
    print*, "  receivers read: ",irec
  enddo

  ! user output
  print*, "done adjoint sources"
  print*, "please, check output files in directory: ",trim(ADJ_PATH)

  end program

!-------------------------------------------------------------------------------------------------

  subroutine create_adjoint_trace()

  use SU_adjoint_var
  implicit none

  ! MODIFY THIS according to your misfit function

  ! waveform misfit:
  adj(:) = syn(:) - dat(:)

  end subroutine
