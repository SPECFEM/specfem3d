program adj_seismogram

  implicit none

  include "mpif.h"

  integer :: myrank,ier

  integer :: ios,irec,i,j,NSTEP,NPROC
  double precision :: DT
  real(kind=4),dimension(:),allocatable :: syn,dat,adj
  integer(kind=4) :: r4head(60)
  integer(kind=4) :: header4(1)
  integer(kind=2) :: header2(2)
  equivalence(header2,header4)
  character(len=512) :: filename,procname,arg

  call MPI_Init(ier)
  call MPI_Comm_Rank(MPI_COMM_WORLD,myrank,ier)
  call MPI_Comm_Size(MPI_COMM_WORLD,NPROC,ier)

  write(procname,"(i4)") myrank

  !! input parameters
  if( iargc() .ne. 2 ) stop 'Usage: ./xadj NSTEP DT'
  j=1;  call getarg(j, arg); read(arg,*,iostat=ios) NSTEP;   if (ios /= 0) stop 'Error reading NSTEP'
  j=2;  call getarg(j, arg); read(arg,*,iostat=ios)    DT;   if (ios /= 0) stop 'Error reading DT'

  !!!! read NSTEP from seismograms
  !!!filename=trim(adjustl(procname))//"_dx_SU"
  !!!open(111,file="../in_out_files/OUTPUT_FILES/"//trim(filename),access='direct',recl=240,iostat = ios)
  !!!read(111,rec=1,iostat=ios) r4head
  !!!close(111)
  !!!header4=r4head(29)
  !!!NSTEP=header2(2)
  !!!header4=r4head(30)
  !!!DT=header2(1)*1.0d-6
  !!!print*, 'irec=',r4head(1)
  !!!print*, 'xs=',r4head(19)
  !!!print*, 'zs=',r4head(20)
  !!!print*, 'xr=',r4head(21)
  !!!print*, 'zr=',r4head(22)
  !!!print*, 'NSTEP=',NSTEP
  !!!print*, "DT=",DT

  allocate(syn(NSTEP),dat(NSTEP),adj(NSTEP))

! read 'syn', 'dat' and then write 'adj'
  ! x-component
  filename=trim(adjustl(procname))//"_dx_SU"
  open(111,file="../in_out_files/SEM/syn/"//trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
  open(112,file="../in_out_files/SEM/dat/"//trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
  open(113,file="../in_out_files/SEM/"//trim(filename)//".adj",access='direct',recl=240+4*NSTEP,iostat = ios)
  irec=1
  do while(ios==0)
     read(111,rec=irec,iostat=ios) r4head,syn
              if (ios /= 0) cycle
     read(112,rec=irec,iostat=ios) r4head,dat
              if (ios /= 0) cycle
     write(113,rec=irec,iostat=ios) r4head,syn-dat
              if (ios /= 0) cycle
     irec=irec+1
  enddo
  close(111)
  close(112)
  close(113)
  ! y-component
  filename=trim(adjustl(procname))//"_dy_SU"
  open(111,file="../in_out_files/SEM/syn/"//trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
  open(112,file="../in_out_files/SEM/dat/"//trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
  open(113,file="../in_out_files/SEM/"//trim(filename)//".adj",access='direct',recl=240+4*NSTEP,iostat = ios)
  irec=1
  do while(ios==0)
     read(111,rec=irec,iostat=ios) r4head,syn
              if (ios /= 0) cycle
     read(112,rec=irec,iostat=ios) r4head,dat
              if (ios /= 0) cycle
     write(113,rec=irec,iostat=ios) r4head,syn-dat
              if (ios /= 0) cycle
     irec=irec+1
  enddo
  close(111)
  close(112)
  close(113)
  ! z-component
  filename=trim(adjustl(procname))//"_dz_SU"
  open(111,file="../in_out_files/SEM/syn/"//trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
  open(112,file="../in_out_files/SEM/dat/"//trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
  open(113,file="../in_out_files/SEM/"//trim(filename)//".adj",access='direct',recl=240+4*NSTEP,iostat = ios)
  irec=1
  do while(ios==0)
     read(111,rec=irec,iostat=ios) r4head,syn
              if (ios /= 0) cycle
     read(112,rec=irec,iostat=ios) r4head,dat
              if (ios /= 0) cycle
     write(113,rec=irec,iostat=ios) r4head,syn-dat
              if (ios /= 0) cycle
     irec=irec+1
  enddo
  close(111)
  close(112)
  close(113)

  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  call MPI_Finalize(ier)

end program adj_seismogram
