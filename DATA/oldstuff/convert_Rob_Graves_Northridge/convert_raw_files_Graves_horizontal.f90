
  program convert_raw_files_Wald

! rotates the horizontal components of Rob Graves' data for Northridge

! with Intel's ifc compiler, need to link with "-Vaxlib"

  implicit none

  character(len=15) arg1,arg2,arg3

  integer i,iargc,numarg,iangle1,iangle2,ianglebak,length1,length2,itarget_angle

! number of values per line in Rob Graves' data files
  integer, parameter :: NVAL = 6

  integer nstep,hour,minutes,nlines,iline,it,ival,NSTEP_OUTPUT

  double precision a1(NVAL),a2(NVAL)
  double precision seconds,dt,initial_time,rotation_angle

  double precision, dimension(:), allocatable :: time,amplitude1,amplitude2,amplitudebak,amplitude_east,amplitude_north

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0

  numarg = iargc()

  if(numarg /= 3) stop 'need exactly 3 arguments'

! get the 3 arguments
  call getarg(1,arg1)
  call getarg(2,arg2)
  call getarg(3,arg3)

! length of arguments
  length1 = len_trim(arg1)
  length2 = len_trim(arg2)

! get digits at the end of name of stations, name is always 000, 090, 175 etc.
! convert it to integer by writing to a file and reading back
  open(unit=11,file='____tutu',status='unknown')
  write(11,*) arg1(length1-2:length1)
  write(11,*) arg2(length2-2:length2)
  close(11)

  open(unit=11,file='____tutu',status='old')
  read(11,*) iangle1
  read(11,*) iangle2
  close(11)

!! DK DK particular case of 360 degree angle
  if(iangle1 == 360) iangle1 = 0
  if(iangle2 == 360) iangle2 = 0

  write(*,*)
  write(*,*) arg1,arg2,arg3

! read both data files
  open(unit=11,file=arg1,status='old')
  open(unit=12,file=arg2,status='old')

! skip title
  read(11,*)
  read(12,*)

! read timing info
  read(11,*) nstep,dt,hour,minutes,seconds
  read(12,*)

! compute number of lines of NVAL values
  nlines = nstep / NVAL

! compute number of time steps (can be different because of integer division)
  NSTEP_OUTPUT = nlines * NVAL

! allocate arrays
  allocate(time(NSTEP_OUTPUT))
  allocate(amplitude1(NSTEP_OUTPUT))
  allocate(amplitude2(NSTEP_OUTPUT))
  allocate(amplitudebak(NSTEP_OUTPUT))
  allocate(amplitude_east(NSTEP_OUTPUT))
  allocate(amplitude_north(NSTEP_OUTPUT))

! initialize time step number
  it = 0

! initial time
  initial_time = 3600 * hour + 60 * minutes + seconds

  print *,'initial time = ',initial_time

! loop on all the lines
  do iline = 1,nlines

    read(11,*) (a1(ival),ival=1,NVAL)
    read(12,*) (a2(ival),ival=1,NVAL)

! unit is cm/s/s, convert to m/s/s (i.e., SI)

    do ival = 1,NVAL
      it = it + 1
      time(it) = (it-1)*dt + initial_time
      amplitude1(it) = a1(ival)/100.d0
      amplitude2(it) = a2(ival)/100.d0
    enddo

  enddo

  close(11)
  close(12)

! swap components if needed, if angles are inverted on command line
  itarget_angle = iangle1 - 90
  if(itarget_angle < 0) itarget_angle = itarget_angle + 360
  if(iangle2 /= itarget_angle) then
    print *,'swapping components **************'
    ianglebak = iangle1
    amplitudebak = amplitude1
    iangle1 = iangle2
    amplitude1 = amplitude2
    iangle2 = ianglebak
    amplitude2 = amplitudebak
  else
    print *,'not swapping components'
  endif

! now check that angles are always correct
  write(*,*) iangle1,iangle2
  itarget_angle = iangle1 - 90
  if(itarget_angle < 0) itarget_angle = itarget_angle + 360
  if(iangle2 /= itarget_angle) stop 'difference of angles is not 90 degrees'

!! DK DK rotation from Rob Graves' azimuth to SEM East/North/up convention
  rotation_angle = - (90 - iangle1) * PI / 180.d0
print *,'rotating radians ',rotation_angle
  do it=1,NSTEP_OUTPUT
    amplitude_east(it) = cos(rotation_angle)*amplitude1(it) + sin(rotation_angle)*amplitude2(it)
    amplitude_north(it) = - sin(rotation_angle)*amplitude1(it) + cos(rotation_angle)*amplitude2(it)
  enddo

!! DK DK write final output files
  open(unit=12,file=arg3(1:len_trim(arg3))//'.DW.BHE.dataa',status='unknown')
  do it=1,NSTEP_OUTPUT
    write(12,*) sngl(time(it)),sngl(amplitude_east(it))
  enddo
  close(12)

  open(unit=12,file=arg3(1:len_trim(arg3))//'.DW.BHN.dataa',status='unknown')
  do it=1,NSTEP_OUTPUT
    write(12,*) sngl(time(it)),sngl(amplitude_north(it))
  enddo
  close(12)

  end program convert_raw_files_Wald

