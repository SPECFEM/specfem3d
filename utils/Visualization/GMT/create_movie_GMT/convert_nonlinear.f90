program convert_nonlinear

  ! reads in the output of shakemaps
  ! and convert it to non-linear scale

  implicit none

  character(len=256) :: infile,outfile,ch_power,ch_zmax
  real*8 :: power,zmax,zmin
  logical :: calculate_zmax
  integer :: ios,i,j
  real*8 :: lon,lat,zv



  call getarg(1,infile)
  call getarg(2,outfile)
  call getarg(3,ch_power)
  call getarg(4,ch_zmax)

  if (trim(infile) == '' .or. trim(outfile) == '' ) then
     stop 'Usage: convert_nonlinear infile outfile power [zmax]'
  endif
  if (trim(ch_power) == '') then
     power = 0.25
  else
     read(ch_power,*) power
  endif
  if (trim(ch_zmax) == '') then
     calculate_zmax = .true.
  else
     calculate_zmax = .false.
     read(ch_zmax,*) zmax
  endif


  open(11,file = infile,status = 'old',iostat = ios)
  if (ios /= 0) stop 'Error open input xyz file'
  open(12,file = outfile,status='unknown',iostat = ios)
  if (ios /= 0) stop 'Error open output xyz file'
  i = 0
  if (calculate_zmax) then
     do while (ios == 0)
        read(11,*,iostat=ios) lon,lat,zv
        if (ios /= 0) exit
        i = i + 1
        if (i==1) then
           zmax = zv
        else if (zmax < zv) then
           zmax = zv
        endif
     enddo
     print *, 'The maximum value of the dataset is ', sngl(zmax)
     close(11)
     open(11,file = infile,status = 'old',iostat = ios)
     if (ios /= 0) stop 'Error open input xyz file'
  endif

  do while (ios == 0)
     read(11,*,iostat=ios) lon,lat,zv
     if (ios /= 0) exit
     write(12,'(f18.12,f16.12,e19.12)',iostat=ios) lon,lat,(zv/zmax)**power
     if (ios /= 0) stop 'Error writing'
  enddo
  close(12)
  close(11)

end program convert_nonlinear
