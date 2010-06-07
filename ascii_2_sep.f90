
  program ascii_2_sep

! to convert ASCII seismograms to SEP binary format

  implicit none

! Parameters to be defined to set the dimension of the arrays:
! - nr = maximum number of receivers.
! - ntime = maximum number of points of each seismogram.
  integer, parameter :: nr = 101
  integer, parameter :: ntime = 512

  real(kind=4), dimension(ntime,nr) :: sy

  integer :: it,ir,i

! read seismogram component in ASCII text format
  print *,'Reading seismograms in text format'

  open(unit=11,file='U_file.txt',status='unknown')
  do ir=1,nr
    do it=1,ntime
      read(11,*) sy(it,ir)

! invert sign of vertical component if needed depending on the reference frame convention
!     sy(it,ir) = - sy(it,ir)

    enddo
  enddo
  close(11)

! write seismogram component in binary SEP format
  print *,'Saving seismograms in SEP format'

  open(unit=11,file='U_file.sep',status='unknown',access='direct',recl=ntime*nr*4)
  write(11,rec=1) (sy(i,1),i=1,ntime*nr)
  close(11)

  end program ascii_2_sep

