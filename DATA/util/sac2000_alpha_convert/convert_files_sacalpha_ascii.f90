
  program sac_convert

!! DK DK convert sac2000 ALPHA files to regular ASCII (e.g. for Gnuplot)

  implicit none

  integer it,nstep,istep

  double precision dt,t0,val1,val2,val3,val4,val5
  double precision dummy1,dummy2,dummy3,dummy4

! read header, ignore the rest of the lines

  read(5,*) dt
  read(5,*) t0

  do it=1,13
    read(5,*)
  enddo

  read(5,*) dummy1,dummy2,dummy3,dummy4,nstep

  do it=1,14
    read(5,*)
  enddo

!! DK DK in sac2000 ALPHA format, there are 5 ASCII values per line
  it = 0
  do istep = 1,nstep/5
    read(5,*) val1,val2,val3,val4,val5

    it = it + 1
    write(*,*) sngl(t0 + (it-1)*dt),sngl(val1)

    it = it + 1
    write(*,*) sngl(t0 + (it-1)*dt),sngl(val2)

    it = it + 1
    write(*,*) sngl(t0 + (it-1)*dt),sngl(val3)

    it = it + 1
    write(*,*) sngl(t0 + (it-1)*dt),sngl(val4)

    it = it + 1
    write(*,*) sngl(t0 + (it-1)*dt),sngl(val5)

  enddo

  end program sac_convert

