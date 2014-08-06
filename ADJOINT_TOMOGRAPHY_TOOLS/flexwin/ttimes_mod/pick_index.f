      program getindex

      character filename*125,ktname*8,phase*8 
      real x(100000), b, dt
      integer npts, ioerror,nchar,nphase

      call getarg(1,filename)
      call getarg(2,phase)
      nphase=lnblnk(phase)
      call rsac1(filename,x,npts,b,dt,100000,ioerror)

c     check kt0
      call getkhv('kt0',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '0'
      endif
	
c     check kt1
      call getkhv('kt1',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '1'
      endif
	
c     check kt2
      call getkhv('kt2',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '2'
      endif
	
c     check kt3
      call getkhv('kt3',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '3'
      endif
	
c     check kt4
      call getkhv('kt4',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '4'
      endif
	
c     check kt5
      call getkhv('kt5',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '5'
      endif
	
c     check kt6
      call getkhv('kt6',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '6'
      endif
	
c     check kt7
      call getkhv('kt7',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '7'
      endif
	
c     check kt8
      call getkhv('kt8',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '8'
      endif
	
c     check kt9
      call getkhv('kt9',ktname,ioerror)
      nchar=lnblnk(ktname)
      if(nchar .eq. nphase .and. ktname .eq. phase ) then
	write(*,*) '9'
      endif
	
c     none match
      write (*,*) 'Cannot find a match for ',phase

      end
