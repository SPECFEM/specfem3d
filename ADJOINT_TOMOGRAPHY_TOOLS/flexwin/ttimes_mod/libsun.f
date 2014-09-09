      subroutine warn(msg)
      character*(*) msg
      write(*,100) msg
 100  format(1x,a)
      return
      end
      subroutine tnoua(ia)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine tnoua writes the character string ia to the standard 
c   output without the trailing newline (allowing user input on the 
c   same line).  Programmed on 17 September 1980 by R. Buland.
c
      save
      character*(*) ia
      write(*,100)ia
 100  format(a,$)
      return
      end
      subroutine dasign(lu,mode,ia,len)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine dasign opens (connects) logical unit lu to the disk file
c   named by the character string ia with mode mode.  If iabs(mode) = 1,
c   then open the file for reading.  If iabs(mode) = 2, then open the
c   file for writing.  If iabs(mode) = 3, then open a scratch file for
c   writing.  If mode > 0, then the file is formatted.  If mode < 0,
c   then the file is unformatted.  All files opened by dasign are
c   assumed to be direct access.  Programmed on 3 December 1979 by
c   R. Buland.
c
      save
      character*(*) ia
      logical exst
c
      if(mode.ge.0) nf=1
      if(mode.lt.0) nf=2
      ns=iabs(mode)
      if(ns.le.0.or.ns.gt.3) ns=3
      go to (1,2),nf
 1    go to (11,12,13),ns
 11   open(lu,file=ia,status='old',form='formatted',
     1 access='direct',recl=len)
      return
 12   inquire(file=ia,exist=exst)
      if(exst) go to 11
 13   open(lu,file=ia,status='new',form='formatted',
     1 access='direct',recl=len)
      return
 2    go to (21,22,23),ns
 21   open(lu,file=ia,status='old',form='unformatted',access='direct',
     1 recl=len)
      return
 22   inquire(file=ia,exist=exst)
      if(exst) go to 21
 23   open(lu,file=ia,status='new',form='unformatted',access='direct',
     1 recl=len)
      return
      end
      subroutine vexit(ierr)
      call exit(ierr)
      end
