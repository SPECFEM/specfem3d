
	program txt_2_sep

c Parameters to be defined to dimension the arrays:
c - nlmax = maximum number of layers.
c - nr = maximum number of receivers.
c - ntime = maximum number of points of each seismogram.
      parameter (nr=56)
      parameter (nrtot=56+56)

      parameter (ntime=512)

      dimension sy(ntime,nrtot)

c ecriture au format text deplacement 
 
      print *,'Reading seismograms in text format'
 
      open(unit=11,file='U_file.txt',status='unknown')
	do ir=1,nr
	do it=1,ntime
      	read(11,*) sy(it,ir+nr)

c eventuellement inverser le signe pour compatibilite 
      	sy(it,ir+nr) = + sy(it,ir+nr)

	enddo
	enddo
      close(11)
 
c ecriture au format binaire deplacement

	print *,'Saving seismograms in SEP format'

      open(unit=11,file='U_file',status='unknown',form='unformatted')
      write(11) sy
      close(11)

      stop                                                                      
      end                                                                       

