
	program comp_traces

	implicit none

      character*60 namefile

	integer ntspec, ntbouchon, nrec
	parameter(ntspec = 750, ntbouchon = 512)
	parameter(nrec = 112)

	real tmax_spec,tmax_bouchon
	parameter(tmax_spec = 48.75, tmax_bouchon = 50.)

	real factor
	parameter(factor = 1.e14)

	real dt_spec,dt_bouchon
	integer irec,it

c sismos extraits en surface (profils)
	real sisuz_bouchon(ntbouchon,nrec)
	real sisuz_specfem(ntspec,nrec)

	print *,'**** Lecture sismos ****'
	print *,' ' 

	dt_spec = tmax_spec / real(ntspec)

c DK DK jan 99 apparamment bouchon prend un pas de plus
c	dt_bouchon = tmax_bouchon / real(ntbouchon - 1)
	dt_bouchon = tmax_bouchon / real(ntbouchon)

	print *,' Time step specfem = ',dt_spec
	print *,' Time step bouchon = ',dt_bouchon
	print *,' Number of time steps specfem = ',ntspec
	print *,' Number of time steps bouchon = ',ntbouchon
	print *,' '

c--- lecture au format binaire sismogrammes
      open(unit=11,file='Uy_file_Xcste',status='unknown',
     .    access='direct',recl=ntspec*nrec*4)
      read(11,rec=1) sisuz_specfem
      close(11)

c--- lecture au format binaire sismogrammes
      open(unit=11,file='Ux_file',status='unknown',
     .      form='unformatted')
      read(11) sisuz_bouchon
      close(11)

c ### boucle sur tous les recepteurs
	do irec=nrec/2+1,nrec

	print *,'Recepteur #',irec

c--- sauvegarde data gnuplot
      write(namefile,60) irec
 60   format('dataxspec',i3.3,'.gnu')
      open(unit=17,file=namefile,status='unknown')
      do it=1,ntspec
! DK DK DK it - 1 change en it pour test
            write(17,*) (it)*dt_spec,sisuz_specfem(it,irec)*(-factor)
      enddo
      close(17)

c--- sauvegarde data gnuplot
      write(namefile,70) irec
 70   format('dataxbouchon',i3.3,'.gnu')
      open(unit=17,file=namefile,status='unknown')
      do it=1,ntbouchon
            write(17,*) (it-1)*dt_bouchon,sisuz_bouchon(it,irec)*factor
     .		* 4400. * 1.025e9 *2.
      enddo
      close(17)

	enddo

c petit script de visu pour gnuplot
      open(unit=17,file='plotall_x.gnu',status='unknown')
	write(17,*) 'set term x11'
	write(17,*) 'set title "Horizontal component"'
	write(17,*) 'set xrange [0:35]'
      do irec=nrec/2+1,nrec
      	write(17,100) irec,irec
		write(17,*) 'pause -1 "Hit any key..."'
	enddo
      close(17)

 100	format('plot "dataxspec',i3.3,'.gnu" w l 1, "dataxbouchon',i3.3
     .		,'.gnu" w l 2')

	write(*,*) 'Termine !!!'

      END
