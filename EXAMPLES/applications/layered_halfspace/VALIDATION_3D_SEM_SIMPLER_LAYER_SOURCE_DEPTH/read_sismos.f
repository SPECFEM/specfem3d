      program lit_sismos
c
c=======================================================================
c
c lecture sismos en sortie de specfem 3D
c      -----------------
c
c ======================================================================
c

      include 'common_grid.h'
      include 'common_solve.h'

      integer ispec,idummy,i

! DK DK DK tester sismogrammes
      integer nreceptsmax
      parameter(nreceptsmax = 200)
      integer nreceptsX,nreceptsY
      real yreceivX(nreceptsmax),xreceivY(nreceptsmax)
      real xdummy1,xdummy2,xdummy3
      real sisuxX(nseis,nreceptsmax)
      real sisuyX(nseis,nreceptsmax)
      real sisuzX(nseis,nreceptsmax)
      real sisuxY(nseis,nreceptsmax)
      real sisuyY(nseis,nreceptsmax)
      real sisuzY(nseis,nreceptsmax)

      print *
      print *,'***************************************'
      print *,'**** lecture sismos Specfem 3D     ****'
      print *,'***************************************'
      print *
      call system('date')
      print *

      print *
      print *,' Time step = ',dtinc
      print *,' Number of time step = ',ncycl
      print *
      print *,' Nb echantillons pour seismogrammes = ',nseis
      print *,' Sous echantillonnage seismogrammes = ',isamp
      print *

      write(*,*)' Lecture des positions des recepteurs...'

c plan de coupe X = cste
      open(unit=20,file='receiversXcste.dat',status='old')
      read(20,*) nreceptsX
      if(nreceptsX .gt. nreceptsmax) stop 'receivers line too big'
      do ispec=1,nreceptsX
       read(20,*) idummy,xdummy1,yreceivX(ispec),xdummy3
       write(*,*) 'Recepteur Xcste numero ',ispec,' en Y = ',yreceivX(ispec)
      enddo
      close(20)

c plan de coupe Y = cste
      open(unit=20,file='receiversYcste.dat',status='old')
      read(20,*) nreceptsY
      if(nreceptsY .gt. nreceptsmax) stop 'receivers line too big'
      do ispec=1,nreceptsY
       read(20,*) idummy,xreceivY(ispec),xdummy2,xdummy3
       write(*,*) 'Recepteur Ycste numero ',ispec,' en X = ',xreceivY(ispec)
      enddo
      close(20)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      write(*,*) 'Lecture des sismos...'
      open(unit=27,file='sismos.bin',status='old',
     .            form='unformatted')
      read(27) sisuxX
      read(27) sisuyX
      read(27) sisuzX
      read(27) sisuxY
      read(27) sisuyY
      read(27) sisuzY
      close(27)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!----

      write(*,*)
      write(*,*) 'Sauvegarde traces format SEP...'

! ecriture au format binaire deplacement ligne X=cste
      open(unit=11,file='Ux_file_Xcste',status='unknown',
     .      access='direct',recl=nseis*nreceptsmax*4)
      write(11,rec=1) (sisuxX(i,1),i=1,nseis*nreceptsmax)
      close(11)
      open(unit=11,file='Uy_file_Xcste',status='unknown',
     .      access='direct',recl=nseis*nreceptsmax*4)
      write(11,rec=1) (sisuyX(i,1),i=1,nseis*nreceptsmax)
      close(11)
      open(unit=11,file='Uz_file_Xcste',status='unknown',
     .      access='direct',recl=nseis*nreceptsmax*4)
      write(11,rec=1) (sisuzX(i,1),i=1,nseis*nreceptsmax)
      close(11)

! ecriture au format binaire deplacement ligne Y=cste
      open(unit=11,file='Ux_file_Ycste',status='unknown',
     .      access='direct',recl=nseis*nreceptsmax*4)
      write(11,rec=1) (sisuxY(i,1),i=1,nseis*nreceptsmax)
      close(11)
      open(unit=11,file='Uy_file_Ycste',status='unknown',
     .      access='direct',recl=nseis*nreceptsmax*4)
      write(11,rec=1) (sisuyY(i,1),i=1,nseis*nreceptsmax)
      close(11)
      open(unit=11,file='Uz_file_Ycste',status='unknown',
     .      access='direct',recl=nseis*nreceptsmax*4)
      write(11,rec=1) (sisuzY(i,1),i=1,nseis*nreceptsmax)
      close(11)

      end

