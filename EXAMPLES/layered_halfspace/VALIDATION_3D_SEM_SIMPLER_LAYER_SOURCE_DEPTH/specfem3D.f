
      program specfem3D
c
c=======================================================================
c
c     "s p e c f e m 3 D" : 3-D spectral elements code
c      -----------------
c
c ======================================================================
c

c DK DK DK
c DK DK DK  march99  Portage en OpenMP pour Dec Alpha IPGP
c DK DK DK

      include 'common_grid.h'
      include 'common_solve.h'

      real zero,two
      parameter(zero=0., two=2.)

      real bigval
      integer bigvalint
      parameter(bigval = 1.e20, bigvalint = 1000000000)

      integer k,ix,iy,iz,ipoin,ispec,iglobnum
      integer iglobnum1,iglobnum2,iglobnum3

      integer it,itmod,indice,istep,ispecloc
      integer ipoinsource,ixgllsource,iygllsource,izgllsource,ispecsource

      integer inbarrays_npoin,inbarrays_nspec
      real taillemem

      real hdgll,ricker
      external hdgll,ricker

      real a1(nxgll,nxgll,nxgll)
      real a2(nxgll,nxgll,nxgll)
      real a3(nxgll,nxgll,nxgll)
      real a4(nxgll,nxgll,nxgll)
      real a5(nxgll,nxgll,nxgll)
      real a6(nxgll,nxgll,nxgll)
      real a7(nxgll,nxgll,nxgll)
      real a8(nxgll,nxgll,nxgll)
      real a9(nxgll,nxgll,nxgll)

! DK DK DK pour source explosive
      real a11(nxgll,nxgll,nxgll)
      real a12(nxgll,nxgll,nxgll)
      real a13(nxgll,nxgll,nxgll)
      real ricktime
! DK DK DK pour source explosive

! petits tableaux temporaires pour la parallelisation
      real
     . tempx1(nxgll,nxgll,nxgll),
     . tempx3(nxgll,nxgll,nxgll),
     . tempx5(nxgll,nxgll,nxgll),
     . tempy1(nxgll,nxgll,nxgll),
     . tempy3(nxgll,nxgll,nxgll),
     . tempy5(nxgll,nxgll,nxgll),
     . tempz1(nxgll,nxgll,nxgll),
     . tempz3(nxgll,nxgll,nxgll),
     . tempz5(nxgll,nxgll,nxgll),
     . Uxoldloc(nxgll,nxgll,nxgll),
     . Uyoldloc(nxgll,nxgll,nxgll),
     . Uzoldloc(nxgll,nxgll,nxgll)

      real xix(nxgll,nxgll,nxgll,nspec)
      real xiy(nxgll,nxgll,nxgll,nspec)
      real xiz(nxgll,nxgll,nxgll,nspec)
      real etax(nxgll,nxgll,nxgll,nspec)
      real etay(nxgll,nxgll,nxgll,nspec)
      real etaz(nxgll,nxgll,nxgll,nspec)
      real gammax(nxgll,nxgll,nxgll,nspec)
      real gammay(nxgll,nxgll,nxgll,nspec)
      real gammaz(nxgll,nxgll,nxgll,nspec)

      real rlambda(nxgll,nxgll,nxgll,nspec)
      real rmu(nxgll,nxgll,nxgll,nspec)

      integer ibool(nxgll,nxgll,nxgll,nspec)

      real displx(npoin)
      real disply(npoin)
      real displz(npoin)
      real velocx(npoin)
      real velocy(npoin)
      real velocz(npoin)
      real accelx(npoin)
      real accely(npoin)
      real accelz(npoin)

      real rmass(npoin)

      real Uxnewloc,Uynewloc,Uznewloc

      real derxxl,deryxl,derzxl,derxyl,deryyl,derzyl,derxzl,deryzl,
     .   derzzl,xixl,etayl,gammaxl,gammayl,gammazl,dvolul,rlambdal,rmul,
     .   xizl,etazl,xiyl,etaxl,rkappal

      real sumderxyl,sumderxzl,sumderyzl
      real sumderyyzzl,sumderxxzzl,sumderxxyyl

      real a1loc,a2loc,a3loc,a4loc,a5loc,a6loc,a7loc,a8loc,a9loc

      real tempx1l,tempx2l,tempx3l
      real tempy1l,tempy2l,tempy3l
      real tempz1l,tempz2l,tempz3l

      real tempx4l,tempx6l
      real tempy4l,tempy6l
      real tempz4l,tempz6l

! DK DK DK DK pour interpolation dans les plans de coupe
      integer nxdisp,nzdisp
      parameter(nxdisp=600,nzdisp=600)
      integer indicedispX(nxdisp,nzdisp)
      integer indicedispY(nxdisp,nzdisp)
      integer indicedispZ(nxdisp,nzdisp)
      real deltaxdispX,deltazdispX
      real deltaxdispY,deltazdispY
      real deltaxdispZ,deltazdispZ
! DK DK DK DK pour interpolation dans les plans de coupe

	integer ispecabs

! DK DK DK bord absorbant du fond
      integer nspecabsfondmax
      parameter(nspecabsfondmax = (nx/4)*(ny/4))

      integer nspecabsfond

      integer numspecabsfond(nspecabsfondmax)
      real dampPfond(nspecabsfondmax,nxgll,nygll)
      real dampSfond(nspecabsfondmax,nxgll,nygll)
! DK DK DK bord absorbant du fond

! DK DK DK bord absorbant X=Xmin
      integer nspecabsxminmax
      parameter(nspecabsxminmax = 1200)

      integer nspecabsxmin

      integer numspecabsxmin(nspecabsxminmax)
      real dampPxmin(nspecabsxminmax,nygll,nzgll)
      real dampSxmin(nspecabsxminmax,nygll,nzgll)

! DK DK DK bord absorbant X=Xmax
      integer nspecabsxmaxmax
      parameter(nspecabsxmaxmax = 1200)

      integer nspecabsxmax

      integer numspecabsxmax(nspecabsxmaxmax)
      real dampPxmax(nspecabsxmaxmax,nygll,nzgll)
      real dampSxmax(nspecabsxmaxmax,nygll,nzgll)

! DK DK DK bord absorbant Y=Ymin
      integer nspecabsyminmax
      parameter(nspecabsyminmax = 1200)

      integer nspecabsymin

      integer numspecabsymin(nspecabsyminmax)
      real dampPymin(nspecabsyminmax,nxgll,nzgll)
      real dampSymin(nspecabsyminmax,nxgll,nzgll)

! DK DK DK bord absorbant Y=Ymax
      integer nspecabsymaxmax
      parameter(nspecabsymaxmax = 1200)

      integer nspecabsymax

      integer numspecabsymax(nspecabsymaxmax)
      real dampPymax(nspecabsymaxmax,nxgll,nzgll)
      real dampSymax(nspecabsymaxmax,nxgll,nzgll)

! DK DK DK tester plans de coupe AVS
      integer nspecsidemax
      parameter(nspecsidemax = 14000)

      integer ngnodside
      parameter(ngnodside=4)

      integer nspecsideX,nspecsideY,nspecsideZ
      integer idummy,icell,numat,inodloc,inode,nbdonnees
      real donnee1,donnee2,donnee3,donnee4,donnee5,donnee6,donnee7,donnee8
      character *60 namefile

      integer icutXibool(4,nspecsidemax)
      real cutXx(4,nspecsidemax)
      real cutXy(4,nspecsidemax)
      real cutXz(4,nspecsidemax)
      real modvitX(4,nspecsidemax)

      integer icutYibool(4,nspecsidemax)
      real cutYx(4,nspecsidemax)
      real cutYy(4,nspecsidemax)
      real cutYz(4,nspecsidemax)
      real modvitY(4,nspecsidemax)

      integer icutZibool(4,nspecsidemax)
      real cutZx(4,nspecsidemax)
      real cutZy(4,nspecsidemax)
      real cutZz(4,nspecsidemax)
      real modvitZ(4,nspecsidemax)

      real maxnormZ(4,nspecsidemax)
! DK DK DK tester plans de coupe AVS

! DK DK DK tester sismogrammes
      integer nreceptsmax
      parameter(nreceptsmax = 200)
      integer nreceptsX,nreceptsY
      integer iboolreceivX(nreceptsmax),iboolreceivY(nreceptsmax)
      real xdummy1,xdummy2,xdummy3
      real sisuxX(nseis,nreceptsmax)
      real sisuyX(nseis,nreceptsmax)
      real sisuzX(nseis,nreceptsmax)
      real sisuxY(nseis,nreceptsmax)
      real sisuyY(nseis,nreceptsmax)
      real sisuzY(nseis,nreceptsmax)
! DK DK DK tester sismogrammes

      real deltat,deltatover2,deltatsqover2

      real xixmin,etaymin,gammaxmin,gammaymin,gammazmin,xizmin,etazmin
      real xixmax,etaymax,gammaxmax,gammaymax,gammazmax,xizmax,etazmax
      real xiymin,etaxmin
      real xiymax,etaxmax
      real dvolumin,dvolumax
      real rlambdamin,rmumin,rlambdamax,rmumax
      real Unormmax

      integer iboolmin,iboolmax

      print *
      print *,'***************************************'
      print *,'**** Specfem Curvilinear 3D Solver ****'
      print *,'***************************************'
      print *
      call system('date')
      print *

      print *
      print *,'Sur Dec Alpha :'
      print *

      print *,' nspec = ',nspec
      print *,' ndime = ',ndime
      print *,' nxgll = ',nxgll
      print *
      print *,' npoin = ',npoin
      print *
      print *,' Time step = ',dtinc
      print *,' Number of time step = ',ncycl
      print *
      print *,' Frequence affichage = ',itaff
      print *
      print *,' Nb echantillons pour seismogrammes = ',nseis
      print *,' Sous echantillonnage seismogrammes = ',isamp
      print *

c affichage de la taille memoire des tableaux
      inbarrays_npoin = 10
      inbarrays_nspec = 12
      taillemem = dble(4*(inbarrays_npoin*npoin +
     .    inbarrays_nspec*nspec*nxgll*nxgll*nxgll))/dble(1024*1024)
      print *,' Nb de tableaux de taille npoin = ',inbarrays_npoin
      print *,' Nb de tableaux de taille nspec = ',inbarrays_nspec
      print *,' '
      print *,' Taille memoire totale des ',
     .      inbarrays_npoin+inbarrays_nspec,
     .      ' tableaux = ',taillemem,' Megabytes'
      print *,' '

c lecture metrique

! DK DK DK pour onde plane
!c---- zmesh
!      print *,'Lecture array zmesh pour onde plane'
!      open(unit=27,file='zmesh.bin',status='old',
!     .                  form='unformatted')
!      read(27) zmesh
!      close(27)
! DK DK DK pour onde plane

      print *,'debut lecture metrique'

c---- xix
      print *,'Lecture array xix'
      open(unit=27,file='xix.bin',status='old',form='unformatted')
      read(27) xix
      close(27)

c---- xiy
      print *,'Lecture array xiy'
      open(unit=27,file='xiy.bin',status='old',form='unformatted')
      read(27) xiy
      close(27)

c---- xiz
      print *,'Lecture array xiz'
      open(unit=27,file='xiz.bin',status='old',form='unformatted')
      read(27) xiz
      close(27)

c---- etax
      print *,'Lecture array etax'
      open(unit=27,file='etax.bin',status='old',form='unformatted')
      read(27) etax
      close(27)

c---- etay
      print *,'Lecture array etay'
      open(unit=27,file='etay.bin',status='old',form='unformatted')
      read(27) etay
      close(27)

c---- etaz
      print *,'Lecture array etaz'
      open(unit=27,file='etaz.bin',status='old',form='unformatted')
      read(27) etaz
      close(27)

c---- gammax
      print *,'Lecture array gammax'
      open(unit=27,file='gammax.bin',status='old',form='unformatted')
      read(27) gammax
      close(27)

c---- gammay
      print *,'Lecture array gammay'
      open(unit=27,file='gammay.bin',status='old',form='unformatted')
      read(27) gammay
      close(27)

c---- gammaz
      print *,'Lecture array gammaz'
      open(unit=27,file='gammaz.bin',status='old',form='unformatted')
      read(27) gammaz
      close(27)

      print *,'fin lecture metrique'

c---- calcul des mins et maxs
      xixmin = bigval
      xiymin = bigval
      xizmin = bigval
      etaxmin = bigval
      etaymin = bigval
      etazmin = bigval
      gammaxmin = bigval
      gammaymin = bigval
      gammazmin = bigval
      dvolumin  = bigval
      xixmax = - bigval
      xiymax = - bigval
      xizmax = - bigval
      etaxmax = - bigval
      etaymax = - bigval
      etazmax = - bigval
      gammaxmax = - bigval
      gammaymax = - bigval
      gammazmax = - bigval
      dvolumax  = - bigval
      do ispec=1,nspec
      do iz=1,nzgll
      do iy=1,nygll
      do ix=1,nxgll
            xixmin = amin1(xixmin,xix(ix,iy,iz,ispec))
            xiymin = amin1(xiymin,xiy(ix,iy,iz,ispec))
            xizmin = amin1(xizmin,xiz(ix,iy,iz,ispec))
            etaxmin = amin1(etaxmin,etax(ix,iy,iz,ispec))
            etaymin = amin1(etaymin,etay(ix,iy,iz,ispec))
            etazmin = amin1(etazmin,etaz(ix,iy,iz,ispec))
            gammaxmin = amin1(gammaxmin,gammax(ix,iy,iz,ispec))
            gammaymin = amin1(gammaymin,gammay(ix,iy,iz,ispec))
            gammazmin = amin1(gammazmin,gammaz(ix,iy,iz,ispec))
            xixmax = amax1(xixmax,xix(ix,iy,iz,ispec))
            xiymax = amax1(xiymax,xiy(ix,iy,iz,ispec))
            xizmax = amax1(xizmax,xiz(ix,iy,iz,ispec))
            etaxmax = amax1(etaxmax,etax(ix,iy,iz,ispec))
            etaymax = amax1(etaymax,etay(ix,iy,iz,ispec))
            etazmax = amax1(etazmax,etaz(ix,iy,iz,ispec))
            gammaxmax = amax1(gammaxmax,gammax(ix,iy,iz,ispec))
            gammaymax = amax1(gammaymax,gammay(ix,iy,iz,ispec))
            gammazmax = amax1(gammazmax,gammaz(ix,iy,iz,ispec))

c DK DK march 99 recalculer dvolu au lieu de le stocker
            xixl = xix(ix,iy,iz,ispec)
            xiyl = xiy(ix,iy,iz,ispec)
            xizl = xiz(ix,iy,iz,ispec)
            etaxl = etax(ix,iy,iz,ispec)
            etayl = etay(ix,iy,iz,ispec)
            etazl = etaz(ix,iy,iz,ispec)
            gammaxl = gammax(ix,iy,iz,ispec)
            gammayl = gammay(ix,iy,iz,ispec)
            gammazl = gammaz(ix,iy,iz,ispec)
            dvolul = 1./(xixl*(etayl*gammazl-etazl*gammayl)
     .               - xiyl*(etaxl*gammazl-etazl*gammaxl)
     .               + xizl*(etaxl*gammayl-etayl*gammaxl))

            dvolumin = amin1(dvolumin,dvolul)
            dvolumax = amax1(dvolumax,dvolul)
      enddo
      enddo
      enddo
      enddo

      print *,' Metrique :'
      print *
      print *,' xix = ',xixmin,xixmax
      print *,' xiy = ',xiymin,xiymax
      print *,' xiz = ',xizmin,xizmax
      print *,' etax = ',etaxmin,etaxmax
      print *,' etay = ',etaymin,etaymax
      print *,' etaz = ',etazmin,etazmax
      print *,' gammax = ',gammaxmin,gammaxmax
      print *,' gammay = ',gammaymin,gammaymax
      print *,' gammaz = ',gammazmin,gammazmax
      print *
      print *,' dvolu = ',dvolumin,dvolumax
      print *

c lecture modele de vitesse

      print *,'debut lecture modele de vitesse'

c---- rlambda
      print *,'Lecture array rlambda'
      open(unit=27,file='rlambda.bin',status='old',form='unformatted')
      read(27) rlambda
      close(27)

c---- rmu
      print *,'Lecture array rmu'
      open(unit=27,file='rmu.bin',status='old',form='unformatted')
      read(27) rmu
      close(27)

      print *,'fin lecture modele de vitesse'

c---- calcul des mins et maxs
      rlambdamin = bigval
      rmumin = bigval
      rlambdamax = - bigval
      rmumax = - bigval
      do ispec=1,nspec
      do iz=1,nzgll
      do iy=1,nygll
      do ix=1,nxgll
            rlambdamin = amin1(rlambdamin,rlambda(ix,iy,iz,ispec))
            rmumin = amin1(rmumin,rmu(ix,iy,iz,ispec))
            rlambdamax = amax1(rlambdamax,rlambda(ix,iy,iz,ispec))
            rmumax = amax1(rmumax,rmu(ix,iy,iz,ispec))
      enddo
      enddo
      enddo
      enddo

      print *,' Modele de vitesse :'
      print *,' '
      print *,' rlambda = ',rlambdamin,rlambdamax
      print *,' rmu = ',rmumin,rmumax
      print *

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c lecture a1 a a9

      print *,'debut lecture tableaux a1 a a9'

      open(unit=27,file='a1_a9.bin',status='old',
     .                  form='unformatted')
      read(27) a1
      read(27) a2
      read(27) a3
      read(27) a4
      read(27) a5
      read(27) a6
      read(27) a7
      read(27) a8
      read(27) a9
      close(27)

      print *,'fin lecture tableaux a1 a a9'

c lecture a11 a a13

      print *,'debut lecture tableaux a11 a a13'

      open(unit=27,file='a11_a13.bin',status='old',
     .                  form='unformatted')
      read(27) a11
      read(27) a12
      read(27) a13
      close(27)

      print *,'fin lecture tableaux a11 a a13'


c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c read the mass matrix

      print *
      print *,'debut lecture matrice de masse'

      open(unit=27,file='rmass.bin',status='old',
     .                  form='unformatted')
      read(27) rmass
      close(27)

      print *,'fin lecture matrice de masse'
      print *

c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c read the global numbering

      print *
      print *,'debut lecture tableau ibool'

      open(unit=27,file='ibool.bin',status='old',
     .                  form='unformatted')
      read(27) ibool
      close(27)

      print *,'fin lecture tableau ibool'

      iboolmin = + bigvalint
      iboolmax = - bigvalint
      do ispec=1,nspec
      do ix=1,nxgll
            do iy=1,nygll
                  do iz=1,nzgll
            iboolmin = min(iboolmin,ibool(ix,iy,iz,ispec))
            iboolmax = max(iboolmax,ibool(ix,iy,iz,ispec))
                  enddo
            enddo
      enddo
      enddo

      print *,'minval maxval ibool = ',iboolmin,iboolmax

      if(iboolmin .ne. 1 .or. iboolmax .ne. npoin)
     .      stop 'incorrect global numbering read'

      print *

c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      print *,'Creating gnuplot file for source time function'

      open(unit=11,file='sources',status='unknown')
      do it=1,ncycl
      write(11,*) it*dtinc,ricker(it*dtinc)
      enddo
      close(11)

c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c---- lecture de la position du point source
      open(unit=54,file='positsource.dat',status='old')
      read(54,*) ipoinsource
      read(54,*) ixgllsource,iygllsource,izgllsource,ispecsource
      close(54)

! lecture des indices des points des plans de coupe
      open(unit=54,file='pointscoupeX.dat',status='old')
      read(54,*) deltaxdispX,deltazdispX
      do ix=1,nxdisp
      do iz=1,nzdisp
            read(54,*) indicedispX(ix,iz)
      enddo
      enddo
      close(54)

! lecture des indices des points des plans de coupe
      open(unit=54,file='pointscoupeY.dat',status='old')
      read(54,*) deltaxdispY,deltazdispY
      do ix=1,nxdisp
      do iz=1,nzdisp
            read(54,*) indicedispY(ix,iz)
      enddo
      enddo
      close(54)

! lecture des indices des points des plans de coupe
      open(unit=54,file='pointscoupeZ.dat',status='old')
      read(54,*) deltaxdispZ,deltazdispZ
      do ix=1,nxdisp
      do iz=1,nzdisp
            read(54,*) indicedispZ(ix,iz)
      enddo
      enddo
      close(54)

c pour la rapidite precalculer l'inverse de la matrice de masse
!$omp parallel do
      do ipoin=1,npoin
            rmass(ipoin) = 1.0 / rmass(ipoin)
      enddo

c initialiser les tableaux a zero
!$omp parallel do
      do ipoin=1,npoin
            displx(ipoin) = zero
            disply(ipoin) = zero
            displz(ipoin) = zero
            velocx(ipoin) = zero
            velocy(ipoin) = zero
            velocz(ipoin) = zero
            accelx(ipoin) = zero
            accely(ipoin) = zero
            accelz(ipoin) = zero
      enddo

! DK DK DK tester enregistrement des sismogrammes

      write(*,*)' Lecture des positions des recepteurs...'

c plan de coupe X = cste
      open(unit=20,file='receiversXcste.dat',status='old')
      read(20,*) nreceptsX
      if(nreceptsX .gt. nreceptsmax) stop 'receivers line too big'
      do ispec=1,nreceptsX
       read(20,*) iboolreceivX(ispec),xdummy1,xdummy2,xdummy3
       write(*,*) 'Recepteur Xcste numero ',ispec,' en ibool = ',iboolreceivX(ispec)
      enddo
      close(20)

c plan de coupe Y = cste
      open(unit=20,file='receiversYcste.dat',status='old')
      read(20,*) nreceptsY
      if(nreceptsY .gt. nreceptsmax) stop 'receivers line too big'
      do ispec=1,nreceptsY
       read(20,*) iboolreceivY(ispec),xdummy1,xdummy2,xdummy3
       write(*,*) 'Recepteur Ycste numero ',ispec,' en ibool = ',iboolreceivY(ispec)
      enddo
      close(20)

! DK DK DK tester enregistrement des sismogrammes

! DK DK DK lecture des donnees pour bords absorbants

      if(bords_abs) then

      write(*,*) 'Lecture des donnees pour bord absorbant du fond...'

c--- bord absorbant du fond

      open(unit=20,file='bord_abs_fond.dat',status='old')
      read(20,*) nspecabsfond
      if(nspecabsfond .gt. nspecabsfondmax) stop 'bord abs fond too big'
      do ispecabs=1,nspecabsfond
            read(20,*) numspecabsfond(ispecabs)
            do iy=1,nygll
                  do ix=1,nxgll
                    read(20,*) dampPfond(ispecabs,ix,iy),dampSfond(ispecabs,ix,iy)
                  enddo
            enddo
      enddo

      write(*,*) 'Il y a ',nspecabsfond,' elements absorbants sur le fond'

c--- bord absorbant X=Xmin

      write(*,*) 'Lecture des donnees pour bord absorbant X=Xmin...'

      open(unit=20,file='bord_abs_X=Xmin.dat',status='old')
      read(20,*) nspecabsxmin
      if(nspecabsxmin .gt. nspecabsxminmax) stop 'bord abs X=Xmin big'
      do ispecabs=1,nspecabsxmin
            read(20,*) numspecabsxmin(ispecabs)
            do iy=1,nygll
                  do iz=1,nzgll
                    read(20,*) dampPxmin(ispecabs,iy,iz),dampSxmin(ispecabs,iy,iz)
                  enddo
            enddo
      enddo

      write(*,*) 'Il y a ',nspecabsxmin,' elements absorbants sur X=Xmin'

c--- bord absorbant X=Xmax

      write(*,*) 'Lecture des donnees pour bord absorbant X=Xmax...'

      open(unit=20,file='bord_abs_X=Xmax.dat',status='old')
      read(20,*) nspecabsxmax
      if(nspecabsxmax .gt. nspecabsxmaxmax) stop 'bord abs X=Xmax big'
      do ispecabs=1,nspecabsxmax
            read(20,*) numspecabsxmax(ispecabs)
            do iy=1,nygll
                  do iz=1,nzgll
                    read(20,*) dampPxmax(ispecabs,iy,iz),dampSxmax(ispecabs,iy,iz)
                  enddo
            enddo
      enddo

      write(*,*) 'Il y a ',nspecabsxmax,' elements absorbants sur X=Xmax'

c--- bord absorbant Y=Ymin

      write(*,*) 'Lecture des donnees pour bord absorbant Y=Ymin...'

      open(unit=20,file='bord_abs_Y=Ymin.dat',status='old')
      read(20,*) nspecabsymin
      if(nspecabsymin .gt. nspecabsyminmax) stop 'bord abs Y=Ymin big'
      do ispecabs=1,nspecabsymin
            read(20,*) numspecabsymin(ispecabs)
            do ix=1,nxgll
                  do iz=1,nzgll
                    read(20,*) dampPymin(ispecabs,ix,iz),dampSymin(ispecabs,ix,iz)
                  enddo
            enddo
      enddo

      write(*,*) 'Il y a ',nspecabsymin,' elements absorbants sur Y=Ymin'

c--- bord absorbant Y=Ymax

      write(*,*) 'Lecture des donnees pour bord absorbant Y=Ymax...'

      open(unit=20,file='bord_abs_Y=Ymax.dat',status='old')
      read(20,*) nspecabsymax
      if(nspecabsymax .gt. nspecabsymaxmax) stop 'bord abs Y=Ymax big'
      do ispecabs=1,nspecabsymax
            read(20,*) numspecabsymax(ispecabs)
            do ix=1,nxgll
                  do iz=1,nzgll
                    read(20,*) dampPymax(ispecabs,ix,iz),dampSymax(ispecabs,ix,iz)
                  enddo
            enddo
      enddo

      write(*,*) 'Il y a ',nspecabsymax,' elements absorbants sur Y=Ymax'

      endif

! DK DK DK lecture des donnees pour bords absorbants

! DK DK DK tester plans de coupe AVS

      write(*,*) 'Lecture des donnees pour plans de coupe format AVS...'

c plan de coupe X = cste
      open(unit=20,file='AVScutplaneXcste.dat',status='old')
      read(20,*) nspecsideX
      if(nspecsideX .gt. nspecsidemax) stop 'cut plane too big'
      do ispec=1,nspecsideX
       read(20,*) idummy,
     .  icutXibool(1,ispec),icutXibool(2,ispec),
     .  icutXibool(3,ispec),icutXibool(4,ispec),
     .  cutXx(1,ispec),cutXx(2,ispec),
     .  cutXx(3,ispec),cutXx(4,ispec),
     .  cutXy(1,ispec),cutXy(2,ispec),
     .  cutXy(3,ispec),cutXy(4,ispec),
     .  cutXz(1,ispec),cutXz(2,ispec),
     .  cutXz(3,ispec),cutXz(4,ispec),
     .  modvitX(1,ispec),modvitX(2,ispec),
     .  modvitX(3,ispec),modvitX(4,ispec)
      enddo
      close(20)

c plan de coupe Y = cste
      open(unit=20,file='AVScutplaneYcste.dat',status='old')
      read(20,*) nspecsideY
      if(nspecsideY .gt. nspecsidemax) stop 'cut plane too big'
      do ispec=1,nspecsideY
       read(20,*) idummy,
     .  icutYibool(1,ispec),icutYibool(2,ispec),
     .  icutYibool(3,ispec),icutYibool(4,ispec),
     .  cutYx(1,ispec),cutYx(2,ispec),
     .  cutYx(3,ispec),cutYx(4,ispec),
     .  cutYy(1,ispec),cutYy(2,ispec),
     .  cutYy(3,ispec),cutYy(4,ispec),
     .  cutYz(1,ispec),cutYz(2,ispec),
     .  cutYz(3,ispec),cutYz(4,ispec),
     .  modvitY(1,ispec),modvitY(2,ispec),
     .  modvitY(3,ispec),modvitY(4,ispec)
      enddo
      close(20)

c plan de coupe Z = ztopo (surface)
      open(unit=20,file='AVScutplaneZcste.dat',status='old')
      read(20,*) nspecsideZ
      if(nspecsideZ .gt. nspecsidemax) stop 'cut plane too big'
      do ispec=1,nspecsideZ
       read(20,*) idummy,
     .  icutZibool(1,ispec),icutZibool(2,ispec),
     .  icutZibool(3,ispec),icutZibool(4,ispec),
     .  cutZx(1,ispec),cutZx(2,ispec),
     .  cutZx(3,ispec),cutZx(4,ispec),
     .  cutZy(1,ispec),cutZy(2,ispec),
     .  cutZy(3,ispec),cutZy(4,ispec),
     .  cutZz(1,ispec),cutZz(2,ispec),
     .  cutZz(3,ispec),cutZz(4,ispec),
     .  modvitZ(1,ispec),modvitZ(2,ispec),
     .  modvitZ(3,ispec),modvitZ(4,ispec)
      enddo
      close(20)

c initialiser le stockage du max de la norme en surface
      do ispec=1,nspecsideZ
      do inodloc = 1,ngnodside
            maxnormZ(inodloc,ispec) = - 1.
      enddo
      enddo

! DK DK DK tester plans de coupe AVS

c definition de constantes pour le schema en temps
      deltat = dtinc
      deltatover2 = dtinc / 2.0
      deltatsqover2 = dtinc * dtinc / 2.0

c
c----          s t a r t   t i m e   i t e r a t i o n s
c

      itmod = 10

      print *
      print *,'Starting time iteration loop...'
      print *

      call system('date')
      call system('rm -f startdate.txt ; date > startdate.txt')

c
c initialiser les predictors au depart
c
!$omp parallel do
      do ipoin=1,npoin
            displx(ipoin) = displx(ipoin) + deltat * velocx(ipoin) +
     .            deltatsqover2 * accelx(ipoin)
            disply(ipoin) = disply(ipoin) + deltat * velocy(ipoin) +
     .            deltatsqover2 * accely(ipoin)
            displz(ipoin) = displz(ipoin) + deltat * velocz(ipoin) +
     .            deltatsqover2 * accelz(ipoin)

            velocx(ipoin) = velocx(ipoin) + deltatover2 * accelx(ipoin)
            velocy(ipoin) = velocy(ipoin) + deltatover2 * accely(ipoin)
            velocz(ipoin) = velocz(ipoin) + deltatover2 * accelz(ipoin)

            accelx(ipoin) = zero
            accely(ipoin) = zero
            accelz(ipoin) = zero
      enddo

c
c grande boucle serial sur les pas de temps
c
      do it=1,ncycl

      write(iout,100) it,it*dtinc,ncycl

      if(mod(it,itmod).eq.0.or.it.eq.2) then
           Unormmax = - bigval
           do ipoin=1,npoin
               Unormmax = amax1(Unormmax,
     . sqrt(abs(displx(ipoin)**2 + disply(ipoin)**2 + displz(ipoin)**2)))
           enddo
           print *,'Max norme vecteur U = ',Unormmax
      endif

c--- tache de calcul pour les elements spectraux

c grande boucle parallele sur tous les elements spectraux
!$omp parallel do default(none)
!$omp& private(ispec,ix,iy,iz,tempx1l,tempx2l,tempx3l)
!$omp& private(iglobnum,iglobnum1,iglobnum2,iglobnum3)
!$omp& private(tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l)
!$omp& private(k,a1loc,a2loc,a3loc,xixl,xiyl,xizl,etaxl,etayl,etazl)
!$omp& private(gammaxl,gammayl,gammazl,derxxl,deryxl,derzxl)
!$omp& private(derxyl,deryyl,derzyl,derxzl,deryzl,derzzl)
!$omp& private(dvolul,rlambdal,rmul,rkappal)
!$omp& private(tempx4l,tempx6l,tempy4l,tempy6l,tempz4l,tempz6l)
!$omp& private(a4loc,a5loc,a6loc,a7loc,a8loc,a9loc,ispecloc)
!$omp& private(Uxnewloc,Uynewloc,Uznewloc),
!$omp& private(sumderxyl,sumderxzl,sumderyzl)
!$omp& private(sumderyyzzl,sumderxxzzl,sumderxxyyl)
!$omp& shared(a1,a2,a3,a4,a5,a6,a7,a8,a9),
!$omp& shared(gammax,gammay,gammaz,xix,xiy,xiz,etax,etay,etaz),
!$omp& shared(displx,disply,displz),
!$omp& shared(velocx,velocy,velocz)
!$omp& shared(accelx,accely,accelz)
!$omp& shared(rlambda,rmu,ibool)
!$omp& shared(nspecabsfond,numspecabsfond,dampSfond,dampPfond)
!$omp& shared(nspecabsxmin,numspecabsxmin,dampSxmin,dampPxmin)
!$omp& shared(nspecabsxmax,numspecabsxmax,dampSxmax,dampPxmax)
!$omp& shared(nspecabsymin,numspecabsymin,dampSymin,dampPymin)
!$omp& shared(nspecabsymax,numspecabsymax,dampSymax,dampPymax)
!$omp& private(tempx1,tempx3,tempx5,tempy1,tempy3,tempy5)
!$omp& private(tempz1,tempz3,tempz5,Uxoldloc,Uyoldloc,Uzoldloc)
      do ispec=1,nspec

c faire un "get" pour recuperer les donnees globales en local
      do iz=1,nzgll
      do iy=1,nygll
      do ix=1,nxgll
            iglobnum = ibool(ix,iy,iz,ispec)
            Uxoldloc(ix,iy,iz) = displx(iglobnum)
            Uyoldloc(ix,iy,iz) = disply(iglobnum)
            Uzoldloc(ix,iy,iz) = displz(iglobnum)
      enddo
      enddo
      enddo

      do iz=1,nxgll
      do iy=1,nxgll
      do ix=1,nxgll

c Produits de matrice

            tempx1l = zero
            tempx2l = zero
            tempx3l = zero

            tempy1l = zero
            tempy2l = zero
            tempy3l = zero

            tempz1l = zero
            tempz2l = zero
            tempz3l = zero

            do k=1,nxgll

                  a1loc = a1(ix,k,iz)
                  a2loc = a2(k,iy,iz)
                  a3loc = a3(ix,k,iz)

                  tempx1l = tempx1l + a1loc*Uxoldloc(k,iy,iz)
                  tempx2l = tempx2l + Uxoldloc(ix,k,iz)*a2loc
                  tempx3l = tempx3l + Uxoldloc(ix,iy,k)*a3loc

                  tempy1l = tempy1l + a1loc*Uyoldloc(k,iy,iz)
                  tempy2l = tempy2l + Uyoldloc(ix,k,iz)*a2loc
                  tempy3l = tempy3l + Uyoldloc(ix,iy,k)*a3loc

                  tempz1l = tempz1l + a1loc*Uzoldloc(k,iy,iz)
                  tempz2l = tempz2l + Uzoldloc(ix,k,iz)*a2loc
                  tempz3l = tempz3l + Uzoldloc(ix,iy,k)*a3loc

            enddo

c obtention des derivees locales dans le domaine physique

      xixl     =     xix(ix,iy,iz,ispec)
      xiyl     =     xiy(ix,iy,iz,ispec)
      xizl     =     xiz(ix,iy,iz,ispec)
      etaxl    =    etax(ix,iy,iz,ispec)
      etayl    =    etay(ix,iy,iz,ispec)
      etazl    =    etaz(ix,iy,iz,ispec)
      gammaxl  =  gammax(ix,iy,iz,ispec)
      gammayl  =  gammay(ix,iy,iz,ispec)
      gammazl  =  gammaz(ix,iy,iz,ispec)

      derxxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
      deryxl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
      derzxl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

      derxyl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
      deryyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
      derzyl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

      derxzl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
      deryzl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
      derzzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

c DK DK march 99 recalculer dvolu au lieu de le stocker
            dvolul = 1./(xixl*(etayl*gammazl-etazl*gammayl)
     .               - xiyl*(etaxl*gammazl-etazl*gammaxl)
     .               + xizl*(etaxl*gammayl-etayl*gammaxl))

c obtention des contributions pour chaque terme

      rlambdal = rlambda(ix,iy,iz,ispec)
      rmul     =     rmu(ix,iy,iz,ispec)
      rkappal  = rlambdal + two*rmul

      sumderxyl = deryxl + derxyl
      sumderxzl = derxzl + derzxl
      sumderyzl = deryzl + derzyl
      sumderyyzzl = deryyl + derzzl
      sumderxxzzl = derxxl + derzzl
      sumderxxyyl = derxxl + deryyl

c---- terme en Ux

      tempx1(ix,iy,iz) = dvolul*
     . (rkappal*xixl*derxxl + rlambdal*xixl*sumderyyzzl +
     . rmul*xiyl*sumderxyl + rmul*xizl*sumderxzl)
      tempx3(ix,iy,iz) = dvolul*
     . (rkappal*etaxl*derxxl + rlambdal*etaxl*sumderyyzzl +
     . rmul*etayl*sumderxyl + rmul*etazl*sumderxzl)
      tempx5(ix,iy,iz) = dvolul*
     . (rkappal*gammaxl*derxxl + rlambdal*gammaxl*sumderyyzzl +
     . rmul*gammayl*sumderxyl + rmul*gammazl*sumderxzl)


c---- terme en Uy

      tempy1(ix,iy,iz) = dvolul*
     . ( rkappal*xiyl*deryyl +  rlambdal*xiyl*sumderxxzzl +
     .  rmul*xixl*sumderxyl + rmul*xizl*sumderyzl)
      tempy3(ix,iy,iz) = dvolul*
     . (rkappal*etayl*deryyl + rlambdal*etayl*sumderxxzzl +
     . rmul*etaxl*sumderxyl + rmul*etazl*sumderyzl)
      tempy5(ix,iy,iz) = dvolul*
     . (rkappal*gammayl*deryyl + rlambdal*gammayl*sumderxxzzl +
     . rmul*gammaxl*sumderxyl + rmul*gammazl*sumderyzl)


c---- terme en Uz

      tempz1(ix,iy,iz) = dvolul*
     . (rkappal*xizl*derzzl + rlambdal*xizl*sumderxxyyl +
     . rmul*xixl*sumderxzl + rmul*xiyl*sumderyzl)
      tempz3(ix,iy,iz) = dvolul*
     . (rkappal*etazl*derzzl + rlambdal*etazl*sumderxxyyl +
     . rmul*etaxl*sumderxzl + rmul*etayl*sumderyzl)
      tempz5(ix,iy,iz) = dvolul*
     . (rkappal*gammazl*derzzl + rlambdal*gammazl*sumderxxyyl +
     . rmul*gammaxl*sumderxzl + rmul*gammayl*sumderyzl)

      enddo
      enddo
      enddo

c Produits de matrice

      do iz=1,nxgll
      do iy=1,nxgll
      do ix=1,nxgll

c---- terme en Ux
            tempx2l = zero
            tempx4l = zero
            tempx6l = zero

c---- terme en Uy
            tempy2l = zero
            tempy4l = zero
            tempy6l = zero

c---- terme en Uz
            tempz2l = zero
            tempz4l = zero
            tempz6l = zero

            do k=1,nxgll

                  a4loc = a4(ix,k,iz)
                  a5loc = a5(k,iy,iz)
                  a6loc = a6(ix,k,iz)

                  tempx2l = tempx2l + a4loc*tempx1(k,iy,iz)
                  tempx4l = tempx4l + tempx3(ix,k,iz)*a5loc
                  tempx6l = tempx6l + tempx5(ix,iy,k)*a6loc

                  tempy2l = tempy2l + a4loc*tempy1(k,iy,iz)
                  tempy4l = tempy4l + tempy3(ix,k,iz)*a5loc
                  tempy6l = tempy6l + tempy5(ix,iy,k)*a6loc

                  tempz2l = tempz2l + a4loc*tempz1(k,iy,iz)
                  tempz4l = tempz4l + tempz3(ix,k,iz)*a5loc
                  tempz6l = tempz6l + tempz5(ix,iy,k)*a6loc

            enddo

c introduction des coefficients de l'integration numerique

      a7loc = a7(ix,iy,iz)
      a8loc = a8(ix,iy,iz)
      a9loc = a9(ix,iy,iz)

      Uxnewloc = a7loc*tempx2l + a8loc*tempx4l + a9loc*tempx6l
      Uynewloc = a7loc*tempy2l + a8loc*tempy4l + a9loc*tempy6l
      Uznewloc = a7loc*tempz2l + a8loc*tempz4l + a9loc*tempz6l

! assemblage final des contribution, avec dependance de donnees
! dependance levee par l'utilisation de la clause "atomic" de OpenMP

      iglobnum = ibool(ix,iy,iz,ispec)

!$omp atomic
      accelx(iglobnum) = accelx(iglobnum) + Uxnewloc
!$omp atomic
      accely(iglobnum) = accely(iglobnum) + Uynewloc
!$omp atomic
      accelz(iglobnum) = accelz(iglobnum) + Uznewloc

c---- DK DK DK DK bords abs en traction DK DK DK DK

c bord abs du fond
      if(bords_abs .and. ispec .le. nspecabsfond .and. iz .eq. 1) then
           ispecloc = numspecabsfond(ispec)
           iglobnum = ibool(ix,iy,iz,ispecloc)
!$omp atomic
           accelx(iglobnum) = accelx(iglobnum) +
     .            dampSfond(ispec,ix,iy)*velocx(iglobnum)
!$omp atomic
           accely(iglobnum) = accely(iglobnum) +
     .            dampSfond(ispec,ix,iy)*velocy(iglobnum)
!$omp atomic
           accelz(iglobnum) = accelz(iglobnum) +
     .            dampPfond(ispec,ix,iy)*velocz(iglobnum)
      endif

c bord abs X=Xmin
      if(bords_abs .and. ispec .le. nspecabsxmin .and. ix .eq. 1) then
           ispecloc = numspecabsxmin(ispec)
           iglobnum = ibool(ix,iy,iz,ispecloc)
!$omp atomic
           accelx(iglobnum) = accelx(iglobnum) +
     .            dampPxmin(ispec,iy,iz)*velocx(iglobnum)
!$omp atomic
           accely(iglobnum) = accely(iglobnum) +
     .            dampSxmin(ispec,iy,iz)*velocy(iglobnum)
!$omp atomic
           accelz(iglobnum) = accelz(iglobnum) +
     .            dampSxmin(ispec,iy,iz)*velocz(iglobnum)
      endif

c bord abs X=Xmax
      if(bords_abs .and. ispec .le. nspecabsxmax .and. ix .eq. nxgll) then
           ispecloc = numspecabsxmax(ispec)
           iglobnum = ibool(ix,iy,iz,ispecloc)
!$omp atomic
           accelx(iglobnum) = accelx(iglobnum) +
     .            dampPxmax(ispec,iy,iz)*velocx(iglobnum)
!$omp atomic
           accely(iglobnum) = accely(iglobnum) +
     .            dampSxmax(ispec,iy,iz)*velocy(iglobnum)
!$omp atomic
           accelz(iglobnum) = accelz(iglobnum) +
     .            dampSxmax(ispec,iy,iz)*velocz(iglobnum)
      endif

c bord abs Y=Ymin
      if(bords_abs .and. ispec .le. nspecabsymin .and. iy .eq. 1) then
           ispecloc = numspecabsymin(ispec)
           iglobnum = ibool(ix,iy,iz,ispecloc)
!$omp atomic
           accelx(iglobnum) = accelx(iglobnum) +
     .            dampSymin(ispec,ix,iz)*velocx(iglobnum)
!$omp atomic
           accely(iglobnum) = accely(iglobnum) +
     .            dampPymin(ispec,ix,iz)*velocy(iglobnum)
!$omp atomic
           accelz(iglobnum) = accelz(iglobnum) +
     .            dampSymin(ispec,ix,iz)*velocz(iglobnum)
      endif

c bord abs Y=Ymax
      if(bords_abs .and. ispec .le. nspecabsymax .and. iy .eq. nygll) then
           ispecloc = numspecabsymax(ispec)
           iglobnum = ibool(ix,iy,iz,ispecloc)
!$omp atomic
           accelx(iglobnum) = accelx(iglobnum) +
     .            dampSymax(ispec,ix,iz)*velocx(iglobnum)
!$omp atomic
           accely(iglobnum) = accely(iglobnum) +
     .            dampPymax(ispec,ix,iz)*velocy(iglobnum)
!$omp atomic
           accelz(iglobnum) = accelz(iglobnum) +
     .            dampSymax(ispec,ix,iz)*velocz(iglobnum)
      endif

c---- DK DK DK DK bords abs en traction DK DK DK DK

      enddo
      enddo
      enddo

! fin de la grosse boucle parallele sur tous les elements spectraux
      enddo

      ricktime = ricker(it*dtinc)

c--- DK DK DK ajout de la source explosive
!        do ix=1,nxgll
!        do iz=1,nzgll
!        do iy=1,nygll
!              iglobnum = ibool(ix,iy,iz,ispecsource)
!        accelx(iglobnum) = accelx(iglobnum) + a11(ix,iy,iz)*ricktime
!        accely(iglobnum) = accely(iglobnum) + a12(ix,iy,iz)*ricktime
!        accelz(iglobnum) = accelz(iglobnum) + a13(ix,iy,iz)*ricktime
!        enddo
!        enddo
!        enddo

c ajout du terme source en force
c force suivant Z pour le moment
c le signe est negatif pour compenser le signe moins de rmass
      accelz(ipoinsource) = accelz(ipoinsource) - factor_force * ricktime

!$omp parallel do
      do ipoin=1,npoin

c multiplication par l'inverse de la matrice de masse pour calcul acceleration
            accelx(ipoin) = - accelx(ipoin)*rmass(ipoin)
            accely(ipoin) = - accely(ipoin)*rmass(ipoin)
            accelz(ipoin) = - accelz(ipoin)*rmass(ipoin)

c mise a jour de la vitesse
            velocx(ipoin) = velocx(ipoin) + deltatover2 * accelx(ipoin)
            velocy(ipoin) = velocy(ipoin) + deltatover2 * accely(ipoin)
            velocz(ipoin) = velocz(ipoin) + deltatover2 * accelz(ipoin)

c DK DK predictors deplaces ici test
            displx(ipoin) = displx(ipoin) + deltat * velocx(ipoin) +
     .            deltatsqover2 * accelx(ipoin)
            disply(ipoin) = disply(ipoin) + deltat * velocy(ipoin) +
     .            deltatsqover2 * accely(ipoin)
            displz(ipoin) = displz(ipoin) + deltat * velocz(ipoin) +
     .            deltatsqover2 * accelz(ipoin)

            velocx(ipoin) = velocx(ipoin) + deltatover2 * accelx(ipoin)
            velocy(ipoin) = velocy(ipoin) + deltatover2 * accely(ipoin)
            velocz(ipoin) = velocz(ipoin) + deltatover2 * accelz(ipoin)

            accelx(ipoin) = zero
            accely(ipoin) = zero
            accelz(ipoin) = zero
c DK DK predictors deplaces ici test

      enddo

c stockage sismogrammes
      if (mod(it,isamp).eq.0) then
            istep = it/isamp
            print *,'seismos it istep = ',it,istep

c ligne X=cste
      do ispec=1,nreceptsX
            sisuxX(istep,ispec) = displx(iboolreceivX(ispec))
            sisuyX(istep,ispec) = disply(iboolreceivX(ispec))
            sisuzX(istep,ispec) = displz(iboolreceivX(ispec))
      enddo

c ligne Y=cste
      do ispec=1,nreceptsY
            sisuxY(istep,ispec) = displx(iboolreceivY(ispec))
            sisuyY(istep,ispec) = disply(iboolreceivY(ispec))
            sisuzY(istep,ispec) = displz(iboolreceivY(ispec))
      enddo

c calcul du max de la norme en surface au cours du temps
c      do ispec=1,nspecsideZ
c      do inodloc = 1,ngnodside
c            ipoin   = icutZibool(inodloc,ispec)
c            donnee7 = sqrt(displx(ipoin)**2 + disply(ipoin)**2 + displz(ipoin)**2)
c            maxnormZ(inodloc,ispec) = amax1(maxnormZ(inodloc,ispec),donnee7)
c      enddo
c      enddo

      endif

c affichage des resultats
      if ((mod(it,itaff).eq.0).or.(it.eq.itfirstaff)) then

      write(iout,*)
      write(iout,110) it*dtinc
      write(iout,*)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c affichage postscript

      write(iout,*) 'Dump PostScript'

! coupe suivant X = cste
      indice = 1
      call plotpost(disply,displz,npoin,indicedispX,nxdisp,nzdisp,
     .            deltaxdispX,deltazdispX,it,dtinc,indice)
! coupe suivant Y = cste
      indice = 2
      call plotpost(displx,displz,npoin,indicedispY,nxdisp,nzdisp,
     .            deltaxdispY,deltazdispY,it,dtinc,indice)
! coupe suivant Z = cste
c DK DK march99    indice = 3
c DK DK march99    call plotpost(displx,disply,npoin,indicedispZ,nxdisp,nzdisp,
c DK DK march99   .            deltaxdispZ,deltazdispZ,it,dtinc,indice)

      write(iout,*) 'Fin dump PostScript'

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      write(*,*) 'Sauvegarde des sismos...'

      open(unit=27,file='sismos.bin',status='unknown',
     .            form='unformatted')
      write(27) sisuxX
      write(27) sisuyX
      write(27) sisuzX
      write(27) sisuxY
      write(27) sisuyY
      write(27) sisuzY
      close(27)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      write(*,*) 'Sauvegarde de l''amplitude max en surface...'
      open(unit=27,file='ampli_max.txt',status='unknown')
      write(27,*) nspecsideZ,ngnodside
      do ispec=1,nspecsideZ
      do inodloc = 1,ngnodside
            write(27,*) cutZx(inodloc,ispec),cutZy(inodloc,ispec),
     .                        cutZz(inodloc,ispec),maxnormZ(inodloc,ispec)
      enddo
      enddo
      close(27)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c affichage AVS

      goto 743

      write(*,*)' Ecriture des plans de coupe format AVS...'

c coupe suivant X

      write(namefile,221) it
 221  format('AVScutplaneX',i5.5,'.inp')

      open(unit=20,file=namefile,status='unknown')

c nb de noeuds, de cellules, de donnees par cellule
      nbdonnees = 8
      write(20,*) 4*nspecsideX,nspecsideX,nbdonnees,0,0

c numero et coordonnees des noeuds
      inode = 0
      do ispec=1,nspecsideX
      do inodloc = 1,ngnodside
      inode = inode + 1
      write(20,200) inode,cutXx(inodloc,ispec),cutXy(inodloc,ispec),
     .                  cutXz(inodloc,ispec)
      enddo
      enddo

c numero et coordonnees des cellules, definition materiau
      icell = 0
      numat = 1
      do ispec=1,nspecsideX
          write(20,310) ispec,numat,icell+1,icell+2,icell+3,icell+4
          icell = icell + 4
      enddo

c structure des donnees aux noeuds
c mettre bon nombre de 1 et de labels et modifier la taille du format 210 aussi
      write(20,*) nbdonnees,' 1 1 1 1 1 1 1 1'
      write(20,*) 'Xcoord, meters'
      write(20,*) 'Ycoord, meters'
      write(20,*) 'Zcoord, meters'
      write(20,*) 'Ux, meters'
      write(20,*) 'Uy, meters'
      write(20,*) 'Uz, meters'
      write(20,*) 'Unorm, meters'
      write(20,*) 'lambda, SI'

c donnees aux noeuds
      inode = 0
      do ispec=1,nspecsideX
      do inodloc = 1,ngnodside
            inode   = inode + 1
            ipoin   = icutXibool(inodloc,ispec)
            donnee1 = cutXx(inodloc,ispec)
            donnee2 = cutXy(inodloc,ispec)
            donnee3 = cutXz(inodloc,ispec)
            donnee4 = displx(ipoin)
            donnee5 = disply(ipoin)
            donnee6 = displz(ipoin)
            donnee7 = sqrt(displx(ipoin)**2 + disply(ipoin)**2 + displz(ipoin)**2)
            donnee8 = modvitX(inodloc,ispec)
            write(20,210) inode,donnee1,donnee2,donnee3,donnee4,
     .                        donnee5,donnee6,donnee7,donnee8
      enddo
      enddo

      close(20)

c coupe suivant Y

      write(namefile,222) it
 222  format('AVScutplaneY',i5.5,'.inp')

      open(unit=20,file=namefile,status='unknown')

c nb de noeuds, de cellules, de donnees par cellule
      nbdonnees = 8
      write(20,*) 4*nspecsideY,nspecsideY,nbdonnees,0,0

c numero et coordonnees des noeuds
      inode = 0
      do ispec=1,nspecsideY
      do inodloc = 1,ngnodside
      inode = inode + 1
      write(20,200) inode,cutYx(inodloc,ispec),cutYy(inodloc,ispec),
     .                  cutYz(inodloc,ispec)
      enddo
      enddo

c numero et coordonnees des cellules, definition materiau
      icell = 0
      numat = 1
      do ispec=1,nspecsideY
          write(20,310) ispec,numat,icell+1,icell+2,icell+3,icell+4
          icell = icell + 4
      enddo

c structure des donnees aux noeuds
c mettre bon nombre de 1 et de labels et modifier la taille du format 210 aussi
      write(20,*) nbdonnees,' 1 1 1 1 1 1 1 1'
      write(20,*) 'Xcoord, meters'
      write(20,*) 'Ycoord, meters'
      write(20,*) 'Zcoord, meters'
      write(20,*) 'Ux, meters'
      write(20,*) 'Uy, meters'
      write(20,*) 'Uz, meters'
      write(20,*) 'Unorm, meters'
      write(20,*) 'lambda, SI'

c donnees aux noeuds
      inode = 0
      do ispec=1,nspecsideY
      do inodloc = 1,ngnodside
            inode   = inode + 1
            ipoin   = icutYibool(inodloc,ispec)
            donnee1 = cutYx(inodloc,ispec)
            donnee2 = cutYy(inodloc,ispec)
            donnee3 = cutYz(inodloc,ispec)
            donnee4 = displx(ipoin)
            donnee5 = disply(ipoin)
            donnee6 = displz(ipoin)
            donnee7 = sqrt(displx(ipoin)**2 + disply(ipoin)**2 + displz(ipoin)**2)
            donnee8 = modvitY(inodloc,ispec)
            write(20,210) inode,donnee1,donnee2,donnee3,donnee4,
     .                        donnee5,donnee6,donnee7,donnee8
      enddo
      enddo

      close(20)

c coupe suivant Z

      write(namefile,223) it
 223  format('AVScutplaneZ',i5.5,'.inp')

      open(unit=20,file=namefile,status='unknown')

c nb de noeuds, de cellules, de donnees par cellule
      nbdonnees = 8
      write(20,*) 4*nspecsideZ,nspecsideZ,nbdonnees,0,0

c numero et coordonnees des noeuds
      inode = 0
      do ispec=1,nspecsideZ
      do inodloc = 1,ngnodside
      inode = inode + 1
      write(20,200) inode,cutZx(inodloc,ispec),cutZy(inodloc,ispec),
     .                  cutZz(inodloc,ispec)
      enddo
      enddo

c numero et coordonnees des cellules, definition materiau
      icell = 0
      numat = 1
      do ispec=1,nspecsideZ
          write(20,310) ispec,numat,icell+1,icell+2,icell+3,icell+4
          icell = icell + 4
      enddo

c structure des donnees aux noeuds
c mettre bon nombre de 1 et de labels et modifier la taille du format 210 aussi
      write(20,*) nbdonnees,' 1 1 1 1 1 1 1 1'
      write(20,*) 'Xcoord, meters'
      write(20,*) 'Ycoord, meters'
      write(20,*) 'Zcoord, meters'
      write(20,*) 'Ux, meters'
      write(20,*) 'Uy, meters'
      write(20,*) 'Uz, meters'
      write(20,*) 'Unorm, meters'
      write(20,*) 'lambda, SI'

c donnees aux noeuds
      inode = 0
      do ispec=1,nspecsideZ
      do inodloc = 1,ngnodside
            inode   = inode + 1
            ipoin   = icutZibool(inodloc,ispec)
            donnee1 = cutZx(inodloc,ispec)
            donnee2 = cutZy(inodloc,ispec)
            donnee3 = cutZz(inodloc,ispec)
            donnee4 = displx(ipoin)
            donnee5 = disply(ipoin)
            donnee6 = displz(ipoin)
            donnee7 = sqrt(displx(ipoin)**2 + disply(ipoin)**2 + displz(ipoin)**2)
            donnee8 = modvitZ(inodloc,ispec)
            write(20,210) inode,donnee1,donnee2,donnee3,donnee4,
     .                        donnee5,donnee6,donnee7,donnee8
      enddo
      enddo

      close(20)

 200  format(i8,1x,e12.5,1x,e12.5,1x,e12.5)
 210  format(i8,8(1x,e12.5))
 310  format(i8,1x,i1,' quad ',i5,1x,i5,1x,i5,1x,i5)

      write(*,*)' Fin ecriture des plans de coupe format AVS...'

 743  continue

! DK DK DK tester plans de coupe AVS

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      endif

      enddo

c
c---- end of time iteration loop
c

      write(iout,*) 'Simulation terminee !'

      call system('date')
      call system('rm -f enddate.txt ; date > enddate.txt')

      close(iout)

      stop

 100  format('Iteration numero ',i5,'   t = ',f6.3,' s sur ',i5,' iterations')
 110  format('Sauvegarde deplacement temps t = ',f6.3,' s')

      end

c *******

c calcul du terme temporel de la source pour un Ricker
      real function ricker(t)

      include 'common_grid.h'

      double precision pi
      parameter(pi=3.141592653589793d0)

      real t
      double precision a,td,rickerdble

c Ricker
      a = pi*pi*f0*f0
      td = dble(t)
      rickerdble = - (1.d0-2.d0*a*(td-t0)*(td-t0))*dexp(-a*(td-t0)*(td-t0))

      ricker = sngl(rickerdble)

      return
      end

