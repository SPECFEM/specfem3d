
      program generate
c
c=======================================================================
c
c     "m e s h 3 D" : 3-D mesh generation code
c      -----------
c
c ======================================================================
c

      include 'common_grid.h'

      integer infotime(3)

      double precision gaussalpha,gaussbeta
      parameter(gaussalpha=0.d0,gaussbeta=0.d0)

      double precision zero,one,one_eight
      parameter(zero = 0.d0, one = 1.d0, one_eight = 0.125d0)

      double precision pi
      parameter(pi = 3.14159265d0)

      integer bigvalint
      double precision epsilon,bigval
      parameter(bigval = 1.d20, bigvalint = 1000000000)
      parameter(epsilon = 1.d-9)

! DK DK DK la tolerance est trop faible si nxgll > 12 environ
      real smallvaltol
      parameter(smallvaltol = 0.00001)
!      double precision smallvaltol
!      parameter(smallvaltol = 0.000000001d0)

! petite securite pour la lecture du modele de vitesse externe
      real small_velocity_jump
      parameter(small_velocity_jump = 10.)

c numero de materiau des trois couches et du modele regional 1D
      integer itop,imiddle,ibottom,i12_25,i25_moho,imoho
      parameter(itop=1,imiddle=2,ibottom=3,i12_25=4,i25_moho=5,imoho=6)

c profondeur des couches du modele regional 1D
      double precision depth12,depth25,depthmoho
      parameter(depth12=-12000.d0,depth25=-25000.d0,depthmoho=-30000.d0)

c nb de points pour les fonctions de forme
      integer ngnod
      parameter(ngnod=8)

! DK DK DK DK pour bord absorbant DK DK DK DK
      real rjacobloc,dampP,dampS
! DK DK DK DK pour bord absorbant DK DK DK DK

      integer i,j,k,ix,iy,iz,ip,ia,ispec,imat
      integer ihh,imm,iss

      integer iglobnum,inode,icell,inodloc

      double precision hdgll
      external hdgll

      double precision zsurf
      external zsurf

      real rand
      external rand

! DK DK DK DK pour interpolation dans les plans de coupe
      integer nxdisp,nzdisp
      parameter(nxdisp=600,nzdisp=600)
      integer indicedispX(nxdisp,nzdisp)
      integer indicedispY(nxdisp,nzdisp)
      integer indicedispZ(nxdisp,nzdisp)
      double precision deltaxdisp,deltazdisp,valplancoupe,valtoler,xcutpos
      double precision xpoint,ypoint,zpoint
      integer indcutx,indcutz
! DK DK DK DK pour interpolation dans les plans de coupe

c ****** DK DK DK calcul du jacobien 2D sur les bords ******
      double precision coorgx(4),coorgz(4)
      double precision dershape2D(2,4,nxgll,nxgll)
      integer ip1,ip2,in,l1,l2
      double precision xjac2_11,xjac2_21,xjac2_12,xjac2_22
      double precision s,t,sp,sm,tp,tm
      double precision quart
      parameter(quart = 0.25d0)
c ****** DK DK DK calcul du jacobien 2D sur les bords ******

      double precision rlambda,rmu,cp,cs,denst,Kvol,poiss
      double precision ra1,ra2,rb1,rb2,rc1,rc2
      double precision sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma

      double precision donnee,zmax
      double precision deltax,deltay,deltaz

      double precision xmesh,ymesh,zmesh
      double precision dvolu,gammax,gammay,gammaz

      integer ibool(nxgll,nxgll,nxgll,nspec)

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
      real xixl,xiyl,xizl,hprimeval,etaxl,etayl,etazl
      real gammaxl,gammayl,gammazl
      integer l,m
! DK DK DK pour source explosive

      real dvolustock(nxgll,nxgll,nxgll,nspec)
      real xixstock(nxgll,nxgll,nxgll,nspec)
      real xiystock(nxgll,nxgll,nxgll,nspec)
      real xizstock(nxgll,nxgll,nxgll,nspec)
      real etaxstock(nxgll,nxgll,nxgll,nspec)
      real etaystock(nxgll,nxgll,nxgll,nspec)
      real etazstock(nxgll,nxgll,nxgll,nspec)
      real gammaxstock(nxgll,nxgll,nxgll,nspec)
      real gammaystock(nxgll,nxgll,nxgll,nspec)
      real gammazstock(nxgll,nxgll,nxgll,nspec)
      real rlambdastock(nxgll,nxgll,nxgll,nspec)
      real rmustock(nxgll,nxgll,nxgll,nspec)
      real xmeshstock(nxgll,nxgll,nxgll,nspec)
      real ymeshstock(nxgll,nxgll,nxgll,nspec)
      real zmeshstock(nxgll,nxgll,nxgll,nspec)

      real rmass(npoin)

      integer numat(nspec)

c pour numerotation locale -> globale
! attention ici plantage en simple precision si nxgll > environ 12
! declarer en double precision dans ce cas
!      double precision tol(ntot),xp(ntot),yp(ntot),zp(ntot),work(ntot)
!      double precision xmaxval,xminval,ymaxval,yminval,zmaxval,zminval,xtol
      real tol(ntot),xp(ntot),yp(ntot),zp(ntot),work(ntot)
      real xmaxval,xminval,ymaxval,yminval,zmaxval,zminval,xtol

      integer loc(ntot),ind(ntot),ninseg(ntot),iglob(ntot)
      logical ifseg(ntot)
      integer ieoff,ilocnum,ie,nseg,ioff,iseg,ig,nglob
      integer iboolmin,iboolmax

      double precision xigll(nxgll)
      double precision yigll(nxgll)
      double precision zigll(nxgll)
      double precision wxgll(nxgll)
      double precision wygll(nxgll)
      double precision wzgll(nxgll)

      double precision hprime(nxgll,nxgll)

c stockage des interfaces "tapered"
      double precision xsimu(0:nx,0:ny)
      double precision ysimu(0:nx,0:ny)
      double precision ztoposimu(0:nx,0:ny)
      double precision z18simu(0:nx,0:ny)
      double precision z28simu(0:nx,0:ny)

c stockage de la grille curvi (x, y et z) pour le cube
      double precision x(0:nx,0:ny,0:nz)
      double precision y(0:nx,0:ny,0:nz)
      double precision z(0:nx,0:ny,0:nz)

c stockage des interfaces interpolees en entree
      double precision xinterp(nxinterf,nyinterf)
      double precision yinterp(nxinterf,nyinterf)
      double precision zinterp18(nxinterf,nyinterf)
      double precision zinterp28(nxinterf,nyinterf)
      double precision zinterpkazu(nxinterf,nyinterf)
      double precision zinterptopobathy(nxinterf,nyinterf)

c stockage du modele de vitesse en entree
      integer nbmat
      parameter(nbmat = 10)
      double precision vp(nbmat),vs(nbmat),rho(nbmat)
      integer nxveloc,nyveloc,nzveloc
      parameter(nxveloc=327,nyveloc=401,nzveloc=90)
      real velocmodel(nxveloc,nyveloc,nzveloc)
      real xmin_veloc,xmax_veloc,ymin_veloc,ymax_veloc,zmin_veloc
      parameter(xmin_veloc=417652.,xmax_veloc=580652.)
      parameter(ymin_veloc=978186.,ymax_veloc=1178186.)
      parameter(zmin_veloc=-4500.)
      real deltax_veloc,deltay_veloc,deltaz_veloc
      parameter(deltax_veloc=500.,deltay_veloc=500.,deltaz_veloc=50.)
      real xcoordpt,ycoordpt,zcoordpt
      integer i_veloc,j_veloc,k_veloc

c stockage des fcts de forme
      double precision shape(ngnod,nxgll,nygll,nzgll)
      double precision dershape(ndime,ngnod,nxgll,nygll,nzgll)

c stockage de la topologie de l'element de controle
      integer iaddx(ngnod)
      integer iaddy(ngnod)
      integer iaddz(ngnod)
      double precision xcell(ngnod)
      double precision ycell(ngnod)
      double precision zcell(ngnod)

c stockage des variables normalisees (xi, eta et gamma)
      double precision xig(0:nx)
      double precision etag(0:ny)
      double precision gammag1(0:nz)
      double precision gammag2(0:nz)
      double precision gammag3(0:nz)

      integer i_interf,j_interf
      integer iavsfile,igenlay
      integer inbarrays_nspec,inbarrays_npoin
      integer ipositsource,ispecsource,ixgllsource,iygllsource,izgllsource

      logical idouble
      double precision taillemem,tailledisk
      double precision dist,distmin

      double precision ztopo,z18,z28,zbot
      double precision delta_x,delta_y
      double precision xmin_interf,ymin_interf
      double precision delta_x_interf,delta_y_interf

      double precision dmaxz0,dmaxz1,dmaxz2,dmaxz3
      double precision dminz0,dminz1,dminz2,dminz3
      double precision delta_z0,delta_z1,delta_z2,delta_z3

      double precision
     .      xiymin,xiymax,xizmin,xizmax,etaxmin,etaxmax,etazmin,etazmax,
     .   xixmin,xixmax,etaymin,etaymax,gammaxmin,gammaxmax,
     .   gammaymin,gammaymax,gammazmin,gammazmax

      double precision
     .      yximin,yximax,zximin,zximax,xetamin,xetamax,zetamin,zetamax,
     .   xximin,xximax,yetamin,yetamax,xgammamin,xgammamax,
     .   ygammamin,ygammamax,zgammamin,zgammamax

      integer nspecside,ngnodside,numelemligne,nxd,nxf,nyd,nyf,nzd,nzf
      double precision xtoler

      print *,'**** Generate Specfem Arrays 3D ****'
      print *

      call itime(infotime)
      ihh = infotime(1)
      imm = infotime(2)
      iss = infotime(3)

      print *,'Current time is :',ihh,'h',imm,'m',iss,'s'
      print *

      print *
      print *,' nspec = ',nspec
      print *
      print *,' ndime = ',ndime
      print *,' nxgll = ',nxgll
      print *
      print *,' nx = ',nx
      print *,' ny = ',ny
      print *
      print *,' nz = ',nz
      print *
      print *,' nz1 = ',nz1
      print *,' nz2 = ',nz2
      print *,' nz3 = ',nz3
      print *
      print *,' npoin theorique = ',npoin

      print *
      print *,' nspec sans deraf = ',nspecsansderaf
      print *,' nspec avec deraf = ',nspec
      print *,' ratio sans/avec = ',dble(nspecsansderaf)/dble(nspec)
      print *

      print *
      print *,' nspec zone topo = ',nspechaut
      print *,' nspec zone bedrock = ',nspecbas
      print *,' nspec zone raccord 1 = ',nspecderaf1
      print *,' nspec zone raccord 2 = ',nspecderaf2
      print *

! verifier qu'on peut faire le deraffinement en profondeur
      if(mod(nx,8) .ne. 0)  stop 'nx doit etre un multiple de 8'
      if(mod(ny,8) .ne. 0)  stop 'ny doit etre un multiple de 8'
      if(mod(nz1,4) .ne. 0) stop 'nz1 doit etre un multiple de 4'
      if(nz1 .lt. 12) stop 'nz1 doit etre superieur ou egal a 12'

! verifier l'ordre polynomial
      if(nygll .ne. nxgll) stop 'ordre polynomial non uniforme'
      if(nzgll .ne. nxgll) stop 'ordre polynomial non uniforme'

c affichage de la taille memoire des tableaux
      inbarrays_nspec = 26
      inbarrays_npoin = 1
      idouble = .true.
      idouble = .false.
      if(idouble) inbarrays_nspec = inbarrays_nspec + 5
      taillemem = 4.d0*dble(inbarrays_npoin*npoin +
     .    inbarrays_nspec*nspec*nxgll*nxgll*nxgll)/dble(1024*1024)
      print *
      print *,' Memory size for the grid generation code :'
      print *,' ------------------------------------------'
      print *
      print *,' Nb de tableaux de taille npoin = ',inbarrays_npoin
      print *,' Nb de tableaux de taille nspec = ',inbarrays_nspec
      print *,' Taille memoire totale des ',
     .      inbarrays_npoin+inbarrays_nspec,
     .      ' tableaux = ',taillemem,' Megabytes'
      print *,' '

c affichage de l'espace disk necessaire
      inbarrays_nspec = 9
      inbarrays_npoin = 1
      tailledisk = 4.d0*dble(inbarrays_npoin*npoin +
     .    inbarrays_nspec*nspec*nxgll*nxgll*nxgll)/dble(1024*1024)
      print *,' Disk storage for the output of the grid generation code :'
      print *,' ---------------------------------------------------------'
      print *
      print *,' Nb de tableaux de taille npoin = ',inbarrays_npoin
      print *,' Nb de tableaux de taille nspec = ',inbarrays_nspec
      print *,' Espace disk total des ',
     .      inbarrays_npoin+inbarrays_nspec,
     .      ' tableaux = ',tailledisk,' Megabytes'
      print *
      print *

      call system('date')

c parametres elastiques : definition du modele de vitesse
c taken from Riki's report Table 3

c topo to layer 18
      vp(itop)  = 1800.
      vs(itop)  =  680.
! DK DK DK DK valeur de 0.41 du Poisson's ratio trop elevee numeriquement
! DK DK DK DK valeur de Vs tronquee a 1000. pour l'instant
      vs(itop)  = 1000.
      rho(itop) = 2000.

c layer 18 to layer 28
      vp(imiddle)  = 2800.
      vs(imiddle)  = 1500.
      rho(imiddle) = 2300.

c layer 18 to 12 km
      vp(ibottom)  = 5500.
      vs(ibottom)  = 3000.
      rho(ibottom) = 2500.

c 12 km to 25 km
      vp(i12_25)  = 6150.
      vs(i12_25)  = 3400.
      rho(i12_25) = 2700.

c 25 km to 30 km (Moho)
      vp(i25_moho)  = 6700.
      vs(i25_moho)  = 3700.
      rho(i25_moho) = 2900.

c below 30 km (Moho)
      vp(imoho)  = 7500.
      vs(imoho)  = 4300.
      rho(imoho) = 3200.

! DK DK DK essai sans modele 1D pour test DK DK DK

c DK DK mise d'un modele plus simple

      do imat=1,6
           vp(imat)  = 7500.
           vs(imat)  = 4300.
           rho(imat) = 3200.
      enddo

c topo to layer 18
      vp(itop)  = 2800.
      vs(itop)  = 1500.
      rho(itop) = 2300.

c layer 18 to layer 28
      vp(imiddle)  = 2800.
      vs(imiddle)  = 1500.
      rho(imiddle) = 2300.

! DK DK DK essai sans modele 1D pour test DK DK DK

      print *
      print *,' Parametres elastiques modele Riki :'
      print *

      do imat=1,6

      rmu   = rho(imat)*vs(imat)*vs(imat)
      rlambda  = rho(imat)*vp(imat)*vp(imat)- 2.e0*rmu
      Kvol  = rlambda + 2.d0*rmu/3.d0
      poiss = 0.5d0*(3.d0*Kvol-2.d0*rmu)/(3.d0*Kvol+rmu)
      if (poiss .lt. 0.d0 .or. poiss .ge. 0.50001d0)
     .            stop 'Poisson''s ratio out of range !'

      print *,' Layer ',imat
      print *,' cp = ',vp(imat)
      print *,' cs = ',vs(imat)
      print *,' rho = ',rho(imat)
      print *,' Poisson''s ratio = ',poiss
      print *,' lambda = ',rlambda
      print *,' mu = ',rmu
      print *

      enddo

      print *
      print *,' Regional 1D model :'
      print *
      print *,' Depth first step = ',depth12
      print *,' Depth second step = ',depth25
      print *,' Depth Moho = ',depthmoho
      print *

c
c----    set-up coordinates of the Gauss-Lobatto-Legendre points
c
      call zwgljd(xigll,wxgll,nxgll,gaussalpha,gaussbeta)

c---- if nb of points is odd, the middle abscissa is exactly zero
      if(mod(nxgll,2).ne.0) xigll((nxgll-1)/2+1) = zero

c---- recopier pour autres directions
      do i=1,nxgll
            yigll(i) = xigll(i)
            zigll(i) = xigll(i)
            wygll(i) = wxgll(i)
            wzgll(i) = wxgll(i)
      enddo

c---- verif
      do i=1,nxgll
            print *,'xi ',i,' = ',xigll(i)
      enddo

      print *

      do i=1,nxgll
            print *,'wx ',i,' = ',wxgll(i)
      enddo

c---- compute hprime coefficients (derivatives of Lagrange polynomials)
c---- (works only if nxgll = nygll)
      do ip=1,nxgll
            do i=1,nxgll
                   hprime(ip,i) = hdgll(ip-1,i-1,xigll,nxgll)
            enddo
      enddo

      write(*,*)
      write(*,*) 'Chaque element comporte ',nxgll,
     .      ' points dans chaque direction'
      write(*,*) 'Les fcts de forme ont ',ngnod,' noeuds'

c------------------------------------------------------

c *** afficher limites du modele
      write(*,*)
      write(*,*) 'Limites physiques du modele :'
      write(*,*)
      write(*,*) 'Xmin = ',xmin,'   Xmax = ',xmax
      write(*,*) 'Ymin = ',ymin,'   Ymax = ',ymax
      write(*,*) 'Zmin = ',zmin
      write(*,*)

c ***
c *** generer les fcts de forme 3D 8 noeuds (independant de la geometrie)
c ***

c--- cas d'un element 8 noeuds 3D (brique, Dhatt-Touzot p. 114)

      print *,'Generation des fcts de forme 3D'

      do k=1,nzgll
      do j=1,nygll
      do i=1,nxgll

      ra1 = one + xigll(i)
      ra2 = one - xigll(i)

      rb1 = one + yigll(j)
      rb2 = one - yigll(j)

      rc1 = one + zigll(k)
      rc2 = one - zigll(k)

      shape(1,i,j,k) = one_eight*ra2*rb2*rc2
      shape(2,i,j,k) = one_eight*ra1*rb2*rc2
      shape(3,i,j,k) = one_eight*ra1*rb1*rc2
      shape(4,i,j,k) = one_eight*ra2*rb1*rc2
      shape(5,i,j,k) = one_eight*ra2*rb2*rc1
      shape(6,i,j,k) = one_eight*ra1*rb2*rc1
      shape(7,i,j,k) = one_eight*ra1*rb1*rc1
      shape(8,i,j,k) = one_eight*ra2*rb1*rc1

      dershape(1,1,i,j,k) = - one_eight*rb2*rc2
      dershape(1,2,i,j,k) = one_eight*rb2*rc2
      dershape(1,3,i,j,k) = one_eight*rb1*rc2
      dershape(1,4,i,j,k) = - one_eight*rb1*rc2
      dershape(1,5,i,j,k) = - one_eight*rb2*rc1
      dershape(1,6,i,j,k) = one_eight*rb2*rc1
      dershape(1,7,i,j,k) = one_eight*rb1*rc1
      dershape(1,8,i,j,k) = - one_eight*rb1*rc1

      dershape(2,1,i,j,k) = - one_eight*ra2*rc2
      dershape(2,2,i,j,k) = - one_eight*ra1*rc2
      dershape(2,3,i,j,k) = one_eight*ra1*rc2
      dershape(2,4,i,j,k) = one_eight*ra2*rc2
      dershape(2,5,i,j,k) = - one_eight*ra2*rc1
      dershape(2,6,i,j,k) = - one_eight*ra1*rc1
      dershape(2,7,i,j,k) = one_eight*ra1*rc1
      dershape(2,8,i,j,k) = one_eight*ra2*rc1

      dershape(3,1,i,j,k) = - one_eight*ra2*rb2
      dershape(3,2,i,j,k) = - one_eight*ra1*rb2
      dershape(3,3,i,j,k) = - one_eight*ra1*rb1
      dershape(3,4,i,j,k) = - one_eight*ra2*rb1
      dershape(3,5,i,j,k) = one_eight*ra2*rb2
      dershape(3,6,i,j,k) = one_eight*ra1*rb2
      dershape(3,7,i,j,k) = one_eight*ra1*rb1
      dershape(3,8,i,j,k) = one_eight*ra2*rb1

      enddo
      enddo
      enddo

      print *,'Fin generation des fcts de forme 3D'

c--- verifier les fcts de forme et leur derivees
      print *,'Verification des fcts de forme et des derivees'

      do k=1,nzgll
      do j=1,nygll
      do i=1,nxgll

      sumshape = zero
      sumdershapexi = zero
      sumdershapeeta = zero
      sumdershapegamma = zero

      do ia=1,ngnod
            sumshape = sumshape + shape(ia,i,j,k)
            sumdershapexi = sumdershapexi +
     .                        dershape(1,ia,i,j,k)
            sumdershapeeta = sumdershapeeta +
     .                        dershape(2,ia,i,j,k)
            sumdershapegamma = sumdershapegamma +
     .                        dershape(3,ia,i,j,k)
      enddo

c la somme des shape fcts doit valoir un
c la somme des derivees des shape fcts doit valoir zero
      if(dabs(sumshape-one).gt. epsilon) stop 'erreur shape fcts'
      if(dabs(sumdershapexi).gt. epsilon)
     .            stop 'erreur deriv xi shape fcts'
      if(dabs(sumdershapeeta).gt. epsilon)
     .            stop 'erreur deriv eta shape fcts'
      if(dabs(sumdershapegamma).gt. epsilon)
     .            stop 'erreur deriv gamma shape fcts'

      enddo
      enddo
      enddo

      print *,'Fin verification des fcts de forme'

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0

      iaddx(2) = 1
      iaddy(2) = 0
      iaddz(2) = 0

      iaddx(3) = 1
      iaddy(3) = 1
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 1
      iaddz(4) = 0

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1

      iaddx(6) = 1
      iaddy(6) = 0
      iaddz(6) = 1

      iaddx(7) = 1
      iaddy(7) = 1
      iaddz(7) = 1

      iaddx(8) = 0
      iaddy(8) = 1
      iaddz(8) = 1


c ***
c *** generer les coordonnees des differents macroblocs
c ***

c---
c--- calcul des interpolateurs (lineaires pour le moment)
c---

!**** DK DK DK DK ici mettre fonctions non lineaires pour densification
      do i=0,nx
              xig(i) = dble(i)/dble(nx)
      enddo
      do j=0,ny
              etag(j) = dble(j)/dble(ny)
      enddo
      do k=0,nz1
              gammag1(k) = dble(k)/dble(nz1)
      enddo
      do k=nz1,nz1+nz2
              gammag2(k) = dble(k-nz1)/dble(nz2)
      enddo
      do k=nz1+nz2,nz1+nz2+nz3
              gammag3(k) = dble(k-(nz1+nz2))/dble(nz3)
      enddo

c**** lecture et tapering des interfaces

! remappe les interfaces de Riki et la topo/bathy vers
! une grille reguliere en (x,y) (projection de Lambert)

      print *
      print *,'definition tapering bassin tokyo'
      print *

      print *,'Grille horizontale ',nx,' x ',ny,' elements spectraux en surface'
      print *,'   (soit ',nx*(nxgll-1)+1,' x ',ny*(nygll-1)+1,' points)'
      print *,'Taille physique modele en X = ',(xmax-xmin)/1000.,' km'
      print *,'Taille physique modele en Y = ',(ymax-ymin)/1000.,' km'
      delta_x = (xmax-xmin)/dble(nx)
      delta_y = (ymax-ymin)/dble(ny)
      print *,'Delta_X element spectral = ',delta_x,' m'
      print *,'Delta_Y element spectral = ',delta_y,' m'

      print *
      print *,'Ecart min topo / couche 18 tapered = ',delta_min_topo_18,' m'
      print *,'Ecart min couche 18 / couche 28 tapered = ',delta_min_18_28,' m'
      print *

! lecture resultat
      print *
      print *,'lecture interfaces'
      iavsfile = 34
      if(.not. use_averaged_layers)
     .      open(unit=iavsfile,file='interfs_interp.dat',status='old')
      do i=1,nxinterf
      do j=1,nyinterf
      if(.not. use_averaged_layers) then
            read(iavsfile,*) xinterp(i,j),
     .      yinterp(i,j),
     .      zinterp18(i,j),
     .      zinterp28(i,j),
     .      zinterpkazu(i,j),
     .      zinterptopobathy(i,j)
      else
            xinterp(i,j) = (i-1)*500.
            yinterp(i,j) = (j-1)*500.
            zinterptopobathy(i,j) = 0.
            zinterp18(i,j) = -1500.
            zinterpkazu(i,j) = -1500.
            zinterp28(i,j) = -3000.
      endif
      enddo
      enddo
      print *,'fin lecture interfaces'
      print *
      if(.not. use_averaged_layers) close(iavsfile)

      delta_x_interf = xinterp(2,1) - xinterp(1,1)
      delta_y_interf = yinterp(1,2) - yinterp(1,1)
      print *,'Delta_X input interfaces = ',delta_x_interf
      print *,'Delta_Y input interfaces = ',delta_y_interf

      xmin_interf = xinterp(1,1)
      ymin_interf = yinterp(1,1)
      print *,'Xmin input interfaces = ',xmin_interf
      print *,'Ymin input interfaces = ',ymin_interf

! calcul et verification des hauteurs non dampees sur la grille de simu
      dmaxz1 = - bigval
      dmaxz2 = - bigval
      dmaxz3 = - bigval
      dminz1 = + bigval
      dminz2 = + bigval
      dminz3 = + bigval

      do i=0,nx
      do j=0,ny

! calcul position physique dans la grille de simu
            xsimu(i,j) = xmin*(1.d0-xig(i)) + xmax*xig(i)
            ysimu(i,j) = ymin*(1.d0-etag(j)) + ymax*etag(j)

! en deduire l'indice juste a gauche dans la grille des interfaces
            i_interf = int((xsimu(i,j) - xmin_interf) / delta_x_interf + 1)
            j_interf = int((ysimu(i,j) - ymin_interf) / delta_y_interf + 1)

! verifier si depassement du cadre du modele d'entree
            if(i_interf .gt. nxinterf) i_interf = nxinterf
            if(j_interf .gt. nyinterf) j_interf = nyinterf
            if(i_interf .lt. 1) i_interf = 1
            if(j_interf .lt. 1) j_interf = 1

! interpolation bilineaire ici (attention cas limite a droite)
! DK DK DK DK DK
            ztopo = zinterptopobathy(i_interf,j_interf)

! change here from Kazusa to Layer18 if needed
            z18 = zinterp18(i_interf,j_interf)
            z18 = zinterpkazu(i_interf,j_interf)

            z28 = zinterp28(i_interf,j_interf)
            zbot = zmin

! calculer les interfaces "tapered" pour la simu
            ztoposimu(i,j) = ztopo
            if(ztopo - z18 .ge. delta_min_topo_18) then
                  z18simu(i,j) = z18
            else
                  z18simu(i,j) = ztopo - delta_min_topo_18
            endif
            if(z18simu(i,j) - z28 .ge. delta_min_18_28) then
                  z28simu(i,j) = z28
            else
                  z28simu(i,j) = z18simu(i,j) - delta_min_18_28
            endif

! calculer les min et max des differences de hauteur entre deux couches
            delta_z0 = ztopo
            delta_z1 = ztopo - z18
            delta_z2 = z18 - z28
            delta_z3 = z28 - zbot

            if(delta_z0 .gt. dmaxz0) dmaxz0 = delta_z0
            if(delta_z1 .gt. dmaxz1) dmaxz1 = delta_z1
            if(delta_z2 .gt. dmaxz2) dmaxz2 = delta_z2
            if(delta_z3 .gt. dmaxz3) dmaxz3 = delta_z3

            if(delta_z0 .lt. dminz0) dminz0 = delta_z0
            if(delta_z1 .lt. dminz1) dminz1 = delta_z1
            if(delta_z2 .lt. dminz2) dminz2 = delta_z2
            if(delta_z3 .lt. dminz3) dminz3 = delta_z3

      enddo
      enddo

      print *
      print *,'Hauteur min et max topo et bathy :'
      print *
      print *,' Ztopo max min = ',dmaxz0,dminz0
      print *

      print *
      print *,'Hauteurs min et max entre vraies interfaces sans tapering :'
      print *
      print *,' Ztopo - Z18 max min = ',dmaxz1,dminz1
      print *,'  Z18 - Z28  max min = ',dmaxz2,dminz2
      print *,'  Z28 - Zbot max min = ',dmaxz3,dminz3
      print *

c---
c--- definition du maillage physique X,Y,Z pour le cube
c---

      zmax = - bigval
      do iz=0,nz
      do iy=0,ny
      do ix=0,nx

            x(ix,iy,iz) = xsimu(ix,iy)
            y(ix,iy,iz) = ysimu(ix,iy)

! distinguer les trois couches suivant la direction z

! couche entre le fond et layer28
      if(iz.le.nz1) then
            z(ix,iy,iz) = zmin*(1.d0-gammag1(iz)) +
     .            z28simu(ix,iy)*gammag1(iz)

! couche entre layer28 et layer18
      else if(iz.le.nz1+nz2) then
            z(ix,iy,iz) = z28simu(ix,iy)*(1.d0-gammag2(iz)) +
     .            z18simu(ix,iy)*gammag2(iz)

! couche entre layer18 et la topo/bathy
      else
            z(ix,iy,iz) = z18simu(ix,iy)*(1.d0-gammag3(iz)) +
     .            ztoposimu(ix,iy)*gammag3(iz)
      endif

            zmax = dmax1(zmax,z(ix,iy,iz))

      enddo
      enddo
      enddo

      write(*,*)
      write(*,*) 'Zmax detecte topo = ',zmax
      write(*,*)

!
!--- premiere verification grossiere des caracteristiques du maillage
!

! couche 1 : topo a layer 18
      deltax = dabs(x(1,0,nz1+nz2)   - x(0,0,nz1+nz2))
      deltay = dabs(y(0,1,nz1+nz2)   - y(0,0,nz1+nz2))
      deltaz = dabs(z(0,0,nz1+nz2+1) - z(0,0,nz1+nz2))
      cp = vp(itop)
      cs = vs(itop)
      denst = rho(itop)
      print *,'Couche 1 verif grossiere : (topo to layer18)'
      print *
      print *,'deltax 1 element spectral = ',deltax
      print *,'deltay 1 element spectral = ',deltay
      print *,'deltaz 1 element spectral = ',deltaz
      print *,'nxgll  = ',nxgll
      print *,'cp  = ',cp
      print *,'cs  = ',cs
      print *,'rho = ',denst
      print *,'Pour f0 central = ',f0
      print *,'Nb pts / lambda_P f0 direction X = ',nxgll*cp/(f0*deltax)
      print *,'Nb pts / lambda_P f0 direction Y = ',nxgll*cp/(f0*deltay)
      print *,'Nb pts / lambda_P f0 direction Z = ',nxgll*cp/(f0*deltaz)
      print *
      print *,'Nb pts / lambda_S f0 direction X = ',nxgll*cs/(f0*deltax)
      print *,'Nb pts / lambda_S f0 direction Y = ',nxgll*cs/(f0*deltay)
      print *,'Nb pts / lambda_S f0 direction Z = ',nxgll*cs/(f0*deltaz)
      print *
      print *,'Pour f_max = ',freqcste*f0
      print *,'Nb pts / lambda_P f_max direction X = ',nxgll*cp/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_P f_max direction Y = ',nxgll*cp/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_P f_max direction Z = ',nxgll*cp/(f0*deltaz*freqcste)
      print *
      print *,'Nb pts / lambda_S f_max direction X = ',nxgll*cs/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_S f_max direction Y = ',nxgll*cs/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_S f_max direction Z = ',nxgll*cs/(f0*deltaz*freqcste)
      print *

! couche 2 : layer 18 a layer 28
      deltax = dabs(x(1,0,nz1)   - x(0,0,nz1))
      deltay = dabs(y(0,1,nz1)   - y(0,0,nz1))
      deltaz = dabs(z(0,0,nz1+1) - z(0,0,nz1))
      cp = vp(imiddle)
      cs = vs(imiddle)
      denst = rho(imiddle)
      print *,'Couche 2 verif grossiere : (layer 18 to layer28)'
      print *
      print *,'deltax 1 element spectral = ',deltax
      print *,'deltay 1 element spectral = ',deltay
      print *,'deltaz 1 element spectral = ',deltaz
      print *,'nxgll  = ',nxgll
      print *,'cp  = ',cp
      print *,'cs  = ',cs
      print *,'rho = ',denst
      print *,'Pour f0 central = ',f0
      print *,'Nb pts / lambda_P f0 direction X = ',nxgll*cp/(f0*deltax)
      print *,'Nb pts / lambda_P f0 direction Y = ',nxgll*cp/(f0*deltay)
      print *,'Nb pts / lambda_P f0 direction Z = ',nxgll*cp/(f0*deltaz)
      print *
      print *,'Nb pts / lambda_S f0 direction X = ',nxgll*cs/(f0*deltax)
      print *,'Nb pts / lambda_S f0 direction Y = ',nxgll*cs/(f0*deltay)
      print *,'Nb pts / lambda_S f0 direction Z = ',nxgll*cs/(f0*deltaz)
      print *
      print *,'Pour f_max = ',freqcste*f0
      print *,'Nb pts / lambda_P f_max direction X = ',nxgll*cp/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_P f_max direction Y = ',nxgll*cp/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_P f_max direction Z = ',nxgll*cp/(f0*deltaz*freqcste)
      print *
      print *,'Nb pts / lambda_S f_max direction X = ',nxgll*cs/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_S f_max direction Y = ',nxgll*cs/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_S f_max direction Z = ',nxgll*cs/(f0*deltaz*freqcste)
      print *

! couche 3 : layer 28 a bottom (elements regroupes par quatre)
      deltax = dabs(x(4,0,0)   - x(0,0,0))
      deltay = dabs(y(0,4,0)   - y(0,0,0))
      deltaz = dabs(z(0,0,4) - z(0,0,0))
      cp = vp(ibottom)
      cs = vs(ibottom)
      denst = rho(ibottom)
      print *,'Couche 3 verif grossiere : (layer 28 to bottom)'
      print *
      print *,'deltax 1 element spectral = ',deltax
      print *,'deltay 1 element spectral = ',deltay
      print *,'deltaz 1 element spectral = ',deltaz
      print *,'nxgll  = ',nxgll
      print *,'cp  = ',cp
      print *,'cs  = ',cs
      print *,'rho = ',denst
      print *,'Pour f0 central = ',f0
      print *,'Nb pts / lambda_P f0 direction X = ',nxgll*cp/(f0*deltax)
      print *,'Nb pts / lambda_P f0 direction Y = ',nxgll*cp/(f0*deltay)
      print *,'Nb pts / lambda_P f0 direction Z = ',nxgll*cp/(f0*deltaz)
      print *
      print *,'Nb pts / lambda_S f0 direction X = ',nxgll*cs/(f0*deltax)
      print *,'Nb pts / lambda_S f0 direction Y = ',nxgll*cs/(f0*deltay)
      print *,'Nb pts / lambda_S f0 direction Z = ',nxgll*cs/(f0*deltaz)
      print *
      print *,'Pour f_max = ',freqcste*f0
      print *,'Nb pts / lambda_P f_max direction X = ',nxgll*cp/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_P f_max direction Y = ',nxgll*cp/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_P f_max direction Z = ',nxgll*cp/(f0*deltaz*freqcste)
      print *
      print *,'Nb pts / lambda_S f_max direction X = ',nxgll*cs/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_S f_max direction Y = ',nxgll*cs/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_S f_max direction Z = ',nxgll*cs/(f0*deltaz*freqcste)
      print *

! couche 6 : moho (elements regroupes par quatre)
      deltax = dabs(x(4,0,0)   - x(0,0,0))
      deltay = dabs(y(0,4,0)   - y(0,0,0))
      deltaz = dabs(z(0,0,4) - z(0,0,0))
      cp = vp(imoho)
      cs = vs(imoho)
      denst = rho(imoho)
      print *,'Couche moho verif grossiere :'
      print *
      print *,'deltax 1 element spectral = ',deltax
      print *,'deltay 1 element spectral = ',deltay
      print *,'deltaz 1 element spectral = ',deltaz
      print *,'nxgll  = ',nxgll
      print *,'cp  = ',cp
      print *,'cs  = ',cs
      print *,'rho = ',denst
      print *,'Pour f0 central = ',f0
      print *,'Nb pts / lambda_P f0 direction X = ',nxgll*cp/(f0*deltax)
      print *,'Nb pts / lambda_P f0 direction Y = ',nxgll*cp/(f0*deltay)
      print *,'Nb pts / lambda_P f0 direction Z = ',nxgll*cp/(f0*deltaz)
      print *
      print *,'Nb pts / lambda_S f0 direction X = ',nxgll*cs/(f0*deltax)
      print *,'Nb pts / lambda_S f0 direction Y = ',nxgll*cs/(f0*deltay)
      print *,'Nb pts / lambda_S f0 direction Z = ',nxgll*cs/(f0*deltaz)
      print *
      print *,'Pour f_max = ',freqcste*f0
      print *,'Nb pts / lambda_P f_max direction X = ',nxgll*cp/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_P f_max direction Y = ',nxgll*cp/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_P f_max direction Z = ',nxgll*cp/(f0*deltaz*freqcste)
      print *
      print *,'Nb pts / lambda_S f_max direction X = ',nxgll*cs/(f0*deltax*freqcste)
      print *,'Nb pts / lambda_S f_max direction Y = ',nxgll*cs/(f0*deltay*freqcste)
      print *,'Nb pts / lambda_S f_max direction Z = ',nxgll*cs/(f0*deltaz*freqcste)
      print *

!
! routine sauvegarde fichier AVS
!

      if(draw_interfs_avs) then

      print *,'Entering AVS file generation tapered interfaces...'

! file number for AVS output
      iavsfile = 34

!---- ouverture du fichier AVS (header + donnees)

! --- layer 18

      open(unit=iavsfile,file='avsheader18tapered.dat',status='unknown')
      write(iavsfile,100)
      write(iavsfile,*) 'ndim=2'
      write(iavsfile,*) 'dim1=',nx+1
      write(iavsfile,*) 'dim2=',ny+1
      write(iavsfile,*) 'nspace=3'
      write(iavsfile,*) 'veclen=1'
      write(iavsfile,*) 'data=double'
      write(iavsfile,*) 'field=irregular'
      write(iavsfile,*) 'label=Hauteur'
      write(iavsfile,*) '#'
      write(iavsfile,*) 'coord 1 file=./avsgrid18tapered.dat filetype=ascii skip=1 offset=0 stride=4'
      write(iavsfile,*) 'coord 2 file=./avsgrid18tapered.dat filetype=ascii skip=1 offset=1 stride=4'
      write(iavsfile,*) 'coord 3 file=./avsgrid18tapered.dat filetype=ascii skip=1 offset=2 stride=4'
      write(iavsfile,*) 'variable 1 file=./avsgrid18tapered.dat filetype=ascii skip=1 offset=3 stride=4'
      close(iavsfile)

      open(unit=iavsfile,file='avsgrid18tapered.dat',status='unknown')
      write(iavsfile,*) 'X     Y   Z  Height'
! numero et coordonnees des points du maillage
      do iy=0,ny
            do ix=0,nx
      write(iavsfile,400) xsimu(ix,iy),ysimu(ix,iy),z18simu(ix,iy),z18simu(ix,iy)
      enddo
      enddo
      close(iavsfile)

! --- layer 28

      open(unit=iavsfile,file='avsheader28tapered.dat',status='unknown')
      write(iavsfile,100)
      write(iavsfile,*) 'ndim=2'
      write(iavsfile,*) 'dim1=',nx+1
      write(iavsfile,*) 'dim2=',ny+1
      write(iavsfile,*) 'nspace=3'
      write(iavsfile,*) 'veclen=1'
      write(iavsfile,*) 'data=double'
      write(iavsfile,*) 'field=irregular'
      write(iavsfile,*) 'label=Hauteur'
      write(iavsfile,*) '#'
      write(iavsfile,*) 'coord 1 file=./avsgrid28tapered.dat filetype=ascii skip=1 offset=0 stride=4'
      write(iavsfile,*) 'coord 2 file=./avsgrid28tapered.dat filetype=ascii skip=1 offset=1 stride=4'
      write(iavsfile,*) 'coord 3 file=./avsgrid28tapered.dat filetype=ascii skip=1 offset=2 stride=4'
      write(iavsfile,*) 'variable 1 file=./avsgrid28tapered.dat filetype=ascii skip=1 offset=3 stride=4'
      close(iavsfile)

      open(unit=iavsfile,file='avsgrid28tapered.dat',status='unknown')
      write(iavsfile,*) 'X     Y   Z  Height'
! numero et coordonnees des points du maillage
      do iy=0,ny
            do ix=0,nx
      write(iavsfile,400) xsimu(ix,iy),ysimu(ix,iy),z28simu(ix,iy),z28simu(ix,iy)
      enddo
      enddo
      close(iavsfile)

! --- layer topo

      open(unit=iavsfile,file='avsheadertopotapered.dat',status='unknown')
      write(iavsfile,100)
      write(iavsfile,*) 'ndim=2'
      write(iavsfile,*) 'dim1=',nx+1
      write(iavsfile,*) 'dim2=',ny+1
      write(iavsfile,*) 'nspace=3'
      write(iavsfile,*) 'veclen=1'
      write(iavsfile,*) 'data=double'
      write(iavsfile,*) 'field=irregular'
      write(iavsfile,*) 'label=Hauteur'
      write(iavsfile,*) '#'
      write(iavsfile,*) 'coord 1 file=./avsgridtopotapered.dat filetype=ascii skip=1 offset=0 stride=4'
      write(iavsfile,*) 'coord 2 file=./avsgridtopotapered.dat filetype=ascii skip=1 offset=1 stride=4'
      write(iavsfile,*) 'coord 3 file=./avsgridtopotapered.dat filetype=ascii skip=1 offset=2 stride=4'
      write(iavsfile,*) 'variable 1 file=./avsgridtopotapered.dat filetype=ascii skip=1 offset=3 stride=4'
      close(iavsfile)

      open(unit=iavsfile,file='avsgridtopotapered.dat',status='unknown')
      write(iavsfile,*) 'X     Y   Z  Height'
! numero et coordonnees des points du maillage
      do iy=0,ny
            do ix=0,nx
      write(iavsfile,400) xsimu(ix,iy),ysimu(ix,iy),ztoposimu(ix,iy),ztoposimu(ix,iy)
      enddo
      enddo
      close(iavsfile)

      print *,'End of AVS file generation tapered interfaces...'

      endif

 100  format('# AVS field file')
 400  format(e12.5,1x,e12.5,1x,e12.5,1x,e12.5)

c**** fin lecture et tapering des interfaces

c ***
c *** generer un fichier 'GNUPLOT' pour le controle de la grille ***
c ***

      if(ignuplot) then

      write(*,*)' Ecriture de la grille de base format GNUPLOT...'

      open(unit=20,file='basegrid.GNU',status='unknown')
      iz = 0

c base cubique
      do ix=0,nx
              do iy=0,ny-1
                write(20,115) x(ix,iy,iz),y(ix,iy,iz)
                write(20,115) x(ix,iy+1,iz),y(ix,iy+1,iz)
              write(20,30)
              enddo
      enddo

      do ix=0,nx-1
              do iy=0,ny
                write(20,115) x(ix,iy,iz),y(ix,iy,iz)
                write(20,115) x(ix+1,iy,iz),y(ix+1,iy,iz)
              write(20,30)
              enddo
      enddo

      close(20)

      open(unit=20,file='plotbasegrid.gnu',status='unknown')
      write(20,*) 'set term x11'
      write(20,*) 'set xlabel "X"'
      write(20,*) 'set ylabel "Y"'
      write(20,*) 'set title "Base de la grille"'
      write(20,110)
      write(20,*) 'pause -1 ''Hit any key...'''
      close(20)

 110   format('plot "basegrid.GNU" using 1:2 t '''' w l 1')
 115   format(e12.5,1x,e12.5)

      write(*,*)' Fin ecriture de la grille format GNUPLOT'

      endif

c---
c--- generer la metrique et le jacobien pour les blocs structures
c---

      xximin =  + bigval
      yximin =  + bigval
      zximin =  + bigval
      xetamin = + bigval
      yetamin = + bigval
      zetamin = + bigval
      xgammamin = + bigval
      ygammamin = + bigval
      zgammamin = + bigval
      xximax =  - bigval
      yximax =  - bigval
      zximax =  - bigval
      xetamax = - bigval
      yetamax = - bigval
      zetamax = - bigval
      xgammamax = - bigval
      ygammamax = - bigval
      zgammamax = - bigval

      xixmin =  + bigval
      xiymin =  + bigval
      xizmin =  + bigval
      etaxmin = + bigval
      etaymin = + bigval
      etazmin = + bigval
      gammaxmin = + bigval
      gammaymin = + bigval
      gammazmin = + bigval
      xixmax =  - bigval
      xiymax =  - bigval
      xizmax =  - bigval
      etaxmax = - bigval
      etaymax = - bigval
      etazmax = - bigval
      gammaxmax = - bigval
      gammaymax = - bigval
      gammazmax = - bigval

      print *,'Generation metrique et jacobien pour le cube'

!--- initialisation de la numerotation des elements spectraux
      ispec = 0

!--- generation du bloc des deux couches du bassin (layers 28, 18 et topo)

      print *
      print *,'Generating the two sedimentary layers'
      print *

      igenlay = 0

      do iz=nz1,nz-1

      igenlay = igenlay + 1

       print *,'generating for iz = ',igenlay,' out of ',nz2 + nz3

      do ix=0,nx-1
      do iy=0,ny-1

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)

c couche entre layer18 et topo/bathy
      if(iz .ge. nz1+nz2) then
            numat(ispec) = itop

c couche entre layer28 et layer18
      else
            numat(ispec) = imiddle
      endif

      enddo
      enddo
      enddo

!--- generation du bloc de taille quatre sous le deraffinement

      print *
      print *,'Generating the coarse block in the bedrock'
      print *

      igenlay = 0

      do iz=0,nz1 - 8 - 4,4

      igenlay = igenlay + 1

       print *,'generating for iz = ',igenlay,' out of ',(nz1 - 8)/4

      do ix=0,nx-4,4
      do iy=0,ny-4,4

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+4*iaddx(ia),iy+4*iaddy(ia),iz+4*iaddz(ia))
            ycell(ia) = y(ix+4*iaddx(ia),iy+4*iaddy(ia),iz+4*iaddz(ia))
            zcell(ia) = z(ix+4*iaddx(ia),iy+4*iaddy(ia),iz+4*iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo
      enddo

! $$$$$$$$$$$$$$$$$$$ DEBUT RACCORD $$$$$$$$$$$$$$$$$$$$$$$

!--- generation du raccord de deraffinement

      print *
      print *,'Generating the matching layer'
      print *

! bloc 1 sur 6

      print *,'bloc 1 upper level'

      iz = nz1 - 1

      do ix=0,nx-2,2
      do iy=0,ny-1,1

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0
      if(mod(ix,4) .ne. 0) iaddz(1) = -1

      iaddx(2) = 1
      iaddy(2) = 0
      iaddz(2) = 0

      iaddx(3) = 1
      iaddy(3) = 1
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 1
      iaddz(4) = 0
      if(mod(ix,4) .ne. 0) iaddz(4) = -1

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1

      iaddx(6) = 1
      iaddy(6) = 0
      iaddz(6) = 1

      iaddx(7) = 1
      iaddy(7) = 1
      iaddz(7) = 1

      iaddx(8) = 0
      iaddy(8) = 1
      iaddz(8) = 1

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 2 sur 6

      print *,'bloc 2 upper level'

      iz = nz1 - 2

      do ix=1,nx-1,2
      do iy=0,ny-1,1

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 1

      iaddx(2) = 1
      iaddy(2) = 0
      iaddz(2) = 0
      if(mod(ix-1,4) .ne. 0) iaddz(2) = 1

      iaddx(3) = 1
      iaddy(3) = 1
      iaddz(3) = 0
      if(mod(ix-1,4) .ne. 0) iaddz(3) = 1

      iaddx(4) = 0
      iaddy(4) = 1
      iaddz(4) = 1

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 2

      iaddx(6) = 1
      iaddy(6) = 0
      iaddz(6) = 2

      iaddx(7) = 1
      iaddy(7) = 1
      iaddz(7) = 2

      iaddx(8) = 0
      iaddy(8) = 1
      iaddz(8) = 2

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 3 sur 6

      print *,'bloc 3 upper level'

      iz = nz1 - 2

      do ix=0,nx-2,2
      do iy=0,ny-1,1

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0

      iaddx(2) = 2
      iaddy(2) = 0
      iaddz(2) = 0

      iaddx(3) = 2
      iaddy(3) = 1
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 1
      iaddz(4) = 0

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1
      if(mod(ix,4) .ne. 0) iaddx(5) = 1

      iaddx(6) = 1
      iaddy(6) = 0
      iaddz(6) = 1
      if(mod(ix,4) .ne. 0) iaddx(6) = 2

      iaddx(7) = 1
      iaddy(7) = 1
      iaddz(7) = 1
      if(mod(ix,4) .ne. 0) iaddx(7) = 2

      iaddx(8) = 0
      iaddy(8) = 1
      iaddz(8) = 1
      if(mod(ix,4) .ne. 0) iaddx(8) = 1

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 4 sur 6

      print *,'bloc 4 upper level'

      iz = nz1 - 3

      do ix=0,nx-2,2
      do iy=0,ny-2,2

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0
      if(mod(iy,4) .ne. 0) iaddz(1) = -1

      iaddx(2) = 2
      iaddy(2) = 0
      iaddz(2) = 0
      if(mod(iy,4) .ne. 0) iaddz(2) = -1

      iaddx(3) = 2
      iaddy(3) = 1
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 1
      iaddz(4) = 0

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1

      iaddx(6) = 2
      iaddy(6) = 0
      iaddz(6) = 1

      iaddx(7) = 2
      iaddy(7) = 1
      iaddz(7) = 1

      iaddx(8) = 0
      iaddy(8) = 1
      iaddz(8) = 1

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 5 sur 6

      print *,'bloc 5 upper level'

      iz = nz1 - 4

      do ix=0,nx-2,2
      do iy=1,ny-1,2

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 1

      iaddx(2) = 2
      iaddy(2) = 0
      iaddz(2) = 1

      iaddx(3) = 2
      iaddy(3) = 1
      iaddz(3) = 0
      if(mod(iy-1,4) .ne. 0) iaddz(3) = 1

      iaddx(4) = 0
      iaddy(4) = 1
      iaddz(4) = 0
      if(mod(iy-1,4) .ne. 0) iaddz(4) = 1

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 2

      iaddx(6) = 2
      iaddy(6) = 0
      iaddz(6) = 2

      iaddx(7) = 2
      iaddy(7) = 1
      iaddz(7) = 2

      iaddx(8) = 0
      iaddy(8) = 1
      iaddz(8) = 2

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 6 sur 6

      print *,'bloc 6 upper level'

      iz = nz1 - 4

      do ix=0,nx-2,2
      do iy=0,ny-2,2

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0

      iaddx(2) = 2
      iaddy(2) = 0
      iaddz(2) = 0

      iaddx(3) = 2
      iaddy(3) = 2
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 2
      iaddz(4) = 0

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1
      if(mod(iy,4) .ne. 0) iaddy(5) = 1

      iaddx(6) = 2
      iaddy(6) = 0
      iaddz(6) = 1
      if(mod(iy,4) .ne. 0) iaddy(6) = 1

      iaddx(7) = 2
      iaddy(7) = 1
      iaddz(7) = 1
      if(mod(iy,4) .ne. 0) iaddy(7) = 2

      iaddx(8) = 0
      iaddy(8) = 1
      iaddz(8) = 1
      if(mod(iy,4) .ne. 0) iaddy(8) = 2

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo


! %%%%%% deuxieme couche du bloc de deraffinement



! bloc 1 sur 6

      print *,'bloc 1 lower level'

      iz = nz1 - 1 - 4

      do ix=0,nx-4,4
      do iy=0,ny-2,2

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0
      if(mod(ix,8) .ne. 0) iaddz(1) = -1

      iaddx(2) = 2
      iaddy(2) = 0
      iaddz(2) = 0

      iaddx(3) = 2
      iaddy(3) = 2
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 2
      iaddz(4) = 0
      if(mod(ix,8) .ne. 0) iaddz(4) = -1

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1

      iaddx(6) = 2
      iaddy(6) = 0
      iaddz(6) = 1

      iaddx(7) = 2
      iaddy(7) = 2
      iaddz(7) = 1

      iaddx(8) = 0
      iaddy(8) = 2
      iaddz(8) = 1

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 2 sur 6

      print *,'bloc 2 lower level'

      iz = nz1 - 2 - 4

      do ix=2,nx-2,4
      do iy=0,ny-2,2

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 1

      iaddx(2) = 2
      iaddy(2) = 0
      iaddz(2) = 0
      if(mod(ix-2,8) .ne. 0) iaddz(2) = 1

      iaddx(3) = 2
      iaddy(3) = 2
      iaddz(3) = 0
      if(mod(ix-2,8) .ne. 0) iaddz(3) = 1

      iaddx(4) = 0
      iaddy(4) = 2
      iaddz(4) = 1

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 2

      iaddx(6) = 2
      iaddy(6) = 0
      iaddz(6) = 2

      iaddx(7) = 2
      iaddy(7) = 2
      iaddz(7) = 2

      iaddx(8) = 0
      iaddy(8) = 2
      iaddz(8) = 2

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 3 sur 6

      print *,'bloc 3 lower level'

      iz = nz1 - 2 - 4

      do ix=0,nx-4,4
      do iy=0,ny-2,2

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0

      iaddx(2) = 4
      iaddy(2) = 0
      iaddz(2) = 0

      iaddx(3) = 4
      iaddy(3) = 2
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 2
      iaddz(4) = 0

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1
      if(mod(ix,8) .ne. 0) iaddx(5) = 2

      iaddx(6) = 2
      iaddy(6) = 0
      iaddz(6) = 1
      if(mod(ix,8) .ne. 0) iaddx(6) = 4

      iaddx(7) = 2
      iaddy(7) = 2
      iaddz(7) = 1
      if(mod(ix,8) .ne. 0) iaddx(7) = 4

      iaddx(8) = 0
      iaddy(8) = 2
      iaddz(8) = 1
      if(mod(ix,8) .ne. 0) iaddx(8) = 2

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 4 sur 6

      print *,'bloc 4 lower level'

      iz = nz1 - 3 - 4

      do ix=0,nx-4,4
      do iy=0,ny-4,4

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0
      if(mod(iy,8) .ne. 0) iaddz(1) = -1

      iaddx(2) = 4
      iaddy(2) = 0
      iaddz(2) = 0
      if(mod(iy,8) .ne. 0) iaddz(2) = -1

      iaddx(3) = 4
      iaddy(3) = 2
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 2
      iaddz(4) = 0

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1

      iaddx(6) = 4
      iaddy(6) = 0
      iaddz(6) = 1

      iaddx(7) = 4
      iaddy(7) = 2
      iaddz(7) = 1

      iaddx(8) = 0
      iaddy(8) = 2
      iaddz(8) = 1

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 5 sur 6

      print *,'bloc 5 lower level'

      iz = nz1 - 4 - 4

      do ix=0,nx-4,4
      do iy=2,ny-2,4

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 1

      iaddx(2) = 4
      iaddy(2) = 0
      iaddz(2) = 1

      iaddx(3) = 4
      iaddy(3) = 2
      iaddz(3) = 0
      if(mod(iy-2,8) .ne. 0) iaddz(3) = 1

      iaddx(4) = 0
      iaddy(4) = 2
      iaddz(4) = 0
      if(mod(iy-2,8) .ne. 0) iaddz(4) = 1

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 2

      iaddx(6) = 4
      iaddy(6) = 0
      iaddz(6) = 2

      iaddx(7) = 4
      iaddy(7) = 2
      iaddz(7) = 2

      iaddx(8) = 0
      iaddy(8) = 2
      iaddz(8) = 2

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo



! bloc 6 sur 6

      print *,'bloc 6 lower level'

      iz = nz1 - 4 - 4

      do ix=0,nx-4,4
      do iy=0,ny-4,4

c definition de la topologie de l'element de controle
c (numerotation de Dhatt-Touzot p. 114, brique 8 noeuds 3D)
      iaddx(1) = 0
      iaddy(1) = 0
      iaddz(1) = 0

      iaddx(2) = 4
      iaddy(2) = 0
      iaddz(2) = 0

      iaddx(3) = 4
      iaddy(3) = 4
      iaddz(3) = 0

      iaddx(4) = 0
      iaddy(4) = 4
      iaddz(4) = 0

      iaddx(5) = 0
      iaddy(5) = 0
      iaddz(5) = 1
      if(mod(iy,8) .ne. 0) iaddy(5) = 2

      iaddx(6) = 4
      iaddy(6) = 0
      iaddz(6) = 1
      if(mod(iy,8) .ne. 0) iaddy(6) = 2

      iaddx(7) = 4
      iaddy(7) = 2
      iaddz(7) = 1
      if(mod(iy,8) .ne. 0) iaddy(7) = 4

      iaddx(8) = 0
      iaddy(8) = 2
      iaddz(8) = 1
      if(mod(iy,8) .ne. 0) iaddy(8) = 4

c definition des coordonnees des noeuds de la cellule
      do ia=1,ngnod
            xcell(ia) = x(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            ycell(ia) = y(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
            zcell(ia) = z(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
      enddo

c calcul de la matrice jacobienne a partir des fonctions de forme
      call calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

! marqueur de materiau pour chaque element spectral (indicateur de couche)
      numat(ispec) = ibottom

      enddo
      enddo

! $$$$$$$$$$$$$$$$$$$ FIN RACCORD $$$$$$$$$$$$$$$$$$$$$$$

      print *
      print *,'Valeur finale de ispec = ',ispec
      print *,'Valeur theorique de nspec = ',nspec
      print *

      if(ispec .ne. nspec) stop 'desaccord nombre d''elements spectraux'

! $$$$$$$$$$$$$$$$$$$ GENERATION MODELE DE VITESSE $$$$$$$$$$$$$$$$$$$$$$

      print *
      print *,'Generation du modele de vitesse'
      print *

!
! DK DK DK supprime pour flat layers model
!
!     print *
!     print *,'**** lecture modele de vitesse binaire Riki ****'
!     print *
!
!     print *
!     print *,' nxveloc = ',nxveloc
!     print *,' nyveloc = ',nyveloc
!     print *,' nzveloc = ',nzveloc
!     print *
!
!     print *,'Reading external velocity model'
!     open(unit=27,file='../../VELOCITY_MODEL/velocmodelfiltered.bin',
!    .      status='old',form='unformatted')
!     read(27) velocmodel
!     close(27)

c assignation du modele aux elements spectraux

      do ispec = 1, nspec

      do k=1,nzgll
            do j=1,nygll
                  do i=1,nxgll

! recuperer le modele de vitesse P et S par l'indicateur de couche
      cp = vp(numat(ispec))
      cs = vs(numat(ispec))
      denst = rho(numat(ispec))

! inclure ce modele dans le modele 1D regional de Riki
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth12) then
            cp = vp(i12_25)
            cs = vs(i12_25)
            denst = rho(i12_25)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth25) then
            cp = vp(i25_moho)
            cs = vs(i25_moho)
            denst = rho(i25_moho)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depthmoho) then
            cp = vp(imoho)
            cs = vs(imoho)
            denst = rho(imoho)
      endif

      rlambda = denst*(cp*cp - 2.d0*cs*cs)
      rmu     = denst*cs*cs

      rlambdastock(i,j,k,ispec) = sngl(rlambda)
      rmustock(i,j,k,ispec)     = sngl(rmu)

                  enddo
            enddo
      enddo

      enddo

! $$$$$$$$$$$$$$$$$$$ VERIFICATION GRILLE PHYSIQUE $$$$$$$$$$$$$$$$$$$$$$$

      print *
      print *,'Verification des termes metriques (eventuellement nuls)'
      print *

      print *
      print *,'Termes directs :'
      print *

      print *,'xximin  = ',xximin
      print *,'xximax  = ',xximax
      print *,'yximin  = ',yximin
      print *,'yximax  = ',yximax
      print *,'zximin  = ',zximin
      print *,'zximax  = ',zximax
      print *
      print *,'xetamin = ',xetamin
      print *,'xetamax = ',xetamax
      print *,'yetamin = ',yetamin
      print *,'yetamax = ',yetamax
      print *,'zetamin = ',zetamin
      print *,'zetamax = ',zetamax
      print *
      print *,'xgammamin = ',xgammamin
      print *,'xgammamax = ',xgammamax
      print *,'ygammamin = ',ygammamin
      print *,'ygammamax = ',ygammamax
      print *,'zgammamin = ',zgammamin
      print *,'zgammamax = ',zgammamax
      print *

      print *
      print *,'Termes inverses :'
      print *

      print *,'xixmin  = ',xixmin
      print *,'xixmax  = ',xixmax
      print *,'xiymin  = ',xiymin
      print *,'xiymax  = ',xiymax
      print *,'xizmin  = ',xizmin
      print *,'xizmax  = ',xizmax
      print *
      print *,'etaxmin = ',etaxmin
      print *,'etaxmax = ',etaxmax
      print *,'etaymin = ',etaymin
      print *,'etaymax = ',etaymax
      print *,'etazmin = ',etazmin
      print *,'etazmax = ',etazmax
      print *
      print *,'gammaxmin = ',gammaxmin
      print *,'gammaxmax = ',gammaxmax
      print *,'gammaymin = ',gammaymin
      print *,'gammaymax = ',gammaymax
      print *,'gammazmin = ',gammazmin
      print *,'gammazmax = ',gammazmax
      print *

!
! DK DK DK mettre ici boucle sur ispec et calcul deltax
! deltay deltaz et np pts / lamdba en fonction du materiau numat(ispec)
!

! $$$$$$$$$$$$$$$$$$$ SAUVEGARDE METRIQUE $$$$$$$$$$$$$$$$$$$$$$$

      if(save_binary_files) then

      print *,'stockage des termes de la metrique'

c --- dvolu
c DK DK march99      print *,'Saving array dvolu'
c DK DK march99      open(unit=27,file='dvolu.bin',status='unknown',
c DK DK march99     .            form='unformatted')
c DK DK march99      write(27) dvolustock
c DK DK march99      close(27)

c --- xix
      print *,'Saving array xix'
      open(unit=27,file='xix.bin',status='unknown',
     .            form='unformatted')
      write(27) xixstock
      close(27)

c --- xiy
      print *,'Saving array xiy'
      open(unit=27,file='xiy.bin',status='unknown',
     .            form='unformatted')
      write(27) xiystock
      close(27)

c --- xiz
      print *,'Saving array xiz'
      open(unit=27,file='xiz.bin',status='unknown',
     .            form='unformatted')
      write(27) xizstock
      close(27)

c --- etax
      print *,'Saving array etax'
      open(unit=27,file='etax.bin',status='unknown',
     .            form='unformatted')
      write(27) etaxstock
      close(27)

c --- etay
      print *,'Saving array etay'
      open(unit=27,file='etay.bin',status='unknown',
     .            form='unformatted')
      write(27) etaystock
      close(27)

c --- etaz
      print *,'Saving array etaz'
      open(unit=27,file='etaz.bin',status='unknown',
     .            form='unformatted')
      write(27) etazstock
      close(27)

c --- gammax
      print *,'Saving array gammax'
      open(unit=27,file='gammax.bin',status='unknown',
     .            form='unformatted')
      write(27) gammaxstock
      close(27)

c --- gammay
      print *,'Saving array gammay'
      open(unit=27,file='gammay.bin',status='unknown',
     .            form='unformatted')
      write(27) gammaystock
      close(27)

c --- gammaz
      print *,'Saving array gammaz'
      open(unit=27,file='gammaz.bin',status='unknown',
     .            form='unformatted')
      write(27) gammazstock
      close(27)

c --- xmesh
!      print *,'Saving array xmesh'
!      open(unit=27,file='xmesh.bin',status='unknown',
!     .            form='unformatted')
!      write(27) xmeshstock
!      close(27)

c --- ymesh
!      print *,'Saving array ymesh'
!      open(unit=27,file='ymesh.bin',status='unknown',
!     .            form='unformatted')
!      write(27) ymeshstock
!      close(27)

c --- zmesh
!      print *,'Saving array zmesh'
!      open(unit=27,file='zmesh.bin',status='unknown',
!     .            form='unformatted')
!      write(27) zmeshstock
!      close(27)

      print *,'Fin stockage metrique et jacobien'

      endif

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if(save_binary_files) then

      print *,'stockage modele de vitesse'

c --- rlambda
      print *,'Saving array rlambda'
      open(unit=27,file='rlambda.bin',status='unknown',
     .            form='unformatted')
      write(27) rlambdastock
      close(27)

c --- rmu
      print *,'Saving array rmu'
      open(unit=27,file='rmu.bin',status='unknown',
     .            form='unformatted')
      write(27) rmustock
      close(27)

      print *,'Fin stockage modele de vitesse'

      endif

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c creation tableaux a1..a9

      print *
      print *,'debut creation tableaux a1 a a9'

c---- a1 a a9
      do iz=1,nzgll
            do iy=1,nygll
                  do ix=1,nxgll
            a1(ix,iy,iz) = sngl(hprime(iy,ix))
            a2(ix,iy,iz) = sngl(hprime(ix,iy))
            a3(ix,iy,iz) = sngl(hprime(iy,iz))
            a4(ix,iy,iz) = sngl(wxgll(iy)*hprime(ix,iy))
            a5(ix,iy,iz) = sngl(wxgll(ix)*hprime(iy,ix))
            a6(ix,iy,iz) = sngl(wxgll(iy)*hprime(iz,iy))
            a7(ix,iy,iz) = sngl(wxgll(iy)*wxgll(iz))
            a8(ix,iy,iz) = sngl(wxgll(ix)*wxgll(iz))
            a9(ix,iy,iz) = sngl(wxgll(ix)*wxgll(iy))
                  enddo
            enddo
      enddo

c---- sauvegarde tableaux a1 a a9
      if(save_binary_files) then
      open(unit=27,file='a1_a9.bin',status='unknown',
     .            form='unformatted')
      write(27) a1
      write(27) a2
      write(27) a3
      write(27) a4
      write(27) a5
      write(27) a6
      write(27) a7
      write(27) a8
      write(27) a9
      close(27)
      endif

      print *,'fin creation tableaux a1 a a9'

c %%%%%%%%%%%%%%%%%%%

! stopper si verification de la grille seule demandee
      if(.not. iexec) stop 'Verification seule demandee, terminee'

c %%%%%%%%%%%%%%%%%%%

c creation du tableau ibool

      print *, 'creation tableau ibool (par routine Fischer)'

      do ispec=1,nspec
         ieoff = nxyz*(ispec - 1)
         ilocnum = 0
      do iz = 1,nxgll
      do iy = 1,nxgll
      do ix = 1,nxgll
            ilocnum = ilocnum + 1
            xp(ilocnum + ieoff) = xmeshstock(ix,iy,iz,ispec)
            yp(ilocnum + ieoff) = ymeshstock(ix,iy,iz,ispec)
            zp(ilocnum + ieoff) = zmeshstock(ix,iy,iz,ispec)
      enddo
      enddo
      enddo
      enddo

c $$$$$$$ routine Fischer ici $$$$$$$

c    Establish initial pointers
      do ie=1,nspec
         ieoff = nxyz*(ie - 1)
         do ix=1,nxyz
            loc (ix+ieoff) = ix+ieoff
         enddo
      enddo
c
c    Set up local geometric tolerances
c
      xminval = + sngl(bigval)
      yminval = + sngl(bigval)
      zminval = + sngl(bigval)
      xmaxval = - sngl(bigval)
      ymaxval = - sngl(bigval)
      zmaxval = - sngl(bigval)
      do ie=1,nspec
         ieoff = nxyz*(ie-1) + 1
         do i=1,nxyz
            xmaxval = amax1(xp(ieoff+i-1),xmaxval)
            xminval = amin1(xp(ieoff+i-1),xminval)
            ymaxval = amax1(yp(ieoff+i-1),ymaxval)
            yminval = amin1(yp(ieoff+i-1),yminval)
            zmaxval = amax1(zp(ieoff+i-1),zmaxval)
            zminval = amin1(zp(ieoff+i-1),zminval)
         enddo
         xtol = smallvaltol*
     .       ((xmaxval-xminval)+(ymaxval-yminval)+(zmaxval-zminval))
         do i=1,nxyz
            tol(ieoff+i-1) = xtol
         enddo
      enddo

      do i=1,ntot
            ifseg(i) = .false.
      enddo
      nseg        = 1
      ifseg(1)    = .true.
      ninseg(1)   = ntot

      do j=1,ndime

c       Sort within each segment
         ioff=1
         do iseg=1,nseg
            if  (j.eq.1) then
               call rank (xp(ioff),ind,ninseg(iseg))
            elseif  (j.eq.2) then
               call rank (yp(ioff),ind,ninseg(iseg))
            else
               call rank (zp(ioff),ind,ninseg(iseg))
            endif
            call swap (xp(ioff),work,ind,ninseg(iseg))
            call swap (yp(ioff),work,ind,ninseg(iseg))
            call swap (zp(ioff),work,ind,ninseg(iseg))
            call swap (tol  (ioff),work,ind,ninseg(iseg))
            call iswap (loc (ioff),work,ind,ninseg(iseg))
            ioff=ioff+ninseg(iseg)
         enddo

c       Check for jumps in current coordinate
         if (j.eq.1) then
           do i=2,ntot
           if (abs(xp(i)-xp(i-1)).gt.(tol(i)+tol(i-1)))
     .        ifseg(i)=.true.
           enddo
         elseif (j.eq.2) then
           do i=2,ntot
           if (abs(yp(i)-yp(i-1)).gt.(tol(i)+tol(i-1)))
     .        ifseg(i)=.true.
           enddo
         else
           do i=2,ntot
           if (abs(zp(i)-zp(i-1)).gt.(tol(i)+tol(i-1)))
     .        ifseg(i)=.true.
           enddo
         endif

c       Count up number of different segments
         nseg = 0
         do i=1,ntot
            if (ifseg(i)) then
               nseg = nseg+1
               ninseg(nseg) = 1
            else
               ninseg(nseg) = ninseg(nseg) + 1
            endif
         enddo
      enddo

C     Assign global node numbers (now sorted lexicographically)
      ig = 0
      do i=1,ntot
         if (ifseg(i)) ig=ig+1
         iglob(loc(i)) = ig
      enddo

      nglob = ig

c $$$$$$$ fin routine Fischer ici $$$$$$$

      print *,' result routine brown : nglob = ',nglob
      print *,' npoin theorique = ',npoin

      if(nglob .ne. npoin) stop 'desaccord numerotation globale'

c recuperer resultat a mon format
      do ispec=1,nspec
         ieoff = nxyz*(ispec - 1)
         ilocnum = 0
      do iz = 1,nxgll
      do iy = 1,nxgll
      do ix = 1,nxgll
            ilocnum = ilocnum + 1
            ibool(ix,iy,iz,ispec) = iglob(ilocnum + ieoff)
      enddo
      enddo
      enddo
      enddo

      iboolmin = + bigvalint
      iboolmax = - bigvalint
      do ispec=1,nspec
      do iz=1,nzgll
            do iy=1,nygll
                  do ix=1,nxgll
            iboolmin = min(iboolmin,ibool(ix,iy,iz,ispec))
            iboolmax = max(iboolmax,ibool(ix,iy,iz,ispec))
                  enddo
            enddo
      enddo
      enddo
      if(iboolmin .ne. 1 .or. iboolmax .ne. npoin)
     .      stop 'incorrect global numbering generated'

c ---- sauvegarder ibool
      if(save_binary_files) then
      print *
      print *,'sauvegarde tableau ibool'
      open(unit=27,file='ibool.bin',status='unknown',
     .            form='unformatted')
      write(27) ibool
      close(27)
      print *,'fin sauvegarde tableau ibool'
      endif

c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c compute the mass matrix by summing the contribution of each point

      print *
      print *,' debut creation matrice de masse'

      do i=1,npoin
            rmass(i) = sngl(zero)
      enddo

      do ispec = 1,nspec

      do k=1,nzgll
            do j=1,nygll
                  do i=1,nxgll

! recuperer le modele de densite par l'indicateur de couche
      denst = rho(numat(ispec))

!
! DK DK DK modele externe supprime pour test Bouchon flat layers
!
! coordonnees du point actuel
!      xcoordpt = dble(xmeshstock(i,j,k,ispec)) + xmin
!      ycoordpt = dble(ymeshstock(i,j,k,ispec)) + ymin
!      zcoordpt = dble(zmeshstock(i,j,k,ispec)) + zmin
!
! point equivalent dans le bloc de vitesse
!      i_veloc = int((xcoordpt - xmin_veloc) / deltax_veloc + 1)
!      j_veloc = int((ycoordpt - ymin_veloc) / deltay_veloc + 1)
!      k_veloc = int((zcoordpt - zmin_veloc) / deltaz_veloc + 1)
!
! DK DK DK codage en dur provisoire de l'inversion
! DK DK DK due au signe moins sur la profondeur oublie
!      k_veloc = 90 - k_veloc + 1
!
! verifier si depassement du cadre du modele d'entree
!      if(i_veloc .gt. nxveloc) i_veloc = nxveloc
!      if(j_veloc .gt. nyveloc) j_veloc = nyveloc
!      if(k_veloc .gt. nzveloc) k_veloc = nzveloc
!      if(i_veloc .lt. 1) i_veloc = 1
!      if(j_veloc .lt. 1) j_veloc = 1
!      if(k_veloc .lt. 1) k_veloc = 1
!
! aller chercher le resultat dans le bloc de vitesse
!      cp = velocmodel(i_veloc,j_veloc,k_veloc)
!
! blindage sur la vitesse lue
!      if(cp .lt. vp(itop) - small_velocity_jump .or.
!     .      cp .gt. vp(ibottom) + small_velocity_jump)
!     .            stop 'vitesse incoherente'
!
! en deduire cp, cs et rho pour le milieu homogene dans chaque couche
!      if(cp .le. vp(itop) + small_velocity_jump) then
!            cp = vp(itop)
!            cs = vs(itop)
!            denst = rho(itop)
!      else if(cp .le. vp(imiddle) + small_velocity_jump) then
!            cp = vp(imiddle)
!            cs = vs(imiddle)
!            denst = rho(imiddle)
!      else
!            cp = vp(ibottom)
!            cs = vs(ibottom)
!            denst = rho(ibottom)
!      endif

! inclure ce modele dans le modele 1D regional de Riki
      if((dble(zmeshstock(i,j,k,ispec))+zmin) .lt. depth12)
     .            denst = rho(i12_25)
      if((dble(zmeshstock(i,j,k,ispec))+zmin) .lt. depth25)
     .            denst = rho(i25_moho)
      if((dble(zmeshstock(i,j,k,ispec))+zmin) .lt. depthmoho)
     .            denst = rho(imoho)

      iglobnum = ibool(i,j,k,ispec)

      rmass(iglobnum) = rmass(iglobnum) +
     .  sngl(denst*wxgll(i)*wygll(j)*wzgll(k)*dble(dvolustock(i,j,k,ispec)))

                  enddo
            enddo
      enddo

      enddo

c---- rmass
      if(save_binary_files) then
      print *,'sauvegarde matrice de masse'
      open(unit=27,file='rmass.bin',status='unknown',
     .            form='unformatted')
      write(27) rmass
      close(27)
      print *,'fin sauvegarde matrice de masse'
      endif

      print *, ' fin creation matrice de masse'
      print *

c----
c---- find the position of the source (closest grid point)
c----

      print *
      print *,' recherche position du point source'
      print *

      distmin = bigval

      do ispec = 1,nspec
      do k=1,nzgll
            do j=1,nygll
                  do i=1,nxgll

            dist = dsqrt(
     .        ((xpositsource - xmin ) - dble(xmeshstock(i,j,k,ispec))) ** 2
     .      + ((ypositsource - ymin ) - dble(ymeshstock(i,j,k,ispec))) ** 2
     .      + ((zpositsource - zmin ) - dble(zmeshstock(i,j,k,ispec))) ** 2)

            if(dist .lt. distmin) then
                  distmin = dist
                  ipositsource = ibool(i,j,k,ispec)
                  ispecsource = ispec
                  ixgllsource = i
                  iygllsource = j
                  izgllsource = k
            endif

                  enddo
            enddo
      enddo
      enddo

      print *,'Point number of closest point for the source : ',ipositsource
      print *,'Distance from the original position of the source : ',distmin
      print *
      print *,'Real position of the source that will be used :'
      print *,'     x_s = ',
     . dble(xmeshstock(ixgllsource,iygllsource,izgllsource,ispecsource)) + xmin
      print *,'     y_s = ',
     . dble(ymeshstock(ixgllsource,iygllsource,izgllsource,ispecsource)) + ymin
      print *,'     z_s = ',
     . dble(zmeshstock(ixgllsource,iygllsource,izgllsource,ispecsource)) + zmin
      print *
      print *,'Decalage position source reelle / demandee :'
      print *,'     delta_x_s = ',
     . dble(xmeshstock(ixgllsource,iygllsource,izgllsource,ispecsource))
     .      - (xpositsource - xmin)
      print *,'     delta_y_s = ',
     . dble(ymeshstock(ixgllsource,iygllsource,izgllsource,ispecsource))
     .      - (ypositsource - ymin)
      print *,'     delta_z_s = ',
     . dble(zmeshstock(ixgllsource,iygllsource,izgllsource,ispecsource))
     .      - (zpositsource - zmin)
      print *

c %%%%%%%%%%%%%%%%%%

c creation tableaux a11 a a13 source explosive

      print *
      print *,'debut creation tableaux a11 a a13 source explosive'

! DK DK DK modification position pour source explo supprimee
      goto 765

!---- eviter source au bord d'un element
      if(ixgllsource .eq. 1) ixgllsource = 2
      if(ixgllsource .eq. nxgll) ixgllsource = nxgll - 1
      if(iygllsource .eq. 1) iygllsource = 2
      if(iygllsource .eq. nygll) iygllsource = nygll - 1
      if(izgllsource .eq. 1) izgllsource = 2
      if(izgllsource .eq. nzgll) izgllsource = nzgll - 1
      ipositsource = ibool(ixgllsource,iygllsource,izgllsource,ispecsource)

 765  continue

!---- definir a11, a12 et a13
      do k=1,nxgll
            do l=1,nygll
                  do m=1,nzgll
                        a11(k,l,m) = zero
                        a12(k,l,m) = zero
                        a13(k,l,m) = zero
                  enddo
            enddo
      enddo

      do k=1,nxgll
            do l=1,nygll
                  do m=1,nzgll

!---- dirac (schema en croix)

      if(l .eq. iygllsource .and. m .eq. izgllsource) then
         xixl = xixstock(ixgllsource,iygllsource,izgllsource,ispecsource)
         xiyl = xiystock(ixgllsource,iygllsource,izgllsource,ispecsource)
         xizl = xizstock(ixgllsource,iygllsource,izgllsource,ispecsource)
         hprimeval = hprime(k,ixgllsource)
         a11(k,l,m) = a11(k,l,m) +
     .            (sig11*xixl + sig12*xiyl + sig13*xizl)*hprimeval
         a12(k,l,m) = a12(k,l,m) +
     .            (sig21*xixl + sig22*xiyl + sig23*xizl)*hprimeval
         a13(k,l,m) = a13(k,l,m) +
     .            (sig31*xixl + sig32*xiyl + sig33*xizl)*hprimeval
      endif

!---- dirac (schema en croix)

      if(k .eq. ixgllsource .and. m .eq. izgllsource) then
         etaxl = etaxstock(ixgllsource,iygllsource,izgllsource,ispecsource)
         etayl = etaystock(ixgllsource,iygllsource,izgllsource,ispecsource)
         etazl = etazstock(ixgllsource,iygllsource,izgllsource,ispecsource)
         hprimeval = hprime(l,iygllsource)

         a11(k,l,m) = a11(k,l,m) +
     .            (sig11*etaxl + sig12*etayl + sig13*etazl)*hprimeval
         a12(k,l,m) = a12(k,l,m) +
     .            (sig21*etaxl + sig22*etayl + sig23*etazl)*hprimeval
         a13(k,l,m) = a13(k,l,m) +
     .            (sig31*etaxl + sig32*etayl + sig33*etazl)*hprimeval

!         a11(k,l,m) = a11(k,l,m) + sig0*etaxl*hprimeval
!         a12(k,l,m) = a12(k,l,m) + sig0*etayl*hprimeval
!         a13(k,l,m) = a13(k,l,m) + sig0*etazl*hprimeval

      endif

!---- dirac (schema en croix)

      if(k .eq. ixgllsource .and. l .eq. iygllsource) then
         gammaxl = gammaxstock(ixgllsource,iygllsource,izgllsource,ispecsource)
         gammayl = gammaystock(ixgllsource,iygllsource,izgllsource,ispecsource)
         gammazl = gammazstock(ixgllsource,iygllsource,izgllsource,ispecsource)
         hprimeval = hprime(m,izgllsource)

         a11(k,l,m) = a11(k,l,m) +
     .            (sig11*gammaxl + sig12*gammayl + sig13*gammazl)*hprimeval
         a12(k,l,m) = a12(k,l,m) +
     .            (sig21*gammaxl + sig22*gammayl + sig23*gammazl)*hprimeval
         a13(k,l,m) = a13(k,l,m) +
     .            (sig31*gammaxl + sig32*gammayl + sig33*gammazl)*hprimeval

!         a11(k,l,m) = a11(k,l,m) + sig0*gammaxl*hprimeval
!         a12(k,l,m) = a12(k,l,m) + sig0*gammayl*hprimeval
!         a13(k,l,m) = a13(k,l,m) + sig0*gammazl*hprimeval

      endif

                  enddo
            enddo
      enddo

c---- sauvegarde tableaux a11 a a13
      if(save_binary_files) then
      open(unit=27,file='a11_a13.bin',status='unknown',
     .            form='unformatted')
      write(27) a11
      write(27) a12
      write(27) a13
      close(27)
      endif

      print *,'fin creation tableaux a11 a a13 source explosive'

c %%%%%%%%%%%%%%%%%%

! sauvegarder la position de la source pour le solver
      open(unit=54,file='positsource.dat',status='unknown')
      write(54,*) ipositsource
      write(54,*) ixgllsource,iygllsource,izgllsource,ispecsource
      close(54)

c----
c---- recherche des points pour le display interpole (plans de coupe)
c----

      print *
      print *,' recherche indices points des plans de coupe'
      print *

!----
!----  position du plan de coupe X = cste
!----
      deltaxdisp = (ymax-ymin)/dble(nxdisp-1)
      deltazdisp = (0.  -zmin)/dble(nzdisp-1)

! DK DK DK plan de coupe au milieu du modele (avec origine en zero)
      valplancoupe = (xmax - xmin) / 2.d0
      valtoler = 20.d0

      print *
      print *,'Plan de coupe suivant X = cste'
      print *
      print *,'Position plan coupe valplancoupe = ',valplancoupe + xmin
      print *,'Tolerance absolue valtoler = ',valtoler
      print *,'Resolution display deltaxdisp = ',deltaxdisp
      print *,'Resolution display deltazdisp = ',deltazdisp
      print *

      do ispec = 1,nspec
      do k=1,nzgll
            do j=1,nygll
                  do i=1,nxgll

            xpoint = dble(xmeshstock(i,j,k,ispec))
            ypoint = dble(ymeshstock(i,j,k,ispec))
            zpoint = dble(zmeshstock(i,j,k,ispec))

! voir si le point est dans la bonne tranche
            if(xpoint .gt. valplancoupe - valtoler .and.
     .            xpoint .lt. valplancoupe + valtoler) then

! trouver indice de cette position physique dans le plan de coupe
                  indcutx = nint(ypoint / deltaxdisp) + 1
                  indcutz = nint(zpoint / deltazdisp) + 1
                  if(indcutx .lt. 1) indcutx = 1
                  if(indcutz .lt. 1) indcutz = 1
                  if(indcutx .gt. nxdisp) indcutx = nxdisp
                  if(indcutz .gt. nzdisp) indcutz = nzdisp
                  indicedispX(indcutx,indcutz) = ibool(i,j,k,ispec)

            endif

                  enddo
            enddo
      enddo
      enddo

! sauvegarder les indices des points du plan de coupe
      open(unit=54,file='pointscoupeX.dat',status='unknown')
      write(54,*) deltaxdisp,deltazdisp
      do i=1,nxdisp
      do j=1,nzdisp
            write(54,*) indicedispX(i,j)
      enddo
      enddo
      close(54)


!----
!----  position du plan de coupe Y = cste
!----
      deltaxdisp = (xmax-xmin)/dble(nxdisp-1)
      deltazdisp = (0.  -zmin)/dble(nzdisp-1)

! DK DK DK plan de coupe au milieu du modele (avec origine en zero)
      valplancoupe = (ymax - ymin) / 2.d0
      valtoler = 20.d0

      print *
      print *,'Plan de coupe suivant Y = cste'
      print *
      print *,'Position plan coupe valplancoupe = ',valplancoupe + ymin
      print *,'Tolerance absolue valtoler = ',valtoler
      print *,'Resolution display deltaxdisp = ',deltaxdisp
      print *,'Resolution display deltazdisp = ',deltazdisp
      print *

      do ispec = 1,nspec
      do k=1,nzgll
            do j=1,nygll
                  do i=1,nxgll

            xpoint = dble(xmeshstock(i,j,k,ispec))
            ypoint = dble(ymeshstock(i,j,k,ispec))
            zpoint = dble(zmeshstock(i,j,k,ispec))

! voir si le point est dans la bonne tranche
            if(ypoint .gt. valplancoupe - valtoler .and.
     .            ypoint .lt. valplancoupe + valtoler) then

! trouver indice de cette position physique dans le plan de coupe
                  indcutx = nint(xpoint / deltaxdisp) + 1
                  indcutz = nint(zpoint / deltazdisp) + 1
                  if(indcutx .lt. 1) indcutx = 1
                  if(indcutz .lt. 1) indcutz = 1
                  if(indcutx .gt. nxdisp) indcutx = nxdisp
                  if(indcutz .gt. nzdisp) indcutz = nzdisp
                  indicedispY(indcutx,indcutz) = ibool(i,j,k,ispec)

            endif

                  enddo
            enddo
      enddo
      enddo

! sauvegarder les indices des points du plan de coupe
      open(unit=54,file='pointscoupeY.dat',status='unknown')
      write(54,*) deltaxdisp,deltazdisp
      do i=1,nxdisp
      do j=1,nzdisp
            write(54,*) indicedispY(i,j)
      enddo
      enddo
      close(54)


!----
!----  position du plan de coupe Z = cste
!----
      deltaxdisp = (xmax-xmin)/dble(nxdisp-1)
      deltazdisp = (ymax-ymin)/dble(nzdisp-1)

! DK DK DK plan de coupe au quart du modele sous la surface (origine en zero)
      valplancoupe = 3.d0 * (zmax - zmin) / 4.d0
      valtoler = 80.d0

      print *
      print *,'Plan de coupe suivant Z = cste'
      print *
      print *,'Position plan coupe valplancoupe = ',valplancoupe + zmin
      print *,'Tolerance absolue valtoler = ',valtoler
      print *,'Resolution display deltaxdisp = ',deltaxdisp
      print *,'Resolution display deltazdisp = ',deltazdisp
      print *

      do ispec = 1,nspec
      do k=1,nzgll
            do j=1,nygll
                  do i=1,nxgll

            xpoint = dble(xmeshstock(i,j,k,ispec))
            ypoint = dble(ymeshstock(i,j,k,ispec))
            zpoint = dble(zmeshstock(i,j,k,ispec))

! voir si le point est dans la bonne tranche
            if(zpoint .gt. valplancoupe - valtoler .and.
     .            zpoint .lt. valplancoupe + valtoler) then

! trouver indice de cette position physique dans le plan de coupe
                  indcutx = nint(xpoint / deltaxdisp) + 1
                  indcutz = nint(ypoint / deltazdisp) + 1
                  if(indcutx .lt. 1) indcutx = 1
                  if(indcutz .lt. 1) indcutz = 1
                  if(indcutx .gt. nxdisp) indcutx = nxdisp
                  if(indcutz .gt. nzdisp) indcutz = nzdisp
                  indicedispZ(indcutx,indcutz) = ibool(i,j,k,ispec)

            endif

                  enddo
            enddo
      enddo
      enddo

! sauvegarder les indices des points du plan de coupe
      open(unit=54,file='pointscoupeZ.dat',status='unknown')
      write(54,*) deltaxdisp,deltazdisp
      do i=1,nxdisp
      do j=1,nzdisp
            write(54,*) indicedispZ(i,j)
      enddo
      enddo
      close(54)


c ***
c *** generer un fichier 'GNUPLOT' pour le controle de la grille ***
c ***

      if(ignuplot) then
      write(*,*)' Ecriture de la grille de surface format GNUPLOT...'

      open(unit=20,file='surfgrid.GNU',status='unknown')
      iz = nz
      do ix=0,nx
              do iy=0,ny
                write(20,15) x(ix,iy,iz),y(ix,iy,iz),z(ix,iy,iz)
                write(20,15) x(ix,iy,iz),y(ix,iy,iz),z(ix,iy,iz)
              enddo
              write(20,30)
      enddo
      close(20)

      open(unit=20,file='botgrid.GNU',status='unknown')
      iz = 0
      do ix=0,nx
              do iy=0,ny
                write(20,15) x(ix,iy,iz),y(ix,iy,iz),z(ix,iy,iz)
                write(20,15) x(ix,iy,iz),y(ix,iy,iz),z(ix,iy,iz)
              enddo
              write(20,30)
      enddo
      close(20)

      open(unit=20,file='plotgrid.gnu',status='unknown')
      write(20,*) 'set term x11'
      write(20,*) '# vue en coupe de face X=cste'
      write(20,*) '#set view 90,0'
      write(20,*) '# vue en coupe de face Y=cste'
      write(20,*) '#set view 90,90'
      write(20,*) '# vue classique'
      write(20,*) 'set view 60,30'
      write(20,*) 'set xlabel "X"'
      write(20,*) 'set ylabel "Y"'
      write(20,*) 'set zlabel "Z"'
      write(20,*) 'set title "Surface de la grille"'
      write(20,*) 'set parametric'
      write(20,10)
      write(20,*) 'pause -1 ''Hit any key...'''
      close(20)

 10   format('splot "botgrid.GNU" using 1:2:3 t '''' w l 2, ',
     .            '"surfgrid.GNU" using 1:2:3 t '''' w l 1')
 15   format(e12.5,1x,e12.5,1x,e12.5)
 30   format()

      write(*,*)' Fin ecriture de la grille format GNUPLOT'
      write(*,*)
      endif

c ***
c *** generer un fichier 'AVS' pour le controle de la grille ***
c ***

      if(draw_grid_avs) then

      write(*,*)' Ecriture de la grille 3D format AVS...'

      open(unit=20,file='AVSgrid3D.inp',status='unknown')

c nb de noeuds, de cellules, de donnees par cellule
      write(20,*) 8*nspec,nspec,1,0,0

c numero et coordonnees des noeuds
      inode = 0
      do ispec=1,nspec
      inode = inode + 1
      ix = 1
      iy = nygll
      iz = nzgll
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = 1
      iy = 1
      iz = nzgll
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = nxgll
      iy = 1
      iz = nzgll
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = nxgll
      iy = nygll
      iz = nzgll
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = 1
      iy = nygll
      iz = 1
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = 1
      iy = 1
      iz = 1
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = nxgll
      iy = 1
      iz = 1
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = nxgll
      iy = nygll
      iz = 1
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      enddo

c numero et coordonnees des cellules, definition materiau
      icell = 0
      do ispec=1,nspec
            write(20,300) ispec,numat(ispec),icell+1,icell+2,
     .  icell+3,icell+4,icell+5,icell+6,icell+7,icell+8
            icell = icell + 8
      enddo

c structure des donnees aux noeuds
      write(20,*) '1 1'
      write(20,*) 'fictifnode, bidonnode'

c donnees aux noeuds (hauteur du maillage)
      ix = 1
      iy = 1
      iz = 1
      inode = 0
      do ispec=1,nspec
      donnee = 255.d0*zmeshstock(ix,iy,iz,ispec)/(zmax-zmin)
      if(donnee.lt.1d0) donnee = 1.d0
      if(donnee.gt.255d0) donnee = 255.d0
      do inodloc = 1,ngnod
            inode = inode + 1
            write(20,210) inode,donnee
      enddo
      enddo

      close(20)

      write(*,*)' Fin ecriture de la grille 3D format AVS'
      write(*,*)' '

      endif

 200  format(i8,1x,e12.5,1x,e12.5,1x,e12.5)
 210  format(i8,1x,e12.5)
 300  format(i8,1x,i1,' hex ',i5,1x,i5,1x,i5,1x,i5,1x,i5,1x,i5,1x,i5,1x,i5)


c ***
c *** generer un fichier 'AVS' pour le controle du deraffinement ***
c ***

      write(*,*)' Ecriture de la grille de deraffinement format AVS...'

      open(unit=20,file='AVSgrid2Dderaf.inp',status='unknown')

c compter le nombre de faces a considerer (tester si sur le cote)
      nspecside = 0
      xtoler = 2.
      do ispec=1,nspec
      if(     xmeshstock(nxgll,1,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,nzgll,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,1,nzgll,ispec) .gt. (xmax - xmin) - xtoler)
     .                  nspecside = nspecside + 1
      enddo

c nb de noeuds, de cellules, de donnees par cellule
      write(20,*) 4*nspecside,nspecside,1,0,0

c numero et coordonnees des noeuds
      inode = 0
      do ispec=1,nspec

      if(     xmeshstock(nxgll,1,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,nzgll,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,1,nzgll,ispec) .gt. (xmax - xmin) - xtoler) then

      inode = inode + 1
      ix = nxgll
      iy = 1
      iz = 1
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = nxgll
      iy = nygll
      iz = 1
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = nxgll
      iy = nygll
      iz = nzgll
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)
      inode = inode + 1
      ix = nxgll
      iy = 1
      iz = nzgll
      write(20,200) inode,xmeshstock(ix,iy,iz,ispec),
     .      ymeshstock(ix,iy,iz,ispec),zmeshstock(ix,iy,iz,ispec)

      endif

      enddo

c numero et coordonnees des cellules, definition materiau
      icell = 0
      nspecside = 0
      do ispec=1,nspec

      if(     xmeshstock(nxgll,1,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,nzgll,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,1,nzgll,ispec) .gt. (xmax - xmin) - xtoler) then
          nspecside = nspecside + 1
          write(20,310) nspecside,numat(ispec),icell+1,icell+2,icell+3,icell+4
          icell = icell + 4
      endif

      enddo

c structure des donnees aux noeuds
      write(20,*) '1 1'
      write(20,*) 'fictifnode, bidonnode'

c donnees aux noeuds (rien pour l'instant)
      inode = 0
      ngnodside = 4
      do ispec=1,nspecside
      donnee = 1.d0
      do inodloc = 1,ngnodside
            inode = inode + 1
            write(20,210) inode,donnee
      enddo
      enddo

      close(20)

      write(*,*)' Fin ecriture de la grille de deraffinement format AVS'
      write(*,*)' '

 310  format(i8,1x,i1,' quad ',i5,1x,i5,1x,i5,1x,i5)

c ***
c *** generer les tableaux pour les bords absorbants
c ***

      if(bords_abs) then

      write(*,*)'Ecriture des tableaux pour les bords absorbants...'

c----
c---- definition des derivees des fonctions de forme 2D, element 4 noeuds
c----

       do l2 = 1,nygll

          t  = yigll(l2)

          do l1 = 1,nxgll

             s  = xigll(l1)

             sp = s + one
             sm = s - one
             tp = t + one
             tm = t - one

             dershape2D(1,1,l1,l2) = quart * tm
             dershape2D(1,2,l1,l2) = - quart * tm
             dershape2D(1,3,l1,l2) =  quart * tp
             dershape2D(1,4,l1,l2) = - quart * tp

             dershape2D(2,1,l1,l2) = quart * sm
             dershape2D(2,2,l1,l2) = - quart * sp
             dershape2D(2,3,l1,l2) =  quart * sp
             dershape2D(2,4,l1,l2) = - quart * sm

          enddo
       enddo

c----
c---- bord du fond
c----

c compter le nombre de faces a considerer (tester si sur le fond)
      nspecside = 0
      xtoler = 2.
      do ispec=1,nspec
      if(     zmeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. zmeshstock(nxgll,1,1,ispec) .lt. xtoler
     .  .and. zmeshstock(nxgll,nygll,1,ispec) .lt. xtoler
     .  .and. zmeshstock(1,nygll,1,ispec) .lt. xtoler)
     .                  nspecside = nspecside + 1
      enddo

c dans le cas du fond, on connait le nombre theorique de faces
      if(nspecside .ne. (nx/4)*(ny/4)) stop 'erreur bords absorbants fond'

      open(unit=20,file='bord_abs_fond.dat',status='unknown')

c nb de faces absorbantes
      write(20,*) nspecside

c numero des elements absorbants, et matrices de damping
      do ispec=1,nspec

      if(     zmeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. zmeshstock(nxgll,1,1,ispec) .lt. xtoler
     .  .and. zmeshstock(nxgll,nygll,1,ispec) .lt. xtoler
     .  .and. zmeshstock(1,nygll,1,ispec) .lt. xtoler) then

      write(20,*) ispec

c ***** definition des points de controle pour le jacobien 2D *******
      coorgx(1) = dble(xmeshstock(1,1,nzgll,ispec))
      coorgx(2) = dble(xmeshstock(nxgll,1,nzgll,ispec))
      coorgx(3) = dble(xmeshstock(nxgll,nygll,nzgll,ispec))
      coorgx(4) = dble(xmeshstock(1,nygll,nzgll,ispec))

      coorgz(1) = dble(ymeshstock(1,1,nzgll,ispec))
      coorgz(2) = dble(ymeshstock(nxgll,1,nzgll,ispec))
      coorgz(3) = dble(ymeshstock(nxgll,nygll,nzgll,ispec))
      coorgz(4) = dble(ymeshstock(1,nygll,nzgll,ispec))
c ***** definition des points de controle pour le jacobien 2D *******

            do iy=1,nygll
                  do ix=1,nxgll

c ****** DK DK DK calcul du jacobien 2D sur les bords ******
      ip1 = ix
      ip2 = iy

          xjac2_11 = zero
          xjac2_21 = zero
          xjac2_12 = zero
          xjac2_22 = zero

       do in = 1,4

      xjac2_11 = xjac2_11 + dershape2D(1,in,ip1,ip2)*coorgx(in)
      xjac2_21 = xjac2_21 + dershape2D(1,in,ip1,ip2)*coorgz(in)
      xjac2_12 = xjac2_12 + dershape2D(2,in,ip1,ip2)*coorgx(in)
      xjac2_22 = xjac2_22 + dershape2D(2,in,ip1,ip2)*coorgz(in)

       enddo

       rjacobloc = sngl(xjac2_11*xjac2_22 - xjac2_12*xjac2_21)

      if (rjacobloc .le. zero) stop 'Error : Jacobian undefined !!'

c print *,'jacob comput = ',rjacobloc

c ****** DK DK DK calcul du jacobien 2D sur les bords ******

      rjacobloc = ((xmax - xmin)/(2.*(nx/4)) * (ymax - ymin)/(2.*(ny/4)))

c print *,'jacob theor = ',rjacobloc

      dampP = rho(imoho) * vp(imoho) * sngl(wxgll(ix) * wxgll(iy)) * rjacobloc
      dampS = rho(imoho) * vs(imoho) * sngl(wxgll(ix) * wxgll(iy)) * rjacobloc

c DK DK DK DK mettre a zero les edges communs aux autres bords abs DK DK DK

c bord X=Xmin
      if(     xmeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. xmeshstock(1,nygll,1,ispec) .lt. xtoler
     .  .and. ix .eq. 1) then
            dampP = zero
            dampS = zero
      endif

c bord X=Xmax
      if(     xmeshstock(nxgll,1,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. ix .eq. nxgll) then
            dampP = zero
            dampS = zero
      endif

c bord Y=Ymin
      if(     ymeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. ymeshstock(nxgll,1,1,ispec) .lt. xtoler
     .  .and. iy .eq. 1) then
            dampP = zero
            dampS = zero
      endif

c bord Y=Ymax
      if(     ymeshstock(1,nygll,1,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. ymeshstock(nxgll,nygll,1,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. iy .eq. nygll) then
            dampP = zero
            dampS = zero
      endif

c DK DK DK DK mettre a zero les edges communs aux autres bords abs DK DK DK

      write(20,*) dampP,dampS

                  enddo
            enddo

      endif

      enddo

      close(20)

c----
c---- bord X = Xmin
c----

c compter le nombre de faces a considerer (tester si sur le bord)
      nspecside = 0
      xtoler = 2.
      do ispec=1,nspec
      if(     xmeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. xmeshstock(1,1,nzgll,ispec) .lt. xtoler
     .  .and. xmeshstock(1,nygll,1,ispec) .lt. xtoler
     .  .and. xmeshstock(1,nygll,nzgll,ispec) .lt. xtoler)
     .                  nspecside = nspecside + 1
      enddo

c afficher le nb de faces retenues
      print *,'Bord abs X = Xmin : ',nspecside,' faces'

      open(unit=20,file='bord_abs_X=Xmin.dat',status='unknown')

c nb de faces absorbantes
      write(20,*) nspecside

c numero des elements absorbants, et matrices de damping
      do ispec=1,nspec

      if(     xmeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. xmeshstock(1,1,nzgll,ispec) .lt. xtoler
     .  .and. xmeshstock(1,nygll,1,ispec) .lt. xtoler
     .  .and. xmeshstock(1,nygll,nzgll,ispec) .lt. xtoler) then

      write(20,*) ispec

c ***** definition des points de controle pour le jacobien 2D *******
      coorgx(1) = dble(ymeshstock(1,1,1,ispec))
      coorgx(2) = dble(ymeshstock(1,nygll,1,ispec))
      coorgx(3) = dble(ymeshstock(1,nygll,nzgll,ispec))
      coorgx(4) = dble(ymeshstock(1,1,nzgll,ispec))

      coorgz(1) = dble(zmeshstock(1,1,1,ispec))
      coorgz(2) = dble(zmeshstock(1,nygll,1,ispec))
      coorgz(3) = dble(zmeshstock(1,nygll,nzgll,ispec))
      coorgz(4) = dble(zmeshstock(1,1,nzgll,ispec))
c ***** definition des points de controle pour le jacobien 2D *******

            do iy=1,nygll
                  do iz=1,nzgll

c ****** DK DK DK calcul du jacobien 2D sur les bords ******
      ip1 = iy
      ip2 = iz

          xjac2_11 = zero
          xjac2_21 = zero
          xjac2_12 = zero
          xjac2_22 = zero

       do in = 1,4

      xjac2_11 = xjac2_11 + dershape2D(1,in,ip1,ip2)*coorgx(in)
      xjac2_21 = xjac2_21 + dershape2D(1,in,ip1,ip2)*coorgz(in)
      xjac2_12 = xjac2_12 + dershape2D(2,in,ip1,ip2)*coorgx(in)
      xjac2_22 = xjac2_22 + dershape2D(2,in,ip1,ip2)*coorgz(in)

       enddo

       rjacobloc = sngl(xjac2_11*xjac2_22 - xjac2_12*xjac2_21)

      if (rjacobloc .le. zero) stop 'Error : Jacobian undefined !!'

c ****** DK DK DK calcul du jacobien 2D sur les bords ******

c **** recherche des valeurs des parametres du modele de vitesse
      k=iz
      j=iy
      i=1

! recuperer le modele de vitesse P et S par l'indicateur de couche
      cp = vp(numat(ispec))
      cs = vs(numat(ispec))
      denst = rho(numat(ispec))

! inclure ce modele dans le modele 1D regional de Riki
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth12) then
            cp = vp(i12_25)
            cs = vs(i12_25)
            denst = rho(i12_25)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth25) then
            cp = vp(i25_moho)
            cs = vs(i25_moho)
            denst = rho(i25_moho)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depthmoho) then
            cp = vp(imoho)
            cs = vs(imoho)
            denst = rho(imoho)
      endif
c **** recherche des valeurs des parametres du modele de vitesse

      dampP = denst * cp * sngl(wxgll(iy) * wxgll(iz)) * rjacobloc
      dampS = denst * cs * sngl(wxgll(iy) * wxgll(iz)) * rjacobloc

      write(20,*) dampP,dampS

                  enddo
            enddo

      endif

      enddo

      close(20)

c----
c---- bord X = Xmax
c----

c compter le nombre de faces a considerer (tester si sur le bord)
      nspecside = 0
      xtoler = 2.
      do ispec=1,nspec
      if(     xmeshstock(nxgll,1,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,1,nzgll,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,nzgll,ispec) .gt. (xmax - xmin) - xtoler)
     .                  nspecside = nspecside + 1
      enddo

c afficher le nb de faces retenues
      print *,'Bord abs X = Xmax : ',nspecside,' faces'

      open(unit=20,file='bord_abs_X=Xmax.dat',status='unknown')

c nb de faces absorbantes
      write(20,*) nspecside

c numero des elements absorbants, et matrices de damping
      do ispec=1,nspec

      if(     xmeshstock(nxgll,1,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,1,nzgll,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,nzgll,ispec) .gt. (xmax - xmin) - xtoler)
     .      then

      write(20,*) ispec

c ***** definition des points de controle pour le jacobien 2D *******
      coorgx(1) = dble(ymeshstock(nxgll,1,1,ispec))
      coorgx(2) = dble(ymeshstock(nxgll,nygll,1,ispec))
      coorgx(3) = dble(ymeshstock(nxgll,nygll,nzgll,ispec))
      coorgx(4) = dble(ymeshstock(nxgll,1,nzgll,ispec))

      coorgz(1) = dble(zmeshstock(nxgll,1,1,ispec))
      coorgz(2) = dble(zmeshstock(nxgll,nygll,1,ispec))
      coorgz(3) = dble(zmeshstock(nxgll,nygll,nzgll,ispec))
      coorgz(4) = dble(zmeshstock(nxgll,1,nzgll,ispec))
c ***** definition des points de controle pour le jacobien 2D *******

            do iy=1,nygll
                  do iz=1,nzgll

c ****** DK DK DK calcul du jacobien 2D sur les bords ******
      ip1 = iy
      ip2 = iz

          xjac2_11 = zero
          xjac2_21 = zero
          xjac2_12 = zero
          xjac2_22 = zero

       do in = 1,4

      xjac2_11 = xjac2_11 + dershape2D(1,in,ip1,ip2)*coorgx(in)
      xjac2_21 = xjac2_21 + dershape2D(1,in,ip1,ip2)*coorgz(in)
      xjac2_12 = xjac2_12 + dershape2D(2,in,ip1,ip2)*coorgx(in)
      xjac2_22 = xjac2_22 + dershape2D(2,in,ip1,ip2)*coorgz(in)

       enddo

       rjacobloc = sngl(xjac2_11*xjac2_22 - xjac2_12*xjac2_21)

      if (rjacobloc .le. zero) stop 'Error : Jacobian undefined !!'

c ****** DK DK DK calcul du jacobien 2D sur les bords ******

c **** recherche des valeurs des parametres du modele de vitesse
      k=iz
      j=iy
      i=nxgll

! recuperer le modele de vitesse P et S par l'indicateur de couche
      cp = vp(numat(ispec))
      cs = vs(numat(ispec))
      denst = rho(numat(ispec))

! inclure ce modele dans le modele 1D regional de Riki
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth12) then
            cp = vp(i12_25)
            cs = vs(i12_25)
            denst = rho(i12_25)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth25) then
            cp = vp(i25_moho)
            cs = vs(i25_moho)
            denst = rho(i25_moho)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depthmoho) then
            cp = vp(imoho)
            cs = vs(imoho)
            denst = rho(imoho)
      endif
c **** recherche des valeurs des parametres du modele de vitesse

      dampP = denst * cp * sngl(wxgll(iy) * wxgll(iz)) * rjacobloc
      dampS = denst * cs * sngl(wxgll(iy) * wxgll(iz)) * rjacobloc

      write(20,*) dampP,dampS

                  enddo
            enddo

      endif

      enddo

      close(20)

c----
c---- bord Y = Ymin
c----

c compter le nombre de faces a considerer (tester si sur le bord)
      nspecside = 0
      xtoler = 2.
      do ispec=1,nspec
      if(     ymeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. ymeshstock(1,1,nzgll,ispec) .lt. xtoler
     .  .and. ymeshstock(nxgll,1,1,ispec) .lt. xtoler
     .  .and. ymeshstock(nxgll,1,nzgll,ispec) .lt. xtoler)
     .                  nspecside = nspecside + 1
      enddo

c afficher le nb de faces retenues
      print *,'Bord abs Y = Ymin : ',nspecside,' faces'

      open(unit=20,file='bord_abs_Y=Ymin.dat',status='unknown')

c nb de faces absorbantes
      write(20,*) nspecside

c numero des elements absorbants, et matrices de damping
      do ispec=1,nspec

      if(     ymeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. ymeshstock(1,1,nzgll,ispec) .lt. xtoler
     .  .and. ymeshstock(nxgll,1,1,ispec) .lt. xtoler
     .  .and. ymeshstock(nxgll,1,nzgll,ispec) .lt. xtoler) then

      write(20,*) ispec

c ***** definition des points de controle pour le jacobien 2D *******
      coorgx(1) = dble(xmeshstock(1,1,1,ispec))
      coorgx(2) = dble(xmeshstock(nxgll,1,1,ispec))
      coorgx(3) = dble(xmeshstock(nxgll,1,nzgll,ispec))
      coorgx(4) = dble(xmeshstock(1,1,nzgll,ispec))

      coorgz(1) = dble(zmeshstock(1,1,1,ispec))
      coorgz(2) = dble(zmeshstock(nxgll,1,1,ispec))
      coorgz(3) = dble(zmeshstock(nxgll,1,nzgll,ispec))
      coorgz(4) = dble(zmeshstock(1,1,nzgll,ispec))
c ***** definition des points de controle pour le jacobien 2D *******

            do ix=1,nxgll
                  do iz=1,nzgll

c ****** DK DK DK calcul du jacobien 2D sur les bords ******
      ip1 = ix
      ip2 = iz

          xjac2_11 = zero
          xjac2_21 = zero
          xjac2_12 = zero
          xjac2_22 = zero

       do in = 1,4

      xjac2_11 = xjac2_11 + dershape2D(1,in,ip1,ip2)*coorgx(in)
      xjac2_21 = xjac2_21 + dershape2D(1,in,ip1,ip2)*coorgz(in)
      xjac2_12 = xjac2_12 + dershape2D(2,in,ip1,ip2)*coorgx(in)
      xjac2_22 = xjac2_22 + dershape2D(2,in,ip1,ip2)*coorgz(in)

       enddo

       rjacobloc = sngl(xjac2_11*xjac2_22 - xjac2_12*xjac2_21)

      if (rjacobloc .le. zero) stop 'Error : Jacobian undefined !!'

c ****** DK DK DK calcul du jacobien 2D sur les bords ******

c **** recherche des valeurs des parametres du modele de vitesse
      k=iz
      j=1
      i=ix

! recuperer le modele de vitesse P et S par l'indicateur de couche
      cp = vp(numat(ispec))
      cs = vs(numat(ispec))
      denst = rho(numat(ispec))

! inclure ce modele dans le modele 1D regional de Riki
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth12) then
            cp = vp(i12_25)
            cs = vs(i12_25)
            denst = rho(i12_25)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth25) then
            cp = vp(i25_moho)
            cs = vs(i25_moho)
            denst = rho(i25_moho)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depthmoho) then
            cp = vp(imoho)
            cs = vs(imoho)
            denst = rho(imoho)
      endif
c **** recherche des valeurs des parametres du modele de vitesse

      dampP = denst * cp * sngl(wxgll(ix) * wxgll(iz)) * rjacobloc
      dampS = denst * cs * sngl(wxgll(ix) * wxgll(iz)) * rjacobloc

c DK DK DK DK mettre a zero les edges communs aux autres bords abs DK DK DK

c bord X=Xmin
      if(     xmeshstock(1,1,1,ispec) .lt. xtoler
     .  .and. xmeshstock(1,1,nzgll,ispec) .lt. xtoler
     .  .and. ix .eq. 1) then
            dampP = zero
            dampS = zero
      endif

c bord X=Xmax
      if(     xmeshstock(nxgll,1,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,1,nzgll,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. ix .eq. nxgll) then
            dampP = zero
            dampS = zero
      endif

c DK DK DK DK mettre a zero les edges communs aux autres bords abs DK DK DK

      write(20,*) dampP,dampS

                  enddo
            enddo

      endif

      enddo

      close(20)

c----
c---- bord Y = Ymax
c----

c compter le nombre de faces a considerer (tester si sur le bord)
      nspecside = 0
      xtoler = 2.
      do ispec=1,nspec
      if(     ymeshstock(1,nygll,1,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. ymeshstock(1,nygll,nzgll,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. ymeshstock(nxgll,nygll,1,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. ymeshstock(nxgll,nygll,nzgll,ispec) .gt. (ymax - ymin) - xtoler)
     .                  nspecside = nspecside + 1
      enddo

c afficher le nb de faces retenues
      print *,'Bord abs Y = Ymax : ',nspecside,' faces'

      open(unit=20,file='bord_abs_Y=Ymax.dat',status='unknown')

c nb de faces absorbantes
      write(20,*) nspecside

c numero des elements absorbants, et matrices de damping
      do ispec=1,nspec

      if(     ymeshstock(1,nygll,1,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. ymeshstock(1,nygll,nzgll,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. ymeshstock(nxgll,nygll,1,ispec) .gt. (ymax - ymin) - xtoler
     .  .and. ymeshstock(nxgll,nygll,nzgll,ispec) .gt. (ymax - ymin) - xtoler)
     .      then

      write(20,*) ispec

c ***** definition des points de controle pour le jacobien 2D *******
      coorgx(1) = dble(xmeshstock(1,nygll,1,ispec))
      coorgx(2) = dble(xmeshstock(nxgll,nygll,1,ispec))
      coorgx(3) = dble(xmeshstock(nxgll,nygll,nzgll,ispec))
      coorgx(4) = dble(xmeshstock(1,nygll,nzgll,ispec))

      coorgz(1) = dble(zmeshstock(1,nygll,1,ispec))
      coorgz(2) = dble(zmeshstock(nxgll,nygll,1,ispec))
      coorgz(3) = dble(zmeshstock(nxgll,nygll,nzgll,ispec))
      coorgz(4) = dble(zmeshstock(1,nygll,nzgll,ispec))
c ***** definition des points de controle pour le jacobien 2D *******

            do ix=1,nxgll
                  do iz=1,nzgll

c ****** DK DK DK calcul du jacobien 2D sur les bords ******
      ip1 = ix
      ip2 = iz

          xjac2_11 = zero
          xjac2_21 = zero
          xjac2_12 = zero
          xjac2_22 = zero

       do in = 1,4

      xjac2_11 = xjac2_11 + dershape2D(1,in,ip1,ip2)*coorgx(in)
      xjac2_21 = xjac2_21 + dershape2D(1,in,ip1,ip2)*coorgz(in)
      xjac2_12 = xjac2_12 + dershape2D(2,in,ip1,ip2)*coorgx(in)
      xjac2_22 = xjac2_22 + dershape2D(2,in,ip1,ip2)*coorgz(in)

       enddo

       rjacobloc = sngl(xjac2_11*xjac2_22 - xjac2_12*xjac2_21)

      if (rjacobloc .le. zero) stop 'Error : Jacobian undefined !!'

c ****** DK DK DK calcul du jacobien 2D sur les bords ******

c **** recherche des valeurs des parametres du modele de vitesse
      k=iz
      j=nygll
      i=ix

! recuperer le modele de vitesse P et S par l'indicateur de couche
      cp = vp(numat(ispec))
      cs = vs(numat(ispec))
      denst = rho(numat(ispec))

! inclure ce modele dans le modele 1D regional de Riki
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth12) then
            cp = vp(i12_25)
            cs = vs(i12_25)
            denst = rho(i12_25)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depth25) then
            cp = vp(i25_moho)
            cs = vs(i25_moho)
            denst = rho(i25_moho)
      endif
      if((dble(zmeshstock(i,j,k,ispec)) + zmin) .lt. depthmoho) then
            cp = vp(imoho)
            cs = vs(imoho)
            denst = rho(imoho)
      endif
c **** recherche des valeurs des parametres du modele de vitesse

      dampP = denst * cp * sngl(wxgll(ix) * wxgll(iz)) * rjacobloc
      dampS = denst * cs * sngl(wxgll(ix) * wxgll(iz)) * rjacobloc

c DK DK DK DK mettre a zero les edges communs aux autres bords abs DK DK DK

c bord X=Xmin
      if(     xmeshstock(1,nygll,1,ispec) .lt. xtoler
     .  .and. xmeshstock(1,nygll,nzgll,ispec) .lt. xtoler
     .  .and. ix .eq. 1) then
            dampP = zero
            dampS = zero
      endif

c bord X=Xmax
      if(     xmeshstock(nxgll,nygll,1,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. xmeshstock(nxgll,nygll,nzgll,ispec) .gt. (xmax - xmin) - xtoler
     .  .and. ix .eq. nxgll) then
            dampP = zero
            dampS = zero
      endif

c DK DK DK DK mettre a zero les edges communs aux autres bords abs DK DK DK

      write(20,*) dampP,dampS

                  enddo
            enddo

      endif

      enddo

      close(20)

      write(*,*)'Fin ecriture des tableaux pour les bords absorbants...'
      write(*,*)

      endif

c ***
c *** generer la liste des elements pour les plans de coupe AVS
c ***

      write(*,*)'Recherche des plans de coupe pour AVS...'

c coupe verticale X = cste

c compter le nombre de faces a considerer (tester si sur le cote)
      nspecside = 0
      xtoler = 2.
      xcutpos = (xmax - xmin) / 2.
      do ispec=1,nspec
      if(     abs(xmeshstock(nxgll,1,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(xmeshstock(nxgll,nygll,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(xmeshstock(nxgll,nygll,nzgll,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(xmeshstock(nxgll,1,nzgll,ispec) - xcutpos) .lt. xtoler)
     .                  nspecside = nspecside + 1
      enddo

c on densifie par quatre le display dans ce plan vertical
      open(unit=20,file='AVScutplaneXcste.dat',status='unknown')
      write(20,*) 4*nspecside
      do ispec=1,nspec

      if(     abs(xmeshstock(nxgll,1,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(xmeshstock(nxgll,nygll,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(xmeshstock(nxgll,nygll,nzgll,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(xmeshstock(nxgll,1,nzgll,ispec) - xcutpos) .lt. xtoler) then

            nyd = 1
            nyf = nygll/2
            nzd = 1
            nzf = nzgll/2

                 write(20,*) ispec,
     . ibool(nxgll,nyd,nzd,ispec),ibool(nxgll,nyf,nzd,ispec),
     . ibool(nxgll,nyf,nzf,ispec),ibool(nxgll,nyd,nzf,ispec),
     . xmeshstock(nxgll,nyd,nzd,ispec),xmeshstock(nxgll,nyf,nzd,ispec),
     . xmeshstock(nxgll,nyf,nzf,ispec),xmeshstock(nxgll,nyd,nzf,ispec),
     . ymeshstock(nxgll,nyd,nzd,ispec),ymeshstock(nxgll,nyf,nzd,ispec),
     . ymeshstock(nxgll,nyf,nzf,ispec),ymeshstock(nxgll,nyd,nzf,ispec),
     . zmeshstock(nxgll,nyd,nzd,ispec),zmeshstock(nxgll,nyf,nzd,ispec),
     . zmeshstock(nxgll,nyf,nzf,ispec),zmeshstock(nxgll,nyd,nzf,ispec),
     . rlambdastock(nxgll,nyd,nzd,ispec),rlambdastock(nxgll,nyf,nzd,ispec),
     . rlambdastock(nxgll,nyf,nzf,ispec),rlambdastock(nxgll,nyd,nzf,ispec)

            nyd = nygll/2
            nyf = nygll
            nzd = 1
            nzf = nzgll/2

                 write(20,*) ispec,
     . ibool(nxgll,nyd,nzd,ispec),ibool(nxgll,nyf,nzd,ispec),
     . ibool(nxgll,nyf,nzf,ispec),ibool(nxgll,nyd,nzf,ispec),
     . xmeshstock(nxgll,nyd,nzd,ispec),xmeshstock(nxgll,nyf,nzd,ispec),
     . xmeshstock(nxgll,nyf,nzf,ispec),xmeshstock(nxgll,nyd,nzf,ispec),
     . ymeshstock(nxgll,nyd,nzd,ispec),ymeshstock(nxgll,nyf,nzd,ispec),
     . ymeshstock(nxgll,nyf,nzf,ispec),ymeshstock(nxgll,nyd,nzf,ispec),
     . zmeshstock(nxgll,nyd,nzd,ispec),zmeshstock(nxgll,nyf,nzd,ispec),
     . zmeshstock(nxgll,nyf,nzf,ispec),zmeshstock(nxgll,nyd,nzf,ispec),
     . rlambdastock(nxgll,nyd,nzd,ispec),rlambdastock(nxgll,nyf,nzd,ispec),
     . rlambdastock(nxgll,nyf,nzf,ispec),rlambdastock(nxgll,nyd,nzf,ispec)

            nyd = 1
            nyf = nygll/2
            nzd = nzgll/2
            nzf = nzgll

                 write(20,*) ispec,
     . ibool(nxgll,nyd,nzd,ispec),ibool(nxgll,nyf,nzd,ispec),
     . ibool(nxgll,nyf,nzf,ispec),ibool(nxgll,nyd,nzf,ispec),
     . xmeshstock(nxgll,nyd,nzd,ispec),xmeshstock(nxgll,nyf,nzd,ispec),
     . xmeshstock(nxgll,nyf,nzf,ispec),xmeshstock(nxgll,nyd,nzf,ispec),
     . ymeshstock(nxgll,nyd,nzd,ispec),ymeshstock(nxgll,nyf,nzd,ispec),
     . ymeshstock(nxgll,nyf,nzf,ispec),ymeshstock(nxgll,nyd,nzf,ispec),
     . zmeshstock(nxgll,nyd,nzd,ispec),zmeshstock(nxgll,nyf,nzd,ispec),
     . zmeshstock(nxgll,nyf,nzf,ispec),zmeshstock(nxgll,nyd,nzf,ispec),
     . rlambdastock(nxgll,nyd,nzd,ispec),rlambdastock(nxgll,nyf,nzd,ispec),
     . rlambdastock(nxgll,nyf,nzf,ispec),rlambdastock(nxgll,nyd,nzf,ispec)

            nyd = nygll/2
            nyf = nygll
            nzd = nzgll/2
            nzf = nzgll

                 write(20,*) ispec,
     . ibool(nxgll,nyd,nzd,ispec),ibool(nxgll,nyf,nzd,ispec),
     . ibool(nxgll,nyf,nzf,ispec),ibool(nxgll,nyd,nzf,ispec),
     . xmeshstock(nxgll,nyd,nzd,ispec),xmeshstock(nxgll,nyf,nzd,ispec),
     . xmeshstock(nxgll,nyf,nzf,ispec),xmeshstock(nxgll,nyd,nzf,ispec),
     . ymeshstock(nxgll,nyd,nzd,ispec),ymeshstock(nxgll,nyf,nzd,ispec),
     . ymeshstock(nxgll,nyf,nzf,ispec),ymeshstock(nxgll,nyd,nzf,ispec),
     . zmeshstock(nxgll,nyd,nzd,ispec),zmeshstock(nxgll,nyf,nzd,ispec),
     . zmeshstock(nxgll,nyf,nzf,ispec),zmeshstock(nxgll,nyd,nzf,ispec),
     . rlambdastock(nxgll,nyd,nzd,ispec),rlambdastock(nxgll,nyf,nzd,ispec),
     . rlambdastock(nxgll,nyf,nzf,ispec),rlambdastock(nxgll,nyd,nzf,ispec)

      endif

      enddo
      close(20)

c coupe verticale Y = cste

c compter le nombre de faces a considerer (tester si sur le cote)
      nspecside = 0
      xtoler = 2.
      xcutpos = (ymax - ymin) / 2.
      do ispec=1,nspec
      if(     abs(ymeshstock(1,nygll,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(ymeshstock(nxgll,nygll,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(ymeshstock(nxgll,nygll,nzgll,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(ymeshstock(1,nygll,nzgll,ispec) - xcutpos) .lt. xtoler)
     .                  nspecside = nspecside + 1
      enddo

c on densifie par quatre le display dans ce plan vertical
      open(unit=20,file='AVScutplaneYcste.dat',status='unknown')
      write(20,*) 4*nspecside
      do ispec=1,nspec

      if(     abs(ymeshstock(1,nygll,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(ymeshstock(nxgll,nygll,1,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(ymeshstock(nxgll,nygll,nzgll,ispec) - xcutpos) .lt. xtoler
     .  .and. abs(ymeshstock(1,nygll,nzgll,ispec) - xcutpos) .lt. xtoler) then

            nxd = 1
            nxf = nxgll/2
            nzd = 1
            nzf = nzgll/2

          write(20,*) ispec,
     .  ibool(nxd,nygll,nzd,ispec),ibool(nxf,nygll,nzd,ispec),
     .  ibool(nxf,nygll,nzf,ispec),ibool(nxd,nygll,nzf,ispec),
     .  xmeshstock(nxd,nygll,nzd,ispec),xmeshstock(nxf,nygll,nzd,ispec),
     .  xmeshstock(nxf,nygll,nzf,ispec),xmeshstock(nxd,nygll,nzf,ispec),
     .  ymeshstock(nxd,nygll,nzd,ispec),ymeshstock(nxf,nygll,nzd,ispec),
     .  ymeshstock(nxf,nygll,nzf,ispec),ymeshstock(nxd,nygll,nzf,ispec),
     .  zmeshstock(nxd,nygll,nzd,ispec),zmeshstock(nxf,nygll,nzd,ispec),
     .  zmeshstock(nxf,nygll,nzf,ispec),zmeshstock(nxd,nygll,nzf,ispec),
     .  rlambdastock(nxd,nygll,nzd,ispec),rlambdastock(nxf,nygll,nzd,ispec),
     .  rlambdastock(nxf,nygll,nzf,ispec),rlambdastock(nxd,nygll,nzf,ispec)

            nxd = nxgll / 2
            nxf = nxgll
            nzd = 1
            nzf = nzgll / 2

          write(20,*) ispec,
     .  ibool(nxd,nygll,nzd,ispec),ibool(nxf,nygll,nzd,ispec),
     .  ibool(nxf,nygll,nzf,ispec),ibool(nxd,nygll,nzf,ispec),
     .  xmeshstock(nxd,nygll,nzd,ispec),xmeshstock(nxf,nygll,nzd,ispec),
     .  xmeshstock(nxf,nygll,nzf,ispec),xmeshstock(nxd,nygll,nzf,ispec),
     .  ymeshstock(nxd,nygll,nzd,ispec),ymeshstock(nxf,nygll,nzd,ispec),
     .  ymeshstock(nxf,nygll,nzf,ispec),ymeshstock(nxd,nygll,nzf,ispec),
     .  zmeshstock(nxd,nygll,nzd,ispec),zmeshstock(nxf,nygll,nzd,ispec),
     .  zmeshstock(nxf,nygll,nzf,ispec),zmeshstock(nxd,nygll,nzf,ispec),
     .  rlambdastock(nxd,nygll,nzd,ispec),rlambdastock(nxf,nygll,nzd,ispec),
     .  rlambdastock(nxf,nygll,nzf,ispec),rlambdastock(nxd,nygll,nzf,ispec)

            nxd = 1
            nxf = nxgll / 2
            nzd = nzgll / 2
            nzf = nzgll

          write(20,*) ispec,
     .  ibool(nxd,nygll,nzd,ispec),ibool(nxf,nygll,nzd,ispec),
     .  ibool(nxf,nygll,nzf,ispec),ibool(nxd,nygll,nzf,ispec),
     .  xmeshstock(nxd,nygll,nzd,ispec),xmeshstock(nxf,nygll,nzd,ispec),
     .  xmeshstock(nxf,nygll,nzf,ispec),xmeshstock(nxd,nygll,nzf,ispec),
     .  ymeshstock(nxd,nygll,nzd,ispec),ymeshstock(nxf,nygll,nzd,ispec),
     .  ymeshstock(nxf,nygll,nzf,ispec),ymeshstock(nxd,nygll,nzf,ispec),
     .  zmeshstock(nxd,nygll,nzd,ispec),zmeshstock(nxf,nygll,nzd,ispec),
     .  zmeshstock(nxf,nygll,nzf,ispec),zmeshstock(nxd,nygll,nzf,ispec),
     .  rlambdastock(nxd,nygll,nzd,ispec),rlambdastock(nxf,nygll,nzd,ispec),
     .  rlambdastock(nxf,nygll,nzf,ispec),rlambdastock(nxd,nygll,nzf,ispec)

            nxd = nxgll / 2
            nxf = nxgll
            nzd = nzgll / 2
            nzf = nzgll

          write(20,*) ispec,
     .  ibool(nxd,nygll,nzd,ispec),ibool(nxf,nygll,nzd,ispec),
     .  ibool(nxf,nygll,nzf,ispec),ibool(nxd,nygll,nzf,ispec),
     .  xmeshstock(nxd,nygll,nzd,ispec),xmeshstock(nxf,nygll,nzd,ispec),
     .  xmeshstock(nxf,nygll,nzf,ispec),xmeshstock(nxd,nygll,nzf,ispec),
     .  ymeshstock(nxd,nygll,nzd,ispec),ymeshstock(nxf,nygll,nzd,ispec),
     .  ymeshstock(nxf,nygll,nzf,ispec),ymeshstock(nxd,nygll,nzf,ispec),
     .  zmeshstock(nxd,nygll,nzd,ispec),zmeshstock(nxf,nygll,nzd,ispec),
     .  zmeshstock(nxf,nygll,nzf,ispec),zmeshstock(nxd,nygll,nzf,ispec),
     .  rlambdastock(nxd,nygll,nzd,ispec),rlambdastock(nxf,nygll,nzd,ispec),
     .  rlambdastock(nxf,nygll,nzf,ispec),rlambdastock(nxd,nygll,nzf,ispec)

      endif

      enddo
      close(20)

c display surface Z = Ztopo

c nombre de faces a considerer connu a l'avance (tous les elements de surface)
      nspecside = nx * ny

      open(unit=20,file='AVScutplaneZcste.dat',status='unknown')
      write(20,*) nspecside

c parcourir la couche sedimentaire dans le meme ordre que pour la generation
      ispec = 0
      do iz=nz1,nz-1
      do ix=0,nx-1
      do iy=0,ny-1

            ispec = ispec + 1

c garder les elements de la couche superieure
      if(iz .eq. nz-1)
     .                  write(20,*) ispec,
     .  ibool(1,1,nzgll,ispec),ibool(nxgll,1,nzgll,ispec),
     .  ibool(nxgll,nygll,nzgll,ispec),ibool(1,nygll,nzgll,ispec),
     .  xmeshstock(1,1,nzgll,ispec),xmeshstock(nxgll,1,nzgll,ispec),
     .  xmeshstock(nxgll,nygll,nzgll,ispec),xmeshstock(1,nygll,nzgll,ispec),
     .  ymeshstock(1,1,nzgll,ispec),ymeshstock(nxgll,1,nzgll,ispec),
     .  ymeshstock(nxgll,nygll,nzgll,ispec),ymeshstock(1,nygll,nzgll,ispec),
     .  zmeshstock(1,1,nzgll,ispec),zmeshstock(nxgll,1,nzgll,ispec),
     .  zmeshstock(nxgll,nygll,nzgll,ispec),zmeshstock(1,nygll,nzgll,ispec),
     .  rlambdastock(1,1,nzgll,ispec),rlambdastock(nxgll,1,nzgll,ispec),
     .  rlambdastock(nxgll,nygll,nzgll,ispec),rlambdastock(1,nygll,nzgll,ispec)

      enddo
      enddo
      enddo
      close(20)

c ***
c *** generer la liste des points pour l'enregistrement des sismogrammes
c ***

      write(*,*)'Recherche des positions des recepteurs...'

c----
c coupe verticale X = cste
c----

c nombre de points a considerer connu
      nspecside = ny

c numero d'element ou l'on place la ligne
      numelemligne = nx / 2

      open(unit=20,file='receiversXcste.dat',status='unknown')
      write(20,*) nspecside

c parcourir la couche sedimentaire dans le meme ordre que pour la generation
      ispec = 0
      do iz=nz1,nz-1
      do ix=0,nx-1
      do iy=0,ny-1

            ispec = ispec + 1

c garder les elements de la couche superieure le long du profil
      if(iz .eq. nz-1 .and. (ix+1) .eq. numelemligne)
     .           write(20,*) ibool(nxgll,1,nzgll,ispec),
     . xmeshstock(nxgll,1,nzgll,ispec),
     . ymeshstock(nxgll,1,nzgll,ispec),
     . zmeshstock(nxgll,1,nzgll,ispec)

      enddo
      enddo
      enddo
      close(20)

c----
c coupe verticale Y = cste
c----

c nombre de points a considerer connu
      nspecside = nx

c numero d'element ou l'on place la ligne
      numelemligne = ny / 2

      open(unit=20,file='receiversYcste.dat',status='unknown')
      write(20,*) nspecside

c parcourir la couche sedimentaire dans le meme ordre que pour la generation
      ispec = 0
      do iz=nz1,nz-1
      do ix=0,nx-1
      do iy=0,ny-1

            ispec = ispec + 1

c garder les elements de la couche superieure le long du profil
      if(iz .eq. nz-1 .and. (iy+1) .eq. numelemligne)
     .         write(20,*) ibool(1,nygll,nzgll,ispec),
     .  xmeshstock(1,nygll,nzgll,ispec),
     .  ymeshstock(1,nygll,nzgll,ispec),
     .  zmeshstock(1,nygll,nzgll,ispec)

      enddo
      enddo
      enddo
      close(20)

      call system('date')

      end


c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c calcul de la matrice jacobienne a partir des fonctions de forme

      subroutine calcjac(dvolustock,xixstock,xiystock,etaxstock,
     .   etaystock,xizstock,etazstock,
     .   gammaxstock,gammaystock,gammazstock,
     .   xmeshstock,ymeshstock,zmeshstock,
     .   xmesh,ymesh,zmesh,gammax,gammay,gammaz,dvolu,
     .   xcell,ycell,zcell,shape,dershape,ngnod,ndime,nxgll,nygll,nzgll,
     .   ispec,nspec,xmin,ymin,zmin,
     .   xixmin,xixmax,xiymin,xiymax,xizmin,xizmax,
     .   etaxmin,etaxmax,etaymin,etaymax,etazmin,etazmax,
     .   gammaxmin,gammaxmax,gammaymin,gammaymax,gammazmin,gammazmax,
     .   xximin,xximax,yximin,yximax,zximin,zximax,
     .   xetamin,xetamax,yetamin,yetamax,zetamin,zetamax,
     .   xgammamin,xgammamax,ygammamin,ygammamax,zgammamin,zgammamax)

      implicit none

      double precision
     .      xiymin,xiymax,xizmin,xizmax,etaxmin,etaxmax,etazmin,etazmax,
     .   xixmin,xixmax,etaymin,etaymax,gammaxmin,gammaxmax,
     .   gammaymin,gammaymax,gammazmin,gammazmax
      double precision
     .      yximin,yximax,zximin,zximax,xetamin,xetamax,zetamin,zetamax,
     .   xximin,xximax,yetamin,yetamax,xgammamin,xgammamax,
     .   ygammamin,ygammamax,zgammamin,zgammamax

      integer ngnod,ndime,nxgll,nygll,nzgll,i,j,k,nspec,ispec
      double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
      double precision xmesh,ymesh,zmesh
      double precision dvolu,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
      double precision xmin,ymin,zmin

      double precision shape(ngnod,nxgll,nygll,nzgll)
      double precision dershape(ndime,ngnod,nxgll,nygll,nzgll)

      double precision xcell(ngnod)
      double precision ycell(ngnod)
      double precision zcell(ngnod)

      real dvolustock(nxgll,nxgll,nxgll,nspec)
      real xixstock(nxgll,nxgll,nxgll,nspec)
      real xiystock(nxgll,nxgll,nxgll,nspec)
      real xizstock(nxgll,nxgll,nxgll,nspec)
      real etaxstock(nxgll,nxgll,nxgll,nspec)
      real etaystock(nxgll,nxgll,nxgll,nspec)
      real etazstock(nxgll,nxgll,nxgll,nspec)
      real gammaxstock(nxgll,nxgll,nxgll,nspec)
      real gammaystock(nxgll,nxgll,nxgll,nspec)
      real gammazstock(nxgll,nxgll,nxgll,nspec)
      real xmeshstock(nxgll,nxgll,nxgll,nspec)
      real ymeshstock(nxgll,nxgll,nxgll,nspec)
      real zmeshstock(nxgll,nxgll,nxgll,nspec)

      integer ia

      double precision zero
      parameter(zero = 0.d0)

! incrementation du numero d'element spectral courant
      ispec = ispec + 1
      if(ispec .gt. nspec) stop 'improper spectral element number'

! calcul pour tous les noeuds interieurs a l'element spectral
      do k=1,nzgll
      do j=1,nygll
      do i=1,nxgll

      xxi = zero
      xeta = zero
      xgamma = zero
      yxi = zero
      yeta = zero
      ygamma = zero
      zxi = zero
      zeta = zero
      zgamma = zero
      dvolu = zero
      xmesh = zero
      ymesh = zero
      zmesh = zero

      do ia=1,ngnod
            xxi    = xxi    + dershape(1,ia,i,j,k)*xcell(ia)
            xeta   = xeta   + dershape(2,ia,i,j,k)*xcell(ia)
            xgamma = xgamma + dershape(3,ia,i,j,k)*xcell(ia)

            yxi    = yxi    + dershape(1,ia,i,j,k)*ycell(ia)
            yeta   = yeta   + dershape(2,ia,i,j,k)*ycell(ia)
            ygamma = ygamma + dershape(3,ia,i,j,k)*ycell(ia)

            zxi    = zxi    + dershape(1,ia,i,j,k)*zcell(ia)
            zeta   = zeta   + dershape(2,ia,i,j,k)*zcell(ia)
            zgamma = zgamma + dershape(3,ia,i,j,k)*zcell(ia)

            xmesh  = xmesh  + shape(ia,i,j,k)*xcell(ia)
            ymesh  = ymesh  + shape(ia,i,j,k)*ycell(ia)
            zmesh  = zmesh  + shape(ia,i,j,k)*zcell(ia)
      enddo

      dvolu = xxi*(yeta*zgamma - ygamma*zeta) -
     .        xeta*(yxi*zgamma - ygamma*zxi) +
     .        xgamma*(yxi*zeta - yeta*zxi)

      if(dvolu.le.zero) stop 'Jacobian undefined !!'

c inverser la transformation (Fletcher p. 50 tome 2)
      xix = (yeta*zgamma - ygamma*zeta) / dvolu
      xiy = (xgamma*zeta - xeta*zgamma) / dvolu
      xiz = (xeta*ygamma - xgamma*yeta) / dvolu

      etax = (ygamma*zxi - yxi*zgamma) / dvolu
      etay = (xxi*zgamma - xgamma*zxi) / dvolu
      etaz = (xgamma*yxi - xxi*ygamma) / dvolu

      gammax = (yxi*zeta - yeta*zxi) / dvolu
      gammay = (xeta*zxi - xxi*zeta) / dvolu
      gammaz = (xxi*yeta - xeta*yxi) / dvolu

c verifier les valeurs min et max de certains termes
      xximin  = dmin1(xximin,xxi)
      yximin  = dmin1(yximin,yxi)
      zximin  = dmin1(zximin,zxi)
      xetamin = dmin1(xetamin,xeta)
      yetamin = dmin1(yetamin,yeta)
      zetamin = dmin1(zetamin,zeta)
      xgammamin = dmin1(xgammamin,xgamma)
      ygammamin = dmin1(ygammamin,ygamma)
      zgammamin = dmin1(zgammamin,zgamma)
      xximax  = dmax1(xximax,xxi)
      yximax  = dmax1(yximax,yxi)
      zximax  = dmax1(zximax,zxi)
      xetamax = dmax1(xetamax,xeta)
      yetamax = dmax1(yetamax,yeta)
      zetamax = dmax1(zetamax,zeta)
      xgammamax = dmax1(xgammamax,xgamma)
      ygammamax = dmax1(ygammamax,ygamma)
      zgammamax = dmax1(zgammamax,zgamma)

      xixmin  = dmin1(xixmin,xix)
      xiymin  = dmin1(xiymin,xiy)
      xizmin  = dmin1(xizmin,xiz)
      etaxmin = dmin1(etaxmin,etax)
      etaymin = dmin1(etaymin,etay)
      etazmin = dmin1(etazmin,etaz)
      gammaxmin = dmin1(gammaxmin,gammax)
      gammaymin = dmin1(gammaymin,gammay)
      gammazmin = dmin1(gammazmin,gammaz)
      xixmax  = dmax1(xixmax,xix)
      xiymax  = dmax1(xiymax,xiy)
      xizmax  = dmax1(xizmax,xiz)
      etaxmax = dmax1(etaxmax,etax)
      etaymax = dmax1(etaymax,etay)
      etazmax = dmax1(etazmax,etaz)
      gammaxmax = dmax1(gammaxmax,gammax)
      gammaymax = dmax1(gammaymax,gammay)
      gammazmax = dmax1(gammazmax,gammaz)

c stocker pour sauvegarde
c on rapporte egalement l'origine des coordonnees en zero
      dvolustock(i,j,k,ispec) = sngl(dvolu)
      xixstock(i,j,k,ispec) = sngl(xix)
      xiystock(i,j,k,ispec) = sngl(xiy)
      xizstock(i,j,k,ispec) = sngl(xiz)
      etaxstock(i,j,k,ispec) = sngl(etax)
      etaystock(i,j,k,ispec) = sngl(etay)
      etazstock(i,j,k,ispec) = sngl(etaz)
      gammaxstock(i,j,k,ispec) = sngl(gammax)
      gammaystock(i,j,k,ispec) = sngl(gammay)
      gammazstock(i,j,k,ispec) = sngl(gammaz)
      xmeshstock(i,j,k,ispec) = sngl(xmesh - xmin)
      ymeshstock(i,j,k,ispec) = sngl(ymesh - ymin)
      zmeshstock(i,j,k,ispec) = sngl(zmesh - zmin)

      enddo
      enddo
      enddo

      return
      end


c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


c ---
c --- surface du milieu
c ---

      double precision function zsurf(x,y)
      implicit double precision (a-h,o-z)
      double precision x,y

      angle = 25.

c pente suivant X
      zsurf = 2000. + x*dtan(angle*3.14159265/180.)

c pente suivant Y
      zsurf = 2000. + y*dtan(angle*3.14159265/180.)

c gaussienne asymetrique AGU
      sigmax = 250.
      sigmay = sigmax / 2.
      height = 180.
      x0 = 1040.
      y0 = 1040.
      zsurf = 1050. +
     .      height*dexp(-((x-x0)**2)/(2.*sigmax**2))*
     .      dexp(-((y-y0)**2)/(2.*sigmay**2))

c tests analytiques sans topo
      zsurf = 1050.

      return
      end


c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      subroutine rank(A,IND,N)
C
C     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
C
!      implicit double precision (a-h,o-z)
!      double precision A(1)

      implicit real (a-h,o-z)
      real A(1)

      integer IND(1)
C
      DO 10 J=1,N
         IND(j)=j
   10 continue
C
      if (n.eq.1) return
      L=n/2+1
      ir=n
  100 CONTINUE
         IF (l.gt.1) THEN
            l=l-1
            indx=ind(l)
            q=a(indx)
         ELSE
            indx=ind(ir)
            q=a(indx)
            ind(ir)=ind(1)
            ir=ir-1
            if (ir.eq.1) then
               ind(1)=indx
               return
            endif
         ENDIF
         i=l
         j=l+l
  200    CONTINUE
         IF (J.le.IR) THEN
            IF (J.lt.IR) THEN
               IF ( A(IND(j)).lt.A(IND(j+1)) ) j=j+1
            ENDIF
            IF (q.lt.A(IND(j))) THEN
               IND(I)=IND(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GOTO 200
         ENDIF
         IND(I)=INDX
      GOTO 100
      END
c-----------------------------------------------------------------------
      subroutine swap(a,w,ind,n)
C
C     Use IND to sort array A   (p 233 Num. Rec.), 5/26/93 pff.
C
!      implicit double precision (a-h,o-z)
!      double precision A(1),W(1)

      implicit real (a-h,o-z)
      real A(1),W(1)

      integer IND(1)
C
      DO 10 J=1,N
         W(j)=A(j)
   10 continue
C
      DO 20 J=1,N
         A(j)=W(ind(j))
   20 continue
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine iswap(a,w,ind,n)
C
C     Use IND to sort array A
C
!      implicit double precision (a-h,o-z)

      implicit real (a-h,o-z)

      INTEGER A(1),W(1),IND(1)
C
      DO 10 J=1,N
         W(j)=A(j)
   10 continue
C
      DO 20 J=1,N
         A(j)=W(ind(j))
   20 continue
      RETURN
      END

