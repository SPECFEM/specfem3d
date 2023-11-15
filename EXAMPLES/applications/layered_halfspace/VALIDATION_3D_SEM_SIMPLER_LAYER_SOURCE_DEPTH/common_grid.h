
      implicit none

      integer nx,ny
      integer nxgll

! definition taille du maillage
      parameter (nxgll=6)

      parameter (nx=112, ny=112)

! nb points de grille pour la simulation
! nz1 couche du bas, nz2 couche du milieu, nz3 couche du haut
      integer nz,nz1,nz2,nz3
      parameter(nz1 = 48, nz2 = 1, nz3 = 1)
      parameter(nz = nz1 + nz2 + nz3)

      integer ndime,ndofn,nygll,nzgll,nxyz,ntot
      integer nspec
      integer nspechaut,nspecbas,nspecderaf1,nspecderaf2,nspecsansderaf
      integer npoin
      integer idegpoly

      parameter(ndime=3)
      parameter(ndofn=ndime)

      parameter(nygll=nxgll)
      parameter(nzgll=nxgll)
      parameter(idegpoly=nxgll-1)

! taille physique grille pour la simulation
! obtenue par projection de la zone physique par Lambert
!  (139.15E <-> 140.85E) et (34.9N <-> 36.55N)
      double precision xmin,xmax,ymin,ymax
      parameter(xmin = 0.)
      parameter(xmax = 134000.)
      parameter(ymin = 0.)
      parameter(ymax = 134000.)

! position du fond de la grille (en metres, signe moins pour la profondeur)
      double precision zmin
      parameter(zmin = - 60000.)

! position physique de la source
! actuellement au milieu dans le plan (X,Y), aux deux tiers de la profondeur
      double precision xpositsource,ypositsource,zpositsource

      parameter(xpositsource = xmin + (xmax - xmin) / 2.0)
      parameter(ypositsource = ymin + (ymax - ymin) / 2.0)
      parameter(zpositsource = - 25000.)

! frequence de reference (frequence centrale de la source) et ratio f_max / f0
! pour estimation du nombre de points par longueur d'onde dans le maillage
      double precision f0,t0,factor_force,freqcste
      double precision sig11,sig12,sig13
      double precision sig21,sig22,sig23
      double precision sig31,sig32,sig33
      double precision M0

      parameter(f0 = 0.4)
      parameter(t0 = 2.6)

      parameter(factor_force = 1.e15)

! valeur de M0 en N.m
      parameter(M0 = 3.8e15)

! Tokyo Tsuboi
      parameter(sig11 = -0.1796 *M0 ,sig12 =  -0.891*M0 , sig13 = -0.2949 *M0 )
      parameter(sig22 = -0.03105 *M0 ,sig23 =  0.2834 *M0 )
      parameter(sig33 =  0.2106 *M0 )

! imposer les autres composantes du moment par symetrie
      parameter(sig21 = sig12)
      parameter(sig31 = sig13)
      parameter(sig32 = sig23)

! DK DK source explosive composantes du moment

      parameter(freqcste = 2.5)

! bords absorbants actifs ou non
      logical bords_abs
      parameter(bords_abs = .true.)

! ecart minimum entre les interfaces "tapered" pour generation de la grille
      double precision delta_min_topo_18,delta_min_18_28
      parameter(delta_min_topo_18=650.,delta_min_18_28=900.)

! nb d'elements spectraux dans les differentes zones
      parameter(nspechaut=nx*ny*(nz2+nz3))
      parameter(nspecbas=(nx/4)*(ny/4)*(nz1-8)/4)
      parameter(nspecderaf1=9*(nx/2)*(ny/2))
      parameter(nspecderaf2=9*(nx/4)*(ny/4))

! nombre total d'elements spectraux dans la simu
      parameter(nspec=nspechaut + nspecbas + nspecderaf1 + nspecderaf2)

! nombre total d'elements spectraux si on ne deraffinait pas
! avec facteur 2 dans le bedrock
      parameter(nspecsansderaf=nx*ny*(nz2+nz3+(nz1/2)))
! avec facteur 4 dans le bedrock
!        parameter(nspecsansderaf=nx*ny*(nz2+nz3+(nz1/4)))

      integer isizebloc,isizefacehoriz,isizefacevert,isizearetes,isizetheor
      integer isizesedim,isizebedrock

! valeur theorique du nombre de points
      parameter(isizebloc = 180 * nxgll**3 - 434 * nxgll**2 + 346 * nxgll - 91)
      parameter(isizefacehoriz = 30 * nxgll**2 - 47*(nxgll-2)
     .      - 22*(2-1) - 6*(3-1) -  9*(4-1) - 3*(6-1))
      parameter(isizefacevert  = 42 * nxgll**2 - 71*(nxgll-2)
     .      - 21*(2-1) - 6*(3-1) - 23*(4-1) - 2*(6-1))
      parameter(isizearetes    = 8*(nxgll-1) + 1)
      parameter(isizesedim     =
     .   (nx*(nxgll-1)+1)*(ny*(nygll-1)+1)*((nz2+nz3)*(nzgll-1)))
      parameter(isizebedrock  =
     .   ((nx/4)*(nxgll-1)+1)*((ny/4)*(nygll-1)+1)*(((nz1-8)/4)*(nzgll-1)))
      parameter(isizetheor = (nx/8)*(ny/8)*isizebloc
     .      - (nx/8)*(ny/8-1)*isizefacehoriz
     .      - (nx/8-1)*(ny/8)*isizefacevert
     .      + (nx/8-1)*(ny/8-1)*isizearetes
     .      + isizesedim + isizebedrock)

      parameter(npoin = isizetheor)

c pour routine de numerotation globale
      parameter(nxyz = nxgll*nygll*nzgll)
      parameter(ntot = nxyz*nspec)

! taille grille pour definition des interfaces en entree
      integer nxinterf,nyinterf
      parameter(nxinterf = 366, nyinterf = 445)

! dessiner interfs et grille 3D AVS ou non
      logical draw_interfs_avs,draw_grid_avs
      parameter(draw_interfs_avs = .true.)
      parameter(draw_grid_avs = .false.)

! sauver fichiers grille GNUPLOT ou non
      logical ignuplot
      parameter(ignuplot = .false.)

! save huge binary files for SPECFEM or not
      logical save_binary_files
      parameter(save_binary_files = .true.)

! read the real shape of the layers, or use a mean value
      logical use_averaged_layers
      parameter(use_averaged_layers = .true.)

! execute the grid generation process or just estimate the grid parameters
      logical iexec
      parameter(iexec = .true.)

