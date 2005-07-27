!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine plotpost_norm(coord2D,normdispl2D,it,dt,nspec2D,NEX_XI,NER, &
               x_target_no_rot,z_target_no_rot,nrec,x_target_source_no_rot,z_target_source_no_rot)

  implicit none

  include "constants.h"

!!!!!!!! DK DK XXXXXXXXXXX YYYYYYYYYY UGLY for Carcione copper aniso
  include "carcione_anisotropy.h"

! threshold to cut display to save on size of Postscript file
  double precision, parameter :: cutvect = 1.d0 / 100.d0

! non-linear power law to enhance small amplitudes
  double precision, parameter :: POWER_NON_LINEAR = 0.35d0

  integer it,nspec2D,NEX_XI,NER
  double precision dt,timeval
  double precision coord2D(2,NGLLX,NGLLZ,nspec2D)
  double precision normdispl2D(NGLLX,NGLLZ,nspec2D)

  double precision xmin,zmin,xmax,zmax,height,usoffset
  double precision grayscale_val,grayscale_val1,grayscale_val2,grayscale_val3,grayscale_val4,grayscale_val_non_linear
  double precision x_target_source_no_rot,z_target_source_no_rot

  integer nrec,irec
  double precision, dimension(nrec) :: x_target_no_rot,z_target_no_rot

  character(len=100) name

  double precision convert
  double precision x1,z1,x2,z2,x3,z3,x4,z4
  double precision xoffset,zoffset

  integer i,k,ispec2D,ix,iz

  logical, parameter :: drawmesh = .true.
  logical, parameter :: usletter = .false.
  double precision, parameter :: scalex   = 1.d0
  double precision, parameter :: scalez   = 1.d0
  double precision, parameter :: sizemax  = 1.d0
  double precision, parameter :: orig_x   = 2.4d0
  double precision, parameter :: orig_z   = 2.9d0
  double precision, parameter :: centim = 28.5d0
  double precision, parameter :: angle = 20.d0
  double precision, parameter :: rapport = 0.40d0

! taille de la fenetre de display Postscript en pourcentage de la feuille
  double precision, parameter :: rpercentx = 70.0d0, rpercentz = 77.0d0

  double precision sizex,sizez,rapp_page,dispmax

! papier A4 ou US letter
  if(usletter) then
    usoffset = 1.75d0
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.
    sizex = 29.7d0
    sizez = 21.d0
  endif

! recherche des positions maximales des points de la grille
  xmin = minval(coord2D(1,:,:,:))
  zmin = minval(coord2D(2,:,:,:))
  xmax = maxval(coord2D(1,:,:,:))
  zmax = maxval(coord2D(2,:,:,:))

! calcul de la valeur maximum de la norme du deplacement projete
  dispmax = maxval(normdispl2D)

! rapport taille page/taille domaine physique
  rapp_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) / 100.d0

! hauteur des numeros de domaine en CM
  height = 0.25d0

!
!---- ouverture du fichier PostScript
!
  write(name,222) it
  open(unit=24,file=name,status='unknown')
  222 format('OUTPUT_FILES/vect',i5.5,'.ps')

!
!---- ecriture de l'entete du fichier PostScript
!
  write(24,10)
  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/MK {mark} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/C1 {lineto closepath stroke} def'
  write(24,*) '/C2 {lineto closepath fill stroke} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '/GG {0 setgray ST} def'
  write(24,*) '/GC {Colmesh ST} def'
  write(24,*) '/RF {setrgbcolor fill} def'
  write(24,*) '/SF {setgray fill} def'
  write(24,*) '/SG {setgray} def'
  write(24,*) '/GS {gsave} def'
  write(24,*) '/GR {grestore} def'
  write(24,*) '/SLW {setlinewidth} def'
  write(24,*) '/SCSF {scalefont setfont} def'
  write(24,*) '% differents symboles utiles'
  write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/Cross {GS 1 0 1 setrgbcolor 0.05 CM SLW'
  write(24,*) 'GS 3 3 MR -6. -6. LR ST GR'
  write(24,*) 'GS 3 -3 MR -6. 6. LR ST GR'
  write(24,*) '0.01 CM SLW} def'
  write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(24,*) '/Losange {GS 0 0 1 setrgbcolor 0.03 CM SLW 0 4.2 2 div MR'
  write(24,*) '-3 2 div -4.2 2 div LR 3 2 div -4.2 2 div LR 3 2 div 4.2 2 div LR CP ST'
  write(24,*) 'GR 0.01 CM SLW} def'
  write(24,*) '%'
  write(24,*) '% niveaux de gris pour le modele de vitesse'
  write(24,*) '/BK {setgray fill} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/BK {pop 1 setgray fill} def'
  write(24,*) '%'
  write(24,*) '% magenta pour les vecteurs deplacement'
  write(24,*) '/Colvects {0.01 CM SLW 1. 0. 1. RG} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/Colvects {0.01 CM SLW 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% chartreuse pour le maillage des macroblocs'
  write(24,*) '/Colmesh {0.02 CM SLW 0.5 1. 0. RG} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/Colmesh {0.02 CM SLW 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% cyan pour les sources et recepteurs'
  write(24,*) '/Colreceiv {0. 1. 1. RG} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/Colreceiv {0. setgray} def'
  write(24,*) '%'
  write(24,*) '% macro dessin fleche'
  write(24,*) '/F {MV LR gsave LR ST grestore LR ST} def'
  write(24,*) '% macro dessin contour elements'
  write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(24,*) '%'
  write(24,*) '.01 CM SLW'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.35 CM SCSF'
  write(24,*) '%'
  write(24,*) '/vshift ',-height/2,' CM def'
  write(24,*) '/Rshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop neg vshift MR show } def'
  write(24,*) '/Cshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(24,*) '/fN {/Helvetica-Bold findfont ',sngl(height),' CM SCSF} def'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(24,*) '%'

!
!--- ecriture des legendes du fichier PostScript
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM SCSF'

  write(24,*) '24. CM 1.2 CM MV'
  write(24,610) usoffset,it
  write(24,*) '%'

  write(24,*) '24. CM 1.95 CM MV'
  timeval = it*dt
  if(timeval  >=  1.d-3) then
  write(24,600) usoffset,timeval
  else
  write(24,601) usoffset,timeval
  endif
  write(24,*) '%'
  write(24,*) '24. CM 2.7 CM MV'
  write(24,640) usoffset,dispmax
  write(24,*) '%'
  write(24,*) '24. CM 3.45 CM MV'
!! DK DK UGLY  write(24,620) usoffset,cutvect*100.d0
  write(24,620) usoffset,ANGLE_MESH

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM SCSF'
  write(24,*) '11 CM 1.1 CM MV'
  write(24,*) '(X axis) show'
  write(24,*) '%'
  write(24,*) '1.4 CM 9.5 CM MV'
  write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(24,*) '(Z axis) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.7 CM SCSF'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) sngl(usoffset),' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(Norm of full 3D displacement vector) show'
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) sngl(usoffset),' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(Anisotropic Copper Crystal) show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) sngl(usoffset),' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(Elastic Wave - Spectral Element Method) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) scalex,' ',scalez,' scale'
  write(24,*) '%'

  convert = pi/180.d0

!
!----  draw the norm of the full 3D displacement vector
!

! return if the maximum displacement equals zero (no source)
  if (dispmax == 0.d0) then
    print *,' null displacement : returning !'
    return
  endif

  write(24,*) '0 setgray'

  write(24,*) '%'
  write(24,*) '% norm of displacement vector with gray scale'
  write(24,*) '%'

  do ispec2D=1,nspec2D
  do i=1,NGLLX-1
  do k=1,NGLLZ-1

  x1 =(coord2D(1,i,k,ispec2D)-xmin)*rapp_page
  z1 =(coord2D(2,i,k,ispec2D)-zmin)*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim

  x2 =(coord2D(1,i+1,k,ispec2D)-xmin)*rapp_page
  z2 =(coord2D(2,i+1,k,ispec2D)-zmin)*rapp_page
  x2 = (orig_x+x2) * centim
  z2 = (orig_z+z2) * centim

  x3 =(coord2D(1,i+1,k+1,ispec2D)-xmin)*rapp_page
  z3 =(coord2D(2,i+1,k+1,ispec2D)-zmin)*rapp_page
  x3 = (orig_x+x3) * centim
  z3 = (orig_z+z3) * centim

  x4 =(coord2D(1,i,k+1,ispec2D)-xmin)*rapp_page
  z4 =(coord2D(2,i,k+1,ispec2D)-zmin)*rapp_page
  x4 = (orig_x+x4) * centim
  z4 = (orig_z+z4) * centim

!! DK DK draw grayscale rectangles, white if no signal, black if signal
!! DK DK use average value based on four corners, and non-linear enhancing
  grayscale_val1 = normdispl2D(i,k,ispec2D)/dispmax
  grayscale_val2 = normdispl2D(i+1,k,ispec2D)/dispmax
  grayscale_val3 = normdispl2D(i+1,k+1,ispec2D)/dispmax
  grayscale_val4 = normdispl2D(i,k+1,ispec2D)/dispmax
  grayscale_val = (grayscale_val1 + grayscale_val2 + grayscale_val3 + grayscale_val4) / 4.d0
  if(grayscale_val < 0.d0) grayscale_val = 0.d0
  if(grayscale_val > 1.d0) grayscale_val = 1.d0
  grayscale_val_non_linear = grayscale_val**POWER_NON_LINEAR

!! DK DK only draw if value is above threshold, otherwise ignore
  if(grayscale_val > cutvect) then
    write(24,300) 1.d0 - grayscale_val_non_linear
    write(24,301) x1,z1,x2,z2,x3,z3,x4,z4
  endif

  enddo
  enddo
  enddo


!--- draw edges of the model in red

  write(24,*) '% drawing SEM grid with red lines'

  write(24,*) '1 0 0 setrgbcolor'

  x1 =(xmin-xmin)*rapp_page
  z1 =(zmin-zmin)*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,*) x1,z1,' M'

  x1 =(xmax-xmin)*rapp_page
  z1 =(zmin-zmin)*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,*) x1,z1,' L'

  x1 =(xmax-xmin)*rapp_page
  z1 =(zmax-zmin)*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,*) x1,z1,' L'

  x1 =(xmin-xmin)*rapp_page
  z1 =(zmax-zmin)*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,*) x1,z1,' C1'

! draw spectral element mesh

  if(drawmesh) then

! first draw vertical lines
  do ix = 0,NEX_XI

  xoffset = dble(ix) * (xmax - xmin) / dble(NEX_XI)

  x1 = xoffset*rapp_page
  z1 = 0.
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,100) x1,z1

  x1 = xoffset*rapp_page
  z1 = (zmax-zmin)*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,101) x1,z1

  enddo

! then draw horizontal lines
  do iz = 0,NER

  zoffset = dble(iz) * (zmax - zmin) / dble(NER)

  x1 = 0.
  z1 = zoffset*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,100) x1,z1

  x1 = (xmax-xmin)*rapp_page
  z1 = zoffset*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,101) x1,z1

  enddo

!! DK DK display the receivers with losanges
!! DK DK rotation is not needed because we saved the original target location
  write(24,*) '% drawing receivers with diamonds'

  do irec = 1,nrec
    x1 =(x_target_no_rot(irec)-xmin)*rapp_page
    z1 =(z_target_no_rot(irec)-zmin)*rapp_page
    x1 = (orig_x+x1) * centim
    z1 = (orig_z+z1) * centim
    write(24,*) sngl(x1),' ',sngl(z1),' M Losange'
  enddo

!! DK DK display the source with a cross
!! DK DK rotation is not needed because we saved the original target location
  write(24,*) '% drawing source with a cross'

  x1 =(x_target_source_no_rot-xmin)*rapp_page
  z1 =(z_target_source_no_rot-zmin)*rapp_page
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  write(24,*) sngl(x1),' ',sngl(z1),' M Cross'

  endif

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  close(24)

 10   format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: snaphot',/, &
          '%% Created by: Specfem3D',/, &
          '%% Author: Dimitri Komatitsch',/,'%%')

 100  format(f8.3,1x,f8.3,' M')
 101  format(f8.3,1x,f8.3,' L ST')

 200  format(80(a1))

 300  format(f8.4,' SG')
 301  format(f8.3,1x,f8.3,' M ',f8.3,1x,f8.3,' L ',f8.3,1x,f8.3,' L ',f8.3,1x,f8.3,' C2')

 600  format(f6.3,' neg CM 0 MR (Time =',f6.3,' s) show')
 601  format(f6.3,' neg CM 0 MR (Time =',1pe10.3,' s) show')
 610  format(f6.3,' neg CM 0 MR (Time step = ',i5,') show')
!! DK DK UGLY 620  format(f6.3,' neg CM 0 MR (Cut =',f5.3,' percent) show')
 620  format(f6.3,' neg CM 0 MR (Angle =',f11.3,' deg) show')
 640  format(f6.3,' neg CM 0 MR (Max norm =',1pe10.3,') show')

  end subroutine plotpost_norm

