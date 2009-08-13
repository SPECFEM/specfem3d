
      subroutine plotpost(Uxnew,Uynew,npoin,indicedisp,nxdisp,nzdisp,
     .            deltaxdisp,deltazdisp,it,dtinc,indice)

c
c routine affichage postscript
c

      implicit none

      integer it,indice,npoin,nxdisp,nzdisp
      real dtinc
      real deltaxdisp,deltazdisp

      integer indicedisp(nxdisp,nzdisp)

      real Uxnew(npoin)
      real Uynew(npoin)

      integer ix,iz,longueur,ii
      real cutvect,sizemax,angle,rapport,sizex,sizez,orig_x,orig_z
      real centim,pi,convert,xmax,zmax,xmin,zmin,rapp_page,dispmax
      real d,height,x1,z1,x2,z2,d1,d2,dummy,theta,thetaup,thetadown
      real xa,za,xb,zb

      character*50 stitle
      character*40 namefile

      character*100 name
      character ch1(100),ch2(100)
      equivalence (name,ch1)
      logical first

      if(indice.eq.1) then
            stitle = 'Simu 3D plan X = constante'
      else if(indice.eq.2) then
            stitle = 'Simu 3D plan Y = constante'
      else if(indice.eq.3) then
            stitle = 'Simu 3D plan Z = constante'
      else
            stop 'Bad value indice in plotpost'
      endif

      cutvect = 0.01
      sizemax = 1.
      angle = 20.
      rapport = 0.40

! format europeen A4
      sizex = 21.
      sizez = 29.7

! format US letter
      sizex = 8.5 * 2.54
      sizez = 11. * 2.54

      orig_x = 2.20
      orig_z = 1.40
      centim = 28.5
      pi = 3.14159265
      convert = pi/180.

      print *,stitle

c recherche des positions maximales des points de la grille
      xmax=deltaxdisp*real(nxdisp-1)
      zmax=deltazdisp*real(nzdisp-1)

      write(*,*) 'Max X = ',xmax
      write(*,*) 'Max Z = ',zmax

      write(*,*) 'nxdisp = ',nxdisp
      write(*,*) 'nzdisp = ',nzdisp

c limite du repere physique
      xmin=0.0
      zmin=0.0

c rapport taille page/taille domaine physique
c on prend 80 % de la page en landscape (attention axes permutes)
      rapp_page = 0.80*amin1(sizex/(zmax-zmin),sizez/(xmax-xmin))

c recherche de la valeur maximum de la norme du deplacement
      dispmax = -10.0
      do ix=1,nxdisp
      do iz=1,nzdisp
      if(indicedisp(ix,iz) .gt. 0) then
        d = sqrt(Uxnew(indicedisp(ix,iz))**2 + Uynew(indicedisp(ix,iz))**2)
        dispmax = amax1(d,dispmax)
      endif
      enddo
      enddo
      write(*,*) 'Max norme = ',dispmax

c hauteur des numeros de domaine en CM
      height = 0.25

c ouverture du fichier PostScript
      if(indice.eq.1) then
            write(namefile,221) it
      else if(indice.eq.2) then
            write(namefile,222) it
      else if(indice.eq.3) then
            write(namefile,223) it
      endif

 221  format('vectx',i5.5,'.ps')
 222  format('vecty',i5.5,'.ps')
 223  format('vectz',i5.5,'.ps')
      open(unit=24,file=namefile,status='unknown')

c ecriture de l'entete du fichier PostScript
      write(24,10) stitle
      write(24,*) '/CM {',centim,' mul} def'
      write(24,*) '/LR {rlineto} def'
      write(24,*) '/LT {lineto} def'
      write(24,*) '/L {lineto} def'
      write(24,*) '/MR {rmoveto} def'
      write(24,*) '/MV {moveto} def'
      write(24,*) '/M {moveto} def'
      write(24,*) '/MK {mark} def'
      write(24,*) '/ST {stroke} def'
      write(24,*) '/CP {closepath} def'
      write(24,*) '/RG {setrgbcolor} def'
      write(24,*) '/GF {gsave fill grestore} def'
      write(24,*) '/GG {0 setgray ST} def'
      write(24,*) '/GC {Colmesh ST} def'
      write(24,*) '/RF {setrgbcolor fill} def'
      write(24,*) '/SF {setgray fill} def'
      write(24,*) '/GS {gsave} def'
      write(24,*) '/GR {grestore} def'
      write(24,*) '% differents symboles utiles'
      write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
      write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
      write(24,*) 'CP fill} def'
      write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
      write(24,*) 'CP fill} def'
      write(24,*) '/Cross {GS 0.05 CM setlinewidth'
      write(24,*) 'GS 3 3 MR -6. -6. LR ST GR'
      write(24,*) 'GS 3 -3 MR -6. 6. LR ST GR'
      write(24,*) '0.01 CM setlinewidth} def'
      write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
      write(24,*) '/Losange {GS 0.05 CM setlinewidth 0 4.2 MR'
      write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
      write(24,*) 'GR 0.01 CM setlinewidth} def'
      write(24,*) '%'
      write(24,*) '% decalage des dessins de div et curl en cm'
      write(24,*) '/U { 0 add } def'
      write(24,*) '%'
      write(24,*) '% niveaux de gris pour le modele de vitesse'
      write(24,*) '/BK {setgray fill} def'
      write(24,*) '% version noir et blanc'
      write(24,*) '%/BK {pop 1 setgray fill} def'
      write(24,*) '%'
      write(24,*) '% magenta pour les vecteurs deplacement'
      write(24,*) '/Colvects {0.01 CM setlinewidth 1. 0. 1. RG} def'
      write(24,*) '% version noir et blanc'
      write(24,*) '%/Colvects {0.01 CM setlinewidth 0. setgray} def'
      write(24,*) '%'
      write(24,*) '% chartreuse pour le maillage des macroblocs'
      write(24,*) '/Colmesh {0.02 CM setlinewidth 0.5 1. 0. RG} def'
      write(24,*) '% version noir et blanc'
      write(24,*) '%/Colmesh {0.02 CM setlinewidth 0. setgray} def'
      write(24,*) '%'
      write(24,*) '% cyan pour les sources et recepteurs'
      write(24,*) '/Colreceiv {0. 1. 1. RG} def'
      write(24,*) '% version noir et blanc'
      write(24,*) '%/Colreceiv {0. setgray} def'
      write(24,*) '%'
      write(24,*) '% macro dessin fleche'
      write(24,*) '/F {MV LR gsave LR ST grestore LR ST} def'
      write(24,*) '% macro dessin contour elements'
      write(24,*)
     .     '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
      write(24,*) '%'
      write(24,*) '.01 CM setlinewidth'
      write(24,*) '/Times-Roman findfont'
      write(24,*) '.35 CM scalefont setfont'
      write(24,*) '%'
      write(24,*) '/vshift ',-height/2,' CM def'
      write(24,*) '/Rshow { currentpoint stroke MV'
      write(24,*) 'dup stringwidth pop neg vshift MR show } def'
      write(24,*) '/Cshow { currentpoint stroke MV'
      write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
      write(24,*) '/fN {/Helvetica-Bold findfont ',height,
     .      ' CM scalefont setfont } def'
      write(24,*) '%'
      write(24,*) 'gsave newpath 90 rotate'
      write(24,*) '0 -21. CM translate 1. 1. scale'
      write(24,*) '%'

c
c--- ecriture des legendes du fichier PostScript
c
      write(24,*) '0 setgray'
      write(24,*) '/Times-Roman findfont'
      write(24,*) '.5 CM scalefont setfont'

      write(name,610) it
      write(24,*) '22.6 CM 1.2 CM MV'
      write(24,*) name
      write(24,*) '%'
      write(name,600) it*dtinc
      write(24,*) '22.6 CM 1.95 CM MV'
      write(24,*) name
      write(24,*) '%'
      write(name,640) dispmax
      write(24,*) '22.6 CM 2.7 CM MV'
      write(24,*) name
      write(24,*) '%'
      write(name,620) cutvect*100.0
      write(24,*) '22.6 CM 3.45 CM MV'
      write(24,*) name

      write(24,*) '%'
      write(24,*) '/Times-Roman findfont'
      write(24,*) '.6 CM scalefont setfont'
      write(24,*) '11 CM 0.68 CM MV'
      write(24,*) '(X axis) show'
      write(24,*) '%'
      write(24,*) '1.4 CM 9.5 CM MV'
      write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
      write(24,*) '(Y axis) show'
      write(24,*) 'grestore'
      write(24,*) '%'
      write(24,*) '/Times-Roman findfont'
      write(24,*) '.7 CM scalefont setfont'
      write(24,*) '24.35 CM 18.9 CM MV'
      write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
      write(24,*) '(',stitle,') show'
      write(24,*) 'grestore'
      write(24,*) '25.45 CM 18.9 CM MV'
      write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
      write(24,*) '(Elastic Wave 3D - Spectral Element Method) show'
      write(24,*) 'grestore'
      write(24,*) '%'
      write(24,*) '1. 1. scale'
      write(24,*) '%'

c trace les contours de la "boite" interpolee

      x1 =((xmin-xmin)*rapp_page + orig_x)*centim
      z1 =((zmin-zmin)*rapp_page + orig_z)*centim
      write(24,*) x1,z1,' MV'
      x2 =((xmax-xmin)*rapp_page + orig_x)*centim
      z2 =((zmin-zmin)*rapp_page + orig_z)*centim
      write(24,*) x2,z2,' LT'
      x2 =((xmax-xmin)*rapp_page + orig_x)*centim
      z2 =((zmax-zmin)*rapp_page + orig_z)*centim
      write(24,*) x2,z2,' LT'
      x2 =((xmin-xmin)*rapp_page + orig_x)*centim
      z2 =((zmax-zmin)*rapp_page + orig_z)*centim
      write(24,*) x2,z2,' LT'
      write(24,*) 'closepath ST'

c ***************************************

c
c----  draw the normalized displacement field
c

c return if the maximum displacement equals zero (no source)
      if (dispmax.eq.0.) then
            print *,' null displacement : returning !'
      else

            write(24,*) '0 setgray'

c tracer les vecteurs deplacement aux noeuds du maillage

        do ix=1,nxdisp
        do iz=1,nzdisp

      if(indicedisp(ix,iz) .gt. 0) then

        x1 =(deltaxdisp*real(ix-1)-xmin)*rapp_page
        z1 =(deltazdisp*real(iz-1)-zmin)*rapp_page

        x2 = Uxnew(indicedisp(ix,iz))*sizemax/dispmax
        z2 = Uynew(indicedisp(ix,iz))*sizemax/dispmax

        d = sqrt(x2**2 + z2**2)

c ignorer si vecteur trop petit
        if (d.gt.cutvect*sizemax) then

        d1 = d * rapport
        d2 = d1 * cos(angle*convert)

        dummy = x2/d
        if (dummy.gt.0.9999) dummy = 0.9999
        if (dummy.lt.-0.9999) dummy = -0.9999
        theta = acos(dummy)

        if(z2.lt.0.0) then
        theta = 360.0*convert - theta
        endif
        thetaup = theta - angle*convert
        thetadown = theta + angle*convert

c tracer le vecteur proprement dit
        x1 = (orig_x+x1) * centim
        z1 = (orig_z+z1) * centim
        x2 = x2 * centim
        z2 = z2 * centim
        xa = -d2*cos(thetaup)
        za = -d2*sin(thetaup)
        xa = xa * centim
        za = za * centim
        xb = -d2*cos(thetadown)
        zb = -d2*sin(thetadown)
        xb = xb * centim
        zb = zb * centim
        write(name,700) xb,zb,xa,za,x2,z2,x1,z1

c filtrer les blancs inutiles pour diminuer taille fichier PostScript
        longueur = 49
        indice = 1
        first = .false.
        do ii=1,longueur-1
            if(ch1(ii).ne.' '.or.first) then
            if(ch1(ii).ne.' '.or.ch1(ii+1).ne.' ') then
                  ch2(indice) = ch1(ii)
                  indice = indice + 1
                  first = .true.
            endif
            endif
        enddo
        ch2(indice) = ch1(longueur)
        write(24,200) (ch2(ii),ii=1,indice)

      endif

      endif

      enddo
      enddo

      endif

c ***************************************

      write(24,*) '%'
      write(24,*) 'grestore'
      write(24,*) 'showpage'

      close(24)

 10   format('%!PS-Adobe-2.0',/,'%%',/,
     .      '%% Title: ',a50,/,
     .      '%% Creator: Dimitri Komatitsch - Specfem 3D',/,'%%')
 510  format(f5.1,1x,f5.1,' M')
 600  format('(Time =',f6.3,' s) show')
 610  format('(Time step =',i4,') show')
 620  format('(Cut =',f5.2,' \%) show')
 640  format('(Max norm =',e10.3,') show')
 910  format(f5.1,1x,f5.1,' U M')

 200  format(80(a1))
 499  format(f5.1,1x,f5.1,' L')
 501  format('fN (G',i2,' M',i2,') Cshow')
 502  format('fN (',i4,') Cshow')
 601  format(f6.2,1x,f6.2)
 602  format('CP ',f4.2,1x,f4.2,1x,f4.2,' RF')
 604  format('CP ',f4.2,' BK')
 700  format(8(f5.1,1x),'F')

      end
