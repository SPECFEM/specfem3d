      program setbrn
c
c Herewith the new version of setbran.f with separated table and header
c and correct index assignments for direct access (hopefully)
c
      include 'limits.inc'
      save
      character*8 code,phcd
      character*20 modnam
      character*20 modnam2
      double precision zm,pm,pb,pu,taup,xp,taul,px,xt,xl,pux,pt,taut,
     1 coef,xa
      double precision tmp(nsl1,2),xm(nsl1,2),deg,dtol,zmax,zoc,zic,
     1 z0
      dimension ndx2(nsr0,2)
      common/umodc/zm(nsr0,2),pm(nsr0,2),ndex(nsr0,2),mt(2)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/brnc/phcd(jbrn)
      common/msqc/pux(jbrn,2),km(2),midx(jbrn,2)
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data nin,nout,xmin,dtol/1,2,200.,1d-6/
      deg=180d0/3.1415927d0
c
c     write(6,*) "rec length for dasign:"
c     read(5,*) ndasr
c
      call assign(nin,-1,'remodl.hed')
      read(nin)ndasr,modnam,zmax,zoc,zic,kb,(pb(i),i=1,kb(2)),
     1 mt,lt,lbb,lcb,xn,pn,tn
      read(nin)((lbrk(i,nph),i=1,lbb(nph)),(code(i,nph),
     1 i=1,lcb(nph)),(zm(i,nph),pm(i,nph),ndex(i,nph),i=1,mt(nph)),
     2 (lvz(i,nph),taul(i,nph),xl(i,nph),i=1,lt(nph)),nph=1,2)
      call retrns(nin)
      print *,'ndasr =',ndasr,'  modnam = ',modnam
c
      call dasign(nin,-1,'remodl.tbl',ndasr)
      nrec=0
      do 1 nph=1,2
      n1=kb(nph)
      ind=0
      do 2 k=1,n1
 2    xm(k,nph)=0d0
 3    nrec=nrec+1
      read(nin,rec=nrec)z0,n,(tmp(k,1),k=1,n),(tmp(k,2),k=1,n)
      if(ind.gt.0) go to 4
      if(dabs(z0-zoc).le.dtol) go to 4
      j=1
      do 5 k=2,n
      xm(k,nph)=dmax1(xm(k,nph),dabs(tmp(j,2)-tmp(k,2)))
 5    j=k
      if(n+1.eq.n1) xm(n1,nph)=tmp(n,2)
      go to 3
 4    ind=ind+1
      do 6 k=1,n
      taup(k,ind,nph)=tmp(k,1)
 6    xp(k,ind,nph)=tmp(k,2)
      if(ind.lt.3) go to 3
 1    continue
c
      xmin=xn*xmin
c
      call assign(10,2,'setbrn1.lis')
      write(10,*)'kb mt lt lbb lcb',kb,mt,lt,lbb,lcb
      write(10,*)'xn pn tn xmin',xn,pn,tn,xmin
      cn=1./xn
      write(10,209)(i,(lbrk(i,j),code(i,j),j=1,2),i=1,lbb(1))
 209  format(/(1x,2i5,2x,a,i5,2x,a))
      write(10,210)(i,lbrk(i,2),code(i,2),i=lbb(1)+1,lbb(2))
 210  format(1x,i5,15x,i5,2x,a)
      write(10,200,iostat=ios)(i,(zm(i,j),pm(i,j),ndex(i,j),j=1,2),
     1 i=1,mt(1))
 200  format(/(1x,i5,2f12.6,i5,2x,2f12.6,i5))
      write(10,201,iostat=ios)(i,zm(i,2),pm(i,2),ndex(i,2),
     1 i=mt(1)+1,mt(2))
 201  format(1x,i5,31x,2f12.6,i5)
      write(10,217)((nph,i,lvz(i,nph),taul(i,nph),deg*xl(i,nph),
     1 i=1,lt(nph)),nph=1,2)
 217  format(/(1x,3i5,f12.6,f12.2))
      write(10,202)(i,pb(i),cn*xm(i,1),cn*xm(i,2),i=1,kb(1))
 202  format(/(5x,i5,f12.6,2f12.2))
      write(10,203)(i,pb(i),cn*xm(i,2),i=kb(1)+1,kb(2))
 203  format(5x,i5,f12.6,12x,f12.2)
      call retrns(10)
      call assign(10,2,'setbrn2.lis')
c
      do 8 nph=1,2
      n1=kb(nph)
      do 9 i=2,n1
      xm(i,nph)=xm(i-1,nph)+xm(i,nph)
      pu(i,nph)=pb(i)
 9    kuse(i,nph)=-1
      do 8 ind=3,2,-1
      jnd=ind-1
      do 8 i=1,n1
      taup(i,ind,nph)=taup(i,ind,nph)-taup(i,jnd,nph)
 8    xp(i,ind,nph)=xp(i,ind,nph)-xp(i,jnd,nph)
      do 10 nph=1,2
 10   call pdecx(kb(nph),nph,xm,2.)
c
      write(10,*)'ku',ku
      write(10,204)(i,(pu(i,nph),cn*xm(i,nph),cn*(xm(i+1,nph)-
     1 xm(i,nph)),nph=1,2),i=1,ku(1))
 204  format(/(5x,i5,2(f12.6,2f12.2)))
      write(10,205)(i,pu(i,2),cn*xm(i,2),cn*(xm(i+1,2)-xm(i,2)),
     1 i=ku(1)+1,ku(2))
 205  format(5x,i5,36x,f12.6,2f12.2)
      do 207 nph=1,2
 207  write(10,206)(i,pb(i),(taup(i,j,nph),j=1,3),(deg*xp(i,j,nph),
     1 j=1,3),i=1,kb(nph))
 206  format(/(1x,i5,4f10.6,3f10.2))
c
      call layout
c     write(10,214)(i,pb(i),(kuse(i,j),j=1,2),i=1,kb(2))
c214  format(/(5x,i5,f12.6,2i5))
      do 11 nph=1,2
      n1=kb(nph)
      k=0
      do 12 i=1,n1
      if(kuse(i,nph).lt.0) go to 12
      k=k+1
      pu(k,nph)=pb(i)
 12   continue
 11   ku(nph)=k
      call kseq
      call mseq
      write(10,215)(i,(pu(i,j),j=1,2),i=1,ku(1))
 215  format(/(5x,i5,2f12.6))
      write(10,216)(i,pu(i,2),i=ku(1)+1,ku(2))
 216  format(5x,i5,12x,f12.6)
      write(10,208)(i,(nafl(i,j),j=1,3),(indx(i,j),j=1,2),(kndx(i,j),
     1 j=1,2),(fcs(i,j),j=1,3),i=1,nseg)
 208  format(/(1x,8i6,3f6.1))
      write(10,211)(i,(jndx(i,j),j=1,2),(mndx(i,j),j=1,2),(px(i,j),
     1 j=1,2),(deg*xt(i,j),j=1,2),phcd(i),i=1,nbrn)
 211  format(/(1x,i3,4i5,2f12.6,2f10.2,2x,a))
      write(10,218)(i,(midx(i,j),j=1,2),(pux(i,j),j=1,2),
     1 i=1,max0(km(1),km(2)))
 218  format(/(1x,i3,2i5,2f12.6))
      write(10,212,iostat=ios)(i,pt(i),taut(i),deg*xa(i),
     1 cn*(xa(i)-xa(i+1)),(coef(j,i),j=1,5),i=1,nl)
 212  format(/(1x,i4,0p2f12.6,2f10.2,1p5d10.2))
      call retrns(10)
c
      call assign(10,2,'setbrn3.lis')
      do 20 nph=1,2
      mt(nph)=mt(nph)-3
      ku(nph)=ku(nph)-1
 20   km(nph)=km(nph)-1
c     icor=33  -  originally 32 records used as header in setbrn
c                 and 2 records used as header in remodl.
      icor=3
      do 14 nph=1,2
      m1=mt(nph)
      icor=icor-3
      do 14 i=2,m1
 14   ndx2(i,nph)=ndex(i,nph)+icor
      len1=ku(2)+km(2)
      len0=8*len1
      len2=5*nl
      write(10,*)'nseg nbrn mt ku km len len1',nseg,nbrn,mt,ku,km,len0,
     1 len1
      write(10,*)
      nasgr = len0
      write(6,*) 'reclength for direct access', nasgr
c++ 
c     write(6,*) 'enter model name'
c     read(5,*) cmodel
      nb=index(modnam,' ')-1
      if(nb.le.0) nb=len(modnam)
c     cnam1 = cmodel(1:nb)//'.tbl'
c     cnam2 = cmodel(1:nb)//'.hed'
      write(6,*) 'header file  :',modnam(1:nb)//'.hed'
      write(6,*) 'table file   :',modnam(1:nb)//'.tbl'
      call assign(nout,-2,modnam(1:nb)//'.hed')
c++ 
      write(nout) nasgr,nl,len2,xn,pn,tn,mt,nseg,nbrn,ku,km,fcs,nafl,
     1 indx,kndx
      write(nout) pm,zm,ndx2
      write(nout) pu,pux
      write(nout) phcd,px,xt,jndx
      write(nout) pt,taut
      write(nout) coef
      call retrns(nout)
c
      call dasign(nout,-2,modnam(1:nb)//'.tbl',nasgr)
      nrec = 0
      do 16 nph=1,2
      m1=mt(nph)
      n1=ku(nph)
      k1=km(nph)
      write(10,*)'nph m1 n1 k1',nph,m1,n1,k1
      do 16 m=2,m1
      if(ndex(m,nph).eq.ndex(m-1,nph)) go to 16
      read(nin,rec=ndex(m,nph))z0,n,(tmp(k,1),k=1,n),(tmp(k,2),k=1,n)
      write(10,*)'m nph ndex n',m,nph,ndex(m,nph),n
      k=0
      l=1
      do 17 i=1,n
      if(kuse(i,nph).lt.0) go to 17
      if(dabs(pux(l,nph)-pb(i)).gt.dtol) go to 18
      tmp(l,2)=tmp(i,2)
      l=l+1
 18   k=k+1
      tmp(k,1)=tmp(i,1)
 17   continue
      write(10,*)'k l nrec',k,l-1,nrec+1,ndx2(m,nph),sngl(tmp(1,1))
      if(k.ge.n1) go to 19
      k=k+1
      do 21 i=k,n1
 21   tmp(i,1)=0d0
 19   if(l.gt.k1) go to 23
      do 22 i=l,k1
 22   tmp(i,2)=0d0
 23   nrec=nrec+1
      write(nout,rec=nrec)(tmp(i,1),i=1,n1),(tmp(i,2),i=1,k1)
 16   continue
c
      call retrns(10)
      call retrns(nin)
      call retrns(nout)
      call vexit(0)
      end
c
      subroutine pdecx(n1,nph,xm,fac)
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl
      double precision xm(nsl1,2),ptol,pa,pax,plim
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      data ptol/.03/
c
      call collct(1,n1,xm(1,nph),fac*xmin)
      k=0
      plim=.7d0*pu(n1,nph)
      do 1 i=1,n1
      if(xm(i,nph).lt.0d0) go to 1
      if(pu(i,nph).lt.plim) go to 2
      if(pu(i,nph)-pu(k,nph).le.ptol) go to 2
      pa=pu(k,nph)+.75d0*(pu(i,nph)-pu(k,nph))
      pax=1d10
      m=0
      do 3 j=i1,i
      if(dabs(pu(j,nph)-pa).ge.pax) go to 3
      m=j
      pax=dabs(pu(j,nph)-pa)
 3    continue
      if(m.eq.i1.or.m.eq.i) go to 2
      k=k+1
      pu(k,nph)=pu(m,nph)
      xm(k,nph)=0d0
      kuse(m,nph)=1
 2    k=k+1
      i1=i
      pu(k,nph)=pu(i,nph)
      xm(k,nph)=xm(i,nph)
      kuse(i,nph)=1
 1    continue
      ku(nph)=k
      return
      end
c
      subroutine collct(i1,i2,x,xmn)
c
c $$$$$ calls varn $$$$$
c
      save
      double precision x(i2)
      data cn/6371./
c
      is=i1+1
      ie=i2-1
      if(ie.lt.is) return
      k1=i1
      var=0.
      m=0
      do 1 i=is,ie
      dx1=dabs(x(k1)-x(i))-xmn
      dx2=dabs(x(k1)-x(i+1))-xmn
      if(abs(dx2).ge.abs(dx1)) go to 2
      x(i)=-x(i)
      go to 1
 2    if(k1.le.i1) kb=i
      k1=i
      var=var+dx1*dx1
      m=m+1
 1    continue
      dx1=dabs(x(k1)-x(i2))-xmn
      var=var+dx1*dx1
      m=m+1
 7    if(m.le.1) return
      k1=i1
      k2=kb
      ks=kb+1
      nch=0
      do 8 i=ks,i2
      if(x(i).lt.0d0) go to 8
      k0=k1
      k1=k2
      k2=i
      var1=varn(x,k0,k1,k2,k1-1,xmn,var,m,m1)
      var2=varn(x,k0,k1,k2,k1+1,xmn,var,m,m2)
      if(amin1(var1/m1,var2/m2).ge.var/m) go to 6
      nch=nch+1
      x(k1)=-x(k1)
      if(var1/m1-var2/m2)3,4,5
 4    if(m1-m2)3,3,5
 3    k1=k1-1
      x(k1)=dabs(x(k1))
      var=var1
      m=m1
      go to 6
 5    k1=k1+1
      x(k1)=dabs(x(k1))
      var=var2
      m=m2
 6    if(k0.eq.i1) kb=k1
 8    continue
      if(nch.gt.0) go to 7
      return
      end
c
      function varn(x,k0,k1,k2,kt,xmn,var,m,mn)
c
c $$$$$ calls only library routines $$$$$
c
      save
      double precision x(k2)
c
      dx1=dabs(x(k0)-x(k1))-xmn
      dx2=dabs(x(k1)-x(k2))-xmn
      varn=var-dx1*dx1-dx2*dx2
      if(kt.le.k0.or.kt.ge.k2) go to 1
      dx1=dabs(x(k0)-dabs(x(kt)))-xmn
      dx2=dabs(dabs(x(kt))-x(k2))-xmn
      varn=varn+dx1*dx1+dx2*dx2
      mn=m
      return
 1    dx1=dabs(x(k0)-dabs(x(k2)))-xmn
      varn=varn+dx1*dx1
      mn=m-1
      return
      end
c
      subroutine layout
c
c   Layout contains the program for the desired travel-time segments 
c   implemented as calls to the mk_br entry points.  Each call does one 
c   segment (which may have many branches).
c
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision dir(3),cref(3),sref(3)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data dir,cref,sref/1d0,1d0,1d0,1d0,2d0,2d0,2d0,2d0,2d0/
c
c   Initialize variables.
      nseg=0
      nbrn=0
      nl=0
      do 1 j=1,3
      do 1 i=1,jseg
 1    fcs(i,j)=0.
      do 2 i=1,jout
 2    taut(i)=0d0
      do 3 j=1,2
      do 3 i=1,jbrn
 3    xt(i,j)=0d0
c
c   Do all of the segments.
c
c   P (up-going branch).
      print *,'Layout:  do Pup'
      call mkubr(ku(1),   +1)
c   P, Pdiff, and PKP.
      print *,'Layout:  do P and PKP'
      call mkdbr(1,lbb(1),-1,3,1,1,dir)
c   PKiKP.
      print *,'Layout:  do PKiKP'
      call mkrbr(2,       -1,2,1,1,dir)
c   pP.
      print *,'Layout:  do pP'
      call mkdbr(1,lbb(1),+1,3,1,1,dir)
c   sP.
      print *,'Layout:  do sP'
      call mkdbr(1,lbb(1),+2,3,1,1,dir)
c   pPKiKP.
      print *,'Layout:  do pPKiKP'
      call mkrbr(2,       +1,2,1,1,dir)
c   sPKiKP.
      print *,'Layout:  do sPKiKP'
      call mkrbr(2,       +2,2,1,1,dir)
c   PcP.
      print *,'Layout:  do PcP'
      call mkrbr(3,       -1,1,1,1,dir)
c   ScP.
      print *,'Layout:  do ScP'
      call mkrbr(3,       -2,1,2,1,dir)
c   SKP.
      print *,'Layout:  do SKP'
      call mkdbr(1,3,     -2,3,2,1,dir)
c   SKiKP.
      print *,'Layout:  do SKiKP'
      call mkrbr(2,       -2,2,2,1,dir)
c   PKKP.
      print *,'Layout:  do PKKP'
      call mkdbr(1,3,     -1,3,1,1,cref)
c   SKKP.
      print *,'Layout:  do SKKP'
      call mkdbr(1,3,     -2,3,2,1,cref)
c   PP and P'P'.
      print *,'Layout:  do PP, P''P'''
      call mkdbr(1,lbb(1),-1,3,1,1,sref)
c   S (up-going branch).
      print *,'Layout:  do Sup'
      call mkubr(ku(2),   +2)
c   S, Sdiff, and SKS.
      print *,'Layout:  do S and SKS'
      call mkdbr(1,lbb(2),-2,3,2,2,dir)
c   pS
      print *,'Layout:  do pS'
      call mkdbr(1,lbb(1),+1,3,2,2,dir)
c   sS
      print *,'Layout:  do sS'
      call mkdbr(1,lbb(2),+2,3,2,2,dir)
c   ScS
      print *,'Layout:  do ScS'
      call mkrbr(4,       -2,1,2,2,dir)
c   PcS
      print *,'Layout:  do PcS'
      call mkrbr(3,       -1,1,1,2,dir)
c   PKS
      print *,'Layout:  do PKS'
      call mkdbr(1,3,     -1,3,1,2,dir)
c   PKKS
      print *,'Layout:  do PKKS'
      call mkdbr(1,3,     -1,3,1,2,cref)
c   SKKS
      print *,'Layout:  do SKKS'
      call mkdbr(1,3,     -2,3,2,2,cref)
c   SS and S'S'.
      print *,'Layout:  do SS and S''S'''
      call mkdbr(1,lbb(2),-2,3,2,2,sref)
c   SP
      print *,'Layout:  do SP'
      call mkcbr(4,lbb(1),-2,1,2,1,sref)
c   PS
      print *,'Layout:  do PS'
      call mkcbr(4,lbb(1),-1,1,1,2,sref)
      return
      end
c
      subroutine mkdbr(l1,l2,isgn,lyr,nph,kph,fac)
c
c   Mkdbr sets up a simple refracted wave segment.  L1 and l2 point to the 
c   lbrk array of slowness break point pointers.  Note that the P and S 
c   break point arrays don't necessarily line up layer by layer.  This is 
c   not generally a problem as most phases need only worry about the 
c   pointer to the surface slowness for one wave type and a pointer 
c   somewhere in the core (which is constrained to be the same for both 
c   P and S).  Isgn is positive if the wave starts out going up and 
c   negative if the wave starts out going down.  Iabs(isng) is 1 if the 
c   wave starts out as a P wave and 2 if the wave starts out as an S wave.  
c   Lyr gives the number of major layers (mantle, outer core, and inner 
c   core) that the wave penetrates.  Nph and kph give the wave type (1 for 
c   P and 2 for S) on the down-going and up-going legs of the ray path 
c   respectively.  Fac is a three element array giving the number of 
c   repeats of the ray path in each major layer.  This scheme incorporates 
c   turning rays (e.g., P and S), turning rays reflected, but not 
c   converted at the surface (e.g., PP and SS), up-going rays reflected 
c   and/or converted at the surface into turning rays (e.g., pP and sP), 
c   turning rays converted during transmission through an interface (e.g., 
c   SKP and PKS), and rays which turn multiple times while reflecting from 
c   the bottom side of a layer (e.g., PKKP or SKKP).  Mkdbr does not 
c   include up-going rays (to the receiver), rays reflected from the top 
c   side of a discontinuity, or rays which are reflected and converted at 
c   the free surface.  See mkubr, mkrbr, and mkcbr respectively for 
c   routines which handle these types of rays.
c
      save
      include 'limits.inc'
      character*8 code,phcd,ks
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision fac(3)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/brnc/phcd(jbrn)
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data ks/'KKKKKKKK'/
c
c   Remember the programming as part of the final phase construction is 
c   done in depcor.
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
c   Using l1 and l2 to get the breakpoints has some shortcommings, 
c   particularly for converted phases.  It would be more general to 
c   have separate indicies for the breakpoints and the layers covered.
      if(l1.gt.1) kndx(nseg,1)=lbrk(l1-1,nph)
      kndx(nseg,2)=min0(min0(lbrk(l2,nph),lbrk(l2,kph)),
     1 lbrk(l2,iabs(isgn)))
      print *,'Mkdbr:  l1 l2 isgn lyr nph kph =',l1,l2,isgn,lyr,nph,kph
      print *,'Mkdbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 14 m=1,lyr
      fcs(nseg,m)=fac(m)
 14   xfc=amax1(xfc,fcs(nseg,m))
c
c   Set up the required slownesses, taus and distances.
c
      j=kndx(nseg,1)
      lz1=1
      lz2=1
c   Loop over the layers of interest.
      do 1 i=l1,l2
c   Be sure that the phase cuts off at the right place.
      l=min0(lbrk(i,nph),kndx(nseg,2))
c   Skip all total internal reflections.
      if(code(i,nph)(1:1).eq.'r'.or.j.ge.l) go to 1
c   Set the starting branch pointer.
      nbrn=nbrn+1
      nt=nl+1
      jndx(nbrn,1)=nt
c   Copy in the desired slownesses.
      do 2 k=j,l
      nl=nl+1
      pt(nl)=pb(k)
c   Add up the tau contributions.
      do 2 m=1,lyr
 2    taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
c   Take care of branch end pointers and slownesses.
      mndx(nbrn,1)=j
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(j)
      px(nbrn,2)=pb(l)
c   Add up distance contributions for the branch end points only.
      do 3 m=1,lyr
      xt(nbrn,1)=xt(nbrn,1)+fac(m)*(xp(j,m,nph)+xp(j,m,kph))
 3    xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
c   Take care of the contribution due to low velocity zones for the 
c   down-going leg(s).
      if(j.ne.lvz(lz1,nph)) go to 9
      do 11 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,nph)
 11   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,nph)
      taut(nt)=taut(nt)+fac(1)*taul(lz1,nph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz1,nph)
      lz1=lz1+1
c   Take care of the contributions due to low velocity zones for the 
c   up-going leg(s).
 9    if(j.ne.lvz(lz2,kph)) go to 10
      do 12 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,kph)
 12   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,kph)
      taut(nt)=taut(nt)+fac(1)*taul(lz2,kph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz2,kph)
      lz2=lz2+1
c   Decimate the slownesses if the branch is oversampled in distance.
 10   call pdect(jndx(nbrn,1),nl,j,iabs(isgn),xfc)
c   Set up the interpolation.
      call tauspl(jndx(nbrn,1),nl,pt,coef)
c   Remember the final branch end slowness value.
      jndx(nbrn,2)=nl
c
c   Take care of the branch name.  First, set up a default.
      phcd(nbrn)=code(i,nph)(2:2)//code(i,kph)(3:)
      if(idint(fac(1)+.5d0).gt.1) go to 5
      if(idint(fac(2)+.5d0).le.1) go to 4
c   Re-do the name if the ray is reflected from the underside of the 
c   core-mantle boundary.
      ind=idint(fac(2)-.5d0)
      phcd(nbrn)=code(i,nph)(2:2)//ks(1:ind)//code(i,kph)(3:)
      go to 4
c   Re-do the name if the ray is reflected from the surface.
 5    if(code(i,nph)(3:3).eq.' ') phcd(nbrn)=code(i,nph)(2:2)//
     1 code(i,kph)(2:)
      if(code(i,nph)(3:3).ne.' '.and.code(i,nph)(3:3).ne.'K')
     1  phcd(nbrn)=code(i,nph)(2:3)//code(i,kph)(2:)
      if(code(i,nph)(3:3).eq.'K') phcd(nbrn)=code(i,nph)(2:2)//''''//
     1 code(i,kph)(2:2)//''''//code(i,kph)(5:)
c   Take care .
 4    ind=max0(index(phcd(nbrn),'KSab'),index(phcd(nbrn),'S''ab'))
      if(phcd(nbrn)(1:1).eq.'S'.and.ind.gt.0) phcd(nbrn)(ind+2:ind+3)=
     1 'ac'
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
 1    j=l
      indx(nseg,2)=nl
      return
c
c   Mkubr handles up-going P and S.  L1 and isgn are as for mkdbr (except 
c   that l1 actually plays the role of l2 with the beginning break point 
c   assumed to be zero).  The other arguments are not needed.
c
      entry mkubr(l1,isgn)
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=0
      nafl(nseg,3)=0
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      l=kb(iabs(isgn))
      kndx(nseg,2)=l
      print *,'Mkubr:  l1 isgn =',l1,isgn
      print *,'Mkubr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      nbrn=nbrn+1
      jndx(nbrn,1)=nl+1
      do 6 k=1,l1
      nl=nl+1
      pt(nl)=pu(k,isgn)
 6    xa(nl)=0d0
      mndx(nbrn,1)=1
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(1)
      px(nbrn,2)=pb(l)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
      phcd(nbrn)=code(1,iabs(isgn))(2:2)
      indx(nseg,2)=nl
      return
c
c   Mkrbr handles reflected phases possibly with a conversion such as 
c   PcP, PcS, and PkiKP.  Arguments are as for mkdbr (except that l1 
c   actually plays the role of l2 with the beginning break point assumed 
c   to be zero).
c
      entry mkrbr(l1,isgn,lyr,nph,kph,fac)
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      l=min0(lbrk(l1,nph),lbrk(l1,kph))
      kndx(nseg,2)=l
      print *,'Mkrbr:  l1 isgn lyr nph kph =',l1,isgn,lyr,nph,kph
      print *,'Mkrbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 15 m=1,lyr
      fcs(nseg,m)=fac(m)
 15   xfc=amax1(xfc,fcs(nseg,m))
      if(lyr.ge.2) xfc=2.
c
      nbrn=nbrn+1
      jndx(nbrn,1)=nl+1
      do 7 k=1,l
      nl=nl+1
      pt(nl)=pb(k)
      do 7 m=1,lyr
 7    taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
      mndx(nbrn,1)=1
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(1)
      px(nbrn,2)=pb(l)
      do 8 m=1,lyr
 8    xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
      call pdect(jndx(nbrn,1),nl,1,iabs(isgn),xfc)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
      if(lyr.eq.1) phcd(nbrn)=code(l1,nph)(2:2)//'c'//code(l1,kph)(2:2)
      if(lyr.eq.2) phcd(nbrn)=code(l1,nph)(2:2)//code(l1,kph)(3:)
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
      indx(nseg,2)=nl
      return
c
c   Mkcbr handles phases reflected and converted at the surface such as 
c   PS and SP.  Arguments are as for mkdbr.
c
      entry mkcbr(l1,l2,isgn,lyr,nph,kph,fac)
      if(nph.gt.0.and.kph.gt.0.and.nph.ne.kph) go to 29
      print *,'Mkcbr:  bad call - nph kph =',nph,kph
      call vexit(1)
 29   nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      if(l1.gt.1) kndx(nseg,1)=min0(lbrk(l1,nph),lbrk(l1,kph))
      kndx(nseg,2)=min0(min0(lbrk(l2,nph),lbrk(l2,kph)),
     1 lbrk(l2,iabs(isgn)))
      print *,'Mkcbr:  l1 l2 isgn lyr nph kph =',l1,l2,isgn,lyr,nph,kph
      print *,'Mkcbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 16 m=1,lyr
      fcs(nseg,m)=fac(m)
 16   xfc=amax1(xfc,fcs(nseg,m))
c
      j=kndx(nseg,1)
      lz1=1
      lz2=1
      ik=l1
c
      print *,'Mkcbr:  start loop'
      do 17 in=l1,l2
 31   l=min0(lbrk(in,nph),kndx(nseg,2))
      if(code(in,nph)(1:1).eq.'r'.or.j.ge.l) go to 17
      l=min0(lbrk(ik,kph),kndx(nseg,2))
      if(code(ik,kph)(1:1).ne.'r'.and.j.lt.l.or.ik.ge.l2) go to 28
      j=max0(j,l)
      ik=ik+1
      go to 31
c
 28   if(lbrk(in,nph).le.lbrk(ik,kph)) go to 26
      l=min0(lbrk(ik,kph),kndx(nseg,2))
      print *,'kph:  kph ik j l code =',kph,ik,j,l,' ',code(ik,kph)
      isw=2
      go to 27
 26   l=min0(lbrk(in,nph),kndx(nseg,2))
      print *,'nph:  nph in j l code =',nph,in,j,l,' ',code(in,nph)
      isw=1
c
 27   nbrn=nbrn+1
      nt=nl+1
      jndx(nbrn,1)=nt
      do 18 k=j,l
      nl=nl+1
      pt(nl)=pb(k)
      do 18 m=1,lyr
 18   taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
      mndx(nbrn,1)=j
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(j)
      px(nbrn,2)=pb(l)
      do 19 m=1,lyr
      xt(nbrn,1)=xt(nbrn,1)+fac(m)*(xp(j,m,nph)+xp(j,m,kph))
 19   xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
      if(j.ne.lvz(lz1,nph)) go to 20
      do 21 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,nph)
 21   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,nph)
      taut(nt)=taut(nt)+fac(1)*taul(lz1,nph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz1,nph)
      lz1=lz1+1
 20   if(j.ne.lvz(lz2,kph)) go to 22
      do 23 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,kph)
 23   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,kph)
      taut(nt)=taut(nt)+fac(1)*taul(lz2,kph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz2,kph)
      lz2=lz2+1
 22   call pdect(jndx(nbrn,1),nl,j,iabs(isgn),xfc)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
c
      if(code(in,nph)(3:3).eq.' ') phcd(nbrn)=code(in,nph)(2:2)//
     1 code(ik,kph)(2:)
      if(code(in,nph)(3:3).ne.' '.and.code(in,nph)(3:3).ne.'K')
     1  phcd(nbrn)=code(in,nph)(2:3)//code(ik,kph)(2:)
      if(code(in,nph)(3:3).eq.'K') phcd(nbrn)=code(in,nph)(2:2)//''''//
     1 code(ik,kph)(2:2)//''''//code(ik,kph)(5:)
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
      print *,'phcd:  in ik phcd =',in,ik,' ',phcd(nbrn)
      if(isw.le.1) go to 17
      ik=ik+1
      j=max0(j,l)
      go to 31
 17   j=max0(j,l)
      indx(nseg,2)=nl
      return
      end
c
      subroutine pdect(i1,i2,j1,nph,fac)
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision h1,h2,hh
      dimension ib(2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
c
      xmn=fac*xmin
      isg=1
      do 1 i=1,2
      ib(i,1)=i1
 1    ib(i,2)=i2
      ii=i1+1
      ie=i2-1
      xa(i1)=xt(nbrn,1)
      do 2 i=ii,ie
      h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
 2    xa(i)=-(h2*taut(i-1)-(h2+h1)*taut(i)+h1*taut(i+1))/hh
      xa(i2)=xt(nbrn,2)
      do 3 i=ii,ie
      if((xa(i+1)-xa(i))*(xa(i)-xa(i-1)).gt.0d0) go to 3
      isg=2
      ib(1,2)=i-2
      ib(2,1)=i+2
 3    continue
      do 4 it=1,isg
 4    call collct(ib(it,1),ib(it,2),xa,xmn)
      k=i1-1
      j=j1
      do 5 i=i1,i2
      if(xa(i).lt.0d0) go to 5
      k=k+1
      pt(k)=pt(i)
      taut(k)=taut(i)
      xa(k)=xa(i)
      kuse(j,nph)=1
 5    j=j+1
      if(k.eq.nl) return
      ii=k+1
      do 6 i=ii,nl
 6    taut(i)=0d0
      nl=k
      return
      end
c
      subroutine kseq
c
c   Kseq makes a correspondence between model slownesses in array pb and 
c   the subset of the same slownesses used for sampling tau which are 
c   stored in pu (separate sets for P and S).  The net result is to 
c   translate the kndx pointers to critical slowness values (bounding 
c   branches actually implemented) from pointing into pb to pointing 
c   into pu.
c   
c
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl
      dimension kl(2),kk(jseg,2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      data kl/0,0/
c
c   Compile a sorted list of unique kndx values in the first column of 
c   kk.
c
      do 1 i=1,nseg
      nph=iabs(nafl(i,1))
      k=kl(nph)
      do 2 j=1,2
      if(k.le.0) go to 3
      do 4 m=1,k
      if(kk(m,nph,1)-kndx(i,j))4,2,5
 4    continue
 3    k=k+1
      kk(k,nph,1)=kndx(i,j)
      go to 2
 5    do 6 l=k,m,-1
 6    kk(l+1,nph,1)=kk(l,nph,1)
      k=k+1
      kk(m,nph,1)=kndx(i,j)
 2    continue
 1    kl(nph)=k
c
c   Make the correspondence between pb and pu for each kndx and save it 
c   in the second column of kk.
c
      do 7 nph=1,2
      n1=ku(nph)
      k=1
      ki=kk(k,nph,1)
      do 8 i=1,n1
      if(pu(i,nph)-pb(ki))8,9,10
 9    kk(k,nph,2)=i
      if(k.ge.kl(nph)) go to 7
      k=k+1
      ki=kk(k,nph,1)
 8    continue
 10   write(*,100)ki,pb(ki),nph
 100  format(1x,'Kseq:  pb(',i3,') =',f7.4,' not found in pu(*,',i1,
     1 ').')
      call vexit(1)
 7    continue
c
c   Replace each kndx pb pointer with the corresponding pu pointer.
c
      do 11 i=1,nseg
      nph=iabs(nafl(i,1))
      k=kl(nph)
      do 11 j=1,2
      do 12 m=1,k
      if(kk(m,nph,1)-kndx(i,j))12,11,13
 12   continue
 13   write(*,101)kndx(i,j)
 101  format(1x,'Kseq:  kndx value',i4,' not translated.')
      call vexit(1)
 11   kndx(i,j)=kk(m,nph,2)
      return
      end
c
      subroutine mseq
c         partial reordering of tables
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pux
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/msqc/pux(jbrn,2),km(2),midx(jbrn,2)
      data km/0,0/
c
      is=1
      do 1 i=1,nbrn
 8    if(jndx(i,2).le.indx(is,2)) go to 7
      is=is+1
      go to 8
 7    nph=iabs(nafl(is,1))
      k=km(nph)
      do 2 j=1,2
      if(k.le.0) go to 3
      do 4 m=1,k
      if(midx(m,nph)-mndx(i,j))4,2,5
 4    continue
 3    k=k+1
      midx(k,nph)=mndx(i,j)
      pux(k,nph)=px(i,j)
      go to 2
 5    do 6 l=k,m,-1
      midx(l+1,nph)=midx(l,nph)
 6    pux(l+1,nph)=pux(l,nph)
      k=k+1
      midx(m,nph)=mndx(i,j)
      pux(m,nph)=px(i,j)
 2    continue
 1    km(nph)=k
      return
      end
