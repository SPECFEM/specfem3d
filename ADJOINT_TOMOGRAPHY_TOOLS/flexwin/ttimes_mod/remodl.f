      program remodl
      save
      include 'limits.inc'
      character*8 code
      character*20 modnam
c     logical ldum
      dimension mt(2),kb(2),lt(2),lvz(nlvz0,2)
      double precision taup(nsl1,2),xp(nsl1,2),ttau,tx,pb,pm,zm,pmj,pmi,
     1 zmj,zmi,zmax,zic,zoc,taul(nlvz0,2),xl(nlvz0,2)
      common/pgridc/delx(2),dpmax(2),drmax(2),pp(nsl1,2),xx(nsl1,2),
     1 rr(nsl1,2)
      common/zgridc/pb(nsl1),pm(ndp1,2),zm(ndp1,2),zm0(ndp1,2),
     1 ndex(ndp1,2),n
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      common/emodc/r(nmd0),vp(nmd0),vs(nmd0)
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
      common/crtptc/ucrt(ncp0),nc
      common/roughc/ric,roc
      data nout,delx,dpmax,drmax,tol,dtol,xtol,dmax/1,2*200.,2*.01,
     1 2*75.,1e-6,1d-6,1.,800./
c
      write (*,*) 'DEBUG - calling rough'
      call rough(50.,nk,ncr,modnam)
      a0=r(1)
      call assign(10,2,'remodl1.lis')
      write(10,200)nk,ncr,a0,(i,r(i),vp(i),vs(i),i=1,nk)
 200  format(/1x,2i5,f10.2/(i5,f10.2,2f10.4))
      pn=vs(1)
      xn=1/r(1)
      tn=pn*xn
      do 1 i=1,nk
      rn=r(i)*xn
      z(i)=alog(rn)
      u(i,1)=rn*pn/vp(i)
      if(i.le.ncr) u(i,2)=rn*pn/vs(i)
 1    if(i.gt.ncr) u(i,2)=u(i,1)
      write(10,201)(i,r(i),a0*z(i),(u(i,j),j=1,2),i=1,nk)
 201  format(/(1x,i5,2f10.2,2f12.6))
      write(10,202)delx,dpmax,drmax
 202  format(/1x,'delx =',2f10.2/1x,'dpmax -',2f12.6/1x,'drmax =',
     1 2f10.2)
      delx(1)=xn*delx(1)
      delx(2)=xn*delx(2)
      xtol=xtol/a0
c
      write (*,*) 'DEBUG - calling findcp'
      call findcp(ncr)
      write(10,203)(i,ucrt(i),i=1,nc)
 203  format(/1x,'critical points'/(1x,i5,f12.6))
c     call assign(11,-2,'crtpt.dat')
c     write(11)nc,(ucrt(i),i=1,nc)
c     call retrns(11)
c
      write(*,*) 'DEBUG Calling pgrid with', 1, n1, xtol
      call pgrid(1,n1,xtol)
c     write(10,204)pp(1,1),xx(1,1),rr(1,1)
c204  format(/5x,'1',f12.6,12x,f12.6,10x,f10.2)
c     write(10,205,iostat=ios)(i,pp(i,1),(pp(i-1,1)-pp(i,1)),xx(i,1),
c    1 a0*(xx(i,1)-xx(i-1,1)),rr(i,1),rr(i-1,1)-rr(i,1),i=2,n1)
c205  format(1x,i5,3f12.6,3f10.2)
      write(*,*) 'DEBUG Calling pgrid with', 2, n2, xtol
      call pgrid(2,n2,xtol)
c     write(10,204)pp(1,2),xx(1,2),rr(1,2)
c     write(10,205,iostat=ios)(i,pp(i,2),(pp(i-1,2)-pp(i,2)),xx(i,2),
c    1 a0*(xx(i,2)-xx(i-1,2)),rr(i,2),rr(i-1,2)-rr(i,2),i=2,n2)
      write(*,*) 'DEBUG Done calling pgrid with', 2, n2, xtol
c
      n=0
      j1=1
      k1=1
      do 2 i=1,nc
      ii=i+1
      if(pp(1,1).ge.ucrt(i)-tol) go to 3
      j=1
      go to 4
 3    do 5 j=j1,n1
      if(abs(pp(j,1)-ucrt(ii)).le.tol) go to 4
 5    continue
      j=n1
 4    do 6 k=k1,n2
      if(abs(pp(k,2)-ucrt(ii)).le.tol) go to 7
 6    continue
      k=n2
 7    if(j-j1.le.k-k1) go to 8
      do 9 l=j1,j
      n=n+1
 9    pb(n)=pp(l,1)
      go to 10
 8    do 11 l=k1,k
      n=n+1
 11   pb(n)=pp(l,2)
 10   j1=j+1
 2    k1=k+1
      call retrns(10)
      call assign(10,2,'remodl2.lis')
      write(10,211)pb(1)
 211  format(/5x,'1',f12.6)
      write(10,212,iostat=ios)(i,pb(i),(pb(i-1)-pb(i)),i=2,n)
 212  format(1x,i5,0pf12.6,1pe12.2)
      call efe8(n,pb)
      n1=n
      n=n-1
c
      call zgrid(1,mt(1))
      mm=mt(1)
      write(10,206,iostat=ios)(i,pm(i,1),zm(i,1),a0*(zm0(i,1)-
     1 zm0(i+1,1)),a0*(zm(i,1)-zm(i+1,1)),ndex(i,1),i=1,mm)
 206  format(/(1x,i5,0p2f12.6,f10.2,1pe12.4,i5))
      write(10,208,iostat=ios)mm+1,pm(mm+1,1),zm(mm+1,1),ndex(mm+1,1)
 208  format(1x,i5,0p2f12.6,22x,i5)
      call zgrid(2,mt(2))
      mm=mt(2)
      write(10,206,iostat=ios)(i,pm(i,2),zm(i,2),a0*(zm0(i,2)-
     1 zm0(i+1,2)),a0*(zm(i,2)-zm(i+1,2)),ndex(i,2),i=1,mm)
      write(10,208,iostat=ios)mm+1,pm(mm+1,2),zm(mm+1,2),ndex(mm+1,2)
c
c   Set up break pointers.
      call brkpts(mt(1),1)
      write(10,209)(lbrk(i,1),i=1,lbb(1))
 209  format(/(1x,20i5))
      write(10,210)(i,code(i,1),i=1,lcb(1))
 210  format(/(1x,i5,2x,a))
      kb(1)=lbrk(lbb(1),1)
      call brkpts(mt(2),2)
      write(10,209)(lbrk(i,2),i=1,lbb(2))
      write(10,210)(i,code(i,2),i=1,lcb(2))
      kb(2)=lbrk(lbb(2),2)
      write(10,*)
      write(10,*)'n1 kb',n1,kb
c
      ndasr  = 8*(1+2*kb(2))+4
      call dasign(nout,-2,'remodl.tbl',ndasr)
      write(10,*) 'reclength for dasign:', ndasr
      n1=kb(2)
      zmax=alog((a0-dmax)*xn)
      zic=alog(ric*xn)
      zoc=alog(roc*xn)
      zlim=zmax
      write(10,*)
      write(10,*)'zmax zoc zic',zmax,zoc,zic
      write(10,*)
c
c   Loop over phases.
c
      nrec=0
      do 17 nph=1,2
      j=1
      lz=0
      n1=kb(nph)
      do 16 i=1,n1
      taup(i,nph)=0d0
 16   xp(i,nph)=0d0
      taup(n1,nph)=tn*1d-6
      xp(n1,nph)=xn*1d-6
      n=n1-1
      mm=mt(nph)+1
      ndex(1,nph)=-1
      zmi=zm(1,nph)
      pmi=pm(1,nph)
c
c   Loop over model slownesses.
c
      do 18 i=2,mm
      zmj=zmi
      zmi=zm(i,nph)
      pmj=pmi
      pmi=pm(i,nph)
      if(dabs(zmj-zmi).le.0d0) go to 19
c
c   Collect the tau and x integrals.
      do 20 k=1,n
      if(pmi.lt.pb(k)) go to 21
      call tauint(pb(k),pmj,pmi,zmj,zmi,ttau,tx)
c     print *,'int',i,k,pb(k),pmj,pmi,dz,ttau,tx
      taup(k,nph)=taup(k,nph)+ttau
 20   xp(k,nph)=xp(k,nph)+tx
      go to 22
 21   n=k-1
 22   if(n.le.1) go to 23
      if(pb(n).eq.pb(n-1)) n=n-1
 23   if(zmj.lt.zlim) go to 18
      j=j+1
      if(zmj.ge.zmax) nrec=nrec+1
      zm(j,nph)=zmi
      pm(j,nph)=pmi
      ndex(j,nph)=nrec
      if(zmj.lt.zmax) go to 18
      write(10,*)'lev1',j,n,nph,sngl(pm(j,nph)),sngl(zm(j,nph))
      write(nout,rec=nrec)zmi,n,(taup(k,nph),k=1,n),(xp(k,nph),k=1,n)
      go to 18
 19   if(dabs(zmi-zoc).gt.dtol.and.dabs(zmi-zic).gt.dtol) go to 25
      if(dabs(zmi-zm(j,nph)).le.dtol) go to 18
      j=j+1
      nrec=nrec+1
      zm(j,nph)=zmi
      pm(j,nph)=pmi
      ndex(j,nph)=nrec
      write(10,*)'lev2',j,n,nph,sngl(pm(j,nph)),sngl(zm(j,nph))
      write(nout,rec=nrec)zmi,n1,(taup(k,nph),k=1,n1),(xp(k,nph),k=1,n1)
      go to 26
 25   if(zmi.lt.zmax) go to 26
      if(dabs(zmi-zm(j-1,nph)).gt.dtol) j=j+1
      zm(j,nph)=zmi
      pm(j,nph)=pmi
      ndex(j,nph)=ndex(j-1,nph)
 26   if(pmi.le.pmj) go to 18
      if(lz.le.0) go to 27
      if(lvz(lz,nph).eq.n) go to 18
 27   lz=lz+1
      lvz(lz,nph)=n
      taul(lz,nph)=taup(n,nph)
      xl(lz,nph)=xp(n,nph)
      write(10,*)'lvz ',lz,n,nph,sngl(taul(lz,nph)),sngl(xl(lz,nph))
 18   continue
      j=j+1
      nrec=nrec+1
      zm(j,nph)=zmi
      pm(j,nph)=pmi
      ndex(j,nph)=nrec
      write(10,*)'lev3',j,n,nph,sngl(pm(j,nph)),sngl(zm(j,nph))
      write(10,*)
      write(nout,rec=nrec)zmi,n1,(taup(k,nph),k=1,n1),(xp(k,nph),k=1,n1)
      mt(nph)=j
      lt(nph)=lz
      if(lz.le.1) go to 34
      call efe4(lz,lvz(1,nph))
      call efe8(lz,taul(1,nph))
      call efe8(lz,xl(1,nph))
c
 34   if(nph.ge.2) go to 17
      do 30 i=1,mm
      if(zm(i,1).lt.zmax) go to 31
 30   continue
      i=mm
 31   plim=pm(i,1)
      write(10,*)'i zmax plim zm(i,1) =',i,zmax,plim,zm(i,1)
      do 32 i=1,mm
      if(pm(i,2).le.plim) go to 33
 32   continue
      i=mm
 33   zlim=zm(i,2)
      write(10,*)'i plim zlim =',i,plim,zlim
c
 17   continue
      call retrns(nout)
c
      call assign(nout,-2,'remodl.hed')
      write(nout)ndasr,modnam,zmax,zoc,zic,kb,(pb(i),i=1,n1),
     1 mt,lt,lbb,lcb,xn,pn,tn
      write(nout)((lbrk(i,nph),i=1,lbb(nph)),(code(i,nph),
     1 i=1,lcb(nph)),(zm(i,nph),pm(i,nph),ndex(i,nph),i=1,mt(nph)),
     2 (lvz(i,nph),taul(i,nph),xl(i,nph),i=1,lt(nph)),nph=1,2)
      call retrns(nout)
      call retrns(10)
      call vexit(0)
      end
c
      subroutine rough(dr,n,ncr,modnam)
c
c $$$$$ calls emdld and emdlv $$$$$
c
c   Rough provides a rough interpolation of the earth model available
c   through routine emdlv.  The model radii, compressional velocity, and
c   shear velocity are provided in arrays r, vp, and vs respectively.
c   Between first order discontinuities available through routine emdld,
c   the radii are equally spaced as close to spacing dr as possible.  The
c   number of radii used is returned in variable n and the index of the
c   core-mantle radius is returned in ncr.  Note that emdld returns the
c   radii of the discontinuities from the center of the earth out while
c   rough returns the model from the surface in.  Also note that in the
c   model returned by rough, each discontinuity will be represented by
c   two model ajacent model points with the same radius.
c
      save
      include 'limits.inc'
      character*(*) modnam
      dimension rd(30)
      common/emodc/r(nmd0),vp(nmd0),vs(nmd0)
      common/roughc/ric,roc
      data tol,vtol/1e-6,2e-5/
c
c   Get the radii of model discontinuities.
      call emdld(np,rd,modnam)
      write(*,*) 'DEBUG',np,rd,modnam
c   Save the radii of the inner core-outer core and core-mantle boundaries
c   respectively.
      ric=rd(1)
      roc=rd(2)
c
c   Begin the interpolation.
c
      n=0
      i=np-1
      r1=rd(np)
c   Loop over each layer (between two discontinuities).
 1    r0=r1
      r1=.001
      if(i.gt.0) r1=rd(i)
      l=(r0-r1)/dr-.5
      dx=(r0-r1)/(l+1)
c   Set the outer most point of the layer.
      n=n+1
      r(n)=r0
      call emdlv(r0*(1.-tol),vp(n),vs(n))
c   Check for continuity across an apparant discontinuity.
      if(n.le.1) go to 2
c   If vp is close to continuous, force it.
      if(abs(vp(n-1)-vp(n)).gt.vtol*vp(n)) go to 6
      vp(n-1)=.5*(vp(n-1)+vp(n))
      vp(n)=vp(n-1)
c   If vs is close to continuous, force it.
 6    if(abs(vs(n-1)-vs(n)).gt.vtol*vs(n)) go to 2
      vs(n-1)=.5*(vs(n-1)+vs(n))
      vs(n)=vs(n-1)
c   Make P and S velocity the same if we are in a fluid.
 2    if(i.le.1) vs(n)=vp(n)
c   Interpolate the model throughout the layer.
      if(l.le.0) go to 4
      do 3 j=1,l
      n=n+1
      r(n)=r0-j*dx
      call emdlv(r(n),vp(n),vs(n))
c   Make P and S velocity the same if we are in a fluid.
 3    if(i.le.1) vs(n)=vp(n)
c   Set the inner most point of the layer.
 4    n=n+1
      r(n)=r1
      call emdlv(r1*(1.+tol),vp(n),vs(n))
c   Make P and S velocity the same if we are in a fluid.
      if(i.le.1) vs(n)=vp(n)
c   Set the index to the core-mantle radius.
      if(i.eq.2) ncr=n
      i=i-1
c   If there is another layer, go do it.
      if(i.ge.0) go to 1
      return
      end
c
      subroutine findcp(ncr)
c
c $$$$$ calls crtpt $$$$$
c
c        find critical points
      save
      include 'limits.inc'
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
      common/crtptc/ucrt(ncp0),nc
      data tol/1e-6/
c
      ifl=1
      kfl=1
      j=1
      call crtpt(1,1)
      call crtpt(1,2)
      do 5 i=2,nk
      if(abs(z(j)-z(i)).gt.tol) go to 6
      call crtpt(j,1)
      call crtpt(i,1)
      if(j.le.ncr) call crtpt(j,2)
      if(i.le.ncr) call crtpt(i,2)
      go to 5
 6    go to (7,8),ifl
 7    if(u(i,1).le.u(j,1)) go to 2
      ifl=2
      call crtpt(j,1)
      go to 2
 8    if(u(i,1).ge.u(j,1)) go to 2
      ifl=1
      call crtpt(j,1)
 2    if(i.gt.ncr) go to 5
      go to (3,4),kfl
 3    if(u(i,2).le.u(j,2)) go to 5
      kfl=2
      call crtpt(j,2)
      go to 5
 4    if(u(i,2).ge.u(j,2)) go to 5
      kfl=1
      call crtpt(j,2)
 5    j=i
      ucrt(nc+1)=0.
      return
      end
c
      subroutine crtpt(k,nph)
c
c $$$$$ calls no other routine $$$$$
c
c   For each critical point (slowness corresponding to a first order
c   discontinuity) save slowness u(k,nph) in array ucrt and sort it
c   into descending order.  Note that k indexes an equivalent depth and
c   nph indexes either P or V slownesses.
c
      save
      include 'limits.inc'
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
      common/crtptc/ucrt(ncp0),nc
      data nc/0/,tol/1e-6/
c   Eliminate duplicates.
      if(nc.le.0) go to 3
      do 2 i=1,nc
      if(abs(ucrt(i)-u(k,nph)).gt.tol) go to 2
      write(10,*)'duplicate critical value eliminated:',u(k,nph)
      return
 2    continue
c
 3    nc=nc+1
      ucrt(nc)=u(k,nph)
      if(nc.le.1) return
      j=nc
      do 1 i=2,nc
      if(ucrt(j-1).ge.ucrt(j)) return
      utmp=ucrt(j-1)
      ucrt(j-1)=ucrt(j)
      ucrt(j)=utmp
 1    j=j-1
      return
      end
c
      subroutine pgrid(nph,n0,xtol)
      save 
      parameter(mxtmp=200)
      include 'limits.inc'
      dimension ps(mxtmp),xs(mxtmp),pb(3),xb(3)
      common/pgridc/delx(2),dpmax(2),drmax(2),pp(nsl1,2),xx(nsl1,2),
     1 rr(nsl1,2)
      common/emodc/r(nmd0),vp(nmd0),vs(nmd0)
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
      common/crtptc/ucrt(ncp0),nc
      data tol/1e-5/
c
      a0=1./xn
      do 1 i=1,nc
      if(abs(u(1,nph)-ucrt(i)).le.tol) go to 2
 1    continue
 2    ic=i
c
      n=0
      j=ic+1
      do 20 i=ic,nc 
      write(10,210)'i',i,ucrt(i),ucrt(j)
 210  format(10x,a,i5,2f12.6)
c
      ifl=1
      ps(1)=ucrt(i)
      xs(1)=xmod(ucrt(i),nph,1,nr,r0)
      x0=xmod(ucrt(j),nph,-1,nr,r1)
      l=max0(int(abs(x0-xs(1))/delx(nph)+.8),1)
      du=(ps(1)-ucrt(j))/(l*l)
c     write(10,*)'overflow???  x0 x1 l du nr r1',
c    1 xs(1),x0,l+1,du,nr,r1
      l=l-1
      do 3 k=1,l
      ps(k+1)=ps(1)-k*k*du
 3    xs(k+1)=xmod(ps(k+1),nph,1,nr,r2)
      l=l+2
      ps(l)=ucrt(j)
      xs(l)=x0
c     write(10,*)'***** i l',i,l
c
 19   x0=xs(2)-xs(1)
      do 4 k=3,l
      if((xs(k)-xs(k-1))*x0.le.0.) go to 5
 4    continue
      go to 6
c
 5    ifl=ifl+1
      r2=r1
      write(10,210)'caustic'
      k=k-1
      kk=k-2
      do 7 m=1,3
      kk=kk+1
      pb(m)=ps(kk)
 7    xb(m)=xs(kk)
      sgn=sign(1.,.5*(xb(1)+xb(3))-xb(2))
 11   u0=.5*(pb(1)+pb(2))
      x0=xmod(u0,nph,1,nr,r1)
      if(sgn*(xb(2)-x0).lt.0.) go to 8
      pb(3)=pb(2)
      xb(3)=xb(2)
      pb(2)=u0
      xb(2)=x0
      go to 9
 8    pb(1)=u0
      xb(1)=x0
      u0=.5*(pb(2)+pb(3))
      x0=xmod(u0,nph,1,nr,r1)
      if(sgn*(xb(2)-x0).lt.0.) go to 10
      pb(1)=pb(2)
      xb(1)=xb(2)
      pb(2)=u0
      xb(2)=x0
      go to 9
 10   pb(3)=u0
      xb(3)=x0
 9    if(abs(xb(3)-xb(1)).gt.xtol) go to 11
      ps(k)=pb(2)
      xs(k)=xb(2)
      lsav=l
      l=k
c
 6    if(i.ne.ic) go to 21
      n=n+1
      pp(n,nph)=ps(1)
      xx(n,nph)=xs(1)
      rr(n,nph)=r0
      write(10,201)'first ',n,pp(n,nph),xx(n,nph),rr(n,nph)
 21   k=max0(int(abs(xs(l)-xs(1))/delx(nph)+.8),1)
      dx=(xs(l)-xs(1))/k
      x0=xs(1)
      rsav=r0
      mm=2
 12   x0=x0+dx
      n=n+1
      if(abs(x0-xs(l)).gt.tol) go to 24
      pp(n,nph)=ps(l)
      xx(n,nph)=xs(l)
      rr(n,nph)=r1
      go to 22
 24   do 14 kk=mm,l
      if((x0-xs(kk))*(x0-xs(kk-1)).le.0.) go to 15
 14   continue
      kk=l
 15   mm=kk
      call finrng(x0,nph,n,ps(kk),ps(kk-1),xs(kk),xs(kk-1),xtol,nr)
      write(10,201)'sol   ',n,pp(n,nph),xx(n,nph),rr(n,nph),nr,
     1 pp(n-1,nph)-pp(n,nph),a0*(xx(n,nph)-xx(n-1,nph)),rsav-rr(n,nph)
 201  format(1x,a6,i5,2f10.6,f9.2,i5,f10.6,2f9.2)
 22   if(abs(pp(n,nph)-pp(n-1,nph)).le.dpmax(nph)) go to 16
      ll=max0(int(abs(ps(l)-pp(n-1,nph))/dpmax(nph)+.99),1)
      pp(n,nph)=pp(n-1,nph)+(ps(l)-pp(n-1,nph))/ll
      xx(n,nph)=xmod(pp(n,nph),nph,1,nr,rr(n,nph))
      write(10,201)' dpmax',n,pp(n,nph),xx(n,nph),rr(n,nph),nr,
     1 pp(n-1,nph)-pp(n,nph),a0*(xx(n,nph)-xx(n-1,nph)),rsav-rr(n,nph)
      k=max0(int(abs(xs(l)-xx(n,nph))/delx(nph)+.8),1)
      dx=(xs(l)-xx(n,nph))/k
      x0=xx(n,nph)
      mm=2
 16   if(abs(rr(n,nph)-rsav).le.drmax(nph)) go to 23
      rnew=rsav-drmax(nph)
 26   if(rnew.le.r(nr)) go to 25
      nr=nr-1
      go to 26
 25   if(nr.lt.nk) du=abs(pp(n-1,nph)-u(nr,nph)*(rnew/r(nr))**
     1 (alog(u(nr+1,nph)/u(nr,nph))/alog(r(nr+1)/r(nr))))
      if(nr.ge.nk) du=abs(pp(n-1,nph)-u(nk,nph)*rnew/r(nk))
      ll=max0(int(abs(ps(l)-pp(n-1,nph))/du+.99),1)
      pp(n,nph)=pp(n-1,nph)+(ps(l)-pp(n-1,nph))/ll
      xx(n,nph)=xmod(pp(n,nph),nph,1,nr,rr(n,nph))
      write(10,201)' drmax',n,pp(n,nph),xx(n,nph),rr(n,nph),nr,
     1 pp(n-1,nph)-pp(n,nph),a0*(xx(n,nph)-xx(n-1,nph)),rsav-rr(n,nph)
      k=max0(int(abs(xs(l)-xx(n,nph))/delx(nph)+.8),1)
      dx=(xs(l)-xx(n,nph))/k
      x0=xx(n,nph)
      mm=2
 23   rsav=rr(n,nph)
      if(abs(x0-xs(l)).gt.tol) go to 12
c
 17   write(10,201)'end   ',n,pp(n,nph),xx(n,nph),rr(n,nph),nr,
     1 pp(n-1,nph)-pp(n,nph),a0*(xx(n,nph)-xx(n-1,nph)),rr(n-1,nph)-
     2 rr(n,nph)
      ifl=ifl-1
      if(ifl.le.0) go to 20
      k=0
      do 18 m=l,lsav
      k=k+1
      ps(k)=ps(m)
 18   xs(k)=xs(m)
c     write(10,*)'+++++ i l lsav k',i,l,lsav,k
      l=lsav-l+1
      r0=r1
      r1=r2
      go to 19
c
 20   j=j+1
      n0=n
      return
      end
c
      function xmod(pk,nph,ipart,nr,rb)
c        partial tau integrals
      save 
      include 'limits.inc'
      common/emodc/r(nmd0),vp(nmd0),vs(nmd0)
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
c
      if(pk.gt.0.) go to 4
      x=-1.5707963
      nr=nk
      rb=0.
      go to 3
 4    p2=pk*pk
      x=0.
      j=1
      do 1 i=2,nk
      if(pk.gt.u(i,nph)) go to 2
      if(u(j,nph).ne.u(i,nph)) x=x+(z(i)-z(j))*
     1 (acos(pk/u(j,nph))-acos(pk/u(i,nph)))/alog(u(j,nph)/u(i,nph))
      if(pk.ne.u(i,nph).or.ipart.ge.0) go to 1
      nr=i
      rb=r(i)
      go to 3
 1    j=i
      nr=nk
      rb=r(nk)*pk/u(nk,nph)
      zb=alog(rb*xn)
      if(u(nk,nph).gt.pk) x=x+(zb-z(nk))*acos(pk/u(nk,nph))/
     1 alog(u(nk,nph)/pk)
      go to 3
 2    nr=j
      rb=r(j)*(pk/u(j,nph))**(alog(r(i)/r(j))/alog(u(i,nph)/u(j,nph)))
      zb=alog(rb*xn)
      if(u(j,nph).gt.pk) x=x+(zb-z(j))*acos(pk/u(j,nph))/
     1 alog(u(j,nph)/pk)
 3    xmod=-2.*x
      return
      end
c
      subroutine finrng(xg,nph,n,p0,p1,x0,x1,xtol,nr)
c        iteration on range
c
c $$$$$ calls xmod $$$$$
c
c   Function find0 returns x0 for which f(x0) = 0 where f(x) is a
c   function supplied by the user (specified external in the calling
c   routine).  X1 and x2 are the starting trial points.  Function
c   evaluations are made until x0 is determined to a relative precision
c   of eps by a process of inverse iterative interpolation using
c   Aitken's method.
c
c                                                     -rpb
      save 
      include 'limits.inc'
      character*54 msg
      common/pgridc/delx(2),dpmax(2),drmax(2),pp(nsl1,2),xx(nsl1,2),
     1 rr(nsl1,2)
      common/findc/x(20),y(20),m
c
      m=1
      if((x1-xg)*(x0-xg).le.0.) go to 7
      stop 'Root not bracketed.'
c   If the limits are already close together, there is no point in
c   iterating.
 7    if(abs(x0-x1).gt.xtol) go to 3
      pp(n,nph)=.5*(p0+p1)
      xx(n,nph)=xmod(pp(n,nph),nph,1,nr,rr(n,nph))
      return
c   Set up iteration with first two trial points, x1 and x2.
 3    y(1)=x0-xg
      if(abs(y(1)).gt.xtol) go to 4
      pp(n,nph)=p0
      xx(n,nph)=xmod(pp(n,nph),nph,1,nr,rr(n,nph))
      return
 4    x(1)=p0
      yi=x1-xg
      if(abs(yi).gt.xtol) go to 5
      pp(n,nph)=p1
      xx(n,nph)=xmod(pp(n,nph),nph,1,nr,rr(n,nph))
      return
 5    if(y(1).gt.yi) go to 8
      ps0=p0
      xs0=y(1)
      ps1=p1
      xs1=yi
      go to 9
 8    ps0=p1
      xs0=yi
      ps1=p0
      xs1=y(1)
 9    xi=(p0*yi-p1*y(1))/(yi-y(1))
c   Iterate.
      do 1 m=2,20
      if((xi-ps0)*(xi-ps1).le.0.) go to 10
      xi=.5*(ps0+ps1)
c   Save the current best guess of the zero.
 10   y(m)=yi
      x(m)=xi
c   Start iteration at the current best guess of the zero.
      yi=xmod(xi,nph,1,nr,r0)-xg
c   Check for convergence.
      if(abs(yi).le.xtol) go to 6
      if(yi.gt.0.) go to 11
      ps0=xi
      xs0=yi
      go to 2
 11   ps1=xi
      xs1=yi
 2    do 1 j=1,m
 1    xi=(x(j)*yi-xi*y(j))/(yi-y(j))
      write(msg,100)n,nph,xg
 100  format('Iteration did not converge:  n, nph, xg -',i4,i2,f7.4)
c     call abort(msg)
c     Following two lines are changed because intrinsic function
c     abort actually does not take any arguments and the gnu f77
c     compiler complains bitterly about it!
      write(*,*)msg
      call abort
c   Return the final best guess of the zero.
 6    pp(n,nph)=xi
      xx(n,nph)=yi+xg
      rr(n,nph)=r0
      return
      end
c
      subroutine zgrid(nph,m0)
c         depth grid
      save 
      include 'limits.inc'
      double precision pb,pm,zm
      common/zgridc/pb(nsl1),pm(ndp1,2),zm(ndp1,2),zm0(ndp1,2),
     1 ndex(ndp1,2),n
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
      data tol/1e-6/
c
      n1=n+1
      do 1 i=n1,1,-1
      if(abs(u(1,nph)-sngl(pb(i))).le.tol) go to 2
 1    continue
 2    n1=i
      write(10,*)'zgrid',n1,u(1,nph),sngl(pb(n1))
      l=1
      i=n1-1
      j=1
      pm(1,nph)=pb(n1)
      zm(1,nph)=z(1)
      zm0(1,nph)=1.
      ndex(1,nph)=n1
 9    if(u(j+1,nph).le.pb(i)) go to 10
      if(u(j+1,nph).gt.u(j,nph)) go to 11
      j=j+1
      go to 9
 11   i=i+2
 12   if(u(j+1,nph)-pb(i)) 15,16,10
 15   j=j+1
      go to 12
 16   j=j+1
 10   l=l+1
      pm(l,nph)=pb(i)
      zm(l,nph)=findep(sngl(pb(i)),j,nph)
      zm0(l,nph)=dexp(zm(l,nph))
      ndex(l,nph)=i
      if(i.le.2) go to 14
      i=i-1
      go to 9
 14   m=l
      m1=m+1
      pm(m1,nph)=pb(1)
      zm(m1,nph)=-1d6
      zm0(m1,nph)=0.
      ndex(m1,nph)=1
      m0=m
      return
      end
c
      function findep(u0,k,nph)
c
c $$$$$ calls emdlv $$$$$
c
c   Function findep returns the equivalent depth for which the model
c   slowness is u0.  If nph = 1, P slownesses are searched.  If nph = 2,
c   S slownesses are searched.  K is taken as the index into equivalent
c   depth near where the desired slowness should be found.  Function
c   evaluations are made until u0 is fit to a relative precision of aep
c   by a process of inverse iterative interpolation using Aitken's method.
c
      save
      include 'limits.inc'
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
      common/roughc/ric,roc
      common/findc/x(20),y(20),m
      data aep/1e-6/
c
      m=1
      a0=1./xn
      if(abs(z(k)-z(k+1)).le.-aep*z(k)) go to 14
      kph=nph
      x1=exp(z(k+1))
      x2=exp(z(k))
      if(a0*x1.lt.roc) kph=1
c   Set up iteration with first two trial points, x1 and x2.
      if(k+1.ge.nk) go to 3
      if(abs(z(k+1)-z(k+2)).gt.-aep*z(k+1)) go to 3
      call emdlv(a0*x1*(1.+aep),vp,vs)
      go to 4
 3    call emdlv(a0*x1,vp,vs)
 4    if(kph.eq.1) y(1)=pn*x1/vp-u0
      if(kph.eq.2) y(1)=pn*x1/vs-u0
      if(abs(y(1)).gt.aep*u0) go to 7
      findep=z(k+1)
      return
 7    x(1)=x1
      if(k.le.1) go to 5
      if(abs(z(k-1)-z(k)).gt.-aep*z(k)) go to 5
      call emdlv(a0*x2*(1.-aep),vp,vs)
      go to 6
 5    call emdlv(a0*x2,vp,vs)
 6    if(kph.eq.1) yi=pn*x2/vp-u0
      if(kph.eq.2) yi=pn*x2/vs-u0
      if(abs(yi).le.aep*u0.or.yi.eq.y(1)) go to 14
      xi=(x1*yi-x2*y(1))/(yi-y(1))
      if(abs(xi-x1).le.aep*amax1(abs(xi),aep).or.abs(xi-x2).le.
     1 aep*amax1(abs(xi),aep)) go to 13
c   Iterate.
      do 1 m=2,20
c   Save the current best guess of the zero.
      y(m)=yi
      x(m)=xi
c   Start iteration at the current best guess of the zero.
      call emdlv(a0*xi,vp,vs)
      if(kph.eq.1) yi=pn*xi/vp-u0
      if(kph.eq.2) yi=pn*xi/vs-u0
      do 2 j=1,m
      if(yi.eq.y(j)) go to 13
 2    xi=(x(j)*yi-xi*y(j))/(yi-y(j))
c   Check for convergence.
 1    if(abs(xi-x(m)).le.aep*amax1(abs(xi),aep)) go to 13
      m=0
c   Return the final best guess of the zero.
 13   findep=alog(amax1(xi,aep))
      return
 14   findep=z(k)
      return
      end
c
      subroutine brkpts(m,nph)
c
c $$$$$ calls efec, efe4, efe8, and phcod $$$$$
c
c        sets up discontinuity information
      save
      include 'limits.inc'
      character*8 code,tmp
      double precision pb,pm,zm
      common/zgridc/pb(nsl1),pm(ndp1,2),zm(ndp1,2),zm0(ndp1,2),
     1 ndex(ndp1,2),n
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      data dtol/1d-6/
c
      lb=1
      lbrk(1,nph)=ndex(1,nph)
      lc=0
      i=1
      isw=1
c   Search for discontinuities.
 17   if(dabs(zm(i,nph)-zm(i+1,nph)).le.-zm(i+1,nph)*dtol) go to 18
      if(i.ge.m) go to 22
      i=i+1
c   No discontinuity.
      go to (29,30),isw
c   Have we hit a high slowness zone?
 29   if(pm(i,nph).le.pm(i-1,nph)) go to 17
c   Yes, flag it.
      isw=2
c   If the high slowness zone is topped with a discontinuity, go back
c   to the usual processing.
      if(lbrk(lb,nph).eq.ndex(i-1,nph)) go to 17
c   Otherwise, mark the grazing ray to the top of the zone.
      lb=lb+1
      lbrk(lb,nph)=ndex(i-1,nph)
      call phcod(lc,nph,sngl(zm(i-1,nph)),-1)
      go to 17
c   We are already in a high slowness zone.  See if we have hit bottom.
 30   if(pm(i,nph).ge.pm(i-1,nph)) go to 17
c   Yes we have.  Reset the high slowness zone flag.
      isw=1
      go to 17
c   Discontinuity!  See what kind.
 18   if(pm(i+1,nph).gt.pm(i,nph)) go to 19
c   Velocity increase.  Flag the bottom of the step.
      lb=lb+1
      lbrk(lb,nph)=ndex(i,nph)
      call phcod(lc,nph,sngl(zm(i,nph)),1)
      if(lbrk(lb,nph).lt.lbrk(lb-1,nph)) go to 20
      lb=lb-1
      lc=lc-1
      code(lc,nph)=code(lc+1,nph)
c   Find the top of the discontinuity.
 20   if(i.ge.m) go to 22
      i=i+1
      if(dabs(zm(i,nph)-zm(i+1,nph)).le.-zm(i+1,nph)*dtol) go to 20
c   Flag the top of the step.
      lb=lb+1
      lbrk(lb,nph)=ndex(i,nph)
      if(lbrk(lb,nph).lt.lbrk(lb-1,nph)) go to 17
      lb=lb-1
      lc=lc-1
      go to 17
c   Velocity decrease. Flag the top of the step.
 19   lb=lb+1
      lbrk(lb,nph)=ndex(i,nph)
      call phcod(lc,nph,sngl(zm(i,nph)),-1)
c   Find the bottom of the discontinuity.
 21   if(i.ge.m) go to 22
      i=i+1
      if(dabs(zm(i,nph)-zm(i+1,nph))+zm(i+1,nph)*dtol)21,21,17
c   We have hit the bottom of the model.
 22   call phcod(lc,nph,sngl(zm(m,nph)),-1)
      call efe4(lb,lbrk(1,nph))
      call efec(lc,code(1,nph),tmp)
      lbb(nph)=lb
      lcb(nph)=lc
      return
      end
c
      subroutine phcod(lc,nph,z0,kfl)
c         set up phase codes
      save 
      include 'limits.inc'
      character*8 code
      character*1 tag,suf
      character*2 pre
      character*4 dep
      character*10 buf
      common/brkptc/code(nbr1,2)
      common/xmodc/z(nmd0),u(nmd0,2),xn,tn,pn,nk
      common/roughc/ric,roc
      data ln/8/
c
      d0=(1.-exp(z0))/xn
      r0=1./xn-d0
      idep=d0+.5
      if(lc.gt.0) go to 1
      lc=1
      pre=' '
      suf=' '
      if(nph.le.1) tag='P'
      if(nph.ge.2) tag='S'
      code(lc,nph)='t'//tag//'g'
      go to 2
 1    lc=lc+1
      if(idep.gt.70) go to 3
      code(lc,nph)='t'//tag//'b'
      if(code(lc-1,nph)(3:3).eq.'b') code(lc-1,nph)(3:3)='g'
      if(code(lc-2,nph)(3:3).eq.'b') code(lc-2,nph)(3:3)='g'
      go to 2
 3    if(pre.eq.' ') code(lc,nph)='t'//tag
      if(pre.ne.' ') code(lc,nph)='t'//tag//pre//suf//tag
      if(code(lc-2,nph)(3:3).ne.'g'.and.code(lc-2,nph)(3:3).ne.'b')
     1  go to 9
      code(lc,nph)(3:3)='n'
      code(lc-1,nph)='r'//tag//'m'//tag
 9    j=0
      do 4 i=1,ln
      if(code(lc,nph)(i:i).eq.' ') go to 4
      j=j+1
      code(lc,nph)(j:j)=code(lc,nph)(i:i)
 4    continue
      if(j.lt.ln) code(lc,nph)(j+1:ln)=' '
      if(code(lc,nph)(2:).eq.'PKP') code(lc,nph)(2:)='PKPab'
      if(code(lc,nph)(2:).eq.'PKIKP') code(lc,nph)(2:)='PKPdf'
      if(code(lc,nph)(2:).eq.'SKS') code(lc,nph)(2:)='SKSab'
      if(code(lc,nph)(2:).eq.'SKIKS') code(lc,nph)(2:)='SKSdf'
      if(abs(r0-roc).gt.20.) go to 5
      pre='K'
      if(kfl.le.0) return
      lc=lc+1
      code(lc,nph)='r'//tag//'c'//tag
      return
 5    if(abs(r0-ric).gt.20.) go to 2
      pre='KI'
      suf='K'
      if(kfl.le.0) return
      lc=lc+1
      code(lc,nph)='r'//tag//'KiK'//tag
      return
 2    if(kfl.le.0) return
      lc=lc+1
      write(dep,100)idep
 100  format(i4)
      buf=tag//pre//'d'//dep//suf//tag
      code(lc,nph)='r'
      j=1
      do 8 i=1,10
      if(buf(i:i).eq.' '.or.j.ge.ln) go to 8
      j=j+1
      code(lc,nph)(j:j)=buf(i:i)
 8    continue
      if(j.lt.ln) code(lc,nph)(j+1:ln)=' '
      return
      end
c
      subroutine efe4(n,na)
c
c $$$$$ calls no other routine $$$$$
c
c   Integer array na(n) is transposed end-for-end.
c
      save
      dimension na(n)
c
      if(n.le.1) return
      n2=n/2
      j=n
      do 1 i=1,n2
      nb=na(i)
      na(i)=na(j)
      na(j)=nb
 1    j=j-1
      return
      end
c
      subroutine efe8(n,da)
c
c $$$$$ calls no other routine $$$$$
c
c   Double precision array da(n) is transposed end-for-end.
c
      save
      double precision da(n)
      double precision db
c
      if(n.le.1) return
      n2=n/2
      j=n
      do 2 i=1,n2
      db=da(i)
      da(i)=da(j)
      da(j)=db
 2    j=j-1
      return
      end
c
      subroutine efec(n,ia,ib)
c
c $$$$$ calls no other routine $$$$$
c
c   Character array ia(n) is transposed end-for-end.  Ib must be a
c   character variable with the same character width as each element of
c   array ia.
c
      save
      character*(*) ia(n),ib
c
      if(n.le.1) return
      n2=n/2
      j=n
      do 3 i=1,n2
      ib=ia(i)
      ia(i)=ia(j)
      ia(j)=ib
 3    j=j-1
      return
      end
