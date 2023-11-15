
      program dwm_michel

c code Michel DWM 3D layer-cake
c DK : portage en double precision Feb 99

c  COMPUTES A SURFACE SEISMIC PROFILE FOR A POINT SOURCE
c  IN A LAYERED MEDIUM.
c  THE CODE USES THE DISCRETE WAVENUMBER METHOD (Bouchon,
c  BSSA, 71, 959-971, 1981) AND THE METHOD OF REFLECTIVITY
c  AND TRANSMISSIVITY MATRICES OF Kennett (Muller, Journal
c  of Geophysics, 58, 153-174, 1985).

      implicit double precision(a-h,o-z)

c Parameters to be defined to dimension the arrays:
c - nlmax = maximum number of layers.
c - nrmax = maximum number of recivers.
c - ntmax = maximum number of points of each seismogram.
      parameter (nlmax=100)
      parameter (nrmax=1000)
      parameter (ntmax=1024)

      parameter (nlmx=nlmax+5)
      dimension th0(nlmx),alpha0(nlmx),beta0(nlmx),dens0(nlmx),
     $  qp0(nlmx),qs0(nlmx),th(nlmx),alp(nlmx),bet(nlmx),
     $  dens(nlmx),qp(nlmx),qs(nlmx),
     $  r(nrmax),
     $  yy(ntmax),y(ntmax+ntmax+3),sy(ntmax,nrmax)

      complex*16 alpha(nlmx),beta(nlmx),aqp(nlmx),aqs(nlmx),emu(nlmx),
     $  emu2(nlmx),elam(nlmx),wa2(nlmx),wb2(nlmx),cth(nlmx),
     $  as(2),bs(2),su(2),sd(2),
     $  al(2),bl(2),
     $  wza(nlmx),wzb(nlmx),
     $  woa(nlmx),wob(nlmx),
     $  rd(2,2,nlmx),ru(2,2,nlmx),td(2,2,nlmx),tu(2,2,nlmx),
     $  mb(2,2,nlmx),mt(2,2,nlmx),nb(2,2,nlmx),nt(2,2,nlmx),
     $  mtt(2,2,nlmx),ntt(2,2,nlmx),
     $  quu(2,2,nlmx),qud(2,2,nlmx),qdu(2,2,nlmx),qdd(2,2,nlmx),
     $  tup(2,2,nlmx,nlmx),tdw(2,2,nlmx,nlmx),
     $  f(2,2,nlmx),g(2,2,nlmx),ee(2,2,nlmx),c(2,2),d(2,2),e(2,2),
     $  source(ntmax/2),u(nrmax,ntmax/2)

      complex*16 ai,omega,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,
     $  c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c35,
     $  c36,c0e,c0f,
     $  cu,cu2,cc,
     $  d1d,d2d,d1u,d2u,dd,xlnf

      complex*16 q1,q2,uu1,uu2

      real*8 arg,bj0(nrmax),bj1(nrmax),by0,by1,a20

      logical down,surf

      ai=(0.d0,1.d0)
      pi=3.141592653589793d0
      pi2=pi+pi
      eps=1.d-4

c ** INPUT DATA:

c input file:
      open(15,file="reyk.da",form="formatted")

c number of layers (including the half-space):
      read(15,*) nl0

c thickness, P-wave velocity, S-wave velocity, density, QP, QS
c for each layer (the value of thickness for the half-space does
c not matter):
      do 55 l=1,nl0
   55 read(15,*) th0(l),alpha0(l),beta0(l),dens0(l),qp0(l),qs0(l)

c source depth:
      read(15,*) z0

c number of receivers, horizontal offsets from the source of the
c closest (r0) and furthest (r1) receivers:
      read(15,*) nr,r0,r1

c type of source: explosion (isource=1) or vertical force (isource=2)
c type of receiver: vertical (icomp=1), horizontal (icomp=2), or
c pressure (icomp=3) sensors:
      read(15,*) isource,icomp

c number of points in time of each ouput trace, length of the desired
c time window for each trace, central period of the Ricker source time
c function, starting time of each trace (t0=0 implies that the traces
c start at the origin time of the source)
c (p.s. if one desires another source time dependence than a Ricker,
c  one only has to change the expression for "source(if)" ).
      read(15,*) ntime,tl,tricker,t0

c periodicity length (xl should be such that arrivals from the next
c closest source arrive at the receivers after time t0+tl). This
c periodicity length is a feature of the discrete wavenumber method.
c m is a large integer number which represents the maximum allowed
c of wavenumber terms for the wavenumber sums. The actual number of
c wavenumbers considered is determined by a convergence criterion
c for each frequency and is written in the output. m should be larger
c than any of these numbers for the results to be accurate.
      read(15,*) xl,m

c down=.true. if only the downgoing wavefield emitted by the source
c is to be considered.
c down=.false. if the total source wavefield is desired.
      down=.false.

c surf=.false. if the free surface effect is to be removed (that is
c surface reflections are not calculated, no surface waves, no
c ground roll,... ). The surface is treated as a transparent (not a
c reflecting) boundary.
c surf=.true. implies that all free surface effects are calculated.
      surf=.true.

      do 65 l=1,nl0
   65 write(6,*) 'layer ',l,':',
     .     th0(l),alpha0(l),beta0(l),dens0(l),qp0(l),qs0(l)
      write(6,*) z0
      write(6,*) nr,r0,r1
      write(6,*) isource,icomp
      write(6,*) ntime,tl,tricker,t0
      write(6,*) xl,m

      pil=pi2/xl
      dt=tl/dble(ntime)
      nfreq=ntime/2
      dfreq=1.d0/tl
      q=tl/2.d0
      aw=-pi/q

c r(ir) is the horizontal offset of receiver ir:
      dr=0.d0
      if(nr.ne.1) dr=(r1-r0)/dble(nr-1)
      do 1 ir=1,nr
    1 r(ir)=r0+dr*dble(ir-1)

      lr=2
      zl=0.d0
      do 2 l=1,nl0
      zl=zl+th0(l)
      if(zl.ge.z0) go to 3
    2 continue
    3 if(z0.ne.0.) go to 33
      ls=l+2
      nl=nl0+3
      go to 34
   33 ls=l+3
      nl=nl0+4
   34 nl1=nl-1
      th(1)=0.d0
      alp(1)=340.d0
      bet(1)=1.d0
      dens(1)=1.3d-3
      qp(1)=1.d6
      qs(1)=1.d6
      do 4 l=2,nl
      if(l.eq.2) ll=1
      if(l.gt.2.and.l.lt.ls) ll=l-2
      if(l.eq.3.and.l.eq.ls) ll=l-2
      if(l.gt.3.and.l.eq.ls) ll=l-3
      if(l.gt.ls.and.ls.eq.3) ll=l-3
      if(l.gt.ls.and.ls.gt.3) ll=l-4
      th(l)=th0(ll)
      if(l.eq.lr) th(l)=0.d0
      if(l.eq.(ls-1).and.ls.ne.3) th(l)=z0-(zl-th0(ll))
      if(l.eq.ls) th(l)=0.d0
      if(l.eq.(ls+1).and.ls.eq.3) th(l)=th0(ll)
      if(l.eq.(ls+1).and.ls.ne.3) th(l)=zl-z0
      alp(l)=alpha0(ll)
      bet(l)=beta0(ll)
      dens(l)=dens0(ll)
      qp(l)=qp0(ll)
    4 qs(l)=qs0(ll)
      th(nl)=0.d0

      do 8 l=1,nl
    8 cth(l)=-ai*th(l)

      write(6,*) lr,ls
      do 1010 l=1,nl
 1010 write(6,*) 'structure layer ',l,':',
     .      th(l),alp(l),bet(l),dens(l),qp(l),qs(l)

c The introduction of anelastic attenuation produces some
c dispersion of the velocities, althought very light. freq0
c is the frequency of reference for the velocities.
c For more... see the discussion in Aki and Richards on Q.
      freq0=20.d0
      do 5 l=1,nl
      aqp(l)=(1.d0+ai/(qp(l)+qp(l)))/(1.d0+.25d0/qp(l)**2)*alp(l)
    5 aqs(l)=(1.d0+ai/(qs(l)+qs(l)))/(1.d0+.25d0/qs(l)**2)*bet(l)

      freq=0.d0

c** Begin loop on frequencies
      do 100 if=1,nfreq
      rw=pi2*freq
      omega=dcmplx(rw,aw)

c  source is the Fourier transform of a Ricker:
      source(if)=-omega**2*cdexp(-(omega/pi2*tricker)**2)
     $  *cdexp(ai*omega*t0)

      zom=dsqrt(rw**2+aw**2)/pi2
      if(if.eq.1) phi=-pi/2.d0
      if(if.ne.1) phi=datan(aw/rw)
      xlnf=(ai*phi+dlog(zom)-dlog(freq0))/pi
      do 10 l=1,nl
      alpha(l)=aqp(l)/(1.d0-xlnf/qp(l))
      if(freq.eq.0.d0) alpha(l)=alp(l)
      beta(l)=aqs(l)/(1.d0-xlnf/qs(l))
      if(freq.eq.0.d0) beta(l)=bet(l)
      emu(l)=dens(l)*beta(l)**2
      emu2(l)=emu(l)+emu(l)
      elam(l)=dens(l)*alpha(l)**2-emu2(l)

      wa2(l)=(omega/alpha(l))**2
   10 wb2(l)=(omega/beta(l))**2

      c0e=pi2/xl
      c0f=1.d0/(2.d0*xl*dens(ls)*omega**2)

      do 899 ir=1,nr
  899 u(if,ir)=0.d0
      uu1=0.d0
      uu2=0.d0

      do 50 k=1,m
      ak=pil*dble(k-1)
      ak2=ak**2

      do 20 ir=1,nr
      a1=ak*r(ir)
      arg=abs(a1)
      call ff01ad(bj0(ir),by0,arg,0)
      call ff02ad(bj1(ir),by1,arg,0)
      if(a1.lt.0.) bj1(ir)=-bj1(ir)
   20 continue

      do 30 l=1,nl
      c1=wa2(l)-ak2
      wza(l)=cdsqrt(c1)
      q1=wza(l)
      if(dimag(q1).gt.0.) wza(l)=-wza(l)
      c2=wb2(l)-ak2
      wzb(l)=cdsqrt(c2)
      q2=wzb(l)
      if(dimag(q2).gt.0.) wzb(l)=-wzb(l)
      woa(l)=wza(l)/omega
      wob(l)=wzb(l)/omega
   30 continue

c --------------------------------------------------------------------------
c
c Calcul des matrices de reflexion et de transmission
c
c --------------------------------------------------------------------------
      cu=ak/omega
      cu2=cu**2
      do 301 l=1,nl1
      l1=l+1
      cc=emu2(l)-emu2(l1)
      c1=cc*cu2
      c2=c1-dens(l)
      c3=c1+dens(l1)
      c4=c1-dens(l)+dens(l1)
      c5=c2*c2
      c6=c3*c3
      c7=c4*c4*cu2
      a1=dens(l)*dens(l1)
      c8=woa(l)*wob(l)
      c9=woa(l)*wob(l1)
      c10=woa(l1)*wob(l)
      c11=woa(l1)*wob(l1)
      c14=a1*c9
      c15=a1*c10
      c16=cc*c1*c8*c11
      c17=c5*c11
      c18=c6*c8
      d1d=c7+c17+c15
      d2d=c16+c18+c14
      d1u=c7+c18+c14
      d2u=c16+c17+c15
      c19=c3*wob(l)-c2*wob(l1)
      c20=c3*woa(l)-c2*woa(l1)
      dd=d1d+d2d
      rd(1,1,l1)=(d2d-d1d)/dd
      ru(1,1,l1)=(d2u-d1u)/dd
      c21=(cu+cu)*woa(l)
      c22=(cu+cu)*wob(l)
      c23=(cu+cu)*woa(l1)
      c24=(cu+cu)*wob(l1)
      c25=(c4*c3+cc*c2*c11)/dd
      rd(2,1,l1)=-c21*c25
      c35=(c4*c2+cc*c3*c8)/dd
      ru(2,1,l1)=c23*c35
      c26=dens(l)/dd
      td(1,1,l1)=(c26+c26)*woa(l)*c19
      td(2,1,l1)=-c26*c21*(c4+cc*c10)
      c27=(a1+a1)*(c10-c9)
      rd(2,2,l1)=(d2d-d1d+c27)/dd
      rd(1,2,l1)=c22*c25
      td(2,2,l1)=(c26+c26)*wob(l)*c20
      td(1,2,l1)=c26*c22*(c4+cc*c9)
      c36=dens(l1)/dd
      tu(1,1,l1)=(c36+c36)*woa(l1)*c19
      tu(2,1,l1)=-c36*c23*(c4+cc*c9)
      ru(2,2,l1)=(d2u-d1u-c27)/dd
      ru(1,2,l1)=-c24*c35
      tu(2,2,l1)=(c36+c36)*wob(l1)*c20
      tu(1,2,l1)=c36*c24*(c4+cc*c10)
  301 continue

      if(surf) go to 500
      do 499 i=1,2
      do 499 j=1,2
      ru(i,j,2)=0.d0
  499 tu(i,j,2)=0.d0
      tu(1,1,2)=1.d0
      tu(2,2,2)=1.d0
  500 continue

      do 304 i=1,2
      do 304 j=1,2
      mt(i,j,nl)=0.d0
      mb(i,j,nl1)=rd(i,j,nl)
      nb(i,j,1)=0.d0
      nt(i,j,2)=ru(i,j,2)
      g(i,j,1)=tu(i,j,2)
  304 f(i,j,nl)=td(i,j,nl)
      do 303 l=2,nl1
      c1=cth(l)*wza(l)
      c2=cth(l)*wzb(l)
      ee(1,1,l)=cdexp(c1)
      ee(2,2,l)=cdexp(c2)
      ee(1,2,l)=0.d0
  303 ee(2,1,l)=0.d0
      do 306 l=nl1,2,-1
      l1=l-1
      c1=ee(1,1,l)
      c2=ee(2,2,l)
      c3=c1*c2
      mt(1,1,l)=mb(1,1,l)*c1*c1
      mt(1,2,l)=mb(1,2,l)*c3
      mt(2,1,l)=mb(2,1,l)*c3
      mt(2,2,l)=mb(2,2,l)*c2*c2
      do 308 i=1,2
      do 308 j=1,2
      c(i,j)=0.d0
      do 308 ij=1,2
  308 c(i,j)=c(i,j)+mt(i,ij,l)*ru(ij,j,l)
      e(1,1)=1.d0-c(1,1)
      e(2,2)=1.d0-c(2,2)
      e(1,2)=-c(1,2)
      e(2,1)=-c(2,1)
      call inv2(e)
      do 310 i=1,2
      do 310 j=1,2
      c(i,j)=0.d0
      do 310 ij=1,2
  310 c(i,j)=c(i,j)+tu(i,ij,l)*e(ij,j)
      do 312 i=1,2
      do 312 j=1,2
      e(i,j)=0.d0
      do 312 ij=1,2
  312 e(i,j)=e(i,j)+c(i,ij)*mt(ij,j,l)
      do 314 i=1,2
      do 314 j=1,2
      mb(i,j,l1)=rd(i,j,l)
      do 314 ij=1,2
  314 mb(i,j,l1)=mb(i,j,l1)+e(i,ij)*td(ij,j,l)
  306 continue
      do 316 l=2,nl1
      l1=l+1
      c1=ee(1,1,l)
      c2=ee(2,2,l)
      c3=c1*c2
      nb(1,1,l)=nt(1,1,l)*c1*c1
      nb(1,2,l)=nt(1,2,l)*c3
      nb(2,1,l)=nt(2,1,l)*c3
      nb(2,2,l)=nt(2,2,l)*c2*c2
      do 318 i=1,2
      do 318 j=1,2
      c(i,j)=0.d0
      do 318 ij=1,2
  318 c(i,j)=c(i,j)+nb(i,ij,l)*rd(ij,j,l1)
      e(1,1)=1.d0-c(1,1)
      e(2,2)=1.d0-c(2,2)
      e(1,2)=-c(1,2)
      e(2,1)=-c(2,1)
      call inv2(e)
      do 320 i=1,2
      do 320 j=1,2
      c(i,j)=0.d0
      do 320 ij=1,2
  320 c(i,j)=c(i,j)+td(i,ij,l1)*e(ij,j)
      do 322 i=1,2
      do 322 j=1,2
      e(i,j)=0.d0
      do 322 ij=1,2
  322 e(i,j)=e(i,j)+c(i,ij)*nb(ij,j,l)
      do 324 i=1,2
      do 324 j=1,2
      nt(i,j,l1)=ru(i,j,l1)
      do 324 ij=1,2
  324 nt(i,j,l1)=nt(i,j,l1)+e(i,ij)*tu(ij,j,l1)
  316 continue
      do 350 l=2,nl1
      do 3061 i=1,2
      ntt(i,1,l)=nt(i,1,l)*ee(1,1,l)
      ntt(i,2,l)=nt(i,2,l)*ee(2,2,l)
      mtt(i,1,l)=mb(i,1,l)*ee(1,1,l)
 3061 mtt(i,2,l)=mb(i,2,l)*ee(2,2,l)
      do 451 i=1,2
      do 451 j=1,2
      c(i,j)=0.d0
      do 451 ij=1,2
  451 c(i,j)=c(i,j)+mt(i,ij,l)*nt(ij,j,l)
      e(1,1)=1.d0-c(1,1)
      e(1,2)=-c(1,2)
      e(2,1)=-c(2,1)
      e(2,2)=1.d0-c(2,2)
      call inv2(e)
      c(1,1)=ee(1,1,l)*mb(1,1,l)
      c(1,2)=ee(1,1,l)*mb(1,2,l)
      c(2,1)=ee(2,2,l)*mb(2,1,l)
      c(2,2)=ee(2,2,l)*mb(2,2,l)
      do 452 i=1,2
      do 452 j=1,2
      quu(i,j,l)=e(i,j)
      qud(i,j,l)=0.d0
      do 452 ij=1,2
  452 qud(i,j,l)=qud(i,j,l)+e(i,ij)*c(ij,j)
      do 461 i=1,2
      do 461 j=1,2
      c(i,j)=0.d0
      do 461 ij=1,2
  461 c(i,j)=c(i,j)+nb(i,ij,l)*mb(ij,j,l)
      e(1,1)=1.d0-c(1,1)
      e(1,2)=-c(1,2)
      e(2,1)=-c(2,1)
      e(2,2)=1.d0-c(2,2)
      call inv2(e)
      c(1,1)=ee(1,1,l)*nt(1,1,l)
      c(1,2)=ee(1,1,l)*nt(1,2,l)
      c(2,1)=ee(2,2,l)*nt(2,1,l)
      c(2,2)=ee(2,2,l)*nt(2,2,l)
      do 462 i=1,2
      do 462 j=1,2
      qdu(i,j,l)=0.d0
      qdd(i,j,l)=e(i,j)
      do 462 ij=1,2
  462 qdu(i,j,l)=qdu(i,j,l)+e(i,ij)*c(ij,j)
  350 continue
      do 3062 i=1,2
      do 3062 j=1,2
      ntt(i,j,1)=0.d0
 3062 mtt(i,j,nl)=0.d0
      do 370 i=1,2
      do 370 j=1,2
      quu(i,j,1)=0.d0
      quu(i,j,nl)=0.d0
      qud(i,j,1)=mb(i,j,1)
      qud(i,j,nl)=0.d0
      qdd(i,j,1)=0.d0
      qdd(i,j,nl)=0.d0
      qdu(i,j,1)=0.d0
  370 qdu(i,j,nl)=nt(i,j,nl)
      do 371 i=1,2
      qdd(i,i,1)=1.d0
  371 quu(i,i,nl)=1.d0
      do 380 l=1,nl1
      l1=l+1
      do 381 i=1,2
      do 381 j=1,2
      c(i,j)=0.d0
      d(i,j)=0.d0
      do 381 ij=1,2
      c(i,j)=c(i,j)-rd(i,ij,l1)*nb(ij,j,l)
  381 d(i,j)=d(i,j)-ru(i,ij,l1)*mt(ij,j,l1)
      c(1,1)=1.d0+c(1,1)
      c(2,2)=1.d0+c(2,2)
      d(1,1)=1.d0+d(1,1)
      d(2,2)=1.d0+d(2,2)
      call inv2(c)
      call inv2(d)
      do 382 i=1,2
      do 382 j=1,2
      g(i,j,l)=0.d0
      f(i,j,l)=0.d0
      do 382 ij=1,2
      g(i,j,l)=g(i,j,l)+c(i,ij)*tu(ij,j,l1)
  382 f(i,j,l)=f(i,j,l)+d(i,ij)*td(ij,j,l1)
  380 continue
      l=ls
      if(lr.gt.ls) go to 555
      l1=l-1
      lu=lr
      do 5100 i=1,2
      do 5100 j=1,2
 5100 tup(i,j,l,lu)=0.d0
      tup(1,1,l,lu)=1.0d0
      tup(2,2,l,lu)=1.0d0
      do 512 ll=lu,l1
      if(ll.eq.lu) go to 514
      tup(1,1,l,lu)=tup(1,1,l,lu)*ee(1,1,ll)
      tup(1,2,l,lu)=tup(1,2,l,lu)*ee(2,2,ll)
      tup(2,1,l,lu)=tup(2,1,l,lu)*ee(1,1,ll)
      tup(2,2,l,lu)=tup(2,2,l,lu)*ee(2,2,ll)
  514 do 511 i=1,2
      do 511 j=1,2
      d(i,j)=0.d0
      do 511 ij=1,2
  511 d(i,j)=d(i,j)+tup(i,ij,l,lu)*g(ij,j,ll)
      do 515 i=1,2
      do 515 j=1,2
  515 tup(i,j,l,lu)=d(i,j)
  512 continue
  510 continue
      go to 556
  555 l1=l+1
      ld=lr
      do 5200 i=1,2
      do 5200 j=1,2
 5200 tdw(i,j,l,ld)=0.d0
      tdw(1,1,l,ld)=1.0d0
      tdw(2,2,l,ld)=1.0d0
      do 522 ll=l,ld-1
      if(ll.eq.l) go to 524
      tdw(1,1,l,ld)=ee(1,1,ll)*tdw(1,1,l,ld)
      tdw(1,2,l,ld)=ee(1,1,ll)*tdw(1,2,l,ld)
      tdw(2,1,l,ld)=ee(2,2,ll)*tdw(2,1,l,ld)
      tdw(2,2,l,ld)=ee(2,2,ll)*tdw(2,2,l,ld)
  524 do 521 i=1,2
      do 521 j=1,2
      d(i,j)=0.d0
      do 521 ij=1,2
  521 d(i,j)=d(i,j)+f(i,ij,ll)*tdw(ij,j,l,ld)
      do 525 i=1,2
      do 525 j=1,2
  525 tdw(i,j,l,ld)=d(i,j)
  522 continue
  520 continue
  556 continue
c ------------------------------------------------------------------------
c
c
c ------------------------------------------------------------------------

      if(isource.ne.1) go to 600
      sd(1)=c0e*ak/(ai*wza(ls))
      sd(2)=0.d0
      su(1)=sd(1)
      su(2)=0.d0
      go to 601
  600 if(isource.ne.2) go to 601
      sd(1)=c0f*ak
      sd(2)=c0f*ak2/wzb(ls)
      su(1)=-sd(1)
      su(2)=sd(2)
  601 continue

      if(down) su(1)=0.d0
      if(down) su(2)=0.d0

      l=ls
      do 801 ip=1,2
      bs(ip)=0.d0
      as(ip)=0.d0
      do 801 jp=1,2
      bs(ip)=bs(ip)+quu(ip,jp,l)*su(jp)+qud(ip,jp,l)*
     $  sd(jp)
  801 as(ip)=as(ip)+qdu(ip,jp,l)*su(jp)+qdd(ip,jp,l)*
     $  sd(jp)
      lj=ls
      li=lr
      if(lr.lt.ls) go to 557
      do 802 ip=1,2
      al(ip)=0.d0
      do 802 jp=1,2
  802 al(ip)=al(ip)+tdw(ip,jp,lj,li)*as(jp)
      do 804 ip=1,2
      bl(ip)=0.d0
      do 804 jp=1,2
  804 bl(ip)=bl(ip)+mtt(ip,jp,li)*al(jp)
      go to 558
  557 do 806 ip=1,2
      bl(ip)=0.d0
      do 806 jp=1,2
  806 bl(ip)=bl(ip)+tup(ip,jp,lj,li)*bs(jp)
      do 808 ip=1,2
      al(ip)=0.d0
      do 808 jp=1,2
  808 al(ip)=al(ip)+ntt(ip,jp,li)*bl(jp)
  558 continue

      c1=ai*wza(lr)*(bl(1)-al(1))-ai*ak*(bl(2)+al(2))
      c2=-ak*(bl(1)+al(1))-wzb(lr)*(bl(2)-al(2))
      c3=-wa2(lr)*(bl(1)+al(1))

      q1=c1
      q2=c2
      uu1=uu1+c1
      uu2=uu2+c2

      if(icomp.eq.2) c1=c2
      if(icomp.eq.3) c1=c3
      do 1000 ir=1,nr
      if(icomp.eq.1.or.icomp.eq.3) a20=bj0(ir)
      if(icomp.eq.2) a20=bj1(ir)
 1000 u(if,ir)=u(if,ir)+c1*a20

      if(k.eq.1) go to 50
      a11=cdabs(q1)
      a12=cdabs(q2)
      a21=cdabs(uu1)*eps
      a22=cdabs(uu2)*eps
      if(a11.le.a21.and.a12.le.a22) go to 21
   50 continue

   21 continue

c writes the current frequency index and the number of wavenumbers
c considered in the wavenumber series:
      write(6,*) if,k

  100 freq=freq+dfreq
c ** End of loop on frequencies


c ** Calculation of the time domain solution:

      n3=ntime+ntime+3
      tex1=-aw*dt
      tex1=dexp(tex1)
      ex7=1.0d0
      do 140 i=1,ntime
      yy(i)=ex7
  140 ex7=ex7*tex1
      do 410 ir=1,nr
      do 517 i=1,n3
  517 y(i)=0.d0
      do 518 i=1,nfreq
c Convolution of the impulse response with the source time
c function (=spectral multiplication):
      q1=u(i,ir)*source(i)
      i2=i+i
      i1=i2-1
      y(i1)=dreal(q1)
      y(i2)=dimag(q1)
      y(n3-i1)=-y(i2)
  518 y(n3-i2)=y(i1)
c FFT:
      call four1(y,ntime,1)
      do 411 i=1,ntime
c yy(i) corrects the FFT result for the time attenuation which results
c from the imaginary part of the frequency:
      sy(i,ir)=y(i+i-1)*yy(i)
  411 continue
  410 continue

c output:
c - ir  :index of receiver considered
c - sy  :seismogram at receiver ir

c ** Plot the seismograms

c sy(i,ir) is the seismogram at receiver ir.
c ...

c ecriture au format text deplacement

      print *,'Saving seismograms in text format'

      open(unit=11,file='U_file.txt',status='unknown')
      do ir=1,nr
c      write(11,*) 'Receiver ',ir
      do it=1,ntime
            write(11,*) sngl(sy(it,ir))
      enddo
      enddo
      close(11)

      stop
      end


      subroutine inv2(a)

c inverser matrice 2x2

      implicit double precision(a-h,o-z)

      complex*16 a(2,2),b(2,2),det
      do 1 i=1,2
      do 1 j=1,2
    1 b(i,j)=a(i,j)
      det=b(1,1)*b(2,2)-b(1,2)*b(2,1)
      a(1,1)=b(2,2)/det
      a(2,2)=b(1,1)/det
      a(1,2)=-b(1,2)/det
      a(2,1)=-b(2,1)/det
      return
      end


c BESSEL FUNCTIONS:

C/     ADD NAME=FF01AD          HSL     F77     DOUBLE
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS FF01AD
      SUBROUTINE FF01AD(VJ0,VY0,XD,N)
C  STANDARD FORTRAN 66(A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION VJ0,VY0,X,Y,Z,Q1,Q2,Q3,FX,X1,X2,X3,
     1                 X4,XD,XLG,A,B,C,D,E
      DIMENSION A(73),B(18),C(19),D(18),E(18)
      EQUIVALENCE (A(1),B(1)),(A(19),C(1)),(A(38),D(1)),(A(56),E(1))
      DATA XLG /1.0D+70/
      DATA B(1),B(2),B(3),B(4),B(5),B(6),B(7),B(8),B(9),B(10),B(11),
     1     B(12),B(13),B(14),B(15),B(16),B(17),B(18)    /
     1   -.17D-18                  , .1222D-16             ,
     2   -.75885D-15               , .4125321D-13          ,
     3   -.194383469D-11           , .7848696314D-10       ,
     4   -.267925353056D-8         , .7608163592419D-7     ,
     5   -.176194690776215D-5      , .3246032882100508D-4  ,
     6   -.46062616620627505D-3    , .48191800694676045D-2 ,
     7   -.34893769411408885D-1    , .15806710233209726D0  ,
     8   -.37009499387264978D0     , .26517861320333681D0  ,
     9   -.87234423528522213D-2    , .31545594294978024D0  /
      DATA C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),
     1     C(12),C(13),C(14),C(15),C(16),C(17),C(18),C(19)    /
     A   -.1D-19                   , .39D-18               ,
     B   -.2698D-16                , .164349D-14           ,
     C   -.8747341D-13             , .402633082D-11        ,
     D   -.15837552542D-9          , .524879478733D-8      ,
     E   -.14407233274019D-6       , .32065325376548D-5    ,
     F   -.5632079141056987D-4     , .75311359325777423D-3 ,
     G   -.72879624795520792D-2    , .47196689595763387D-1 ,
     H   -.17730201278114358D0     , .26156734625504664D0  ,
     I    .17903431407718266D0     ,-.27447430552974527D0  ,
     J   -.66292226406569883D-1     /
      DATA D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),D(10),D(11),
     1     D(12),D(13),D(14),D(15),D(16),D(17),D(18)    /
     K   -.1D-19                   , .2D-19                ,
     L   -.11D-18                  , .55D-18               ,
     M   -.288D-17                 , .1631D-16             ,
     N   -.10012D-15               , .67481D-15            ,
     O   -.506903D-14              , .4326596D-13          ,
     O   -.43045789D-12            , .516826239D-11        ,
     P   -.7864091377D-10          , .163064646352D-8      ,
     Q   -.5170594537606D-7        , .307518478751947D-5   ,
     R   -.53652204681321174D-3    , .19989206986950373D1 /
      DATA E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10),E(11),
     1     E(12),E(13),E(14),E(15),E(16),E(17),E(18)   /
     S    .1D-19                   ,-.3D-19                ,
     T    .13D-18                  ,-.62D-18               ,
     U    .311D-17                 ,-.1669D-16             ,
     V    .9662D-16                ,-.60999D-15            ,
     W    .425523D-14              ,-.3336328D-13          ,
     X    .30061451D-12            ,-.320674742D-11        ,
     Y    .4220121905D-10          ,-.72719159369D-9       ,
     Z    .1797245724797D-7        ,-.74144984110606D-6    ,
     1    .683851994261165D-4      ,-.31111709210674018D-1 /
      X=XD
      Y=DABS(X)
      Z=Y*.125D0
      IF(Z .LE.1.0D0)GO TO 10
      Z=1.0D0/Z
      X2=4.0D0*Z*Z-2.0D0
      N1=38
      N2=55
      GO TO 70
   10 IF(Z .EQ. 0.0D0)GO TO  78
      X2=4.0D0*Z*Z-2.0D0
      N1=1
      N2=18
   70 DO 80 J=1,2
      Q3=0.0D0
      Q2=0.0D0
      DO 40  I=N1,N2
      Q1=Q2
      Q2=Q3
   40 Q3=X2*Q2-Q1+A(I)
      FX=(Q3-Q1)*.5D0
      IF(N1-19)50,51,52
   50 VJ0=FX
      IF(N .LE. 0)GO TO 75
      N1=19
      N2=37
      GO TO 80
   52 IF(N1.EQ.56)GO TO 53
      X1=FX
      N1=56
      N2=73
   80 CONTINUE
   78 VJ0=1.0D0
      VY0=-XLG
      GO TO 75
   51 VY0=.6366197723675813D0*DLOG(Y)*VJ0+FX
      GO TO 75
   53 X2=DCOS(Y-.7853981633974483D0)
      X3=DSIN(Y-.7853981633974483D0)
      X4=.7978845608028654D0/DSQRT(Y)
      FX=FX*Z
      VJ0=X4*(X1*X2-FX*X3)
      VY0=X4*(FX*X2+X1*X3)
   75 RETURN
      END


C/     ADD NAME=FF02AD          HSL     F77     DOUBLE
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS FF02AD
      SUBROUTINE FF02AD(VJ1,VY1,XD,N)
C  STANDARD FORTRAN 66(A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION VJ1,VY1,X,Y,Z,Q1,Q2,Q3,FX,X1,X2,X3,X4,
     1                 XD,XLG,A,B,C,D,E
      DIMENSION A(72),B(18),C(18),D(18),E(18)
      EQUIVALENCE (A(1),B(1)),(A(19),C(1)),(A(37),D(1)),(A(55),E(1))
      DATA XLG/1.0D+70 /
      DATA B(1),B(2),B(3),B(4),B(5),B(6),B(7),B(8),B(9),B(10),B(11),
     1     B(12),B(13),B(14),B(15),B(16),B(17),B(18)   /
     1   -.4D-19                   , .295D-17              ,
     2   -.19554D-15               , .1138572D-13          ,
     3   -.57774042D-12            , .2528123664D-10       ,
     4   -.94242129816D-9          , .2949707007278D-7     ,
     5   -.76175878054003D-6       , .1588701923993213D-4  ,
     6   -.26044438934858068D-3    , .32402701826838575D-2 ,
     7   -.29175524806154208D-1    , .17770911723972828D0  ,
     8   -.66144393413454325D0     , .12879940988576776D1  ,
     9   -.11918011605412169D1     , .12967175412105298D1  /
      DATA C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),
     1     C(12),C(13),C(14),C(15),C(16),C(17),C(18)   /
     A    .9D-19                   ,-.658D-17              ,
     B    .42773D-15               ,-.2440949D-13          ,
     C    .121143321D-11           ,-.5172121473D-10       ,
     D    .187547032473D-8         ,-.5688440039919D-7     ,
     E    .141662436449235D-5      ,-.283046401495148D-4   ,
     F    .44047862986709951D-3    ,-.51316411610610848D-2 ,
     G    .42319180353336904D-1    ,-.22662499155675492D0  ,
     H    .67561578077218767D0     ,-.76729636288664594D0  ,
     I   -.12869738438135000D0     , .40608211771868508D-1 /
      DATA D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),D(10),D(11),
     1     D(12),D(13),D(14),D(15),D(16),D(17),D(18)   /
     J    .1D-19                   ,-.2D-19                ,
     K    .12D-18                  ,-.58D-18               ,
     L    .305D-17                 ,-.1731D-16             ,
     M    .10668D-15               ,-.72212D-15            ,
     N    .545267D-14              ,-.4684224D-13          ,
     O    .46991955D-12            ,-.570486364D-11        ,
     P    .881689866D-10           ,-.187189074911D-8      ,
     Q    .6177633960644D-7        ,-.398728430048891D-5   ,
     R    .89898983308594085D-3    , .20018060817200274D1  /
      DATA E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10),E(11),
     1     E(12),E(13),E(14),E(15),E(16),E(17),E(18)   /
     S   -.1D-19                   , .3D-19                ,
     T   -.14D-18                  , .65D-18               ,
     U   -.328D-17                 , .1768D-16             ,
     V   -.10269D-15               , .65083D-15            ,
     W   -.456125D-14              , .3596777D-13          ,
     X   -.32643157D-12            , .351521879D-11        ,
     Y   -.4686363688D-10          , .82291933277D-9       ,
     Z   -.2095978138408D-7        , .91386152579555D-6    ,
     1   -.9627723549157079D-4     , .93555574139070650D-1 /
      X=XD
      Y=DABS(X)
      Z=Y*.125D0
      IF(Z.LE.1.0D0)GO TO 10
      Z=1.0D0/Z
      X2=4.0D0*Z*Z-2.0D0
      N1=37
      N2=54
      GO TO 70
   10 IF(Z .LE. 0.0D0)GO TO 78
      X2=4.0D0*Z*Z-2.0D0
      N1=1
      N2=18
   70 DO 80 J=1,2
      Q3=0.0D0
      Q2=0.0D0
      DO 40 I=N1,N2
      Q1=Q2
      Q2=Q3
      Q3=X2*Q2-Q1+A(I)
   40 CONTINUE
      FX=(Q3-Q1)*.5D0
      IF(N1-19)50,51,52
   50 VJ1=FX*Z
      IF(N.LE.0)GO TO 75
      N1=19
      N2=36
      GO TO 80
   52 IF(N1.EQ.55)GO TO 53
      X1=FX
      N1=55
      N2=72
   80 CONTINUE
   78 VJ1=0.0D0
      VY1=-XLG
      GO TO 75
   51 VY1=.6366197723675813D0*(DLOG(Y)*VJ1-1.0D0/Y)+FX*Z
      GO TO 75
   53 X2=DCOS(Y-2.356194490192345D0)
      X3=DSIN(Y-2.356194490192345D0)
      X4=.7978845608028654D0/DSQRT(Y)
      FX=FX*Z
      VJ1=X4*(X1*X2-FX*X3)
      VY1=X4*(FX*X2+X1*X3)
   75 RETURN
      END


c  FFT
      subroutine four1(data,n,isign)

      implicit double precision(a-h,o-z)

      dimension data(2050)
      ip0=2
      ip3=ip0*n
      i3rev=1
      do 50 i3=1,ip3,ip0
      if(i3-i3rev) 10,20,20
   10 tempr=data(i3)
      tempi=data(i3+1)
      data(i3)=data(i3rev)
      data(i3+1)=data(i3rev+1)
      data(i3rev)=tempr
      data(i3rev+1)=tempi
   20 ip1=ip3/2
   30 if(i3rev-ip1) 50,50,40
   40 i3rev=i3rev-ip1
      ip1=ip1/2
      if(ip1-ip0) 50,30,30
   50 i3rev=i3rev+ip1
      ip1=ip0
   60 if(ip1-ip3) 70,100,100
   70 ip2=ip1*2
      theta=3.141592653589793d0 + 3.141592653589793d0
      theta=theta/dble(isign*ip2/ip0)
      sinth=dsin(theta/2.0d0)
      wstpr=-2.0d0*sinth*sinth
      wstpi=dsin(theta)
      wr=1.d0
      wi=0.d0
      do 90 i1=1,ip1,ip0
      do 80 i3=i1,ip3,ip2
      i2a=i3
      i2b=i2a+ip1
      tempr=wr*data(i2b)-wi*data(i2b+1)
      tempi=wr*data(i2b+1)+wi*data(i2b)
      data(i2b)=data(i2a)-tempr
      data(i2b+1)=data(i2a+1)-tempi
      data(i2a)=data(i2a)+tempr
   80 data(i2a+1)=data(i2a+1)+tempi
      tempr=wr
      wr=wr*wstpr-wi*wstpi+wr
   90 wi=wi*wstpr+tempr*wstpi+wi
      ip1=ip2
      go to 60
  100 return
      end



