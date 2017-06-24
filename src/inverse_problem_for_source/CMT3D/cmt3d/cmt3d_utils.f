c This subroutine takes the value of a moment tensor
c and outputs the strike-dip-rake of the best-fitting 
c double couple, and their moment. 

	 subroutine mij(m,epsilon,s1,d1,r1,s2,d2,r2,m0,m00,iflag)
	 implicit  REAL*8(A-H,O-Z)

	 dimension a(3,3),d(3),v(3,3),ae(3),an(3)
	 dimension ae1(3),an1(3),an2(3),ae2(3)
	 real*8  m(6), m0, m00
	 common/mpa/ft,fd,fl,ft1,fd1,fl1,azt,aiht,azp,aihp,az0,aih0,
     *d,clvd

         if (iflag == 0) then ! Hiroo's convention
            a(1,1)=m(1)
            a(2,2)=m(2)
            a(3,3)=m(3)
            a(1,2)=m(4)
            a(1,3)=m(5)
            a(2,3)=m(6)
            a(2,1)=a(1,2)
            a(3,1)=a(1,3)
            a(3,2)=a(2,3)
         else                 ! Harvard's convention
            a(1,1)=m(2)
            a(2,2)=m(3)
            a(3,3)=m(1)
            a(1,2)=m(6)
            a(1,3)=-m(4)
            a(2,3)=-m(5)
            a(2,1)=a(1,2)
            a(3,1)=a(1,3)
            a(3,2)=a(2,3)
         endif

c        calculate the full moment
         m00 =0.
         do i = 1, 3
            do j = 1, 3
               m00 = m00 + a(i,j) ** 2
            enddo
         enddo
         m00 = sqrt(m00 /2 )
c       take eigen values and eigen vectors
	 call jacobi(a,3,3,d,v,nrot)
11	 format(1x,i3,3(f5.2,2x))

c       subtract isotropic part
         sum = 0.
         do i = 1, 3
            sum = sum + d(i)
         enddo
         do i = 1, 3
            d(i) = d(i) - sum/3
         enddo
c      calculate the devitoric part
	 amax=d(1)
	 amin=d(1)
	 absmax=abs(d(1))
	 absmin=abs(d(1))
	 imax = 1
	 imin = 1
	 jmax = 1
	 jmin = 1
	 do i=2,3
            if(d(i).ge.amax) then
               amax=d(i)
               imax=i
            end if
            if(d(i).le.amin) then
               amin=d(i)
               imin=i
            end if
            if(abs(d(i)).ge.absmax) then
               absmax=abs(d(i))
               jmax = i
            endif
            if(abs(d(i)).le.absmin) then
               absmin=abs(d(i))
               jmin = i
            endif
         enddo
         do i=1,3
            if(i.ne.imax.and.i.ne.imin) k=i
         enddo

c        calculate moment and normalize tensor
	 m0=sqrt((amax*amx+amin*amin+d(k)*d(k))/2.)
	 mxx=mxx/m0
	 myy=myy/m0
	 mzz=mzz/m0
	 mxy=mxy/m0
	 mxz=mxz/m0
	 myz=myz/m0

c        calculate epsilon
	 epsilon = - d(jmin) / d(jmax)

c        obtain strike dip and rake
	do 990 ii=1,3
990	d(ii)=d(ii)/m0
	dt=d(imax)
	dp=d(imin)
	d0=d(k)
	d(1)=dt
	d(2)=dp
	d(3)=d0
	in0=k
	 do 22 i=1,3
   	 ae(i)=(v(i,imax)+v(i,imin))/sqrt(2.0)
   	 an(i)=(v(i,imax)-v(i,imin))/sqrt(2.0)
22	 continue
	 aer=sqrt(ae(1)**2+ae(2)**2+ae(3)**2)
	 anr=sqrt(an(1)**2+an(2)**2+an(3)**2)
	 do 23 i=1,3
	 ae(i)=ae(i)/aer
23	 an(i)=an(i)/anr
	aaa=1.0e-06
	if(abs(an(3)).lt.aaa) an(3)=0.
	if(abs(ae(3)).lt.aaa) ae(3)=0.
	if(an(3).le.0.) then
	do 81 i=1,3
	an1(i)=an(i)
 81	ae1(i)=ae(i)
	if(ae(3).le.0.) then
	do 82 i=1,3
	an2(i)=ae(i)
 82	ae2(i)=an(i)
	else
	do 83 i=1,3
	an2(i)=-ae(i)
 83	ae2(i)=-an(i)
	end if
	else
	do 84 i=1,3
	an1(i)=-an(i)
 84	ae1(i)=-ae(i)
	do 85 i=1,3
	do 86 j=1,3
 86	v(j,i)=-v(j,i)
 85 	continue
	if(ae(3).le.0.) then
	do 87 i=1,3
	an2(i)=ae(i)
 87	ae2(i)=an(i)
	else
	do 88 i=1,3
	an2(i)=-ae(i)
 88	ae2(i)=-an(i)
	end if
	end if
	call tdl(an1,ae1,ft,fd,fl)
	call tdl(an2,ae2,ft1,fd1,fl1)
        s1 = 360-ft
        d1 = fd
        r1 = 180-fl
        s2 = 360-ft1
        d2 = fd1
        r2 = 180-fl1

c      calculate the TPN axes for the moment
	do 44 i=1,3
	an1(i)=v(i,in0)
	an(i)=v(i,imax)
 44	ae(i)=v(i,imin)
	call azih(an,azt,aiht)
	call azih(ae,azp,aihp)
	call azih(an1,az0,aih0)
c	write(*,*) 'az,ih,for t-axis:  ',azt,aiht
c	write(*,*) 'az,ih,for p-axis:  ',azp,aihp
c	write(*,*) 'az,ih,for 0-axis:  ',az0,aih0
        return
        end

      subroutine jacobi(a,n,np,d,v,nrot)
      implicit  REAL*8(A-H,O-Z)
      parameter (nmax=100)
      dimension a(np,np),d(np),v(np,np),b(nmax),z(nmax)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     *         .and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      stop '50 iterations should never happen'
      return
      end

	subroutine tdl(an,bn,ft,fd,fl)
	implicit  REAL*8(A-H,O-Z)
	dimension an(3),bn(3)
	xn=an(1)
	yn=an(2)
	zn=an(3)
	xe=bn(1)
	ye=bn(2)
	ze=bn(3)
	aaa=1.0e-06
	con=57.2957795
	if(abs(zn).lt.aaa) then
	fd=90.
	axn=abs(xn)
	write(*,*) 'axn  ',axn
	if(axn.gt.1.0) axn=1.0
	ft=asin(axn)*con
	st=-xn
	ct=yn
c	write(*,*) 'st,ct  ',st,ct
cif(abs(st).lt.aaa) st=0.
cif(abs(ct).lt.aaa) ct=0.
	if(st.ge.0..and.ct.lt.0) ft=180.-ft
	if(st.lt.0..and.ct.le.0) ft=180.+ft
	if(st.lt.0..and.ct.gt.0) ft=360.-ft
	fl=asin(abs(ze))*con
	sl=-ze
	if(abs(xn).lt.aaa) then
	cl=xe/yn
	else
	cl=-ye/xn
	end if
cif(abs(sl).lt.aaa) sl=0.
cif(abs(cl).lt.aaa) cl=0.
	if(sl.ge.0..and.cl.lt.0) fl=180.-fl
	if(sl.lt.0..and.cl.le.0) fl=fl-180.
	if(sl.lt.0..and.cl.gt.0) fl=-fl
	else
	if(-zn.gt.1.0) zn=-1.0
	fdh=acos(-zn)
	fd=fdh*con
cwrite(*,*) 'fd  ',fd
	sd=sin(fdh)
	if (sd.eq.0)then
	write(*,*) 'sd=0 ! stop'
	stop
	end if
	st=-xn/sd
	ct=yn/sd
	sx=abs(st)
	if(sx.gt.1.0) sx=1.0
	ft=asin(sx)*con
cif(abs(st).lt.aaa) st=0.
cif(abs(ct).lt.aaa) ct=0.
	if(st.ge.0..and.ct.lt.0) ft=180.-ft
	if(st.lt.0..and.ct.le.0) ft=180.+ft
	if(st.lt.0..and.ct.gt.0) ft=360.-ft
cwrite(*,*) 'ft  ',ft
	sl=-ze/sd
	sx=abs(sl)
	if(sx.gt.1.0) sx=1.0
	fl=asin(sx)*con
	if(st.eq.0) then
	cl=xe/ct
	else
	xxx=yn*zn*ze/sd/sd+ye
	cl=-sd*xxx/xn
	if(ct.eq.0) cl=ye/st
	end if
cif(abs(sl).lt.aaa) sl=0.
cif(abs(cl).lt.aaa) cl=0.
	if(sl.ge.0..and.cl.lt.0) fl=180.-fl
	if(sl.lt.0..and.cl.le.0) fl=fl-180.
	if(sl.lt.0..and.cl.gt.0) fl=-fl
cwrite(*,*) 'fl  ',fl
	end if
	return
	end

	subroutine azih(aa,az,aih)
	implicit  REAL*8(A-H,O-Z)
	dimension aa(3)
	con=57.2957795
	gx=aa(1)
	gy=aa(2)
	gz=aa(3)
	if(gz.lt.0.) then
	gx=-gx
	gy=-gy
	gz=-gz
	end if
	if(abs(gx).lt.0.00001) then
	az=270.
	if(gy.ge.0.) az=90.
	else
 100	az=atan(gy/gx)*con
	if(gx.gt.0..and.gy.lt.0.) az=az+360.
	if(gx.lt.0.) az=az+180.
	end if
	if(gz.gt.1.0) gz=1.0
	 aih=asin(gz)*con
	 if(gz.lt.0.00001) aih=0.
cwrite(*,*)'gz:  ',gz,aih
	if(az.lt.0.) az=az+360.
	return
	end

c*****************************************************************************

	subroutine sdr2moment(sphif,sdlt,slmda,moment,mrr,mtt,mpp,mrt,mrp,mtp)
	
	implicit none
c       
c       this subroutine converts eqk fault parameters to moment tensors
c       for definitions see hiroo''s class notes
c       
	real*8 sphif,sdlt,slmda,moment
	real*8 mrr,mtt,mpp,mrt,mrp,mtp
	real*8 pi,phif,dlt,lmda
c       
	pi=3.1415926
	phif=sphif*pi/180
	dlt=sdlt*pi/180
	lmda=slmda*pi/180
	mrr=(sin(2*dlt)*sin(lmda))*moment
	mtt=(-(sin(phif))**2*sin(2*dlt)*sin(lmda)-sin(2*phif)*cos(lmda)*sin(dlt))*moment
	mpp=(-(cos(phif))**2*sin(2*dlt)*sin(lmda)+sin(2*phif)*cos(lmda)*sin(dlt))*moment
	mrt=-(cos(phif)*cos(dlt)*cos(lmda)+sin(phif)*sin(lmda)*cos(2*dlt))*moment
	mrp=-(-sin(phif)*cos(dlt)*cos(lmda)+cos(phif)*sin(lmda)*cos(2*dlt))*moment
      	mtp=(-0.5*sin(2*phif)*sin(2*dlt)*sin(lmda)-cos(2*phif)*cos(lmda)*sin(dlt))*moment
        return
	end
