cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine convspc
     &     ( max_nstation,maxnfreq,n_station,
     &       time_series_length,n_frequency,omega_imag,
     &       station_displacement,
     &       source_depth,source_lat,source_lon,station_lat,station_lon,
     &       work_spc,work_time,sac_file )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cc DK DK compute displacement or velocity seismograms
      logical DISPLACEMENT
      parameter (DISPLACEMENT = .false.)

c variables for input/output
      integer max_nstation,maxnfreq,n_frequency,n_station
      real*8 time_series_length,omega_imag
      real*8 source_depth,source_lat,source_lon
      real*8 station_lat(max_nstation),station_lon(max_nstation)
      complex*16 station_displacement(3,max_nstation,maxnfreq)
      real work_time(32*maxnfreq)
      complex*16 work_spc(16*maxnfreq)
      character*80 sac_file(max_nstation)
c other variables
      integer i_station,icomp,i_frequency,ismooth,n1,m1,nn,isign,i
      integer lnblnk
      real*8 t
      character*8 cext(3)
      character*80 outfile

      data ismooth / 8 /
      data cext / '.bhz.sac','.bhr.sac','.bht.sac' /

      do 1000 i_station=1,n_station
        do 900 icomp=1,3
          work_spc(1) = dcmplx(0.d0)

cc DK DK following a suggestion by Takeuchi, I changed the code here to compute
cc DK DK displacement instead of velocity
          if (DISPLACEMENT) then
c compute displacement seismograms (divide the spectrum by i\omega)
            do i_frequency=1,n_frequency
              work_spc(i_frequency+1)
     &              = station_displacement(icomp,i_station,i_frequency)
     &    / dcmplx(0.d0, 2.d0 * 3.1415926535897932d0 * dble(i_frequency)
     &                    / time_series_length)
            enddo

          else
c compute velocity seismograms
            do i_frequency=1,n_frequency
              work_spc(i_frequency+1)
     &              = station_displacement(icomp,i_station,i_frequency)
            enddo
          endif

          do 210 i_frequency=n_frequency+2,ismooth*n_frequency+1
            work_spc(i_frequency) = dcmplx(0.d0)
  210     continue
          do 220 i_frequency=1,ismooth*n_frequency-1
            n1 = ismooth*n_frequency + i_frequency + 1
            m1 = ismooth*n_frequency - i_frequency + 1
            work_spc(n1) = dconjg( work_spc(m1) )
  220     continue
          nn = 2 * ismooth*n_frequency
          isign = 1
          call four1( work_spc,nn,isign )
          do 230 i=-nn+1,nn
            t = time_series_length
     &                           * dble(i-1) / dble(nn)
            work_time(i+nn)
     &              = real(
     &                  dble( work_spc( mod(i-1+nn,nn)+1 ) )
     &                  * dexp( omega_imag * t )
     &                )
c amplitude correction
c --- for FFT
          work_time(i+nn) = work_time(i+nn) / real( time_series_length )
c --- km -> m
            work_time(i+nn) = work_time(i+nn) * 1.e3
  230     continue
          outfile
     &            = sac_file(i_station)(1:lnblnk(sac_file(i_station)))
     &              //cext(icomp)
          call sacconv( 2*time_series_length,2*nn,
     &                   -time_series_length,work_time,
     &                    source_depth,source_lat,source_lon,
     &                    station_lat(i_station),station_lon(i_station),
     &                    icomp,outfile )
  900   continue
 1000 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sacconv( tlen,np,b0,y,d0,theta0,phi0,theta,phi,
     &                          icomp,output )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'header.h'
      include 'equiv.h'
c
      integer np,icomp
      real*8 tlen,b0,d0
      real y(*)
      real*8 theta0,phi0,theta,phi
c variables for files
      character*80 output
c
      call inithdr
      call compprm( tlen,np,b0,
     &                    d0,theta0,phi0,theta,phi,
     &                    begin,delta,npts,stla,stlo,evdp,evla,evlo )
      call wsact0( begin,delta,npts,stla,stlo,evdp,evla,evlo,icomp,
     &                   y,output,
     &                   depmin,depmax,depmen,ennd,
     &                   dist,az,baz,gcarc,cmpaz,cmpinc,
     &                   nvhdr,ninf,nhst,iftype,leven,lovrok,lcalda )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compprm( tlen,np,beg0,
     &                      d0,theta0,phi0,theta,phi,
     &                      begin,delta,npts,stla,stlo,evdp,evla,evlo )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
      integer np,npts
      real*8 tlen,beg0,d0,theta0,phi0,theta,phi
      real begin,delta,stla,stlo,evdp,evla,evlo
c
      begin = real(beg0)
      npts = np
      delta = tlen / real(npts)
      stla = real( theta )
      stlo = real( phi )
      evdp = real( d0 )
      evla = real( theta0 )
      evlo = real( phi0 )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inithdr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'header.h'
      integer i
c
      do 100 i=1,mfhdr
        fhdr(i) = -12345.0
  100 continue
      do 110 i=1,mnhdr
        nhdr(i) = -12345
  110 continue
      do 120 i=1,mihdr
        ihdr(i) = -12345
  120 continue
      do 130 i=1,mlhdr
        lhdr(i) = 0
  130 continue
      do 140 i=1,mkhdr
        khdr(i) = '-12345  '
  140 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wsact0( begin,delta,npts,stla,stlo,evdp,evla,evlo,
     &                      icmpt,data,filen,
     &                      depmin,depmax,depmen,ennd,
     &                      dist,az,baz,gcarc,cmpaz,cmpinc,
     &                      nvhdr,ninf,nhst,iftype,leven,lovrok,lcalda )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Generating SAC-binary file for evenly sampled data.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'header.h'
c required variables
      integer npts,icmpt
      real begin,delta,stla,stlo,evdp,evla,evlo,data(*)
      character*80 filen
c computed variables
      real depmin,depmax,depmen,ennd,dist,az,baz,gcarc,cmpaz,cmpinc
      integer nvhdr,ninf,nhst,iftype,leven,lovrok,lcalda
c other variables
      integer i,ier,filenchk,sacbw
c
c     icmpt
c       3: Transverse Component
c       otherwise: not supported by this version
c
c
c
      depmax = data(1)
      depmin = data(2)
      depmen = 0.d0
      do 100 i=1,npts
        if ( data(i).gt.depmax ) depmax = data(i)
        if ( data(i).lt.depmin ) depmin = data(i)
        depmen = depmen + data(i)
  100 continue
      depmen = depmen / real(npts)
c
      ennd = begin + delta * real(npts-1)
c
      az = 0.0
      baz = 0.0
      dist = 0.0
      call distaz1( evla,evlo,stla,stlo,1,dist,az,baz,gcarc,ier )
c
      if ( icmpt.eq.1 ) then
        cmpaz = 0.0
        cmpinc = 0.0
      else
      if ( icmpt.eq.2 ) then
        cmpaz = baz + 180.0
        if ( cmpaz.ge.360.0 ) cmpaz = cmpaz - 360.0
        cmpinc = 90.0
      else
      if ( icmpt.eq.3 ) then
        cmpaz = baz + 90.0
        if ( cmpaz.ge.360.0 ) cmpaz = cmpaz - 360.0
        cmpinc = 90.0
      else
        stop 'Illegal component type (wsact0)'
      endif
      endif
      endif
      nvhdr = 6
      ninf = 0
      nhst = 0
      iftype = 1
      leven = 1
      lovrok = 1
      lcalda = 1
c
      ier = filenchk( 80,filen )
      ier = sacbw( mfhdr,mnhdr,mihdr,mlhdr,mkhdr,
     &                   fhdr,nhdr,ihdr,lhdr,khdr,
     &                   npts,data,filen )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine distaz1(the,phe,ths,phs,ns,dist,az,baz,xdeg,nerr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*=====================================================================
* PURPOSE:  To compute the distance and azimuth between locations.
*=====================================================================
* INPUT ARGUMENTS:
*    THE:     Event latitude in decimal degrees, North positive. [r]
*    PHE:     Event longitude, East positive. [r]
*    THS:     Array of station latitudes. [r]
*    PHS:     Array of station longitudes. [r]
*    NS:      Length of THS and PHS. [i]
*=====================================================================
* OUTPUT ARGUMENTS:
*    DIST:    Array of epicentral distances in km. [r]
*    AZ:      Array of azimuths in degrees. [r]
*    BAZ:     Array of back azimuths. [r]
*    XDEG:    Array of great circle arc lengths. [r]
*    NERR:    Error flag:
*             =    0   No error.
*             = 0904   Calculation failed internal consistency checks.
*=====================================================================
* MODULE/LEVEL:  DFM/4
*=====================================================================
* GLOBAL INPUT:
*    MACH:
*=====================================================================
* SUBROUTINES CALLED:
*    SACLIB:  SETMSG, APCMSG
*=====================================================================
* LOCAL VARIABLES:
*=====================================================================
* KNOWN ERRORS:
* - Problem with equation for distance. See discussion below.
*=====================================================================
      real pi,todeg,torad
      parameter (PI=3.141592654)
      parameter (TODEG=57.29577950)
      parameter (TORAD=1./TODEG)
c
      integer ns,nerr,i
      real ths(ns), phs(ns)
      real dist(ns), az(ns), baz(ns), xdeg(ns)
      logical laz,lbaz,ldist,lxdeg
      real the,phe,ec2,onemec2,eps,temp,therad,pherad,thg
      real d,e,f,c,a,b,g,h,thsrad,phsrad
      real d1,e1,f1,c1,a1,b1,g1,h1,sc,sd,ss,t1,p1,t2,p2,el
      real costhi,costhk,sinthi,sinthk,tanthi,tanthk,al,dl
      real a12top,a12bot,a12,cosa12,sina12
      real e2,e3,c0,c2,c4,v1,v2,z1,z2,x2,y2
      real e1p1,sqrte1p1,u1bot,u1,u2top,u2bot,u2,b0,du,pdist
      real rad,fl,twopideg,degtokm
      real c00,c01,c02,c03,c21,c22,c23,c42,c43
* PROCEDURE:
c
* - Calculations are based upon the reference spheroid of 1968 and
*   are defined by the major radius (RAD) and the flattening (FL).
c
      data rad/6378.160/,fl/0.00335293/
      data twopideg/360./
      data c00,c01,c02,c03/1.,0.25,-4.6875e-02,1.953125e-02/
      data c21,c22,c23/-0.125,3.125e-02,-1.46484375e-02/
      data c42,c43/-3.90625e-03,2.9296875e-03/
      data degtokm/111.3199/
c
* - Initialize.
c
      nerr=0
      ec2=2.*fl-fl*fl
      onemec2=1.-ec2
      eps=1.+ec2/onemec2
c
* - Check which output items are required.
c
      laz=.true.
      if(az(1).lt.0.)laz=.false.
      lbaz=.true.
      if(baz(1).lt.0.)lbaz=.false.
      ldist=.true.
      if(dist(1).lt.0.)ldist=.false.
      lxdeg=.true.
      if(xdeg(1).lt.0.)lxdeg=.false.
c
* - Convert event location to radians.
*   (Equations are unstable for latidudes of exactly 0 degrees.)
c
      temp=the
      if(temp.eq.0.)temp=1.0e-08
      therad=torad*temp
      pherad=torad*phe
c
* - Must convert from geographic to geocentric coordinates in order
*   to use the spherical trig equations.  This requires a latitude
*   correction given by: 1-EC2=1-2*FL+FL*FL
c
      thg=atan(onemec2*tan(therad))
      d=sin(pherad)
      e=-cos(pherad)
      f=-cos(thg)
      c=sin(thg)
      a= f*e
      b=-f*d
      g=-c*e
      h=c*d
c
* - Loop on stations:
c
      do 5000 i=1,ns
c
* -- Convert to radians.
        temp=ths(i)
        if(temp.eq.0.)temp=1.0e-08
        thsrad=torad*temp
        phsrad=torad*phs(i)
c
* -- Calculate some trig constants.
        thg=atan(onemec2*tan(thsrad))
        d1=sin(phsrad)
        e1=-cos(phsrad)
        f1=-cos(thg)
        c1=sin(thg)
        a1=f1*e1
        b1=-f1*d1
        g1=-c1*e1
        h1=c1*d1
        sc=a*a1+b*b1+c*c1
c
* - Spherical trig relationships used to compute angles.
c
        if(lxdeg)then
          sd=0.5*sqrt(((a-a1)**2+(b-b1)**2+(c-c1)**2)*((a+a1)**2
     #       +(b+b1)**2+(c+c1)**2))
          xdeg(i)=atan2(sd,sc)*todeg
          if(xdeg(i).lt.0.)xdeg(i)=xdeg(i)+twopideg
        endif
        if(laz)then
          ss = ((a1-d)**2+(b1-e)**2+c1**2-2.)
          sc = ((a1-g)**2+(b1-h)**2+(c1-f)**2-2.)
          az(i)=atan2(ss,sc)*todeg
          if(az(i).lt.0.)az(i)=az(i)+twopideg
        endif
        if(lbaz)then
          ss=((a-d1)**2+(b-e1)**2+c**2-2.)
          sc=((a-g1)**2+(b-h1)**2+(c-f1)**2-2.)
          baz(i)=atan2(ss,sc)*todeg
          if(baz(i).lt.0.)baz(i)=baz(i)+twopideg
        endif
c
* - Now compute the distance between the two points using Rudoe's
*   formula given in GEODESY, section 2.15(b).
*   (There is some numerical problem with the following formulae.
*   If the station is in the southern hemisphere and the event in
*   in the northern, these equations give the longer, not the
*   shorter distance between the two locations.  Since the equations
*   are fairly messy, the simplist solution is to reverse the
*   meanings of the two locations for this case.)
        if(ldist)then
          if(thsrad.gt.0.)then
            t1=thsrad
            p1=phsrad
            t2=therad
            p2=pherad
          else
            t1=therad
            p1=pherad
            t2=thsrad
            p2=phsrad
          endif
          el=ec2/onemec2
          e1=1.+el
          costhi=cos(t1)
          costhk=cos(t2)
          sinthi=sin(t1)
          sinthk=sin(t2)
          tanthi=sinthi/costhi
          tanthk=sinthk/costhk
          al=tanthi/(e1*tanthk)+
     #       ec2*sqrt((e1+tanthi**2)/(e1+tanthk**2))
          dl=p1-p2
          a12top=sin(dl)
          a12bot=(al-cos(dl))*sinthk
          a12=atan2(a12top,a12bot)
          cosa12=cos(a12)
          sina12=sin(a12)
          e1=el*((costhk*cosa12)**2+sinthk**2)
          e2=e1*e1
          e3=e1*e2
          c0=c00+c01*e1+c02*e2+c03*e3
          c2=c21*e1+c22*e2+c23*e3
          c4=c42*e2+c43*e3
          v1=rad/sqrt(1.-ec2*sinthk**2)
          v2=rad/sqrt(1.-ec2*sinthi**2)
          z1=v1*(1.-ec2)*sinthk
          z2=v2*(1.-ec2)*sinthi
          x2=v2*costhi*cos(dl)
          y2=v2*costhi*sin(dl)
          e1p1=e1+1.
          sqrte1p1=sqrt(e1p1)
          u1bot=sqrte1p1*cosa12
          u1=atan2(tanthk,u1bot)
          u2top=v1*sinthk+e1p1*(z2-z1)
          u2bot=sqrte1p1*(x2*cosa12-y2*sinthk*sina12)
          u2=atan2(u2top,u2bot)
          b0=v1*sqrt(1.+el*(costhk*cosa12)**2)/e1p1
          du=u2 -u1
          pdist=b0*(c2*(sin(2.*u2)-sin(2.*u1))+
     #       c4*(sin(4.*u2)-sin(4.*u1)))
          dist(i)=abs(b0*c0*du+pdist)
          if(lxdeg .and. (abs(dist(i)-degtokm*xdeg(i))).gt.100.)then
            nerr=0904
c            call setmsg('ERROR',nerr)
c            call apimsg(i)
          endif
        endif
 5000   continue
c
 8888 return
c
*=====================================================================
* MODIFICATION HISTORY:
*    830603:  Fixed bug with negative station latiudes.
*    810000:  Original version.
*=====================================================================
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------------------------------------------------
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8 DATA(*)
      INTEGER NN,ISIGN
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      REAL*8 TEMPR,TEMPI
      INTEGER N,I,J,M,MMAX,ISTEP

      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
