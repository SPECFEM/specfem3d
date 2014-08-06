      subroutine distaz(elat,elon,slat,slon,azm,bzm,ddg,dkm)
c----------------------------------------------------------------------
c
c   Computes distance and azimuth between two points on the surface
c   of the earth.
c
c   Input via call:
c   elat, elon -=- coordinates of 1st point (decimal degrees)
c   slat, slon -=- coordinates of 2nd point (decimal degrees)
c
c   Output via call:
c   azm -=- event to station azimuth in decimal degrees
c   bzm -=- station to event azimuth in decimal degrees
c   ddg -=- epicentral distance in decimal degrees
c   dkm -=- epicentral distance in kilometers
c
c   Sign convention: North and East = positive
c		     South and West = negative
c
c   Note: Single precision version.
c
c----------------------------------------------------------------------
      c1=57.29578
      c2=1.570796
      c3=.9932313
      c4=.0174533
      a=c2-atan(c3*tan(c4*elat))
      b=-elon
      c=c2-atan(c3*tan(c4*slat))
      d=-slon
      delo=b-d
      if (delo .lt. -180.) delo=360.+delo
      if (delo .gt.  180.) delo=delo-360.
      delo=delo*c4
      del=acosf(cos(a)*cos(c)+sin(a)*sin(c)*cos(delo))
      ddg=c1*del
      dkm=(6371.227*(1.+.0033785*(1./3.-cos((a+c)/2.)**2)))*del
      e=acosf((cos(c)-cos(a)*cos(del))/(sin(a)*sin(del)))
      s=acosf((cos(a)-cos(c)*cos(del))/(sin(c)*sin(del)))
      if (delo .lt. 0.) then
        azm=360.-c1*e
        bzm=c1*s
      else
        azm=c1*e
        bzm=360.-c1*s
      endif
      return
      end
      function acosf(arg)
      if (abs(arg) .gt. 1.e-8) goto 1
      acosf=1.570796
      goto 3
1     if (arg .gt. 0.) goto 2
      if (arg .lt. -1.) arg=-.9999999
      acosf=3.141593-atan((sqrt(1.-arg**2))/abs(arg))
      goto 3
2     if (arg .gt. 1.) arg=.9999999
      acosf=atan((sqrt(1.-arg**2))/arg)
3     return
      end
