c
c-----------------------------------------------------------------------
c
      double precision function gammaf (x)
c
c=======================================================================
c
c     G a m m a f :
c     -----------
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
c
c---- remove above card for single precision operations
c
      parameter(half=0.5d0,one=1.d0,two=2.d0,four=4.d0)
c
      pi   = four*datan(one)
      gammaf = one
c
      if (x.eq.-half) gammaf = -two*dsqrt(pi)
      if (x.eq. half) gammaf =  dsqrt(pi)
      if (x.eq. one ) gammaf =  one
      if (x.eq. two ) gammaf =  one
      if (x.eq. 1.5d0) gammaf =  dsqrt(pi)/2.d0
      if (x.eq. 2.5d0) gammaf =  1.5d0*dsqrt(pi)/2.d0
      if (x.eq. 3.5d0) gammaf =  2.5d0*1.5d0*dsqrt(pi)/2.d0
      if (x.eq. 3.d0 ) gammaf =  2.d0
      if (x.eq. 4.d0 ) gammaf = 6.d0
      if (x.eq. 5.d0 ) gammaf = 24.d0
      if (x.eq. 6.d0 ) gammaf = 120.d0
c
      return
      end
