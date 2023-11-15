      subroutine jacg (xjac,np,alpha,beta)
c
c=======================================================================
c
c     J a c g : Compute np Gauss points, which are the zeros of the
c               Jacobi polynomial with parameter alpha and beta.
c
c=======================================================================
c
c     Note :
c     ----
c          .Alpha and Beta determines the specific type of gauss points.
c                  .alpha = beta =  0.0  ->  Legendre points
c                  .alpha = beta = -0.5  ->  Chebyshev points
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
c
c---- remove above card for single precision operations
c
      double precision xjac(*)
c
      data kstop/10/
      data eps/1.0d-12/
c
c-----------------------------------------------------------------------
c
      xlast = 0.d0
      n   = np-1
      dth = 4.d0*datan(1.d0)/(2.d0*dble(n)+2.d0)
      do 40 j=1,np
         if (j.eq.1) then
            x = dcos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
         else
            x1 = dcos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
            x2 = xlast
            x  = (x1+x2)/2.d0
         endif
         do 30 k=1,kstop
            call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
            recsum = 0.d0
            jm = j-1
            do 29 i=1,jm
               recsum = recsum+1.d0/(x-xjac(np-i+1))
 29         continue
            delx = -p/(pd-recsum*p)
            x    = x+delx
            if (abs(delx) .lt. eps) goto 31
 30      continue
 31      continue
         xjac(np-j+1) = x
         xlast        = x
 40   continue
      do 200 i=1,np
         xmin = 2.d0
         do 100 j=i,np
            if (xjac(j).lt.xmin) then
               xmin = xjac(j)
               jmin = j
            endif
 100     continue
         if (jmin.ne.i) then
            swap = xjac(i)
            xjac(i) = xjac(jmin)
            xjac(jmin) = swap
         endif
 200  continue
      return
      end
