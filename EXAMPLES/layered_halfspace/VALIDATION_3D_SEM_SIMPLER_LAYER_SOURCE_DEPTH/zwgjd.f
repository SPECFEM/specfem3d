      subroutine zwgjd (z,w,np,alpha,beta)
c
c=======================================================================
c
c     Z w g j d : Generate np Gauss Jacobi points and weights
c                 associated with Jacobi polynomial of degree n = np-1
c
c=======================================================================
c
c     Note : Coefficients alpha and beta must be greater than -1.
c     ----
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
c
c---- remove above card for single precision operations
c
      parameter(one=1.d0,two=2.d0)
c
      dimension z(*),w(*)
c
c-----------------------------------------------------------------------
c
      n     = np-1
      apb   = alpha+beta

      if (np.le.0) stop 'Minimum number of Gauss points is 1'
      if ((alpha.le.-one).or.(beta.le.-one))
     .   stop 'Alpha and Beta must be greater than -1'

      if (np.eq.1) then
         z(1) = (beta-alpha)/(apb+two)
         w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two)
     .          * two**(apb+one)
         return
      endif

      call jacg (z,np,alpha,beta)

      np1   = n+1
      np2   = n+2
      dnp1  = dble(np1)
      dnp2  = dble(np2)
      fac1  = dnp1+alpha+beta+one
      fac2  = fac1+dnp1
      fac3  = fac2+one
      fnorm = pnormj(np1,alpha,beta)
      rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
      do 100 i=1,np
         call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
         w(i) = -rcoef/(p*pdm1)
 100  continue
      return
      end
