      subroutine zwgljd (z,w,np,alpha,beta)
c
c=======================================================================
c
c     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
c     -----------   weights associated with Jacobi polynomials of degree
c                   n = np-1.
c
c=======================================================================
c
c     Note : alpha and beta coefficients must be greater than -1.
c     ----
c            Legendre polynomials are special case of Jacobi polynomials
c            just by setting alpha and beta to 0.
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
c
c---- remove above card for single precision operations
c
      parameter(one=1.d0,two=2.d0)

      dimension z(*), w(*)
c
c-----------------------------------------------------------------------
c
      n     = np-1
      nm1   = n-1

      if (np.le.1) stop 'Minimum number of Gauss-Lobatto points is 2'

      if ((alpha.le.-one).or.(beta.le.-one))
     .      write (*,*) 'Alpha and Beta must be greater than -1'

      if (nm1.gt.0) then
         alpg  = alpha+one
         betg  = beta+one
         call zwgjd (z(2),w(2),nm1,alpg,betg)
      endif
      z(1)  = -one
      z(np) =  one
      do  110 i=2,np-1
         w(i) = w(i)/(one-z(i)**2)
  110 continue
      call jacobf (p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
      w(1)  = endw1 (n,alpha,beta)/(two*pd)
      call jacobf (p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
      w(np) = endw2 (n,alpha,beta)/(two*pd)
c
      return
      end
