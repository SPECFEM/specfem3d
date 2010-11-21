      double precision function pnormj (n,alpha,beta)
c
c=======================================================================
c
c     P n o r m j
c     -----------
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
c
c---- remove above card for single precision operation
c
      one   = 1.d0
      two   = 2.d0
      dn    = dble(n)
      const = alpha+beta+one
      if (n.le.1) then
         prod   = gammaf(dn+alpha)*gammaf(dn+beta)
         prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
         pnormj = prod * two**const/(two*dn+const)
         return
      endif
      prod  = gammaf(alpha+one)*gammaf(beta+one)
      prod  = prod/(two*(one+const)*gammaf(const+one))
      prod  = prod*(one+alpha)*(two+alpha)
      prod  = prod*(one+beta)*(two+beta)
      do 100 i=3,n
         dindx = dble(i)
         frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
         prod  = prod*frac
 100  continue
      pnormj = prod * two**const/(two*dn+const)
c
      return
      end
