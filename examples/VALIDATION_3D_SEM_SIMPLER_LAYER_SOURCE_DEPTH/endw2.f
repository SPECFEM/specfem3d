c
c-----------------------------------------------------------------------
c
      double precision function endw2 (n,alpha,beta)
c
c=======================================================================
c
c     E n d w 2 :
c     ---------
c
c=======================================================================
c
      implicit double precision (a-h,o-z)
c
c---- remove above card for single precision operations
c
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
c
c-----------------------------------------------------------------------
c
      apb   = alpha+beta
      if (n.eq.0) then
         endw2 = zero
         return
      endif
      f1   = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
      f1   = f1*(apb+two)*two**(apb+two)/two
      if (n.eq.1) then
         endw2 = f1
         return
      endif
      fint1 = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
      fint1 = fint1*two**(apb+two)
      fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
      fint2 = fint2*two**(apb+three)
      f2    = (two*(alpha+two)*fint1 - (apb+four)*fint2)
     $        * (apb+three)/four
      if (n.eq.2) then
         endw2 = f2
         return
      endif
      do 100 i=3,n
         di   = dble(i-1)
         abn  = alpha+beta+di
         abnn = abn+di
         a1   =  -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
         a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
         a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
         f3   =  -(a2*f2+a1*f1)/a3
         f1   = f2
         f2   = f3
 100  continue
      endw2  = f3
      return
      end
