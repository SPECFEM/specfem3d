      subroutine jacobf (poly,pder,polym1,pderm1,polym2,pderm2,n,alp,
     .                   bet,x)
c
c=======================================================================
c
c     J a c o b f : Compute the Jacobi polynomial and its derivative
c     -----------   of degree n at x.
c
c=======================================================================
      implicit double precision (a-h,o-z)
c
c---- remove above card for single precision operations
c
      apb  = alp+bet
      poly = 1.d0
      pder = 0.d0
c
      if (n .eq. 0) return
c
      polyl = poly
      pderl = pder
      poly  = (alp-bet+(apb+2.d0)*x)/2.d0
      pder  = (apb+2.d0)/2.d0
      if (n .eq. 1) return
      do 20 k=2,n
         dk = dble(k)
         a1 = 2.d0*dk*(dk+apb)*(2.d0*dk+apb-2.d0)
         a2 = (2.d0*dk+apb-1.d0)*(alp**2-bet**2)
         b3 = (2.d0*dk+apb-2.d0)
         a3 = b3*(b3+1.d0)*(b3+2.d0)
         a4 = 2.d0*(dk+alp-1.d0)*(dk+bet-1.d0)*(2.d0*dk+apb)
         polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
         pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
         psave  = polyl
         pdsave = pderl
         polyl  = poly
         poly   = polyn
         pderl  = pder
         pder   = pdern
 20   continue
      polym1 = polyl
      pderm1 = pderl
      polym2 = psave
      pderm2 = pdsave
c
      return
      end
