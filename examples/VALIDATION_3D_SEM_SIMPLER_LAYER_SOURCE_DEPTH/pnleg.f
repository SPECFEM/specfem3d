      double precision FUNCTION PNLEG (Z,N)
C-------------------------------------------------------------
C
C     Compute the value of the Nth order Legendre polynomial at Z.
C     Based on the recursion formula for the Legendre polynomials.
C
C-------------------------------------------------------------
      implicit double precision (A-H,O-Z)

      P1   = 1.d0
      P2   = Z
      P3   = P2
      DO 10 K = 1, N-1
         FK  = dble(K)
         P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
         P1  = P2
         P2  = P3
 10   CONTINUE
      PNLEG = P3
      RETURN
      END
