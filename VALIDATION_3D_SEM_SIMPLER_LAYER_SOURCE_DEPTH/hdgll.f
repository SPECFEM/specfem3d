        double precision FUNCTION HDGLL (I,j,ZGLL,NZ)
C-------------------------------------------------------------
C
C     Compute the value of the derivative of the I-th
C     Lagrangian interpolant through the
C     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j).
C
C-------------------------------------------------------------

      implicit double precision (A-H,O-Z)

        integer i,j,nz
        double precision zgll(0:nz-1)

          idegpoly = nz - 1
        if ((i.eq.0).and.(j.eq.0)) then
                hdgll = - dble(idegpoly)*(dble(idegpoly)+1.d0)/4.d0
        else if ((i.eq.idegpoly).and.(j.eq.idegpoly)) then
                hdgll = dble(idegpoly)*(dble(idegpoly)+1.d0)/4.d0
        else if (i.eq.j) then
                hdgll = 0.d0
        else
               rlegendre1 = pnleg(zgll(j),idegpoly)
               rlegendre2 = pndleg(zgll(j),idegpoly)
               rlegendre3 = pnleg(zgll(i),idegpoly)
        hdgll = rlegendre1 / (rlegendre3*(zgll(j)-zgll(i)))
     .  + (1.d0-zgll(j)*zgll(j))*rlegendre2/(dble(idegpoly)*
     .  (dble(idegpoly)+1.d0)*rlegendre3*
     .  (zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
        endif

        return
        end
