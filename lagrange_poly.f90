!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine lagrange_any(xi,NGLL,xigll,h,hprime)

! subroutine to compute the Lagrange interpolants based upon the GLL points
! and their first derivatives at any point xi in [-1,1]

  implicit none

  integer NGLL
  double precision xi,xigll(NGLL),h(NGLL),hprime(NGLL)

  integer dgr,i,j
  double precision prod1,prod2

  do dgr=1,NGLL

  prod1 = 1.0d0
  prod2 = 1.0d0
  do i=1,NGLL
    if(i /= dgr) then
      prod1 = prod1*(xi-xigll(i))
      prod2 = prod2*(xigll(dgr)-xigll(i))
    endif
  enddo
  h(dgr)=prod1/prod2

  hprime(dgr)=0.0d0
  do i=1,NGLL
    if(i /= dgr) then
      prod1=1.0d0
      do j=1,NGLL
        if(j /= dgr .and. j /= i) prod1 = prod1*(xi-xigll(j))
      enddo
      hprime(dgr) = hprime(dgr)+prod1
    endif
  enddo
  hprime(dgr) = hprime(dgr)/prod2

  enddo

  end subroutine lagrange_any

!
!=====================================================================
!

! subroutine to compute the derivative of the Lagrange interpolants
! at the GLL points at any given GLL point

  double precision function lagrange_deriv_GLL(I,j,ZGLL,NZ)

!------------------------------------------------------------------------
!
!     Compute the value of the derivative of the I-th
!     Lagrange interpolant through the
!     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j)
!
!------------------------------------------------------------------------

  implicit none

  integer i,j,nz
  double precision zgll(0:nz-1)

  integer degpoly

  double precision, external :: pnleg,pndleg

  degpoly = nz - 1
  if (i == 0 .and. j == 0) then
    lagrange_deriv_GLL = - dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0
  else if (i == degpoly .and. j == degpoly) then
    lagrange_deriv_GLL = dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0
  else if (i == j) then
    lagrange_deriv_GLL = 0.d0
  else
    lagrange_deriv_GLL = pnleg(zgll(j),degpoly) / &
      (pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))) &
      + (1.d0-zgll(j)*zgll(j))*pndleg(zgll(j),degpoly) / (dble(degpoly)* &
      (dble(degpoly)+1.d0)*pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
  endif

  end function lagrange_deriv_GLL

