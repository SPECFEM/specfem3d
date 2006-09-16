!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
         hprime_xx,hprime_yy,hprime_zz, &
         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz)

  implicit none

  include "constants.h"

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! function for calculating derivatives of Lagrange polynomials
  double precision, external :: lagrange_deriv_GLL

  integer i,j,k,i1,i2,j1,j2,k1,k2

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly ZERO
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then

! calculate derivatives of the Lagrange polynomials
! and precalculate some products in double precision
! hprime(i,j) = h'_i(xigll_j) by definition of the derivation matrix
  do i1=1,NGLLX
    do i2=1,NGLLX
      hprime_xx(i1,i2) = sngl(lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX))
      hprimewgll_xx(i1,i2) = sngl(lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)*wxgll(i2))
    enddo
  enddo

  do j1=1,NGLLY
    do j2=1,NGLLY
      hprime_yy(j1,j2) = sngl(lagrange_deriv_GLL(j1-1,j2-1,yigll,NGLLY))
      hprimewgll_yy(j1,j2) = sngl(lagrange_deriv_GLL(j1-1,j2-1,yigll,NGLLY)*wygll(j2))
    enddo
  enddo

  do k1=1,NGLLZ
    do k2=1,NGLLZ
      hprime_zz(k1,k2) = sngl(lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ))
      hprimewgll_zz(k1,k2) = sngl(lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ)*wzgll(k2))
    enddo
  enddo

  do i=1,NGLLX
    do j=1,NGLLY
      wgllwgll_xy(i,j) = sngl(wxgll(i)*wygll(j))
    enddo
  enddo

  do i=1,NGLLX
    do k=1,NGLLZ
      wgllwgll_xz(i,k) = sngl(wxgll(i)*wzgll(k))
    enddo
  enddo

  do j=1,NGLLY
    do k=1,NGLLZ
      wgllwgll_yz(j,k) = sngl(wygll(j)*wzgll(k))
    enddo
  enddo

  else  ! double precision version

! calculate derivatives of the Lagrange polynomials
! and precalculate some products in double precision
! hprime(i,j) = h'_i(xigll_j) by definition of the derivation matrix
  do i1=1,NGLLX
    do i2=1,NGLLX
      hprime_xx(i1,i2) = lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)
      hprimewgll_xx(i1,i2) = lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)*wxgll(i2)
    enddo
  enddo

  do j1=1,NGLLY
    do j2=1,NGLLY
      hprime_yy(j1,j2) = lagrange_deriv_GLL(j1-1,j2-1,yigll,NGLLY)
      hprimewgll_yy(j1,j2) = lagrange_deriv_GLL(j1-1,j2-1,yigll,NGLLY)*wygll(j2)
    enddo
  enddo

  do k1=1,NGLLZ
    do k2=1,NGLLZ
      hprime_zz(k1,k2) = lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ)
      hprimewgll_zz(k1,k2) = lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ)*wzgll(k2)
    enddo
  enddo

  do i=1,NGLLX
    do j=1,NGLLY
      wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
    enddo
  enddo

  do i=1,NGLLX
    do k=1,NGLLZ
      wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo

  do j=1,NGLLY
    do k=1,NGLLZ
      wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
    enddo
  enddo

  endif

  end subroutine define_derivation_matrices

