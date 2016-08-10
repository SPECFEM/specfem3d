module wave2d_define_der_matrices

  use wave2d_constants

contains

  subroutine define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,wgllwgll_xz)

  implicit none

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz

! array with all the weights in the square
  double precision, dimension(NGLLX,NGLLZ) :: wgll_square

! function for calculating derivatives of Lagrange polynomials
  double precision, external :: lagrange_deriv_GLL

  integer i1,i2,k1,k2,i,k

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if (mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = 0.d0
  if (mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = 0.d0

! calculate derivatives of the Lagrange polynomials
! and precalculate some products in double precision
! hprime(i,j) = h'_i(xigll_j) by definition of the derivative matrix
  do i1=1,NGLLX
    do i2=1,NGLLX
      hprime_xx(i1,i2) = lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)
    enddo
  enddo

  do k1=1,NGLLZ
    do k2=1,NGLLZ
      hprime_zz(k1,k2) = lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ)
    enddo
  enddo

  do i=1,NGLLX
    do k=1,NGLLZ
      wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo

end subroutine define_derivative_matrices

end module wave2d_define_der_matrices

