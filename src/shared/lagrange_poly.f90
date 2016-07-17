!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine lagrange_any(xi,NGLL,xigll,h,hprime)

! subroutine to compute the Lagrange interpolants based upon the GLL points
! and their first derivatives at any point xi in [-1,1]

  implicit none

  double precision,intent(in) :: xi

  integer,intent(in) :: NGLL
  double precision,dimension(NGLL),intent(in) :: xigll
  double precision,dimension(NGLL),intent(out) :: h,hprime

  ! local parameters
  integer :: dgr,i,j
  double precision :: prod1,prod2,prod3
  double precision :: prod2_inv
  double precision :: sum
  double precision :: x0,x

! note: this routine is hit pretty hard by the mesher, optimizing the loops here will be beneficial

  do dgr = 1,NGLL

    prod1 = 1.0d0
    prod2 = 1.0d0

    ! lagrangian interpolants
    x0 = xigll(dgr)
    do i = 1,NGLL
      if (i /= dgr) then
        x = xigll(i)
        prod1 = prod1*(xi-x)
        prod2 = prod2*(x0-x)
      endif
    enddo

    ! takes inverse to avoid additional divisions
    ! (multiplications are cheaper than divisions)
    prod2_inv = 1.d0/prod2

    h(dgr) = prod1 * prod2_inv

    ! first derivatives
    sum = 0.0d0
    do i = 1,NGLL
      if (i /= dgr) then
        prod3 = 1.0d0
        do j = 1,NGLL
          if (j /= dgr .and. j /= i) prod3 = prod3*(xi-xigll(j))
        enddo
        sum = sum + prod3
      endif
    enddo

    hprime(dgr) = sum * prod2_inv

  enddo

  end subroutine lagrange_any

!
!=====================================================================
!

! subroutine to compute the derivative of the Lagrange interpolants
! at the GLL points at any given GLL point

  double precision function lagrange_deriv_GLL(i,j,ZGLL,NZ)

!------------------------------------------------------------------------
!
!     Compute the value of the derivative of the I-th
!     Lagrange interpolant through the
!     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j)
!
!------------------------------------------------------------------------

  implicit none

  integer :: i,j,nz
  double precision :: zgll(0:nz-1)

  ! local parameters
  integer :: degpoly

  double precision, external :: pnleg,pndleg

  degpoly = nz - 1
  if (i == 0 .and. j == 0) then
    lagrange_deriv_GLL = - dble(degpoly)*(dble(degpoly)+1.d0) * 0.25d0  ! / 4.d0
  else if (i == degpoly .and. j == degpoly) then
    lagrange_deriv_GLL = dble(degpoly)*(dble(degpoly)+1.d0) * 0.25d0  ! / 4.d0
  else if (i == j) then
    lagrange_deriv_GLL = 0.d0
  else
    lagrange_deriv_GLL = pnleg(zgll(j),degpoly) / &
      (pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))) &
      + (1.d0-zgll(j)*zgll(j))*pndleg(zgll(j),degpoly) / (dble(degpoly)* &
      (dble(degpoly)+1.d0)*pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
  endif

  end function lagrange_deriv_GLL

