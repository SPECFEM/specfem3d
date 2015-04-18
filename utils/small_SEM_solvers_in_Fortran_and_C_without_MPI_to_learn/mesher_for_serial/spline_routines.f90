!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! compute spline coefficients

  subroutine spline_construction(xpoint,ypoint,npoint,tangent_first_point,tangent_last_point,spline_coefficients)

  implicit none

! tangent to the spline imposed at the first and last points
  double precision, intent(in) :: tangent_first_point,tangent_last_point

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients output by the routine
  double precision, dimension(npoint), intent(out) :: spline_coefficients

  integer :: i

  double precision, dimension(:), allocatable :: temporary_array

  allocate(temporary_array(npoint))

  spline_coefficients(1) = - 1.d0 / 2.d0

  temporary_array(1) = (3.d0/(xpoint(2)-xpoint(1)))*((ypoint(2)-ypoint(1))/(xpoint(2)-xpoint(1))-tangent_first_point)

  do i = 2,npoint-1

    spline_coefficients(i) = ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))-1.d0) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

    temporary_array(i) = (6.d0*((ypoint(i+1)-ypoint(i))/(xpoint(i+1)-xpoint(i)) &
       - (ypoint(i)-ypoint(i-1))/(xpoint(i)-xpoint(i-1)))/(xpoint(i+1)-xpoint(i-1)) &
       - (xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*temporary_array(i-1)) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

  enddo

  spline_coefficients(npoint) = ((3.d0/(xpoint(npoint)-xpoint(npoint-1))) &
      * (tangent_last_point-(ypoint(npoint)-ypoint(npoint-1))/(xpoint(npoint)-xpoint(npoint-1))) &
      - 1.d0/2.d0*temporary_array(npoint-1))/(1.d0/2.d0*spline_coefficients(npoint-1)+1.d0)

  do i = npoint-1,1,-1
    spline_coefficients(i) = spline_coefficients(i)*spline_coefficients(i+1) + temporary_array(i)
  enddo

  deallocate(temporary_array)

  end subroutine spline_construction

! --------------

! evaluate a spline

  subroutine spline_evaluation(xpoint,ypoint,spline_coefficients,npoint,x_evaluate_spline,y_spline_obtained)

  implicit none

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients to use
  double precision, dimension(npoint), intent(in) :: spline_coefficients

! abscissa at which we need to evaluate the value of the spline
  double precision, intent(in):: x_evaluate_spline

! ordinate evaluated by the routine for the spline at this abscissa
  double precision, intent(out):: y_spline_obtained

  integer :: index_loop,index_lower,index_higher

  double precision :: coef1,coef2

! initialize to the whole interval
  index_lower = 1
  index_higher = npoint

! determine the right interval to use, by dichotomy
  do while (index_higher - index_lower > 1)
! compute the middle of the interval
    index_loop = (index_higher + index_lower) / 2
    if(xpoint(index_loop) > x_evaluate_spline) then
      index_higher = index_loop
    else
      index_lower = index_loop
    endif
  enddo

! test that the interval obtained does not have a size of zero
! (this could happen for instance in the case of duplicates in the input list of points)
  if(xpoint(index_higher) == xpoint(index_lower)) stop 'incorrect interval found in spline evaluation'

  coef1 = (xpoint(index_higher) - x_evaluate_spline) / (xpoint(index_higher) - xpoint(index_lower))
  coef2 = (x_evaluate_spline - xpoint(index_lower)) / (xpoint(index_higher) - xpoint(index_lower))

  y_spline_obtained = coef1*ypoint(index_lower) + coef2*ypoint(index_higher) + &
        ((coef1**3 - coef1)*spline_coefficients(index_lower) + &
         (coef2**3 - coef2)*spline_coefficients(index_higher))*((xpoint(index_higher) - xpoint(index_lower))**2)/6.d0

  end subroutine spline_evaluation

