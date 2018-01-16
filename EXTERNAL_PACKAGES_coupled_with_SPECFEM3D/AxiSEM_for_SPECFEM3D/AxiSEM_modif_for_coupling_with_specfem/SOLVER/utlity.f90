!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
module utlity

  use global_parameters
  implicit none

  public :: compute_coordinates, scoord, zcoord, rcoord, thetacoord
  public :: dblreldiff_small, reldiff_small
  public :: dblereldiff, reldiff
  public :: dbleabsreldiff, absreldiff
  public :: to_lower
  private

contains

!-----------------------------------------------------------------------------------------
pure logical function dblreldiff_small(x1,x2)

  real(kind=dp), intent(in) :: x1,x2

  dblreldiff_small = .false.

  if (x1 /= zero) then
     if (abs((x1-x2)/x1) <= smallval_dble) dblreldiff_small = .true.
  else if (x2 /= zero) then
     if (abs((x1-x2)/x2) <= smallval_dble) dblreldiff_small = .true.
  else
     dblreldiff_small = .true.
  endif

end function dblreldiff_small
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure logical function reldiff_small(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2
  real(kind=realkind)             ::  smallval1

  if (realkind == sp) smallval1 = smallval_sngl
  if (realkind == dp) smallval1 = smallval_dble

  reldiff_small = .false.

  if (x1 /= zero) then
     if (abs((x1-x2)/x1) <= smallval1) reldiff_small = .true.
  else if (x2 /= zero) then
     if (abs((x1-x2)/x2) <= smallval1) reldiff_small = .true.
  else
     reldiff_small = .true.
  endif

end function reldiff_small
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=realkind) function reldiff(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2

  if (x1 /= zero) then
     reldiff=(x1-x2)/x1
  else if (x2 /= zero) then
     reldiff=(x1-x2)/x2
  else
     reldiff=zero
  endif

end function reldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function dblereldiff(x1,x2)

  real(kind=dp), intent(in) :: x1,x2

  if (x1 /= zero) then
     dblereldiff=(x1-x2)/x1
  else if (x2 /= zero) then
     dblereldiff=(x1-x2)/x2
  else
     dblereldiff=zero
  endif

end function dblereldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=realkind) function absreldiff(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2

  if (x1 /= zero) then
     absreldiff=abs((x1-x2)/x1)
  else if (x2 /= zero) then
     absreldiff=abs((x1-x2)/x2)
  else
     absreldiff=zero
  endif

end function absreldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function dbleabsreldiff(x1,x2)

  real(kind=dp), intent(in) :: x1,x2

  if (x1 /= zero) then
     dbleabsreldiff=abs((x1-x2)/x1)
  else if (x2 /= zero) then
     dbleabsreldiff=abs((x1-x2)/x2)
  else
     dbleabsreldiff=zero
  endif

end function dbleabsreldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine compute_coordinates(s,z,r,theta,ielem,ipol,jpol)
! < Given the elemental grid point index, outputs s,z,r,theta coordinate [m,rad].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).

  use data_mesh, only: min_distance_dim
  use data_mesh, only: lnods, crd_nodes, axis
  use data_spec, only: xi_k, eta
  use analytic_mapping, only: mapping

  real(kind=dp), intent(out)    :: s,z,r,theta
  integer, intent(in)           :: ielem,ipol,jpol
  integer                       :: ipt,inode
  real(kind=dp)                 :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  enddo

  ! Fill global coordinate array
  if ( axis(ielem) ) then

     s= mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
     z= mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)

  else
     s= mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
     z= mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  endif

  ! Eliminate roundoff errors
  if (abs(s) < min_distance_dim) s=zero
  if (abs(z) < min_distance_dim) z=zero

  r = dsqrt(s**2+z**2)
  theta = datan(s/(z+epsi))
  if ( zero > theta ) theta = pi + theta
  if (theta == zero .and. z < 0) theta = pi

end subroutine compute_coordinates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function scoord(ipol,jpol,ielem)
! < Given the elemental grid point index, outputs the s coordinate [m].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).

  use data_mesh, only: min_distance_dim
  use data_mesh, only: lnods, crd_nodes, axis
  use data_spec, only: xi_k, eta
  use analytic_mapping, only: mapping

  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  enddo

  ! Fill global coordinate array
  if ( axis(ielem) ) then
     scoord = mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
  else
     scoord = mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
  endif

  ! Eliminate roundoff errors
  if (abs(scoord) < min_distance_dim) scoord=zero

end function scoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function zcoord(ipol,jpol,ielem)
! < Given the elemental grid point index, outputs the z coordinate [m].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).

  use data_mesh, only: min_distance_dim
  use data_mesh, only: lnods, crd_nodes, axis
  use data_spec, only: xi_k, eta
  use analytic_mapping, only: mapping

  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  enddo

  ! Fill global coordinate array
  if ( axis(ielem) ) then
     zcoord = mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)
  else
     zcoord = mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  endif

  ! Eliminate roundoff errors
  if (abs(zcoord) < min_distance_dim) zcoord=zero

end function zcoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function rcoord(ipol,jpol,ielem)
! < Given the elemental grid point index, outputs the radius coordinate [m].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).

  use data_mesh, only: min_distance_dim
  use data_mesh, only: lnods, crd_nodes, axis
  use data_spec, only: xi_k, eta
  use analytic_mapping, only: mapping

  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2),s,z

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  enddo

  ! Fill global coordinate array
  if ( axis(ielem) ) then
     s = mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)
  else
     s = mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  endif

  rcoord = sqrt(s**2 + z**2)
  ! Eliminate roundoff errors
  if (abs(rcoord) < min_distance_dim) rcoord=zero

end function rcoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function thetacoord(ipol,jpol,ielem)
! < Given the elemental grid point index, outputs the theta coordinate [rad].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).

  use data_mesh, only: min_distance_dim
  use data_mesh, only: lnods, crd_nodes,axis
  use data_spec, only: xi_k, eta
  use analytic_mapping, only: mapping

  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2),s,z

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  enddo

  ! Fill global coordinate array
  if ( axis(ielem) ) then
     s = mapping(xi_k(ipol), eta(jpol), nodes_crd, 1, ielem)
     z = mapping(xi_k(ipol), eta(jpol), nodes_crd, 2, ielem)
  else
     s = mapping(eta(ipol), eta(jpol), nodes_crd, 1, ielem)
     z = mapping(eta(ipol), eta(jpol), nodes_crd, 2, ielem)
  endif

  thetacoord = datan(s/(z+epsi))
  if ( zero > thetacoord ) thetacoord = pi + thetacoord
  if (thetacoord == zero .and. z < 0) thetacoord = pi

end function thetacoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function to_lower(strIn) result(strOut)
! < Converts string to lowercase, adapted from http://www.star <= ac.uk/~cgp/Fortran.html
    implicit none

    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer                      :: i,j

    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j >= iachar("A") .and. j <= iachar("Z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))+32)
        else
            strOut(i:i) = strIn(i:i)
        endif
    enddo

end function to_lower
!-----------------------------------------------------------------------------------------

end module utlity
!=========================================================================================
