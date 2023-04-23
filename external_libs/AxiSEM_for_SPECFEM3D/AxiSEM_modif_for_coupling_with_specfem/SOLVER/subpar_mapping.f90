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
module subpar_mapping
! < Module used to compute the subparametric mapping that defines the mesh.

  use global_parameters

  implicit none

  public :: mapping_subpar
  public :: compute_partial_d_subpar
  private

contains

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function mapping_subpar(xil,etal,nodes_crd,iaxis)
! < This routines computes the coordinates along the iaxis axis of the image of
!! any point in the reference domain in the physical domain.
!
! 7 - - - 6 - - - 5
! |       ^       |
! |   eta |       |
! |       |       |
! 8        --->   4
! |        xi     |
! |               |
! |               |
! 1 - - - 2 - - - 3 .
!
! iaxis = 1 : along the cylindrical radius axis
! iaxis = 2 : along the vertical(rotation) axis

  integer, intent(in)       :: iaxis
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  integer                   :: inode
  real(kind=dp)             :: shp(8)

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8(xil,etal,shp)

  mapping_subpar = zero

  do inode = 1, 8
     mapping_subpar = mapping_subpar + shp(inode)*nodes_crd(inode,iaxis)
  enddo

end function mapping_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine compute_partial_d_subpar(dsdxi, dzdxi, dsdeta, dzdeta, xil, etal, nodes_crd)

  real(kind=dp), intent(out) :: dsdxi,dzdxi,dsdeta,dzdeta
  real(kind=dp), intent(in)  :: xil,etal,nodes_crd(8,2)
  integer                    :: inode
  real(kind=dp)              :: shpder(8,2)

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  dsdxi = zero ; dzdxi = zero ; dsdeta = zero ; dzdeta = zero

  do inode = 1, 8
     dsdxi  =  dsdxi + nodes_crd(inode,1)*shpder(inode,1)
     dzdeta = dzdeta + nodes_crd(inode,2)*shpder(inode,2)
     dsdeta = dsdeta + nodes_crd(inode,1)*shpder(inode,2)
     dzdxi  =  dzdxi + nodes_crd(inode,2)*shpder(inode,1)
  enddo

end subroutine compute_partial_d_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine shp8(xil,etal,shp)
!
! < This routine computes and returns the quadratic
!! shape functions axixiociated with a 8-nodes serendip
!! element for a given point of coordinates (xi,eta).
!
! Topology is defined as follows
!
! 7 - - - 6 - - - 5
! |       ^       |
! |   eta |       |
! |       |       |
! 8        --->   4
! |        xi     |
! |               |
! |               |
! 1 - - - 2 - - - 3
!

  real(kind=dp), intent(in)  :: xil, etal
  real(kind=dp), intent(out) :: shp(8)
  real(kind=dp)              :: xip, xim, etap, etam, xixi, etaeta

  shp(:) = zero


  xip    = one +  xil
  xim    = one -  xil
  etap   = one + etal
  etam   = one - etal
  xixi   =  xil *  xil
  etaeta = etal * etal

  ! Corners first:
  shp(1) = quart * xim * etam * (xim + etam - three)
  shp(3) = quart * xip * etam * (xip + etam - three)
  shp(5) = quart * xip * etap * (xip + etap - three)
  shp(7) = quart * xim * etap * (xim + etap - three)

  ! Then midpoints:
  shp(2) = half  * etam * (one -   xixi)
  shp(4) = half  *  xip * (one - etaeta)
  shp(6) = half  * etap * (one -   xixi)
  shp(8) = half  *  xim * (one - etaeta)

end subroutine shp8
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine shp8der(xil,etal,shpder)
!
! < This routine computes and returns the derivatives
!! of the shape functions axixiociated with a 8-nodes serendip
!! element for a given point of coordinates (xi,eta).
!! shpder(:,1) : derivative wrt xi
!! shpder(:,2) : derivative wrt eta
!
! Topology is defined as follows
!
! 7 - - - 6 - - - 5
! |       ^       |
! |   eta |       |
! |       |       |
! 8        --->   4
! |        xi     |
! |               |
! |               |
! 1 - - - 2 - - - 3
!
!

  real(kind=dp), intent(in)  :: xil, etal
  real(kind=dp), intent(out) :: shpder(8,2)
  real(kind=dp)              :: xip, xim, etap, etam, xixi, etaeta

  shpder(:,:) = zero

  xip    = one +  xil
  xim    = one -  xil
  etap   = one + etal
  etam   = one - etal
  xixi   =  xil *  xil
  etaeta = etal * etal

  ! Corners first:
  shpder(1,1) = -quart * etam * ( xim + xim + etam - three)
  shpder(1,2) = -quart *  xim * (etam + xim + etam - three)
  shpder(3,1) =  quart * etam * ( xip + xip + etam - three)
  shpder(3,2) = -quart *  xip * (etam + xip + etam - three)
  shpder(5,1) =  quart * etap * ( xip + xip + etap - three)
  shpder(5,2) =  quart *  xip * (etap + xip + etap - three)
  shpder(7,1) = -quart * etap * ( xim + xim + etap - three)
  shpder(7,2) =  quart *  xim * (etap + xim + etap - three)

  ! Then midside points :
  shpder(2,1) = -one  * xil * etam
  shpder(2,2) = -half * (one - xixi)
  shpder(4,1) =  half * (one - etaeta)
  shpder(4,2) = -one  * etal * xip
  shpder(6,1) = -one  * xil * etap
  shpder(6,2) =  half * (one - xixi)
  shpder(8,1) = -half * (one - etaeta)
  shpder(8,2) = -one  * etal * xim

end subroutine shp8der
!-----------------------------------------------------------------------------------------

end module subpar_mapping
!=========================================================================================
