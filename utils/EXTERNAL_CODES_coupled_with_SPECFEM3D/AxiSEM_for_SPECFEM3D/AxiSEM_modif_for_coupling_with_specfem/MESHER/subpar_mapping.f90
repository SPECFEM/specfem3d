!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
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
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!=======================
module subpar_mapping  
!=======================
!
! Module used to compute the  
! subparametric mapping that defines the mesh. 

  use global_parameters

  implicit none 

  public :: jacobian_subpar, alpha_subpar, beta_subpar
  public :: gamma_subpar, delta_subpar, epsilon_subpar, zeta_subpar
  public :: alphak_subpar, betak_subpar, gammak_subpar
  public :: deltak_subpar, epsilonk_subpar, zetak_subpar
  public :: jacobian_srf_subpar,  quadfunc_map_subpar, grad_quadfunc_map_subpar
  public :: mgrad_pointwise_subpar, mgrad_pointwisek_subpar
  public :: mapping_subpar, s_over_oneplusxi_axis_subpar
  public :: compute_partial_d_subpar
  private

contains

!-----------------------------------------------------------------------------------------
real(kind=dp)    function mapping_subpar(xil,etal,nodes_crd,iaxis)
! This routines computes the coordinates along the iaxis axis of the image of
! any point in the reference domain in the physical domain.
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
  
  integer          :: iaxis
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer          :: inode
  real(kind=dp)    :: shp(8)

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8(xil,etal,shp)
  
  mapping_subpar = zero

  do inode = 1, 8
     mapping_subpar = mapping_subpar + shp(inode)*nodes_crd(inode,iaxis)
  end do

end function mapping_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function quadfunc_map_subpar(p, s, z, nodes_crd)
! This routines computes the quadratic functional (s-s(xi,eta))**2 + (z-z(xi,eta))**2

  real(kind=dp)    :: p(2),xil,etal,s,z,nodes_crd(8,2)
  
  xil  = p(1)
  etal = p(2)

  quadfunc_map_subpar = (s-mapping_subpar(xil,etal,nodes_crd,1))**2 &
                      + (z-mapping_subpar(xil,etal,nodes_crd,2))**2

end function quadfunc_map_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine grad_quadfunc_map_subpar(grd, p, s, z, nodes_crd)
! This routine returns the gradient of the quadratic functional associated with
! the mapping.

  real(kind=dp)    :: grd(2),p(2),xil, etal, s,z, nodes_crd(8,2)
  real(kind=dp)    :: shpder(8,2),a,d,b,c
  integer          :: inode

  xil  = p(1)
  etal = p(2)

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)

  a = zero ; b = zero ; c = zero ; d = zero

     do inode = 1, 8
        a = a + nodes_crd(inode,1)*shpder(inode,1)
        d = d + nodes_crd(inode,2)*shpder(inode,2)
        b = b + nodes_crd(inode,1)*shpder(inode,2)
        c = c + nodes_crd(inode,2)*shpder(inode,1)
     end do

  grd(1) = -((s-mapping_subpar(xil,etal,nodes_crd,1))*a&
            +(z-mapping_subpar(xil,etal,nodes_crd,2))*c)  
  grd(2) = -((s-mapping_subpar(xil,etal,nodes_crd,1))*b&
            +(z-mapping_subpar(xil,etal,nodes_crd,2))*d)

end subroutine grad_quadfunc_map_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_partial_d_subpar(dsdxi, dzdxi, dsdeta, dzdeta, xil, etal, nodes_crd)

  real(kind=dp)   , intent(out) :: dsdxi,dzdxi,dsdeta,dzdeta
  real(kind=dp)   , intent(in) :: xil,etal,nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2)

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  dsdxi = zero ; dzdxi = zero ; dsdeta = zero ; dzdeta = zero

  do inode = 1, 8
     dsdxi  =  dsdxi + nodes_crd(inode,1)*shpder(inode,1)
     dzdeta = dzdeta + nodes_crd(inode,2)*shpder(inode,2)
     dsdeta = dsdeta + nodes_crd(inode,1)*shpder(inode,2)
     dzdxi  =  dzdxi + nodes_crd(inode,2)*shpder(inode,1)
  end do

end subroutine compute_partial_d_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function s_over_oneplusxi_axis_subpar(xil, etal, nodes_crd)
! This routine returns the value of the quantity
!  
!              s/(1+xi) 
!
! when the associated element lies along the axis of 
! symmetry. Again, in this routine, we assume that the 
! spectral element is a 8-node Serendipity element :
! 
! 7 - - - 6 - - - 5
! |       ^       |       Control points 1,8, and 7 
! |   eta |       | belong to the axis of symmetry. 
! |       |       |        
! 8        --->   4
! |        xi     |
! |               |
! |               |
! 1 - - - 2 - - - 3 .
 
  
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: shp(2:6)
  integer :: inode 
 
  s_over_oneplusxi_axis_subpar = zero

  shp(2) =  half*(one-xil)*(one-etal)
  shp(3) = quart*(one-etal)*(xil-etal-one)
  shp(4) =  half*(one-etal**2)
  shp(5) = quart*(one+etal)*(xil+etal-one) 
  shp(6) =  half*(one-xil)*(one+etal)

  do inode = 2, 6
     s_over_oneplusxi_axis_subpar = s_over_oneplusxi_axis_subpar &
                                   +nodes_crd(inode,1)*shp(inode) 
  end do

end function s_over_oneplusxi_axis_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function jacobian_subpar(xil, etal, nodes_crd)
! This routines the value of the Jacobian (that is, 
! the determinant of the Jacobian matrix), for any point
! inside a given element. IT ASSUMES 8 nodes 2D isoparametric
! formulation of the geometrical transformation and therefore
! requires the knowledge of the coordinated of the 8 control
! points, which are defined as follows :
!
!     7 - - - 6 - - - 5
!     |       ^       |
!     |   eta |       |
!     |       |       |
!     8        --->   4
!     |        xi     |
!     |               |
!     |               |
!     1 - - - 2 - - - 3 .
!
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  
  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)          
     d = d + nodes_crd(inode,2)*shpder(inode,2)          
     b = b + nodes_crd(inode,1)*shpder(inode,2)          
     c = c + nodes_crd(inode,2)*shpder(inode,1)          
  end do 

  jacobian_subpar = a*d - b*c 


end function jacobian_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function jacobian_srf_subpar(xil, crdedge)
! This routine computes the Jacobian of the transformation
! that maps [-1,+1] into a portion of the boundary of domain.  
!
!         xi
!        ---->
! 1 - - - 2 - - - 3 .
!
  real(kind=dp)    :: xil, crdedge(3,2)
  real(kind=dp)    :: dsdxi,dzdxi,s1,s2,s3,z1,z2,z3 

  s1 = crdedge(1,1) ; s2 = crdedge(2,1) ; s3 = crdedge(3,1) 
  z1 = crdedge(1,2) ; z2 = crdedge(2,2) ; z3 = crdedge(3,2) 

  dsdxi = s1*(xil-half) + s2*(-two*xil) + s3*(xil+half)
  dzdxi = z1*(xil-half) + z2*(-two*xil) + z3*(xil+half)

  jacobian_srf_subpar = dsqrt(dsdxi**2+dzdxi**2)


end function jacobian_srf_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function alphak_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of 
!
!    alphak =  ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is 
! defined by the isoparametric transformation involving eight control
! nodes.J is the determinant of the Jacobian matrix of the
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer          :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)

  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = one/(a*d - b*c)

  alphak_subpar  = -inv_jacob*a*b

end function alphak_subpar 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function betak_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!     betak =   ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. beta is defined within an element, and s(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the 
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer          :: inode
  real(kind=dp)    :: shpder(8,2), a, d, b, c, inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)

  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = one/(a*d - b*c)

  betak_subpar   = inv_jacob*(a**2)

end function betak_subpar 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function gammak_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!     gammak =  ( ds/deta ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. gamma is defined within an element, and s(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.
!
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  
  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  gammak_subpar   = inv_jacob*(b**2)

end function gammak_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function deltak_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!     deltak =  - ( dz/dxi ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. delta is defined within an element, and (s,z)(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.
  
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  
  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  deltak_subpar   = -inv_jacob*c*d

end function deltak_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function epsilonk_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!   epsilonk =  ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. epsilon is defined within an element, and (s,z)(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)

  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  epsilonk_subpar   =  inv_jacob*(c**2)

end function epsilonk_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function zetak_subpar(xil,etal,nodes_crd)
!
! This routines returns the value of
!
!   zetak =  ( dz/deta ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. zeta is defined within an element, and (s,z)(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode 
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob
 
  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
 
  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  zetak_subpar   =  inv_jacob*(d**2)

end function zetak_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function alpha_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of 
!
!    alpha =  s(xi,eta) * ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is 
! defined by the isoparametric transformation involving eight control
! nodes.J is the determinant of the Jacobian matrix of the
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer          :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  
  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  alpha_subpar  = -mapping_subpar(xil,etal,nodes_crd,1)*inv_jacob*a*b

end function alpha_subpar 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function  beta_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!     beta =  s(xi,eta) * ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. beta is defined within an element, and s(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the 
! transformation.
!
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)

  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  beta_subpar   = mapping_subpar(xil,etal,nodes_crd,1)*inv_jacob*(a**2)

end function beta_subpar 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function  gamma_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!     gamma =  s(xi,eta) * ( ds/deta ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. gamma is defined within an element, and s(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.
!
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)

  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  gamma_subpar   = mapping_subpar(xil,etal,nodes_crd,1)*inv_jacob*(b**2)

end function gamma_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function delta_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!     delta =  - s(xi,eta) * ( dz/dxi ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. delta is defined within an element, and (s,z)(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  
  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  delta_subpar   = -mapping_subpar(xil,etal,nodes_crd,1)*inv_jacob*c*d

end function delta_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function epsilon_subpar(xil, etal, nodes_crd)
!
! This routines returns the value of
!
!   epsilon =  s(xi,eta) * ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. epsilon is defined within an element, and (s,z)(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  
  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  epsilon_subpar   =  mapping_subpar(xil,etal,nodes_crd,1)*inv_jacob*(c**2)

end function epsilon_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function zeta_subpar(xil,etal,nodes_crd)
!
! This routines returns the value of
!
!   zeta =  s(xi,eta) * ( dz/deta ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. zeta is defined within an element, and (s,z)(xi,eta) is
! defined by the isoparametric transformation involving eight control
! nodes. J is the determinant of the Jacobian matrix of the
! transformation.

  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer :: inode 
  real(kind=dp)    :: shpder(8,2),a,d,b,c,inv_jacob

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)

  a = zero
  b = zero
  c = zero
  d = zero

  do inode = 1, 8
     a = a + nodes_crd(inode,1)*shpder(inode,1)
     d = d + nodes_crd(inode,2)*shpder(inode,2)
     b = b + nodes_crd(inode,1)*shpder(inode,2)
     c = c + nodes_crd(inode,2)*shpder(inode,1)
  end do

  inv_jacob  = 1./(a*d - b*c)

  zeta_subpar   =  mapping_subpar(xil,etal,nodes_crd,1)*inv_jacob*(d**2)

end function zeta_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine mgrad_pointwise_subpar(mg, xil, etal, nodes_crd)
!
! This routines returns the following matrix:
!                      +                     +
!                      |(ds/dxi)  | (ds/deta)|
!    mg =  s(xi,eta) * | ---------|--------- |(xi,eta)
!                      |(dz/dxi ) | (dz/deta)|
!                      +                     +
! This 2*2 matrix is needed when defining and storing
! gradient/divergence related arrays.


  real(kind=dp)    :: mg(2,2)
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer          :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c

  mg(:,:) = zero

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  a = zero ; b = zero ; c = zero ; d = zero

     do inode = 1, 8
        a = a + nodes_crd(inode,1)*shpder(inode,1)
        d = d + nodes_crd(inode,2)*shpder(inode,2)
        b = b + nodes_crd(inode,1)*shpder(inode,2)
        c = c + nodes_crd(inode,2)*shpder(inode,1)
     end do

  mg(1,1)  = mapping_subpar(xil,etal,nodes_crd,1)*a
  mg(1,2)  = mapping_subpar(xil,etal,nodes_crd,1)*b 
  mg(2,1)  = mapping_subpar(xil,etal,nodes_crd,1)*c
  mg(2,2)  = mapping_subpar(xil,etal,nodes_crd,1)*d

end subroutine mgrad_pointwise_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine mgrad_pointwisek_subpar(mg, xil, etal, nodes_crd)
!
! This routines returns the following matrix:
!            +                     +
!            |(ds/dxi)  | (ds/deta)|
!    mg =    | ---------|--------- |(xi,eta)
!            |(dz/dxi ) | (dz/deta)|
!            +                     +
! This 2*2 matrix is needed when defining and storing
! gradient/divergence related arrays.

  real(kind=dp)    :: mg(2,2)
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  integer          :: inode
  real(kind=dp)    :: shpder(8,2),a,d,b,c

  mg(:,:) = zero

  ! Compute the appropriate derivatives of the shape
  ! functions

  call shp8der(xil,etal,shpder)
  a = zero ; b = zero ; c = zero ; d = zero

     do inode = 1, 8
        a = a + nodes_crd(inode,1)*shpder(inode,1)
        d = d + nodes_crd(inode,2)*shpder(inode,2)
        b = b + nodes_crd(inode,1)*shpder(inode,2)
        c = c + nodes_crd(inode,2)*shpder(inode,1)
     end do

  mg(1,1)  = a
  mg(1,2)  = b 
  mg(2,1)  = c
  mg(2,2)  = d

end subroutine mgrad_pointwisek_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine shp8(xil,etal,shp)
!
! This routine computes and returns the quadratic
! shape functions axixiociated with a 8-nodes serendip
! element for a given point of coordinates (xi,eta).
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

  real(kind=dp)    :: xil, etal
  real(kind=dp)    :: shp(8)
  real(kind=dp)    :: xip,xim,etap,etam,xixi,etaeta

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
subroutine shp8der(xil,etal,shpder)
!
! This routine computes and returns the derivatives
! of the shape functions axixiociated with a 8-nodes serendip
! element for a given point of coordinates (xi,eta).
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
! shpder(:,1) : derivative wrt xi
! shpder(:,2) : derivative wrt eta

  real(kind=dp)    :: xil, etal
  real(kind=dp)    :: shpder(8,2)
  real(kind=dp)    :: xip,xim,etap,etam,xixi,etaeta
 
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


!==========================
end module subpar_mapping
!==========================
