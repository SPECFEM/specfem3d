!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
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

!==========================
module analytic_mapping
!==========================

  use global_parameters
  use analytic_spheroid_mapping
  use analytic_semi_mapping
  use subpar_mapping
  use data_mesh

  implicit none

  public :: jacobian_anal
  public :: alpha_anal,  beta_anal,  gamma_anal,  delta_anal,  epsilon_anal,  zeta_anal
  public :: alphak_anal, betak_anal, gammak_anal, deltak_anal, epsilonk_anal, zetak_anal
  public :: jacobian_srf_anal,  quadfunc_map_anal, grad_quadfunc_map_anal
  public :: mgrad_pointwise_anal, mgrad_pointwisek_anal
  public :: mapping_anal, s_over_oneplusxi_axis_anal

  public :: Ms_z_eta_s_xi,   Ms_z_eta_s_eta
  public :: Ms_z_xi_s_eta,   Ms_z_xi_s_xi
  public :: Ms_z_eta_s_xi_k, Ms_z_eta_s_eta_k
  public :: Ms_z_xi_s_eta_k, Ms_z_xi_s_xi_k

  public :: compute_partial_derivatives

  private

contains

!-----------------------------------------------------------------------------------------
real(kind=dp) function mapping_anal(xil, etal, nodes_crd, iaxis, ielem0)
! This routine computes the coordinates along the iaxis axis
! of the image of any point in the reference domain in the physical domain
! using the implicit assumption that the domain is spheroidal.

  integer,intent(in)          :: iaxis,ielem0
  real(kind=dp)   ,intent(in) :: xil, etal, nodes_crd(8,2)

  if (eltype(ielem0) == 'curved') &
     mapping_anal = map_spheroid(xil,etal,nodes_crd,iaxis)
  if (eltype(ielem0) == 'linear') &
     mapping_anal = mapping_subpar(xil,etal,nodes_crd,iaxis)
  if (eltype(ielem0) == 'semino') &
     mapping_anal = map_semino(xil,etal,nodes_crd,iaxis)
  if (eltype(ielem0) == 'semiso') &
     mapping_anal = map_semiso(xil,etal,nodes_crd,iaxis)

end function mapping_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function quadfunc_map_anal(p, s, z, nodes_crd, ielem0)
! This routines computes the quadratic functional
! (s-s(xi,eta))**2 + (z-z(xi,eta))**2

  integer :: ielem0
  real(kind=dp)    :: p(2), xil,etal,s,z, nodes_crd(8,2)

  xil  = p(1)
  etal = p(2)

  quadfunc_map_anal = (s-mapping_anal(xil,etal,nodes_crd,1,ielem0))**2 &
                    + (z-mapping_anal(xil,etal,nodes_crd,2,ielem0))**2

end function quadfunc_map_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine grad_quadfunc_map_anal(grd, p, s, z, nodes_crd, ielem0)
! This routine returns the gradient of the quadratic functional
! associated with the mapping.

  integer :: ielem0
  real(kind=dp)    :: grd(2),p(2),xil, etal, s,z, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdxi,dsdeta,dzdeta

  xil  = p(1)
  etal = p(2)

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,&
                                   nodes_crd,ielem0)

  grd(1) = -((s-mapping_anal(xil,etal,nodes_crd,1,ielem0))*dsdxi&
            +(z-mapping_anal(xil,etal,nodes_crd,2,ielem0))*dzdxi)
  grd(2) = -((s-mapping_anal(xil,etal,nodes_crd,1,ielem0))*dsdeta&
            +(z-mapping_anal(xil,etal,nodes_crd,2,ielem0))*dzdeta)

end subroutine grad_quadfunc_map_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function s_over_oneplusxi_axis_anal(xil, etal, nodes_crd, ielem0)
! This routine returns the value of the quantity
!
!              s/(1+xi)
!
! when the associated element lies along the axis of
! symmetry, in the case of an analytical transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdxi,dsdeta,dzdeta

  if ( xil == -one ) then
     ! Apply L'Hopital's rule
     call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                      nodes_crd,ielem0)
     s_over_oneplusxi_axis_anal = dsdxi
  else
     s_over_oneplusxi_axis_anal = mapping_anal(xil,etal,nodes_crd,1,ielem0) / &
                                  (one+xil)
  endif

end function s_over_oneplusxi_axis_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function jacobian_anal(xil, etal, nodes_crd, ielem0)
! This function returns the value of the jacobian of the
! analytical mapping between the reference square [-1,1]^2 and
! the deformed element in the spheroid.

  integer,intent(in)          :: ielem0
  real(kind=dp)   ,intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)               :: dsdxi,dzdxi,dsdeta,dzdeta

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  jacobian_anal = dsdxi*dzdeta-dsdeta*dzdxi

end function jacobian_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function jacobian_srf_anal(xil, crdedge)
! This routine computes the Jacobian of the transformation
! that maps [-1,+1] into a portion of the boundary of domain.
!
!         xi(or eta)
!        ---->
! 1 - - - 2 - - - 3 .

  real(kind=dp)    :: xil, crdedge(3,2)
  real(kind=dp)    :: dsdxi,dzdxi,s1,s2,s3,z1,z2,z3
  real(kind=dp)    :: thetabar,deltatheta,a,b
  real(kind=dp)    :: arg,dist

  s1 = crdedge(1,1) ; s2 = crdedge(2,1) ; s3 = crdedge(3,1)
  z1 = crdedge(1,2) ; z2 = crdedge(2,2) ; z3 = crdedge(3,2)

  call compute_parameters_srf(s1,s3,z1,z3,a,b,deltatheta,thetabar)

  if (dabs(deltatheta) > 1.0d-10 ) then
     arg = xil*deltatheta + thetabar
     dsdxi = -deltatheta*a*dsin(arg)
     dzdxi =  deltatheta*b*dcos(arg)
     jacobian_srf_anal = dsqrt(dsdxi**2+dzdxi**2)
  else
     dist = dsqrt((z3-z1)**2+(s3-s1)**2)
     jacobian_srf_anal = half*dist
  endif

end function jacobian_srf_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function alphak_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    alphak =  ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  alphak_anal  = -inv_jacob*dsdxi*dsdeta

end function alphak_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function betak_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    betak =  ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  betak_anal  = inv_jacob*dsdxi**2

  end function betak_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function gammak_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    gammak =  ( ds/deta ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  gammak_anal  = inv_jacob*dsdeta**2

end function gammak_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function deltak_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    deltak = -( dz/dxi ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  deltak_anal  = -inv_jacob*dzdxi*dzdeta

end function deltak_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function epsilonk_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    epsilonk = ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  epsilonk_anal  = inv_jacob*dzdxi**2

end function epsilonk_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function zetak_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    zetak_anal = ( dz/deta ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  zetak_anal  = inv_jacob*dzdeta**2

end function zetak_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function alpha_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    alpha = s(xi,eta) * ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  alpha_anal  = -inv_jacob*dsdxi*dsdeta*mapping_anal(xil,etal,nodes_crd,1,ielem0)

  end function alpha_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function beta_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    beta =  s(xi,eta) * ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  beta_anal  = inv_jacob*dsdxi**2*mapping_anal(xil,etal,nodes_crd,1,ielem0)

  end function beta_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function gamma_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    gamma = s(xi,eta) * ( ds/deta ) * ( ds/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  gamma_anal  = inv_jacob*dsdeta**2*mapping_anal(xil,etal,nodes_crd,1,ielem0)

end function gamma_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function delta_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    delta = -s(xi,eta) * ( dz/dxi ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  delta_anal  = -inv_jacob*dzdxi*dzdeta*mapping_anal(xil,etal,nodes_crd,1,ielem0)

end function delta_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function epsilon_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    epsilon = s(xi,eta) * ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  epsilon_anal  = inv_jacob*dzdxi**2*mapping_anal(xil,etal,nodes_crd,1,ielem0)

end function epsilon_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function zeta_anal(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    zeta_anal = s(xi,eta) * ( dz/deta ) * ( dz/deta) / J(xi,eta),
!
! a quantity that is needed in the calculation of the laplacian
! operator. alpha is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,&
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  zeta_anal  = inv_jacob*dzdeta**2*mapping_anal(xil,etal,nodes_crd,1,ielem0)

  end function zeta_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_eta_s_xi(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    Ms_z_eta_s_xi = s(xi,eta) / J(xi,eta) * ( ds/dxi ) * ( dz/deta)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the FIRST TERM OF dsdz_0
!          in the THIRD TERM OF dzds_0
! It is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_xi=inv_jacob*dsdxi*dzdeta*mapping_anal(xil,etal,nodes_crd,1,ielem0)

end function Ms_z_eta_s_xi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_eta_s_eta(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    Ms_z_eta_s_eta = - s(xi,eta) / J(xi,eta) * ( ds/deta ) * ( dz/deta)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the SECOND TERM OF dsdz_0
!          in the SECOND TERM OF dzds_0

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_eta=-inv_jacob*dsdeta*dzdeta*mapping_anal(xil,etal,nodes_crd, &
                                                       1,ielem0)
end function Ms_z_eta_s_eta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_xi_s_eta(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    Ms_z_xi_s_eta = s(xi,eta) / J(xi,eta) * ( ds/deta ) * ( dz/xi)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the THIRD TERM OF dsdz_0
!          in the FIRST TERM OF dzds_0

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_xi_s_eta=inv_jacob*dsdeta*dzdxi*mapping_anal(xil,etal,nodes_crd,1,ielem0)

end function Ms_z_xi_s_eta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_xi_s_xi(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    Ms_z_xi_s_xi = - s(xi,eta) / J(xi,eta) * ( ds/dxi ) * ( dz/xi)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the FOURTH TERM OF dsdz_0
!          in the FOURTH TERM OF dzds_0

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_xi_s_xi = -inv_jacob*dsdxi*dzdxi*mapping_anal(xil,etal,nodes_crd,1,ielem0)
end function Ms_z_xi_s_xi
!-----------------------------------------------------------------------------------------


!*******************************************************************
!**********BEGIN axial part of M_* definitions**********************
!*******************************************************************
!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_eta_s_xi_k(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    Ms_z_eta_s_xi_k = 1 / J(xi,eta) * ( ds/dxi ) * ( dz/deta)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the FIRST TERM OF dsdz_0
!          in the THIRD TERM OF dzds_0
! It is defined within an element, and s(xi,eta) is
! defined by the analytic transformation .J is the determinant of
! the Jacobian matrix of the transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_xi_k = inv_jacob*dsdxi*dzdeta

end function Ms_z_eta_s_xi_k
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_eta_s_eta_k(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    Ms_z_eta_s_eta_k = - 1 / J(xi,eta) * ( ds/deta ) * ( dz/deta)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the SECOND TERM OF dsdz_0
!          in the SECOND TERM OF dzds_0

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_eta_k = -inv_jacob*dsdeta*dzdeta

end function Ms_z_eta_s_eta_k
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_xi_s_eta_k(xil, etal, nodes_crd, ielem0)
! This routines returns the value of
!
!    Ms_z_xi_s_eta_k = 1 / J(xi,eta) * ( ds/deta ) * ( dz/xi)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the THIRD TERM OF dsdz_0
!          in the FIRST TERM OF dzds_0

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_xi_s_eta_k = inv_jacob*dsdeta*dzdxi

end function Ms_z_xi_s_eta_k
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function Ms_z_xi_s_xi_k(xil,etal,nodes_crd,ielem0)
! This routines returns the value of
!
!    Ms_z_xi_s_xi = - 1 / J(xi,eta) * ( ds/dxi ) * ( dz/xi)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the FOURTH TERM OF dsdz_0
!          in the FOURTH TERM OF dzds_0

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_xi_s_xi_k = -inv_jacob*dsdxi*dzdxi

end function Ms_z_xi_s_xi_k
!-----------------------------------------------------------------------------------------
! NOTE: M_xi e.g. does not need to be specifically defined for the axis
!       since there is no factor s(), the axis case is therefore equal to
!       the non-axial case

!%********************************************************************
!% **********END axial part of M_* definitions************************
!%********************************************************************

!-----------------------------------------------------------------------------------------
subroutine mgrad_pointwise_anal(mg, xil, etal, nodes_crd, ielem0)
! This routines returns the following matrix:
!                      +                     +
!                      |(ds/dxi)  | (ds/deta)|
!    mg =  s(xi,eta) * | ---------|--------- |(xi,eta)
!                      |(dz/dxi ) | (dz/deta)|
!                      +                     +
!       This 2*2 matrix is needed when defining and storing
!gradient/divergence related arrays.

  integer :: ielem0
  real(kind=dp)    :: mg(2,2)
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)
  real(kind=dp)    :: dsdxi,dsdeta,dzdxi,dzdeta
  real(kind=dp)    :: sloc

  mg(:,:) = zero

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)

  sloc = mapping_anal(xil,etal,nodes_crd,1,ielem0)

  mg(1,1)  = sloc*dsdxi
  mg(1,2)  = sloc*dsdeta
  mg(2,1)  = sloc*dzdxi
  mg(2,2)  = sloc*dzdeta

end subroutine mgrad_pointwise_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine mgrad_pointwisek_anal(mg, xil, etal, nodes_crd, ielem0)
! This routines returns the following matrix:
!           +                     +
!           |(ds/dxi)  | (ds/deta)|
!    mg =   | ---------|--------- |(xi,eta)
!           |(dz/dxi ) | (dz/deta)|
!           +                     +
! This 2*2 matrix is needed when defining and storing
! gradient/divergence related arrays.

  integer,intent(in)           :: ielem0
  real(kind=dp)   ,intent(in)  :: xil, etal, nodes_crd(8,2)
  real(kind=dp)   ,intent(out) :: mg(2,2)
  real(kind=dp)                :: dsdxi,dsdeta,dzdxi,dzdeta

  mg(:,:) = zero

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,&
                                   nodes_crd,ielem0)
  mg(1,1) = dsdxi
  mg(1,2) = dsdeta
  mg(2,1) = dzdxi
  mg(2,2) = dzdeta

  end subroutine mgrad_pointwisek_anal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_partial_derivatives(dsdxi, dzdxi, dsdeta, dzdeta, xil, etal, &
                                         nodes_crd,ielem0)
! This routine returns the analytical values of the partial derivatives
! of the analytic spheroidal mapping.
!
  integer,intent(in)            :: ielem0
  real(kind=dp), intent(in)     :: xil, etal, nodes_crd(8,2)
  real(kind=dp), intent(out)    :: dsdxi, dzdxi, dsdeta, dzdeta

  if (eltype(ielem0) == 'curved') &
  call compute_partial_d_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,nodes_crd)
  if (eltype(ielem0) == 'linear') &
  call compute_partial_d_subpar(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,nodes_crd)
  if (eltype(ielem0) == 'semino') &
  call compute_partial_d_semino(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,nodes_crd)
  if (eltype(ielem0) == 'semiso') &
  call compute_partial_d_semiso(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,nodes_crd)

end subroutine compute_partial_derivatives
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_parameters(nodes_crd, a1, a2, b1, b2, deltatheta, thetabar)
  real(kind=dp)    :: nodes_crd(8,2)
  real(kind=dp)    :: a1,a2,b1,b2,deltatheta,thetabar,theta3,theta1
  real(kind=dp)    ::  s1,z1,s3,z3,s5,z5,s7,z7

  s1 = nodes_crd(1,1)
  z1 = nodes_crd(1,2)

  s3 = nodes_crd(3,1)
  z3 = nodes_crd(3,2)

  s5 = nodes_crd(5,1)
  z5 = nodes_crd(5,2)

  s7 = nodes_crd(7,1)
  z7 = nodes_crd(7,2)

  theta1 = datan(s1/(z1+epsi))
  theta3 = datan(s3/(z3+epsi))

  if ( zero > theta1 ) theta1 = pi + theta1
  if (theta1 == zero .and. z1 < 0) theta1 = pi
  if ( zero > theta3 ) theta3 = pi + theta3
  if (theta3 == zero .and. z3 < 0) theta3 = pi

  a1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(z3**2-z1**2))
  a2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(z5**2-z7**2))

  b1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(s1**2-s3**2))
  b2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(s7**2-s5**2))

  deltatheta = half*(theta3-theta1)
  thetabar   = half*(theta3+theta1)

end subroutine compute_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_parameters_new(nodes_crd, a1, a2, b1, b2, deltatheta1, &
                                  thetabar1, deltatheta2, thetabar2)

  real(kind=dp)    :: nodes_crd(8,2)
  real(kind=dp)    :: a1,a2,b1,b2,deltatheta1,thetabar1,deltatheta2,thetabar2
  real(kind=dp)    :: theta3,theta1,theta5,theta7
  real(kind=dp)    ::  s1,z1,s3,z3,s5,z5,s7,z7

  s1 = nodes_crd(1,1)
  z1 = nodes_crd(1,2)

  s3 = nodes_crd(3,1)
  z3 = nodes_crd(3,2)

  s5 = nodes_crd(5,1)
  z5 = nodes_crd(5,2)

  s7 = nodes_crd(7,1)
  z7 = nodes_crd(7,2)

  theta1 = datan(s1/(z1+epsi))
  theta3 = datan(s3/(z3+epsi))
  theta7 = datan(s7/(z7+epsi))
  theta5 = datan(s5/(z5+epsi))

  if ( zero > theta1 ) theta1 = pi + theta1
  if (theta1 == zero .and. z1 < 0) theta1 = pi
  if ( zero > theta3 ) theta3 = pi + theta3
  if (theta3 == zero .and. z3 < 0) theta3 = pi
  if ( zero > theta5 ) theta5 = pi + theta5
  if (theta5 == zero .and. z5 < 0) theta5 = pi
  if ( zero > theta7 ) theta7 = pi + theta7
  if (theta7 == zero .and. z7 < 0) theta7 = pi

  a1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(z3**2-z1**2))
  a2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(z5**2-z7**2))

  b1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(s1**2-s3**2))
  b2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(s7**2-s5**2))

  deltatheta1 = half*(theta3-theta1)
  thetabar1   = half*(theta3+theta1)
  deltatheta2 = half*(theta5-theta7)
  thetabar2   = half*(theta5+theta7)

end subroutine compute_parameters_new
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_parameters_srf(s1, s3, z1, z3, a, b, deltatheta, thetabar)

  real(kind=dp)   ,intent(out) :: a,b,deltatheta,thetabar
  real(kind=dp)    :: theta3,theta1
  real(kind=dp)   ,intent(in) ::  s1,z1,s3,z3

  a= zero
  b = zero
  deltatheta = zero
  thetabar = zero

  if (z1/=z3) a = dsqrt(dabs(((s1*z3)**2 - (s3*z1)**2)/(z3**2-z1**2)))

  if (s1/=s3) b = dsqrt(dabs(((s1*z3)**2 - (s3*z1)**2)/(s1**2-s3**2)))

  theta1 = datan(s1*b/(z1*a+epsi)) ; theta3 = datan(s3*b/(z3*a+epsi))

  if ( zero > theta1 ) theta1 = pi + theta1
  if (theta1 == zero .and. z1 < 0) theta1 = pi
  if ( zero > theta3 ) theta3 = pi + theta3
  if (theta3 == zero .and. z3 < 0) theta3 = pi

  deltatheta = half*(theta3-theta1)
  thetabar   = half*(theta3+theta1)

end subroutine compute_parameters_srf
!-----------------------------------------------------------------------------------------
!
!=============================
  end module analytic_mapping
!=============================
