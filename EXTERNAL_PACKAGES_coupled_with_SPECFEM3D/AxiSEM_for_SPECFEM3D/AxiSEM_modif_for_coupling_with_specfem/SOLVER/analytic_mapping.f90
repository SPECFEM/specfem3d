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
module analytic_mapping

  use global_parameters
  use analytic_spheroid_mapping
  use analytic_semi_mapping
  use subpar_mapping
  use data_mesh

  implicit none

  public :: jacobian
  public :: alpha,  beta,  gamma1,  delta,  epsilon1,  zeta
  public :: alphak, betak, gammak, deltak, epsilonk, zetak
  public :: mapping, s_over_oneplusxi_axis

  public :: Ms_z_eta_s_xi,   Ms_z_eta_s_eta
  public :: Ms_z_xi_s_eta,   Ms_z_xi_s_xi
  public :: Ms_z_eta_s_xi_k, Ms_z_eta_s_eta_k
  public :: Ms_z_xi_s_eta_k, Ms_z_xi_s_xi_k

  public :: compute_partial_derivatives

  private

contains

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function mapping(xil, etal, nodes_crd, iaxis, ielem0)
! < This routine computes the coordinates along the iaxis axis
!! of the image of any point in the reference domain in the physical domain
!! using the implicit assumption that the domain is spheroidal.

  integer, intent(in)          :: iaxis,ielem0
  real(kind=dp), intent(in)    :: xil, etal, nodes_crd(8,2)

  if (eltype(ielem0) == 'curved') &
     mapping = map_spheroid(xil,etal,nodes_crd,iaxis)
  if (eltype(ielem0) == 'linear') &
     mapping = mapping_subpar(xil,etal,nodes_crd,iaxis)
  if (eltype(ielem0) == 'semino') &
     mapping = map_semino(xil,etal,nodes_crd,iaxis)
  if (eltype(ielem0) == 'semiso') &
     mapping = map_semiso(xil,etal,nodes_crd,iaxis)

end function mapping
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function s_over_oneplusxi_axis(xil, etal, nodes_crd, ielem0)
! < This routine returns the value of the quantity
!!
!!              s/(1+xi)
!!
!! when the associated element lies along the axis of
!! symmetry, in the case of an analytical transformation.

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdxi, dsdeta, dzdeta

  if ( xil == -one ) then
     ! Apply L'Hopital's rule
     call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                      nodes_crd,ielem0)
     s_over_oneplusxi_axis = dsdxi
  else
     s_over_oneplusxi_axis = mapping(xil,etal,nodes_crd,1,ielem0) / &
                                  (one+xil)
  endif

end function s_over_oneplusxi_axis
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function jacobian(xil, etal, nodes_crd, ielem0)
! < This function returns the value of the jacobian of the
!! analytical mapping between the reference square [-1,1]^2 and
!! the deformed element in the spheroid.

  integer,intent(in)       :: ielem0
  real(kind=dp),intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)            :: dsdxi,dzdxi,dsdeta,dzdeta

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  jacobian = dsdxi*dzdeta-dsdeta*dzdxi

end function jacobian
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function alphak(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    alphak =  ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  alphak  = -inv_jacob*dsdxi*dsdeta

end function alphak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function betak(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    betak =  ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  betak  = inv_jacob*dsdxi**2

  end function betak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function gammak(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    gammak =  ( ds/deta ) * ( ds/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  gammak  = inv_jacob*dsdeta**2

end function gammak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function deltak(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    deltak = -( dz/dxi ) * ( dz/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  deltak  = -inv_jacob*dzdxi*dzdeta

end function deltak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function epsilonk(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    epsilonk = ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  epsilonk  = inv_jacob*dzdxi**2

end function epsilonk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function zetak(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    zetak = ( dz/deta ) * ( dz/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  zetak  = inv_jacob*dzdeta**2

end function zetak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function alpha(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    alpha = s(xi,eta) * ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  alpha  = -inv_jacob*dsdxi*dsdeta*mapping(xil,etal,nodes_crd,1,ielem0)

  end function alpha
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function beta(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    beta =  s(xi,eta) * ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  beta  = inv_jacob*dsdxi**2*mapping(xil,etal,nodes_crd,1,ielem0)

  end function beta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function gamma1(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    gamma = s(xi,eta) * ( ds/deta ) * ( ds/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  gamma1  = inv_jacob*dsdeta**2*mapping(xil,etal,nodes_crd,1,ielem0)

end function gamma1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function delta(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    delta = -s(xi,eta) * ( dz/dxi ) * ( dz/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  delta  = -inv_jacob*dzdxi*dzdeta*mapping(xil,etal,nodes_crd,1,ielem0)

end function delta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function epsilon1(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    epsilon = s(xi,eta) * ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  epsilon1  = inv_jacob*dzdxi**2*mapping(xil,etal,nodes_crd,1,ielem0)

end function epsilon1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function zeta(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    zeta = s(xi,eta) * ( dz/deta ) * ( dz/deta) / J(xi,eta),
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator. alpha is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.
!
  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  zeta  = inv_jacob*dzdeta**2*mapping(xil,etal,nodes_crd,1,ielem0)

  end function zeta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function Ms_z_eta_s_xi(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    Ms_z_eta_s_xi = s(xi,eta) / J(xi,eta) * ( ds/dxi ) * ( dz/deta)
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator in the FIRST TERM OF dsdz_0
!!          in the THIRD TERM OF dzds_0
!! It is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_xi=inv_jacob*dsdxi*dzdeta*mapping(xil,etal,nodes_crd,1,ielem0)

end function Ms_z_eta_s_xi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function Ms_z_eta_s_eta(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!    Ms_z_eta_s_eta = - s(xi,eta) / J(xi,eta) * ( ds/deta ) * ( dz/deta)
!
! a quantity that is needed in the calculation of the laplacian
! operator in the SECOND TERM OF dsdz_0
!          in the SECOND TERM OF dzds_0

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_eta=-inv_jacob*dsdeta*dzdeta*mapping(xil,etal,nodes_crd, &
                                                       1,ielem0)
end function Ms_z_eta_s_eta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function Ms_z_xi_s_eta(xil, etal, nodes_crd, ielem0)
!! This routines returns the value of
!!
!!    Ms_z_xi_s_eta = s(xi,eta) / J(xi,eta) * ( ds/deta ) * ( dz/xi)
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator in the THIRD TERM OF dsdz_0
!!          in the FIRST TERM OF dzds_0

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_xi_s_eta=inv_jacob*dsdeta*dzdxi*mapping(xil,etal,nodes_crd,1,ielem0)

end function Ms_z_xi_s_eta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function Ms_z_xi_s_xi(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    Ms_z_xi_s_xi = - s(xi,eta) / J(xi,eta) * ( ds/dxi ) * ( dz/xi)
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator in the FOURTH TERM OF dsdz_0
!!          in the FOURTH TERM OF dzds_0

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi, dzdeta, dzdxi, dsdeta, inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_xi_s_xi = -inv_jacob*dsdxi*dzdxi*mapping(xil,etal,nodes_crd,1,ielem0)
end function Ms_z_xi_s_xi
!-----------------------------------------------------------------------------------------


!*******************************************************************
!**********BEGIN axial part of M_* definitions**********************
!*******************************************************************
!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function Ms_z_eta_s_xi_k(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    Ms_z_eta_s_xi_k = 1 / J(xi,eta) * ( ds/dxi ) * ( dz/deta)
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator in the FIRST TERM OF dsdz_0
!!          in the THIRD TERM OF dzds_0
!! It is defined within an element, and s(xi,eta) is
!! defined by the analytic transformation .J is the determinant of
!! the Jacobian matrix of the transformation.

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal,nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_xi_k = inv_jacob*dsdxi*dzdeta

end function Ms_z_eta_s_xi_k
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function Ms_z_eta_s_eta_k(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    Ms_z_eta_s_eta_k = - 1 / J(xi,eta) * ( ds/deta ) * ( dz/deta)
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator in the SECOND TERM OF dsdz_0
!!          in the SECOND TERM OF dzds_0

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                   nodes_crd,ielem0)
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
  Ms_z_eta_s_eta_k = -inv_jacob*dsdeta*dzdeta

end function Ms_z_eta_s_eta_k
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function Ms_z_xi_s_eta_k(xil, etal, nodes_crd, ielem0)
! < This routines returns the value of
!!
!!    Ms_z_xi_s_eta_k = 1 / J(xi,eta) * ( ds/deta ) * ( dz/xi)
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator in the THIRD TERM OF dsdz_0
!!          in the FIRST TERM OF dzds_0

integer, intent(in)       :: ielem0
real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
real(kind=dp)             :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xil,etal, &
                                 nodes_crd,ielem0)
inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
Ms_z_xi_s_eta_k = inv_jacob*dsdeta*dzdxi

end function Ms_z_xi_s_eta_k
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function Ms_z_xi_s_xi_k(xil,etal,nodes_crd,ielem0)
! < This routines returns the value of
!!
!!    Ms_z_xi_s_xi = - 1 / J(xi,eta) * ( ds/dxi ) * ( dz/xi)
!!
!! a quantity that is needed in the calculation of the laplacian
!! operator in the FOURTH TERM OF dsdz_0
!!          in the FOURTH TERM OF dzds_0

  integer, intent(in)       :: ielem0
  real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)
  real(kind=dp)             :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob

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
pure subroutine compute_partial_derivatives(dsdxi, dzdxi, dsdeta, dzdeta, xil, etal, &
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
pure subroutine compute_parameters(nodes_crd, a1, a2, b1, b2, deltatheta, thetabar)
  real(kind=dp), intent(in)    :: nodes_crd(8,2)
  real(kind=dp), intent(out)   :: a1, a2, b1, b2, deltatheta, thetabar
  real(kind=dp)                :: theta3, theta1
  real(kind=dp)                :: s1, z1, s3, z3, s5, z5, s7, z7

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
pure subroutine compute_parameters_new(nodes_crd, a1, a2, b1, b2, deltatheta1, &
                                      thetabar1, deltatheta2, thetabar2)

  real(kind=dp), intent(in)  :: nodes_crd(8,2)
  real(kind=dp), intent(out) :: a1, a2, b1, b2, deltatheta1, thetabar1, deltatheta2, thetabar2
  real(kind=dp)              :: theta3, theta1, theta5, theta7
  real(kind=dp)              :: s1, z1, s3, z3, s5, z5, s7, z7

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
pure subroutine compute_parameters_srf(s1, s3, z1, z3, a, b, deltatheta, thetabar)

  real(kind=dp), intent(out) :: a, b, deltatheta, thetabar
  real(kind=dp)              :: theta3, theta1
  real(kind=dp), intent(in)  ::  s1, z1, s3, z3

  a= zero
  b = zero
  deltatheta = zero
  thetabar = zero

  if (z1 /= z3) a = dsqrt(dabs(((s1*z3)**2 - (s3*z1)**2)/(z3**2-z1**2)))

  if (s1 /= s3) b = dsqrt(dabs(((s1*z3)**2 - (s3*z1)**2)/(s1**2-s3**2)))

  theta1 = datan(s1*b/(z1*a+epsi)) ; theta3 = datan(s3*b/(z3*a+epsi))

  if ( zero > theta1 ) theta1 = pi + theta1
  if (theta1 == zero .and. z1 < 0) theta1 = pi
  if ( zero > theta3 ) theta3 = pi + theta3
  if (theta3 == zero .and. z3 < 0) theta3 = pi

  deltatheta = half*(theta3-theta1)
  thetabar   = half*(theta3+theta1)

end subroutine compute_parameters_srf
!-----------------------------------------------------------------------------------------

end module analytic_mapping
!=========================================================================================
