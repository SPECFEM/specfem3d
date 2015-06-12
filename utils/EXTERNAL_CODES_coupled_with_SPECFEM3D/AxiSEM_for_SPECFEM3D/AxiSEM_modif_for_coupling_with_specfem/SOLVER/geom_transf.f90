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

!=====================
  module geom_transf
!=====================

  use data_mesh
  use subpar_mapping
  use analytic_mapping

  implicit none
  public :: jacobian

  public :: alpha,  beta,  gamma1, delta,  epsilon1, zeta
  public :: alphak, betak, gammak, deltak, epsilonk, zetak

  public :: jacobian_srf,  quadfunc_map, grad_quadfunc_map
  public :: mgrad_pointwise, mgrad_pointwisek
  public :: mapping, s_over_oneplusxi_axis

  logical, parameter :: ana_map = .true. ! We employ analytical mapping here.

  private
!
  contains
!//////////////////////////////////////////////////////////
!
!dk mapping----------------------------------------------------------
  real(kind=dp)    function mapping(xil,etal,nodes_crd,iaxis,ielem0)
!
  integer          :: iaxis,ielem0
  real(kind=dp)    :: xil,etal,nodes_crd(8,2)!,dumbdummy

  if (     ana_map) mapping = mapping_anal(xil,etal,nodes_crd,iaxis,ielem0)
  if (.not.ana_map) mapping = mapping_subpar(xil,etal,nodes_crd,iaxis)

!  if(ielem0==1 ) then
!     dumbdummy=mapping_anal(xi,eta,nodes_crd,iaxis,ielem0)
!     write(6,*)'IELGEOM:',ana_map,xi,eta,dumbdummy
!  endif

  if ( iaxis == 1 .and. dabs(mapping/router) < 1.d-12 ) mapping = 0.d0

  end function mapping
!
!--------------------------------------------------------------------
!
!dk quadfunc_map--------------------------------------------
  real(kind=dp)    function quadfunc_map(p,s,z,nodes_crd,ielem0)
!
!        This routines computes the
!quadratic functional (s-s(xi,eta))**2 + (z-z(xi,eta))**2
!
  integer :: ielem0
  real(kind=dp)    :: p(2), s,z, nodes_crd(8,2)
!
  if (     ana_map) quadfunc_map = quadfunc_map_anal(p,s,z,nodes_crd,ielem0)
  if (.not.ana_map) quadfunc_map = quadfunc_map_subpar(p,s,z,nodes_crd)

  end function quadfunc_map
!
!-----------------------------------------------------------------
!
!dk grad_quadfunc_map------------------------------------------
  subroutine grad_quadfunc_map(grd,p,s,z,nodes_crd,ielem0)
!
!       This routine returns the gradient of the quadratic
!functional associated with the mapping.
!
  integer :: ielem0
  real(kind=dp)    :: grd(2),p(2), s,z, nodes_crd(8,2)

  if (     ana_map) call grad_quadfunc_map_anal(grd,p,s,z,nodes_crd,ielem0)
  if (.not.ana_map) call grad_quadfunc_map_subpar(grd,p,s,z,nodes_crd)

  end subroutine grad_quadfunc_map
!
!--------------------------------------------------------------
!
!dk s_over_oneplusxi_axis--------------------------------------------
  real(kind=dp)    function s_over_oneplusxi_axis(xil,etal,nodes_crd,ielem0)
!
! This routine returns the value of the quantity
!
!              s/(1+xi)
!
! when the associated element lies along the axis of
! symmetry, in the case of an analytical transformation.

  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map)s_over_oneplusxi_axis=s_over_oneplusxi_axis_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map)s_over_oneplusxi_axis=s_over_oneplusxi_axis_subpar(xil,etal,nodes_crd)

  end function s_over_oneplusxi_axis
!
!--------------------------------------------------------------------


!=========================================================================
! -----TARJE-----------
!=========================================================================
!!$!dk one_over_oneplusxi_axis--------------------------------------------
!!$  real(kind=dp)    function one_over_oneplusxi_axis(xi,eta,nodes_crd,ielem0)
!!$!
!!$! This routine returns the value of the quantity
!!$!
!!$!              1/(1+xi)
!!$!
!!$! when the associated element lies along the axis of
!!$! symmetry, in the case of an analytical transformation.
!!$
!!$  integer :: ielem0
!!$  real(kind=dp)    :: xi, eta, nodes_crd(8,2)
!!$! WRONG WRONG WRONG: STILL HAVE TO DEFINE SUBPARAM FOR ONE_OVER_....!!!!!!!!!!!!!!!!!!!!!
!!$  if (.not.ana_map)one_over_oneplusxi_axis=one_over_oneplusxi_axis_subpar(xi,eta,nodes_crd)
!!$  if (     ana_map)one_over_oneplusxi_axis=one_over_oneplusxi_axis_anal(xi,eta,nodes_crd,ielem0)
!!$
!!$  end function one_over_oneplusxi_axis
!--------------------------------------------------------------------
!=========================================================================
! -----END TARJE-----------
!=========================================================================


!
!dk jacobian----------------------------------------------------
  real(kind=dp)    function jacobian(xil, etal, nodes_crd, ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) jacobian = jacobian_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) jacobian = jacobian_subpar(xil,etal,nodes_crd)

  end function jacobian
!----------------------------------------------------------------
!
!dk jacobian_srf------------------------------------------------------------
  real(kind=dp)    function jacobian_srf(xil,crdedge,ielem0)
!
!       This routine computes the Jacobian of the transformation
!that maps [-1,+1] into a portion of the boundary of domain.
!
  integer :: ielem0
  real(kind=dp)    :: xil, crdedge(3,2)

  if (     ana_map) then
     if (eltype(ielem0) /= 'linear') &
         jacobian_srf = jacobian_srf_anal(xil,crdedge)
     if (eltype(ielem0) == 'linear') &
         jacobian_srf = jacobian_srf_subpar(xil,crdedge)
  endif
  if (.not.ana_map) jacobian_srf = jacobian_srf_subpar(xil,crdedge)

!
  end function jacobian_srf
!---------------------------------------------------------------------------
!
!dk alphak------------------------------------------------------
  real(kind=dp)    function alphak(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) alphak = alphak_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) alphak = alphak_subpar(xil,etal,nodes_crd)

  end function alphak
!---------------------------------------------------------------
!
!dk betak------------------------------------------------------
  real(kind=dp)    function betak(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) betak = betak_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) betak = betak_subpar(xil,etal,nodes_crd)

  end function betak
!---------------------------------------------------------------
!
!dk gammak------------------------------------------------------
  real(kind=dp)    function gammak(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) gammak = gammak_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) gammak = gammak_subpar(xil,etal,nodes_crd)

  end function gammak
!---------------------------------------------------------------
!
!dk deltak------------------------------------------------------
  real(kind=dp)    function deltak(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) deltak = deltak_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) deltak = deltak_subpar(xil,etal,nodes_crd)

  end function deltak
!---------------------------------------------------------------
!
!dk epsilonk----------------------------------------------------
  real(kind=dp)    function epsilonk(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) epsilonk = epsilonk_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) epsilonk = epsilonk_subpar(xil,etal,nodes_crd)

  end function epsilonk
!---------------------------------------------------------------
!
!dk zetak------------------------------------------------------
  real(kind=dp)    function zetak(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) zetak = zetak_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) zetak = zetak_subpar(xil,etal,nodes_crd)

  end function zetak
!---------------------------------------------------------------
!
!dk alpha------------------------------------------------------
  real(kind=dp)    function alpha(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) alpha = alpha_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) alpha = alpha_subpar(xil,etal,nodes_crd)

  end function alpha
!---------------------------------------------------------------
!
!dk beta------------------------------------------------------
  real(kind=dp)    function beta(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) beta = beta_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) beta = beta_subpar(xil,etal,nodes_crd)

  end function beta
!---------------------------------------------------------------
!
!dk gamma1------------------------------------------------------
  real(kind=dp)    function gamma1(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) gamma1 = gamma_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) gamma1 = gamma_subpar(xil,etal,nodes_crd)

  end function gamma1
!---------------------------------------------------------------
!
!dk delta------------------------------------------------------
  real(kind=dp)    function delta(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) delta = delta_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) delta = delta_subpar(xil,etal,nodes_crd)

  end function delta
!---------------------------------------------------------------
!
!dk epsilon1----------------------------------------------------
  real(kind=dp)    function epsilon1(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if (     ana_map) epsilon1 = epsilon_anal(xil,etal,nodes_crd,ielem0)
  if (.not.ana_map) epsilon1 = epsilon_subpar(xil,etal,nodes_crd)

  end function epsilon1
!---------------------------------------------------------------
!
!dk zeta------------------------------------------------------
  real(kind=dp)    function zeta(xil,etal,nodes_crd,ielem0)
!
  integer :: ielem0
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

 if (     ana_map) zeta = zeta_anal(xil,etal,nodes_crd,ielem0)
 if (.not.ana_map) zeta = zeta_subpar(xil,etal,nodes_crd)

  end function zeta
!---------------------------------------------------------------
!
!dk mgrad_pointwise------------------------------------------------------
  subroutine mgrad_pointwise(mg,xil,etal,nodes_crd,ielem0)
!
! This routines returns the following matrix:
!                      +                     +
!                      |(ds/dxi)  | (ds/deta)|
!    mg =  s(xi,eta) * | ---------|--------- |(xi,eta)
!                      |(dz/dxi ) | (dz/deta)|
!                      +                     +
!       This 2*2 matrix is needed when defining and storing
!gradient/divergence related arrays.

  implicit none

  integer :: ielem0
  real(kind=dp)    :: mg(2,2)
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if(     ana_map) call mgrad_pointwise_anal(mg,xil,etal,nodes_crd,ielem0)
  if(.not.ana_map) call mgrad_pointwise_subpar(mg,xil,etal,nodes_crd)

  end subroutine mgrad_pointwise
!---------------------------------------------------------------------------
!
!dk mgrad_pointwisek------------------------------------------------------
  subroutine mgrad_pointwisek(mg,xil,etal,nodes_crd,ielem0)
!
! This routines returns the following matrix:
!          +                     +
!          |(ds/dxi)  | (ds/deta)|
!    mg =  | ---------|--------- |(xi,eta)
!          |(dz/dxi ) | (dz/deta)|
!          +                     +
!       This 2*2 matrix is needed when defining and storing
!gradient/divergence related arrays.

  implicit none

  integer :: ielem0
  real(kind=dp)    :: mg(2,2)
  real(kind=dp)    :: xil, etal, nodes_crd(8,2)

  if(     ana_map) call mgrad_pointwisek_anal(mg,xil,etal,nodes_crd,ielem0)
  if(.not.ana_map) call mgrad_pointwisek_subpar(mg,xil,etal,nodes_crd)

  end subroutine mgrad_pointwisek
!---------------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////
!
!=========================
  end module geom_transf
!=========================
