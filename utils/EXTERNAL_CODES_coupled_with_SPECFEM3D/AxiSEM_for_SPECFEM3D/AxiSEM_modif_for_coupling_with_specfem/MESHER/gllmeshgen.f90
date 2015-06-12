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

module gllmeshgen

  use data_mesh
  use analytic_mapping 
  use data_spec
  use data_gllmesh 
  use splib

  implicit none

  public :: test_mapping, create_gllmesh
  public :: extract_fluid_solid_submeshes
  private

  contains

!-----------------------------------------------------------------------------------------
subroutine create_gllmesh

  use data_diag

  real(kind=dp)    :: crd_nodes(8,2)

  integer           :: iel, jpol, ipol
  real(kind=dp)     :: stest

  ! Generate collocation points in the two directions of space
  allocate(sgll(0:npol,0:npol,neltot),zgll(0:npol,0:npol,neltot))

  ! QUADRATURE POINTS and weights
  allocate (eta(0:npol))
  allocate (dxi(0:npol))
  allocate (wt(0:npol))
  allocate(xi_k(0:npol), wt_axial_k(0:npol))

  call zemngl2(npol,xi_k)                       ! Gauss-Jacobi(0,1) quadrature
  call get_welegl_axial(npol,xi_k,wt_axial_k,2) !

  ! In the z-direction and in the s-direction for any other element

  call zelegl(npol,eta,dxi)                 ! Gauss-Lobatto Points
  call get_welegl(npol,eta,wt)              !

  !$omp parallel do shared(sgll, zgll) private(crd_nodes, jpol, ipol, stest)
  do iel = 1, neltot
  
     ! define dummy coordinate arrays
     crd_nodes(:,:) = 0.d0
     crd_nodes(1,1) = sg(lnodesg(1,iel))
     crd_nodes(1,2) = zg(lnodesg(1,iel))
     crd_nodes(3,1) = sg(lnodesg(2,iel))
     crd_nodes(3,2) = zg(lnodesg(2,iel))
     crd_nodes(5,1) = sg(lnodesg(3,iel))
     crd_nodes(5,2) = zg(lnodesg(3,iel))
     crd_nodes(7,1) = sg(lnodesg(4,iel))
     crd_nodes(7,2) = zg(lnodesg(4,iel))
     
     crd_nodes(2,:) = .5d0 * ( crd_nodes(1,:) + crd_nodes(3,:) )  ! midpoints are necessary
     crd_nodes(4,:) = .5d0 * ( crd_nodes(3,:) + crd_nodes(5,:) )  ! for subparametric mapping
     crd_nodes(6,:) = .5d0 * ( crd_nodes(5,:) + crd_nodes(7,:) )  ! (Serendipity elements).
     crd_nodes(8,:) = .5d0 * ( crd_nodes(7,:) + crd_nodes(1,:) )  !  

     stest = minval(sg(lnodesg(1:4,iel)))
     if ( stest < smallval_dble ) then 
        do jpol = 0, npol
           do ipol = 0, npol
              sgll(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,1,iel)
              zgll(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,2,iel)
           end do
        end do  
     else
        do jpol = 0, npol
           do ipol = 0, npol
              sgll(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,1,iel)
              zgll(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,2,iel)
           end do
        end do
     end if
  end do
  !$omp end parallel do
end subroutine create_gllmesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine extract_fluid_solid_submeshes
  integer :: iel_fluid,iel
  integer :: iel_solid
  
  allocate(sgll_fluid(0:npol,0:npol,neltot_fluid))
  sgll_fluid(:,:,:) = 0.d0 
  allocate(zgll_fluid(0:npol,0:npol,neltot_fluid))
  zgll_fluid(:,:,:) = 0.d0 
  
  do iel_fluid = 1, neltot_fluid
     iel = ielem_fluid(iel_fluid)
     sgll_fluid(:,:,iel_fluid) = sgll(:,:,iel)
     zgll_fluid(:,:,iel_fluid) = zgll(:,:,iel)
  end do
 
  allocate(sgll_solid(0:npol,0:npol,neltot_solid))
  sgll_solid(:,:,:) = 0.d0 
  allocate(zgll_solid(0:npol,0:npol,neltot_solid))
  zgll_solid(:,:,:) = 0.d0 
  do iel_solid = 1, neltot_solid
     iel = ielem_solid(iel_solid)
     sgll_solid(:,:,iel_solid) = sgll(:,:,iel)
     zgll_solid(:,:,iel_solid) = zgll(:,:,iel)
  end do
 
end subroutine extract_fluid_solid_submeshes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mapping
  real(kind=dp)     :: crd_nodes(8,2)
  integer           :: npoin
  integer           :: iel,jpol,ipol
  real(kind=dp)     ::  stest
  real(kind=dp)   , dimension(:,:,:), allocatable :: sglltmp,zglltmp

  write(6,*)'Test mapping...'; call flush(6)

  ! Generate collocation points in the  two directions of space
  npoin = neltot*(npol+1)**2
  allocate(sglltmp(0:npol,0:npol,neltot),zglltmp(0:npol,0:npol,neltot))

  ! QUADRATURE POINTS 
  
  call zemngl2(npol,xi_k)                       ! Gauss-Jacobi(0,1) quadrature
  call get_welegl_axial(npol,xi_k,wt_axial_k,2) !

  ! In the z-direction and in the s-direction for any other element

  call ZELEGL(npol,eta,dxi)                 ! Gauss-Lobatto Points
  call get_welegl(npol,eta,wt)              !
  
  do iel = 1, neltot

     ! define dummy coordinate arrays
     crd_nodes(:,:) = 0. 
     crd_nodes(1,1) = sg(lnodesg(1,iel)) ; crd_nodes(1,2) = zg(lnodesg(1,iel))
     crd_nodes(3,1) = sg(lnodesg(2,iel)) ; crd_nodes(3,2) = zg(lnodesg(2,iel))
     crd_nodes(5,1) = sg(lnodesg(3,iel)) ; crd_nodes(5,2) = zg(lnodesg(3,iel))
     crd_nodes(7,1) = sg(lnodesg(4,iel)) ; crd_nodes(7,2) = zg(lnodesg(4,iel))
     crd_nodes(2,:) = .5d0 * ( crd_nodes(1,:) + crd_nodes(3,:) )  ! midpoints are necessary
     crd_nodes(4,:) = .5d0 * ( crd_nodes(3,:) + crd_nodes(5,:) )  ! for subparametric mapping
     crd_nodes(6,:) = .5d0 * ( crd_nodes(5,:) + crd_nodes(7,:) )  ! (Serendipity elements).
     crd_nodes(8,:) = .5d0 * ( crd_nodes(7,:) + crd_nodes(1,:) )  !  

     stest = minval(sg(lnodesg(1:4,iel)))
     if ( stest < smallval_dble ) then 
        do jpol = 0, npol
           do ipol = 0, npol
              sglltmp(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,1,iel)
              zglltmp(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,2,iel)
           end do
        end do  
     else
        do jpol = 0, npol
           do ipol = 0, npol
              sglltmp(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,1,iel)
              zglltmp(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,2,iel)
           end do
        end do
     end if
  end do

  open(21,file='colloc_grid.dat')
  open(23,file='mesh.dat')
  do iel = 1,neltot
     do jpol = 0, npol
        do ipol = 0, npol
           write(21,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
        end do
        write(21,*)
     end do
                                                   
     do ipol = 0, npol
        do jpol = 0, npol
           write(21,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
        end do
        write(21,*)
     end do

     do jpol = 0,npol,npol
        do ipol = 0, npol
           write(23,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
        end do
        write(23,*)
     end do

     do ipol = 0,npol,npol
        do jpol = 0, npol
           write(23,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
        end do
        write(23,*)
     end do
  end do
  close(23)
  close(21)

  deallocate(sglltmp, zglltmp, eta, dxi, wt, xi_k, wt_axial_k)

end subroutine test_mapping
!-----------------------------------------------------------------------------------------

end module gllmeshgen
