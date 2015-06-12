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

!> Variables concerned with elemental & spectral aspects only
!! (e.g. GLL points, Lagrange interpolant derivatives, quadrature weights)
!===================
 module data_spec
!===================

use global_parameters
implicit none
public 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  real(kind=dp), allocatable, dimension(:) :: xi_k, eta   ! Allocated in splib
  real(kind=dp), allocatable, dimension(:) :: dxi         ! " 
  real(kind=dp), allocatable, dimension(:) :: wt          !Quadrature weights
  real(kind=dp), allocatable, dimension(:) :: wt_axial_k  !Quad. wgts for the   
                                                     !gaus jacobi(0,1) integration

! Lagrange interpolant derivatives
! 3rd index: 1 - non-axial, 2 - axial
! 4th index: 1 - \partial_\xi, 2 - \partial_\eta
  real(kind=realkind), allocatable, dimension(:,:,:,:) :: shp_deri_k
  real(kind=realkind), allocatable, dimension(:,:)     :: G1, G1T, G2, G2T
  real(kind=realkind), allocatable, dimension(:)       :: G0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_spec
!=======================
