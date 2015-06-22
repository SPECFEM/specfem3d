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

  subroutine get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

  implicit none

  include "constants.h"

  integer myrank

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  integer i,j,k,ia

! location of the nodes of the 3D quadrilateral elements
  double precision xi,eta,gamma
  double precision l1xi,l2xi,l3xi,l1eta,l2eta,l3eta,l1gamma,l2gamma,l3gamma
  double precision l1pxi,l2pxi,l3pxi,l1peta,l2peta,l3peta,l1pgamma,l2pgamma,l3pgamma

! for checking the 3D shape functions
  double precision sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma

! check that the parameter file is correct
  if(NGNOD /= 27) call exit_MPI(myrank,'elements should have 27 control nodes')

! generate the 3D shape functions and their derivatives (27 nodes)
  do i=1,NGLLX

  xi=xigll(i)

  l1xi=HALF*xi*(xi-ONE)
  l2xi=ONE-xi**2
  l3xi=HALF*xi*(xi+ONE)

  l1pxi=xi-HALF
  l2pxi=-TWO*xi
  l3pxi=xi+HALF

  do j=1,NGLLY

    eta=yigll(j)

    l1eta=HALF*eta*(eta-ONE)
    l2eta=ONE-eta**2
    l3eta=HALF*eta*(eta+ONE)

    l1peta=eta-HALF
    l2peta=-TWO*eta
    l3peta=eta+HALF

    do k=1,NGLLZ

      gamma=zigll(k)

      l1gamma=HALF*gamma*(gamma-ONE)
      l2gamma=ONE-gamma**2
      l3gamma=HALF*gamma*(gamma+ONE)

      l1pgamma=gamma-HALF
      l2pgamma=-TWO*gamma
      l3pgamma=gamma+HALF

!     corner nodes

      shape3D(1,i,j,k)=l1xi*l1eta*l1gamma
      shape3D(2,i,j,k)=l3xi*l1eta*l1gamma
      shape3D(3,i,j,k)=l3xi*l3eta*l1gamma
      shape3D(4,i,j,k)=l1xi*l3eta*l1gamma
      shape3D(5,i,j,k)=l1xi*l1eta*l3gamma
      shape3D(6,i,j,k)=l3xi*l1eta*l3gamma
      shape3D(7,i,j,k)=l3xi*l3eta*l3gamma
      shape3D(8,i,j,k)=l1xi*l3eta*l3gamma

      dershape3D(1,1,i,j,k)=l1pxi*l1eta*l1gamma
      dershape3D(1,2,i,j,k)=l3pxi*l1eta*l1gamma
      dershape3D(1,3,i,j,k)=l3pxi*l3eta*l1gamma
      dershape3D(1,4,i,j,k)=l1pxi*l3eta*l1gamma
      dershape3D(1,5,i,j,k)=l1pxi*l1eta*l3gamma
      dershape3D(1,6,i,j,k)=l3pxi*l1eta*l3gamma
      dershape3D(1,7,i,j,k)=l3pxi*l3eta*l3gamma
      dershape3D(1,8,i,j,k)=l1pxi*l3eta*l3gamma

      dershape3D(2,1,i,j,k)=l1xi*l1peta*l1gamma
      dershape3D(2,2,i,j,k)=l3xi*l1peta*l1gamma
      dershape3D(2,3,i,j,k)=l3xi*l3peta*l1gamma
      dershape3D(2,4,i,j,k)=l1xi*l3peta*l1gamma
      dershape3D(2,5,i,j,k)=l1xi*l1peta*l3gamma
      dershape3D(2,6,i,j,k)=l3xi*l1peta*l3gamma
      dershape3D(2,7,i,j,k)=l3xi*l3peta*l3gamma
      dershape3D(2,8,i,j,k)=l1xi*l3peta*l3gamma

      dershape3D(3,1,i,j,k)=l1xi*l1eta*l1pgamma
      dershape3D(3,2,i,j,k)=l3xi*l1eta*l1pgamma
      dershape3D(3,3,i,j,k)=l3xi*l3eta*l1pgamma
      dershape3D(3,4,i,j,k)=l1xi*l3eta*l1pgamma
      dershape3D(3,5,i,j,k)=l1xi*l1eta*l3pgamma
      dershape3D(3,6,i,j,k)=l3xi*l1eta*l3pgamma
      dershape3D(3,7,i,j,k)=l3xi*l3eta*l3pgamma
      dershape3D(3,8,i,j,k)=l1xi*l3eta*l3pgamma

!     midside nodes

      shape3D(9,i,j,k)=l2xi*l1eta*l1gamma
      shape3D(10,i,j,k)=l3xi*l2eta*l1gamma
      shape3D(11,i,j,k)=l2xi*l3eta*l1gamma
      shape3D(12,i,j,k)=l1xi*l2eta*l1gamma
      shape3D(13,i,j,k)=l1xi*l1eta*l2gamma
      shape3D(14,i,j,k)=l3xi*l1eta*l2gamma
      shape3D(15,i,j,k)=l3xi*l3eta*l2gamma
      shape3D(16,i,j,k)=l1xi*l3eta*l2gamma
      shape3D(17,i,j,k)=l2xi*l1eta*l3gamma
      shape3D(18,i,j,k)=l3xi*l2eta*l3gamma
      shape3D(19,i,j,k)=l2xi*l3eta*l3gamma
      shape3D(20,i,j,k)=l1xi*l2eta*l3gamma

      dershape3D(1,9,i,j,k)=l2pxi*l1eta*l1gamma
      dershape3D(1,10,i,j,k)=l3pxi*l2eta*l1gamma
      dershape3D(1,11,i,j,k)=l2pxi*l3eta*l1gamma
      dershape3D(1,12,i,j,k)=l1pxi*l2eta*l1gamma
      dershape3D(1,13,i,j,k)=l1pxi*l1eta*l2gamma
      dershape3D(1,14,i,j,k)=l3pxi*l1eta*l2gamma
      dershape3D(1,15,i,j,k)=l3pxi*l3eta*l2gamma
      dershape3D(1,16,i,j,k)=l1pxi*l3eta*l2gamma
      dershape3D(1,17,i,j,k)=l2pxi*l1eta*l3gamma
      dershape3D(1,18,i,j,k)=l3pxi*l2eta*l3gamma
      dershape3D(1,19,i,j,k)=l2pxi*l3eta*l3gamma
      dershape3D(1,20,i,j,k)=l1pxi*l2eta*l3gamma

      dershape3D(2,9,i,j,k)=l2xi*l1peta*l1gamma
      dershape3D(2,10,i,j,k)=l3xi*l2peta*l1gamma
      dershape3D(2,11,i,j,k)=l2xi*l3peta*l1gamma
      dershape3D(2,12,i,j,k)=l1xi*l2peta*l1gamma
      dershape3D(2,13,i,j,k)=l1xi*l1peta*l2gamma
      dershape3D(2,14,i,j,k)=l3xi*l1peta*l2gamma
      dershape3D(2,15,i,j,k)=l3xi*l3peta*l2gamma
      dershape3D(2,16,i,j,k)=l1xi*l3peta*l2gamma
      dershape3D(2,17,i,j,k)=l2xi*l1peta*l3gamma
      dershape3D(2,18,i,j,k)=l3xi*l2peta*l3gamma
      dershape3D(2,19,i,j,k)=l2xi*l3peta*l3gamma
      dershape3D(2,20,i,j,k)=l1xi*l2peta*l3gamma

      dershape3D(3,9,i,j,k)=l2xi*l1eta*l1pgamma
      dershape3D(3,10,i,j,k)=l3xi*l2eta*l1pgamma
      dershape3D(3,11,i,j,k)=l2xi*l3eta*l1pgamma
      dershape3D(3,12,i,j,k)=l1xi*l2eta*l1pgamma
      dershape3D(3,13,i,j,k)=l1xi*l1eta*l2pgamma
      dershape3D(3,14,i,j,k)=l3xi*l1eta*l2pgamma
      dershape3D(3,15,i,j,k)=l3xi*l3eta*l2pgamma
      dershape3D(3,16,i,j,k)=l1xi*l3eta*l2pgamma
      dershape3D(3,17,i,j,k)=l2xi*l1eta*l3pgamma
      dershape3D(3,18,i,j,k)=l3xi*l2eta*l3pgamma
      dershape3D(3,19,i,j,k)=l2xi*l3eta*l3pgamma
      dershape3D(3,20,i,j,k)=l1xi*l2eta*l3pgamma

!     side center nodes

      shape3D(21,i,j,k)=l2xi*l2eta*l1gamma
      shape3D(22,i,j,k)=l2xi*l1eta*l2gamma
      shape3D(23,i,j,k)=l3xi*l2eta*l2gamma
      shape3D(24,i,j,k)=l2xi*l3eta*l2gamma
      shape3D(25,i,j,k)=l1xi*l2eta*l2gamma
      shape3D(26,i,j,k)=l2xi*l2eta*l3gamma

      dershape3D(1,21,i,j,k)=l2pxi*l2eta*l1gamma
      dershape3D(1,22,i,j,k)=l2pxi*l1eta*l2gamma
      dershape3D(1,23,i,j,k)=l3pxi*l2eta*l2gamma
      dershape3D(1,24,i,j,k)=l2pxi*l3eta*l2gamma
      dershape3D(1,25,i,j,k)=l1pxi*l2eta*l2gamma
      dershape3D(1,26,i,j,k)=l2pxi*l2eta*l3gamma

      dershape3D(2,21,i,j,k)=l2xi*l2peta*l1gamma
      dershape3D(2,22,i,j,k)=l2xi*l1peta*l2gamma
      dershape3D(2,23,i,j,k)=l3xi*l2peta*l2gamma
      dershape3D(2,24,i,j,k)=l2xi*l3peta*l2gamma
      dershape3D(2,25,i,j,k)=l1xi*l2peta*l2gamma
      dershape3D(2,26,i,j,k)=l2xi*l2peta*l3gamma

      dershape3D(3,21,i,j,k)=l2xi*l2eta*l1pgamma
      dershape3D(3,22,i,j,k)=l2xi*l1eta*l2pgamma
      dershape3D(3,23,i,j,k)=l3xi*l2eta*l2pgamma
      dershape3D(3,24,i,j,k)=l2xi*l3eta*l2pgamma
      dershape3D(3,25,i,j,k)=l1xi*l2eta*l2pgamma
      dershape3D(3,26,i,j,k)=l2xi*l2eta*l3pgamma

!     center node

      shape3D(27,i,j,k)=l2xi*l2eta*l2gamma

      dershape3D(1,27,i,j,k)=l2pxi*l2eta*l2gamma
      dershape3D(2,27,i,j,k)=l2xi*l2peta*l2gamma
      dershape3D(3,27,i,j,k)=l2xi*l2eta*l2pgamma

    enddo
  enddo
  enddo

! check the shape functions
  do i=1,NGLLX
    do j=1,NGLLY
      do k=1,NGLLZ

      sumshape=ZERO

      sumdershapexi=ZERO
      sumdershapeeta=ZERO
      sumdershapegamma=ZERO

      do ia=1,NGNOD

        sumshape=sumshape+shape3D(ia,i,j,k)

        sumdershapexi=sumdershapexi+dershape3D(1,ia,i,j,k)
        sumdershapeeta=sumdershapeeta+dershape3D(2,ia,i,j,k)
        sumdershapegamma=sumdershapegamma+dershape3D(3,ia,i,j,k)

      enddo

!     the sum of the shape functions should be 1
      if(abs(sumshape-ONE) > TINYVAL) call exit_MPI(myrank,'error in 3D shape functions')

!     the sum of the derivatives of the shape functions should be 0
      if(abs(sumdershapexi) > TINYVAL) &
        call exit_MPI(myrank,'error in xi derivatives of 3D shape function')

      if(abs(sumdershapeeta) > TINYVAL) &
        call exit_MPI(myrank,'error in eta derivatives of 3D shape function')

      if(abs(sumdershapegamma) > TINYVAL) &
        call exit_MPI(myrank,'error in gamma derivatives of 3D shape function')

      enddo
    enddo
  enddo

  end subroutine get_shape3D

