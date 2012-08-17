!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

  subroutine get_jacobian_boundary_face(myrank,nspec, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
                        ispec,iface,jacobian2Dw_face,normal_face,NGLLA,NGLLB)

! returns jacobian2Dw_face and normal_face (pointing outwards of element)

  implicit none

  include "constants.h"

  integer nspec,myrank,nglob

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL) :: xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob)

! face information
  integer :: iface,ispec,NGLLA,NGLLB
  real(kind=CUSTOM_REAL) jacobian2Dw_face(NGLLA,NGLLB)
  real(kind=CUSTOM_REAL) normal_face(NDIM,NGLLA,NGLLB)

  double precision dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ)
  double precision dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ)
  double precision dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  double precision, dimension(NGLLX,NGLLY) :: wgllwgll_xy
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  double precision, dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! local parameters
! face corners
  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)

! check that the parameter file is correct
  if(NGNOD /= 8) call exit_MPI(myrank,'elements should have 8 control nodes')
  if(NGNOD2D /= 4) call exit_MPI(myrank,'surface elements should have 4 control nodes')

  select case ( iface )
  ! on reference face: xmin
  case(1)
    xelm(1)=xstore_dummy( ibool(1,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(1,NGLLY,1,ispec) )
    yelm(2)=ystore_dummy( ibool(1,NGLLY,1,ispec) )
    zelm(2)=zstore_dummy( ibool(1,NGLLY,1,ispec) )
    xelm(3)=xstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,1,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,1,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,1,NGLLZ,ispec) )

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_x,wgllwgll_yz, &
                  jacobian2Dw_face,normal_face,NGLLY,NGLLZ)

! on boundary: xmax
  case(2)
    xelm(1)=xstore_dummy( ibool(NGLLX,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(NGLLX,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(NGLLX,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_x,wgllwgll_yz, &
                  jacobian2Dw_face,normal_face,NGLLY,NGLLZ)

! on boundary: ymin
  case(3)
    xelm(1)=xstore_dummy( ibool(1,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,1,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,1,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,1,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,1,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,1,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,1,NGLLZ,ispec) )

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_y,wgllwgll_xz, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLZ)

! on boundary: ymax
  case(4)
    xelm(1)=xstore_dummy( ibool(1,NGLLY,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,NGLLY,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,NGLLY,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_y, wgllwgll_xz, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLZ)


! on boundary: bottom
  case(5)
    xelm(1)=xstore_dummy( ibool(1,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,1,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,1,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,1,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(4)=xstore_dummy( ibool(1,NGLLY,1,ispec) )
    yelm(4)=ystore_dummy( ibool(1,NGLLY,1,ispec) )
    zelm(4)=zstore_dummy( ibool(1,NGLLY,1,ispec) )

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm,&
                  dershape2D_bottom,wgllwgll_xy, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLY)

! on boundary: top
  case(6)
    xelm(1)=xstore_dummy( ibool(1,1,NGLLZ,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,NGLLZ,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,NGLLZ,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm,&
                  dershape2D_top, wgllwgll_xy, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLY)

  case default
    stop 'error 2D jacobian'
  end select

  end subroutine get_jacobian_boundary_face


! -------------------------------------------------------

  subroutine compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                                dershape2D,wgllwgll, &
                                jacobian2Dw_face,normal_face,NGLLA,NGLLB)

  implicit none

  include "constants.h"

! generic routine that accepts any polynomial degree in each direction
! returns 2D jacobian and normal for this face only

  integer NGLLA,NGLLB,myrank

  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)
  double precision wgllwgll(NGLLA,NGLLB)

  real(kind=CUSTOM_REAL) jacobian2Dw_face(NGLLA,NGLLB)
  real(kind=CUSTOM_REAL) normal_face(NDIM,NGLLA,NGLLB)

  integer i,j,ia
  double precision xxi,xeta,yxi,yeta,zxi,zeta
  double precision unx,uny,unz,jacobian

  do j=1,NGLLB
    do i=1,NGLLA

    xxi=ZERO
    xeta=ZERO
    yxi=ZERO
    yeta=ZERO
    zxi=ZERO
    zeta=ZERO
    do ia=1,NGNOD2D
      xxi=xxi+dershape2D(1,ia,i,j)*xelm(ia)
      xeta=xeta+dershape2D(2,ia,i,j)*xelm(ia)
      yxi=yxi+dershape2D(1,ia,i,j)*yelm(ia)
      yeta=yeta+dershape2D(2,ia,i,j)*yelm(ia)
      zxi=zxi+dershape2D(1,ia,i,j)*zelm(ia)
      zeta=zeta+dershape2D(2,ia,i,j)*zelm(ia)
    enddo

!   calculate the unnormalized normal to the boundary
    unx=yxi*zeta-yeta*zxi
    uny=zxi*xeta-zeta*xxi
    unz=xxi*yeta-xeta*yxi
    jacobian=dsqrt(unx**2+uny**2+unz**2)
    if(jacobian == ZERO) call exit_MPI(myrank,'2D Jacobian undefined')

!   normalize normal vector and store weighted surface jacobian

! distinguish if single or double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      jacobian2Dw_face(i,j) = sngl(jacobian * wgllwgll(i,j) )
      normal_face(1,i,j)=sngl(unx/jacobian)
      normal_face(2,i,j)=sngl(uny/jacobian)
      normal_face(3,i,j)=sngl(unz/jacobian)
    else
      jacobian2Dw_face(i,j) = jacobian * wgllwgll(i,j)
      normal_face(1,i,j)=unx/jacobian
      normal_face(2,i,j)=uny/jacobian
      normal_face(3,i,j)=unz/jacobian
    endif

    enddo
  enddo

  end subroutine compute_jacobian_2D_face


! This subroutine recompute the 3D jacobian for one element
! based upon 125 GLL points
! Hejun Zhu, Oct 16, 2009

! input: myrank,
!        xstore,ystore,zstore ----- input position
!        xigll,yigll,zigll ----- gll points position
!        ispec,nspec       ----- element number
!        ACTUALLY_STORE_ARRAYS   ------ save array or not

! output: xixstore,xiystore,xizstore,
!         etaxstore,etaystore,etazstore,
!         gammaxstore,gammaystore,gammazstore ------ parameters used for calculating jacobian


  subroutine recalc_jacobian_gll2D(myrank,xstore,ystore,zstore, &
                                  xigll,yigll,wgllwgll,NGLLA,NGLLB, &
                                  ispec,nspec,jacobian2Dw_face,normal_face)

  implicit none

  include "constants.h"

  ! input parameter
  integer::myrank,ispec,nspec
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer :: NGLLA,NGLLB
  double precision, dimension(NGLLA):: xigll
  double precision, dimension(NGLLB):: yigll
  double precision:: wgllwgll(NGLLA,NGLLB)

  real(kind=CUSTOM_REAL) jacobian2Dw_face(NGLLA,NGLLB)
  real(kind=CUSTOM_REAL) normal_face(NDIM,NGLLA,NGLLB)

  ! other parameters for this subroutine
  integer:: i,j,k,i1,j1,k1
  double precision:: xxi,xeta,yxi,yeta,zxi,zeta
  double precision:: xi,eta
  double precision,dimension(NGLLA):: hxir,hpxir
  double precision,dimension(NGLLB):: hetar,hpetar
  double precision:: hlagrange,hlagrange_xi,hlagrange_eta
  double precision:: jacobian
  double precision:: unx,uny,unz



  ! test parameters which can be deleted
  double precision:: xmesh,ymesh,zmesh
  double precision:: sumshape,sumdershapexi,sumdershapeeta

  ! first go over all gll points on face
  k=1
  do j=1,NGLLB
    do i=1,NGLLA

      xxi = 0.0
      xeta = 0.0
      yxi = 0.0
      yeta = 0.0
      zxi = 0.0
      zeta = 0.0

      xi = xigll(i)
      eta = yigll(j)

      ! calculate lagrange polynomial and its derivative
      call lagrange_any(xi,NGLLA,xigll,hxir,hpxir)
      call lagrange_any(eta,NGLLB,yigll,hetar,hpetar)

      ! test parameters
      sumshape = 0.0
      sumdershapexi = 0.0
      sumdershapeeta = 0.0
      xmesh = 0.0
      ymesh = 0.0
      zmesh = 0.0

      k1=1
      do j1 = 1,NGLLB
        do i1 = 1,NGLLA
         hlagrange = hxir(i1)*hetar(j1)
         hlagrange_xi = hpxir(i1)*hetar(j1)
         hlagrange_eta = hxir(i1)*hpetar(j1)


         xxi = xxi + xstore(i1,j1,k1,ispec)*hlagrange_xi
         xeta = xeta + xstore(i1,j1,k1,ispec)*hlagrange_eta

         yxi = yxi + ystore(i1,j1,k1,ispec)*hlagrange_xi
         yeta = yeta + ystore(i1,j1,k1,ispec)*hlagrange_eta

         zxi = zxi + zstore(i1,j1,k1,ispec)*hlagrange_xi
         zeta = zeta + zstore(i1,j1,k1,ispec)*hlagrange_eta

         ! test the lagrange polynomial and its derivate
         xmesh = xmesh + xstore(i1,j1,k1,ispec)*hlagrange
         ymesh = ymesh + ystore(i1,j1,k1,ispec)*hlagrange
         zmesh = zmesh + zstore(i1,j1,k1,ispec)*hlagrange
         sumshape = sumshape + hlagrange
         sumdershapexi = sumdershapexi + hlagrange_xi
         sumdershapeeta = sumdershapeeta + hlagrange_eta

         end do
      end do

      ! Check the lagrange polynomial and its derivative
      if (xmesh /=xstore(i,j,k,ispec).or.ymesh/=ystore(i,j,k,ispec).or.zmesh/=zstore(i,j,k,ispec)) then
        call exit_MPI(myrank,'new mesh positions are wrong in recalc_jacobian_gall3D.f90')
      end if
      if(abs(sumshape-one) >  TINYVAL) then
        call exit_MPI(myrank,'error shape functions in recalc_jacobian_gll3D.f90')
      end if
      if(abs(sumdershapexi) >  TINYVAL) then
        call exit_MPI(myrank,'error derivative xi shape functions in recalc_jacobian_gll3D.f90')
      end if
      if(abs(sumdershapeeta) >  TINYVAL) then
        call exit_MPI(myrank,'error derivative eta shape functions in recalc_jacobian_gll3D.f90')
      end if

!   calculate the unnormalized normal to the boundary
      unx=yxi*zeta-yeta*zxi
      uny=zxi*xeta-zeta*xxi
      unz=xxi*yeta-xeta*yxi
      jacobian=dsqrt(unx**2+uny**2+unz**2)
      if(jacobian <= ZERO) call exit_MPI(myrank,'2D Jacobian undefined')

!   normalize normal vector and store weighted surface jacobian

! distinguish if single or double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        jacobian2Dw_face(i,j) = sngl(jacobian * wgllwgll(i,j) )
        normal_face(1,i,j)=sngl(unx/jacobian)
        normal_face(2,i,j)=sngl(uny/jacobian)
        normal_face(3,i,j)=sngl(unz/jacobian)
      else
        jacobian2Dw_face(i,j) = jacobian * wgllwgll(i,j)
        normal_face(1,i,j)=unx/jacobian
        normal_face(2,i,j)=uny/jacobian
        normal_face(3,i,j)=unz/jacobian
      endif

    enddo
  enddo

  end subroutine recalc_jacobian_gll2D

