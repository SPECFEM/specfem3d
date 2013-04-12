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

  subroutine store_boundaries(myrank,iboun,nspec,&
                              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                              NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer nspec,myrank
  integer NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX

  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP)

  logical iboun(6,nspec)

  ! global element numbering
  integer ispec

  ! counters to keep track of number of elements on each of the boundaries
  integer ispecb1,ispecb2,ispecb3,ispecb4,ispecb5,ispecb6

  ! initializes
  ispecb1 = 0
  ispecb2 = 0
  ispecb3 = 0
  ispecb4 = 0
  ispecb5 = 0
  ispecb6 = 0

  ! determine if the element falls on a boundary
  do ispec=1,nspec

    ! on boundary: xmin
    if(iboun(1,ispec)) then
      ispecb1=ispecb1+1
      if( ispecb1 > NSPEC2DMAX_XMIN_XMAX ) stop 'error NSPEC2DMAX_XMIN_XMAX too small'
      ibelm_xmin(ispecb1)=ispec
    endif

    ! on boundary: xmax
    if(iboun(2,ispec)) then
      ispecb2=ispecb2+1
      if( ispecb2 > NSPEC2DMAX_XMIN_XMAX ) stop 'error NSPEC2DMAX_XMIN_XMAX too small'
      ibelm_xmax(ispecb2)=ispec
    endif

    ! on boundary: ymin
    if(iboun(3,ispec)) then
      ispecb3=ispecb3+1
      if( ispecb3 > NSPEC2DMAX_YMIN_YMAX ) stop 'error NSPEC2DMAX_YMIN_YMAX too small'
      ibelm_ymin(ispecb3)=ispec
    endif

    ! on boundary: ymax
    if(iboun(4,ispec)) then
      ispecb4=ispecb4+1
      if( ispecb4 > NSPEC2DMAX_YMIN_YMAX ) stop 'error NSPEC2DMAX_YMIN_YMAX too small'
      ibelm_ymax(ispecb4)=ispec
    endif

    ! on boundary: bottom
    if(iboun(5,ispec)) then
      ispecb5=ispecb5+1
      if( ispecb5 > NSPEC2D_BOTTOM ) stop 'error NSPEC2D_BOTTOM too small'
      ibelm_bottom(ispecb5)=ispec
    endif

    ! on boundary: top
    if(iboun(6,ispec)) then
      ispecb6=ispecb6+1
      if( ispecb6 > NSPEC2D_TOP ) stop 'error NSPEC2D_TOP too small'
      ibelm_top(ispecb6)=ispec
    endif

  enddo

  ! check theoretical value of elements at the bottom
  if(ispecb5 /= NSPEC2D_BOTTOM) call exit_MPI(myrank,'ispecb5 should equal NSPEC2D_BOTTOM')

  ! check theoretical value of elements at the top
  if(ispecb6 /= NSPEC2D_TOP) call exit_MPI(myrank,'ispecb6 should equal NSPEC2D_TOP')

  nspec2D_xmin = ispecb1
  nspec2D_xmax = ispecb2
  nspec2D_ymin = ispecb3
  nspec2D_ymax = ispecb4

  end subroutine store_boundaries

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
                            dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            jacobian2D_xmin,jacobian2D_xmax, &
                            jacobian2D_ymin,jacobian2D_ymax, &
                            jacobian2D_bottom,jacobian2D_top, &
                            normal_xmin,normal_xmax, &
                            normal_ymin,normal_ymax, &
                            normal_bottom,normal_top, &
                            NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer nspec,myrank
  integer NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX

  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP)

  logical iboun(6,nspec)

  double precision xstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision ystore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision zstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)

  real(kind=CUSTOM_REAL) jacobian2D_xmin(NGLLY_M,NGLLZ_M,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_xmax(NGLLY_M,NGLLZ_M,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymin(NGLLX_M,NGLLZ_M,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymax(NGLLX_M,NGLLZ_M,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_bottom(NGLLX_M,NGLLY_M,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) jacobian2D_top(NGLLX_M,NGLLY_M,NSPEC2D_TOP)

  real(kind=CUSTOM_REAL) normal_xmin(NDIM,NGLLY_M,NGLLZ_M,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_xmax(NDIM,NGLLY_M,NGLLZ_M,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_ymin(NDIM,NGLLX_M,NGLLZ_M,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_ymax(NDIM,NGLLX_M,NGLLZ_M,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_bottom(NDIM,NGLLX_M,NGLLY_M,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) normal_top(NDIM,NGLLX_M,NGLLY_M,NSPEC2D_TOP)

  double precision dershape2D_x(NDIM2D,NGNOD2D_FOUR_CORNERS,NGLLY_M,NGLLZ_M)
  double precision dershape2D_y(NDIM2D,NGNOD2D_FOUR_CORNERS,NGLLX_M,NGLLZ_M)
  double precision dershape2D_bottom(NDIM2D,NGNOD2D_FOUR_CORNERS,NGLLX_M,NGLLY_M)
  double precision dershape2D_top(NDIM2D,NGNOD2D_FOUR_CORNERS,NGLLX_M,NGLLY_M)

! global element numbering
  integer ispec

! counters to keep track of number of elements on each of the boundaries
  integer ispecb1,ispecb2,ispecb3,ispecb4,ispecb5,ispecb6

  double precision xelm(NGNOD2D_FOUR_CORNERS),yelm(NGNOD2D_FOUR_CORNERS),zelm(NGNOD2D_FOUR_CORNERS)

  ispecb1 = 0
  ispecb2 = 0
  ispecb3 = 0
  ispecb4 = 0
  ispecb5 = 0
  ispecb6 = 0

  do ispec=1,nspec

! determine if the element falls on a boundary

! on boundary: xmin

  if(iboun(1,ispec)) then

    ispecb1=ispecb1+1
    ibelm_xmin(ispecb1)=ispec

!   specify the 4 nodes for the 2-D boundary element
    xelm(1)=xstore(1,1,1,ispec)
    yelm(1)=ystore(1,1,1,ispec)
    zelm(1)=zstore(1,1,1,ispec)
    xelm(2)=xstore(1,NGLLY_M,1,ispec)
    yelm(2)=ystore(1,NGLLY_M,1,ispec)
    zelm(2)=zstore(1,NGLLY_M,1,ispec)
    xelm(3)=xstore(1,NGLLY_M,NGLLZ_M,ispec)
    yelm(3)=ystore(1,NGLLY_M,NGLLZ_M,ispec)
    zelm(3)=zstore(1,NGLLY_M,NGLLZ_M,ispec)
    xelm(4)=xstore(1,1,NGLLZ_M,ispec)
    yelm(4)=ystore(1,1,NGLLZ_M,ispec)
    zelm(4)=zstore(1,1,NGLLZ_M,ispec)

    call compute_jacobian_2D(myrank,ispecb1,xelm,yelm,zelm,dershape2D_x, &
                  jacobian2D_xmin,normal_xmin,NGLLY_M,NGLLZ_M,NSPEC2DMAX_XMIN_XMAX)

  endif

! on boundary: xmax

  if(iboun(2,ispec)) then

    ispecb2=ispecb2+1
    ibelm_xmax(ispecb2)=ispec

!   specify the 4 nodes for the 2-D boundary element
    xelm(1)=xstore(NGLLX_M,1,1,ispec)
    yelm(1)=ystore(NGLLX_M,1,1,ispec)
    zelm(1)=zstore(NGLLX_M,1,1,ispec)
    xelm(2)=xstore(NGLLX_M,NGLLY_M,1,ispec)
    yelm(2)=ystore(NGLLX_M,NGLLY_M,1,ispec)
    zelm(2)=zstore(NGLLX_M,NGLLY_M,1,ispec)
    xelm(3)=xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(3)=ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(3)=zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(4)=xstore(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(4)=ystore(NGLLX_M,1,NGLLZ_M,ispec)
    zelm(4)=zstore(NGLLX_M,1,NGLLZ_M,ispec)

    call compute_jacobian_2D(myrank,ispecb2,xelm,yelm,zelm,dershape2D_x, &
                  jacobian2D_xmax,normal_xmax,NGLLY_M,NGLLZ_M,NSPEC2DMAX_XMIN_XMAX)

  endif

! on boundary: ymin

  if(iboun(3,ispec)) then

    ispecb3=ispecb3+1
    ibelm_ymin(ispecb3)=ispec

!   specify the 4 nodes for the 2-D boundary element
    xelm(1)=xstore(1,1,1,ispec)
    yelm(1)=ystore(1,1,1,ispec)
    zelm(1)=zstore(1,1,1,ispec)
    xelm(2)=xstore(NGLLX_M,1,1,ispec)
    yelm(2)=ystore(NGLLX_M,1,1,ispec)
    zelm(2)=zstore(NGLLX_M,1,1,ispec)
    xelm(3)=xstore(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(3)=ystore(NGLLX_M,1,NGLLZ_M,ispec)
    zelm(3)=zstore(NGLLX_M,1,NGLLZ_M,ispec)
    xelm(4)=xstore(1,1,NGLLZ_M,ispec)
    yelm(4)=ystore(1,1,NGLLZ_M,ispec)
    zelm(4)=zstore(1,1,NGLLZ_M,ispec)

    call compute_jacobian_2D(myrank,ispecb3,xelm,yelm,zelm,dershape2D_y, &
                  jacobian2D_ymin,normal_ymin,NGLLX_M,NGLLZ_M,NSPEC2DMAX_YMIN_YMAX)

  endif

! on boundary: ymax

  if(iboun(4,ispec)) then

    ispecb4=ispecb4+1
    ibelm_ymax(ispecb4)=ispec

!   specify the 4 nodes for the 2-D boundary element
    xelm(1)=xstore(1,NGLLY_M,1,ispec)
    yelm(1)=ystore(1,NGLLY_M,1,ispec)
    zelm(1)=zstore(1,NGLLY_M,1,ispec)
    xelm(2)=xstore(NGLLX_M,NGLLY_M,1,ispec)
    yelm(2)=ystore(NGLLX_M,NGLLY_M,1,ispec)
    zelm(2)=zstore(NGLLX_M,NGLLY_M,1,ispec)
    xelm(3)=xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(3)=ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(3)=zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(4)=xstore(1,NGLLY_M,NGLLZ_M,ispec)
    yelm(4)=ystore(1,NGLLY_M,NGLLZ_M,ispec)
    zelm(4)=zstore(1,NGLLY_M,NGLLZ_M,ispec)

    call compute_jacobian_2D(myrank,ispecb4,xelm,yelm,zelm,dershape2D_y, &
                  jacobian2D_ymax,normal_ymax,NGLLX_M,NGLLZ_M,NSPEC2DMAX_YMIN_YMAX)

  endif

! on boundary: bottom

  if(iboun(5,ispec)) then

    ispecb5=ispecb5+1
    ibelm_bottom(ispecb5)=ispec

    xelm(1)=xstore(1,1,1,ispec)
    yelm(1)=ystore(1,1,1,ispec)
    zelm(1)=zstore(1,1,1,ispec)
    xelm(2)=xstore(NGLLX_M,1,1,ispec)
    yelm(2)=ystore(NGLLX_M,1,1,ispec)
    zelm(2)=zstore(NGLLX_M,1,1,ispec)
    xelm(3)=xstore(NGLLX_M,NGLLY_M,1,ispec)
    yelm(3)=ystore(NGLLX_M,NGLLY_M,1,ispec)
    zelm(3)=zstore(NGLLX_M,NGLLY_M,1,ispec)
    xelm(4)=xstore(1,NGLLY_M,1,ispec)
    yelm(4)=ystore(1,NGLLY_M,1,ispec)
    zelm(4)=zstore(1,NGLLY_M,1,ispec)

    call compute_jacobian_2D(myrank,ispecb5,xelm,yelm,zelm,dershape2D_bottom, &
                  jacobian2D_bottom,normal_bottom,NGLLX_M,NGLLY_M,NSPEC2D_BOTTOM)

  endif

! on boundary: top

  if(iboun(6,ispec)) then

    ispecb6=ispecb6+1
    ibelm_top(ispecb6)=ispec

    xelm(1)=xstore(1,1,NGLLZ_M,ispec)
    yelm(1)=ystore(1,1,NGLLZ_M,ispec)
    zelm(1)=zstore(1,1,NGLLZ_M,ispec)
    xelm(2)=xstore(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(2)=ystore(NGLLX_M,1,NGLLZ_M,ispec)
    zelm(2)=zstore(NGLLX_M,1,NGLLZ_M,ispec)
    xelm(3)=xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(3)=ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(3)=zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(4)=xstore(1,NGLLY_M,NGLLZ_M,ispec)
    yelm(4)=ystore(1,NGLLY_M,NGLLZ_M,ispec)
    zelm(4)=zstore(1,NGLLY_M,NGLLZ_M,ispec)

    call compute_jacobian_2D(myrank,ispecb6,xelm,yelm,zelm,dershape2D_top, &
                  jacobian2D_top,normal_top,NGLLX_M,NGLLY_M,NSPEC2D_TOP)

  endif

  enddo

! check theoretical value of elements at the bottom
  if(ispecb5 /= NSPEC2D_BOTTOM) call exit_MPI(myrank,'ispecb5 should equal NSPEC2D_BOTTOM')

! check theoretical value of elements at the top
  if(ispecb6 /= NSPEC2D_TOP) call exit_MPI(myrank,'ispecb6 should equal NSPEC2D_TOP')

  nspec2D_xmin = ispecb1
  nspec2D_xmax = ispecb2
  nspec2D_ymin = ispecb3
  nspec2D_ymax = ispecb4

  end subroutine get_jacobian_boundaries

! -------------------------------------------------------

  subroutine compute_jacobian_2D(myrank,ispecb,xelm,yelm,zelm,dershape2D,jacobian2D,normal,NGLLA,NGLLB,NSPEC2DMAX_AB)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

! generic routine that accepts any polynomial degree in each direction

  integer ispecb,NGLLA,NGLLB,NSPEC2DMAX_AB,myrank

  double precision xelm(NGNOD2D_FOUR_CORNERS),yelm(NGNOD2D_FOUR_CORNERS),zelm(NGNOD2D_FOUR_CORNERS)
  double precision dershape2D(NDIM2D,NGNOD2D_FOUR_CORNERS,NGLLA,NGLLB)

  real(kind=CUSTOM_REAL) jacobian2D(NGLLA,NGLLB,NSPEC2DMAX_AB)
  real(kind=CUSTOM_REAL) normal(3,NGLLA,NGLLB,NSPEC2DMAX_AB)

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
    do ia=1,NGNOD2D_FOUR_CORNERS
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

!   normalize normal vector and store surface jacobian

! distinguish if single or double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      jacobian2D(i,j,ispecb)=sngl(jacobian)
      normal(1,i,j,ispecb)=sngl(unx/jacobian)
      normal(2,i,j,ispecb)=sngl(uny/jacobian)
      normal(3,i,j,ispecb)=sngl(unz/jacobian)
    else
      jacobian2D(i,j,ispecb)=jacobian
      normal(1,i,j,ispecb)=unx/jacobian
      normal(2,i,j,ispecb)=uny/jacobian
      normal(3,i,j,ispecb)=unz/jacobian
    endif

    enddo
  enddo

  end subroutine compute_jacobian_2D

