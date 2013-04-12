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

  subroutine get_flags_boundaries(nspec,iproc_xi,iproc_eta,ispec,idoubling, &
             xstore,ystore,zstore,iboun,iMPIcut_xi,iMPIcut_eta, &
             NPROC_XI,NPROC_ETA, &
             UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer nspec
  integer ispec,idoubling
  integer NPROC_XI,NPROC_ETA

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK

  logical iboun(6,nspec)
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  double precision xstore(NGLLX_M,NGLLY_M,NGLLZ_M)
  double precision ystore(NGLLX_M,NGLLY_M,NGLLZ_M)
  double precision zstore(NGLLX_M,NGLLY_M,NGLLZ_M)

! use iproc_xi and iproc_eta to determine MPI cut planes along xi and eta
  integer iproc_xi,iproc_eta

  double precision target,sizeslice,TOLERANCE_METERS
  double precision xelm(8),yelm(8),zelm(8)

! find the coordinates of the eight corner nodes of the element
  xelm(1)=xstore(1,1,1)
  yelm(1)=ystore(1,1,1)
  zelm(1)=zstore(1,1,1)
  xelm(2)=xstore(NGLLX_M,1,1)
  yelm(2)=ystore(NGLLX_M,1,1)
  zelm(2)=zstore(NGLLX_M,1,1)
  xelm(3)=xstore(NGLLX_M,NGLLY_M,1)
  yelm(3)=ystore(NGLLX_M,NGLLY_M,1)
  zelm(3)=zstore(NGLLX_M,NGLLY_M,1)
  xelm(4)=xstore(1,NGLLY_M,1)
  yelm(4)=ystore(1,NGLLY_M,1)
  zelm(4)=zstore(1,NGLLY_M,1)
  xelm(5)=xstore(1,1,NGLLZ_M)
  yelm(5)=ystore(1,1,NGLLZ_M)
  zelm(5)=zstore(1,1,NGLLZ_M)
  xelm(6)=xstore(NGLLX_M,1,NGLLZ_M)
  yelm(6)=ystore(NGLLX_M,1,NGLLZ_M)
  zelm(6)=zstore(NGLLX_M,1,NGLLZ_M)
  xelm(7)=xstore(NGLLX_M,NGLLY_M,NGLLZ_M)
  yelm(7)=ystore(NGLLX_M,NGLLY_M,NGLLZ_M)
  zelm(7)=zstore(NGLLX_M,NGLLY_M,NGLLZ_M)
  xelm(8)=xstore(1,NGLLY_M,NGLLZ_M)
  yelm(8)=ystore(1,NGLLY_M,NGLLZ_M)
  zelm(8)=zstore(1,NGLLY_M,NGLLZ_M)

! compute geometrical tolerance small compared to size of model to detect edges
  TOLERANCE_METERS = dabs(UTM_X_MAX - UTM_X_MIN) / 100000.

! ****************************************************
!     determine if the element falls on a boundary
! ****************************************************

  iboun(:,ispec)=.false.

! on boundary 1: x=xmin
  target= UTM_X_MIN + TOLERANCE_METERS
  if(xelm(1)<target .and. xelm(4)<target .and. xelm(5)<target .and. xelm(8)<target) iboun(1,ispec)=.true.

! on boundary 2: xmax
  target= UTM_X_MAX - TOLERANCE_METERS
  if(xelm(2)>target .and. xelm(3)>target .and. xelm(6)>target .and. xelm(7)>target) iboun(2,ispec)=.true.

! on boundary 3: ymin
  target= UTM_Y_MIN + TOLERANCE_METERS
  if(yelm(1)<target .and. yelm(2)<target .and. yelm(5)<target .and. yelm(6)<target) iboun(3,ispec)=.true.

! on boundary 4: ymax
  target= UTM_Y_MAX - TOLERANCE_METERS
  if(yelm(3)>target .and. yelm(4)>target .and. yelm(7)>target .and. yelm(8)>target) iboun(4,ispec)=.true.

! on boundary 5: bottom
  target = Z_DEPTH_BLOCK + TOLERANCE_METERS
  if(zelm(1)<target .and. zelm(2)<target .and. zelm(3)<target .and. zelm(4)<target) iboun(5,ispec)=.true.

! on boundary 6: top
  if(idoubling == IFLAG_ONE_LAYER_TOPOGRAPHY) iboun(6,ispec)=.true.

! *******************************************************************
!     determine if the element falls on an MPI cut plane along xi
! *******************************************************************

! detect the MPI cut planes along xi in the cubed sphere

  iMPIcut_xi(:,ispec)=.false.

! angular size of a slice along xi
  sizeslice = (UTM_X_MAX-UTM_X_MIN) / NPROC_XI

! left cut-plane in the current slice along X = constant (Xmin of this slice)
! and add geometrical tolerance

  target = UTM_X_MIN + iproc_xi*sizeslice + TOLERANCE_METERS
  if(xelm(1)<target .and. xelm(4)<target .and. xelm(5)<target .and. xelm(8)<target) &
    iMPIcut_xi(1,ispec)=.true.

! right cut-plane in the current slice along X = constant (Xmax of this slice)
! and add geometrical tolerance

  target = UTM_X_MIN + (iproc_xi+1)*sizeslice - TOLERANCE_METERS
  if(xelm(2)>target .and. xelm(3)>target .and. xelm(6)>target .and. xelm(7)>target) &
    iMPIcut_xi(2,ispec)=.true.

! ********************************************************************
!     determine if the element falls on an MPI cut plane along eta
! ********************************************************************

  iMPIcut_eta(:,ispec)=.false.

! angular size of a slice along eta
  sizeslice = (UTM_Y_MAX-UTM_Y_MIN) / NPROC_ETA

! left cut-plane in the current slice along Y = constant (Ymin of this slice)
! and add geometrical tolerance

  target = UTM_Y_MIN + iproc_eta*sizeslice + TOLERANCE_METERS
  if(yelm(1)<target .and. yelm(2)<target .and. yelm(5)<target .and. yelm(6)<target) &
    iMPIcut_eta(1,ispec)=.true.

! right cut-plane in the current slice along Y = constant (Ymax of this slice)
! and add geometrical tolerance

  target = UTM_Y_MIN + (iproc_eta+1)*sizeslice - TOLERANCE_METERS
  if(yelm(3)>target .and. yelm(4)>target .and. yelm(7)>target .and. yelm(8)>target) &
    iMPIcut_eta(2,ispec)=.true.

  end subroutine get_flags_boundaries

