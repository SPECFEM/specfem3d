!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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
                                  UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,NEX_XI,NEX_ETA)

  use constants
  use constants_meshfem3D, only: NGLLX_M,NGLLY_M,NGLLZ_M,IFLAG_ONE_LAYER_TOPOGRAPHY

  implicit none

  integer,intent(in) :: nspec
  integer,intent(in) :: ispec,idoubling
  integer,intent(in) :: NPROC_XI,NPROC_ETA

  double precision,intent(in) :: UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK

  logical,intent(inout) :: iboun(6,nspec)
  logical,intent(inout) :: iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  double precision,intent(in) :: xstore(NGLLX_M,NGLLY_M,NGLLZ_M)
  double precision,intent(in) :: ystore(NGLLX_M,NGLLY_M,NGLLZ_M)
  double precision,intent(in) :: zstore(NGLLX_M,NGLLY_M,NGLLZ_M)

! use iproc_xi and iproc_eta to determine MPI cut planes along xi and eta
  integer,intent(in) :: iproc_xi,iproc_eta

  integer,intent(in) :: NEX_XI,NEX_ETA

  ! local parameters
  double precision :: target_val,sizeslice,TOLERANCE_METERS
  double precision :: dx,dy
  double precision :: xelm(NGNOD_EIGHT_CORNERS),yelm(NGNOD_EIGHT_CORNERS),zelm(NGNOD_EIGHT_CORNERS)

! find the coordinates of the eight corner nodes of the element
  xelm(1) = xstore(1,1,1)
  yelm(1) = ystore(1,1,1)
  zelm(1) = zstore(1,1,1)
  xelm(2) = xstore(NGLLX_M,1,1)
  yelm(2) = ystore(NGLLX_M,1,1)
  zelm(2) = zstore(NGLLX_M,1,1)
  xelm(3) = xstore(NGLLX_M,NGLLY_M,1)
  yelm(3) = ystore(NGLLX_M,NGLLY_M,1)
  zelm(3) = zstore(NGLLX_M,NGLLY_M,1)
  xelm(4) = xstore(1,NGLLY_M,1)
  yelm(4) = ystore(1,NGLLY_M,1)
  zelm(4) = zstore(1,NGLLY_M,1)
  xelm(5) = xstore(1,1,NGLLZ_M)
  yelm(5) = ystore(1,1,NGLLZ_M)
  zelm(5) = zstore(1,1,NGLLZ_M)
  xelm(6) = xstore(NGLLX_M,1,NGLLZ_M)
  yelm(6) = ystore(NGLLX_M,1,NGLLZ_M)
  zelm(6) = zstore(NGLLX_M,1,NGLLZ_M)
  xelm(7) = xstore(NGLLX_M,NGLLY_M,NGLLZ_M)
  yelm(7) = ystore(NGLLX_M,NGLLY_M,NGLLZ_M)
  zelm(7) = zstore(NGLLX_M,NGLLY_M,NGLLZ_M)
  xelm(8) = xstore(1,NGLLY_M,NGLLZ_M)
  yelm(8) = ystore(1,NGLLY_M,NGLLZ_M)
  zelm(8) = zstore(1,NGLLY_M,NGLLZ_M)

! compute geometrical tolerance small compared to size of model and size of element to detect edges
  dx = dabs(UTM_X_MAX - UTM_X_MIN)
  dy = dabs(UTM_Y_MAX - UTM_Y_MIN)
  TOLERANCE_METERS = min( dx / 100000. , dy / 100000. )

  ! in case of very large meshes
  if (NEX_XI * NGLLX_M > 100000) &
    TOLERANCE_METERS = min( TOLERANCE_METERS, dx / (NEX_XI * NGLLX_M) )
  if (NEX_ETA * NGLLY_M > 100000) &
    TOLERANCE_METERS = min( TOLERANCE_METERS, dy / (NEX_ETA * NGLLY_M) )

! ****************************************************
!     determine if the element falls on a boundary
! ****************************************************

  iboun(:,ispec)=.false.

! on boundary 1: x=xmin
  target_val = UTM_X_MIN + TOLERANCE_METERS
  if (xelm(1) < target_val .and. xelm(4) < target_val .and. xelm(5) < target_val .and. xelm(8) < target_val) iboun(1,ispec)=.true.

! on boundary 2: xmax
  target_val = UTM_X_MAX - TOLERANCE_METERS
  if (xelm(2) > target_val .and. xelm(3) > target_val .and. xelm(6) > target_val .and. xelm(7) > target_val) iboun(2,ispec)=.true.

! on boundary 3: ymin
  target_val = UTM_Y_MIN + TOLERANCE_METERS
  if (yelm(1) < target_val .and. yelm(2) < target_val .and. yelm(5) < target_val .and. yelm(6) < target_val) iboun(3,ispec)=.true.

! on boundary 4: ymax
  target_val = UTM_Y_MAX - TOLERANCE_METERS
  if (yelm(3) > target_val .and. yelm(4) > target_val .and. yelm(7) > target_val .and. yelm(8) > target_val) iboun(4,ispec)=.true.

! on boundary 5: bottom
  target_val = Z_DEPTH_BLOCK + TOLERANCE_METERS
  if (zelm(1) < target_val .and. zelm(2) < target_val .and. zelm(3) < target_val .and. zelm(4) < target_val) iboun(5,ispec)=.true.

! on boundary 6: top
  if (idoubling == IFLAG_ONE_LAYER_TOPOGRAPHY) iboun(6,ispec)=.true.

! *******************************************************************
!     determine if the element falls on an MPI cut plane along xi
! *******************************************************************

! detect the MPI cut planes along xi in the cubed sphere

  iMPIcut_xi(:,ispec)=.false.

! angular size of a slice along xi
  sizeslice = (UTM_X_MAX-UTM_X_MIN) / NPROC_XI

! left cut-plane in the current slice along X = constant (Xmin of this slice)
! and add geometrical tolerance

  target_val = UTM_X_MIN + iproc_xi*sizeslice + TOLERANCE_METERS
  if (xelm(1) < target_val .and. xelm(4) < target_val .and. xelm(5) < target_val .and. xelm(8) < target_val) &
    iMPIcut_xi(1,ispec)=.true.

! right cut-plane in the current slice along X = constant (Xmax of this slice)
! and add geometrical tolerance

  target_val = UTM_X_MIN + (iproc_xi+1)*sizeslice - TOLERANCE_METERS
  if (xelm(2) > target_val .and. xelm(3) > target_val .and. xelm(6) > target_val .and. xelm(7) > target_val) &
    iMPIcut_xi(2,ispec)=.true.

! ********************************************************************
!     determine if the element falls on an MPI cut plane along eta
! ********************************************************************

  iMPIcut_eta(:,ispec)=.false.

! angular size of a slice along eta
  sizeslice = (UTM_Y_MAX-UTM_Y_MIN) / NPROC_ETA

! left cut-plane in the current slice along Y = constant (Ymin of this slice)
! and add geometrical tolerance

  target_val = UTM_Y_MIN + iproc_eta*sizeslice + TOLERANCE_METERS
  if (yelm(1) < target_val .and. yelm(2) < target_val .and. yelm(5) < target_val .and. yelm(6) < target_val) &
    iMPIcut_eta(1,ispec)=.true.

! right cut-plane in the current slice along Y = constant (Ymax of this slice)
! and add geometrical tolerance

  target_val = UTM_Y_MIN + (iproc_eta+1)*sizeslice - TOLERANCE_METERS
  if (yelm(3) > target_val .and. yelm(4) > target_val .and. yelm(7) > target_val .and. yelm(8) > target_val) &
    iMPIcut_eta(2,ispec)=.true.

  end subroutine get_flags_boundaries

