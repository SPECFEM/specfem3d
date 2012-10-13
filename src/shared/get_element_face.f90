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

  subroutine get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface_id )

! returns iface_id of face in reference element, determined by corner locations xcoord/ycoord/zcoord;

  implicit none

  include "constants.h"

  integer :: ispec,nspec,nglob,iface_id

! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord

! index array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! global point locations
  real(kind=CUSTOM_REAL) :: xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob)

! local parameters
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord_face,ycoord_face,zcoord_face
  real(kind=CUSTOM_REAL) :: midpoint_faces(NDIM,6),midpoint(NDIM),midpoint_distances(6)

! corners indices of reference cube faces
  ! shapes of arrays below
  integer,dimension(2),parameter :: face_shape = (/3,4/)
  integer,dimension(3),parameter :: all_faces_shape = (/3,4,6/)

  ! xmin
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),face_shape)
  ! xmax
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
       reshape((/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),face_shape)
  ! ymin
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
       reshape((/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),face_shape)
  ! ymax
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
       reshape((/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),face_shape)
  ! bottom
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),face_shape)
  ! top
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
       reshape((/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),face_shape)
  ! all faces
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
       reshape((/ iface1_corner_ijk,iface2_corner_ijk, &
                  iface3_corner_ijk,iface4_corner_ijk, &
                  iface5_corner_ijk,iface6_corner_ijk /),all_faces_shape)

! face orientation
  integer  :: ifa,icorner,i,j,k,iglob,iloc(1)

! initializes
  iface_id = -1

! gets face midpoint by its corners
  midpoint(:) = 0.0
  do icorner=1,NGNOD2D_FOUR_CORNERS
    midpoint(1) = midpoint(1) + xcoord(icorner)
    midpoint(2) = midpoint(2) + ycoord(icorner)
    midpoint(3) = midpoint(3) + zcoord(icorner)
  enddo
  midpoint(:) = midpoint(:) / 4.0

! determines element face by minimum distance of midpoints
  midpoint_faces(:,:) = 0.0
  do ifa=1,6

    ! face corners
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      i = iface_all_corner_ijk(1,icorner,ifa)
      j = iface_all_corner_ijk(2,icorner,ifa)
      k = iface_all_corner_ijk(3,icorner,ifa)

      ! coordinates
      iglob = ibool(i,j,k,ispec)
      xcoord_face(icorner) = xstore_dummy(iglob)
      ycoord_face(icorner) = ystore_dummy(iglob)
      zcoord_face(icorner) = zstore_dummy(iglob)

      ! face midpoint coordinates
      midpoint_faces(1,ifa) =  midpoint_faces(1,ifa) + xcoord_face(icorner)
      midpoint_faces(2,ifa) =  midpoint_faces(2,ifa) + ycoord_face(icorner)
      midpoint_faces(3,ifa) =  midpoint_faces(3,ifa) + zcoord_face(icorner)
    enddo

    midpoint_faces(:,ifa) = midpoint_faces(:,ifa) / 4.0

    ! distance
    midpoint_distances(ifa) = (midpoint(1)-midpoint_faces(1,ifa))**2 &
                            + (midpoint(2)-midpoint_faces(2,ifa))**2 &
                            + (midpoint(3)-midpoint_faces(3,ifa))**2
  enddo

! gets closest point, which determines face
  iloc = minloc(midpoint_distances)

  ! checks that found midpoint is close enough
  if( midpoint_distances(iloc(1)) > 1.e-4 * &
          (   (xcoord(1)-xcoord(2))**2 &
            + (ycoord(1)-ycoord(2))**2 &
            + (zcoord(1)-zcoord(2))**2 ) ) then
    print*,'error element face midpoint distance:',midpoint_distances(iloc(1)), &
          (xcoord(1)-xcoord(2))**2+(ycoord(1)-ycoord(2))**2+(zcoord(1)-zcoord(2))**2
    ! corner locations
    do icorner=1,NGNOD2D_FOUR_CORNERS
      i = iface_all_corner_ijk(1,icorner,iloc(1))
      j = iface_all_corner_ijk(2,icorner,iloc(1))
      k = iface_all_corner_ijk(3,icorner,iloc(1))
      iglob = ibool(i,j,k,ispec)
      print*,'error corner:',icorner,'xyz:',xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
    enddo
    ! stop
    print*,'xcoord:', xcoord(:)
    print*,'ycoord:', ycoord(:)
    print*,'zcoord:', zcoord(:)
    stop 'error element face midpoint'
  else
    iface_id = iloc(1)

  endif

  end subroutine get_element_face_id

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_element_face_gll_indices(iface,ijk_face,NGLLA,NGLLB )

! returns local indices in ijk_face for specified face

  implicit none

  include "constants.h"

  integer :: iface !,nspec,nglob

! gll point indices i,j,k for face, format corresponds to ijk_face(1,*) = i, ijk_face(2,*) = j, ijk_face(3,*) = k
  integer :: NGLLA,NGLLB
  integer,dimension(3,NGLLA,NGLLB) :: ijk_face

  integer :: i,j,k
  integer :: ngll,i_gll,j_gll,k_gll

! sets i,j,k indices of GLL points on boundary face
  ngll = 0
  select case( iface )

  ! reference xmin face
  case(1)
    if( NGLLA /= NGLLY .or. NGLLB /= NGLLZ ) stop 'error absorbing face 1 indexing'
    i_gll = 1
    do k=1,NGLLZ
      do j=1,NGLLY
        ngll = ngll + 1
        ijk_face(1,j,k) = i_gll
        ijk_face(2,j,k) = j
        ijk_face(3,j,k) = k
      enddo
    enddo

  ! reference xmax face
  case(2)
    if( NGLLA /= NGLLY .or. NGLLB /= NGLLZ ) stop 'error absorbing face 2 indexing'
    i_gll = NGLLX
    do k=1,NGLLZ
      do j=1,NGLLY
        ngll = ngll + 1
        ijk_face(1,j,k) = i_gll
        ijk_face(2,j,k) = j
        ijk_face(3,j,k) = k
      enddo
    enddo

  ! reference ymin face
  case(3)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLZ ) stop 'error absorbing face 3 indexing'
    j_gll = 1
    do k=1,NGLLZ
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,k) = i
        ijk_face(2,i,k) = j_gll
        ijk_face(3,i,k) = k
      enddo
    enddo

  ! reference ymax face
  case(4)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLZ ) stop 'error absorbing face 4 indexing'
    j_gll = NGLLY
    do k=1,NGLLZ
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,k) = i
        ijk_face(2,i,k) = j_gll
        ijk_face(3,i,k) = k
      enddo
    enddo

  ! reference bottom face
  case(5)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLY ) stop 'error absorbing face 5 indexing'
    k_gll = 1
    do j=1,NGLLY
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,j) = i
        ijk_face(2,i,j) = j
        ijk_face(3,i,j) = k_gll
      enddo
    enddo

  ! reference bottom face
  case(6)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLY ) stop 'error absorbing face 6 indexing'
    k_gll = NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,j) = i
        ijk_face(2,i,j) = j
        ijk_face(3,i,j) = k_gll
      enddo
    enddo

  case default
    stop 'error element face not found'

  end select

  ! checks number of gll points set on face
  if( ngll /= NGLLA*NGLLB ) then
    print*,'error element face ngll:',ngll,NGLLA,NGLLB
    stop 'error element face ngll'
  endif

end subroutine get_element_face_gll_indices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                ibool,nspec,nglob, &
                                xstore_dummy,ystore_dummy,zstore_dummy, &
                                normal)

! only changes direction of normal to point outwards of element

  implicit none

  include "constants.h"

  integer :: ispec,iface,nspec,nglob

! face corner locations
!! DK DK Oct 2012: in principle we should use NGNOD2D instead of NGNOD2D_FOUR_CORNERS when
!! DK DK Oct 2012: computing the normal in the case of HEX27 elements, to be more precise;
!! DK DK Oct 2012: but the code below would be a bit difficult to modify therefore we keep
!! DK DK Oct 2012: using NGNOD2D_FOUR_CORNERS only for now.
!! DK DK Oct 2012: If the face is flat (for instance along an absorbing edge) then the above
!! DK DK Oct 2012: code is still exact even for HEX27, but if there is bathymetry for instance
!! DK DK Oct 2012: along a curved fluid-solid interface then the code below is not as precise
!! DK DK Oct 2012: as it should for HEX27.
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord

! index array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! global point locations
  real(kind=CUSTOM_REAL),dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! face normal
  real(kind=CUSTOM_REAL),dimension(NDIM) :: normal

! local parameters
  real(kind=CUSTOM_REAL) :: face_n(NDIM),tmp,v_tmp(NDIM)
  integer :: iglob

! determines initial orientation given by three corners on the face
  ! cross-product of vectors from corner 1 to corner 2 and from corner 1 to corner 3
  face_n(1) =   (ycoord(2)-ycoord(1))*(zcoord(3)-zcoord(1)) - (zcoord(2)-zcoord(1))*(ycoord(3)-ycoord(1))
  face_n(2) = - (xcoord(2)-xcoord(1))*(zcoord(3)-zcoord(1)) + (zcoord(2)-zcoord(1))*(xcoord(3)-xcoord(1))
  face_n(3) =   (xcoord(2)-xcoord(1))*(ycoord(3)-ycoord(1)) - (ycoord(2)-ycoord(1))*(xcoord(3)-xcoord(1))
  tmp = sqrt( face_n(1)*face_n(1) + face_n(2)*face_n(2) + face_n(3)*face_n(3) )
  if( abs(tmp) < TINYVAL ) then
    print*,'error get face normal: length',tmp
    print*,'normal:',face_n(:)
    call exit_mpi(0,'error get element face normal')
  endif
  face_n(:) = face_n(:)/tmp

! checks that this normal direction is outwards of element:
  ! takes additional corner out of face plane and determines scalar product (dot product) to normal
  iglob = 0
  select case( iface )
  case(1) ! opposite to xmin face
    iglob = ibool(NGLLX,1,1,ispec)
  case(2) ! opposite to xmax face
    iglob = ibool(1,1,1,ispec)
  case(3) ! opposite to ymin face
    iglob = ibool(1,NGLLY,1,ispec)
  case(4) ! opposite to ymax face
    iglob = ibool(1,1,1,ispec)
  case(5) ! opposite to bottom
    iglob = ibool(1,1,NGLLZ,ispec)
  case(6) ! opposite to top
    iglob = ibool(1,1,1,ispec)
  case default
    call exit_mpi(0,'error get element face iface value')
  end select
  ! vector from corner 1 to this opposite one
  v_tmp(1) = xstore_dummy(iglob) - xcoord(1)
  v_tmp(2) = ystore_dummy(iglob) - ycoord(1)
  v_tmp(3) = zstore_dummy(iglob) - zcoord(1)

  ! scalar product (dot product)
  tmp = v_tmp(1)*face_n(1) + v_tmp(2)*face_n(2) + v_tmp(3)*face_n(3)

  ! makes sure normal points outwards, that is points away from this additional corner and scalar product (dot product) is negative
  if( tmp > 0.0_CUSTOM_REAL ) then
    face_n(:) = - face_n(:)
  endif

  ! in case given normal has zero length, sets it to computed face normal
  ! note: to avoid floating-point exception we use dble()
  !         values of normal(:) could be very small, almost zero, and lead to underflow
  tmp = sngl(dble(normal(1))**2 + dble(normal(2))**2 + dble(normal(3))**2)
  if( tmp < TINYVAL ) then
    normal(:) = face_n(:)
    return
  endif

  ! otherwise determines orientation of normal and flips direction such that normal points outwards
  tmp = face_n(1)*normal(1) + face_n(2)*normal(2) + face_n(3)*normal(3)
  if( tmp < 0.0_CUSTOM_REAL ) then
    !swap
    normal(:) = - normal(:)
  endif

  end subroutine get_element_face_normal

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_element_face_normal_idirect(ispec,iface,xcoord,ycoord,zcoord, &
                                ibool,nspec,nglob, &
                                xstore_dummy,ystore_dummy,zstore_dummy, &
                                normal,idirect)

! returns direction of normal:
!   idirect = 1 to point outwards of/away from element
!   idirect = 2 to point into element

  implicit none

  include "constants.h"

  integer :: ispec,iface,nspec,nglob

! face corner locations
!! DK DK Oct 2012: in principle we should use NGNOD2D instead of NGNOD2D_FOUR_CORNERS when
!! DK DK Oct 2012: computing the normal in the case of HEX27 elements, to be more precise;
!! DK DK Oct 2012: but the code below would be a bit difficult to modify therefore we keep
!! DK DK Oct 2012: using NGNOD2D_FOUR_CORNERS only for now.
!! DK DK Oct 2012: If the face is flat (for instance along an absorbing edge) then the above
!! DK DK Oct 2012: code is still exact even for HEX27, but if there is bathymetry for instance
!! DK DK Oct 2012: along a curved fluid-solid interface then the code below is not as precise
!! DK DK Oct 2012: as it should for HEX27.
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord

! index array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! global point locations
  real(kind=CUSTOM_REAL) :: xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob)

! face normal
  real(kind=CUSTOM_REAL),dimension(NDIM) :: normal

! direction type
  integer, intent(out) :: idirect

! local parameters
  real(kind=CUSTOM_REAL) :: face_n(3),tmp,v_tmp(3)
  integer :: iglob

  ! determines initial orientation given by three corners on the face
  ! cross-product of vectors from corner 1 to corner 2 and from corner 1 to corner 3
  face_n(1) =   (ycoord(2)-ycoord(1))*(zcoord(3)-zcoord(1)) - (zcoord(2)-zcoord(1))*(ycoord(3)-ycoord(1))
  face_n(2) = - (xcoord(2)-xcoord(1))*(zcoord(3)-zcoord(1)) + (zcoord(2)-zcoord(1))*(xcoord(3)-xcoord(1))
  face_n(3) =   (xcoord(2)-xcoord(1))*(ycoord(3)-ycoord(1)) - (ycoord(2)-ycoord(1))*(xcoord(3)-xcoord(1))
  tmp = sqrt( face_n(1)**2 + face_n(2)**2 + face_n(3)**2 )
  if( abs(tmp) < TINYVAL ) then
    print*,'error get face normal: length',tmp
    print*,'normal:',face_n(:)
    call exit_mpi(0,'error get element face normal')
  endif
  face_n(:) = face_n(:)/tmp

  ! checks that this normal direction is outwards of element:
  ! takes additional corner out of face plane and determines scalar product (dot product) to normal
  iglob = 0
  select case( iface )
  case(1) ! opposite to xmin face
    iglob = ibool(NGLLX,1,1,ispec)
  case(2) ! opposite to xmax face
    iglob = ibool(1,1,1,ispec)
  case(3) ! opposite to ymin face
    iglob = ibool(1,NGLLY,1,ispec)
  case(4) ! opposite to ymax face
    iglob = ibool(1,1,1,ispec)
  case(5) ! opposite to bottom
    iglob = ibool(1,1,NGLLZ,ispec)
  case(6) ! opposite to top
    iglob = ibool(1,1,1,ispec)
  case default
    call exit_mpi(0,'error get element face iface value')
  end select
  ! vector from corner 1 to this opposite one
  v_tmp(1) = xstore_dummy(iglob) - xcoord(1)
  v_tmp(2) = ystore_dummy(iglob) - ycoord(1)
  v_tmp(3) = zstore_dummy(iglob) - zcoord(1)

  ! scalar product (dot product)
  tmp = v_tmp(1)*face_n(1) + v_tmp(2)*face_n(2) + v_tmp(3)*face_n(3)

  ! makes sure normal points outwards, that is points away from this additional corner and scalar product (dot product) is negative
  if( tmp > 0.0 ) then
    face_n(:) = - face_n(:)
  endif

  ! in case given normal has zero length, exit
  if( ( normal(1)**2 + normal(2)**2 + normal(3)**2 ) < TINYVAL ) then
    print*,'problem: given normal is zero'
    idirect = 0
    return
  endif

! otherwise determines orientation of normal
  tmp = face_n(1)*normal(1) + face_n(2)*normal(2) + face_n(3)*normal(3)
  if( tmp < 0.0 ) then
    ! points into element
    idirect = 2
  else
    ! points away from element/ outwards
    idirect = 1
  endif

  end subroutine get_element_face_normal_idirect

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_element_corners(ispec,iface_ref,xcoord,ycoord,zcoord,iglob_corners_ref, &
                                ibool,nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
                                iface_all_corner_ijk)

  implicit none

  include "constants.h"

  integer,intent(in) :: ispec,iface_ref,nspec,nglob

  ! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS),intent(out) :: xcoord,ycoord,zcoord
  integer,dimension(NGNOD2D_FOUR_CORNERS),intent(out):: iglob_corners_ref

  ! index array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  ! global point locations
  real(kind=CUSTOM_REAL) :: xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob)

  ! assumes NGNOD2D_FOUR_CORNERS == 4
!! DK DK Oct 2012: in principle we should use NGNOD2D instead of NGNOD2D_FOUR_CORNERS when
!! DK DK Oct 2012: computing the normal in the case of HEX27 elements, to be more precise;
!! DK DK Oct 2012: but the code below would be a bit difficult to modify therefore we keep
!! DK DK Oct 2012: using NGNOD2D_FOUR_CORNERS only for now.
!! DK DK Oct 2012: If the face is flat (for instance along an absorbing edge) then the above
!! DK DK Oct 2012: code is still exact even for HEX27, but if there is bathymetry for instance
!! DK DK Oct 2012: along a curved fluid-solid interface then the code below is not as precise
!! DK DK Oct 2012: as it should for HEX27.
  integer,dimension(3,4,6) :: iface_all_corner_ijk

  ! local parameters
  integer :: icorner,i,j,k

  ! loops over corners
  do icorner = 1,NGNOD2D_FOUR_CORNERS
    i = iface_all_corner_ijk(1,icorner,iface_ref)
    j = iface_all_corner_ijk(2,icorner,iface_ref)
    k = iface_all_corner_ijk(3,icorner,iface_ref)

    ! global reference indices
    iglob_corners_ref(icorner) = ibool(i,j,k,ispec)

    ! reference corner coordinates
    xcoord(icorner) = xstore_dummy(iglob_corners_ref(icorner))
    ycoord(icorner) = ystore_dummy(iglob_corners_ref(icorner))
    zcoord(icorner) = zstore_dummy(iglob_corners_ref(icorner))
  enddo

  end subroutine get_element_corners

