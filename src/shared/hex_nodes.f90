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

  subroutine usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer,intent(in) :: NGNOD
  integer, dimension(NGNOD),intent(out) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! check that the parameter file is correct
  if (NGNOD /= 8 .and. NGNOD /= 27) &
    stop 'volume elements should have 8 or 27 control nodes'

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=2
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=2
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=2

  if (NGNOD == 27) then

    ! note: put further initialization into subroutine to avoid compilation errors
    !       in case NGNOD == 8
    call usual_hex_nodes_27(NGNOD,iaddx,iaddy,iaddz)

  endif

  end subroutine usual_hex_nodes

!
!-------------------------------------------------------------------------------------------------
!

  subroutine usual_hex_nodes_27(NGNOD,iaddx,iaddy,iaddz)

  implicit none

  integer,intent(in) :: NGNOD
  integer, dimension(NGNOD),intent(inout) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

  ! midside nodes (nodes located in the middle of an edge)

  iaddx(9) = 1
  iaddy(9) = 0
  iaddz(9) = 0

  iaddx(10) = 2
  iaddy(10) = 1
  iaddz(10) = 0

  iaddx(11) = 1
  iaddy(11) = 2
  iaddz(11) = 0

  iaddx(12) = 0
  iaddy(12) = 1
  iaddz(12) = 0

  iaddx(13) = 0
  iaddy(13) = 0
  iaddz(13) = 1

  iaddx(14) = 2
  iaddy(14) = 0
  iaddz(14) = 1

  iaddx(15) = 2
  iaddy(15) = 2
  iaddz(15) = 1

  iaddx(16) = 0
  iaddy(16) = 2
  iaddz(16) = 1

  iaddx(17) = 1
  iaddy(17) = 0
  iaddz(17) = 2

  iaddx(18) = 2
  iaddy(18) = 1
  iaddz(18) = 2

  iaddx(19) = 1
  iaddy(19) = 2
  iaddz(19) = 2

  iaddx(20) = 0
  iaddy(20) = 1
  iaddz(20) = 2

  ! side center nodes (nodes located in the middle of a face)

  iaddx(21) = 1
  iaddy(21) = 1
  iaddz(21) = 0

  iaddx(22) = 1
  iaddy(22) = 0
  iaddz(22) = 1

  iaddx(23) = 2
  iaddy(23) = 1
  iaddz(23) = 1

  iaddx(24) = 1
  iaddy(24) = 2
  iaddz(24) = 1

  iaddx(25) = 0
  iaddy(25) = 1
  iaddz(25) = 1

  iaddx(26) = 1
  iaddy(26) = 1
  iaddz(26) = 2

  ! center node (barycenter of the eight corners)

  iaddx(27) = 1
  iaddy(27) = 1
  iaddz(27) = 1

  end subroutine usual_hex_nodes_27

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes1(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=4
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=4
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=2
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=4
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=4
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=2
  iaddy(8)=4
  iaddz(8)=2

  end subroutine unusual_hex_nodes1

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes1p(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=4
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=4
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

  end subroutine unusual_hex_nodes1p

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes2(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=2

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=2

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=4

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=4

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=4

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=4

  end subroutine unusual_hex_nodes2

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes2p(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=-2

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=-2

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

  end subroutine unusual_hex_nodes2p

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes3(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

  end subroutine unusual_hex_nodes3

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes4(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=2

  end subroutine unusual_hex_nodes4

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes4p(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=2
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=2
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

  end subroutine unusual_hex_nodes4p

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes6(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=2
  iaddz(3)=-2

  iaddx(4)=0
  iaddy(4)=2
  iaddz(4)=-2

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=2

  end subroutine unusual_hex_nodes6

!
!-------------------------------------------------------------------------------------------------
!

  subroutine unusual_hex_nodes6p(iaddx,iaddy,iaddz)

  use constants

  implicit none

  integer, dimension(NGNOD_EIGHT_CORNERS) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=2
  iaddz(3)=2

  iaddx(4)=0
  iaddy(4)=2
  iaddz(4)=2

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=4

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=4

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=4

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=4

  end subroutine unusual_hex_nodes6p

!
!-------------------------------------------------------------------------------------------------
!


  subroutine hex_nodes_anchor_ijk_NGLL(NGNOD,anchor_iax,anchor_iay,anchor_iaz,NGLLX,NGLLY,NGLLZ)

! gets control point indices
!
! to get coordinates of control points (corners,midpoints) for an element ispec, they can be use as:
!do ia = 1,NGNOD
!  iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
!  xelm(ia) = dble(xstore(iglob))
!  yelm(ia) = dble(ystore(iglob))
!  zelm(ia) = dble(zstore(iglob))
!enddo

  implicit none

  integer, intent(in) :: NGNOD
  integer, dimension(NGNOD), intent(out) :: anchor_iax,anchor_iay,anchor_iaz
  integer, intent(in) :: NGLLX,NGLLY,NGLLZ

  ! local parameters
  integer :: ia
  integer :: iax,iay,iaz

  ! topology of the control/anchor points of the surface element
  integer :: iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

  ! define topology of the control element
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

  ! define (i,j,k) indices of the control/anchor points of the elements
  do ia = 1,NGNOD
    ! control point index
    iax = 0
    if (iaddx(ia) == 0) then
      iax = 1
    else if (iaddx(ia) == 1) then
      iax = (NGLLX+1)/2
    else if (iaddx(ia) == 2) then
      iax = NGLLX
    else
      stop 'incorrect value of iaddx'
    endif
    anchor_iax(ia) = iax

    iay = 0
    if (iaddy(ia) == 0) then
      iay = 1
    else if (iaddy(ia) == 1) then
      iay = (NGLLY+1)/2
    else if (iaddy(ia) == 2) then
      iay = NGLLY
    else
      stop 'incorrect value of iaddy'
    endif
    anchor_iay(ia) = iay

    iaz = 0
    if (iaddz(ia) == 0) then
      iaz = 1
    else if (iaddz(ia) == 1) then
      iaz = (NGLLZ+1)/2
    else if (iaddz(ia) == 2) then
      iaz = NGLLZ
    else
      stop 'incorrect value of iaddr'
    endif
    anchor_iaz(ia) = iaz
  enddo

  end subroutine hex_nodes_anchor_ijk_NGLL

