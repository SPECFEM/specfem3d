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

  subroutine usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer :: NGNOD
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! check that the parameter file is correct
  if(NGNOD /= 8 .and. NGNOD /= 27) &
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

  if(NGNOD == 27) then

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

  integer :: NGNOD
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

  implicit none

  include "constants.h"

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

