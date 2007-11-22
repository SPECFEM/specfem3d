!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

  subroutine usual_hex_nodes(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! check that the parameter file is correct
  if(NGNOD /= 8) stop 'elements should have 8 control nodes'

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

  end subroutine usual_hex_nodes

  subroutine unusual_hex_nodes1(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes1p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes2(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes2p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes3(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes4(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes4p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes6(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

  subroutine unusual_hex_nodes6p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

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

