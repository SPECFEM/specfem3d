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

!----
!---- Stubs for parallel routines. Used by the serial version.
!----


  subroutine stop_all()
  stop 'error, program ended in exit_MPI'
  end subroutine stop_all

!
!----
!

  double precision function wtime()
  wtime = 0.d0
  end function wtime

!
!----
!

  subroutine sync_all()
  end subroutine sync_all

!
!----
!

  subroutine bcast_all_i(buffer, count)

  integer count
  integer, dimension(count) :: buffer

  end subroutine bcast_all_i

!
!----
!

  subroutine bcast_all_cr(buffer, count)

  include "constants.h"

  integer count
  real(kind=CUSTOM_REAL), dimension(count) :: buffer

  end subroutine bcast_all_cr

!
!----
!

  subroutine bcast_all_dp(buffer, count)

  integer count
  double precision, dimension(count) :: buffer

  end subroutine bcast_all_dp

!
!----
!

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

  integer sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_i

!
!----
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

  integer sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_dp

!
!----
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

  include "constants.h"

  integer sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_cr

!
!----
!

  subroutine gather_all_all_cr(sendbuf, recvbuf, counts,NPROC)

  implicit none

  include "constants.h"

  integer  NPROC,counts
  real(kind=CUSTOM_REAL), dimension(counts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(counts,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_all_cr

!
!----
!

 subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  implicit none

  include "constants.h"

  integer sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  recvbuf(:) = sendbuf(:)
  
  end subroutine gatherv_all_cr

!
!----
!


  subroutine init()
  end subroutine init

!
!----
!

  subroutine finalize()
  end subroutine finalize


!
!----
!

  subroutine world_size(size)

  implicit none

  integer size

  size = 1

  end subroutine world_size

!
!----
!

  subroutine world_rank(rank)

  implicit none

  integer rank

  rank = 0

  end subroutine world_rank

!
!----
!

  subroutine min_all_dp(sendbuf, recvbuf)

  implicit none

  double precision sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_dp

!
!----
!

  subroutine max_all_dp(sendbuf, recvbuf)

  implicit none

  double precision sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_dp

!
!----
!

  subroutine max_all_cr(sendbuf, recvbuf)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_cr

!
!----
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_all_cr
  
!
!----
!

  subroutine max_all_all_dp(sendbuf, recvbuf)

  implicit none

  double precision :: sendbuf, recvbuf

  recvbuf = sendbuf
  
  end subroutine max_all_all_dp
  

!
!----
!
!
!  subroutine min_all_all_dp(sendbuf, recvbuf)
!
!  implicit none
!
!  double precision :: sendbuf, recvbuf
!
!  recvbuf = sendbuf
!  
!  end subroutine min_all_all_dp
!
!----
!

  subroutine min_all_cr(sendbuf, recvbuf)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_cr

!
!----
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_all_cr

!
!----
!

  subroutine max_all_i(sendbuf, recvbuf)

  implicit none
  integer :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_i

!
!----
!

  subroutine min_all_i(sendbuf, recvbuf)

  implicit none
  integer:: sendbuf, recvbuf

  recvbuf = sendbuf
  
  end subroutine min_all_i

!
!----
!


  subroutine sum_all_dp(sendbuf, recvbuf)

  implicit none

  double precision sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_dp

!
!----
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_cr

!
!----
!

  subroutine sum_all_i(sendbuf, recvbuf)

  implicit none

  integer sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_i


!
!----
!
  subroutine sum_all_all_i(sendbuf, recvbuf)

  implicit none

  integer sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_all_i

!
!----
!

  subroutine any_all_l(sendbuf, recvbuf)

  implicit none

  logical sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine any_all_l

!
!----
!

  subroutine sendrecv_all_cr(sendbuf, sendcount, dest, sendtag, &
                             recvbuf, recvcount, source, recvtag)

  implicit none

  include "constants.h"

  integer sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  stop 'sendrecv_all_cr not implemented for serial code'

  end subroutine sendrecv_all_cr

!
!----
!

  integer function proc_null()
  proc_null = 0
  end function proc_null

!
!----
!

  subroutine issend_cr(sendbuf, sendcount, dest, sendtag, req)

  implicit none

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  
  stop 'issend_cr not implemented for serial code'

  end subroutine issend_cr

!
!----
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  implicit none


  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  stop 'irecv_cr not implemented for serial code'

  end subroutine irecv_cr

!
!----
!

  subroutine issend_i(sendbuf, sendcount, dest, sendtag, req)

  implicit none

  integer sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  stop 'issend_i not implemented for serial code'

  end subroutine issend_i

!
!----
!

  subroutine irecv_i(recvbuf, recvcount, dest, recvtag, req)

  implicit none

  integer recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf

  stop 'irecv_i not implemented for serial code'

  end subroutine irecv_i


!
!----
!

  subroutine recv_i(recvbuf, recvcount, dest, recvtag )

  implicit none
  
  !integer recvbuf,recvcount,dest,recvtag
  integer dest,recvtag
  integer recvcount
  integer,dimension(recvcount):: recvbuf
  
  stop 'recv_i not implemented for serial code'

  end subroutine recv_i

!
!----
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  implicit none
  
  integer recvcount,dest,recvtag
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf

  stop 'recvv_cr not implemented for serial code'

  end subroutine recvv_cr


!
!----
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  implicit none

  integer sendbuf,sendcount,dest,sendtag
  
  stop 'send_i not implemented for serial code'

  end subroutine send_i


!
!----
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  implicit none

  integer sendcount,dest,sendtag
  real(kind=CUSTOM_REAL),dimension(sendcount) :: sendbuf

  stop 'sendv_cr not implemented for serial code'

  end subroutine sendv_cr
!
!----
!

  subroutine wait_req(req)

  implicit none

  integer :: req

  end subroutine wait_req
  