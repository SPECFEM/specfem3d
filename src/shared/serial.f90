!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  implicit none
  real :: ct

  ! note: for simplicity, we take cpu_time which returns the elapsed CPU time in seconds
  !          (instead of wall clock time for parallel MPI function)
  call cpu_time(ct)

  wtime = ct

  end function wtime

!
!----
!

  subroutine synchronize_all()
  end subroutine synchronize_all

!
!----
!

  subroutine bcast_all_i(buffer, count)

  use unused_mod
  implicit none

  integer count
  integer, dimension(count) :: buffer

  unused_i4 = buffer(1)

  end subroutine bcast_all_i

!
!----
!

  subroutine bcast_all_cr(buffer, count)

  use unused_mod
  use constants

  implicit none

  integer count
  real(kind=CUSTOM_REAL), dimension(count) :: buffer

  unused_cr = buffer(1)

  end subroutine bcast_all_cr

!
!----
!

  subroutine bcast_all_dp(buffer, count)

  use unused_mod
  implicit none

  integer count
  double precision, dimension(count) :: buffer

  unused_dp = buffer(1)

  end subroutine bcast_all_dp

!
!----
!

  subroutine bcast_all_r(buffer, count)

  use unused_mod
  implicit none

  integer count
  real, dimension(count) :: buffer

  unused_r = buffer(1)

  end subroutine bcast_all_r


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

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  implicit none

  integer NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  recvbuf(0) = sendbuf

  end subroutine gather_all_singlei

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

  use constants

  implicit none

  integer sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_cr

!
!----
!

  subroutine gather_all_all_cr(sendbuf, recvbuf, counts,NPROC)

  use constants

  implicit none

  integer  NPROC,counts
  real(kind=CUSTOM_REAL), dimension(counts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(counts,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_all_cr

!
!----
!

 subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use unused_mod
  use constants

  implicit none

  integer sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  recvbuf(:) = sendbuf(:)

  unused_i4 = recvcount(1)
  unused_i4 = recvoffset(1)

  end subroutine gatherv_all_cr

!
!----
!

  subroutine init()

  use constants, only: NUMBER_OF_SIMULTANEOUS_RUNS

  if(NUMBER_OF_SIMULTANEOUS_RUNS <= 0) stop 'NUMBER_OF_SIMULTANEOUS_RUNS <= 0 makes no sense'

  if(NUMBER_OF_SIMULTANEOUS_RUNS > 1) stop 'serial runs require NUMBER_OF_SIMULTANEOUS_RUNS == 1'

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

  use constants

  implicit none

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_cr


!
!----
!

  subroutine min_all_cr(sendbuf, recvbuf)

  use constants

  implicit none

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_cr

!
!----
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  use constants

  implicit none

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_all_cr

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


  subroutine max_all_i(sendbuf, recvbuf)

  implicit none
  integer :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_i

!
!----
!

  subroutine max_allreduce_i(buffer,count)

  use unused_mod
  implicit none

  integer :: count
  integer,dimension(count),intent(inout) :: buffer

  unused_i4 = buffer(1)

  end subroutine max_allreduce_i

!
!----
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  use constants

  implicit none

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

  subroutine min_all_i(sendbuf, recvbuf)

  implicit none
  integer:: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_i

!
!----
!

  subroutine maxloc_all_dp(sendbuf, recvbuf)

  implicit none
  double precision, dimension(2) :: sendbuf,recvbuf

  recvbuf(1) = sendbuf(1)  ! maximum value
  recvbuf(2) = sendbuf(2)  ! index of myrank 0

  end subroutine maxloc_all_dp

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

  use constants

  implicit none

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_cr

!
!----
!

  subroutine sum_all_1Darray_dp(sendbuf, recvbuf, nx)

  implicit none

  integer :: nx
  double precision, dimension(nx) :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_1Darray_dp

!
!----
!

  subroutine sum_all_all_cr(sendbuf, recvbuf)

  use constants

  implicit none

  real(kind=CUSTOM_REAL) sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_all_cr

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

  use unused_mod
  use constants

  implicit none

  integer sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  stop 'sendrecv_all_cr not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_i4 = source
  unused_i4 = recvtag
  unused_cr = sendbuf(1)
  unused_cr = recvbuf(1)

  end subroutine sendrecv_all_cr

!
!----
!

  integer function proc_null()

  implicit none

  proc_null = 0

  end function proc_null

!
!----
!

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use unused_mod
  use constants

  implicit none

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  stop 'isend_cr not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_i4 = req
  unused_cr = sendbuf(1)

  end subroutine isend_cr

!
!----
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  use unused_mod
  use constants

  implicit none

  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  stop 'irecv_cr not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_i4 = req
  unused_cr = recvbuf(1)

  end subroutine irecv_cr

!
!----
!

  subroutine isend_i(sendbuf, sendcount, dest, sendtag, req)

  use unused_mod

  implicit none

  integer sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  stop 'isend_i not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_i4 = req
  unused_i4 = sendbuf(1)

  end subroutine isend_i

!
!----
!

  subroutine irecv_i(recvbuf, recvcount, dest, recvtag, req)

  use unused_mod

  implicit none

  integer recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf

  stop 'irecv_i not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_i4 = req
  unused_i4 = recvbuf(1)

  end subroutine irecv_i


!
!----
!

  subroutine recv_i(recvbuf, recvcount, dest, recvtag )

  use unused_mod

  implicit none

  integer dest,recvtag
  integer recvcount
  integer,dimension(recvcount):: recvbuf

  stop 'recv_i not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_i4 = recvbuf(1)

  end subroutine recv_i

!
!----
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  use unused_mod
  use constants

  implicit none

  integer recvcount,dest,recvtag
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf

  stop 'recvv_cr not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_cr = recvbuf(1)

  end subroutine recvv_cr


!
!----
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  use unused_mod

  implicit none

  integer dest,sendtag
  integer sendcount
  integer, dimension(sendcount) :: sendbuf

  stop 'send_i not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_i4 = sendbuf(1)

  end subroutine send_i


!
!----
!

  subroutine send_i_t(sendbuf,sendcount,dest)

  use unused_mod

  implicit none

  integer :: dest,sendcount
  integer, dimension(sendcount) :: sendbuf

  stop 'send_i_t not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendbuf(1)

  end subroutine send_i_t

!
!----
!


  subroutine recv_i_t(recvbuf,recvcount,source)

  use unused_mod

  implicit none

  integer :: source,recvcount
  integer, dimension(recvcount) :: recvbuf

  stop 'recv_i_t not implemented for serial code'
  unused_i4 = source
  unused_i4 = recvbuf(1)

  end subroutine recv_i_t


!
!----
!
!
!
!  subroutine send_dp_t(sendbuf,sendcount,dest)
!
!  implicit none
!
!  integer :: dest,sendcount
!  double precision, dimension(sendcount) :: sendbuf
!
!  stop 'send_dp_t not implemented for serial code'
!
!  end subroutine send_dp_t
!
!
!----
!

!  subroutine recv_dp_t(recvbuf,recvcount,source)
!
!  implicit none
!
!  integer :: recvcount,source
!  double precision, dimension(recvcount) :: recvbuf
!
!  stop 'recv_dp_t not implemented for serial code'
!
!  end subroutine recv_dp_t
!
!
!
!----
!  the following two subroutines are needed by locate_receivers.f90
  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  use unused_mod

  implicit none

  integer dest,sendtag
  integer sendcount
  double precision,dimension(sendcount):: sendbuf

  stop 'send_dp not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_dp = sendbuf(1)

  end subroutine send_dp
!
!----
!
  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  use unused_mod

  implicit none

  integer dest,recvtag
  integer recvcount
  double precision,dimension(recvcount):: recvbuf

  stop 'recv_dp not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_dp = recvbuf(1)

  end subroutine recv_dp

!
!----
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  use unused_mod

  use constants

  implicit none

  integer sendcount,dest,sendtag
  real(kind=CUSTOM_REAL),dimension(sendcount) :: sendbuf

  stop 'sendv_cr not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_cr = sendbuf(1)

  end subroutine sendv_cr
!
!----
!

  subroutine wait_req(req)

  use unused_mod

  implicit none

  integer :: req

  unused_i4 = req

  end subroutine wait_req

!
!----
!

  subroutine world_get_comm(comm)

  implicit none

  integer,intent(out) :: comm

  comm = 0

  end subroutine world_get_comm

!
!----
!

  subroutine world_duplicate(comm)

  implicit none

  integer,intent(out) :: comm

  comm = 0

  end subroutine world_duplicate

!
!----
!

  subroutine world_split()

  use constants

  implicit none

  mygroup = 0

  end subroutine world_split

!
!----
!

  subroutine world_unsplit()

  end subroutine world_unsplit

