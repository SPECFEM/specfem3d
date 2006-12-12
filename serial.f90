!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
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

  subroutine sum_all_dp(sendbuf, recvbuf)

  implicit none

  double precision sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_dp

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
