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
!---- Parallel routines.  All MPI calls belong in this file!
!----


  subroutine stop_all()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer ier

! stop all the MPI processes, and exit
! on some machines, MPI_FINALIZE needs to be called before MPI_ABORT
  call MPI_FINALIZE(ier)
  call MPI_ABORT(ier)
  stop 'error, program ended in exit_MPI'

  end subroutine stop_all

!
!----
!

  double precision function wtime()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  wtime = MPI_WTIME()

  end function wtime

!
!----
!

  subroutine sync_all()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer ier

  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  end subroutine sync_all

!
!----
!

  subroutine bcast_all_i(buffer, count)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer count
  integer, dimension(count) :: buffer

  integer ier

  call MPI_BCAST(buffer,count,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_i

!
!----
!

  subroutine bcast_all_dp(buffer, count)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer count
  double precision, dimension(count) :: buffer

  integer ier

  call MPI_BCAST(buffer,count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_dp

!
!----
!

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_i

!
!----
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_dp

!
!----
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,CUSTOM_MPI_TYPE, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_cr

!
!----
!

  subroutine init()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer ier

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

  end subroutine init

!
!----
!

  subroutine finalize()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer ier

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

  end subroutine finalize

!
!----
!

  subroutine world_size(size)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer size
  integer ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ier)

  end subroutine world_size

!
!----
!

  subroutine world_rank(rank)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer rank
  integer ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)

  end subroutine world_rank

!
!----
!

  subroutine min_all_dp(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_dp

!
!----
!

  subroutine max_all_dp(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_dp

!
!----
!

  subroutine max_all_cr(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_cr

!
!----
!

  subroutine sum_all_dp(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_dp

!
!----
!

  subroutine sum_all_i(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_i

!
!----
!

  subroutine sendrecv_all_cr(sendbuf, sendcount, dest, sendtag, &
                             recvbuf, recvcount, source, recvtag)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

  integer ier

  call MPI_SENDRECV(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                    recvbuf,recvcount,CUSTOM_MPI_TYPE,source,recvtag, &
                    MPI_COMM_WORLD,msg_status,ier)

  end subroutine sendrecv_all_cr

!
!----
!

  integer function proc_null()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  proc_null = MPI_PROC_NULL

  end function proc_null
