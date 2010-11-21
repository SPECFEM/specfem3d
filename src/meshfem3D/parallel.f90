!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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
!---- Parallel routines.  All MPI calls belong in this file!
!----


  subroutine stop_all()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer ier

! stop all the MPI processes, and exit
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
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

  integer ier!,size

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

  ! print*,'init'
  !call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ier)
  !print*,'size = ',size

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

  subroutine maxloc_all_dp(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  double precision, dimension(2) :: sendbuf,recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION, &
                  MPI_MAXLOC,MPI_COMM_WORLD,ier)

  end subroutine maxloc_all_dp


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

subroutine send_i(sendbuf,sendcount,dest)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: dest,sendcount,ier
  integer :: tag = 100
  integer, dimension(sendcount) :: sendbuf

  ! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,tag, &
       MPI_COMM_WORLD,msg_status,ier)

end subroutine send_i

!
!----
!

subroutine recv_i(recvbuf,recvcount,source)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: source,recvcount,ier
  integer :: tag = 100
  integer, dimension(recvcount) :: recvbuf

  ! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,source,tag, &
       MPI_COMM_WORLD,msg_status,ier)

end subroutine recv_i


!
!----
!

subroutine send_dp(sendbuf,sendcount,dest)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: dest,sendcount,ier
  integer :: tag = 100
  double precision, dimension(sendcount) :: sendbuf

  ! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,tag, &
       MPI_COMM_WORLD,msg_status,ier)

end subroutine send_dp

!
!----
!

subroutine recv_dp(recvbuf,recvcount,source)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: recvcount,source,ier
  integer :: tag = 100
  double precision, dimension(recvcount) :: recvbuf

  ! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,source,tag, &
       MPI_COMM_WORLD,msg_status,ier)

end subroutine recv_dp


!
!----
!

  integer function proc_null()

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  proc_null = MPI_PROC_NULL

  end function proc_null
