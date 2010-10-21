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

  subroutine bcast_all_cr(buffer, count)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  
  include "constants.h"  
  include "precision.h"

  integer count
  real(kind=CUSTOM_REAL), dimension(count) :: buffer

  integer ier

  call MPI_BCAST(buffer,count,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_cr

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

  subroutine gather_all_all_cr(sendbuf, recvbuf, counts, NPROC)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer NPROC,counts
  real(kind=CUSTOM_REAL), dimension(counts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(counts,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_ALLGATHER(sendbuf,counts,CUSTOM_MPI_TYPE,recvbuf,counts,CUSTOM_MPI_TYPE, &
                 MPI_COMM_WORLD,ier)

  end subroutine gather_all_all_cr

!
!----
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  integer ier

  call MPI_GATHERV(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,recvoffset,CUSTOM_MPI_TYPE, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gatherv_all_cr

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

  subroutine min_all_cr(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  include "constants.h"
  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_cr


!
!----
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  include "constants.h"
  include "precision.h"

  real(kind=CUSTOM_REAL):: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MIN,MPI_COMM_WORLD,ier)

  end subroutine min_all_all_cr

!
!----
!
!
!
!  subroutine min_all_all_dp(sendbuf, recvbuf)
!
!  implicit none
!
!! standard include of the MPI library
!  include 'mpif.h'
!  include "constants.h"
!  include "precision.h"
!
!  double precision :: sendbuf, recvbuf
!  integer ier
!
!  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
!                  MPI_MIN,MPI_COMM_WORLD,ier)
!
!  end subroutine min_all_all_dp
!
!
!----
!

  subroutine max_all_i(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_i

!
!----
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  include "constants.h"
  include "precision.h"

  real(kind=CUSTOM_REAL):: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MAX,MPI_COMM_WORLD,ier)

  end subroutine max_all_all_cr


!
!----
!


  subroutine max_all_all_dp(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  include "constants.h"
  include "precision.h"

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_MAX,MPI_COMM_WORLD,ier)

  end subroutine max_all_all_dp


!
!----
!

  subroutine min_all_i(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_i

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

  subroutine sum_all_cr(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  include "constants.h"
  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_cr

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

  subroutine sum_all_all_i(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_SUM,MPI_COMM_WORLD,ier)

  end subroutine sum_all_all_i

!
!----
!

  subroutine any_all_l(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  logical sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL, &
                  MPI_LOR,MPI_COMM_WORLD,ier)

  end subroutine any_all_l

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

!
!----
!

  subroutine issend_cr(sendbuf, sendcount, dest, sendtag, req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISSEND(sendbuf(1),sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine issend_cr

!
!----
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_IRECV(recvbuf(1),recvcount,CUSTOM_MPI_TYPE,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_cr

!
!----
!

  subroutine issend_i(sendbuf, sendcount, dest, sendtag, req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISSEND(sendbuf(1),sendcount,MPI_INTEGER,dest,sendtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine issend_i

!
!----
!

  subroutine irecv_i(recvbuf, recvcount, dest, recvtag, req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf
  integer ier

  call MPI_IRECV(recvbuf(1),recvcount,MPI_INTEGER,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_i


!
!----
!

  subroutine recv_i(recvbuf, recvcount, dest, recvtag )

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  
  integer dest,recvtag
  integer recvcount
  !integer recvbuf
  integer,dimension(recvcount):: recvbuf
  integer req(MPI_STATUS_SIZE)
  integer ier
  
  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,MPI_COMM_WORLD,req,ier)

  end subroutine recv_i

!
!----
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  
  include "constants.h"
  include "precision.h"
  
  integer recvcount,dest,recvtag
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf
  integer req(MPI_STATUS_SIZE)
  integer ier
  
  call MPI_RECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,MPI_COMM_WORLD,req,ier)


  end subroutine recvv_cr


!
!----
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  
  !integer sendbuf,sendcount,dest,sendtag
  integer dest,sendtag
  integer sendcount
  integer,dimension(sendcount):: sendbuf
  integer ier
  
  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_i


!
!----
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  implicit none

! standard include of the MPI library
  include 'mpif.h'
  
  include "constants.h"
  include "precision.h"

  integer sendcount,dest,sendtag
  real(kind=CUSTOM_REAL),dimension(sendcount) :: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine sendv_cr
!
!----
!

  subroutine wait_req(req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: req

  integer, dimension(MPI_STATUS_SIZE) :: req_mpi_status

  integer :: ier

  call mpi_wait(req,req_mpi_status,ier)

  end subroutine wait_req
