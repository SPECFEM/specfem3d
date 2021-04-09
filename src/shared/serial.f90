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

!----
!---- Stubs for parallel routines. Used by the serial version.
!----

  subroutine abort_mpi()
  stop 'error, program ended in exit_MPI'
  end subroutine abort_mpi

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

  subroutine synchronize_all_comm(comm)

  implicit none

  integer,intent(in) :: comm

  ! local parameters
  integer :: idummy

  ! synchronizes MPI processes
  idummy = comm

  end subroutine synchronize_all_comm

!
!----
!

  subroutine bcast_all_i(buffer, countval)

  implicit none

  integer :: countval
  integer, dimension(countval) :: buffer
  integer(kind=4) :: unused_i4

  unused_i4 = buffer(1)

  end subroutine bcast_all_i

!
!----
!

  subroutine bcast_all_cr(buffer, countval)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: countval
  real(kind=CUSTOM_REAL), dimension(countval) :: buffer
  real(kind=CUSTOM_REAL) :: unused_cr

  unused_cr = buffer(1)

  end subroutine bcast_all_cr

!
!----
!

  subroutine bcast_all_singlecr(buffer)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: buffer
  real(kind=CUSTOM_REAL) :: unused_cr

  unused_cr = buffer

  end subroutine bcast_all_singlecr

!
!----
!

  subroutine bcast_all_dp(buffer, countval)

  implicit none

  integer :: countval
  double precision, dimension(countval) :: buffer
  double precision :: unused_dp

  unused_dp = buffer(1)

  end subroutine bcast_all_dp

!
!----
!

  subroutine bcast_all_singledp(buffer)

  implicit none

  double precision :: buffer
  double precision :: unused_dp

  unused_dp = buffer

  end subroutine bcast_all_singledp

!
!----
!

  subroutine bcast_all_r(buffer, countval)

  implicit none

  integer :: countval
  real, dimension(countval) :: buffer
  real :: unused_r

  unused_r = buffer(1)

  end subroutine bcast_all_r

!
!---
!

  subroutine bcast_all_ch_array(buffer,countval,STRING_LEN)

    implicit none

    integer :: countval,STRING_LEN

    character(len=STRING_LEN), dimension(countval) :: buffer
    character(len=STRING_LEN) :: unused_ch

    unused_ch=buffer(1)

  end subroutine bcast_all_ch_array

!
!----
!

  subroutine bcast_all_i_for_database(buffer, countval)

  implicit none

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  integer :: buffer
  integer(kind=4) :: unused_i4

  unused_i4 = countval

  unused_i4 = buffer

  end subroutine bcast_all_i_for_database

!
!----
!

  subroutine bcast_all_l_for_database(buffer, countval)

  implicit none

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  logical :: buffer
  integer(kind=4) :: unused_i4
  logical :: unused_l

  unused_i4 = countval

  unused_l = buffer

  end subroutine bcast_all_l_for_database

!
!----
!

  subroutine bcast_all_cr_for_database(buffer, countval)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real(kind=CUSTOM_REAL) :: buffer
  integer(kind=4) :: unused_i4
  real(kind=CUSTOM_REAL) :: unused_cr

  unused_i4 = countval

  unused_cr = buffer

  end subroutine bcast_all_cr_for_database

!
!----
!

  subroutine bcast_all_dp_for_database(buffer, countval)

  implicit none

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  double precision :: buffer
  integer(kind=4) :: unused_i4
  double precision :: unused_dp

  unused_i4 = countval

  unused_dp = buffer

  end subroutine bcast_all_dp_for_database

!
!----
!
! unused so far...
!
!  subroutine bcast_all_r_for_database(buffer, countval)
!
!  implicit none
!
!  integer countval
!  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
!  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
!  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
!  real :: buffer
!  integer(kind=4) :: unused_i4
!  real :: unused_r
!
!  unused_i4 = countval
!
!  unused_r = buffer
!
!  end subroutine bcast_all_r_for_database
!
!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

  implicit none

  integer :: buffer,idummy

  idummy = buffer

  end subroutine bcast_all_singlei

!
!----
!

  subroutine bcast_all_singlel(buffer)

  implicit none

  logical :: buffer,ldummy

  ldummy = buffer

  end subroutine bcast_all_singlel

!
!----
!

  subroutine bcast_all_string(buffer)

  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN) :: buffer,stringdummy

  stringdummy = buffer

  end subroutine bcast_all_string

!
!----
!

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_i

!
!----
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  recvbuf(0) = sendbuf

  end subroutine gather_all_singlei

!
!----
!

  subroutine gather_all_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_all_i

!
!----
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

  integer :: sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_dp

!
!----
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_cr

!
!----
!

  subroutine gather_all_all_cr(sendbuf, recvbuf, counts,NPROC)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: NPROC,counts
  real(kind=CUSTOM_REAL), dimension(counts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(counts,0:NPROC-1) :: recvbuf

  recvbuf(:,0) = sendbuf(:)

  end subroutine gather_all_all_cr

!
!----
!

  subroutine gatherv_all_i(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  implicit none

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcounttot) :: recvbuf

  integer(kind=4) :: unused_i4

  recvbuf(:) = sendbuf(:)

  unused_i4 = recvcount(1)
  unused_i4 = recvoffset(1)

  end subroutine gatherv_all_i

!
!----
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  integer(kind=4) :: unused_i4

  recvbuf(:) = sendbuf(:)

  unused_i4 = recvcount(1)
  unused_i4 = recvoffset(1)

  end subroutine gatherv_all_cr

!
!----
!

  subroutine all_gather_all_i(sendbuf, recvbuf, NPROC)

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(NPROC) :: recvbuf

  recvbuf(1) = sendbuf

  end subroutine all_gather_all_i

!
!----
!

  subroutine all_gather_all_r(sendbuf, sendcnt, recvbuf, recvcnt, recvoffset, dim1, NPROC)

  implicit none

  integer :: sendcnt, dim1, NPROC

  real, dimension(sendcnt) :: sendbuf
  real, dimension(dim1, NPROC) :: recvbuf

  integer, dimension(NPROC) :: recvoffset, recvcnt

  integer(kind=4) :: unused_i4

  recvbuf(1:sendcnt,1) = sendbuf(:)

  unused_i4 = recvcnt(1)
  unused_i4 = recvoffset(1)

  end subroutine all_gather_all_r

!
!----
!

  subroutine all_gather_all_ch(sendbuf, sendcnt, recvbuf, recvcnt, recvoffset, dim1, dim2, NPROC)

  implicit none

  integer :: sendcnt, dim1, dim2, NPROC

  character(len=dim2), dimension(sendcnt) :: sendbuf
  character(len=dim2), dimension(dim1, NPROC) :: recvbuf

  integer, dimension(NPROC) :: recvoffset, recvcnt

  integer(kind=4) :: unused_i4

  recvbuf(1:sendcnt,1) = sendbuf(:)

  unused_i4 = recvcnt(1)
  unused_i4 = recvoffset(1)

  end subroutine all_gather_all_ch


!
!----
!

  subroutine init_mpi()

  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  integer :: myrank
  logical :: BROADCAST_AFTER_READ

  ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS is read, thus read the parameter file
  myrank = 0
  BROADCAST_AFTER_READ = .false.
  call read_parameter_file(BROADCAST_AFTER_READ)

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 0) stop 'NUMBER_OF_SIMULTANEOUS_RUNS <= 0 makes no sense'

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) stop 'serial runs require NUMBER_OF_SIMULTANEOUS_RUNS == 1'

  end subroutine init_mpi

!
!----
!

  subroutine finalize_mpi()
  end subroutine finalize_mpi

!
!----
!

  subroutine world_size(size)

  implicit none

  integer :: size

  size = 1

  end subroutine world_size

!
!----
!

  subroutine world_size_comm(sizeval,comm)

  implicit none

  integer,intent(out) :: sizeval
  integer,intent(in) :: comm

  ! local parameters
  integer :: idummy

  sizeval = 1
  idummy = comm

  end subroutine world_size_comm

!
!----
!

  subroutine world_rank(rank)

  implicit none

  integer :: rank

  rank = 0

  end subroutine world_rank


!
!----
!

  subroutine world_rank_comm(rank,comm)

  implicit none

  integer,intent(out) :: rank
  integer,intent(in) :: comm

  ! local parameters
  integer :: idummy

  rank = 0
  idummy = comm

  end subroutine world_rank_comm

!
!----
!

  subroutine min_all_dp(sendbuf, recvbuf)

  implicit none

  double precision :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_dp

!
!----
!

  subroutine max_all_dp(sendbuf, recvbuf)

  implicit none

  double precision :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_dp

!
!----
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine max_all_cr


!
!----
!

  subroutine min_all_cr(sendbuf, recvbuf)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_cr

!
!----
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_all_cr


!
!----
!

  subroutine min_all_all_dp(sendbuf, recvbuf)

  implicit none

  double precision :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine min_all_all_dp

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

  subroutine max_allreduce_i(buffer,countval)

  implicit none

  integer :: countval
  integer,dimension(countval),intent(inout) :: buffer

  integer(kind=4) :: unused_i4

  unused_i4 = buffer(1)

  end subroutine max_allreduce_i

!
!----
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

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
  integer :: sendbuf, recvbuf

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

  double precision :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_dp

!
!----
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_cr

!
!----
!

  subroutine sum_all_1Darray_dp(sendbuf, recvbuf, nx)

  implicit none

  integer :: nx
  double precision, dimension(nx) :: sendbuf, recvbuf

  recvbuf(:) = sendbuf(:)

  end subroutine sum_all_1Darray_dp

!
!----
!

  subroutine any_all_1Darray_l(sendbuf, recvbuf, nx)

  implicit none

  integer :: nx
  logical, dimension(nx) :: sendbuf, recvbuf

  recvbuf(:) = sendbuf(:)

  end subroutine any_all_1Darray_l

!
!----
!

  subroutine sum_all_all_cr(sendbuf, recvbuf)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_all_cr

!
!----
!

  subroutine sum_all_i(sendbuf, recvbuf)

  implicit none

  integer :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_i


!
!----
!
  subroutine sum_all_all_i(sendbuf, recvbuf)

  implicit none

  integer :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine sum_all_all_i

!
!----
!

  subroutine any_all_l(sendbuf, recvbuf)

  implicit none

  logical :: sendbuf, recvbuf

  recvbuf = sendbuf

  end subroutine any_all_l

!
!----
!

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer(kind=4) :: unused_i4
  real(kind=CUSTOM_REAL) :: unused_cr

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

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer(kind=4) :: unused_i4
  real(kind=CUSTOM_REAL) :: unused_cr

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

  implicit none

  integer :: sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  integer(kind=4) :: unused_i4

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

  implicit none

  integer :: recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf

  integer(kind=4) :: unused_i4

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

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  integer,dimension(recvcount):: recvbuf

  integer(kind=4) :: unused_i4

  stop 'recv_i not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_i4 = recvbuf(1)

  end subroutine recv_i

!
!----
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: recvcount,dest,recvtag
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf

  integer(kind=4) :: unused_i4
  real(kind=CUSTOM_REAL) :: unused_cr

  stop 'recvv_cr not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_cr = recvbuf(1)

  end subroutine recvv_cr

!
!----
!

  subroutine recv_r(recvbuf, recvcount, dest, recvtag )

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  real,dimension(recvcount) :: recvbuf

  integer(kind=4) :: unused_i4
  real :: unused_r

  stop 'recv_r not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_r = recvbuf(1)

  end subroutine recv_r


!
!----
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  integer, dimension(sendcount) :: sendbuf

  integer(kind=4) :: unused_i4

  stop 'send_i not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_i4 = sendbuf(1)

  end subroutine send_i


!
!----
!

  subroutine send_i_t(sendbuf,sendcount,dest)

  implicit none

  integer :: dest,sendcount
  integer, dimension(sendcount) :: sendbuf

  integer(kind=4) :: unused_i4

  stop 'send_i_t not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendbuf(1)

  end subroutine send_i_t

!
!----
!


  subroutine recv_i_t(recvbuf,recvcount,source)

  implicit none

  integer :: source,recvcount
  integer, dimension(recvcount) :: recvbuf

  integer(kind=4) :: unused_i4

  stop 'recv_i_t not implemented for serial code'
  unused_i4 = source
  unused_i4 = recvbuf(1)

  end subroutine recv_i_t

!
!----
!

  subroutine send_r(sendbuf, sendcount, dest, sendtag)

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  real,dimension(sendcount):: sendbuf

  integer(kind=4) :: unused_i4
  real :: unused_r

  stop 'send_i_t not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_r = sendbuf(1)

  end subroutine send_r

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

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  double precision,dimension(sendcount):: sendbuf

  integer(kind=4) :: unused_i4
  double precision :: unused_dp

  stop 'send_dp not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_dp = sendbuf(1)

  end subroutine send_dp
!
!----
!
  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision,dimension(recvcount) :: recvbuf

  integer(kind=4) :: unused_i4
  double precision :: unused_dp

  stop 'recv_dp not implemented for serial code'
  unused_i4 = dest
  unused_i4 = recvtag
  unused_dp = recvbuf(1)

  end subroutine recv_dp

!
!----
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: sendcount,dest,sendtag
  real(kind=CUSTOM_REAL),dimension(sendcount) :: sendbuf

  integer(kind=4) :: unused_i4
  real(kind=CUSTOM_REAL) :: unused_cr

  stop 'sendv_cr not implemented for serial code'
  unused_i4 = dest
  unused_i4 = sendtag
  unused_cr = sendbuf(1)

  end subroutine sendv_cr
!
!----
!

  subroutine wait_req(req)

  implicit none

  integer :: req

  integer(kind=4) :: unused_i4

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

  use constants, only: mygroup

  implicit none

  mygroup = 0

  end subroutine world_split

!
!----
!

  subroutine world_unsplit()

  end subroutine world_unsplit

!
!----
!

  subroutine bcast_all_l_array(buffer, countval)

  implicit none
  integer  :: countval
  logical, dimension(countval) :: buffer
  logical :: unused_l

  unused_l = buffer(1)

  end subroutine bcast_all_l_array

