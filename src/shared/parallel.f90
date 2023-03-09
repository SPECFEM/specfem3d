!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

! Dimitri Komatitsch, July 2014, CNRS Marseille, France:
! added the ability to run several calculations (several earthquakes)
! in an embarrassingly-parallel fashion from within the same run;
! this can be useful when using a very large supercomputer to compute
! many earthquakes in a catalog, in which case it can be better from
! a batch job submission point of view to start fewer and much larger jobs,
! each of them computing several earthquakes in parallel.
! To turn that option on, set parameter NUMBER_OF_SIMULTANEOUS_RUNS to a value greater than 1 in the Par_file.
! To implement that, we create NUMBER_OF_SIMULTANEOUS_RUNS MPI sub-communicators,
! each of them being labeled "my_local_mpi_comm_world", and we use them
! in all the routines in "src/shared/parallel.f90", except in MPI_ABORT() because in that case
! we need to kill the entire run.
! When that option is on, of course the number of processor cores used to start
! the code in the batch system must be a multiple of NUMBER_OF_SIMULTANEOUS_RUNS,
! all the individual runs must use the same number of processor cores,
! which as usual is NPROC in the input file DATA/Par_file,
! and thus the total number of processor cores to request from the batch system
! should be NUMBER_OF_SIMULTANEOUS_RUNS * NPROC.
! All the runs to perform must be placed in directories called run0001, run0002, run0003 and so on
! (with exactly four digits).

!-------------------------------------------------------------------------------------------------
!
! Parallel routines.  All MPI calls belong in this file!
!
!-------------------------------------------------------------------------------------------------

module my_mpi

! main parameter module for specfem simulations

  use mpi

  implicit none

  integer :: my_local_mpi_comm_world, my_local_mpi_comm_for_bcast

end module my_mpi

!-------------------------------------------------------------------------------------------------
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi()

  use my_mpi
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: myrank,ier

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

  ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read before calling world_split()
  ! thus read the parameter file
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
  if (myrank == 0) then
    call open_parameter_file_from_main_only(ier)
    ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read
    call read_value_integer(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NUMBER_OF_SIMULTANEOUS_RUNS'
    call read_value_logical(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter BROADCAST_SAME_MESH_AND_MODEL'
    ! close parameter file
    call close_parameter_file()
  endif

  ! broadcast parameters read from main to all processes
  my_local_mpi_comm_world = MPI_COMM_WORLD
  call bcast_all_singlei(NUMBER_OF_SIMULTANEOUS_RUNS)
  call bcast_all_singlel(BROADCAST_SAME_MESH_AND_MODEL)

  ! create sub-communicators if needed, if running more than one earthquake from the same job
  call world_split()

  end subroutine init_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_mpi()

  use my_mpi

  implicit none

  integer :: ier

! close sub-communicators if needed, if running more than one earthquake from the same job
  call world_unsplit()

  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0) stop 'Error finalizing MPI'

  end subroutine finalize_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine abort_mpi()

  use my_mpi
  use constants, only: MAX_STRING_LEN,mygroup
  use shared_input_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer :: my_local_rank,my_global_rank,ier
  logical :: run_file_exists
  character(len=MAX_STRING_LEN) :: filename

  ! get my local rank and my global rank (in the case of simultaneous jobs, for which we split
  ! the MPI communicator, they will be different; otherwise they are the same)
  call world_rank(my_local_rank)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_global_rank,ier)

  ! write a stamp file to disk to let the user know that the run failed
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    ! notifies which run directory failed
    write(filename,"('run_with_directory_',i4.4,'_failed')") mygroup + 1
    inquire(file=trim(filename), exist=run_file_exists)
    if (run_file_exists) then
      open(unit=9765,file=trim(filename),status='old',position='append',action='write',iostat=ier)
    else
      open(unit=9765,file=trim(filename),status='new',action='write',iostat=ier)
    endif
    if (ier == 0) then
      write(9765,*) 'run ',mygroup+1,' with local rank ',my_local_rank,' and global rank ',my_global_rank,' failed'
      close(9765)
    endif

    ! notifies which rank failed
    write(filename,"('run_with_local_rank_',i8.8,'_and_global_rank_',i8.8,'_failed')") my_local_rank,my_global_rank
    open(unit=9765,file=trim(filename),status='unknown',action='write')
    write(9765,*) 'run with local rank ',my_local_rank,' and global rank ',my_global_rank,' failed'
    close(9765)
  else
    ! note: we already output an OUTPUT_FILES/error_message***.txt file for each failed rank (single runs)
    ! debug
    !write(filename,"('run_with_local_rank_',i8.8,'_failed')") my_local_rank
    !open(unit=9765,file=filename,status='unknown',action='write')
    !write(9765,*) 'run with local rank ',my_local_rank,' failed'
    !close(9765)
  endif

  ! note: MPI_ABORT does not return, it makes the program exit with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  stop 'error, program ended in abort_mpi'

  end subroutine abort_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all()

  use my_mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(my_local_mpi_comm_world,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

  end subroutine synchronize_all


!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all_comm(comm)

  use my_mpi

  implicit none

  integer,intent(in) :: comm

  ! local parameters
  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(comm,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes for specified communicator'

  end subroutine synchronize_all_comm

!
!-------------------------------------------------------------------------------------------------
!

  double precision function wtime()

  use my_mpi

  implicit none

  wtime = MPI_WTIME()

  end function wtime

!
!-------------------------------------------------------------------------------------------------
!

!  integer function null_process()
!  end function null_process

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine test_request(request,flag_result_test)
!  end subroutine test_request

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wait_req(req)

  use my_mpi

  implicit none

  integer :: req

  integer :: ier

  call mpi_wait(req,MPI_STATUS_IGNORE,ier)

  end subroutine wait_req

!
!-------------------------------------------------------------------------------------------------
!

  logical function is_valid_comm(comm)

  use my_mpi

  implicit none

  integer, intent(in) :: comm

  ! tests if communicator is valid
  if (comm == MPI_COMM_NULL) then
    is_valid_comm = .false.
  else
    is_valid_comm = .true.
  endif

  end function is_valid_comm

!-------------------------------------------------------------------------------------------------
!
! MPI broadcasting helper
!
!-------------------------------------------------------------------------------------------------

  subroutine bcast_all_i(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  integer, dimension(countval) :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

  use my_mpi

  implicit none

  integer :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlel(buffer)

  use my_mpi

  implicit none

  logical :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr(buffer, countval)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: countval
  real(kind=CUSTOM_REAL), dimension(countval) :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlecr(buffer)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singlecr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_r(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  real, dimension(countval) :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return

  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  double precision, dimension(countval) :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singledp(buffer)

  use my_mpi

  implicit none

  double precision :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singledp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch_array(buffer,countval,STRING_LEN)

  use my_mpi

  implicit none

  integer :: countval, STRING_LEN

  character(len=STRING_LEN), dimension(countval) :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return

  call MPI_BCAST(buffer,STRING_LEN*countval,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_ch_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l_array(buffer, countval)

  use my_mpi

  implicit none

  integer :: countval
  logical, dimension(countval) :: buffer
  integer :: ier

  ! checks if anything to do
  if (countval == 0) return

  call MPI_BCAST(buffer,countval,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_l_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_string(buffer)

  use my_mpi
  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,MAX_STRING_LEN,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_string

!
!---- broadcast to send the mesh and model to other simultaneous runs
!

  subroutine bcast_all_i_for_database(buffer, countval)

  use my_mpi
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  integer :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return
  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_i_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l_for_database(buffer, countval)

  use my_mpi
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  logical :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return
  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_LOGICAL,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_l_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr_for_database(buffer, countval)

  use my_mpi
  use constants, only: CUSTOM_REAL
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  include "precision.h"

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real(kind=CUSTOM_REAL) :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return
  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_cr_for_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp_for_database(buffer, countval)

  use my_mpi
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  double precision :: buffer

  integer :: ier

  ! checks if anything to do
  if (countval == 0) return
  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_dp_for_database

!
!-------------------------------------------------------------------------------------------------
!
! unused so far...
!
!  subroutine bcast_all_r_for_database(buffer, countval)
!
!  use my_mpi
!  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL
!
!  implicit none
!
!  integer :: countval
!  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
!  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
!  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
!  real :: buffer
!
!  integer :: ier
!
!  ! checks if anything to do
!  if (countval == 0) return
!  if (.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return
!
!  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_for_bcast,ier)
!
!  end subroutine bcast_all_r_for_database
!
!-------------------------------------------------------------------------------------------------
!
! MPI math helper
!
!-------------------------------------------------------------------------------------------------

  subroutine min_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_allreduce_i(buffer,countval)

  use my_mpi

  implicit none

  integer :: countval
  integer,dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  integer,dimension(countval) :: send

  ! checks if anything to do
  if (countval == 0) return

  ! seems not to be supported on all kind of MPI implementations...
  !! DK DK: yes, I confirm, using MPI_IN_PLACE is tricky
  !! DK DK (see the answer at http://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
  !! DK DK      for how to use it right)
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)
  if (ier /= 0) stop 'Allreduce to get max values failed.'

  end subroutine max_allreduce_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_allreduce_singlei(val,recvval)

  use my_mpi

  implicit none

  integer,intent(in) :: val
  integer,intent(out) :: recvval

  ! local parameters
  integer :: ier

  call MPI_ALLREDUCE(val, recvval, 1, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)
  if (ier /= 0) stop 'Allreduce to get single max value failed.'

  end subroutine max_allreduce_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MIN,my_local_mpi_comm_world,ier)

  end subroutine min_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MIN,my_local_mpi_comm_world,ier)

  end subroutine min_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine max_allreduce_cr(sendbuf, recvbuf)
!  end subroutine max_allreduce_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine maxloc_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision, dimension(2) :: sendbuf,recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,my_local_mpi_comm_world,ier)

  end subroutine maxloc_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine any_all_l(sendbuf, recvbuf)

  use my_mpi

  implicit none

  logical :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL,MPI_LOR,my_local_mpi_comm_world,ier)

  end subroutine any_all_l

!
!-------------------------------------------------------------------------------------------------
!

! MPI summations

  subroutine sum_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_dp

!
!-------------------------------------------------------------------------------------------------
!


  subroutine sum_all_1Darray_dp(sendbuf, recvbuf, nx)

  use my_mpi

  implicit none

  integer :: nx
  double precision, dimension(nx) :: sendbuf, recvbuf
  integer :: ier

  ! checks if anything to do
  if (nx == 0) return

  call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_1Darray_dp

!
!-------------------------------------------------------------------------------------------------
!


  subroutine any_all_1Darray_l(sendbuf, recvbuf, nx)

  use my_mpi

  implicit none

  integer :: nx
  logical, dimension(nx) :: sendbuf, recvbuf
  integer :: ier

  ! checks if anything to do
  if (nx == 0) return

  call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_LOGICAL,MPI_LOR,0,my_local_mpi_comm_world,ier)

  end subroutine any_all_1Darray_l

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sum_all_3Darray_dp(sendbuf, recvbuf, nx,ny,nz)
!  end subroutine sum_all_3Darray_dp


!-------------------------------------------------------------------------------------------------
!
! Send/Receive MPI
!
!-------------------------------------------------------------------------------------------------

! asynchronuous send/receive

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_ISEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine isend_i(sendbuf, sendcount, dest, sendtag, req)

  use my_mpi

  implicit none

  integer :: sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_ISEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_i

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine isend_dp(sendbuf, sendcount, dest, sendtag, req)
!  end subroutine isend_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_IRECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine irecv_dp(recvbuf, recvcount, dest, recvtag, req)
!  end subroutine irecv_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_i(recvbuf, recvcount, dest, recvtag, req)

  use my_mpi

  implicit none

  integer :: recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf
  integer :: ier

  call MPI_IRECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_i

!
!-------------------------------------------------------------------------------------------------
!

! synchronuous/blocking send/receive

  subroutine recv_i(recvbuf, recvcount, dest, recvtag )

  use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  integer,dimension(recvcount):: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_r(recvbuf, recvcount, dest, recvtag )

  use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  real,dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_REAL,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_r


!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  use my_mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision,dimension(recvcount):: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: recvcount,dest,recvtag
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recvv_cr


!
!-------------------------------------------------------------------------------------------------
!

!  subroutine recv_singlei(recvbuf, dest, recvtag)
!  end subroutine recv_singlei

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine recv_singlel(recvbuf, dest, recvtag)
!  end subroutine recv_singlel

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine send_ch(sendbuf, sendcount, dest, sendtag)
!  end subroutine send_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  integer,dimension(sendcount):: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_i_t(sendbuf,sendcount,dest)

  use my_mpi

  implicit none

  integer :: dest,sendcount,ier
  integer :: tag = 100
  integer, dimension(sendcount) :: sendbuf

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,tag,my_local_mpi_comm_world,ier)

  end subroutine send_i_t

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_i_t(recvbuf,recvcount,source)

  use my_mpi

  implicit none

  integer :: source,recvcount,ier
  integer :: tag = 100
  integer, dimension(recvcount) :: recvbuf

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,source,tag,my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i_t

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_r(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  real,dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_REAL,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_r


!
!-------------------------------------------------------------------------------------------------
!
!

  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  double precision,dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcount,dest,sendtag
  real(kind=CUSTOM_REAL),dimension(sendcount) :: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine sendv_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sendrecv_cr(sendbuf, sendcount, dest, sendtag, &
!                         recvbuf, recvcount, source, recvtag)
!  end subroutine sendrecv_cr

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine sendrecv_dp(sendbuf, sendcount, dest, sendtag, &
!                         recvbuf, recvcount, source, recvtag)
!  end subroutine sendrecv_dp


!-------------------------------------------------------------------------------------------------
!
! MPI gather helper
!
!-------------------------------------------------------------------------------------------------

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_ALLGATHER(sendbuf,sendcnt,MPI_INTEGER, &
                     recvbuf,recvcount,MPI_INTEGER, &
                     my_local_mpi_comm_world,ier)

  end subroutine gather_all_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  use my_mpi

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,CUSTOM_MPI_TYPE, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_dp

!
!-------------------------------------------------------------------------------------------------
!

! unused so far...
!
!  subroutine gather_all_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)
!
!  use my_mpi
!
!  implicit none
!
!  integer :: sendcnt, recvcount, NPROC
!  double precision, dimension(sendcnt) :: sendbuf
!  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf
!
!  integer :: ier
!
!  call MPI_ALLGATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
!                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
!                  my_local_mpi_comm_world,ier)
!
!  end subroutine gather_all_all_dp
!

!
!-------------------------------------------------------------------------------------------------
!


  subroutine gather_all_all_cr(sendbuf, recvbuf, counts, NPROC)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: NPROC,counts
  real(kind=CUSTOM_REAL), dimension(counts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(counts,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_ALLGATHER(sendbuf,counts,CUSTOM_MPI_TYPE,recvbuf,counts,CUSTOM_MPI_TYPE,my_local_mpi_comm_world,ier)

  end subroutine gather_all_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_i(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use my_mpi

  implicit none

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,MPI_INTEGER,recvbuf,recvcount,recvoffset,MPI_INTEGER, &
                   0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_i


!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                   recvbuf,recvcount,recvoffset,CUSTOM_MPI_TYPE, &
                   0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_cr


!
!-------------------------------------------------------------------------------------------------
!

  subroutine all_gather_all_i(sendbuf, recvbuf, NPROC)

  use constants
  use my_mpi

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(NPROC) :: recvbuf

  integer :: ier

  call MPI_Allgather(sendbuf,1,MPI_INTEGER, &
                     recvbuf,1,MPI_INTEGER, &
                     my_local_mpi_comm_world,ier)

  end subroutine all_gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine all_gather_all_r(sendbuf, sendcnt, recvbuf, recvcnt, recvoffset, dim1, NPROC)

  use constants
  use my_mpi

  implicit none

  integer :: sendcnt, dim1, NPROC

  real, dimension(sendcnt) :: sendbuf
  real, dimension(dim1, NPROC) :: recvbuf

  integer, dimension(NPROC) :: recvoffset, recvcnt

  integer :: ier

  call MPI_Allgatherv(sendbuf,sendcnt,MPI_REAL, &
                      recvbuf,recvcnt,recvoffset,MPI_REAL, &
                      my_local_mpi_comm_world,ier)

  end subroutine all_gather_all_r


!
!-------------------------------------------------------------------------------------------------
!

  subroutine all_gather_all_ch(sendbuf, sendcnt, recvbuf, recvcnt, recvoffset, dim1, dim2, NPROC)

  use constants
  use my_mpi

  implicit none

  integer :: sendcnt, dim1, dim2, NPROC

  character(len=dim2), dimension(sendcnt) :: sendbuf
  character(len=dim2), dimension(dim1, NPROC) :: recvbuf

  integer, dimension(NPROC) :: recvoffset, recvcnt

  integer :: ier

  call MPI_Allgatherv(sendbuf,sendcnt,MPI_CHARACTER, &
                      recvbuf,recvcnt,recvoffset,MPI_CHARACTER, &
                      my_local_mpi_comm_world,ier)

  end subroutine all_gather_all_ch


! not used yet...
!  subroutine gatherv_all_dp(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)
!
!  use my_mpi
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  include "precision.h"
!
!  integer :: sendcnt,recvcounttot,NPROC
!  integer, dimension(NPROC) :: recvcount,recvoffset
!  double precision, dimension(sendcnt) :: sendbuf
!  double precision, dimension(recvcounttot) :: recvbuf
!
!  integer :: ier
!
!  call MPI_GATHERV(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
!                  recvbuf,recvcount,recvoffset,MPI_DOUBLE_PRECISION, &
!                  0,my_local_mpi_comm_world,ier)
!
!  end subroutine gatherv_all_dp

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine gatherv_all_r(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)
!  end subroutine gatherv_all_r

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine scatter_all_singlei(sendbuf, recvbuf, NPROC)
!  end subroutine scatter_all_singlei


!-------------------------------------------------------------------------------------------------
!
! MPI world helper
!
!-------------------------------------------------------------------------------------------------

  subroutine world_size(sizeval)

  use my_mpi

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(my_local_mpi_comm_world,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_size_comm(sizeval,comm)

  use my_mpi

  implicit none

  integer,intent(out) :: sizeval
  integer,intent(in) :: comm

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(comm,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

  end subroutine world_size_comm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

  use my_mpi

  implicit none

  integer,intent(out) :: rank

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(my_local_mpi_comm_world,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

  end subroutine world_rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank_comm(rank,comm)

  use my_mpi

  implicit none

  integer,intent(out) :: rank
  integer,intent(in) :: comm

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(comm,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

  end subroutine world_rank_comm


!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_duplicate(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm
  integer :: ier

  call MPI_COMM_DUP(my_local_mpi_comm_world,comm,ier)
  if (ier /= 0) stop 'error duplicating my_local_mpi_comm_world communicator'

  end subroutine world_duplicate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm

  comm = my_local_mpi_comm_world

  end subroutine world_get_comm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm_self(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm

  comm = MPI_COMM_SELF

  end subroutine world_get_comm_self


!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_comm_free(comm)

  use my_mpi

  implicit none

  integer,intent(inout) :: comm

  ! local parameters
  integer :: ier

  call MPI_Comm_free(comm,ier)
  if (ier /= 0 ) stop 'Error freeing MPI communicator'

  end subroutine world_comm_free

!
!-------------------------------------------------------------------------------------------------
!

!  subroutine world_get_info_null(info)
!  end subroutine world_get_info_null

!
!-------------------------------------------------------------------------------------------------
!

! create sub-communicators if needed, if running more than one earthquake from the same job.
! create a sub-communicator for each independent run;
! if there is a single run to do, then just copy the default communicator to the new one
  subroutine world_split()

  use my_mpi
  use constants, only: MAX_STRING_LEN,OUTPUT_FILES, &
    IMAIN,ISTANDARD_OUTPUT,mygroup,I_should_read_the_database
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: sizeval,myrank,ier,key,my_group_for_bcast,my_local_rank_for_bcast,NPROC

  character(len=MAX_STRING_LEN) :: path_to_add

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 0) stop 'NUMBER_OF_SIMULTANEOUS_RUNS <= 0 makes no sense'

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mod(sizeval,NUMBER_OF_SIMULTANEOUS_RUNS) /= 0) then
    if (myrank == 0) then
      print *,'Error: the number of MPI processes ',sizeval, &
              ' is not a multiple of NUMBER_OF_SIMULTANEOUS_RUNS = ',NUMBER_OF_SIMULTANEOUS_RUNS
      print *
      print *,'make sure to launch program with NPROC * NUMBER_OF_SIMULTANEOUS_RUNS processes'
      print *,'for example: NPROC = 1 and NUMBER_OF_SIMULTANEOUS_RUNS = 4'
      print *,' > mpirun -np 4 ./bin/xspecfem3D'
      print *
    endif
    stop 'the number of MPI processes is not a multiple of NUMBER_OF_SIMULTANEOUS_RUNS'
  endif

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. IMAIN == ISTANDARD_OUTPUT) &
    stop 'must not have IMAIN == ISTANDARD_OUTPUT when NUMBER_OF_SIMULTANEOUS_RUNS > 1 otherwise output to screen is mingled'

  if (NUMBER_OF_SIMULTANEOUS_RUNS == 1) then

    my_local_mpi_comm_world = MPI_COMM_WORLD

! no broadcast of the mesh and model databases to other runs in that case
    my_group_for_bcast = 0
    my_local_mpi_comm_for_bcast = MPI_COMM_NULL

  else

!--- create a subcommunicator for each independent run

    NPROC = sizeval / NUMBER_OF_SIMULTANEOUS_RUNS

!   create the different groups of processes, one for each independent run
    mygroup = myrank / NPROC
    key = myrank
    if (mygroup < 0 .or. mygroup > NUMBER_OF_SIMULTANEOUS_RUNS-1) stop 'invalid value of mygroup'

!   build the sub-communicators
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, mygroup, key, my_local_mpi_comm_world, ier)
    if (ier /= 0) stop 'error while trying to create the sub-communicators'

!   add the right directory for that run
!   (group numbers start at zero, but directory names start at run0001, thus we add one)
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    OUTPUT_FILES = path_to_add(1:len_trim(path_to_add))//OUTPUT_FILES(1:len_trim(OUTPUT_FILES))

!--- create a subcommunicator to broadcast the identical mesh and model databases if needed
    if (BROADCAST_SAME_MESH_AND_MODEL) then

      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
!     to broadcast the model, split along similar ranks per run instead
      my_group_for_bcast = mod(myrank,NPROC)
      key = myrank
      if (my_group_for_bcast < 0 .or. my_group_for_bcast > NPROC-1) stop 'invalid value of my_group_for_bcast'

!     build the sub-communicators
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, my_group_for_bcast, key, my_local_mpi_comm_for_bcast, ier)
      if (ier /= 0) stop 'error while trying to create the sub-communicators'

!     see if that process will need to read the mesh and model database and then broadcast it to others
      call MPI_COMM_RANK(my_local_mpi_comm_for_bcast,my_local_rank_for_bcast,ier)
      if (my_local_rank_for_bcast > 0) I_should_read_the_database = .false.

    else

! no broadcast of the mesh and model databases to other runs in that case
      my_group_for_bcast = 0
      my_local_mpi_comm_for_bcast = MPI_COMM_NULL

    endif

  endif

  end subroutine world_split

!
!-------------------------------------------------------------------------------------------------
!

! close sub-communicators if needed, if running more than one earthquake from the same job.
  subroutine world_unsplit()

  use my_mpi
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: ier

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    call MPI_COMM_FREE(my_local_mpi_comm_world,ier)
    if (BROADCAST_SAME_MESH_AND_MODEL) call MPI_COMM_FREE(my_local_mpi_comm_for_bcast,ier)
  endif

  end subroutine world_unsplit

