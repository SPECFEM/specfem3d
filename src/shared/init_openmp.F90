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


  subroutine init_openmp()

! outputs OpenMP support info

#ifdef USE_OPENMP
  use constants, only: myrank,IMAIN
#endif

  implicit none

#ifdef USE_OPENMP
  ! local parameters
  integer :: thread_id,num_threads
  integer :: num_procs,max_threads
  logical :: is_dynamic,is_nested
  ! OpenMP functions
  integer,external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer,external :: OMP_GET_NUM_PROCS,OMP_GET_MAX_THREADS
  logical,external :: OMP_GET_DYNAMIC,OMP_GET_NESTED

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(myrank) &
!$OMP PRIVATE(thread_id,num_threads,num_procs,max_threads,is_dynamic,is_nested)
  ! gets thread number
  thread_id = OMP_GET_THREAD_NUM()

  ! gets total number of threads for this MPI process
  num_threads = OMP_GET_NUM_THREADS()

  ! OpenMP main thread only
  if (thread_id == 0) then
    ! gets additional environment info
    num_procs = OMP_GET_NUM_PROCS()
    max_threads = OMP_GET_MAX_THREADS()
    is_dynamic = OMP_GET_DYNAMIC()
    is_nested = OMP_GET_NESTED()

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'OpenMP information:'
      write(IMAIN,*) '  number of threads (per MPI process) = ', num_threads
      write(IMAIN,*)
      write(IMAIN,*) '  number of processors available      = ', num_procs
      write(IMAIN,*) '  maximum number of threads available = ', num_procs
      write(IMAIN,*) '  dynamic thread adjustement          = ', is_dynamic
      write(IMAIN,*) '  nested parallelism                  = ', is_nested
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif
!$OMP END PARALLEL

#else
  ! nothing to do..
  return
#endif

  end subroutine init_openmp

