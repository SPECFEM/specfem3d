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


  subroutine heap_sort( N, array )

! heap sort algorithm
! sorts integer array (in increasing order, like 1 - 5 - 6 - 9 - 12 - 13 - 14 -...)

  implicit none
  integer,intent(in) :: N
  integer,dimension(N),intent(inout) :: array

  ! local parameters
  integer :: tmp
  integer :: i

  ! checks if anything to do
  if (N < 2 ) return

  ! builds heap
  do i = N/2, 1, -1
    call heap_sort_siftdown(N,array,i,N)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    tmp = array(1)
    array(1) = array(i)
    array(i) = tmp
    call heap_sort_siftdown(N,array,1,i-1)
  enddo

  contains

!
!----
!

  subroutine heap_sort_siftdown(N,array,start,bottom)

  implicit none

  integer,intent(in):: N
  integer,dimension(N),intent(inout) :: array
  integer :: start,bottom

  ! local parameters
  integer :: i,j
  integer :: tmp

  i = start
  tmp = array(i)
  j = 2*i
  do while( j <= bottom )
    ! chooses larger value first in this section
    if (j < bottom) then
      if (array(j) <= array(j+1) ) j = j + 1
    endif

    ! checks if section already smaller than initial value
    if (array(j) < tmp ) exit

    array(i) = array(j)
    i = j
    j = 2*i
  enddo

  array(i) = tmp
  return

  end subroutine heap_sort_siftdown

  end subroutine heap_sort

!
!-------------------------------------------------------------------------------------------------
!

!    - Implementation of a Heap Sort Routine
!    Input
!      n = Input
!         Length of arrays
!      X = Input/Output
!         Vector to be sorted
!         dimension(n)
!      Y = Output
!         Sorted Indices of vector X
!
!      Example:
!         X = [ 4 3 1 2 ] on Input
!         Y = [ 1 2 3 4 ] Computed Internally (in order)
!
!         X = [ 1 2 3 4 ] on Output
!         Y = [ 3 4 2 1 ] on Output
!
  subroutine heap_sort_local(N, X, Y)

  implicit none
  integer, intent(in) :: N
  double precision, dimension(N), intent(inout) :: X
  integer, dimension(N), intent(out) :: Y

  ! local parameters
  double precision :: tmp
  integer :: itmp
  integer :: i

  do i = 1,N
     Y(i) = i
  enddo

  ! checks if anything to do
  if (N < 2) return

  ! builds heap
  do i = N/2, 1, -1
    call heap_sort_siftdown(i, N)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    tmp = X(1)
    X(1) = X(i)
    X(i) = tmp
    itmp = Y(1)
    Y(1) = Y(i)
    Y(i) = itmp

    call heap_sort_siftdown(1, i - 1)
  enddo

!
!----
!

  contains

    subroutine heap_sort_siftdown(start, bottom)

    implicit none

    integer, intent(in) :: start, bottom

    ! local parameters
    integer :: i, j
    double precision :: xtmp
    integer :: ytmp

    i = start
    xtmp = X(i)
    ytmp = Y(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
      if (j < bottom) then
        if (X(j) <= X(j+1)) j = j + 1
      endif

      ! checks if section already smaller than initial value
      if (X(j) < xtmp) exit

      X(i) = X(j)
      Y(i) = Y(j)
      i = j
      j = 2 * i
    enddo

    X(i) = xtmp
    Y(i) = ytmp

    end subroutine heap_sort_siftdown

  end subroutine heap_sort_local
