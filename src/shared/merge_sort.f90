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

!-----------------------------------------------------------------
! merge sort
!
! see https://en.wikipedia.org/wiki/Merge_sort
!
! Fortran modified version for double precision
! originally from: https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
!
! uses insertion sort for small arrays to increase efficiency
!
! example:
!   ..
!   integer, parameter :: N = 8
!   double precision, dimension(N) :: A = (/ 1., 5., 2., 7., 3., 9., 4., 6. /)
!   double precision, dimension ((N+1)/2) :: T
!   call merge_sort(A,N,T)
!
!-------------------------------------------------------------------

  recursive subroutine merge_sort(A,N,T)

  implicit none
  integer, intent(in) :: N
  double precision, dimension(N), intent(inout) :: A     ! result, sorted in increasing order
  double precision, dimension((N+1)/2), intent(out) :: T ! temporary array

  integer :: NA,NB
  double precision :: val

  integer,parameter :: threshold_size = 16 ! fallback to insertion sort

  ! trivial case: single entry
  if (N < 2) return

  ! trivial case: 2 entries
  if (N == 2) then
    if (A(1) > A(2)) then
      val = A(1)
      A(1) = A(2)
      A(2) = val
    endif
    return
  endif

  ! small array case: insertion sort has less operations for worst case
  if (N <= threshold_size) then
    call InsertionSort(A,N)
    return
  endif

  NA = (N+1)/2
  NB = N-NA

  call MergeSort(A,NA,T)
  call MergeSort(A(NA+1),NB,T)

  if (A(NA) > A(NA+1)) then
    T(1:NA) = A(1:NA)
    call Merge(T,NA,A(NA+1),NB,A,N)
  endif
  return

  end subroutine merge_sort

!
!-------------------------------------------------------------------------------------------------
!

  subroutine merge(A,NA,B,NB,C,NC)

  implicit none
  integer, intent(in)    :: NA,NB,NC         ! Normal usage: NA+NB = NC
  double precision, intent(inout) :: A(NA)   ! B overlays C(NA+1:NC)
  double precision, intent(in)    :: B(NB)
  double precision, intent(inout) :: C(NC)

  integer :: i,j,k

  i = 1
  j = 1
  k = 1
  do while(i <= NA .and. j <= NB)
    if (A(i) <= B(j)) then
       C(k) = A(i)
       i = i+1
    else
       C(k) = B(j)
       j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= NA)
    C(k) = A(i)
    i = i+1
    k = k+1
  enddo

  return

  end subroutine merge

!
!-------------------------------------------------------------------------------------------------
!


! insertion sort: good locality and efficiency for small n
!
! see: https://en.wikipedia.org/wiki/Insertion_sort
!
! modified double precision version,
! originally from: https://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort

  subroutine insertion_sort(A,N)

  implicit none
  integer,intent(in) :: N
  double precision,intent(inout) :: A(N)

  double precision :: x
  integer :: i, j

  do i = 2, N
    x = A(i)
    j = i-1
    do while (j >= 1)
      if (A(j) <= x) exit
      A(j+1) = A(j)
      j = j-1
    enddo
    A(j+1) = x
  enddo

  end subroutine insertion_sort


