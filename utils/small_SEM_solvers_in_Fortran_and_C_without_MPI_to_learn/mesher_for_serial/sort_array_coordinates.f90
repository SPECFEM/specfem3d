!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! subroutines to sort MPI buffers to assemble between chunks

  subroutine sort_array_coordinates(npointot,x,y,z,ibool,iglob,locval,ifseg, &
                                    nglob,ninseg)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

  implicit none

  include "constants.h"

  integer, intent(in) :: npointot
  double precision, dimension(npointot), intent(inout) :: x, y, z
  integer, dimension(npointot), intent(inout) :: ibool
  integer, dimension(npointot), intent(out) :: iglob, locval, ninseg
  logical, dimension(npointot), intent(out) :: ifseg
  integer, intent(out) :: nglob
  double precision xtol

  ! local parameters
  integer :: i, j
  integer :: nseg, ioff, iseg, ig

  ! establish initial pointers
  do i=1,npointot
    locval(i) = i
  enddo

! define a tolerance, normalized radius is 1., so let's use a small value
  xtol = SMALLVALTOL

  ifseg(:) = .false.

  nseg = 1
  ifseg(1) = .true.
  ninseg(1) = npointot

  do j=1,NDIM
    ! sort within each segment
    ioff = 1
    do iseg=1,nseg
      if (j == 1) then

        call heap_sort_multi(ninseg(iseg), x(ioff), y(ioff), z(ioff), ibool(ioff), locval(ioff))

      else if (j == 2) then

        call heap_sort_multi(ninseg(iseg), y(ioff), x(ioff), z(ioff), ibool(ioff), locval(ioff))

      else

        call heap_sort_multi(ninseg(iseg), z(ioff), x(ioff), y(ioff), ibool(ioff), locval(ioff))

      endif

      ioff = ioff + ninseg(iseg)
    enddo

    ! check for jumps in current coordinate
    if (j == 1) then
      do i=2,npointot
        if (dabs(x(i) - x(i-1)) > xtol) ifseg(i) = .true.
      enddo
    else if (j == 2) then
      do i=2,npointot
        if (dabs(y(i) - y(i-1)) > xtol) ifseg(i) = .true.
      enddo
    else
      do i=2,npointot
        if (dabs(z(i) - z(i-1)) > xtol) ifseg(i) = .true.
      enddo
    endif

    ! count up number of different segments
    nseg = 0
    do i=1,npointot
      if (ifseg(i)) then
        nseg = nseg + 1
        ninseg(nseg) = 1
      else
        ninseg(nseg) = ninseg(nseg) + 1
      endif
    enddo
  enddo

  ! assign global node numbers (now sorted lexicographically)
  ig = 0
  do i=1,npointot
    if (ifseg(i)) ig = ig + 1
    iglob(locval(i)) = ig
  enddo

  nglob = ig

  end subroutine sort_array_coordinates

! -------------------- library for sorting routine ------------------

  subroutine heap_sort_multi(N, dx, dy, dz, ia, ib)

  implicit none
  integer, intent(in) :: N
  double precision, dimension(N), intent(inout) :: dx
  double precision, dimension(N), intent(inout) :: dy
  double precision, dimension(N), intent(inout) :: dz
  integer, dimension(N), intent(inout) :: ia
  integer, dimension(N), intent(inout) :: ib

  integer :: i

  ! checks if anything to do
  if (N < 2) return

  ! builds heap
  do i = N/2, 1, -1
    call heap_sort_siftdown(i, n)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    call dswap(dx, 1, i)
    call dswap(dy, 1, i)
    call dswap(dz, 1, i)
    call iswap(ia, 1, i)
    call iswap(ib, 1, i)
    call heap_sort_siftdown(1, i - 1)
  enddo

  contains

    subroutine dswap(A, i, j)

    double precision, dimension(:), intent(inout) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j

    double precision :: tmp

    tmp = A(i)
    A(i) = A(j)
    A(j) = tmp

    end subroutine

    subroutine iswap(A, i, j)

    integer, dimension(:), intent(inout) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j

    integer :: tmp

    tmp = A(i)
    A(i) = A(j)
    A(j) = tmp

    end subroutine

    subroutine heap_sort_siftdown(start, bottom)

    integer, intent(in) :: start
    integer, intent(in) :: bottom

    integer :: i, j
    double precision :: xtmp, ytmp, ztmp
    integer :: atmp, btmp

    i = start
    xtmp = dx(i)
    ytmp = dy(i)
    ztmp = dz(i)
    atmp = ia(i)
    btmp = ib(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
      if (j < bottom) then
        if (dx(j) <= dx(j+1)) j = j + 1
      endif

      ! checks if section already smaller than initial value
      if (dx(j) < xtmp) exit

      dx(i) = dx(j)
      dy(i) = dy(j)
      dz(i) = dz(j)
      ia(i) = ia(j)
      ib(i) = ib(j)
      i = j
      j = 2 * i
    enddo

    dx(i) = xtmp
    dy(i) = ytmp
    dz(i) = ztmp
    ia(i) = atmp
    ib(i) = btmp

    end subroutine heap_sort_siftdown

  end subroutine heap_sort_multi
