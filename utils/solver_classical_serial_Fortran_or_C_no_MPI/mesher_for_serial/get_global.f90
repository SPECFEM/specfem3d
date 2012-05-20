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

  subroutine get_global(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  include "constants.h"

! parameters
  integer, intent(in) :: npointot,nspec
  double precision, intent(in) :: xp(npointot),yp(npointot),zp(npointot)

  integer, intent(out) :: iglob(npointot),loc(npointot)
  logical, intent(out) :: ifseg(npointot)
  integer, intent(out) :: nglob

! variables
  integer ispec,i,j
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

! dynamically allocate arrays
  allocate(ind(npointot))
  allocate(ninseg(npointot))
  allocate(iwork(npointot))
  allocate(work(npointot))

! establish initial pointers
  do ispec=1,nspec
    ieoff=NGLLX * NGLLY * NGLLZ * (ispec-1)
    do ilocnum=1,NGLLX * NGLLY * NGLLZ
      loc(ilocnum+ieoff)=ilocnum+ieoff
    enddo
  enddo

  ifseg(:)=.false.

  nseg=1
  ifseg(1)=.true.
  ninseg(1)=npointot

do j=1,NDIM

    ! sort within each segment
    ioff=1
    do iseg=1,nseg
        if(j == 1) then
            call rank(xp(ioff),ind,ninseg(iseg))
        else if(j == 2) then
            call rank(yp(ioff),ind,ninseg(iseg))
        else
            call rank(zp(ioff),ind,ninseg(iseg))
        endif
        call swap_all(loc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))
        ioff=ioff+ninseg(iseg)
    enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
    if(j == 1) then
        do i=2,npointot
            if(dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
        enddo
    else if(j == 2) then
        do i=2,npointot
           if(dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
        enddo
    else
        do i=2,npointot
            if(dabs(zp(i)-zp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
        enddo
    endif

! count up number of different segments
    nseg=0
    do i=1,npointot
        if(ifseg(i)) then
        nseg=nseg+1
        ninseg(nseg)=1
        else
        ninseg(nseg)=ninseg(nseg)+1
        endif
    enddo
enddo

! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npointot
    if(ifseg(i)) ig=ig+1
    iglob(loc(i))=ig
  enddo

  nglob=ig

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine get_global

! -----------------------------------

! sorting routines put in same file to allow for inlining

  subroutine rank(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j=1,n
   IND(j)=j
  enddo

  if (n == 1) return

  L=n/2+1
  ir=n
  100 CONTINUE
   IF (l>1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   ENDIF
   i=l
   j=l+l
  200    CONTINUE
   IF (J <= IR) THEN
      IF (J<IR) THEN
         IF ( A(IND(j))<A(IND(j+1)) ) j=j+1
      ENDIF
      IF (q<A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      ENDIF
   goto 200
   ENDIF
   IND(I)=INDX
  goto 100
  end subroutine rank

! ------------------------------------------------------------------

  subroutine swap_all(IA,A,B,C,IW,W,ind,n)
!
! swap arrays IA, A, B and C according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IW(n)
  double precision A(n),B(n),C(n),W(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
    IA(i)=IW(ind(i))
    A(i)=W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
    B(i)=W(ind(i))
  enddo

  W(:) = C(:)

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all

