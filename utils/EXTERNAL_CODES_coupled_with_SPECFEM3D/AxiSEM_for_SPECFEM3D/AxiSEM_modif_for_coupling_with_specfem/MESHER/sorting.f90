!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

module sorting
  use global_parameters, only: sp, dp
  implicit none

  private
  public    :: mergesort_3

contains

!-----------------------------------------------------------------------------------------
subroutine smerge(a, b, c, ind1, ind2, ind)
! merge two sorted arrays a and b into the array c, such that c is also sorted
! (starting with smallest elements)
! carry on index arrays
 
  integer                    :: na, nb
  real(kind=dp), intent(in)  :: a(:)
  real(kind=dp), intent(in)  :: b(:)

  real(kind=dp), intent(out) :: c(size(a)+size(b))
  
  integer, dimension(size(a)), intent(in)           :: ind1
  integer, dimension(size(b)), intent(in)           :: ind2
  integer, dimension(size(b)+size(a)), intent(out)  :: ind
 
  integer                       :: i, j, k
  
  na = size(a)
  nb = size(b)

  i = 1
  j = 1
  k = 1

  do while(i <= na .and. j <= nb)
     if (a(i) <= b(j)) then
        c(k) = a(i)
        ind(k) = ind1(i)
        i = i + 1
     else
        c(k) = b(j)
        ind(k) = ind2(j)
        j = j + 1
     endif
     k = k + 1
  enddo

  do while (i <= na)
     c(k) = a(i)
     ind(k) = ind1(i)
     i = i + 1
     k = k + 1
  enddo

  do while (j <= nb)
     c(k) = b(j)
     ind(k) = ind2(j)
     j = j + 1
     k = k + 1
  enddo

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmerge(a, b, c, ind1, ind2, ind, p)
! split two sorted arrays a and b into p sections such, that the elements of
! each pair of splits form a continuous sequence in a merged and sorted array c,
! then merge in trivially parallel fashion
! see S. Odeh et al (2012), Merge Path - Parallel Merging Made Simple
! carry on index arrays
  
  !$ use omp_lib     
 
  integer                    :: na, nb
  real(kind=dp), intent(in)  :: a(:)
  real(kind=dp), intent(in)  :: b(:)
  real(kind=dp), intent(out) :: c(size(a)+size(b))
  
  integer, dimension(size(a)), intent(in)           :: ind1
  integer, dimension(size(b)), intent(in)           :: ind2
  integer, dimension(size(b)+size(a)), intent(out)  :: ind

  integer, intent(in)           :: p
 
  integer                       :: i, diagnum(0:p)
  integer                       :: ia(0:p), ib(0:p), iamin, iamax, sa, sb
  !$ integer                    :: nthreads
  
  
  na = size(a)
  nb = size(b)

  ia(0) = 0
  ib(0) = 0
  diagnum(0) = 0

  if (p > 1) then
     ! this loop is not parallel, but scales with O(p * log(n))
     do i=1, p
        ! find intersection of merge path with i'th diagonal
        if (i < p) then
           diagnum(i) = i * (na + nb) / p
           
           iamin = 0
           iamax = diagnum(i)

           if (diagnum(i) > na) iamax = na
           if (diagnum(i) > nb) iamin = diagnum(i) - nb

           ! binary search
           do
              sa = (iamin + iamax + 1) / 2
              sb = diagnum(i) - sa + 1
              if (a(sa) < b(sb)) then
                  iamin = sa
              else
                  iamax = sa - 1
              endif
              if (iamin == iamax) exit
           enddo

           ia(i) = iamax 
           ib(i) = diagnum(i) - iamax
        else 
           ! last split is trivial
           ia(i) = na
           ib(i) = nb 
           diagnum(i) = na + nb
        endif

     enddo
     
     !$ nthreads = min(OMP_get_max_threads(), p)
     !$ if (nthreads > 1) then
        !$ call omp_set_num_threads(nthreads)
        !$omp parallel shared (a, b, c, ia, ib, diagnum)
        
        ! this loop is embarrassingly parallel
        !$omp do 
        !$ do i=1, p
        !$    ! merge splitted arrays
        !$    call smerge(a(ia(i-1)+1:ia(i)), b(ib(i-1)+1:ib(i)), &
        !$                c(diagnum(i-1)+1:diagnum(i)), &
        !$                ind1(ia(i-1)+1:ia(i)), ind2(ib(i-1)+1:ib(i)), &
        !$                ind(diagnum(i-1)+1:diagnum(i)))
        !$ enddo
        !$omp end do
        !$omp end parallel
     !$ else
     do i=1, p
        ! merge splitted arrays
        call smerge(a(ia(i-1)+1:ia(i)), b(ib(i-1)+1:ib(i)), c(diagnum(i-1)+1:diagnum(i)), &
                    ind1(ia(i-1)+1:ia(i)), ind2(ib(i-1)+1:ib(i)), &
                    ind(diagnum(i-1)+1:diagnum(i)))
     enddo
     !$ endif
  else
     call smerge(a, b, c, ind1, ind2, ind)
  endif

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
recursive subroutine pmsort(a, ind, p)
! parallel mergesort parallalizing both splitting (in the largest power of 2 smaller 
! or equal p) and merging (exactly p threads) using merge path

  real(kind=dp), intent(inout)      :: a(:)
  integer, intent(inout)            :: ind(size(a))
  integer, intent(in)               :: p
  integer                           :: na, n1, n2
  real(kind=dp)                     :: buf
  integer                           :: ibuf
  real(kind=dp), dimension((size(a)+1)/2)            :: t1
  real(kind=dp), dimension(size(a) - (size(a)+1)/2)  :: t2
  integer, dimension((size(a)+1)/2)                  :: ind1
  integer, dimension(size(a) - (size(a)+1)/2)        :: ind2
 
  na = size(a)

  if (na < 2) then
     return
  elseif (na == 2) then
     if (a(1) > a(2)) then
        buf = a(1)
        a(1) = a(2)
        a(2) = buf
        ibuf = ind(1)
        ind(1) = ind(2)
        ind(2) = ibuf
     endif
     return
  endif

  n1 = (na + 1) / 2
  n2 = na - n1

  if (p < 2) then
     call pmsort(a(1:n1), ind(1:n1), 1)
     call pmsort(a(n1+1:), ind(n1+1:), 1) 
     if (a(n1) > a(n1+1)) then
        ! merge needs memcopies of first inputs!
        t1 = a(1:n1)
        ind1 = ind(1:n1)
        call smerge(t1, a(n1+1:), a, ind1, ind(n1+1:), ind)
     endif
  else
     !$omp parallel sections shared(a, n1, n2) 
     !$omp section
     call pmsort(a(1:n1), ind(1:n1), p/2)
     !$omp section
     call pmsort(a(n1+1:), ind(n1+1:), p/2)
     !$omp end parallel sections
     
     if (a(n1) > a(n1+1)) then
        ! pmerge needs memcopies of both inputs!
        t1 = a(1:n1)
        t2 = a(n1+1:)
        ind1 = ind(1:n1)
        ind2 = ind(n1+1:)
        call pmerge(t1, t2, a, ind1, ind2, ind, p)
     endif
  endif

end subroutine 
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine mergesort_3(a, b, il, il2, p)
    
  use data_time
  use clocks_mod
  !$ use omp_lib     
  
  real(kind=dp), intent(inout) :: a(:)
  real(kind=dp), intent(inout), optional :: b(size(a))
  integer, intent(inout), optional       :: il(size(a))
  integer, intent(inout), optional       :: il2(size(a))
  integer, intent(in), optional          :: p
  
  real(kind=dp), allocatable   :: bbuff(:)
  integer, allocatable         :: ibuff(:)
  integer, allocatable         :: i2buff(:)
  integer                      :: na, i, p_loc
  integer, allocatable         :: ind(:)
  !$ integer                   :: nthreads
 
  if (present(p)) then
     p_loc = p
  else
     p_loc = 1
  endif

  na = size(a)
  allocate(ind(na))
  
  do i=1, na
     ind(i) = i
  enddo
  
  iclock01 = tick()
  call pmsort(a, ind, p_loc) 
  iclock01 = tick(id=idold01, since=iclock01)

  iclock05 = tick()
  if (present(il)) then
     allocate(ibuff(na))
     ibuff(:) = il(:)
  endif

  if (present(il2)) then
     allocate(i2buff(na))
     i2buff(:) = il2(:)
  endif

  if (present(b)) then
     allocate(bbuff(na))
     bbuff(:) = b(:)
  endif
  
  !$ nthreads = min(OMP_get_max_threads(), p_loc)
  !$ call omp_set_num_threads(nthreads)
  !$omp parallel shared (il, ibuff, il2, i2buff, b, bbuff)
  
  if (present(il)) then
     !$omp do
     do i=1, na
        il(i) = ibuff(ind(i))
     enddo
     !$omp end do
  endif
  
  if (present(il2)) then
     !$omp do
     do i=1, na
        il2(i) = i2buff(ind(i))
     enddo
     !$omp end do
  endif
  
  if (present(b)) then
     !$omp do
     do i=1, na
        b(i) = bbuff(ind(i))
     enddo
     !$omp end do
  endif
  !$omp end parallel
  
  iclock05 = tick(id=idold05, since=iclock05)

end subroutine

end module
