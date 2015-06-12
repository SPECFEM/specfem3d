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

module numbering

  use data_gllmesh
  use data_numbering
  use data_spec
  use data_mesh
  use data_grid
  use data_diag

  implicit none
  public :: define_global_global_numbering
  public :: define_global_flobal_numbering
  public :: define_global_slobal_numbering
  public :: get_global
  private

contains

!--------------------------------------------------------------------------
subroutine define_global_global_numbering

  use data_time
  use clocks_mod
  
  real(kind=dp)   , dimension(:), allocatable   :: sgtmp, zgtmp
  logical, dimension(:), allocatable            :: ifseg
  integer, dimension(:), allocatable            :: loc
  integer   :: npointot

  ngllcube = (npol + 1)**2 
  npointot = neltot * (npol+1)**2

  if (dump_mesh_info_screen) then
   write(6,*) 
   write(6,*) 'NPOINTOT GLOBAL IS ' , npointot
  end if

  allocate(sgtmp(npointot))
  sgtmp = pack(sgll, .true.)

  allocate(zgtmp(npointot))
  zgtmp = pack(zgll, .true.)

  allocate(iglob(npointot))
  iglob(:) = 0
  allocate(loc(npointot))
  loc(:) = 0
  allocate(ifseg(npointot))

  iclock07 = tick()
  call get_global(neltot, sgtmp, zgtmp, iglob, loc, ifseg, nglobglob, npointot, &
                  ngllcube)
  iclock07 = tick(id=idold07, since=iclock07)

  deallocate(ifseg)
  deallocate(loc)
  deallocate(sgtmp)
  deallocate(zgtmp)

  if (dump_mesh_info_screen) write(6,*) 'NGLOBGLOB IS ' , NGLOBGLOB

end subroutine define_global_global_numbering
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine define_global_flobal_numbering
  
  use data_time
  use clocks_mod
  
  real(kind=dp)   , dimension(:), allocatable :: sgtmp,zgtmp
  integer, dimension(:), allocatable :: loc_fluid
  logical, dimension(:), allocatable ::   ifseg
  integer :: npointot

  npointot = neltot_fluid * (npol+1)**2

  if (dump_mesh_info_screen) then 
     write(6,*) 
     write(6,*) 'NPOINTOT FLOBAL IS ' , npointot
  end if

  allocate(sgtmp(npointot))
  sgtmp = pack(sgll_fluid, .true.)
  deallocate(sgll_fluid)

  allocate(zgtmp(npointot))
  zgtmp = pack(zgll_fluid, .true.)
  deallocate(zgll_fluid)

  allocate(iglob_fluid(npointot))
  iglob_fluid(:) = 0
  allocate(loc_fluid(npointot))
  loc_fluid(:) = 0
  allocate(ifseg(npointot))

  iclock07 = tick()
  call get_global(neltot_fluid, sgtmp, zgtmp, iglob_fluid, loc_fluid, ifseg, &
                  nglobflob, npointot, NGLLcube)
  iclock07 = tick(id=idold07, since=iclock07)

  deallocate(ifseg)
  deallocate(loc_fluid)
  deallocate(zgtmp)
  deallocate(sgtmp)

  if (dump_mesh_info_screen) write(6,*) 'NGLOBFLOB IS ' , NGLOBFLOB

end subroutine define_global_flobal_numbering
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine define_global_slobal_numbering

  use data_time
  use clocks_mod
  
  real(kind=dp)   , dimension(:), allocatable   :: sgtmp, zgtmp
  integer, dimension(:), allocatable            :: loc_solid
  logical, dimension(:), allocatable            :: ifseg

  integer :: npointot

  npointot = neltot_solid * (npol+1)**2

  if (dump_mesh_info_screen) then 
     write(6,*) 
     write(6,*) 'NPOINTOT SLOBAL IS ' , npointot
  end if

  allocate(sgtmp(npointot))
  sgtmp = pack(sgll_solid, .true.)
  deallocate(sgll_solid) ! not needed anymore 
  
  allocate(zgtmp(npointot))
  zgtmp = pack(zgll_solid, .true.)
  deallocate(zgll_solid) ! not needed anymore

  allocate(iglob_solid(npointot))
  iglob_solid(:) = 0
  allocate(loc_solid(npointot))
  loc_solid(:) = 0
  allocate(ifseg(npointot))

  iclock07 = tick()
  call get_global(neltot_solid, sgtmp, zgtmp, iglob_solid, loc_solid, ifseg, nglobslob, &
                  npointot, NGLLcube)
  iclock07 = tick(id=idold07, since=iclock07)

  deallocate(ifseg)
  deallocate(loc_solid)
  deallocate(zgtmp)
  deallocate(sgtmp)

  if (dump_mesh_info_screen) write(6,*) 'NGLOBSLOB IS ' , NGLOBSLOB 

end subroutine define_global_slobal_numbering
!-------------------------------------------------------------------------


!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

subroutine get_global(nspec2, xp, yp, iglob2, loc2, ifseg2, nglob2, npointot2, &
                      NGLLCUBE2, nmax_threads_in)

  ! this routine MUST be in real(kind=dp)    to avoid sensitivity
  ! to roundoff errors in the coordinates of the points

  ! non-structured global numbering software provided by Paul F. Fischer

  ! leave sorting subroutines in same source file to allow for inlining

  use sorting,      only: mergesort_3
  !$ use omp_lib     

  integer, intent(in)                :: nspec2, npointot2, NGLLCUBE2
  real(kind=dp)   , intent(inout)    :: xp(npointot2), yp(npointot2)
  integer, intent(out)               :: iglob2(npointot2), loc2(npointot2), nglob2
  logical, intent(out)               :: ifseg2(npointot2)
  integer, intent(in), optional      :: nmax_threads_in

  integer                            :: nmax_threads

  integer                            :: ioffs(npointot2)

  integer                            :: ispec, i, nthreads
  integer                            :: ieoff, ilocnum, nseg, ioff, iseg, ig

  integer, dimension(:), allocatable :: ind, ninseg

  real(kind=dp), parameter           :: SMALLVALTOL = 1.d-08

  if (present(nmax_threads_in)) then
     nmax_threads = nmax_threads_in
  else
     nmax_threads = 1
     !$ nmax_threads = omp_get_max_threads()
  endif

  ! establish initial pointers
  do ispec=1, nspec2
     ieoff = NGLLCUBE2 * (ispec - 1)
     do ilocnum=1, NGLLCUBE2
        loc2(ilocnum + ieoff) = ilocnum + ieoff
     enddo
  enddo

  ifseg2 = .false.

  allocate(ind(npointot2))
  allocate(ninseg(npointot2))

  ninseg = 0
  ind = 0

  nseg = 1
  ifseg2(1) = .true.
  ninseg(1) = npointot2

  call mergesort_3(xp, b=yp, il=loc2, p=nmax_threads)

  ! check for jumps in current coordinate
  ! compare the coordinates of the points within a small tolerance
  do i=2, npointot2
     if (dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg2(i) = .true.
  enddo

  ! count up number of different segments
  nseg = 0
  do i=1, npointot2
     if(ifseg2(i)) then
        nseg = nseg + 1
        ninseg(nseg) = 1
     else
        ninseg(nseg) = ninseg(nseg) + 1
     endif

  enddo
  
  ! sort within each segment
  ioff = 1

  ioffs(1) = 1
  do iseg=2, nseg
     ioffs(iseg) = ioffs(iseg-1) + ninseg(iseg-1)
  end do
  
  !$ nthreads = min(min(omp_get_max_threads(),8), nmax_threads)
  !$ call omp_set_num_threads(nthreads)
  !$ print *, 'Using ', nthreads, ' threads for sorting in get_global!'
  !$ call flush(6)
  !$omp parallel do shared(xp,yp,loc2,ninseg) private(ind,ioff)
  do iseg=1, nseg
     ioff = ioffs(iseg)
     call rank_y(yp(ioff), ind, ninseg(iseg))
     call swapall(loc2(ioff), xp(ioff), yp(ioff), ind, ninseg(iseg))
  enddo
  !$omp end parallel do        

  ! check for jumps in current coordinate
  ! compare the coordinates of the points within a small tolerance
  do i=2,npointot2
     if (dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg2(i) = .true.
  enddo

  deallocate(ind)
  deallocate(ninseg)

  ! assign global node numbers (now sorted lexicographically)
  ig = 0
  do i=1, npointot2
     if(ifseg2(i)) ig = ig + 1
     iglob2(loc2(i)) = ig
  enddo
  nglob2 = ig

end subroutine get_global
!-------------------------------------------------------------------------

! sorting routines put in same file to allow for inlining
!-------------------------------------------------------------------------
subroutine rank_y(A,IND,N)
  !
  ! Use Heap Sort (Numerical Recipes)
  !

  integer n
  real(kind=dp)    A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  real(kind=dp)    q

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
200 CONTINUE
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

end subroutine rank_y
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine swapall(IA,A,B,ind,n)
  !
  ! swap arrays IA, A, B and C according to addressing in array IND
  !

  integer, intent(in)             :: n
  integer, intent(in)             :: IND(n)
  integer, intent(inout)          :: IA(n)
  real(kind=dp)   ,intent(inout)  :: A(n),B(n)
  real(kind=dp)                   :: W(n), W2(n)
  integer                         :: IW(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)
  W2(:) = B(:)

  do i=1,n
     IA(i) = IW(ind(i))
     A(i)  = W(ind(i))
     B(i)  = W2(ind(i))
  enddo

end subroutine swapall
!-------------------------------------------------------------------------

!=========================
end module numbering
!=========================
