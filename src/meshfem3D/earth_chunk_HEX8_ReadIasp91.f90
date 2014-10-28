!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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
!===========================================================================
!

  subroutine earth_chunk_HEX8_ReadIasp91(vp,vs,rho,rb,n)

  implicit none

  integer i,j,n,iunit,nlay,nco(n),ifanis
  double precision vp(n,4),vs(n,4),rho(n,4),rb(n)
  real fref,vph,vsh,qm,qk,eta
  character text*80,cnlay*2,form*11
  do i=1,n
     !qm(i)=0.d0
        !qk(i)=0.d0
     rb(i)=0.d0
     !iflso(i)=0
     nco(i)=0
     do j=1,4
        rho(i,j)=0.d0
        vp(i,j)=0.d0
        !vph(i,j)=0.d0
        vs(i,j)=0.d0
        !vsh(i,j)=0.d0
        !eta(i,j)=0.d0
     enddo
  enddo
  iunit=26
  open(unit=iunit,file='iasp91',status='old')


1 read(iunit,'(a72)') text
  if (text(1:1)=='#') then
     goto 1
  endif
  backspace iunit


  read(iunit,'(i2)') nlay                ! Number of layers

  write(cnlay,'(i2)') nlay               !
  form='('//cnlay//'i2)'                 ! Number of polynomal
  read(iunit,form) (nco(i),i=1,nlay)     ! coefficients for each layer

  read(iunit,*) fref               ! reference frequency of Qs in Hertz
  read(iunit,*) ifanis             ! Transversal isotropic? 1=y, else=n
  read(iunit,'(1x/1x/)')


  do i = 1, nlay

     !read(iunit,*) rb(i-1),rho(i,1),vpv(i,1),vph(i,1),vsv(i,1),vsh(i,1),qm(i),qk(i),eta(i,1)
     read(iunit,*) rb(i),rho(i,1),vp(i,1),vph,vs(i,1),vsh,qm,qk,eta
     !write(*,*) i,rb(i),rho(i,1)
     do j = 2, nco(i)
        read(iunit,*) rho(i,j),vp(i,j),vph,vs(i,j),vsh,eta
        !write(*,*) i,j,rho(i,j)
     enddo
     read(iunit,'(1x)')
  enddo
  i = nlay+1
  read(iunit,*) rb(i)
  j = 1
  rho(i,j) =  rho(i-1,j)
  vp(i,j) = vp(i-1,j)
  vs(i,j) = vs(i-1,j)

  return

  end subroutine earth_chunk_HEX8_ReadIasp91

!
!===========================================================================
!

  subroutine Read_dsm_model(model_file,vp,vs,rho,rb,n)

  implicit none

  integer i,n,iunit,nco(n)
  double precision vp(n,4),vs(n,4),rho(n,4),rb(n),eta(4),vrmin,vrmax
  real vph(4),vsh(4),qm,qk
  integer nzone

  character(len=250) model_file

  rb    = 0.d0
  rho   = 0.d0
  vp    = 0.d0
  vs    = 0.d0
  nco   = 0
  iunit = 26

  open(unit=iunit,file=trim(model_file),status='old',action='read')

  read(iunit,*) nzone

  do i=1, nzone
     read(iunit,*) vrmin, vrmax, &
          rho(i,1), rho(i,2), rho(i,3), rho(i,4), &
          vp(i,1), vp(i,2), vp(i,3), vp(i,4), &
          vph(1), vph(2), vph(3), vph(4), &
          vs(i,1), vs(i,2), vs(i,3), vs(i,4), &
          vsh(1), vsh(2), vsh(3), vsh(4), &
          eta(1), eta(2), eta(3), eta(4),&
          qm, qk
          rb(i)=vrmin
  enddo

  i        = nzone+1
  rb(i)    = vrmax
  vp(i,:)  = vp(i-1,:)
  vs(i,:)  = vs(i-1,:)
  rho(i,:) = rho(i-1,:)

  close(iunit)

  end subroutine Read_dsm_model

!
!===========================================================================
!

subroutine Lyfnd(r,rb,n,i)
  implicit none
  integer i,n
  double precision r,rb(n)
  i=1
  do while (r > rb(i) )
     i = i + 1
  enddo
  i = i - 1
  return
end subroutine Lyfnd

function IsNewLayer(x,r,n)
  implicit none
  integer IsNewLayer,n,i
  double precision x,r(n)
  IsNewLayer = 0
  ! ce test fonctionne que si les mailles sont suffisament petites !! ATTENTION
  do i = 1, n-1
     !write(*,*) x,r(i),x-r(i)
     if (abs(x-r(i))<1.d-10) then
        !write(*,*) 'its layer'
        IsNewLayer = 1
        return
     endif
  enddo
end function IsNewLayer


subroutine StorePoint(z,k,zc)
  implicit none
  integer k
  double precision z(*),zc
  if (k==0) then
     k = k + 1
     z(k) = zc
     return
  else
     if (z(k)==zc) then
        return
     else
        k = k + 1
        z(k) = zc
     endif
  endif
end subroutine StorePoint

subroutine StorePointZ(z,k,zc,NoInter)
  implicit none
  integer k
  double precision z(*),zc
  logical NoInter
  if (k==0) then
     k = k + 1
     z(k) = zc
     return
  else
     if (z(k)==zc.and.NoInter) then
        return
     else
        k = k + 1
        z(k) = zc
     endif
  endif
end subroutine StorePointZ

!!$subroutine FindNum(i,X,X0,n)
!!$  implicit none
!!$  integer i,n
!!$  real X(n),X0
!!$  do i=1,n
!!$     if (X0==X(n)) return
!!$  enddo
!!$  write(*,*) ' warrnig FindnNum ',X0
!!$  stop
!!$end subroutine FindNum
 subroutine CalGridProf(ProfForGemini,Niveau_elm,zlayer,nlayer,NEX_GAMMA,Z_DEPTH_BLOCK)

  implicit none
  integer NEX_GAMMA,nlayer,nbbloc(100000),Niveau_elm(0:NEX_GAMMA-1)
  double precision ProfForGemini(0:NEX_GAMMA-1,3),zlayer(nlayer)
  double precision Z_DEPTH_BLOCK,zpoint(100000),zz(100000)
  double precision epsillon
  integer nb, n, i,j,k,ilayer,ilay,nd,niveau
  double precision p, pas, longeur
  logical test

  epsillon=1d-3
   nbbloc(:)=0
   ! point de depart
   zpoint(1)=zlayer(nlayer) - Z_DEPTH_BLOCK
   write(*,*) zlayer(nlayer) ,  Z_DEPTH_BLOCK
   !! niveau de depart
   call FindLayer_in_mesh_chunk(ilayer,zlayer,zpoint(1),nlayer)
   write(*,*) '              INITIALISATION calcul du niveau de depart : '
   write(*,*)
   write(*,*) 'zlayer : ', zlayer
   write(*,*) 'premier point : '   , zpoint(1),ilayer
    write(*,*)

  !!


  !! on compte le nombre d'elements par niveau
  i = 1
  k = ilayer - 1
  nb = 0
  do while (zpoint(i)<zlayer(nlayer))
    i = i + 1
    k = k + 1
    zpoint(i) = zlayer(k)
  enddo

  nb = i
  nd = i-1
  longeur = zlayer(nlayer) - zpoint(1)


  do i=1,nb-1

     pas = zpoint(i+1) - zpoint(i)
     p = NEX_GAMMA * pas / longeur

    if (p < 0.8d0) then
        n = 1
    else
        n = max(int(p),2)
    endif
    !write(*,*) 'n :',n

    nbbloc(i)=n
  enddo






  do j=1,nb-1
    write(*,*) j,nbbloc(j)
  enddo


  !! on elimine les blocs en trop
   write(*,*) 'SUM ',sum(nbbloc)

   nb = sum(nbbloc)


   do while ( nb > NEX_GAMMA)

      k  =  1
      test = .true.

    do  while (test)

         j =  maxval(nbbloc)
         ! on cherche l'indice du max

         if (j == nbbloc(k) ) then
            nbbloc(k ) = nbbloc(k) -1
            test = .false.
         endif

         k = k + 1

      enddo

      nb = sum(nbbloc)
      write(*,*) 'nb, ',nb,NEX_GAMMA
   enddo

  !!
  longeur = zlayer(nlayer) - zpoint(1)
  k=1
  zz(k)=zpoint(1)
  !zpoint(nb+1)=zlayer(nlayer)
  do i=1,nd
     !write(*,*) i,nbbloc(i)
     pas = (zpoint(i+1) - zpoint(i)) / nbbloc(i)
     write(*,*) i,nbbloc(i),pas
     do while (zz(k) < zpoint(i+1) - epsillon)
        k = k + 1
        zz(k) = zz(k-1) + pas
        write(*,*) zz(k), zpoint(i+1)
     enddo
  enddo




   do ilay=1,NEX_GAMMA

      ProfForGemini(ilay-1,1)  =  zz(ilay)
      ProfForGemini(ilay-1,2)  =  zz(ilay+1)
      ProfForGemini(ilay-1,3)  = 0.5d0 * (zz(ilay) + zz(ilay+1))

      call FindLayer_in_mesh_chunk(niveau,zlayer, ProfForGemini(ilay-1,3),nlayer)
      Niveau_elm(ilay-1)=niveau
      write(*,'(i5,2f15.3,i10)') ilay,zz(ilay),zz(ilay+1),niveau
   enddo

   do ilay=0,NEX_GAMMA-1
     !write(*,*) ProfForGemini(ilay,1), ProfForGemini(ilay,2), ProfForGemini(ilay,3)
   enddo


 end subroutine CalGridProf

 subroutine  FindLayer_for_earth_chunk_mesh(i,z,r,n)
   implicit none
   integer i,n
   double precision z(n),r

   if (r>z(n) .or. r<z(1)) then
    write(*,*) 'STOP :: point ouside grid'
    stop
   endif
   i = 1
   do while (r > z(i))
     i = i + 1
   enddo


 end subroutine FindLayer_for_earth_chunk_mesh





!=====================================================================

  subroutine get_global1(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,UTM_X_MIN,UTM_X_MAX)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  !include "constants.h"

  integer npointot
  integer nspec,nglob
  integer iglob(npointot),loc(npointot)
  logical ifseg(npointot)
  double precision xp(npointot),yp(npointot),zp(npointot)
  double precision UTM_X_MIN,UTM_X_MAX

  integer ispec,i,j,NGLLCUBE1
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
  double precision SMALLVALTOL

  NGLLCUBE1=8
! define geometrical tolerance based upon typical size of the model
  SMALLVALTOL = 1.d-10 * dabs(UTM_X_MAX - UTM_X_MIN)
  write(*,*) dabs(UTM_X_MAX - UTM_X_MIN)
  write(*,*) ' SMALLVALTOL  ',SMALLVALTOL
! dynamically allocate arrays
  allocate(ind(npointot))
  allocate(ninseg(npointot))
  allocate(iwork(npointot))
  allocate(work(npointot))

! establish initial pointers
  do ispec=1,nspec
    ieoff=NGLLCUBE1*(ispec-1)
    do ilocnum=1,NGLLCUBE1
      loc(ilocnum+ieoff)=ilocnum+ieoff
    enddo
  enddo

  ifseg(:)=.false.

  nseg=1
  ifseg(1)=.true.
  ninseg(1)=npointot

  do j=1,3 !,NDIM

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

!! verif VM
!   open(27,file='zp_rank.txt')

!  do i=1,npointot

!     write(27,*) zp(i)
!  enddo
!  close(27)

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine get_global1


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
   endif
   i=l
   j=l+l
  200    CONTINUE
   IF (J <= IR) THEN
      IF (J<IR) THEN
         IF ( A(IND(j))<A(IND(j+1)) ) j=j+1
      endif
      IF (q<A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      endif
   goto 200
   endif
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


!!$!=====================================================================
!!$
!!$! 3D shape functions for 8-node element
!!$
!!$  subroutine get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ)
!!$
!!$  implicit none
!!$
!!$  !include "constants.h"
!!$
!!$  integer myrank,NGNOD,NGLLX,NGLLY,NGLLZ
!!$  integer, parameter :: NDIM=3
!!$
!!$! Gauss-Lobatto-Legendre points of integration
!!$  double precision xigll(NGLLX)
!!$  double precision yigll(NGLLY)
!!$  double precision zigll(NGLLZ)
!!$
!!$! 3D shape functions and their derivatives
!!$  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
!!$  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
!!$
!!$  integer i,j,k,ia
!!$
!!$! location of the nodes of the 3D quadrilateral elements
!!$  double precision xi,eta,gamma
!!$  double precision ra1,ra2,rb1,rb2,rc1,rc2
!!$
!!$! for checking the 3D shape functions
!!$  double precision sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma
!!$
!!$  double precision, parameter :: ONE_EIGHTH = 0.125d0, ZERO = 0.d0, one=1.d0,TINYVAL = 1.d-9
!!$
!!$
!!$! check that the parameter file is correct
!!$  !myrank=0
!!$  if(NGNOD /= 8) call exit_MPI(myrank,'elements should have 8 control nodes')
!!$
!!$! ***
!!$! *** create 3D shape functions and jacobian
!!$! ***
!!$
!!$!--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
!!$
!!$  do i=1,NGLLX
!!$  do j=1,NGLLY
!!$  do k=1,NGLLZ
!!$
!!$  xi = xigll(i)
!!$  eta = yigll(j)
!!$  gamma = zigll(k)
!!$
!!$  ra1 = one + xi
!!$  ra2 = one - xi
!!$
!!$  rb1 = one + eta
!!$  rb2 = one - eta
!!$
!!$  rc1 = one + gamma
!!$  rc2 = one - gamma
!!$
!!$  shape3D(1,i,j,k) = ONE_EIGHTH*ra2*rb2*rc2
!!$  shape3D(2,i,j,k) = ONE_EIGHTH*ra1*rb2*rc2
!!$  shape3D(3,i,j,k) = ONE_EIGHTH*ra1*rb1*rc2
!!$  shape3D(4,i,j,k) = ONE_EIGHTH*ra2*rb1*rc2
!!$  shape3D(5,i,j,k) = ONE_EIGHTH*ra2*rb2*rc1
!!$  shape3D(6,i,j,k) = ONE_EIGHTH*ra1*rb2*rc1
!!$  shape3D(7,i,j,k) = ONE_EIGHTH*ra1*rb1*rc1
!!$  shape3D(8,i,j,k) = ONE_EIGHTH*ra2*rb1*rc1
!!$
!!$  dershape3D(1,1,i,j,k) = - ONE_EIGHTH*rb2*rc2
!!$  dershape3D(1,2,i,j,k) = ONE_EIGHTH*rb2*rc2
!!$  dershape3D(1,3,i,j,k) = ONE_EIGHTH*rb1*rc2
!!$  dershape3D(1,4,i,j,k) = - ONE_EIGHTH*rb1*rc2
!!$  dershape3D(1,5,i,j,k) = - ONE_EIGHTH*rb2*rc1
!!$  dershape3D(1,6,i,j,k) = ONE_EIGHTH*rb2*rc1
!!$  dershape3D(1,7,i,j,k) = ONE_EIGHTH*rb1*rc1
!!$  dershape3D(1,8,i,j,k) = - ONE_EIGHTH*rb1*rc1
!!$
!!$  dershape3D(2,1,i,j,k) = - ONE_EIGHTH*ra2*rc2
!!$  dershape3D(2,2,i,j,k) = - ONE_EIGHTH*ra1*rc2
!!$  dershape3D(2,3,i,j,k) = ONE_EIGHTH*ra1*rc2
!!$  dershape3D(2,4,i,j,k) = ONE_EIGHTH*ra2*rc2
!!$  dershape3D(2,5,i,j,k) = - ONE_EIGHTH*ra2*rc1
!!$  dershape3D(2,6,i,j,k) = - ONE_EIGHTH*ra1*rc1
!!$  dershape3D(2,7,i,j,k) = ONE_EIGHTH*ra1*rc1
!!$  dershape3D(2,8,i,j,k) = ONE_EIGHTH*ra2*rc1
!!$
!!$  dershape3D(3,1,i,j,k) = - ONE_EIGHTH*ra2*rb2
!!$  dershape3D(3,2,i,j,k) = - ONE_EIGHTH*ra1*rb2
!!$  dershape3D(3,3,i,j,k) = - ONE_EIGHTH*ra1*rb1
!!$  dershape3D(3,4,i,j,k) = - ONE_EIGHTH*ra2*rb1
!!$  dershape3D(3,5,i,j,k) = ONE_EIGHTH*ra2*rb2
!!$  dershape3D(3,6,i,j,k) = ONE_EIGHTH*ra1*rb2
!!$  dershape3D(3,7,i,j,k) = ONE_EIGHTH*ra1*rb1
!!$  dershape3D(3,8,i,j,k) = ONE_EIGHTH*ra2*rb1
!!$
!!$  enddo
!!$  enddo
!!$  enddo
!!$
!!$!--- check the shape functions and their derivatives
!!$
!!$  do i=1,NGLLX
!!$  do j=1,NGLLY
!!$  do k=1,NGLLZ
!!$
!!$  sumshape = ZERO
!!$  sumdershapexi = ZERO
!!$  sumdershapeeta = ZERO
!!$  sumdershapegamma = ZERO
!!$
!!$  do ia=1,NGNOD
!!$    sumshape = sumshape + shape3D(ia,i,j,k)
!!$    sumdershapexi = sumdershapexi + dershape3D(1,ia,i,j,k)
!!$    sumdershapeeta = sumdershapeeta + dershape3D(2,ia,i,j,k)
!!$    sumdershapegamma = sumdershapegamma + dershape3D(3,ia,i,j,k)
!!$  enddo
!!$
!!$! sum of shape functions should be one
!!$! sum of derivative of shape functions should be zero
!!$  if(abs(sumshape-one) >  TINYVAL) call exit_MPI(myrank,'error shape functions')
!!$  if(abs(sumdershapexi) >  TINYVAL) call exit_MPI(myrank,'error derivative xi shape functions')
!!$  if(abs(sumdershapeeta) >  TINYVAL) call exit_MPI(myrank,'error derivative eta shape functions')
!!$  if(abs(sumdershapegamma) >  TINYVAL) call exit_MPI(myrank,'error derivative gamma shape functions')
!!$
!!$  enddo
!!$  enddo
!!$  enddo
!!$
!!$  end subroutine get_shape3D !

!
!!! subroutine exit_MPI(myrank,error_msg)

!!!  implicit none
!!!  integer myrank
!!!  character (len=*) error_msg
!!! write(*,*) error_msg
!!!  stop
!!!end subroutine exit_MPI

!
subroutine calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
  implicit none
  integer NGNOD,NGLLX,NGLLY,NGLLZ
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD),xmesh,ymesh,zmesh
  double precision, parameter :: ZERO = 0.d0
  integer ia,i,j,k

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           ! compute mesh coordinates
           xmesh = ZERO
           ymesh = ZERO
           zmesh = ZERO
           do ia=1,NGNOD
              xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
              ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
              zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
           enddo
           xstore(i,j,k) = xmesh
           ystore(i,j,k) = ymesh
           zstore(i,j,k) = zmesh
        enddo
     enddo
  enddo
end subroutine calc_gll_points
