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

  subroutine compute_gradient_in_acoustic(ispec,scalar_field,vector_field_element)

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  use specfem_par, only: NGLOB_AB

  implicit none

  integer,intent(in) :: ispec
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: scalar_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out) :: vector_field_element

  ! no need to test the two others because NGLLX == NGLLY = NGLLZ in unstructured meshes
  if (NGLLX == 5 .or. NGLLX == 6 .or. NGLLX == 7) then
    call compute_gradient_in_acoustic_fast_Deville(ispec,scalar_field,vector_field_element)
  else
    call compute_gradient_in_acoustic_generic_slow(ispec,scalar_field,vector_field_element)
  endif

  end subroutine compute_gradient_in_acoustic

!
!------------------------------------------------------
!

  subroutine compute_gradient_in_acoustic_generic_slow(ispec,scalar_field,vector_field_element)

! calculates gradient of given acoustic scalar (potential) field on all GLL points in one, single element
! note:
!   displacement s = (rho)^{-1} \del \chi
!   velocity     v = (rho)^{-1} \del \ddot \chi
!
!  in case of gravity:
!   displacement s = \del \chi
!   velocity     v = \del \ddot \chi
! returns: (1/rho) times gradient vector field (vector_field_element) in specified element
!             or in gravity case, just gradient vector field

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  use specfem_par, only: NGLOB_AB,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,xix_regular,irregular_element_number, &
                         ibool,rhostore,hprime_xx,hprime_yy,hprime_zz ! ,GRAVITY

  implicit none

  integer,intent(in) :: ispec
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: scalar_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out) :: vector_field_element

! local parameters
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: temp1l,temp2l,temp3l
  real(kind=CUSTOM_REAL) :: rho_invl
  integer :: i,j,k,l,ispec_irreg

! double loop over GLL points to compute and store gradients
  vector_field_element(:,:,:,:) = 0._CUSTOM_REAL

  ispec_irreg = irregular_element_number(ispec)

  if (ispec_irreg /= 0) then
    ! irregular element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          ! derivative along x
          temp1l = 0._CUSTOM_REAL
          do l = 1,NGLLX
            temp1l = temp1l + scalar_field(ibool(l,j,k,ispec))*hprime_xx(i,l)
          enddo

          ! derivative along y
          temp2l = 0._CUSTOM_REAL
          do l = 1,NGLLY
            temp2l = temp2l + scalar_field(ibool(i,l,k,ispec))*hprime_yy(j,l)
          enddo

          ! derivative along z
          temp3l = 0._CUSTOM_REAL
          do l = 1,NGLLZ
            temp3l = temp3l + scalar_field(ibool(i,j,l,ispec))*hprime_zz(k,l)
          enddo

          ! Daniel Peter: TODO - check gravity case here
!! DK DK NEVER put an "if" statement inside a critical loop, since that prevents vectorization
!! DK DK and thus drastically slows down the code; put it outside the loop, and duplicate the content of the loop
          !if (GRAVITY) then
          ! rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !else
            rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !endif

          xixl = xixstore(i,j,k,ispec_irreg)
          xiyl = xiystore(i,j,k,ispec_irreg)
          xizl = xizstore(i,j,k,ispec_irreg)
          etaxl = etaxstore(i,j,k,ispec_irreg)
          etayl = etaystore(i,j,k,ispec_irreg)
          etazl = etazstore(i,j,k,ispec_irreg)
          gammaxl = gammaxstore(i,j,k,ispec_irreg)
          gammayl = gammaystore(i,j,k,ispec_irreg)
          gammazl = gammazstore(i,j,k,ispec_irreg)

          ! derivatives of acoustic scalar potential field on GLL points
          vector_field_element(1,i,j,k) = (temp1l*xixl + temp2l*etaxl + temp3l*gammaxl) * rho_invl
          vector_field_element(2,i,j,k) = (temp1l*xiyl + temp2l*etayl + temp3l*gammayl) * rho_invl
          vector_field_element(3,i,j,k) = (temp1l*xizl + temp2l*etazl + temp3l*gammazl) * rho_invl

        enddo
      enddo
    enddo

  else
    ! regular element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          ! derivative along x
          temp1l = 0._CUSTOM_REAL
          do l = 1,NGLLX
            temp1l = temp1l + scalar_field(ibool(l,j,k,ispec))*hprime_xx(i,l)
          enddo

          ! derivative along y
          temp2l = 0._CUSTOM_REAL
          do l = 1,NGLLY
            temp2l = temp2l + scalar_field(ibool(i,l,k,ispec))*hprime_yy(j,l)
          enddo

          ! derivative along z
          temp3l = 0._CUSTOM_REAL
          do l = 1,NGLLZ
            temp3l = temp3l + scalar_field(ibool(i,j,l,ispec))*hprime_zz(k,l)
          enddo

          ! Daniel Peter: TODO - check gravity case here
!! DK DK NEVER put an "if" statement inside a critical loop, since that prevents vectorization
!! DK DK and thus drastically slows down the code; put it outside the loop, and duplicate the content of the loop
          !if (GRAVITY) then
          ! rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !else
            rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !endif

          ! derivatives of acoustic scalar potential field on GLL points
          vector_field_element(1,i,j,k) = temp1l * xix_regular * rho_invl
          vector_field_element(2,i,j,k) = temp2l * xix_regular * rho_invl
          vector_field_element(3,i,j,k) = temp3l * xix_regular * rho_invl

        enddo
      enddo
    enddo

  endif

  end subroutine compute_gradient_in_acoustic_generic_slow

!
!------------------------------------------------------
!

  subroutine compute_gradient_in_acoustic_fast_Deville(ispec,scalar_field,vector_field_element)

! calculates gradient of given acoustic scalar (potential) field on all GLL points in one, single element
! note:
!   displacement s = (rho)^{-1} \del \chi
!   velocity     v = (rho)^{-1} \del \ddot \chi
!
!  in case of gravity:
!   displacement s = \del \chi
!   velocity     v = \del \ddot \chi
! returns: (1/rho) times gradient vector field (vector_field_element) in specified element
!             or in gravity case, just gradient vector field

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,m1,m2
  use specfem_par, only: NGLOB_AB,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,xix_regular,irregular_element_number, &
                         ibool,rhostore,hprime_xx,hprime_xxT ! ,GRAVITY

  implicit none

  integer,intent(in) :: ispec
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: scalar_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out) :: vector_field_element

  ! local parameters
  ! note: declaring arrays in this subroutine here will allocate them generally on the stack
  !       (intel by default; not for gfortran though, it always uses heap memory).
  !       stack memory access is faster, thus please let these declarations here for local element arrays...
  !
  ! arrays for elemental computations inside a given spectral element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: field_local
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1,temp2,temp3

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: rho_invl
  integer :: i,j,k,ispec_irreg

  ! double loop over GLL points to compute and store gradients
  vector_field_element(:,:,:,:) = 0._CUSTOM_REAL

  ! gets value of the field inside the element and make it local
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        field_local(i,j,k) = scalar_field(ibool(i,j,k,ispec))
      enddo
    enddo
  enddo

  ! derivatives along x, y and z

  ! subroutines adapted from Deville, Fischer and Mund, High-order methods
  ! for incompressible fluid flow, Cambridge University Press (2002),
  ! pages 386 and 389 and Figure 8.3.1

  ! computes 1. matrix multiplication for temp1
  ! computes 2. matrix multiplication for temp2
  ! computes 3. matrix multiplication for temp3
  if (NGLLX == 5) then
    call mxm5_single(hprime_xx,m1,field_local,temp1,m2)
    call mxm5_3dmat_single(field_local,m1,hprime_xxT,m1,temp2,NGLLX)
    call mxm5_single(field_local,m2,hprime_xxT,temp3,m1)
  else if (NGLLX == 6) then
    call mxm6_single(hprime_xx,m1,field_local,temp1,m2)
    call mxm6_3dmat_single(field_local,m1,hprime_xxT,m1,temp2,NGLLX)
    call mxm6_single(field_local,m2,hprime_xxT,temp3,m1)
  else if (NGLLX == 7) then
    call mxm7_single(hprime_xx,m1,field_local,temp1,m2)
    call mxm7_3dmat_single(field_local,m1,hprime_xxT,m1,temp2,NGLLX)
    call mxm7_single(field_local,m2,hprime_xxT,temp3,m1)
  endif

  ispec_irreg = irregular_element_number(ispec)

  if (ispec_irreg /= 0) then
    ! irregular element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          ! Daniel Peter: TODO - check gravity case here
!! DK DK NEVER put an "if" statement inside a critical loop, since that prevents vectorization
!! DK DK and thus drastically slows down the code; put it outside the loop, and duplicate the content of the loop
          !if (GRAVITY) then
          ! rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !else
            rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !endif

          xixl = xixstore(i,j,k,ispec_irreg)
          xiyl = xiystore(i,j,k,ispec_irreg)
          xizl = xizstore(i,j,k,ispec_irreg)
          etaxl = etaxstore(i,j,k,ispec_irreg)
          etayl = etaystore(i,j,k,ispec_irreg)
          etazl = etazstore(i,j,k,ispec_irreg)
          gammaxl = gammaxstore(i,j,k,ispec_irreg)
          gammayl = gammaystore(i,j,k,ispec_irreg)
          gammazl = gammazstore(i,j,k,ispec_irreg)

          ! derivatives of acoustic scalar potential field on GLL points
          vector_field_element(1,i,j,k) = (temp1(i,j,k)*xixl + temp2(i,j,k)*etaxl + temp3(i,j,k)*gammaxl) * rho_invl
          vector_field_element(2,i,j,k) = (temp1(i,j,k)*xiyl + temp2(i,j,k)*etayl + temp3(i,j,k)*gammayl) * rho_invl
          vector_field_element(3,i,j,k) = (temp1(i,j,k)*xizl + temp2(i,j,k)*etazl + temp3(i,j,k)*gammazl) * rho_invl

        enddo
      enddo
    enddo

  else
    ! regular element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          ! Daniel Peter: TODO - check gravity case here
!! DK DK NEVER put an "if" statement inside a critical loop, since that prevents vectorization
!! DK DK and thus drastically slows down the code; put it outside the loop, and duplicate the content of the loop
          !if (GRAVITY) then
          ! rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !else
            rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
          !endif

          ! derivatives of acoustic scalar potential field on GLL points
          vector_field_element(1,i,j,k) = temp1(i,j,k) * xix_regular * rho_invl
          vector_field_element(2,i,j,k) = temp2(i,j,k) * xix_regular * rho_invl
          vector_field_element(3,i,j,k) = temp3(i,j,k) * xix_regular * rho_invl

        enddo
      enddo
    enddo

  endif

  contains

!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices (5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications is in general slower
!
! please leave the routines here to help compilers inline the code

  subroutine mxm5_single(A,n1,B,C,n3)

! two-dimensional arrays (25,5)/(5,25)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_single

  !-------------

  subroutine mxm6_single(A,n1,B,C,n3)

! two-dimensional arrays (36,6)/(6,36)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,6),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(6,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j) &
              + A(i,6) * B(6,j)
    enddo
  enddo

  end subroutine mxm6_single

  !-------------

  subroutine mxm7_single(A,n1,B,C,n3)

! two-dimensional arrays (49,7)/(7,49)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,7),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(7,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j) &
              + A(i,6) * B(6,j) &
              + A(i,7) * B(7,j)
    enddo
  enddo

  end subroutine mxm7_single

!--------------------------------------------------------------------------------------------

  subroutine mxm5_3dmat_single(A,n1,B,n2,C,n3)

! three-dimensional arrays (5,5,5) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3dmat_single

  !-------------

  subroutine mxm6_3dmat_single(A,n1,B,n2,C,n3)

! three-dimensional arrays (6,6,6) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,6,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(6,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j) &
                  + A(i,6,k) * B(6,j)
      enddo
    enddo
  enddo

  end subroutine mxm6_3dmat_single

  !-------------

  subroutine mxm7_3dmat_single(A,n1,B,n2,C,n3)

! three-dimensional arrays (7,7,7) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,7,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(7,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j) &
                  + A(i,6,k) * B(6,j) &
                  + A(i,7,k) * B(7,j)
      enddo
    enddo
  enddo

  end subroutine mxm7_3dmat_single

  end subroutine compute_gradient_in_acoustic_fast_Deville

