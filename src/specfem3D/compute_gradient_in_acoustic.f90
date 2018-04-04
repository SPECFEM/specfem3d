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

  subroutine compute_gradient_in_acoustic(ispec, &
                        scalar_field, vector_field_element)

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

  use constants
  use specfem_par, only: NGLOB_AB,xix,xiy,xiz,etax,etay,etaz, &
                          gammax,gammay,gammaz,xix_regular,irregular_element_number, &
                          ibool,rhostore,GRAVITY,hprime_xx,hprime_yy,hprime_zz

  implicit none

  integer,intent(in) :: ispec
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: scalar_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out) :: vector_field_element

! local parameters
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) temp1l,temp2l,temp3l
  real(kind=CUSTOM_REAL) rho_invl
  integer :: i,j,k,l,ispec_irreg

! double loop over GLL points to compute and store gradients
  vector_field_element(:,:,:,:) = 0._CUSTOM_REAL

  ispec_irreg = irregular_element_number(ispec)

  do k= 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        ! derivative along x
        temp1l = ZERO
        do l = 1,NGLLX
          temp1l = temp1l + scalar_field(ibool(l,j,k,ispec))*hprime_xx(i,l)
        enddo

        ! derivative along y
        temp2l = ZERO
        do l = 1,NGLLZ
          temp2l = temp2l + scalar_field(ibool(i,l,k,ispec))*hprime_yy(j,l)
        enddo

        ! derivative along z
        temp3l = ZERO
        do l = 1,NGLLZ
          temp3l = temp3l + scalar_field(ibool(i,j,l,ispec))*hprime_zz(k,l)
        enddo

        ! daniel: TODO - check gravity case here
        if (GRAVITY) then
          rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
        else
          rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)
        endif

        if (ispec_irreg /= 0 ) then !irregular element

          xixl = xix(i,j,k,ispec_irreg)
          xiyl = xiy(i,j,k,ispec_irreg)
          xizl = xiz(i,j,k,ispec_irreg)
          etaxl = etax(i,j,k,ispec_irreg)
          etayl = etay(i,j,k,ispec_irreg)
          etazl = etaz(i,j,k,ispec_irreg)
          gammaxl = gammax(i,j,k,ispec_irreg)
          gammayl = gammay(i,j,k,ispec_irreg)
          gammazl = gammaz(i,j,k,ispec_irreg)

          ! derivatives of acoustic scalar potential field on GLL points
          vector_field_element(1,i,j,k) = (temp1l*xixl + temp2l*etaxl + temp3l*gammaxl) * rho_invl
          vector_field_element(2,i,j,k) = (temp1l*xiyl + temp2l*etayl + temp3l*gammayl) * rho_invl
          vector_field_element(3,i,j,k) = (temp1l*xizl + temp2l*etazl + temp3l*gammazl) * rho_invl

        else !regular element

          ! derivatives of acoustic scalar potential field on GLL points
          vector_field_element(1,i,j,k) = temp1l * xix_regular * rho_invl
          vector_field_element(2,i,j,k) = temp2l * xix_regular * rho_invl
          vector_field_element(3,i,j,k) = temp3l * xix_regular * rho_invl

        endif
      enddo
    enddo
  enddo

  end subroutine compute_gradient_in_acoustic


