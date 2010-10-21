!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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


subroutine compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                        scalar_field, vector_field_element,&
                        hprime_xx,hprime_yy,hprime_zz, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        ibool,rhostore)

! calculates gradient of given acoustic scalar (potential) field on all GLL points in one, single element
! note: 
!   displacement s = (rho)^{-1} \del \chi
!   velocity          v = (rho)^{-1} \del \ddot \chi
!
! returns: (1/rho) times gradient vector field (vector_field_element) in specified element 

  implicit none
  include 'constants.h'

  integer,intent(in) :: ispec,NSPEC_AB,NGLOB_AB
  
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: scalar_field
  
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out) :: vector_field_element
  
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB):: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhostore

! array with derivatives of Lagrange polynomials 
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! local parameters  
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) temp1l,temp2l,temp3l 
  real(kind=CUSTOM_REAL) rho_invl  
  integer :: i,j,k,l

! double loop over GLL points to compute and store gradients
  vector_field_element(:,:,:,:) = 0._CUSTOM_REAL
  
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
        
        xixl = xix(i,j,k,ispec)
        xiyl = xiy(i,j,k,ispec)
        xizl = xiz(i,j,k,ispec)
        etaxl = etax(i,j,k,ispec)
        etayl = etay(i,j,k,ispec)
        etazl = etaz(i,j,k,ispec)
        gammaxl = gammax(i,j,k,ispec)
        gammayl = gammay(i,j,k,ispec)
        gammazl = gammaz(i,j,k,ispec)
        
        rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)                              
        
        ! derivatives of acoustic scalar potential field on GLL points
        vector_field_element(1,i,j,k) = (temp1l*xixl + temp2l*etaxl + temp3l*gammaxl) * rho_invl
        vector_field_element(2,i,j,k) = (temp1l*xiyl + temp2l*etayl + temp3l*gammayl) * rho_invl
        vector_field_element(3,i,j,k) = (temp1l*xizl + temp2l*etazl + temp3l*gammazl) * rho_invl
                
      enddo
    enddo
  enddo

end subroutine compute_gradient


