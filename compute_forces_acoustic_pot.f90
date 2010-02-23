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

! for acoustic solver

  subroutine compute_forces_acoustic_pot( iphase, NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic )

! computes forces for acoustic elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot  
!
  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,TINYVAL_SNGL
  use PML_par,only:PML,ispec_is_PML_inum
  implicit none
  !include "constants.h"
  integer :: NSPEC_AB,NGLOB_AB

! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: &
        potential_acoustic,potential_dot_dot_acoustic

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        rhostore,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! communication overlap
!  logical, dimension(NSPEC_AB) :: ispec_is_inner
!  logical :: phase_is_inner
  
!  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

  integer :: iphase
  integer :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic
  integer, dimension(num_phase_ispec_acoustic,2) :: phase_ispec_inner_acoustic

! local variables
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1,temp2,temp3
  real(kind=CUSTOM_REAL) temp1l,temp2l,temp3l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) rho_invl
  
  integer :: ispec,iglob,i,j,k,l,ispec_p,num_elements

  if( iphase == 1 ) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

! loop over spectral elements
  do ispec_p = 1,num_elements

    !if ( (ispec_is_inner(ispec) .eqv. phase_is_inner) ) then

      ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

      ! only elements outside PML, inside "regular" domain
      if( PML ) then
        if( ispec_is_PML_inum(ispec) > 0 ) then
         cycle
        endif
      endif
      
!      if( ispec_is_acoustic(ispec) ) then

        ! gets values for element
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              chi_elem(i,j,k) = potential_acoustic(ibool(i,j,k,ispec))
            enddo
          enddo
        enddo
        ! would check if anything to do, but might lower accuracy of computation
        !if( maxval( abs( chi_elem ) ) < TINYVAL_SNGL ) cycle

        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX

              ! density (reciproc)
              rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec) 
              
              ! derivative along x, y, z
              ! first double loop over GLL points to compute and store gradients
              ! we can merge the loops because NGLLX == NGLLY == NGLLZ
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              temp3l = 0._CUSTOM_REAL
              do l = 1,NGLLX
                temp1l = temp1l + chi_elem(l,j,k)*hprime_xx(i,l)
                temp2l = temp2l + chi_elem(i,l,k)*hprime_yy(j,l)
                temp3l = temp3l + chi_elem(i,j,l)*hprime_zz(k,l)
              enddo 

              ! get derivatives of potential with respect to x, y and z
              xixl = xix(i,j,k,ispec)
              xiyl = xiy(i,j,k,ispec)
              xizl = xiz(i,j,k,ispec)
              etaxl = etax(i,j,k,ispec)
              etayl = etay(i,j,k,ispec)
              etazl = etaz(i,j,k,ispec)
              gammaxl = gammax(i,j,k,ispec)
              gammayl = gammay(i,j,k,ispec)
              gammazl = gammaz(i,j,k,ispec)
              jacobianl = jacobian(i,j,k,ispec)

              ! derivatives of potential
              dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l
              dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l
              dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l

              ! for acoustic medium
              ! also add GLL integration weights
              temp1(i,j,k) = rho_invl * wgllwgll_yz(j,k) * jacobianl* &
                            (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
              temp2(i,j,k) = rho_invl * wgllwgll_xz(i,k) * jacobianl* &
                            (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
              temp3(i,j,k) = rho_invl * wgllwgll_xy(i,j) * jacobianl* &
                            (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)
            enddo
          enddo
        enddo

        ! second double-loop over GLL to compute all the terms
        do k = 1,NGLLZ
          do j = 1,NGLLZ
            do i = 1,NGLLX

              ! along x,y,z direction
              ! and assemble the contributions
              !!! can merge these loops because NGLLX = NGLLY = NGLLZ   
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              temp3l = 0._CUSTOM_REAL
              do l=1,NGLLX
                temp1l = temp1l + temp1(l,j,k) * hprimewgll_xx(l,i)
                temp2l = temp2l + temp2(i,l,k) * hprimewgll_yy(l,j)
                temp3l = temp3l + temp3(i,j,l) * hprimewgll_zz(l,k)
              enddo

              ! sum contributions from each element to the global values              
              iglob = ibool(i,j,k,ispec)
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                  - ( temp1l + temp2l + temp3l )

            enddo
          enddo 
        enddo

!      endif ! end of test if acoustic element
!    endif ! ispec_is_inner
    
  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_acoustic_pot

