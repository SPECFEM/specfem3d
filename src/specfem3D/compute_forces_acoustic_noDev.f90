!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

  subroutine compute_forces_acoustic_noDev(iphase,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool,deltat, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic)

! computes forces for acoustic elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot
!
  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,TINYVAL_SNGL,ABSORB_USE_PML,ABSORBING_CONDITIONS,PML_CONDITIONS
  use pml_par

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

! time step
  real(kind=CUSTOM_REAL) :: deltat

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

  integer :: iphase
  integer :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic
  integer, dimension(num_phase_ispec_acoustic,2) :: phase_ispec_inner_acoustic

! C-PML absorbing boundary conditions
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer, dimension(nspec2D_xmin) :: ibelm_xmin
  integer, dimension(nspec2D_xmax) :: ibelm_xmax
  integer, dimension(nspec2D_ymin) :: ibelm_ymin
  integer, dimension(nspec2D_ymax) :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM) :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP) :: ibelm_top

! local variables
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
       tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1,temp2,temp3
  real(kind=CUSTOM_REAL) :: temp1l,temp2l,temp3l
  real(kind=CUSTOM_REAL) :: temp1l_new,temp2l_new,temp3l_new

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: dpotentialdxl_new,dpotentialdyl_new,dpotentialdzl_new
  real(kind=CUSTOM_REAL) :: rho_invl

  integer :: ispec,ispec2D,iglob,i,j,k,l,ispec_p,num_elements

  ! local C-PML absorbing boundary conditions parameters
  integer :: ispec_CPML

  if( iphase == 1 ) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! gets values for element
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          chi_elem(i,j,k) = potential_acoustic(ibool(i,j,k,ispec))
        enddo
      enddo
    enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

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

          if (PML_CONDITIONS) then
             ! do not merge this second line with the first using an ".and." statement
             ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
             if(CPML_mask_ibool(ispec)) then
                temp1l_new = temp1l
                temp2l_new = temp2l
                temp3l_new = temp3l

                do l=1,NGLLX
                   hp1 = hprime_xx(l,i)
                   iglob = ibool(l,j,k,ispec)
                   temp1l_new = temp1l_new + deltat*potential_dot_acoustic(iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
                   hp2 = hprime_yy(l,j)
                   iglob = ibool(i,l,k,ispec)
                   temp2l_new = temp2l_new + deltat*potential_dot_acoustic(iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
                   hp3 = hprime_zz(l,k)
                   iglob = ibool(i,j,l,ispec)
                   temp3l_new = temp3l_new + deltat*potential_dot_acoustic(iglob)*hp3
                enddo
             endif
          endif

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

          ! stores derivatives of ux, uy and uz with respect to x, y and z
          if (PML_CONDITIONS) then
             ! do not merge this second line with the first using an ".and." statement
             ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
             if(CPML_mask_ibool(ispec)) then
                ispec_CPML = spec_to_CPML(ispec)

                PML_dpotential_dxl(i,j,k,ispec_CPML) = dpotentialdxl
                PML_dpotential_dyl(i,j,k,ispec_CPML) = dpotentialdyl
                PML_dpotential_dzl(i,j,k,ispec_CPML) = dpotentialdzl

                dpotentialdxl_new = xixl*temp1l_new + etaxl*temp2l_new + gammaxl*temp3l_new
                dpotentialdyl_new = xiyl*temp1l_new + etayl*temp2l_new + gammayl*temp3l_new
                dpotentialdzl_new = xizl*temp1l_new + etazl*temp2l_new + gammazl*temp3l_new

                PML_dpotential_dxl_new(i,j,k,ispec_CPML) = dpotentialdxl_new
                PML_dpotential_dyl_new(i,j,k,ispec_CPML) = dpotentialdyl_new
                PML_dpotential_dzl_new(i,j,k,ispec_CPML) = dpotentialdzl_new
             endif
          endif

          ! density (reciproc)
          rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

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

    if (PML_CONDITIONS) then
       ! do not merge this second line with the first using an ".and." statement
       ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
       if(CPML_mask_ibool(ispec)) then
          ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
          call pml_compute_memory_variables(ispec,ispec_CPML,deltat,jacobianl,tempx1,tempy1,tempz1,tempx2,tempy2,tempz2, &
               tempx3,tempy3,tempz3,sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz, &
               sigma_yx,sigma_zx,sigma_zy,lambdal,mul,lambdalplus2mul,xixl,xiyl,xizl, &
               etaxl,etayl,etazl,gammaxl,gammayl,gammazl)

          ! calculates contribution from each C-PML element to update acceleration
          call pml_compute_accel_contribution(ispec,ispec_CPML,deltat,jacobianl,accel_elastic_CPML)
       endif
    endif

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

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - ( temp1l + temp2l + temp3l )

          ! updates potential_dot_dot_acoustic with contribution from each C-PML element
          if (PML_CONDITIONS) then
             ! do not merge this second line with the first using an ".and." statement
             ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
             if(CPML_mask_ibool(ispec)) then
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                     potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML)
             endif
          endif

        enddo
      enddo
    enddo

  enddo ! end of loop over all spectral elements

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

  end subroutine compute_forces_acoustic_noDev

