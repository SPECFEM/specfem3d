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

  subroutine compute_forces_acoustic_Dev(iphase,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic,backward_simulation)

! computes forces for acoustic elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot
!
  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,TINYVAL_SNGL, &
                        m1,m2,NGLLCUBE,PML_CONDITIONS
  use pml_par, only: is_CPML, spec_to_CPML, NSPEC_CPML, &
                     PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl,&
                     PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old,&
                     potential_dot_dot_acoustic_CPML,rmemory_dpotential_dxl,rmemory_dpotential_dyl,&
                     rmemory_dpotential_dzl,rmemory_potential_acoustic,potential_acoustic_old

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        rhostore,jacobian

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprimewgll_xx,hprimewgll_xxT

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

  integer :: iphase
  integer :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic
  integer, dimension(num_phase_ispec_acoustic,2) :: phase_ispec_inner_acoustic

! CPML adjoint
  logical :: backward_simulation

  integer :: ispec_CPML

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_old,tempx2_old,tempx3_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempx2,newtempx3

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: rho_invl

  integer :: ispec,iglob,i,j,k,ispec_p,num_elements

  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem_old

  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points_old
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points_old
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: E1_m1_m2_5points

  equivalence(chi_elem,B1_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(chi_elem_old,B1_m1_m2_5points_old)
  equivalence(tempx1_old,C1_m1_m2_5points_old)
  equivalence(newtempx1,E1_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points_old
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points_old
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: E1_mxm_m2_m1_5points

  equivalence(chi_elem,A1_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(chi_elem_old,A1_mxm_m2_m1_5points_old)
  equivalence(tempx3_old,C1_mxm_m2_m1_5points_old)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

  if( iphase == 1 ) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! gets values for element
#ifndef FORCE_VECTORIZATION
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          chi_elem(i,j,k) = potential_acoustic(ibool(i,j,k,ispec))
        enddo
      enddo
    enddo
#else
! this will (purposely) give out-of-bound array accesses if run through range checking,
! thus use only for production runs with no bound checking
    do ijk = 1,NGLLCUBE
      chi_elem(ijk,1,1) = potential_acoustic(ibool(ijk,1,1,ispec))
    enddo
#endif

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    do j=1,m2
      do i=1,m1
        C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
                                hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
                                hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
                                hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
                                hprime_xx(i,5)*B1_m1_m2_5points(5,j)
      enddo
    enddo

    do k = 1,NGLLX
      do j=1,m1
        do i=1,m1
          tempx2(i,j,k) = chi_elem(i,1,k)*hprime_xxT(1,j) + &
                          chi_elem(i,2,k)*hprime_xxT(2,j) + &
                          chi_elem(i,3,k)*hprime_xxT(3,j) + &
                          chi_elem(i,4,k)*hprime_xxT(4,j) + &
                          chi_elem(i,5,k)*hprime_xxT(5,j)
        enddo
      enddo
    enddo

    do j=1,m1
      do i=1,m2
        C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                    A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                    A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                    A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                    A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
      enddo
    enddo

    if(PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an ".and." statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false

      if(is_CPML(ispec)) then
      ! gets values for element
#ifndef FORCE_VECTORIZATION
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              chi_elem_old(i,j,k) = potential_acoustic_old(ibool(i,j,k,ispec))
            enddo
          enddo
        enddo
#else
        ! this will (purposely) give out-of-bound array accesses if run through range checking,
        ! thus use only for production runs with no bound checking
        do ijk = 1,NGLLCUBE
          chi_elem_old(ijk,1,1) = potential_acoustic_old(ibool(ijk,1,1,ispec))
        enddo
#endif

        ! subroutines adapted from Deville, Fischer and Mund, High-order methods
        ! for incompressible fluid flow, Cambridge University Press (2002),
        ! pages 386 and 389 and Figure 8.3.1
        do j=1,m2
          do i=1,m1
            C1_m1_m2_5points_old(i,j) = hprime_xx(i,1)*B1_m1_m2_5points_old(1,j) + &
                                        hprime_xx(i,2)*B1_m1_m2_5points_old(2,j) + &
                                        hprime_xx(i,3)*B1_m1_m2_5points_old(3,j) + &
                                        hprime_xx(i,4)*B1_m1_m2_5points_old(4,j) + &
                                        hprime_xx(i,5)*B1_m1_m2_5points_old(5,j)
          enddo
        enddo

        do k = 1,NGLLX
          do j=1,m1
            do i=1,m1
              tempx2_old(i,j,k) = chi_elem_old(i,1,k)*hprime_xxT(1,j) + &
                                  chi_elem_old(i,2,k)*hprime_xxT(2,j) + &
                                  chi_elem_old(i,3,k)*hprime_xxT(3,j) + &
                                  chi_elem_old(i,4,k)*hprime_xxT(4,j) + &
                                  chi_elem_old(i,5,k)*hprime_xxT(5,j)
            enddo
          enddo
        enddo

        do j=1,m1
          do i=1,m2
            C1_mxm_m2_m1_5points_old(i,j) = A1_mxm_m2_m1_5points_old(i,1)*hprime_xxT(1,j) + &
                                            A1_mxm_m2_m1_5points_old(i,2)*hprime_xxT(2,j) + &
                                            A1_mxm_m2_m1_5points_old(i,3)*hprime_xxT(3,j) + &
                                            A1_mxm_m2_m1_5points_old(i,4)*hprime_xxT(4,j) + &
                                            A1_mxm_m2_m1_5points_old(i,5)*hprime_xxT(5,j)
          enddo
        enddo

#ifndef FORCE_VECTORIZATION
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX

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
              PML_dpotential_dxl(i,j,k) = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              PML_dpotential_dyl(i,j,k) = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              PML_dpotential_dzl(i,j,k) = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              PML_dpotential_dxl_old(i,j,k) = xixl*tempx1_old(i,j,k) + etaxl*tempx2_old(i,j,k) + gammaxl*tempx3_old(i,j,k)
              PML_dpotential_dyl_old(i,j,k) = xiyl*tempx1_old(i,j,k) + etayl*tempx2_old(i,j,k) + gammayl*tempx3_old(i,j,k)
              PML_dpotential_dzl_old(i,j,k) = xizl*tempx1_old(i,j,k) + etazl*tempx2_old(i,j,k) + gammazl*tempx3_old(i,j,k)

            enddo
          enddo
        enddo
#else
        do ijk = 1,NGLLCUBE
           ! get derivatives of potential with respect to x, y and z
           xixl = xix(ijk,1,1,ispec)
           xiyl = xiy(ijk,1,1,ispec)
           xizl = xiz(ijk,1,1,ispec)
           etaxl = etax(ijk,1,1,ispec)
           etayl = etay(ijk,1,1,ispec)
           etazl = etaz(ijk,1,1,ispec)
           gammaxl = gammax(ijk,1,1,ispec)
           gammayl = gammay(ijk,1,1,ispec)
           gammazl = gammaz(ijk,1,1,ispec)
           jacobianl = jacobian(ijk,1,1,ispec)

           ! derivatives of potential
           PML_dpotential_dxl(ijk,1,1) = xixl*tempx1(ijk,1,1) + etaxl*tempx2(ijk,1,1) + gammaxl*tempx3(ijk,1,1)
           PML_dpotential_dyl(ijk,1,1) = xiyl*tempx1(ijk,1,1) + etayl*tempx2(ijk,1,1) + gammayl*tempx3(ijk,1,1)
           PML_dpotential_dzl(ijk,1,1) = xizl*tempx1(ijk,1,1) + etazl*tempx2(ijk,1,1) + gammazl*tempx3(ijk,1,1)

           PML_dpotential_dxl_old(ijk,1,1) = xixl*tempx1_old(ijk,1,1) + etaxl*tempx2_old(ijk,1,1) + gammaxl*tempx3_old(ijk,1,1)
           PML_dpotential_dyl_old(ijk,1,1) = xiyl*tempx1_old(ijk,1,1) + etayl*tempx2_old(ijk,1,1) + gammayl*tempx3_old(ijk,1,1)
           PML_dpotential_dzl_old(ijk,1,1) = xizl*tempx1_old(ijk,1,1) + etazl*tempx2_old(ijk,1,1) + gammazl*tempx3_old(ijk,1,1)

        enddo
#endif

      endif
    endif

#ifndef FORCE_VECTORIZATION
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

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
          dpotentialdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
          dpotentialdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
          dpotentialdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

          ! density (reciproc)
          rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

          ! for acoustic medium
          ! also add GLL integration weights
          tempx1(i,j,k) = rho_invl * jacobianl* &
                        (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
          tempx2(i,j,k) = rho_invl * jacobianl* &
                        (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
          tempx3(i,j,k) = rho_invl * jacobianl* &
                        (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)
        enddo
      enddo
    enddo
#else
    do ijk = 1,NGLLCUBE
         ! get derivatives of potential with respect to x, y and z
          xixl = xix(ijk,1,1,ispec)
          xiyl = xiy(ijk,1,1,ispec)
          xizl = xiz(ijk,1,1,ispec)
          etaxl = etax(ijk,1,1,ispec)
          etayl = etay(ijk,1,1,ispec)
          etazl = etaz(ijk,1,1,ispec)
          gammaxl = gammax(ijk,1,1,ispec)
          gammayl = gammay(ijk,1,1,ispec)
          gammazl = gammaz(ijk,1,1,ispec)
          jacobianl = jacobian(ijk,1,1,ispec)

          ! derivatives of potential
          dpotentialdxl = xixl*tempx1(ijk,1,1) + etaxl*tempx2(ijk,1,1) + gammaxl*tempx3(ijk,1,1)
          dpotentialdyl = xiyl*tempx1(ijk,1,1) + etayl*tempx2(ijk,1,1) + gammayl*tempx3(ijk,1,1)
          dpotentialdzl = xizl*tempx1(ijk,1,1) + etazl*tempx2(ijk,1,1) + gammazl*tempx3(ijk,1,1)

          ! density (reciproc)
          rho_invl = 1.0_CUSTOM_REAL / rhostore(ijk,1,1,ispec)

          ! for acoustic medium
          ! also add GLL integration weights
          tempx1(ijk,1,1) = rho_invl * jacobianl* &
                        (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
          tempx2(ijk,1,1) = rho_invl * jacobianl* &
                        (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
          tempx3(ijk,1,1) = rho_invl * jacobianl* &
                        (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)
    enddo
#endif

    if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
       ! do not merge this second line with the first using an ".and." statement
       ! because array is_CPML() is unallocated when PML_CONDITIONS is false
       if(is_CPML(ispec)) then
          tempx1 = 0._CUSTOM_REAL; tempx2 = 0._CUSTOM_REAL; tempx3 = 0._CUSTOM_REAL
          ispec_CPML = spec_to_CPML(ispec)
          ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
          call pml_compute_memory_variables_acoustic(ispec,ispec_CPML,tempx1,tempx2,tempx3,&
                                                     rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl)

          ! calculates contribution from each C-PML element to update acceleration
          call pml_compute_accel_contribution_acoustic(ispec,ispec_CPML,potential_acoustic,&
                                                       potential_dot_acoustic,rmemory_potential_acoustic)
       endif
    endif

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    do j=1,m2
      do i=1,m1
        E1_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C1_m1_m2_5points(1,j) + &
                                hprimewgll_xxT(i,2)*C1_m1_m2_5points(2,j) + &
                                hprimewgll_xxT(i,3)*C1_m1_m2_5points(3,j) + &
                                hprimewgll_xxT(i,4)*C1_m1_m2_5points(4,j) + &
                                hprimewgll_xxT(i,5)*C1_m1_m2_5points(5,j)
      enddo
    enddo

    do k = 1,NGLLX
      do j=1,m1
        do i=1,m1
          newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
                             tempx2(i,2,k)*hprimewgll_xx(2,j) + &
                             tempx2(i,3,k)*hprimewgll_xx(3,j) + &
                             tempx2(i,4,k)*hprimewgll_xx(4,j) + &
                             tempx2(i,5,k)*hprimewgll_xx(5,j)
        enddo
      enddo
    enddo

    do j=1,m1
      do i=1,m2
        E1_mxm_m2_m1_5points(i,j) = C1_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                    C1_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                    C1_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                    C1_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                    C1_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
      enddo
    enddo

    ! second double-loop over GLL to compute all the terms
#ifndef FORCE_VECTORIZATION
    do k = 1,NGLLZ
      do j = 1,NGLLZ
        do i = 1,NGLLX

          ! sum contributions from each element to the global values
          iglob = ibool(i,j,k,ispec)
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - (wgllwgll_yz_3D(i,j,k)*newtempx1(i,j,k) &
                       + wgllwgll_xz_3D(i,j,k)*newtempx2(i,j,k) + wgllwgll_xy_3D(i,j,k)*newtempx3(i,j,k))

        enddo
      enddo
    enddo
#else
! we can force vectorization using a compiler directive here because we know that there is no dependency
! inside a given spectral element, since all the global points of a local elements are different by definition
! (only common points between different elements can be the same)
!DIR$ IVDEP
    do ijk = 1,NGLLCUBE
          ! sum contributions from each element to the global values
       iglob = ibool(ijk,1,1,ispec)
       potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - (wgllwgll_yz_3D(ijk,1,1)*newtempx1(ijk,1,1) &
                       + wgllwgll_xz_3D(ijk,1,1)*newtempx2(ijk,1,1) + wgllwgll_xy_3D(ijk,1,1)*newtempx3(ijk,1,1))
    enddo
#endif

    ! updates potential_dot_dot_acoustic with contribution from each C-PML element
    if (PML_CONDITIONS .and. (.not. backward_simulation)  .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an ".and." statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if(is_CPML(ispec)) then
#ifndef FORCE_VECTORIZATION
        do k = 1,NGLLZ
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,k,ispec)
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                                                  potential_dot_dot_acoustic_CPML(i,j,k)
            enddo
          enddo
        enddo
#else
        do ijk = 1,NGLLCUBE
          ! sum contributions from each element to the global values
          iglob = ibool(ijk,1,1,ispec)
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_CPML(ijk,1,1)
        enddo
#endif
      endif
    endif


  enddo ! end of loop over all spectral elements

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

  end subroutine compute_forces_acoustic_Dev

