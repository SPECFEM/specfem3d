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
!=====================================================================

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


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
  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
                        m1,m2,NGLLCUBE,PML_CONDITIONS

  use pml_par, only: is_CPML, spec_to_CPML, NSPEC_CPML, &
                     PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl,&
                     PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old,&
                     PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new,&
                     potential_dot_dot_acoustic_CPML,rmemory_dpotential_dxl,rmemory_dpotential_dyl,&
                     rmemory_dpotential_dzl,rmemory_potential_acoustic, &
                     potential_acoustic_old,potential_acoustic_new

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: &
        rhostore,jacobian

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprimewgll_xx,hprimewgll_xxT

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

  integer,intent(in) :: iphase
  integer,intent(in) :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic
  integer, dimension(num_phase_ispec_acoustic,2),intent(in) :: phase_ispec_inner_acoustic

  ! CPML adjoint
  logical,intent(in) :: backward_simulation

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempx2,newtempx3

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: rho_invl

  integer :: ispec,iglob,ispec_p,num_elements

  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem

  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: E1_m1_m2_5points

  equivalence(chi_elem,B1_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(newtempx1,E1_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: E1_mxm_m2_m1_5points

  equivalence(chi_elem,A1_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)

  ! CPML
  integer :: ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_old,tempx2_old,tempx3_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_new,tempx2_new,tempx3_new

  ! CPML Deville
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem_old,chi_elem_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points_old,B1_m1_m2_5points_new
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points_old,C1_m1_m2_5points_new

  equivalence(chi_elem_old,B1_m1_m2_5points_old)
  equivalence(tempx1_old,C1_m1_m2_5points_old)

  equivalence(chi_elem_new,B1_m1_m2_5points_new)
  equivalence(tempx1_new,C1_m1_m2_5points_new)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points_old,A1_mxm_m2_m1_5points_new
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points_old,C1_mxm_m2_m1_5points_new

  equivalence(chi_elem_old,A1_mxm_m2_m1_5points_old)
  equivalence(tempx3_old,C1_mxm_m2_m1_5points_old)

  equivalence(chi_elem_new,A1_mxm_m2_m1_5points_new)
  equivalence(tempx3_new,C1_mxm_m2_m1_5points_new)

  integer :: i,j,k
#ifdef FORCE_VECTORIZATION
! this will (purposely) give out-of-bound array accesses if run through range checking,
! thus use only for production runs with no bound checking
  integer :: ijk
#endif


  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! gets values for element
    DO_LOOP_IJK
      chi_elem(INDEX_IJK) = potential_acoustic(ibool(INDEX_IJK,ispec))
    ENDDO_LOOP_IJK

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

    if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an ".and." statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false

      if (is_CPML(ispec)) then
      ! gets values for element
        DO_LOOP_IJK
          chi_elem_old(INDEX_IJK) = potential_acoustic_old(ibool(INDEX_IJK,ispec))
          chi_elem_new(INDEX_IJK) = potential_acoustic_new(ibool(INDEX_IJK,ispec))
        ENDDO_LOOP_IJK

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

            C1_m1_m2_5points_new(i,j) = hprime_xx(i,1)*B1_m1_m2_5points_new(1,j) + &
                                        hprime_xx(i,2)*B1_m1_m2_5points_new(2,j) + &
                                        hprime_xx(i,3)*B1_m1_m2_5points_new(3,j) + &
                                        hprime_xx(i,4)*B1_m1_m2_5points_new(4,j) + &
                                        hprime_xx(i,5)*B1_m1_m2_5points_new(5,j)
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

              tempx2_new(i,j,k) = chi_elem_new(i,1,k)*hprime_xxT(1,j) + &
                                  chi_elem_new(i,2,k)*hprime_xxT(2,j) + &
                                  chi_elem_new(i,3,k)*hprime_xxT(3,j) + &
                                  chi_elem_new(i,4,k)*hprime_xxT(4,j) + &
                                  chi_elem_new(i,5,k)*hprime_xxT(5,j)
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

            C1_mxm_m2_m1_5points_new(i,j) = A1_mxm_m2_m1_5points_new(i,1)*hprime_xxT(1,j) + &
                                            A1_mxm_m2_m1_5points_new(i,2)*hprime_xxT(2,j) + &
                                            A1_mxm_m2_m1_5points_new(i,3)*hprime_xxT(3,j) + &
                                            A1_mxm_m2_m1_5points_new(i,4)*hprime_xxT(4,j) + &
                                            A1_mxm_m2_m1_5points_new(i,5)*hprime_xxT(5,j)
          enddo
        enddo

      endif ! is_CPML
    endif ! PML_CONDITIONS

    DO_LOOP_IJK
      ! get derivatives of potential with respect to x, y and z
      xixl = xix(INDEX_IJK,ispec)
      xiyl = xiy(INDEX_IJK,ispec)
      xizl = xiz(INDEX_IJK,ispec)
      etaxl = etax(INDEX_IJK,ispec)
      etayl = etay(INDEX_IJK,ispec)
      etazl = etaz(INDEX_IJK,ispec)
      gammaxl = gammax(INDEX_IJK,ispec)
      gammayl = gammay(INDEX_IJK,ispec)
      gammazl = gammaz(INDEX_IJK,ispec)
      jacobianl = jacobian(INDEX_IJK,ispec)

      ! derivatives of potential
      dpotentialdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
      dpotentialdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
      dpotentialdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

      ! stores derivatives of ux, uy and uz with respect to x, y and z
      if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
        ! do not merge this second line with the first using an ".and." statement
        ! because array is_CPML() is unallocated when PML_CONDITIONS is false
        if (is_CPML(ispec)) then
          PML_dpotential_dxl(INDEX_IJK) = dpotentialdxl
          PML_dpotential_dyl(INDEX_IJK) = dpotentialdyl
          PML_dpotential_dzl(INDEX_IJK) = dpotentialdzl

          PML_dpotential_dxl_old(INDEX_IJK) = &
            xixl*tempx1_old(INDEX_IJK) + etaxl*tempx2_old(INDEX_IJK) + gammaxl*tempx3_old(INDEX_IJK)
          PML_dpotential_dyl_old(INDEX_IJK) = &
            xiyl*tempx1_old(INDEX_IJK) + etayl*tempx2_old(INDEX_IJK) + gammayl*tempx3_old(INDEX_IJK)
          PML_dpotential_dzl_old(INDEX_IJK) = &
            xizl*tempx1_old(INDEX_IJK) + etazl*tempx2_old(INDEX_IJK) + gammazl*tempx3_old(INDEX_IJK)

          PML_dpotential_dxl_new(INDEX_IJK) = &
            xixl*tempx1_new(INDEX_IJK) + etaxl*tempx2_new(INDEX_IJK) + gammaxl*tempx3_new(INDEX_IJK)
          PML_dpotential_dyl_new(INDEX_IJK) = &
            xiyl*tempx1_new(INDEX_IJK) + etayl*tempx2_new(INDEX_IJK) + gammayl*tempx3_new(INDEX_IJK)
          PML_dpotential_dzl_new(INDEX_IJK) = &
            xizl*tempx1_new(INDEX_IJK) + etazl*tempx2_new(INDEX_IJK) + gammazl*tempx3_new(INDEX_IJK)
        endif
      endif

      ! density (reciproc)
      rho_invl = 1.0_CUSTOM_REAL / rhostore(INDEX_IJK,ispec)

      ! for acoustic medium
      ! also add GLL integration weights
      tempx1(INDEX_IJK) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
      tempx2(INDEX_IJK) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
      tempx3(INDEX_IJK) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)
    ENDDO_LOOP_IJK

    if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an ".and." statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
        ispec_CPML = spec_to_CPML(ispec)

        ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
        call pml_compute_memory_variables_acoustic(ispec,ispec_CPML,tempx1,tempx2,tempx3,&
                                                   rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl)

        ! calculates contribution from each C-PML element to update acceleration
        call pml_compute_accel_contribution_acoustic(ispec,ispec_CPML,potential_acoustic,&
                                                     potential_dot_acoustic,rmemory_potential_acoustic)
      endif
    endif ! PML_CONDITIONS

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
#ifdef FORCE_VECTORIZATION
! we can force vectorization using a compiler directive here because we know that there is no dependency
! inside a given spectral element, since all the global points of a local elements are different by definition
! (only common points between different elements can be the same)
!DIR$ IVDEP
#endif
    DO_LOOP_IJK
      ! sum contributions from each element to the global values
      iglob = ibool(INDEX_IJK,ispec)
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                          - ( wgllwgll_yz_3D(INDEX_IJK)*newtempx1(INDEX_IJK) &
                                            + wgllwgll_xz_3D(INDEX_IJK)*newtempx2(INDEX_IJK) &
                                            + wgllwgll_xy_3D(INDEX_IJK)*newtempx3(INDEX_IJK))
    ENDDO_LOOP_IJK

    ! updates potential_dot_dot_acoustic with contribution from each C-PML element
    if (PML_CONDITIONS .and. (.not. backward_simulation)  .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an ".and." statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_CPML(INDEX_IJK)
        ENDDO_LOOP_IJK
      endif
    endif ! PML_CONDITIONS

  enddo ! end of loop over all spectral elements

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

  end subroutine compute_forces_acoustic_Dev

